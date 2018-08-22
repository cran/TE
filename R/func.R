#' Estimate TE dynamics using mismatch data
#'
#' Given the number of mismatches and element lengths for an LTR retrotransposon family, estimate the age distribution, insertion rate, and deletion rates. 
#'
#' @param mismatch A vector containing the number of mismatches.
#' @param len A vector containing the length of each element. 
#' @param r Mutation rate (substitutions/(million year * site)) used in the calculation. 
#' @param perturb A scalar multiple to perturb the estimated death rate from the null hypothesis estimate. Used to generate the sensitivity analysis.
#' @param rateRange A vector of death rates, an alternative to \code{perturb} for specifying the death rates.
#' @param plotFit Whether to plot the distribution fits.
#' @param plotSensitivity Whether to plot the sensitivity analysis.
#' @param pause Whether to pause after each plot. 
#' @param main The title for the plot.
#'
#' @details \code{EstDynamics} estimates the TE dynamics through fitting a negative binomial fit to the mismatch data, while \code{EstDynamics2} uses a mixture model. For detailed implementation see References.
#' 
#' @return \code{EstDynamics} returns a \code{TEfit} object, containing the following fields, where the unit for time is million years ago (Mya):
#' \item{pvalue}{The p-value for testing H_0: The insertion rate is uniform over time.}
#' \item{ageDist}{A list containing the estimated age distributions.}
#' \item{insRt}{A list containing the estimated insertion rates.}
#' \item{agePeakLoc}{The maximum point (in age) of the age distribution. }
#' \item{insPeakLoc}{The maximum point (in time) of the insertion rate.}
#' \item{estimates}{The parameter estimates from fitting the distributions; see References}
#' \item{sensitivity}{A list containing the results for the sensitivity analysis, with fields \code{time}: time points; \code{delRateRange}: A vector for the range of deletion rates; \code{insRange}: A matrix whose columns contain the insertion rates under different scenarios.}
#' \item{n}{The sample size.}
#' \item{meanLen}{The mean of element length.}
#' \item{meanDiv}{The mean of divergence.}
#' \item{KDE}{A list containing the kernel density estimate for the mismatch data. }
#' \item{logLik}{The log-likelihoods of the parametric fits.}
#' 
#' @references
#' \cite{Dai, X., Wang, H., Dvorak, J., Bennetzen, J., Mueller, H.-G. (2018). "Birth and Death of LTR Retrotransposons in Aegilops tauschii". Genetics}
#' @export
#' 
#' @examples 
#' # Analyze Gypsy family 24 (Nusif)
#' data(AetLTR)
#' dat <- subset(AetLTR, GroupID == 24 & !is.na(Chr))
#' set.seed(1)
#' res1 <- EstDynamics(dat$Mismatch, dat$UngapedLen, plotFit=TRUE, plotSensitivity=FALSE, pause=FALSE)
#'
#' # p-value for testing a uniform insertion rate
#' res1$pvalue
#' 

EstDynamics <- function(mismatch, len, r=1.3e-2, perturb=2, rateRange=NULL, plotFit=FALSE, plotSensitivity=FALSE, pause=plotFit && plotSensitivity, main=sprintf('n = %d', n)) {
  
  stopifnot(length(mismatch) == length(len))
  stopifnot(all(perturb > 0))
  if (any(is.na(mismatch)) || any(is.na(len))) {
    stop('NA values in the input are not allowed')
  }
  
  n <- length(mismatch)
  perturb[perturb < 1] <- 1 / perturb[perturb < 1]
  
  # mutation rate
  meanLen <- mean(len)

## Fit distributions to the mutation data.
  geom1 <- MASS::fitdistr(mismatch, 'geometric')
  geom1$loglik <- sum(stats::dgeom(mismatch, geom1[['estimate']], log=TRUE))
  nb1 <- suppressWarnings({
    MASS::fitdistr(mismatch, 'negative binomial')
  })
  bw <- stats::bw.nrd0(mismatch)
  bw <- max(bw, 0.5)
  den1 <- stats::density(mismatch, bw, kernel='gaussian')

  logLik <- c(nb1 = nb1[['loglik']], geom1 = geom1[['loglik']])
  logLikDiff <- logLik['nb1'] - logLik['geom1']
  pv <- stats::pchisq(2 * logLikDiff, 1, lower.tail=FALSE)

## Recover age distribution.
  alpha1 <- nb1$estimate['size']
  beta1 <- nb1$estimate['size'] / nb1$estimate['mu']

## Time axis
  maxPts <- floor(max(mismatch) * 1.5)
  pts <- seq(0, maxPts)
  ptsDense <- seq(0.1, maxPts, by=0.1)
  timeScale <- pts / meanLen / r / 2
  timeScaleDense <- ptsDense / meanLen / r / 2

  nbDensity <- stats::dnbinom(pts, nb1$estimate['size'], mu=nb1$estimate['mu'])
  lineGeom <- stats::dgeom(pts, geom1$estimate['prob'])
  lineGammaRec <- stats::dgamma(ptsDense, shape=alpha1, rate=beta1) # TODO: use log scale.
  yMax <- max(max(den1[['y']]), max(nbDensity), max(lineGeom), max(lineGammaRec)) * 1.1

  if (plotFit) {
    lwd <- 2
    MASS::truehist(as.integer(mismatch), 
         xlab='Mismatch', 
         main=main, 
         ylim=c(0, yMax), 
         prob=TRUE, 
         col='white', 
         nbins=min(20, length(unique(as.integer(mismatch))) + 1))
         # , xlim=c(0, max(mismatch)))
    ind <- den1$x >= 0
    lines(den1$x[ind], den1$y[ind], lwd=lwd)
    lines(pts, nbDensity, lty=2, lwd=lwd, col='blue')
    # lines(pts, lineGeom, lty=3)
    lines(ptsDense, lineGammaRec, lty=4, lwd=lwd, col='blue')
    legend('topright', 
           legend=c('KDE', 
                    'neg. bin.', 
                    # 'geometric', 
                    'age dist'), 
           lty=c(1, 2, 4), 
           lwd=lwd,
           col=c('black', 'blue', 'blue')) 
           # title='Density fit')
    if (pause)
      readline('Press enter to continue')
  }

## Recover birth rate. The death rate is estimated from H0 (to be improved).
  p1 <- geom1$estimate['prob']
  deathRate <- p1 / (1 - p1)
  names(deathRate) <- 'rate'
  lineExpRec <- stats::dexp(ptsDense, deathRate)
  # plot(pts, lineExpRec, type='l', xlim=c(0, 40))

  birthRateRec <- lineGammaRec / stats::pexp(ptsDense, deathRate, lower.tail=FALSE)
  # Normalize to unit integral
  birthRateRec <- birthRateRec / sum(birthRateRec) * 
                  length(birthRateRec) / diff(range(timeScaleDense)) 
  ageDistRec <- lineGammaRec / sum(lineGammaRec) * 
                  length(lineGammaRec) / diff(range(timeScaleDense))

## Sensitivity analysis.
  if (is.null(rateRange)) {
    rateRange <- signif(deathRate * c(rev(1 / perturb), 1, perturb), 3)
  } else {
    rateRange <- signif(rateRange, 3)
  }
    
## In mutation scale.
  birthRange <- sapply(rateRange, function(rate) {
    birth <- lineGammaRec / (1 - stats::pexp(ptsDense, rate))
    birth[is.infinite(birth)] <- NaN
    birth <- birth / sum(birth, na.rm=TRUE) / (timeScaleDense[2] - timeScaleDense[1])
    birth
  })
  colnames(birthRange) <- rateRange
  
  obj1 <- rainbow::fts(timeScaleDense, birthRange, xname='Mya', yname='Normalized birth rate')
  ind <- max(which(birthRange > 1e-4, arr.ind=TRUE)[, 1])
  if (ind == -Inf)
    ind <- length(timeScaleDense)
  xMaxSen <- timeScaleDense[ind]
  if (plotSensitivity) {
    rainbow::plot.fds(obj1, plot.type='functions', main=sprintf('Sensitivity analysis'), xlim=c(0, xMaxSen), plotlegend=TRUE)
    if (pause)
      readline('Press enter to continue')
  }

  estimates <- c(nb_n       = unname(nb1$estimate['size']), 
                 nb_p      = unname(nb1$estimate['size'] /
                                    (nb1$estimate['size'] +
                                     nb1$estimate['mu'])), 
                 gamma_alpha= unname(alpha1)         , 
                 gamma_beta = unname(beta1) * 2 * r * meanLen, 
                 geom_p     = unname(geom1$estimate) , 
                 exp_lambda = unname(deathRate) * 2 * r * meanLen,
                 exp_half = log(2) / (unname(deathRate) * 2 * r * meanLen))

  res <- list(
    pvalue = unname(pv), 
    # Gamma mode, in My
    ageDist = list(x =timeScaleDense, y=ageDistRec),
    insRt = list(x=timeScaleDense, y=birthRateRec), 
    agePeakLoc = unname((alpha1 - 1) / beta1 / 2 / r / meanLen), 
    insPeakLoc = timeScaleDense[which.max(birthRateRec)],
    estimates = estimates, 
    sensitivity = list(
      time      = timeScaleDense, 
      delRateRange = rateRange, 
      insRange     = birthRange
    ), 
    n = n, 
    meanLen = meanLen, 
    meanDiv = mean(mismatch / len), 
    KDE = den1, 
    logLik = logLik)
  
  class(res) <- 'TEfit'
  
  res
}


#' Print a TEfit or TEfit2 object
#' 
#' @param x A TEfit or TEfit2 object
#' @param ... Not used
#' 
#' @export
print.TEfit <- function(x, ...) {
  
  cat(sprintf('Fit for LTR-RTN family\n\nn=%s, mean divergence=%.2g, mean length=%.4g\n', x[['n']], x[['meanDiv']], x[['meanLen']]))
  cat(sprintf('P-value for testing uniform insertion rate: %.3g\n', x[['pvalue']]))
  
  invisible(x)
}

#' @rdname print.TEfit
#' @export
print.TEfit2 <- function(x, ...) {
  
  cat('Mixture ')
  NextMethod()
}

#' Generate sensitivity plots 
#' 
#' Create sensitivity plots of a few families to investigate different death rate scenarios
#' 
#' @param resList A list of families returned by \code{\link{EstDynamics}}
#' @param col A vector of colors
#' @param xMax The maximum of the x-axis
#' @param markHalfPeak Whether to mark the time points with half-intensity
#' @param famLegend Whether to create legend for families
#' @param rLegend Text for the legend for families
#' @param ... Passed into \code{matplot}
#' @export
#' @examples
#' data(AetLTR)
#' copia3 <- subset(AetLTR, GroupID == 3 & !is.na(Chr))
#' copia9 <- subset(AetLTR, GroupID == 9 & !is.na(Chr))
#' res3 <- EstDynamics(copia3$Mismatch, copia3$UngapedLen)
#' res9 <- EstDynamics(copia9$Mismatch, copia9$UngapedLen)
#' SensitivityPlot(list(`Copia 3`=res3, `Copia 9`=res9))

SensitivityPlot <- function(resList, col, xMax, markHalfPeak=FALSE, famLegend=TRUE, rLegend=names(resList), ...) {

  if (!is.null(names(resList)[1]) && 'pvalue' %in% names(resList)) { # not a list
    resList <- list(resList) 
  }

  # if (any(sapply(resList, function(x) 'TEfit2' %in% class(x)))) {
    # warning("SensitivityPlot is only implemented for 'TEfit' objects. Using fits produced by EstDynamics")
  # }

  if (missing(col)) {
    col <- grDevices::rainbow(length(resList))
  }
  
  # v: A vector
  # len: The desired length after padding in NAs
  padNA <- function(v, len) {
    if (is.vector(v)) {
      c(v, rep(NA, len - length(v)))
    } else {
      rbind(v, matrix(NA, len - nrow(v), ncol(v)))
    }
  }
  timeList <- lapply(resList, function(res) if (inherits(res, 'TEfit2'))  res[['sensitivity2']][['time']] else res[['sensitivity']][['time']])
  maxl <- max(sapply(timeList, length))
  timeMat <- sapply(timeList, padNA, len=maxl)
  birthRangeAll <- sapply(resList, function(res) {
    if (inherits(res, 'TEfit2')) {
      padNA(res[['sensitivity2']][['insRange']], maxl)
    } else {
      padNA(res[['sensitivity']][['insRange']], maxl)
    }
  }, simplify='array') # Reorder to a 3D array

  if (any(is.nan(birthRangeAll))) {
    warning('Some scenario may produce unreasonable results')
  }

  if (missing(xMax)) {
    xMax <- max(timeMat[which(birthRangeAll > 1e-4, arr.ind=TRUE)[, c(1, 3)]])
  }
  yMax <- max(birthRangeAll, na.rm=TRUE)

  baselwd <- 2
  matplot(timeMat, birthRangeAll[, 2, ], xlab='Mya', lty=1, col=col, type='l', xlim=c(0, xMax), ylim=c(0, yMax), lwd=baselwd, ...)
  matplot(timeMat, birthRangeAll[, 1, ], xlab='Mya', lty=2, col=col, add=TRUE, type='l', lwd=baselwd, ...)
  matplot(timeMat, birthRangeAll[, 3, ], xlab='Mya', lty=3, col=col, add=TRUE, type='l', lwd=baselwd + 2, ...)
  legend('topright', c('low', 'median', 'high'), lty=c(2, 1, 3), title='Death rate', lwd=c(baselwd, baselwd, baselwd + 2))

  # indices are (subject, scenario, before/after, x/y)
  peakLoc <- array(NA, c(dim(birthRangeAll)[3], 3, 2, 2))
  dimnames(peakLoc) <- list(names(resList), NULL, c('before', 'after'), c('x', 'y'))

  for (k in 1:dim(birthRangeAll)[3]) {
    halfInd <- apply(birthRangeAll[, , k], 2, findHalfIntensity)

    for (i in 1:2) {
      xx <- timeMat[halfInd[i, ], k]
      yy <- halfInd[3, ]
      peakLoc[k, , i, ] <- cbind(xx, yy)
      if (markHalfPeak) {
        matlines(matrix(c(xx - 0.09,
                          xx + 0.09), nrow=2, 
                        byrow=TRUE), 
                 rbind(yy, yy),
                 pch='-', lwd=c(baselwd, baselwd, baselwd + 2), col=col[k], 
                 lty=c(2, 1, 3))
      }
    } 
  }

  if (length(resList) != 1 && !is.null(col) && famLegend) {
    legend('right', rLegend, col=col, lty=1, lwd=baselwd)
  }

  invisible(list(timeMat=timeMat, 
                 birthRange=birthRangeAll,
                 peakLoc=peakLoc))
}


findHalfIntensity <- function(y) {

  peakLoc <- which.max(y)
  peak <- y[peakLoc]
  half <- peak / 2
  beforePeak <- seq_len(peakLoc)
  afterPeak <- seq(peakLoc, length(y), by=1)
  beforex <- beforePeak[which.min(abs(y[beforePeak] - half))]
  afterx <- afterPeak[which.min(abs(y[afterPeak] - half))]

  c(beforex=beforex, afterx=afterx, halfPeak=half)
}


#' Calcualte the KL divergence of a negative binomial fit to the mismatch data.
#' 
#' @param res A TEfit object.
#'
#'
#' @export
#' @examples
#' # Analyze Gypsy family 24 (Nusif)
#' data(AetLTR)
#' dat <- subset(AetLTR, GroupID == 24 & !is.na(Chr))
#' set.seed(1)
#' res1 <- EstDynamics(dat$Mismatch, dat$UngapedLen, plotFit=TRUE, plotSensitivity=FALSE, pause=FALSE)
#' nbLackOfFitKL(res1)
nbLackOfFitKL <- function(res) {

  pts <- seq(0, 0.2 * res$meanLen) # TODO: modify this

  if ('TEfit2' %in% class(res)) {
    est2 <- res[['estimates2']]
    nbDensity <- dTwoNB(pts, 
                        size1=est2[['size1']],
                        prob1=est2[['prob1']],
                        size2=est2[['size2']],
                        prob2=est2[['prob2']], 
                        p=est2[['p']])
  } else {
    est <- res[['estimates']]
    nbDensity <- stats::dnbinom(pts, est['nb_n'], est['nb_p'])
  }
  
  KDE <- res[['KDE']]
  KDEval <- stats::approx(KDE[['x']], KDE[['y']], pts)[['y']]
  KDEval <- KDEval / sum(KDEval, na.rm=TRUE)
  kl <- sum(nbDensity * log(nbDensity / KDEval), na.rm=TRUE)
  # browser()

  kl
}


dTwoNB <- function(x, size1, prob1, size2, prob2, p) {
  stopifnot(p >= 0 && p <= 1)
  stats::dnbinom(x, size=size1, prob=prob1) * p + 
    stats::dnbinom(x, size=size2, prob=prob2) * (1 - p)
}

#' @param nTrial The number of starting points for searching for the MLE. 
#' @param ... Pass to EstDynamics 
#' @return 
#' This function returns a \code{TEfit2} object, containing all the above fields for \code{TEfit} and the following:
#' \item{estimates2}{The parameter estimates from fitting the mixture distribution.}
#' \item{ageDist2}{The estimated age distribution from fitting the mixture distribution.} 
#' \item{insRt2}{The estimated insertion rate from fitting the mixture distribution.} 
#' \item{agePeakLoc2}{Maximum point(s) for the age distribution.}
#' \item{insPeakLoc2}{Maximum point(s) for the insertion rate.}
#' 
#' @rdname EstDynamics
#' @export
#'
#' @examples
#'
#' # Use a mixture distribution to improve fit
#' res2 <- EstDynamics2(dat$Mismatch, dat$UngapedLen, plotFit=TRUE)
#'
#' # A larger number of trials is recommended to achieve the global MLE
#' \dontrun{
#' res3 <- EstDynamics2(dat$Mismatch, dat$UngapedLen, plotFit=TRUE, nTrial=1000L)
#' }

EstDynamics2 <- function(mismatch, len, r=1.3e-2, nTrial=10L, perturb=2, rateRange=NULL, plotFit=FALSE, plotSensitivity=FALSE, pause=plotFit && plotSensitivity, ...) {

  ddd <- list(...)
  main <- ddd[['main']]
  meanLen <- mean(len)

  # plotFit <- ddd[['plotFit']]
  # if (is.null(plotFit)) {
    # plotFit <- FALSE
  # }
  
  res_ <- EstDynamics(mismatch, len, plotFit=FALSE, plotSensitivity=FALSE, pause=FALSE, ...)
  size0 <- res_[['estimates']]['nb_n']
  prob0 <- res_[['estimates']]['nb_p']
  
  if (nTrial <= 0)
    stop('need at least 1 trial')
    
  loglik <- -Inf
  res <- NULL
  for (i in seq_len(nTrial)) {
    start <- list(size1 = size0 * stats::runif(1, 0.5, 3), 
                  prob1 = prob0 * stats::runif(1, 0.5, 1.5),
                  size2 = size0 * stats::runif(1, 0.5, 3), 
                  prob2 = prob0 * stats::runif(1, 0.5, 1.5),
                  p = stats::runif(1, 0.1, 0.9))
    # In case the starting value does not work.
    tmpRes <- suppressWarnings({
      tryCatch(
        MASS::fitdistr(mismatch, dTwoNB, start=start, 
                       lower=c(0, 0, 0, 0, 0), upper=c(Inf, 1, Inf, 1, 1), 
                       control=list(maxit=500L)),
        error=function(e) NULL
      )
    })
    if (!is.null(tmpRes) && stats::logLik(tmpRes) > loglik) {
      loglik <- unname(stats::logLik(tmpRes))
      res <- tmpRes
    }
  }
  
  if (is.null(res))
    stop('No starting point works')

  est2 <- res[['estimate']]
  mu1 <- est2['size1'] * (1 - est2['prob1']) / est2['prob1']
  alpha1 <- est2['size1']
  beta1 <- est2['size1'] / mu1
  mu2 <- est2['size2'] * (1 - est2['prob2']) / est2['prob2']
  alpha2 <- est2['size2']
  beta2 <- est2['size2'] / mu2

  lof <- nbLackOfFitKL(res_)
  maxMis <- floor(max(mismatch) * 1.5)
  pts <- seq(0, maxMis)
  ptsDense <- seq(0.1, maxMis, by=0.1)
  timeScale <- pts / meanLen / r / 2
  timeScaleDense <- ptsDense / meanLen / r / 2

  den1 <- res_[['KDE']]
  lineNB <- dnbinom(pts, 
                    res_[['estimates']]['nb_n'], 
                    res_[['estimates']]['nb_p'])
  lineGammaRec <- res_[['ageDist']][['y']]
  lineGammaRec <- lineGammaRec * max(res_[['ageDist']][['x']]) / maxMis
  lineNB2 <- dTwoNB(pts, size1=est2['size1'], prob1=est2['prob1'], 
                    size2=est2['size2'], prob2=est2['prob2'], p=est2['p'])
  lineGammaRec2 <- stats::dgamma(ptsDense, alpha1, beta1) * est2['p'] + stats::dgamma(ptsDense, alpha2, beta2) * (1 - est2['p'])
  
  allLines <- c(den1[['y']], lineNB, lineGammaRec, lineNB2, lineGammaRec2)
  yMax <- max(allLines[!is.infinite(allLines)])
    
  if (plotFit) {
    if (is.null(main)) {
      main <- sprintf('n = %d, LoF = %.3f', 
                      length(mismatch), lof)
    }
    hist(mismatch, xlab='Mismatch', 
         main=main,
         ylim=c(0, yMax), probability=TRUE, col='white')
    lwd <- 2
    ind <- den1$x >= 0
    lines(den1$x[ind], den1$y[ind], lwd=lwd)
    lines(pts, lineNB, lty=2, lwd=lwd, col='blue')
    lg <- c('KDE', 'neg. bin.', 'age dist.', '2 neg. bin.', '2 age dist.')
    lty <- c(1, 2, 4, 2, 4)
    col <- c('black', rep('blue', 2), rep('red', 2))
    legend('topright', legend=lg, col=col, lty=lty, lwd=lwd)
    lines(ptsDense, lineGammaRec, lty=4, lwd=lwd, col='blue')
  
    lines(pts, lineNB2, lty=2, col='red', lwd=lwd)
    lines(ptsDense, lineGammaRec2, lty=4, col='red', lwd=lwd)
  }

  ## Sensitivity analysis.
  p1 <- res_$estimate['geom_p']
  deathRate <- p1 / (1 - p1)
  if (is.null(rateRange)) {
    rateRange <- signif(deathRate * c(rev(1 / perturb), 1, perturb), 3)
  } else {
    rateRange <- signif(rateRange, 3)
  }
    
  ## In mutation scale.
  birthRange2 <- sapply(rateRange, function(rate) {
    birth <- lineGammaRec2 / (1 - stats::pexp(ptsDense, rate))
    birth[is.infinite(birth)] <- NaN
    birth <- birth / sum(birth, na.rm=TRUE) / (timeScaleDense[2] - timeScaleDense[1])
    birth
  })
  colnames(birthRange2) <- rateRange
  
  obj1 <- rainbow::fts(timeScaleDense, birthRange2, xname='Mya', yname='Normalized birth rate')
  ind <- max(which(birthRange2 > 1e-4, arr.ind=TRUE)[, 1])
  if (ind == -Inf) {
    ind <- length(timeScaleDense)
  }
  xMaxSen <- timeScaleDense[ind]
  if (plotSensitivity) {
    rainbow::plot.fds(obj1, plot.type='functions', main=sprintf('Sensitivity analysis'), xlim=c(0, xMaxSen), plotlegend=TRUE)
    if (pause) {
      readline('Press enter to continue')
    }
  }
    
  estimates2 <- c(
    alpha1 = unname(alpha1), 
    beta1 = unname(beta1 * 2 * r * res_[['meanLen']]), 
    size1 = unname(est2['size1']), 
    prob1 = unname(est2['prob1']), 
    alpha2 = unname(alpha2), 
    beta2 = unname(beta2 * 2 * r * res_[['meanLen']]), 
    size2 = unname(est2['size2']), 
    prob2 = unname(est2['prob2']),
    p = unname(est2['p']), 
    loglik = loglik)
  
  timeScaleDense <- ptsDense / res_[['meanLen']] / r / 2
  birthRateRec2 <- lineGammaRec2 / stats::pexp(timeScaleDense, res_[['estimates']]['exp_lambda'], lower.tail=FALSE)
  # Normalize to unit integral
  birthRateRec2 <- birthRateRec2 / sum(birthRateRec2) * 
                  length(birthRateRec2) / diff(range(timeScaleDense)) 
  ageDist2Rec <- lineGammaRec2 / sum(lineGammaRec2) * 
                   length(lineGammaRec2) / diff(range(timeScaleDense))

  res <- c(
    list(
      estimates2 = estimates2, 
      ageDist2 = list(x=timeScaleDense, y=ageDist2Rec), 
      insRt2 = list(x=timeScaleDense, y=birthRateRec2), 
      agePeakLoc2 = timeScaleDense[findPeaks(lineGammaRec2)], 
      insPeakLoc2 = timeScaleDense[findPeaks(birthRateRec2)],
      sensitivity2 = list(
        time      = timeScaleDense, 
        delRateRange = rateRange, 
        insRange     = birthRange2
      )
    ), 
    res_)
  class(res) <- c('TEfit2', class(res_))
  res
}


# # Get the recovered (gamma) age distribution, possibly with mixtures.
# getGammaCurve <- function(est, pts=seq(0, 10, length.out=1000L)) {
  # if (all(c('alpha1', 'alpha2', 'beta1', 'beta2', 'p') %in% 
                 # names(est))) {
  # # Two mixture case
    # alpha1 <- est['alpha1']
    # beta1 <- est['beta1']
    # alpha2 <- est['alpha2']
    # beta2 <- est['beta2']
    # p <- est['p']
    
    # den <- stats::dgamma(pts, alpha1, beta1) * p +
           # stats::dgamma(pts, alpha2, beta2) * (1 - p)
    # peaks <- c((alpha1 - 1) / beta1, (alpha2 - 1) / beta2)
  # } else if (all(c('gamma_alpha', 'gamma_beta') %in% names(est))) {
  # # One mixture case; use Mya scale
    # alpha <- est['gamma_alpha']
    # beta <- est['gamma_beta']
    # den <- stats::dgamma(pts, alpha, beta)
    # peaks <- (alpha - 1) / beta
  # } else
    # stop('`est` does not have the required fields')
  
  # list(x=pts, y=den, peaks=peaks)
# }


# #' Find the indices of peaks (local maximas) for y.
# #'
# #' Note that each function may have multiple peaks/local maximas, representing multiple bursts of insertions.
# #'
# #' @param y A vector containing function values at a set of increasing time points
# #' @export
findPeaks <- function(y) {
  yl <- c(-Inf, y[-length(y)])
  yr <- c(y[-1], -Inf)
  
  which(y > yl + 1e-12 & y > yr + 1e-12)
}

#' Plot the age distributions or insertion rates for multiple families.
#' 
#' @param resList A list of TEfit/TEfit2 objects, which can be mixed
#' @param type Whether to plot the insertion rates ('insRt') or the age distributions ('ageDist').
#' @param ... Passed into plotting functions.
#' 
#' @return A list of line data (plotDat) and peak locations (peakDat).
#' @export
#' @examples
#' data(AetLTR)
#' copia3 <- subset(AetLTR, GroupID == 3 & !is.na(Chr))
#' gypsy24 <- subset(AetLTR, GroupID == 24 & !is.na(Chr))
#' res3 <- EstDynamics(copia3$Mismatch, copia3$UngapedLen)
#' res24 <- EstDynamics2(gypsy24$Mismatch, gypsy24$UngapedLen)
#' 
#' # Plot insertion rates
#' PlotFamilies(list(`Copia 3`=res3, `Gypsy 24`=res24))
#' 
#' # Plot age distributions
#' PlotFamilies(list(`Copia 3`=res3, `Gypsy 24`=res24), type='ageDist')
PlotFamilies <- function(resList, type=c('insRt', 'ageDist'), ...) {

  if (is.null(names(resList))) {
    names(resList) <- paste('Group', seq_along(resList))
  }
  type <- match.arg(type)

  titleText <- ifelse(type == 'insRt', 'Insertion Rate', 'Age Distribution')
  ylab <- ifelse(type == 'insRt', 'Normalized Intensity', 'Density')
  pars <- list(lty = 1, pch=17, xlab = 'Mya', ylab = ylab, main = titleText)
  ddd <- list(...)
  pars[names(ddd)] <- ddd
  

  name1 <- type
  name2 <- paste0(type, '2')
  plotDat <- lapply(seq_along(resList), function(i) {

    res <- resList[[i]]
    mixture <- 'TEfit2' %in% class(res)

    if (mixture) {
      dat <- as.data.frame(res[[name2]])
    } else {
      dat <- as.data.frame(res[[name1]])
    }

    # cbind(dat, 
          # GroupID=names(resList)[i], 
          # mix=mixture, 
          # stringsAsFactors=FALSE)
    dat
  })


  pname1 <- sprintf('%sPeakLoc', substr(type, 1, 3))
  pname2 <- paste0(pname1, '2')
  peakDat <- unlist(sapply(seq_along(resList), function(i) {

    res <- resList[[i]]
    mixture <- 'TEfit2' %in% class(res)

    if (mixture) {
      # dat <- data.frame(x=res[[pname2]])
      res[[pname2]]
    } else {
      # dat <- data.frame(x=res[[pname1]])
      res[[pname1]]
    }

    # cbind(dat, 
          # GroupID=names(resList)[i], 
          # mix=mixture, 
          # stringsAsFactors=FALSE)
  }))

  n <- length(resList)
  maxLen <- max(sapply(plotDat, nrow))
  xMat <- yMat <- matrix(NA, maxLen, n)
  for (i in seq_len(n)) {
    xMat[seq_len(nrow(plotDat[[i]])), i] <- plotDat[[i]][['x']]
    yMat[seq_len(nrow(plotDat[[i]])), i] <- plotDat[[i]][['y']]
  }

  do.call(matplot, c(list(xMat, yMat, type='l'), pars))
  do.call(points, c(list(peakDat, rep(0, length(peakDat))), pars))

  # firstLast <- peakDat[peakDat[['x']] %in% range(peakDat[['x']]), ]


  # p <- ggplot2::ggplot(plotDat, aes(x=x, y=y, color=GroupID)) + settings + ggplot2::geom_point(aes(x=x, y=y), data=cbind(peakDat, y=0), size=1, shape=17, color='grey40')
  invisible(list(plotDat = plotDat, peakDat = peakDat))
}



#' Implements the matrix model in Promislow et al (1999)
#'
#' @param mismatch A vector containing the number of mismatches.
#' @param len A vector containing the length of each element. 
#' @param nsolo An integer giving the number of solo elements.
#' @param r Mutation rate (substitutions/(million year * site)) used in the calculation. 
#' @param plotFit Whether to plot the distribution fits.
#' @param main The title for the plot.
#'
#' @details For the method implemented see References.
#' 
#' @return This function returns various parameter estimates described in Promislow et al. (1999), containing the following fields. The unit for time is million years ago (Mya):
#' \item{B}{The constant insertion rate}
#' \item{q}{The constant excision rate}
#' \item{lam}{The population growth rate}
#' \item{R}{The ratio of the number of elements in class j over class j+1, which is constant by assumption}
#' \item{age1}{The age of the system under model 1 (lambda > 1)}
#' \item{age2}{The age of the system under model 2 (an initial burst followed by stasis lambda = 1)}
#' 
#' @references
#' \cite{Promislow, D., Jordan, K. and McDonald, J. "Genomic demography: a life-history analysis of transposable element evolution." Proceedings of the Royal Society of London B: Biological Sciences 266, no. 1428 (1999): 1555-1560.}
#' @export
#' 
#' @examples 
#' # Analyze Gypsy family 24 (Nusif)
#' data(AetLTR)
#' dat <- subset(AetLTR, GroupID == 24 & !is.na(Chr))
#' res1 <- MatrixModel(dat$Mismatch, dat$UngapedLen, nsolo=450, plotFit=TRUE)
MatrixModel <- function(mismatch, len, nsolo, r=1.3e-2, plotFit=FALSE, main=sprintf('n = %d', n)) {
  # dat <- filter(allDat, GroupID == 1)
               # if (dat$GroupID[1] == 31) browser()
  meanLen <- mean(len)
  theta <- nsolo
  n <- length(mismatch)
  EN <- mean(mismatch)
  mu <- r * 2 * mean(len)
  B <- mu / EN
  q <- theta / (theta + n) / EN * mu
  lam <- 1 + B - q
  gamma_a <- log(n) / log(lam)

  # Model 2, burst + stasis. Estimate R using all the complete elements. Not reliable 
  R <- EN / (EN + 1)
  gamma_b <- theta * R / (n * mu * (1 - R))


  if (plotFit) {
    maxPts <- floor(max(mismatch) * 1.5)
    pts <- seq(0, maxPts)
    ptsDense <- seq(0.1, maxPts, by=0.1)
    timeScale <- pts / meanLen / r / 2
    timeScaleDense <- ptsDense / meanLen / r / 2

    MASS::truehist(mismatch / len / r / 2, 
                   xlab='Age (Mya)', 
                   main=main, 
                   # ylim=c(0, yMax), 
                   prob=TRUE, 
                   col='white', 
                   nbins=min(20, length(unique(as.integer(mismatch))) + 1))
    lines(timeScaleDense, stats::dexp(timeScaleDense, B), lwd=2)
  }

  c(B=B, q=q, lam=lam, R=R, age1=gamma_a, age2=gamma_b)
}


#' Implements the master gene model in Marchani et al (2009)
#'
#' @param mismatch A vector containing the number of mismatches.
#' @param len A vector containing the length of each element. 
#' @param r Mutation rate (substitutions/(million year * site)) used in the calculation. 
#' @param plotFit Whether to plot the distribution fits.
#' @param main The title for the plot.
#'
#' @details For the method implemented see References.
#' 
#' @return This function returns various parameter estimates described in Marchani et al (2009), containing the following fields. The unit for time is million years ago (mya):
#' \item{B}{The constant insertion rate}
#' \item{q}{The constant excision rate}
#' \item{lam}{The population growth rate}
#' \item{R}{The ratio of the number of elements in class j over class j+1, which is constant by assumption}
#' \item{age1}{The age of the system under model 1 (lambda > 1)}
#' \item{age2}{The age of the system under model 2 (an initial burst followed by stasis lambda = 1)}
#' 
#' @references
#' \cite{Marchani, Elizabeth E., Jinchuan Xing, David J. Witherspoon, Lynn B. Jorde, and Alan R. Rogers. "Estimating the age of retrotransposon subfamilies using maximum likelihood." Genomics 94, no. 1 (2009): 78-82.}
#' @export
#' 
#' @examples 
#' # Analyze Gypsy family 24 (Nusif)
#' data(AetLTR)
#' dat <- subset(AetLTR, GroupID == 24 & !is.na(Chr))
#' res2 <- MasterGene(dat$Mismatch, dat$UngapedLen, plotFit=TRUE)
MasterGene <- function(mismatch, len, r=1.3e-2, plotFit=FALSE, main=sprintf('n = %d', n)) {

  # ad hoc method: mean age * 2
  meanLen <- mean(len)
  n <- length(mismatch)
  x <- mismatch
  EN <- mean(x)
  mu <- 2 * r * meanLen
  adhoc <- EN / mu * 2

  # MLE which follows formula (3)
  ff <- function(x, T) {
    exp(-meanLen * T) * (meanLen * T)^x
  }

  dLdT <- function(T) {
    - n / T + sum(vapply(seq_len(n), function(i) {
        ff(x[i], T) / stats::integrate(ff, 0, T, x=x[i])$value
      }, 0))
  }

  root <- stats::uniroot(dLdT, c(1e-4, 4 * adhoc * r), tol=.Machine$double.eps^0.5)
  famAge <- root$root / (2 * r)

  if (plotFit) {
    maxPts <- floor(max(mismatch) * 1.5)
    pts <- seq(0, maxPts)
    ptsDense <- seq(0.1, maxPts, by=0.1)
    timeScale <- pts / meanLen / r / 2
    timeScaleDense <- ptsDense / meanLen / r / 2

    MASS::truehist(mismatch / len / r / 2, 
                   xlab='Age (Mya)', 
                   main=main, 
                   # ylim=c(0, yMax), 
                   prob=TRUE, 
                   col='white', 
                   nbins=min(20, length(unique(as.integer(mismatch))) + 1))
    dens <- stats::dunif(timeScaleDense, 0, famAge)
    is.na(dens) <- 0
    lines(timeScaleDense, dens, lwd=2)
  }

  c(adhoc=adhoc, MLE=famAge)
}



