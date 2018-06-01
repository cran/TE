#' TE: Insertion/Deletion Dynamics for Transposable Elements
#'
#' TE package for analyzing insertion/deletion dynamics for transposable elements
#' 
#' Provides functions to estimate the insertion and deletion rates of
#' transposable element (TE) families. The estimation of insertion rate
#' consists of an improved estimate of the age distribution that takes into
#' account random mutations, and an adjustment by the deletion rate. This
#' package includes functions \code{EstDynamics} and \code{EstDynamics2} for
#' analyzing the TE divergence, and visualization functions such as
#' \code{PlotFamilies} and \code{SensitivityPlot}. 
#' This package implements the methods proposed in Dai et al (2018+).
#' @references
#' \cite{Luo, Ming-Cheng, et al. (2017) "Genome sequence of the progenitor of the wheat D genome Aegilops tauschii." Nature 551.7681.}
#'
#' \cite{Dai, X., Wang, H., Dvorak, J., Bennetzen, J., Mueller, H.-G. (2018+, submitted). "Birth and Death of LTR Retrotransposons in Aegilops tauschii"}
#'
#' 
#' @author
#' Xiongtao Dai \email{dai@@ucdavis.edu},
#' Hao Wang
#' Jan Dvorak
#' Jeffrey Bennetzen
#' Hans-Georg Mueller
#'
#' Maintainer: Xiongtao Dai \email{dai@@ucdavis.edu}
#'
#' @docType package
#' @name TE
#' @importFrom rainbow fts plot.fds 
#' @importFrom graphics hist legend lines matlines matplot points
#' @importFrom stats dgamma dnbinom dexp dunif integrate uniroot
#' @importFrom MASS fitdistr
NULL

