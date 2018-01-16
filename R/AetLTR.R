#' LTR retrotransposons in \emph{Aegilops tauschii}
#'
#' This data file contains the LTR retrotransposons in \emph{Ae. tauschii}.
#'
#' @name AetLTR
#' @docType data
#' @format A data frame with 18024 rows and 12 columns. Each row corresponds to a unique LTR retrotransposon, and each column corresponds to a feature of the LTR-RTN. The columns are:
#' \describe{
#' \item{SeqID}{LTR retrotransposon sequence ID}
#' \item{UngapedLen}{Length of each LTR}
#' \item{Mismatch}{Number of mismatches}
#' \item{Distance}{Divergence, as defined by (# of mismatches) / (LTR length)}
#' \item{Chr}{Chromosome number}
#' \item{Start}{Start location in bp}
#' \item{Stop}{Ending location in bp}
#' \item{GroupID}{LTR retrotransposon Family ID}
#' \item{sup}{Super family membership}
#' \item{distToGene}{Distance to the nearest gene in bp}
#' \item{logDist}{Log distance to the nearest gene in bp}
#' \item{recRt5}{Recombination rate}
#' } 
#' @references
#' \cite{Luo, Ming-Cheng, et al. (2017) "Genome sequence of the progenitor of the wheat D genome Aegilops tauschii." Nature 551.7681.}
#'
#' \cite{Dai, X., Wang, H., Dvorak, J., Bennetzen, J., Mueller, H.-G. (2018+). "Birth and Death of LTR Retrotransposons in Aegilops tauschii"}
NULL

