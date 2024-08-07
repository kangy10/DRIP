\name{threeStageParSel}
\alias{threeStageParSel}
\title{Parameter Selection in Three-Stage Image Restoration}
\description{
    Select the bandwidth value for the image restoration method
    implemented in the function \code{\link{threeStage}}
}
\usage{
threeStageParSel(image, bandwidth, edge1, edge2, nboot,
    blur = FALSE)
}
\arguments{
    \item{image}{A square matrix, no  missing value allowed.}
    \item{bandwidth}{Bandwidth values to be chosen from. Each of
    these values need to be an positive integer specifying the
    number of pixels to use in the local smoothing.}
    \item{edge1}{A matrix of 0 and 1 representing the detected
    step edge pixels.}
    \item{edge2}{A matrix of 0 and 1 representing the detected
    roof/valley edge pixels.}
    \item{nboot}{Required when \code{blur} is TRUE. Unused when
    \code{blur} is FALSE. It must be a positive integer
    specifying the number of bootstraps to perform. See
    Qiu and Kang (2015) for details.}
    \item{blur}{TRUE if the image contains blur, FALSE otherwise.
    If TRUE, the hybrid selection method proposed in
    Qiu and Kang (2015) is used. If FALSE, the leave-one-out cross
    validation is used.} 
}
\value{
    A list of the selected bandwidth, and a matrix of
    cross-validation scores with each entry corresponding to a
    choice of bandwidth.}
\examples{
\dontrun{
step.edges <- stepEdge(peppers, bandwidth = 9, thresh = 17) 
roof.edges <- roofEdge(peppers, bandwidth = 6, thresh = 3000,
    edge1 = step.edges) 
set.seed(24)
# Time consuming
parSel <- threeStageParSel(image = peppers, edge1 = step.edges,
    edge2 = roof.edges, bandwidth = 3:4, nboot = 1, blur = TRUE) 
}
}
\author{Yicheng Kang}
\references{
Qiu, P. and Kang, Y. (2015) ``Blind Image Deblurring Using Jump 
Regression Analysis'', \emph{Statistica Sinica}, \bold{25},
879 -- 899, \doi{10.5705/ss.2014.054}.
}