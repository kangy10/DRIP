# This is R source code for function 'threeStage', in the
# R package "image".
# Date: May 7, 2013
# Creator: Yicheng Kang

threeStage <- function(image, bandwidth, edge1, edge2,
                       blur = FALSE, plot = FALSE){
  if (!is.matrix(image) || !is.matrix(edge1) || !is.matrix(edge2)) {
    stop("image, edge1 and edge2 must be a matrix")
  } else {
    n1 <- dim(image)[1]
    n2 <- dim(image)[2]
  }
  if (n1 != n2)
    stop("image data must be a square matrix")
  if (!is.numeric(bandwidth))
    stop("bandwidth must be numeric")
  if (length(bandwidth) > 1) {
    stop('bandwidth must be a positive integer')
  }
#   if (n1 + 2 * bandwidth + 2 > 600)
#     stop("some choice of bandwidth or the resolution of the
# image is too large")
  n1 <- dim(image)[1]
  z <- matrix(as.double(image), ncol = n1)
  k <- as.integer(bandwidth)
  if (nrow(edge1) != n1 || ncol(edge1) != n1) {
    stop("different size in edge1 and image")
  }
  if (nrow(edge2) != n1 || ncol(edge2) != n1) {
    stop("different size in edge2 and image")
  }
  if(!all(edge1 == 0 | edge1 == 1))
    stop("edge1 must be either 0 or 1.")
  if(!all(edge2 == 0 | edge2 == 1))
    stop("edge2 must be either 0 or 1.")
  edge1 <- matrix(as.integer(edge1), ncol = n1)
  edge2 <- matrix(as.integer(edge2), ncol = n2)
  if (blur == TRUE) {
    out <- .Fortran(C_deblur_3stage, n = as.integer(n1 - 1),
                    obsImg = z, bandwidth = k, edge1 = edge1, edge2 = edge2,
                    estImg = z)
  }
  else {
    out <- .Fortran(C_denoise_3stage, n = as.integer(n1 - 1),
                    obsImg = z, bandwidth = k, edge1 = edge1, edge2 = edge2,
                    estImg = z)
  }

  if (plot == FALSE) {
    return(out$estImg)
  }
  else { x <- seq(0, 1, length = n1); y <- x
  image(x, y, out$estImg, zlim = c(0, 255),
        col = gray((0:255)/255))
  return(out$estImg)
  }
}
