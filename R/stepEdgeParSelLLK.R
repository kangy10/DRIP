# This is R source code for function 'stepEdgeParSelLLK' in
# the R package 'image'.
# Creator: Yicheng Kang
# Date: April 29, 2013

stepEdgeParSelLLK <- function(image, bandwidth, thresh, nboot) {
  if (!is.matrix(image))
    stop("image data must be a matrix")
  n1 <- as.integer(dim(image)[1])
  n2 <- as.integer(dim(image)[2])
  if (n1 != n2)
    stop("image data must be a square matrix")
  if (!is.numeric(bandwidth))
    stop("bandwidth must be numeric")
  if (min(bandwidth) < 1)
    stop("bandwidth must be a positive integer")
#   if (n1 + 2 * max(bandwidth) + 2 > 600)
#     stop("some bandwidth is too large or the resolution
# of the image is too high.")
  if (!is.numeric(thresh))
    stop("threshold must be numeric")
  if ((!is.numeric(nboot)) || (length(nboot) > 1) ||
      (as.integer(nboot) < 1))
    stop('nboot must be a positive integer.')
  bandwidth <- as.integer(bandwidth)
  nboot <- as.integer(nboot)
  nband <- length(bandwidth)
  n <- dim(image)[1]
  z <- matrix(as.double(image), ncol = n)
  u <- as.double(thresh)
  nthresh <- length(u)
  out.mat <- array(0, c(nband, nthresh))
  jp.llk <- JPLLK_surface(image = z, bandwidth = 2:7, plot = FALSE)
  fitted <- jp.llk$fitted
  resid <- jp.llk$resid
  imgs.boot <- array(0, c(n, n, nboot))
  edge.boot <- imgs.boot
  for (iboot in 1:nboot) {
    imgs.boot[ , , iboot] <- fitted + matrix(sample(resid, n^2,
                                                    replace = TRUE), ncol = n)
  }
  for (iband in 1:nband) {
    k <- bandwidth[iband]
    metric <- array(0, c(nthresh, nboot))
    foo.boot <- function(iboot) {
      diff.iboot <- diffLLK(imgs.boot[ , , iboot], k)
      return(diff.iboot)
    }
    diff <- diffLLK(z, k)
    #diff.boot <- array(sapply(1:nboot, foo.boot), c(n, n, nboot))
    diff.boot <- vapply(1:nboot, foo.boot, matrix(0, n, n))
    for (ithresh in 1:nthresh){
      edge <- array(as.integer(0), c(n, n))
      edge[diff < u[ithresh]] <- as.integer(0)
      edge[diff >= u[ithresh]] <- as.integer(1)
      edge.boot[diff.boot < u[ithresh]] <- as.integer(0)
      edge.boot[diff.boot >= u[ithresh]] <- as.integer(1)
      foo.dKQ <- function(iboot){
        return(dKQ(edge, edge.boot[ , , iboot]))
      }
      #metric[ithresh, ] <- as.numeric(sapply(1:nboot, foo.dKQ))
      metric[ithresh, ] <- vapply(1:nboot, foo.dKQ, 0)
    }
    dKQ.boot <- apply(metric, 1, mean)
    out.mat[iband, ] <- dKQ.boot
  }
  rownames(out.mat) <- paste('bandwidth=', bandwidth, sep='')
  colnames(out.mat) <- paste('thresh=', thresh, sep='')
  out.mat.min <- min(out.mat)
  for(iband in 1:nband){
    for(ithresh in 1:nthresh){
      if(out.mat[iband, ithresh] == out.mat.min){
        band_sel <- bandwidth[iband]
        thresh_sel <- thresh[ithresh]
      }
    }
  }
  paste('The selected bandwidth is', band_sel)
  paste('The selected threshold is', thresh_sel)
  return(list(output_matrix=out.mat, selected_bandwidth=band_sel,
              selected_threshold=thresh_sel))
}
