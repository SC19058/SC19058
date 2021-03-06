#' calculate the energy statistics from balanced designed samples
#'
#' This fucntion calculate the statistics to test multrivariate test for equal distributions based on energy method
#' @param size the sample size
#' @param group the number of the groups of the samples
#' @param index for permutation test, indicate how to divide the samples into groups
#' @param distance the distance matrix
#' @return energy statistics used for multivariate test for equal distribution
#' @examples
#' \dontrun{ 
#' pdim <- 10;delta <- 0.5;size = 30;group = 4
#' dat <- c(dat <- c(rt(n = size*pdim,df = 4,ncp = delta),
#'                   rt(n = size*pdim*(group - 1),df = 4))
#' dat <- matrix(dat,byrow = TRUE,nrow = total.size)
#' distance <- as.matrix(dist(dat,diag = TRUE,upper = TRUE,
#'                            method = "minkowski",p = 1))
#' index <- 1:size*group
#' d_stat(size,group,index,distance)
#' }
#' @export
d_stat <- function(size,group,index,distance){
  total <- size * group
  ind <- index;ind <- matrix(ind,byrow = TRUE,nrow = group)
  add.w <- numeric(group)
  for (i in 1:group) {
    add.w[i] <- size/2 * mean(distance[ind[i,],ind[i,]])
  }
  w <- sum(add.w)
  s <- numeric(1)
  for(i in 1:(group-1)){
    for(j in (i+1):group){
      s = s + size^2/total *
        ( 2*mean(distance[ind[i,],ind[j,]]) -
            mean(distance[ind[j,],ind[j,]]) -
            mean(distance[ind[i,],ind[i,]]) )
    }
  }
  s/(group - 1)/(w/(total - group))
}
