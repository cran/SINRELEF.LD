especov <- function(thres1,thres2,ncat,r,C){

  sum <- 0
  ncat1 <- ncat[1]

  for (i in 1:ncat1){
    for (j in 1:ncat1){
      tmp <- i * j * C[i,j]
      sum <- sum + tmp
    }
  }

  espcov <- sum
  return(espcov)
}
