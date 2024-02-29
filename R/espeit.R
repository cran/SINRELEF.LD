espeit <- function(thres){

  k <- size(thres)[1]
  cat <- k+1
  p <- matrix(0,cat,1)
  p[1] <- fnor3(thres[1])
  sum <- p[1]
  sum2 <- 0
  sum3 <- 0

  for (i in 2:k){
    tmp <- 0
    tmp <- fnor3(thres[i]) - fnor3(thres[i-1])
    sum <- sum + tmp
    p[i] <- tmp
  }
  p[cat] <- 1 - sum

  for (j in 1:cat){
    sum2 <- sum2 + (j * p[j])
    sum3 <- sum3 + (j * j * p[j])
  }

  esp <- sum2
  var <- sum3 - (esp * esp)

  out <- list('esp'=esp,'var'=var)
  return(out)
}
