ppolyser <- function(r,thx,sigx){
  n <- size(thx)[2]
  sum <- 0

  for (i in 1:n){
    sum <- sum + ordnor(thx[i])
  }
  pr <- (r / sigx) * sum
  return(pr)
}
