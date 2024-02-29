discreB<-function(X, THRES){

  #discretize a continuous variable
  # X     -> continuous data
  # THRES -> Thresholds

  siz <- size(X)
  n <- siz[1]

  siz <- size(THRES)
  k <- siz[1]
  m <- siz[2]

  Y <- matrix(0,n,m)
  for (i in 1:n){
    for (j in 1:m){
      Y[i,j] <- k+1
      for (h in 1:k){
        if (X[i,j] <= THRES[h,j]){
          Y[i,j]=h
          break
        }
      }
    }
  }

  return(Y)

}
