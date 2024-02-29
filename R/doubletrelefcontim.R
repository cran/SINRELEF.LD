doubletrelefcontim<-function(LAM, Ree, display = TRUE){

  ## Relative efficiency of a doublet score base on linear model - multiple

  L <- LAM

  siz <- size(L)
  m <- siz[1]
  r <- siz[2]


  if(m < r){
    L <- transpose(L)
    siz <- size(L)
    m <- siz[1]
    r <- siz[2]
  }

  relefm <- matrix(NA,m,m)

  for (i in 1:m){
    for (j in 1:m){
      relefm[i,j] <- doubletrelefconti(L[i], L[j], Ree[i,j])
    }
  }

  diag(relefm) <- NA

  invisible(relefm)
}
