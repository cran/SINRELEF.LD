replicdiscre<-function(THRES, lam, Ree, PSI, N,display = TRUE){


  ## 1 replic discre UVA Omega

  siz <- size(lam)
  n <- siz[1]
  m <- siz[2]


  th <- transpose(rnorm(N))
  L <- chol(Ree)

  TMP <- matrix(rnorm(N*n), ncol=n)

  E <- TMP %*% L
  EI <- TMP

  ##Note: Z1 is the data matrix with local dependencies according to Ree
  #       Z2 is the data matrix if the item were locally independent

  Z1 <- transpose((lam %*% transpose(th)) + (PSI %*% transpose(E)))
  Z2 <- transpose((lam %*% transpose(th)) + (PSI %*% transpose(EI)))

  X1 <- discreB(Z1, THRES)
  X2 <- discreB(Z2, THRES)

  # Note: the omegas in the simulation are obtained as the squared correlations between
  # the sum scores and the true theta levels. This is possible here because the true theta levels are known.


  sumsco1 <- X1 %*% matrix(1,n,1)
  tmp <- cor(sumsco1,th)
  omrep <- tmp * tmp
  sumsco2 <- X2 %*% matrix(1,n,1)
  tmp2 <- cor(sumsco2,th)
  omrepli <- tmp2 * tmp2
  relefrep <- (omrep * (1-omrepli)) / (omrepli * (1-omrep))



  OUT<-list('omrep'=omrep,'omrepli'=omrepli, 'relefrep'=relefrep)

  invisible(OUT)


}
