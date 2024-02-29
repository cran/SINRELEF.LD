relefdeleteduva<-function(THRES, LAM, Ree, PSI, ncat, display = TRUE){

  ## Omega reliability and relative efficiency with correlated residuals
  ## based on linear model

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

  OUT <- omegalduva(THRES,LAM,Ree,PSI,ncat, display = FALSE)

  omld <- OUT$omld

  infold <- omld/(1-omld)

  if (display==TRUE){
    cat(sprintf('  Omega reliability estimate for the entire set: %7.5f\n\n',omld))
    cat(sprintf('  amount of information for the entire item set: %7.5f\n\n',infold))
  }

  om <- matrix(NA,m,1)
  inf <- matrix(NA,m,1)
  infch <- matrix(NA,m,1)

  items <- 1:m

  for (i in 1:m){

    THRESRED <- THRES[,items[items!=i]]
    siz <- size(THRES)
    m <- siz[1]
    r <- siz[2]
    if (r<m) {THRESRED <- transpose(THRESRED)}

    LAMRED <- transpose(L[items[items!=i]])
    ReeRED <- transpose(Ree[items[items!=i], items[items!=i]])
    PSIRED <- transpose(PSI[items[items!=i], items[items!=i]])

    om[i] <- omegalduva(THRESRED,LAMRED,ReeRED,PSIRED,ncat, display = FALSE)$omld
    inf[i] <- om[i] / (1-om[i])
    infch[i] <- (inf[i] - infold) / infold

    if (display==TRUE){
      cat(sprintf('  Omega reliability estimate if item %2.0f is deleted: %7.5f\n',i,om[i]))
      cat(sprintf('  Relative information change if item %2.0f is deleted: %7.5f\n\n',i,infch[i]))
    }
  }


  OUT<-list('om'=om,'infch'=infch)

  invisible(OUT)


}
