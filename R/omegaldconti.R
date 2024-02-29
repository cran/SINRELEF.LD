omegaldconti<-function(L, Ree, PSI, display = TRUE){

  if (missing(L)){
    stop("The argument L is not optional, please provide Loading column pattern")
  }

  if (missing(Ree)){
    stop("The argument Ree is not optional, please provide Residual correlation matrix")
  }

  if (missing(PSI)){
    stop("The argument PSI is not optional, please prvide residual standard deviations (diagonal)")
  }

  if (display!=0 && display!=1){
    stop("display argument has to be logical (TRUE or FALSE, 0 or 1)")
  }

  ################################# Everything  OK #################################
  ################################# Begin Analysis #################################


  ## Omega reliability and relative efficiency with correlated residuals
  ## based on linear model

  siz <- size(L)
  m <- siz[1]
  r <- siz[2]


  if(m < r){
    L <- transpose(L)
    siz <- size(L)
    m <- siz[1]
    r <- siz[2]
  }

  psivec<- transpose(diag(PSI))
  u <- matrix(1,m,1)
  nume <- transpose(u) %*% L %*% transpose(L) %*% u

  denold <- transpose(psivec) %*% Ree %*% psivec
  denoli <- transpose(psivec) %*% psivec

  omld <- nume / (nume+denold)
  omli <- nume / (nume+denoli)

  infold <- nume / denold
  infoli <- nume /denoli

  reld <- infold / infoli

  if (display==TRUE){

    cat(sprintf('  Omega reliability estimate: %7.5f\n',omld))
    cat(sprintf('  Predicted omega estimate if items were locally independent: %7.5f\n',omli))
    cat(sprintf('  Relative efficiency of the locally dependent scores: %7.5f\n',reld))

  }
  OUT<-list('omld'=omld,'omli'=omli, 'relef'=reld)

  invisible(OUT)


}
