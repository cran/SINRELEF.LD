relefdeleted<-function(LAM, Ree, PSI, display = TRUE){

  if (missing(LAM)){
    stop("The argument LAM is not optional, please provide Loading column pattern")
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

  psivec<- transpose(diag(PSI))
  u <- matrix(1,m,1)
  v <- matrix(1,(m-1),1)
  nume <- transpose(u) %*% L %*% transpose(L) %*% u
  denold <- transpose(psivec) %*% Ree %*% psivec
  omld <- nume / (nume+denold)

  if (display==TRUE){
    cat(sprintf('  Omega reliability estimate for the entire set: %7.5f\n\n',omld))
  }

  # Global CI for Relative information
  #IC <- .90
  #RELEF_global <- rreplicbase(L, Ree, PSI, 1000, 1000, IC, display = FALSE)$RELEF
  #RELEF_ic_low <- quantile(RELEF_global,probs = 1-IC)
  #RELEF_ic_high <- quantile(RELEF_global, probs = IC)

  infold <- nume / denold

  om <- matrix(NA,m,1)
  inf <- matrix(NA,m,1)
  infch <- matrix(NA,m,1)

  #relefic <- matrix(NA,m,2)

  items <- 1:m

  for (i in 1:m){

    tmp1 <- transpose(L[items[items!=i]])
    ptmp1 <- transpose(psivec[items[items!=i]])
    Red1 <- transpose(Ree[items[items!=i], items[items!=i]])

    om[i] <- (transpose(v) %*% tmp1 %*% transpose(tmp1) %*% v) / ((transpose(v) %*% tmp1 %*% transpose(tmp1) %*% v) + (transpose(ptmp1) %*% Red1 %*% ptmp1))
    inf[i] <- (transpose(v) %*% tmp1 %*% transpose(tmp1) %*% v) / (transpose(ptmp1) %*% Red1 %*% ptmp1)
    infch[i] <- (inf[i] - infold) / infold

    if (display==TRUE){
      cat(sprintf('  Omega reliability estimate if item %2.0f is deleted: %7.5f\n',i,om[i]))
      cat(sprintf('  Relative information change if item %2.0f is deleted: %7.5f\n\n',i,infch[i]))
    }
  }


  OUT<-list('om'=om,'infch'=infch)

  invisible(OUT)


}
