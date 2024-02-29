omegalduva<-function(THRES, L, Ree, PSI, ncat, display = TRUE){


  ## Omega reliability and relative efficiency with correlated residuals
  ## based on non-linear UVA model

  siz <- size(L)
  m <- siz[1]
  r <- siz[2]


  if(m < r){
    L <- transpose(L)
    siz <- size(L)
    m <- siz[1]
    r <- siz[2]
  }

  espe <- matrix(0,m,1)
  sds <- matrix(0,m,1)
  I <- diag(m)
  Robsld <- matrix(0,m,m)
  Robsli <- matrix(0,m,m)

  #Step 1 Model-implied inter-item correlation matrices
  Rld <- L %*% transpose(L) + PSI %*% Ree %*% PSI
  Rli <- L %*% transpose(L) + PSI %*% I %*% PSI

  #Step 2. Model-implied items means and standard deviations
  for (i in 1:m){
    out <- espeit(transpose(THRES[,i]))
    espe[i] <- out$esp
    sds[i] <- sqrt(out$var)
  }


  #Step 3. Model-implied total test variance (denominator term in omega)
  SDD <- matrix(0,m,m)
  diag(SDD)=sds

  for (i in 1:m){
    for (j in 1:m){
      if (i==j){
        Robsld[i,j] <- 1
        Robsli[i,j] <- 1
      }
      else {
        Robsld[i,j] <- pmpolycho(THRES[,i], THRES[,j], espe[i], espe[j], sds[i], sds[j], Rld[i,j], ncat)
        Robsli[i,j] <- pmpolycho(THRES[,i], THRES[,j], espe[i], espe[j], sds[i], sds[j], Rli[i,j], ncat)
      }
    }
  }

  Covobsld <- SDD %*% Robsld %*% SDD
  Covobsli <- SDD %*% Robsli %*% SDD
  vartotld <- sum(sum(Covobsld))
  vartotli <- sum(sum(Covobsli))

  #Step 4. Numerator term in Omega. This term is the same for omega-ld and for omega-li

  sumanume <- 0

  for (i in 1:m){
    tmp1 <- ppolyser(L[i], THRES[,i], sds[i])
    tmp2 <- tmp1 * sds[i]
    sumanume <- sumanume + tmp2
  }
  nume <- sumanume^2


  # nonlinear-omega estimates and relative eficiency

  omld <- nume / vartotld
  omli <- nume/vartotli
  relef <- (omld * (1-omli)) / (omli*(1-omld))


  if (display == TRUE){
    cat(sprintf('  Omega reliability estimate: %7.5f\n',omld))
    cat(sprintf('  Predicted omega estimate if items were locally independent: %7.5f\n',omli))
    cat(sprintf('  Relative efficiency of the locally dependent scores: %7.5f\n',relef))
  }

  OUT<-list('omld'=omld,'omli'=omli, 'relef'=relef)

  invisible(OUT)


}
