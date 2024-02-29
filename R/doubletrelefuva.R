doubletrelefuva<-function(thres1, thres2, lam1, lam2, poly12, poly12li, ncat, display = TRUE){

  ## Relative efficiency of a doublet score base on graded model

  K <- c(ncat,ncat)

  # Step 1. Model-implied items means and standard deviations


  out <- espeit(transpose(thres1))
  esp1 <- out$esp
  var1 <- out$var

  out <- espeit(transpose(thres2))
  esp2 <- out$esp
  var2 <- out$var

  sd1 <- sqrt(var1)
  sd2 <- sqrt(var2)

  # Step 2. Model-implied items inter-item product-moment correlations

  pps1 <- ppolyser(lam1, transpose(thres1), sd1)
  pps2 <- ppolyser(lam2, transpose(thres2), sd2)

  # Step 3. Model-implied items inter-item product-moment correlations

  threse1 <- c(-100, transpose(thres1), 5)
  threse2 <- c(-100, transpose(thres2), 5)

  Ct <- contiteor(threse1, threse2, K, poly12)

  Ctli <- contiteor(threse1, threse2, K, poly12li)

  esppm <- especov(threse1, threse2, K, poly12, Ct)
  esppmli <- especov(threse1, threse2, K, poly12li, Ctli)

  r12 <- (esppm - (esp1*esp2)) / (sd1*sd2)
  r12li <- (esppmli - (esp1*esp2)) / (sd1*sd2)

  # Step 4. Doublet-omega coefficient and relative efficiency

  nume <- ((pps1 * sd1) + (pps2 * sd2))^2
  deno <- var1 + var2 + 2 * sd1 * sd2 * r12
  denoli <- var1 + var2 + 2 * sd1 * sd2 * r12li
  omega <- nume / deno
  omegali <- nume / denoli

  relef <- (omega * (1-omegali)) / (omegali * (1-omega))

  if (display == TRUE){

    'Results under the nonlinear UVA FA\n\n'

    cat(sprintf('Doublet correct omega estimate: %7.5f\n\n',omega))
    cat(sprintf('Doublet omega estimate under local independence: %7.5f\n\n',omega))
    cat(sprintf('Doublet score relative efficiency: %7.5f\n\n',omega))

  }

  OUT <- list('omega'=omega,'omegali'=omegali,'relef'=relef)
  invisible(OUT)
}
