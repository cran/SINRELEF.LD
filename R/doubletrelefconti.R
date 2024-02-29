doubletrelefconti<-function(psi1, psi2, ree, display = TRUE){

  ## Relative efficiency of a doublet score base on linear model

  psi <- transpose(c(psi1,psi2))

  Ree <- diag(2)
  Ree[1,2] <- ree
  Ree[2,1] <- ree
  deno <- transpose(psi) %*%Ree %*% psi
  nume <- transpose(psi) %*% psi
  relef <- nume / deno

  if (display == TRUE){
    cat(sprintf('  Doublet score relative efficiency: %7.5f\n',relef))
  }

  invisible(relef)
}
