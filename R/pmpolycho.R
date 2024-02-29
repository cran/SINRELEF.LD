pmpolycho <- function(thres1,thres2,espe1,espe2,sd1,sd2,poly,ncat){

  k <- transpose(c(ncat,ncat))

  threse1 <- transpose(c(-100,transpose(thres1), 5))
  threse2 <- transpose(c(-100,transpose(thres2), 5))

  Ct <- contiteor(threse1, threse2, k, poly)
  esppm <- especov(threse1, threse2, k, poly, Ct)
  pmcorr <- (esppm - (espe1*espe2)) / (sd1 * sd2)

  return(pmcorr)

}
