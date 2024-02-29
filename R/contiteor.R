contiteor<-function(thres1,thres2,ncat,r){

  TMP1 <- normalbi(thres1,thres2,ncat,r)
  TMP2 <- transpose(TMP1)

  ncat1 <- ncat[1]
  tope <- ncat1 + 1
  Ct <- matrix(0,ncat1,ncat1)

  for (i in 1:ncat1){
    for (j in 1:ncat1){
      Ct[i,j] <- TMP2[i,(tope - j)]
    }
  }
  Ct <- transpose(Ct)

}
