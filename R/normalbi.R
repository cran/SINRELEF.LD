normalbi<-function(tha,thb,Ki,r){

  Ka <- Ki[1]
  Kb <- Ki[2]

  C <- matrix(0,Kb,Ka)

  for (j in 1:Kb){
    for (i in 1:Ka){
      ai <- -1*tha[i]
      as <- -1*tha[i+1]
      bi <- -1*thb[j]
      bs <- -1*thb[j+1]
      Ass <- bivnor(as,bs,r)
      Bsi <- 0
      Cis <- 0
      Dii <- 0

      if ((ai != -100) & (bi != -100)){
        Bsi <- bivnor(as,bi,r)
        Cis <- bivnor(ai,bs,r)
        Dii <- bivnor(ai,bi,r)
      }
      if ((ai != -100) & (bi != -100)){
        Cis <- bivnor(ai,bs,r)
        Dii <- bivnor(ai,bi,r)
      }
      if ((ai == -100) & (bi != -100)){
        Bsi <- bivnor(as,bi,r)
        Dii <- bivnor(ai,bi,r)
      }
      T <- Ass - Bsi - Cis + Dii
      C[Kb+1-j,i] <- T
    }
  }

  return(C)

}
