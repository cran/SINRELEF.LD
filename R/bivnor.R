bivnor <- function(ah,ak,r){

  idig <- 15
  b <- 0

  gh <- gauss(-ah) / 2
  gk <- gauss(-ak) / 2

  if (r == 0){
    b <- 4 * gh * gk
    b <- max(b,0)
    b <- min(b,1)
    value <- b
    return(value)
  }

  rr <- (1 + r) * (1 - r)

  if (rr < 0 ){
    #ERROR 1 < |r|
  }

  if (rr == 0){

    if (r < 0){
      if (ah + ak < 0){
        b <- 2 * (gh + gk) - 1
      }
    }
    else {
      if (ah - ak < 0){
        b <- 2 * gk
      }
      else {
        b <- 2 *gh
      }
    }

    b <- max(b,0)
    b <- min(b,1)
    value <- b
    return(value)
  }

  sqr <- sqrt(rr)

  if (idig == 15){
    con <- 2.0 * pi * 0.0000000000000010 / 2
  }
  else {
    con <- pi
    for (i in 1:idig){
      con <- con / 10
    }
  }

  if (ah == 0){
    if (ak ==0){
      b <- 0.25 + 0.5 * asin(r) / pi
      b <- max(b,0)
      b <- min(b,1)
      value <- b
      return(value)
    }
  }

  if (ah ==0){
    if (ak != 0){
      b <- gk
      wh <- -ak
      wk <- (ah / ak - r) / sqr
      gw <- 2 * gk
      is <- 1
    }
  }

  if (ah != 0){
    if (ak == 0){
      b <- gh
      wh <- -ah
      wk <- (ak / ah - r) /sqr
      gw <- 2 * gh
      is <- -1
    }
  }

  if (ah != 0){
    if (ak !=0){
      b <- gh + gk
      if ((ah * ak) < 0){
        b <- b - 0.5
      }
      wh <- -ah
      wk <- (ak / ah -  r) / sqr
      gw <- 2 * gh
      is <- -1
    }
  }

  while (1){

    sgn <- -1
    t <- 0

    if (wk != 0){
      if (abs(wk) ==1){
        t <- wk * gw * (1-gw) / 2
        b <- b + sgn * t
      }
      else {
        if (1 < abs(wk)){
          sgn <- -sgn
          wh <- wh * wk
          g2 <- gauss(wh)
          wk <- 1 / wk

          if (wk < 0){
            b <- b + 0.5
          }

          b <- b - (gw + g2) / 2 + gw *g2
        }

        h2 <- wh * wh
        a2 <- wk * wk
        h4 <- h2 / 2
        ex <- exp(-h4)
        w2 <- h4 * ex
        ap <- 1
        s2 <- ap - ex
        sp <- ap
        s1 <- 0
        sn <- s1
        conex <- abs(con/wk)

        while (1){
          cn <- ap * s2 / (sn + sp)
          s1 <- s1 + cn
          if (abs(cn) <= conex){
            break
          }

          sn <- sp
          sp <- sp + 1
          s2 <- s2 - w2
          w2 <- w2 * h4 / sp
          ap <- -ap * a2

        }

        t <- 0.5 * (atan(wk) - wk * s1) /pi
        b <- b + sgn * t

      }
    }

    if ( 0 <= is){
      break
    }

    if (ak == 0){
      break
    }

    wh <- -ak
    wk <- (ah / ak - r) / sqr
    gw <- 2 * gk
    is <- 1

  }

  b <- max(b, 0)
  b <- min(b, 1)
  value <- b

  return(value)




}
