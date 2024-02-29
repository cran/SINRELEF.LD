fnor3 <- function(x){
  if (x > 0){
    d <- erf(x/sqrt(2)) / 2 + 0.5
  }
  else {
    d <- 0.5 - erf(-x/sqrt(2)) / 2
  }

  if (d > 0.999){ d <- 0.999}
  if (d < 0.0001) {d <- 0.0001}

  return(d)
}
