rreplicbase<-function(LAM, Ree, PSI, N, nreplic, IC, display = TRUE){

  ## Multiple replic omega continuous

  IC <- IC/100

  OMLD <- matrix(0,nreplic,1)
  OMLI <- matrix(0,nreplic,1)
  RELEF <- matrix(0,nreplic,1)

  for (i in 1:nreplic){
    out <- replicbase(LAM, Ree, PSI, N)

    OMLD[i] <- out$omrep
    OMLI[i] <- out$omrepli
    RELEF[i] <- out$relefrep
  }

  low_ci_omld <- quantile(OMLD,probs = 1-IC)
  high_ci_omld <- quantile(OMLD,probs = IC)

  low_ci_omli <- quantile(OMLI, probs = 1- IC)
  high_ci_omli <- quantile(OMLI, probs =IC)

  low_ci_relef <- quantile(RELEF,probs = 1-IC)
  high_ci_relef <- quantile(RELEF,probs = IC)

  if (display==TRUE){
    cat(sprintf('  %2.0f%% confidence interval omega reliability estimate: %7.5f %7.5f\n',IC*100,quantile(OMLD,probs = 1-IC),quantile(OMLD,probs = IC)))
    cat(sprintf('  %2.0f%% confidence interval omega estimate if items were locally independent: %7.5f %7.5f\n',IC*100,quantile(OMLI,probs = 1-IC),quantile(OMLI,probs = IC)))
    cat(sprintf('  %2.0f%% confidence interval relative efficiency of the locally dependent scores: %7.5f %7.5f\n',IC*100,quantile(RELEF,probs = 1-IC),quantile(RELEF,probs = IC)))
  }

  OUT<-list('OMLD'=OMLD,'OMLI'=OMLI, 'RELEF'=RELEF,'low_ci_omld'=low_ci_omld,'high_ci_omld'= high_ci_omld,'low_ci_omli'=low_ci_omli,'high_ci_omli'=high_ci_omli,'low_ci_relef'=low_ci_relef,'high_ci_relef'=high_ci_relef)

  invisible(OUT)
}
