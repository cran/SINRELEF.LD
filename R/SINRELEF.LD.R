SINRELEF.LD<-function(L, PSI, THRES, ncat, model = 'linear', doublet_list, cor_doublet, N, CI = 90, display = TRUE){

  #### GLOBAL CHECKS ####

  if (missing(L)){
    stop("The argument L is not optional, please provide a valid Lambda vector.")
  }
  else {
    siz <- size(L)
    if (min(siz)>1){
      stop("The argument L is a matrix, not a vector. SINRELEF_LD is designed for assessing unifactorial solutions.")
    }
  }

  if (model == 'linear'){
    if (missing(PSI)){
      stop("The argument PSI is not optional for linear model. Please, provide a vector containing the residual variances.")
    }
    else {
      siz <- size(PSI)
      if (siz[1]>1){
        PSI <- transpose(PSI)
      }
    }
  }


  if (model == 'graded'){
    if (missing(THRES)){
      stop("The argument THRES is not optional for graded model. Please, provide a valid threshold vector.")
    }
    if (missing(ncat)){
      stop("The argument ncat is not optional for graded model. Please, provide the number of categories.")
    }
  }

  if (missing(doublet_list)){
    stop("The argument doublet_list is not optional, please provide a valid vector containing all the doublets.")
  }
  else {
    buff1 <- size(doublet_list)[1]
    buff2 <- size(doublet_list)[2]
    if (buff1 == 1 && buff2 == 2){ # ok, just 1 doublet
    }
    if (buff1 == 2 && buff2 == 1){ # ok, just 1 doublet
    }
    if (min(buff1,buff2)==2){
      # not a vector
      if (buff1< buff2){
        #2 vector structure
        doublet_list = transpose(doublet_list)
      }
      # 2 column structure
      doublet_list <- c(t(doublet_list))
    }


  }

  if (missing(cor_doublet)){
    stop("The argument doublet_list is not optional, please provide a valid vector containing the correlation of each doublet included in doublet_list.")
  }

  #if (min(size(doublet_list))/2 != min(size(cor_doublet))){
  #  stop("The doublet_list has to be consistent with the size of the cor_doublet (for each doublet, 1 correlation).")
  #}

  if (missing(N)){
    stop("The argument N is not optional, please provide the sample size used in the original analysis.")
  }

  #### ALL CHECKS OK, BEGIN ANALYSIS ####


  if (model == 'linear'){
    # The PSI diagonal elements are the square roots  of the residual variances
    psivec <- sqrt(PSI)
    PSI <- diag(as.vector(psivec))

    #Ree: Must be obtained from the "with" elements in the standardized model results. STDXY

    siz <- size(L)
    m <- siz[1]
    r <- siz[2]


    if(m < r){
      L <- transpose(L)
      siz <- size(L)
      m <- siz[1]
      r <- siz[2]
    }

    Ree <- diag(1, m, m)

    #Ree
    #number of doublets
    n_doublets <- max(size(cor_doublet))

    j <- 0
    for (i in 1:n_doublets){
      Ree[doublet_list[i+j],doublet_list[i+j+1]] = cor_doublet[i]
      Ree[doublet_list[i+j+1], doublet_list[i+j]] = cor_doublet[i]
      j <- j +1
    }

    OUT <- omegaldconti(L, Ree, PSI, display = FALSE)
    OUT_replic <- rreplicbase(L, Ree, PSI, N, 1000, CI, display = FALSE)

    if (display==TRUE){
      cat(sprintf('Results under the linear FA\n\n'))
      cat(sprintf('Correct omega estimate:\n'))
      cat(sprintf('omld with %2.0f%% C.I.\n',CI))

      cat(sprintf('%3.2f  (%3.2f; %3.2f)',OUT$omld, OUT_replic$low_ci_omld, OUT_replic$high_ci_omld))

      cat(sprintf('\n\nOmega estimate under local independence\n'))
      cat(sprintf('omli with %2.0f%% C.I.\n',CI))
      cat(sprintf('%3.2f  (%3.2f; %3.2f)',OUT$omli, OUT_replic$low_ci_omli, OUT_replic$high_ci_omli))

      cat(sprintf('\n\nScore relative efficiency with the %2.0f%% C.I.\n',CI))
      cat(sprintf('%3.2f  (%3.2f; %3.2f)',OUT$relef, OUT_replic$low_ci_relef, OUT_replic$high_ci_relef))

      cat('\n\n')

      #doublets
      cat('\nBivariate doublet results\n\n')
      cat(sprintf('Doublet     Relative efficiency    %2.0f%% C.I.\n',CI))
    }

    relef_doublet <- matrix(NA,n_doublets)

    j <- 0
    for (i in 1:n_doublets){
      relef <- doubletrelefconti(L[doublet_list[i+j]], L[doublet_list[i+j+1]],Ree[doublet_list[i+j],doublet_list[i+j+1]], display = FALSE)
      L_d <- transpose(c(L[doublet_list[i+j]], L[doublet_list[i+j+1]]))
      Ree_d <- matrix(1,2,2)
      Ree_d[1,2] <- Ree[doublet_list[i+j],doublet_list[i+j+1]]
      Ree_d[2,1] <- Ree[doublet_list[i+j],doublet_list[i+j+1]]
      PSI_d <- diag(c(PSI[doublet_list[i+j], doublet_list[i+j]],PSI[doublet_list[i+j+1], doublet_list[i+j+1]]))
      OUT_replic_d <- rreplicbase(L_d, Ree_d, PSI_d, N, 500, CI , display = FALSE)
      if (display==TRUE){
        cat(sprintf('%2.0f-%2.0f       %3.2f                   (%3.2f; %3.2f)\n',doublet_list[i+j], doublet_list[i+j+1],relef, OUT_replic_d$low_ci_relef, OUT_replic_d$high_ci_relef))
      }
      j <- j + 1
      relef_doublet[i] <- relef
    }


    #Individual item-deletion results

    out_del <- relefdeleted(L, Ree, PSI, display = FALSE)
    omegadel <- out_del$om
    infodel <- out_del$infch

    if (display==TRUE){
      cat('\n')
      cat('\nIndividual item-deletion results\n\n')
      cat('Item     Reliability estimate   Relative info change\n')

      for (i in 1:m){
        cat(sprintf('%3.0f      %5.4f                 % 5.4f',i, omegadel[i], infodel[i]))
        if (omegadel[i] < OUT_replic$low_ci_omld){
          cat('*')
        }
        cat('\n')
      }
      cat('\n')

      }
    }


  #############################################
  ###############    GRADED    ################
  #############################################

  if (model == 'graded'){

    siz <- size(L)
    m <- siz[1]
    r <- siz[2]


    if(m < r){
      L <- transpose(L)
      siz <- size(L)
      m <- siz[1]
      r <- siz[2]
    }

    # PSI: The residual variances are not provided in the output because the total variances are 1. To obtain them:

    psivec <- c()
    for (i in 1:m){
      psivec[i] <- sqrt(1-(L[i]^2))
    }
    PSI <- diag(psivec)

    #Ree
    #number of doublets
    n_doublets <- max(size(cor_doublet))
    Ree <- diag(1, m, m)

    j <- 0
    for (i in 1:n_doublets){
      Ree[doublet_list[i+j],doublet_list[i+j+1]] = cor_doublet[i]
      Ree[doublet_list[i+j+1], doublet_list[i+j]] = cor_doublet[i]
      j <- j +1
    }

    # Check THRES.
    #THRES: MPLUS provides them but with an awkward ordering. The user should provide them as vector, and we are going to re-arrange them

    buff1 <- size(THRES)[1]
    buff2 <- size(THRES)[2]

    if (buff1 == 1){
      THRES <- transpose(THRES)
    }

    if (min(buff1,buff2)==1){
      # vector, re-arrange
      THRES_matrix <- matrix(NA, nrow = ncat-1, ncol = m)
      k <- 1
      for (i in 1:m){
        for (j in 1:(ncat-1)){
          THRES_matrix[j,i] <- THRES[k]
          k <- k+1
        }
      }

      THRES <- THRES_matrix
    }
    else {
      # matrix, check if everything is ok
      if (buff2 < buff1){
        THRES <- transpose(THRES)
      }
      if (buff1 != (ncat-1)){
        stop("The thresholds (THRES) are not consistent with the number of categories (ncat).")
      }
      if (buff2 != m){
        stop("The thresholds (THRES) are not consistent with the number of items.")
      }
    }

    OUT <- omegalduva(THRES, L, Ree, PSI, ncat, display = FALSE)

    OUT_replic <- rreplicdiscre(THRES, L, Ree, PSI, N, 500, CI, display = FALSE)

    if (display==TRUE){
      cat(sprintf('Results under the nonlinear UVA FA\n\n'))
      cat(sprintf('Correct omega estimate:\n'))
      cat(sprintf('omld with %2.0f%% C.I.\n',CI))

      cat(sprintf('%3.2f  (%3.2f; %3.2f)',OUT$omld, OUT_replic$low_ci_omld, OUT_replic$high_ci_omld))

      cat(sprintf('\n\nOmega estimate under local independence\n'))
      cat(sprintf('omli with %2.0f%% C.I.\n',CI))
      cat(sprintf('%3.2f  (%3.2f; %3.2f)',OUT$omli, OUT_replic$low_ci_omli, OUT_replic$high_ci_omli))

      cat(sprintf('\n\nScore relative efficiency with the %2.0f%% C.I.\n',CI))
      cat(sprintf('%3.2f  (%3.2f; %3.2f)',OUT$relef, OUT_replic$low_ci_relef, OUT_replic$high_ci_relef))

      cat('\n\n')

      #doublets
      cat('\nBivariate doublet results\n\n')
      cat(sprintf('Doublet     Relative efficiency    %2.0f%% C.I.\n',CI))
    }

    relef_doublet <- matrix(NA,n_doublets)

    j <- 0
    for (i in 1:n_doublets){


      #Polychoric correlation between the two items in the model
      poly12 <- L[doublet_list[i+j]] * L[doublet_list[i+j+1]] + PSI[doublet_list[i+j], doublet_list[i+j]] * PSI[doublet_list[i+j+1], doublet_list[i+j+1]] * cor_doublet[i]

      #Polychoric correlation between the two items if they were independant
      poly12li <- L[doublet_list[i+j]] * L[doublet_list[i+j+1]]

      L_d <- transpose(c(L[doublet_list[i+j]], L[doublet_list[i+j+1]]))
      Ree_d <- matrix(1,2,2)
      Ree_d[1,2] <- Ree[doublet_list[i+j],doublet_list[i+j+1]]
      Ree_d[2,1] <- Ree[doublet_list[i+j],doublet_list[i+j+1]]
      PSI_d <- diag(c(PSI[doublet_list[i+j], doublet_list[i+j]],PSI[doublet_list[i+j+1], doublet_list[i+j+1]]))

      relef <- omegalduva(cbind(THRES[,doublet_list[i+j]], THRES[,doublet_list[i+j+1]]), transpose(c(L[doublet_list[i+j]], L[doublet_list[i+j+1]])), Ree_d, PSI_d, ncat, display=FALSE)$relef


      THRES_d <- cbind(transpose(THRES[,doublet_list[i+j]]),transpose(THRES[,doublet_list[i+j+1]]))

      OUT_replic_d <- rreplicdiscre(THRES_d, L_d, Ree_d, PSI_d, N, 500, CI, display = FALSE)

      if (display==TRUE){
        cat(sprintf('%2.0f-%2.0f       %3.2f                   (%3.2f; %3.2f)\n',doublet_list[i+j], doublet_list[i+j+1],relef, OUT_replic_d$low_ci_relef, OUT_replic_d$high_ci_relef))
      }
      j <- j + 1

      relef_doublet[i] <- relef
    }


    #Individual item-deletion results

    out_del <- relefdeleteduva(THRES, L, Ree, PSI, ncat, display = FALSE)
    omegadel <- out_del$om
    infodel <- out_del$infch

    if (display==TRUE){
      cat('\n')
      cat('\nIndividual item-deletion results\n\n')
      cat('Item     Reliability estimate   Relative info change\n')

      for (i in 1:m){
        cat(sprintf('%3.0f      %5.4f                 % 5.4f',i, omegadel[i], infodel[i]))
        if (omegadel[i] < OUT_replic$low_ci_omld){
          cat('*')
        }
        cat('\n')
      }
      cat('\n')
    }
  }

  OUT<-list('omld'=OUT$omld,'omli'=OUT$omli, 'relef'=OUT$relef, 'relef_doublet'=relef_doublet,'omega_del'=omegadel, 'r_info_del'=infodel )
  invisible(OUT)

}
