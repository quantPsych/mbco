# Name modification: I added within_lav; 4/27/2018
# This sim function uses within design with focus on lavaan except for prof-lik and mbco
# I added lav_bca.R at the end of this file 4/6/2018
#I modified this function for simulation study single mediator and two-mediators chains for mbco test in psych methods paper: 5/26/2017
# I added BCa 3/30/16
# Simulation function for 3 path indirect effects
# This function accepts sample size n and a list of named parameter vectors
# generateData <- function(n,param){
#   #This function generates data. The model is a micromediational chain with three mediators.
#   #The vector parameters are bs & cp
#   #note: param[] is a data frame and therefore, for some strange reason, if you multiply the vector by the data frame,
#   # it does return single number of type data frame

popModel <- function(param) {
  quant <- param[['quant']][[1]] #get the formula
  varName <- all.vars(quant) #variables names in quant
  ## Single Mediator Chain
  l1 <- paste0("m1~", param[['b1']], "*x")
  l2 <- ""
  l3 <- paste0("y~", param[['b2']], "*m1")
  l4 <- ""
  l5 <- paste0("x~~", param[["s2x"]], "*x")
  l6 <- paste0("m1~~", param[["s2em1"]], "*m1")
  l7 <- ""
  l8 <- ""
  l9 <- paste0("y~~", param[["s2ey"]], "*y")
  ## Tow Mediator Chain
  if (length(varName) == 3) {
    l1 <- paste0("m1~", param[['b1']], "*x")
    l2 <- paste0("m2~", param[['b2']], "*m1")
    l3 <- paste0("y~", param[['b3']], "*m2")
    l4 <- "\n"
    l5 <- paste0("x~~", param[["s2x"]], "*x")
    l6 <- paste0("m1~~", param[["s2em1"]], "*m1")
    l7 <- paste0("m2~~", param[["s2em2"]], "*m2")
    l8 <- "\n"
    l9 <- paste0("y~~", param[["s2ey"]], "*y")
  }
  popmodel <- paste(l1, l2, l3, l4, l5, l6, l7, l8, l9, sep = "\n")
  return(popmodel)
}

dataGen <- function(param, popModel) {
  n <- param[['samp']]
  df <- lavaan::simulateData(popModel, sample.nobs = n)
  return(df)
}

#This function accepts a result of lavaan fit
#Output: a list of two things. 1. a vector of coef estimates 2. covariance matrix of the coef estimates
getCoefs <- function(fit, quant) {
  pEst <- coef(fit) #parameter estimate
  name1 <- all.vars(quant)
  coefRes <- pEst[name1]
  ##Cov
  cov1 <- (vcov(fit)) #covariance of the coef estimates
  covRes <- cov1[name1, name1]
  return(list(coefRes, covRes))
}
#This function calculates CIs based on bootstrap

### This for a new simulation study where we use two seperarte dfs
bootLavanCI3 <- function(fit, boot.ci.type = "perc") {
  res <- parameterEstimates(fit, boot.ci.type = boot.ci.type)
  res1 <- res[res$op == ':=', c("ci.lower", "ci.upper")]
  return(unlist(res1))
}

#This function accepts x as a list of two things: 1. a vector of the coef estimates 2. the covariance matrix of the coef estimates
# argument quant=indirect eff & b &ect. argument type accepts two values corresponding to the methods of calculating CI (MC)
#Output: a vector of confidence interval limits
CI <- function(x, quant, type = "mc") {
  bEst <- x[[1]]
  bCov <- lav_matrix_vech(x[[2]])
  if (type == "mc") {
    res1 <- ci(bEst, bCov, quant = quant, type = type)
    res <- list(res1[[1]], res1[[2]])
  }
  else{
    mu.x <- bEst[1]
    mu.y <- bEst[2]
    se.x <- sqrt(bCov[1])
    se.y <- sqrt(bCov[3])
    res2 <- medci(mu.x, mu.y, se.x, se.y)
    res <- res2[[1]]
  }

  res
}

#This function accepts a vector of confidence interval limits
#Output: returns 1 if it is a false positive (incorrect rejection) or 0 otherwise
falsePositive <- function(x) {
  ifelse(x[1] * x[2] > 0, 1, 0)
}


#################### Main Sim Function

sim <- function(param,
                R = 1000L,
                type = c("mc"),
                df_prof = 1) {
  res <- numeric(0) #final results vector
  popmodel <-
    popModel(param) # generating population model lavaan syntax
  #model1 <- as.character( param[['model']] )
  quant <- param[['quant']][[1]]
  varName <- all.vars(quant) #variables names in quant
  quant_ch <- sub("~", replacement = "", quant)[2]
  ind_cons <- paste0("ind:=", quant_ch) #constraints for Asymptotic

  if (length(varName) == 3) {
    model1 <- "m1 ~ b1*x \n m2 ~ b2*m1 \n y ~ b3*m2"
  }
  else{
    model1 <- "m1 ~ b1*x \n y ~ b2*m1"
  }

  df <-  dataGen(param, popmodel)
  fit <-
    sem(model1,
        data = df,
        constraints = ind_cons) #estimate the sem using lavaan package

  fit_boot <- sem(model1,
                  data = df,
                  se = "boot",
                  bootstrap = R,
                  constraints = ind_cons) #estimate the sem using lavaan package

  ### MC ci
  if ("mc" %in% type) {
    CoefCovList <-
      getCoefs(fit, quant) # list of two things: 1. a vector of coef estimates 2. covariance matrix of the coef estimates
    ciRes <-
      CI(CoefCovList, quant = quant, type = 'mc') #produces a vector of CI
    ciResMC <- ciRes[[1]]
    falsePosMC <- falsePositive(ciResMC) #scalar 0 and 1
    names(ciResMC) <- c("MCL", "MCU")
    names(falsePosMC) <- "MC_Rej"
    res <- c(res, ciResMC, falsePosMC)
  }

  #### DOP, only for a single mediator model
  if ("dop" %in% type) {
    CoefCovList <-
      getCoefs(fit, quant) # list of two things: 1. a vector of coef estimates 2. covariance matrix of the coef estimates
    ciResDOP <-
      CI(CoefCovList, quant = quant, type = 'dop') #produces a vector of CI
    falsePosDOP <- falsePositive(ciResDOP) #scalar 0 and 1
    names(ciResDOP) <- c("dopL", "dopU")
    names(falsePosDOP) <- "dop_Rej"
    res <- c(res, ciResDOP, falsePosDOP)
  }

  # Asymp
  if ("asymp" %in% type) {
    resAsymp <- subset(parameterEstimates(fit), op == ":=")
    ciResAsymp <- resAsymp[c("ci.lower", "ci.upper")]
    falsePosAsymp <- falsePositive(ciResAsymp) #scala0 and 1
    names(ciResAsymp) <- c("AsympL", "AsympU")
    names(falsePosAsymp) <- "Asymp_Rej"
    res <- c(res, ciResAsymp, falsePosAsymp)

  }

  #Bootstrap Percentile
  if ("perc" %in% type) {
    ciResPerc <- bootLavanCI3(fit_boot, boot.ci.type = "perc")
    falsePosPerc <- falsePositive(ciResPerc)
    names(ciResPerc) <- c("PercL", "PercU") # Percentile CI
    names(falsePosPerc)  <- "perc_Rej"
    res <- c(res, ciResPerc, falsePosPerc)
  }

  t1_1m <- t1_2m <- FALSE
  #Bootstrap BC
  if(length(varName) == 2) {
    b1 <- param[['b1']]
    b2 <- param[['b2']]
    t1_1m <- ifelse(b1 * b2 == 0, TRUE, FALSE)
  }

  if(length(varName) == 3) {
    b1 <- param[['b1']]
    b2 <- param[['b2']]
    b2 <- param[['b3']]
    t1_2m <- ifelse(b1 * b2 * b3 == 0, TRUE, FALSE)
  }

  if ("bc" %in% type)
    if (t1_1m || t1_2m) {
      ciResBC <-  bootLavanCI3(fit_boot, boot.ci.type = "bca.simple")
      falsePosBC <- falsePositive(ciResBC)
      names(ciResBC)  <- c("bcL", "bcU")  # BC CI
      names(falsePosBC) <-  "bc_Rej"
      res <- c(res, ciResBC, falsePosBC)
    } else {
      ciResBC <- c(NA, NA)
      falsePosBC <- NA
      names(ciResBC)  <- c("bcL", "bcU")  # BC CI
      names(falsePosBC) <-  "bc_Rej"
      res <- c(res, ciResBC, falsePosBC)
    }

  ## Joing Sig Test
  if ("joint" %in% type) {
    b_est <-
      subset(parameterEstimates(fit), label %in% all.vars(quant))
    joint_rej <- 0
    if (all(b_est[['pvalue']] < .05))
      joint_rej <- 1
    res <- c(res, joint_rej = joint_rej)
  }

  #### BCA Bootstrap
  if ("bca" %in% type)
    if (ifelse(length(varName) == 3, b1 * b2 * b3 == 0, FALSE)  |
        ifelse(length(varName) == 2, b1 * b2 == 0, FALSE)) {
      ciResBCa <- lav_bca(df,
                          model = model1,
                          R = R,
                          quant = quant)
      falsePosBCa <- falsePositive(ciResBCa)
      names(ciResBCa)  <- c("bcaL", "bcaU")  # BC CI
      names(falsePosBCa) <-  "bca_Rej"
      res <- c(res, ciResBCa, falsePosBCa)
    } else {
      ciResBCa <- c(NA, NA)
      falsePosBCa <- NA
      names(ciResBCa)  <- c("bcaL", "bcaU")  # BC CI
      names(falsePosBCa) <-  "bca_Rej"
      res <- c(res, ciResBCa, falsePosBCa)
    }

  #### MBCO
  if ("mbco" %in% type) {
    if (length(varName) == 3) {
      mbco_res <- mbcoSimFun2Med(param, sim_df = df)
      res <- c(res, mbco_res)
    }
    else{
      mbco_res <- mbcoSimFun1Med(param, sim_df = df)
      res <- c(res, mbco_res)
    }
  }

  # Profile Likelihood
  if ("prof" %in% type) {
    if (length(varName) == 3) {
      ciProf <- profLikSim2Med(param, sim_df = df)
      falsePosProf <- falsePositive(ciProf)
      names(falsePosProf) <- "prof_Rej"
      res <- c(res, ciProf, falsePosProf)
    }
    else{
      ciProf <- profLikSim1Med(param, sim_df = df)
      falsePosProf <- falsePositive(ciProf)
      names(falsePosProf) <- "prof_Rej"
      res <- c(res, ciProf, falsePosProf)
    }
  }

  param[['quant']] <- as.character(param[['quant']])
  res <- c(param, res)
  return(res)
}
