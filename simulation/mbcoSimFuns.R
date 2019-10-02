## This document contains mbcoFun1Med and mbcoFun2Med
## Last Modified 4/6/2018


#############################    mbcoSimFun1Med
##

mbcoSimFun1Med <-
  function(param, optim = "NPSOL", sim_df) {
    ## Population Model to Generate Data
    b1 <- param$b1 #param["b1"]
    b2 <- param$b2 #param["b2"]
    #cp <- param$cp #param["cp"]
    n <- param$samp # param["n"]
    s2x <- param$s2x #param["s2x"]
    s2m <- param$s2em1  #param["s2m"]
    s2y <- param$s2ey  #param["s2y"]

    manifests <- c("x", "m1", "y")
    mxOption(NULL, "Default optimizer", optim)  ## I changed optimizer per suggestion
    mxOption(NULL, "Standard Errors" , "No")  ## I changed optimizer per suggestion

    popModel <- mxModel(
      name = "populationModel",
      type = "RAM",
      manifestVars = manifests,
      mxPath(
        from = 'one',
        to = manifests,
        free = TRUE,
        values = 0,
        labels = c("mu_x", "int_m", "int_y")
      ),
      mxPath(
        from = "x",
        to = "m1",
        arrows = 1,
        free = TRUE,
        values = b1,
        labels = "b1"
      ),
      mxPath(
        from = "m1",
        to = "y",
        arrows = 1,
        free = TRUE,
        values = b2,
        labels = "b2"
      ),
      mxPath(
        from = manifests,
        arrows = 2,
        free = TRUE,
        values = c(s2x, s2em1, s2y),
        labels = c("s2x", "s2em1", "s2ey")
      )
    )


    medModel <-
      mxModel(model = popModel,
              name = "Unconstrained Model",
              mxData(sim_df, type = "raw"))

    ## Constrained Model
    medModelconst <- mxModel(
      model = popModel,
      name = "constrained Model",
      mxData(sim_df, type = "raw"),
      mxAlgebra(b1 * b2, name = "ind"),
      mxConstraint(ind == 0, name = "ind0")
    )

    ## Runing Constrained Model
    medModelFit <- mxRun(medModel)
    sumMedModel <- summary(medModelFit)


    medModelconstFit <- mxRun(medModelconst)

    statSum <-
      summary(
        medModelconstFit,
        SaturatedLikelihood = medModelFit,
        SaturatedDoF = sumMedModel$degreesOfFreedom
      )

    mbcoTest <-
      c(
        mbco_chisq = statSum$Chi,
        mbco_pvalue = statSum$p,
        mbco_df = statSum$ChiDoF,
        mbco_reject = ifelse(statSum$p < .05, 1, 0)
      )
    return(mbcoTest)

  }

#############################    mbcoSimFun2Med
##

mbcoSimFun2Med <-
  function(param, optim = "NPSOL", sim_df) {
    ## Population Model to Generate Data
    b1 <- param$b1
    b2 <- param$b2
    b3 <- param$b3
    n <- param$samp
    s2x <- param$s2x
    s2m1 <- param$s2em1
    s2m2 <- param$s2em2
    s2y <- param$s2ey

    manifests <- c("x", "m1", "m2", "y")
    mxOption(NULL, "Default optimizer", optim)  ## I changed optimizer per suggestion
    mxOption(NULL, "Standard Errors" , "No")  ## I changed optimizer per suggestion


    popModel <- mxModel(
      name = "populationModel",
      type = "RAM",
      manifestVars = manifests,
      mxPath(
        from = 'one',
        to = manifests,
        free = FALSE,
        values = 0,
        labels = c("mu_x", "int_m1", "int_m2", "int_y")
      ),
      mxPath(
        from = "x",
        to = "m1",
        arrows = 1,
        free = TRUE,
        values = b1,
        labels = "b1"
      ),
      mxPath(
        from = "m1",
        to = "m2",
        arrows = 1,
        free = TRUE,
        values = b2,
        labels = "b2"
      ),
      mxPath(
        from = "m2",
        to = "y",
        arrows = 1,
        free = TRUE,
        values = b3,
        labels = "b3"
      ),
      mxPath(
        from = manifests,
        arrows = 2,
        free = TRUE,
#        values = c(s2x, s2em1, s2em2, s2ey),
        labels = c("s2x", "s2em1", "s2em2", "s2ey")
      )
    )

    medModel <-
      mxModel(model = popModel,
              name = "Unconstrained Model",
              mxData(sim_df, type = "raw"))

    ## Constrained Model
    medModelconst <- mxModel(
      model = popModel,
      name = "constrained Model",
      mxData(sim_df, type = "raw"),
      mxAlgebra(b1 * b2 * b3, name = "ind"),
      mxConstraint(ind == 0, name = "ind0")
    )
    medModelFit <- mxRun(medModel)

    medModelconstFit <- mxRun(medModelconst)

    sumMedModel <- summary(medModelFit)
    statSum <-
      summary(
        medModelconstFit,
        SaturatedLikelihood = medModelFit,
        SaturatedDoF = sumMedModel$degreesOfFreedom
      )
    mbcoTest <-
      c(
        mbco_chisq = statSum$Chi,
        mbco_pvalue = statSum$p,
        mbco_df = statSum$ChiDoF,
        mbco_reject = ifelse(statSum$p < .05, 1, 0)
      )
    return(c(mbcoTest))
  }
