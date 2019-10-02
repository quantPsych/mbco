## This document contains profLikSim1Med and profLikSim2Med
## Last Modified 4/6/2018


#############################    profLikSim1Med
##

profLikSim1Med <- function(param, optim = "NPSOL", sim_df) {
  # I added sim_df to accept data from sim function 4/6/2018

  manifests = c("x", "m1", "y")
  b1 <- param$b1 #param["b1"]
  b2 <- param$b2 #param["b2"]
  n <- param$samp # param["n"]
  s2x <- param$s2x #param["s2x"]
  s2m <- param$s2em1  #param["s2m"]
  s2y <- param$s2ey  #param["s2y"]

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
      values = c(s2x, s2m, s2y),
      labels = c("s2x", "s2em1", "s2ey")
    )
  )
  popModel <- mxOption(popModel, "Standard Errors"  , "No")

  medModel <-
    mxModel(
      model = popModel,
      name = "Single Mediator Model",
      mxData(sim_df, type = "raw"),
      mxAlgebra(b1 * b2, name = "ind"),
      mxCI("ind")
    )

  medModelFit <- mxRun(medModel, intervals = TRUE)

  stat <- summary(medModelFit)
  res <- unlist(stat$CI, use.names = TRUE)
  res <- res[c(1, 3)]
  names(res) <- c("profL", "profU")

  return(res)
}

#############################    profLikSim2Med
##

profLikSim2Med <- function(param, optim = "NPSOL", sim_df) {
  # I added sim_df to accept data from sim function 4/6/2018

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
      values = c(s2x, s2m1, s2m2, s2y),
      labels = c("s2x", "s2em1", "s2em2", "s2y")
    )
  )


  medModel <-
    mxModel(
      model = popModel,
      name = "Two Mediator Model",
      mxData(sim_df, type = "raw"),
      mxAlgebra(b1 * b2 * b3, name = "ind"),
      mxCI("ind")
    )

  medModelFit <- mxRun(medModel, intervals = TRUE)

  stat <- summary(medModelFit)
  res <- unlist(stat$CI, use.names = TRUE)
  res <- res[c(1, 3)]
  names(res) <- c("profL", "profU")

  return(res)
}
