rsq <- function(model = NULL, name = NULL) {
  if (is.null(model))
    stop("Please provide a model!")
  if (is.null(name))
    stop("Please provide a name of an endogenous variable!")
  df <- model@data@observed
  df <- df[name]
  S <- cov(df)
  res_mat <- model$S$values[name, name]
  R2 <-  1 - diag(res_mat) / diag(S)
  #names(R2) <- name
  return(R2)
}