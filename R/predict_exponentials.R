predict_exponentials = function(x, coeff, description){
  transform = purrr::map(description, ~ x^.x)
  X = do.call(cbind, transform)
  coef_matrix = matrix(coeff, nrow = length(coeff), ncol = 1, byrow = TRUE)

  pred = X %*% coef_matrix
  pred[, 1]
}
