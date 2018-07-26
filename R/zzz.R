#' pick map_at enumerator
#'
#' Allows the user to set which enumerator is being used internally in the functions.
#'
#' @param map_function which function to use, assigns it to an internal object
#'
#' @export
#' @return NULL
set_internal_map <- function(map_function = NULL){
  if (is.null(map_function)) {
    assign("map_function", purrr::map_at, envir = internal_map)
  } else {
    assign("map_function", map_function, envir = internal_map)
  }
}


internal_map <- new.env(hash = TRUE)
assign("map_function", purrr::map_at, envir = internal_map)
