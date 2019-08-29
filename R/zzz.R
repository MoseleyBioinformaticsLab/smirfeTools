#' pick map enumerator
#'
#' Allows the user to set which enumerator is being used internally in the functions.
#'
#' @param map_function which function to use, assigns it to an internal object
#'
#' @export
#' @return NULL
set_internal_map <- function(map_function = NULL){
  if (is.null(map_function)) {
    assign("map_function", purrr::map, envir = internal_map)
  } else {
    assign("map_function", map_function, envir = internal_map)
  }
}


internal_map <- new.env(hash = TRUE)
assign("map_function", purrr::map, envir = internal_map)
has_logger = new.env(hash = TRUE)
assign("logger", FALSE, envir = has_logger)
assign("memory", FALSE, envir = has_logger)

.onLoad <- function(libname, pkgname) {
  tmp_packages = installed.packages()
  if ("logger" %in% rownames(tmp_packages)) {
    assign("logger", TRUE, envir = has_logger)
    sys_info = Sys.info()
    if (grepl("linux", sys_info["sysname"], ignore.case = TRUE)) {
      assign("memory", TRUE, envir = has_logger)
    }
    logger::log_appender(logger::appender_file(paste0("smirfeTools_run_", substring(make.names(Sys.time()), 2), ".log")), namespace = "smirfeTools")
  }
  debugme::debugme()
}