#' log memory usage
#'
#' Logs the amount of memory being used to a log file if it is available, and
#' generating warnings if the amount of RAM hits zero.
#'
#' @export
#' @return NULL
log_memory = function(){
  linux_memory = system("cat /proc/meminfo", intern = TRUE)
  linux_memory = grep("^Mem", linux_memory, value = TRUE)

  memory_values = stringr::str_extract(linux_memory, "[[:digit:]].*")
  memory_numbers = as.numeric(stringr::str_extract(memory_values, "[[:digit:]].* "))
  memory_string = paste0("Memory: ", paste(paste(c("Total: ", "Free: ", "Available: "), memory_values, sep = ""), collapse = ", "))

  if (any(memory_numbers == 0)) {
    if (get("logger", envir = has_logger)) {
      logger::log_warn(memory_string, namespace = "smirfeTools")
    } else {
      warning(memory_string)
    }
  } else {
    if (get("logger", envir = has_logger)) {
      logger::log_info(memory_string, namespace = "smirfeTools")
    }
  }
}

#' log messages
#'
#' If a log_appender is available, logs the given message at the `info` level.
#'
#' @param message_string the string to put in the message
#'
#' @export
#' @return NULL
log_message = function(message_string){
  if (get("logger", envir = has_logger)) {
    logger::log_info(message_string, namespace = "smirfeTools")
  }
}

#' disable logging
#'
#' In certain situations it may be good to disable logging altogether. This allows
#' that.
#'
#' @export
#' @return NULL
disable_logging = function(){
  assign("logger", FALSE, envir = has_logger)
  message("Logging disabled.char")
}