#' log memory usage
#'
#' @param log_obj the log object to use to log the memory
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
      logger::log_error(memory_string, namespace = "smirfeTools")
    } else {
      warning(memory_string)
    }
  } else {
    if (get("logger", envir = has_logger)) {
      logger::log_info(memory_string, namespace = "smirfeTools")
    }
  }
}
