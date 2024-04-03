#' Display debug messages if verbose allows it
#'
#' A function that extend rlang::inform
#' to display a message if the verbose is at "debug"
#' and show content of the variable
#' @param message String to display
#' @param ... Other parameters for rlang::inform
#' @return Nothing
debug_msg <- function(message = NULL, ...) {
    is_debug_mode <- (getOption("BREW3R.r.verbose", "quiet") == "debug")
    if (is_debug_mode) {
        rlang::local_options(rlib_message_verbosity = "verbose")
        rlang::inform(message = message, ...)
        if (exists(message)) {
            methods::show(eval(parse(text = message)))
        }
    }
}
#' Display progression messages if verbose allows it
#'
#' A function that extend rlang::inform
#' to display a message if the verbose is at "debug" or "progression"
#' @param ... Parameters for rlang::inform
#' @return Nothing
progression_msg <- function(...) {
    is_progression_or_debug_mode <-
        (getOption("BREW3R.r.verbose", "quiet") %in% c("debug", "progression"))
    if (is_progression_or_debug_mode) {
        rlang::local_options(rlib_message_verbosity = "verbose")
        rlang::inform(...)
    }
}
