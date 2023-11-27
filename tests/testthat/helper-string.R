# A replacement for stringr::str_count()
str_count <- function(x, pattern, fixed = TRUE) {
  loc <- gregexpr(pattern, text = x, fixed = fixed)
  unlist(lapply(loc, function(x) length(x[x > 0])))
}

# A replacement for stringr::str_detect()
str_detect <- function(x, pattern) {
  grepl(pattern, x = x)
}
