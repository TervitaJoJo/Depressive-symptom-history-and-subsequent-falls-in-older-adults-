# Generated data, imputation caches and logs must stay outside this package.
script_path <- sub("^--file=", "", commandArgs()[grep("^--file=", commandArgs())][1])
package_root <- normalizePath(file.path(dirname(script_path), "../.."), winslash="/", mustWork=TRUE)
assert_external <- function(path) {
  path <- gsub("\\", "/", path, fixed=TRUE)
  if(!grepl("^([A-Za-z]:/|/)",path)) path <- file.path(getwd(),path)
  resolved <- normalizePath(path,winslash="/",mustWork=FALSE)
  resolved <- gsub("\\", "/", resolved, fixed=TRUE)
  a <- tolower(sub("/+$","",resolved)); b <- tolower(sub("/+$","",package_root))
  if(a==b || startsWith(paste0(a,"/"),paste0(b,"/"))) stop("Use a private working directory outside the submission package.")
  invisible(path)
}
