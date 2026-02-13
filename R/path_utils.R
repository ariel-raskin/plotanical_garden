# path_utils.R — Windows-safe path handling

#' Get Windows 8.3 short path (avoids &, spaces, parens issues)
#' Falls through to normalizePath on non-Windows
safe_path <- function(p) {
  if (is.null(p) || !nzchar(p)) return(p)
  p <- normalizePath(p, winslash = "/", mustWork = FALSE)
  if (.Platform$OS.type == "windows" && file.exists(p)) {
    # utils::shortPathName converts to 8.3 format, avoiding special chars
    p <- utils::shortPathName(p)
    p <- gsub("\\\\", "/", p)
  }
  p
}
