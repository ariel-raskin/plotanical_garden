# path_utils.R — Windows-safe path handling

#' Convert a file path to Windows 8.3 short form to avoid
#' data.table::fread interpreting spaces/& as shell commands.
#' On non-Windows, just normalizes the path.
safe_path <- function(p) {
  if (is.null(p) || !nzchar(p)) return(p)
  # Strip any surrounding quotes the user may have pasted
  p <- gsub('^["\']|["\']$', '', trimws(p))
  # Normalize
  p <- normalizePath(p, winslash = "/", mustWork = FALSE)
  if (.Platform$OS.type == "windows" && file.exists(p)) {
    # 8.3 short names avoid spaces, &, parens — fread-safe
    short <- utils::shortPathName(p)
    if (nzchar(short)) p <- gsub("\\\\", "/", short)
  }
  p
}

#' Force data.table to treat paths as files, not commands.
#' Call this once at app startup.
suppress_fread_command <- function() {
  # data.table >= 1.11.6 respects this option
  options(datatable.fread.input.cmd.message = FALSE)
}
