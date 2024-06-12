
#' @title Internal Function for Creating Directories
#'
#' @param path the directory to create
#'
#' @inheritDotParams dir.create
#'
#' @return see [dir.create]
checked_dir_create <- function(path, ...) {
  r <- readline(prompt = sprintf(
    "Try to create path '%s' ? (y/anything else == no) ", path
  ))
  if ((r != "") && switch(r, Y = , y = TRUE, FALSE)) {
    dir.create(path, ...)
  } else {
    invisible(FALSE)
  }
}

#' @title Internal Function for Confirming Overwrite
#'
#' @param path the directory to create
#'
#' @return logical; `TRUE` if the overwrite can proceed, `FALSE` if not
allow_overwrite <- function(path) {
  return(readline(prompt = paste0(
    "Allow overwrite of files in '", path, "'? (y/n) "
  )) == "y")
}

#' Create a local copy of the tutorials.
#'
#' @param path, string; the path to enclosing directory. If this directory does
#'   not exist, will create it
#'
#' @param overwrite, logical; overwrite existing files (corresponding to the
#'   course scripts)?
#'
#' @param what, string: the "scripts" or the "solutions"?
#'
#' @details This function creates a local copy of the R scripts associated with
#' the ICI3D short courses, for students to edit and run during learning
#' activities. These scripts should only require the `ICI3D` package and
#' associated dependencies.
#'
#' The created directory contains several subdirectories, each corresponding to
#' a different session. Within each directory, the exercise files are generally
#' named `NN_short_description.R` (where NN is 00, 01, etc - corresponding to
#' the expected order of completion, with the 00 session being background /
#' preparatory content).
#'
#' If used with default `what` argument, you will get a "scripts" subdirectory,
#' and if you ask for `what = "solutions"` there will be a "solutions"
# subdirectory. Both will contain sessions by name, scripts by order and
# shortname, without and with solutions filled in, respectively.
#'
#' @return string or `NULL`; non-`NULL` indicates successful creation + copy and
#' is the top level root of the course material.
#'
#' @examples
#' require(MTM)
#' tardir <- scripts()
#' list.files(tardir, recursive = TRUE)
#'
#' @export
scripts <- function(
    path = file.path(
      if (.Platform$OS.type == "windows") {
        Sys.getenv("USERPROFILE")
      } else { "~" },
      "Downloads", "ICI3D"
    ),
    overwrite = FALSE,
    what = c("scripts", "solutions")
) {
  stopifnot(
    "'path' must be a string." = is.character(path),
    "'path' must be a single string." = length(path) == 1
  )

  refpath <- match.arg(what, several.ok = FALSE)
  path <- file.path(path, refpath)

  dir_exists <- dir.exists(path)
  if (!dir_exists && !checked_dir_create(path, recursive = TRUE)) {
    warning(sprintf("'%s' does not exist and/or was not created.", path))
    return(invisible())
  }

  # Check user-requested overwrite
  if (
    overwrite && dir_exists && !allow_overwrite(path)
  ) {
    return(NULL)
  } # else perform the overwrite

  srcdir   <- system.file(refpath, package = "ICI3D")
  srcfiles <- list.files(
    srcdir, full.names = TRUE, recursive = FALSE, include.dirs = TRUE
  )

  res <- file.copy(from      = srcfiles,
                   to        = path,
                   overwrite = overwrite,
                   recursive = TRUE)

  # `file.copy` returns `FALSE` if `overwrite == FALSE` and the directory
  # exists, and we don't want that to trigger a warning; hence the
  # `(overwrite || !dir_exists)` here
  if ((overwrite || !dir_exists) && !all(res)) {
    warning("Something may have gone wrong with copying.")
    return(invisible())
  } else {
    return(path)
  }
}
