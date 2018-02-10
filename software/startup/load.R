# Function assigned to `.Rprofile`. Loads packages specified in `.Rprofile$Load`
# and sources directories to objects specified in `.Rprofile$Source_dirs`.
# Session information is then printed in an attractive format. This printed
# output is suitable for use as session input information (though setting
# `show_loaded = TRUE` is recommended in this case)

.Rprofile$load_requirements <- function(Rprofile = .Rprofile, show_system = TRUE, show_attached = TRUE, show_loaded = FALSE) {

  # Check for required packages specified in .Rprofile
  needs_install <- .Rprofile$install_requirements(check_only = TRUE)
  if (length(needs_install)) {
    message(
      'This project requires the following package(s):\n  ',
      paste0(needs_install, collapse = '\n  '),
      '\nPlease install them by running .Rprofile$install_requirements()',
      '\nOnce all requirements have been installed, run .Rprofile$load_requirements()'
    )
    return(invisible())
  }

  .source_directory <- function(from, ext = '[.][Rr]$') {
    invisible(lapply(list.files(from, ext, full.names = T), source))
  }

  # Load libraries and source directories
  suppressPackageStartupMessages(lapply(Rprofile$Load, library, character.only = TRUE, warn.conflicts = FALSE))
  lapply(Rprofile$Source_dirs, .source_directory)

  # Everything below here is just for fun printing
  package_version <- function(x) utils::packageDescription(x)$Version

  repo_version <- function(x) {
    d <- utils::packageDescription(x)
    if (!is.null(d$GithubRef)) {
      return(paste0('(Github.com/', d$GithubUsername, '/', d$GithubRepo, '@', d$GithubRef, ')'))
    } else if (!is.null(d$biocViews)) {
      return(paste0('(Bioc ', Rprofile$Bioc_version, ')'))
    } else {
      return(paste0('(CRAN ', Rprofile$CRAN_date, ')'))
    }
  }

  session <- utils::sessionInfo()

  if (show_system) {
    message(cli::rule(left = crayon::bold(with(session$R.version, paste0('R ', major, '.', minor, ' "', nickname, '"'))), width = 80))
    info <- c(session$running, session$platform, normalizePath(Rprofile$Library_root))
    fields <- c('System', 'Platform', 'Library')
    msg <-
      paste0(
        crayon::magenta(cli::symbol$bullet), ' ',
        crayon::blue(format(fields)), ' ',
        crayon::col_align(info, max(crayon::col_nchar(info)))
      )
    message(paste(msg, collapse = '\n'))
  }

  if (show_attached) {
    message(cli::rule(left = crayon::bold("Packages attached"), width = 80))
    attached <- rev(names(session$otherPkgs))
    versions <- vapply(attached, package_version, character(1))
    repos    <- vapply(attached, repo_version, character(1))
    i <- 1:floor(length(attached) / 2)
    col1 <-
      paste0(
        crayon::green(cli::symbol$tick), " ",
        crayon::blue(format(attached[i])), " ",
        crayon::col_align(versions[i], max(crayon::col_nchar(versions[i]))), " ",
        crayon::col_align(repos[i], max(crayon::col_nchar(repos[i])))
      )
    col2 <-
      paste0(
        crayon::green(cli::symbol$tick), " ",
        crayon::blue(format(attached[-i])), " ",
        crayon::col_align(versions[-i], max(crayon::col_nchar(versions[-i]))), " ",
        crayon::col_align(repos[-i], max(crayon::col_nchar(repos[-i])))
      )

    info <- paste0(col1, "   ", col2)
    message(paste(info, collapse = "\n"))
  }


  if (show_loaded) {
    message(cli::rule(left = crayon::bold("Packages loaded via namespace"), width = 80))
    loaded <- rev(names(session$loadedOnly))
    versions <- vapply(loaded, package_version, character(1))
    repos    <- vapply(loaded, repo_version, character(1))
    i <- 1:floor(length(loaded) / 2)
    col1 <-
      paste0(
        crayon::green(cli::symbol$tick), " ",
        crayon::blue(format(loaded[i])), " ",
        crayon::col_align(versions[i], max(crayon::col_nchar(versions[i]))), " ",
        crayon::col_align(repos[i], max(crayon::col_nchar(repos[i])))
      )
    col2 <-
      paste0(
        crayon::green(cli::symbol$tick), " ",
        crayon::blue(format(loaded[-i])), " ",
        crayon::col_align(versions[-i], max(crayon::col_nchar(versions[-i]))), " ",
        crayon::col_align(repos[-i], max(crayon::col_nchar(repos[-i])))
      )

    info <- paste0(col1, "   ", col2)
    message(paste(info, collapse = "\n"))
  }
  invisible()
}

if (.Rprofile$Auto_load) {
  cat('\014')
  message('Copyright (C) The R Foundation for Statistical Computing')
  message('R is free software and comes with ABSOLUTELY NO WARRANTY')
  message(
    '  Help:         help(), and ?<function>\n',
    '  Demo:         demo()\n',
    '  Quit:         q()\n',
    '  Licence:      licence()\n',
    '  Contributors: contributors()\n',
    '  Citing R:     citation()\n'
  )
  .Rprofile$load_requirements()
}
