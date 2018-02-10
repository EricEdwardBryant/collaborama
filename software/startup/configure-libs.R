.Rprofile$Library_here <- with(.Rprofile, file.path(Library_root, R_version, R.version$platform))

.Rprofile$Libraries <-
  list(
    CRAN         = with(.Rprofile, file.path(Library_here, 'CRAN', CRAN_date)),
    Bioconductor = with(.Rprofile, file.path(Library_here, 'Bioconductor', Bioc_version)),
    GitHub       = with(.Rprofile, file.path(Library_here, 'GitHub'))
  )

# Create directories if not present
invisible(lapply(.Rprofile$Libraries, dir.create, showWarnings = F, recursive = T))

.libPaths(unlist(.Rprofile$Libraries)) # Configures libraries
