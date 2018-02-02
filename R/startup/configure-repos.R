.Rprofile$Repositories <- c(
  CRAN      = file.path('https://mran.microsoft.com/snapshot', .Rprofile$CRAN_date),
  BioCsoft  = file.path('https://bioconductor.org/packages',   .Rprofile$Bioc_version, 'bioc'),
  BioCann   = file.path('https://bioconductor.org/packages',   .Rprofile$Bioc_version, 'data/annotation'),
  BioCexp   = file.path('https://bioconductor.org/packages',   .Rprofile$Bioc_version, 'data/experiment'),
  BioCextra = file.path('https://bioconductor.org/packages',   .Rprofile$Bioc_version, 'extra')
)

options(repos = .Rprofile$Repositories, keep.source.pkgs = TRUE)
