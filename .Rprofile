# ==== Profile options ========================================================
.Rprofile <-
  list(
    # Version specification
    R_version    = '3.4.3',      # Not strict - just warns if mismatch
    Bioc_version = '3.5',        # Bioconductor has fixed releases
    CRAN_date    = '2017-11-20', # Date of CRAN snapshot (uses MRAN repository)
    Library_root = 'software/library',  # where to install packages
    # Packages to install from CRAN, Bioconductor and GitHub
    CRAN = c(
      'tidyverse', # Collection of data manipulation and visualization packages
      'devtools',  # Used to install packages from GitHub
      'fuzzyjoin', # Usefull genome_join
      'pbapply',   # Long loops with progress bars and parallel processing
      'rmarkdown', # For reproducible reporting
      'caTools'    # Required by rmarkdown in RStudio
    ),
    Bioconductor = c(
      'BiocInstaller',               # Install packages from Bioconductor
      'BSgenome',                    # Biostrings genome object infrastructure
      'BSgenome.Hsapiens.UCSC.hg38', # Human genome reference GRCh38
      'BSgenome.Hsapiens.UCSC.hg19'  # Human genome reference GRCh37
    ),
    GitHub = c(
      'EricEdwardBryant/iSTOP@a2fd5d3',      # The original iSTOP tool
      'EricEdwardBryant/mutagenesis@0fe642a' # Variant effect prediciton
    ),
    # Packages to load (in order) when using .Rprofile$load_requirements
    Load = c(
      'methods',
      'stats',
      'utils',
      'tidyverse',
      'iSTOP',
      'mutagenesis'
    ),
    # Directories to source when using .Rprofile$load_requirements
    Source_dirs = 'R',
    # Should loading occur automatically?
    Auto_load = interactive()  # only if interactive (i.e. not in reports)
  )
# =============================================================================

source('software/startup/configure-libs.R')  # as per .Rprofile$Library_root
source('software/startup/configure-repos.R') # As per .Rprofile
source('software/startup/install.R')         # .Rprofile$install_requirements()
source('software/startup/load.R')            # .Rprofile$load_requirements()
source('software/startup/check-r-version.R') # As per .Rprofile$R_version
