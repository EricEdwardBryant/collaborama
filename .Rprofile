# ==== Profile options ========================================================
.Rprofile <-
  list(
    # Version specification
    R_version    = '3.4.3',      # Not strict - just warns if mismatch
    Bioc_version = '3.5',        # Bioconductor has fixed releases
    CRAN_date    = '2017-11-20', # Date of CRAN snapshot (uses MRAN repository)
    Library_root = 'R/library',  # where to install packages
    # Packages to install from CRAN, Bioconductor and GitHub
    CRAN = c(
      
      'tidyverse', # Collection of data manipulation and visualization packages
      'devtools',  # Used to install packages from GitHub
      'fuzzyjoin',
      'pbapply',
      'rmarkdown', # For reproducible reporting
      'caTools',   # Required by rmarkdown in RStudio
      'R.utils'
    ),
    Bioconductor = c(
      'BiocInstaller',               # Install packages from Bioconductor
      'BSgenome',                    # Biostrings genome object infrastructure
      #'Biostrings',                  # DNA / RNA / Protein sequence manipulation
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
    # Directories to source as objects when using .Rprofile$load_requirements
    Source_dirs = c(
      data_import = 'R/data-import',
      analyses    = 'R/analyses',
      figures     = 'R/figures',
      fxn         = 'R/fxn'
    ),
    # Should loading occur automatically?
    Auto_load = interactive()  # only if interactive (i.e. not in reports)
  )
# =============================================================================

source('R/startup/check-r-version.R') # Warning if not .Rprofile$R_version
source('R/startup/configure-libs.R')  # As specified by .Rprofile$Library_root
source('R/startup/configure-repos.R') # As specified in .Rprofile
source('R/startup/install.R')         # Create .Rprofile$install_requirements()
source('R/startup/load.R')            # Create .Rprofile$load_requirements()
