# Collection of functions to manage package versioning and installation
# Main package installation script will install all packages listed in package manifest
# as well as specific versions of GitHub packages. Upon completion, manifest is updated
# to include all currently installed packages.

.Rprofile$install_requirements <- function(Rprofile = .Rprofile, check_only = FALSE) {

  # Vectorized check for package installation status
  .installed <- function(pkgs) pkgs %in% rownames(utils::installed.packages())

  # Vectorized extraction of GitHub package commit SHA (i.e. version)
  .github_sha <- function(pkgs) sapply(pkgs, function(x) if (.installed(x)) utils::packageDescription(x)$GithubRef else 'not installed')

  # Install a package from GitHub given user name, package name, and specific commit version only if necessary. Requires devtools
  .install_from_github <- function(url, libs, .check_only = check_only) {
    user <- gsub('/.*$', '', url)
    pkgs <- gsub('^.*?/|@.*?$', '', url)
    sha  <- gsub('^.*?@', '', url)

    not_installed <- !.installed(pkgs)
    wrong_version <- sha != .github_sha(pkgs)
    needs_install <- url[not_installed | wrong_version]
    if (!.check_only) {
      lapply(needs_install, function(pkg) {
        withr::with_libpaths(unique(unlist(c(libs$GitHub, libs))), devtools::install_github(pkg, dependencies = T, upgrade_dependencies = F))
      })
    }
    return(needs_install)
  }

  # Install from CRAN style repositories
  .install_from_repos <- function(pkgs, lib, .check_only = check_only) {
    not_installed <- !.installed(pkgs)
    needs_install <- pkgs[not_installed]
    if (!.check_only) {
      lapply(needs_install, utils::install.packages, lib = lib)
    }
    return(needs_install)
  }

  CRAN <- .install_from_repos(Rprofile$CRAN, lib = Rprofile$Libraries$CRAN)
  Bioc <- .install_from_repos(Rprofile$Bioconductor, lib = Rprofile$Libraries$Bioconductor)
  GitH <- .install_from_github(Rprofile$GitHub, libs = Rprofile$Libraries)
  INSTALLED <- c(CRAN, Bioc, GitH)

  if (!check_only & length(INSTALLED)) {
    message('If all packages installed successfully, please run .Rprofile$load_requirements()')
  }

  return(invisible(INSTALLED))
}
