if (.Rprofile$R_version != with(R.version, paste0(major, '.', minor))) {
  warning(
    '\nThe current R version (', with(R.version, paste0(major, '.', minor)),
    ") does not match this project's expected R version (", .Rprofile$R_version,
    '). If you encounter problems, consider installing the expected R version: \n',
    '  Mac:   https://cran.r-project.org/bin/macosx/old/ (see also https://r.research.att.com/RSwitch-1.2.dmg) \n',
    '  Win:   https://cran.r-project.org/bin/windows/base/old/ \n',
    '  Linux: https://cran.r-project.org/bin/linux/ \n')
}
