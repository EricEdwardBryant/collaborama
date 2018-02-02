as.sgEdit <- function(x, edit = 'c') {
  stopifnot(all(str_detect(x, edit)))
  attr(x, 'target') <- edit
  class(x) <- 'sgEdit'
  return(x)
}

sgEdit <- function(..., edit = 'c') as.sgEdit(c(...), edit)

print.sgEdit <- function(x, ...) {
  edit <- str_locate(x, attr(guide, 'target'))[,'start']
  
  lead   <- str_c(str_dup(' ', 8 - edit), str_sub(x, 1, 3))
  upwind <- str_sub(x, 4, edit - 1L) %>% crayon::blue$underline()
  target <- str_sub(x, edit, edit)   %>% crayon::green$bold$underline()
  dnwind <- str_sub(x, edit + 1L, 8) %>% crayon::blue$underline()
  spacer <- str_sub(x, 9, 20)
  pam    <- str_sub(x, 21)           %>% crayon::red()
  cat(str_c(lead, upwind, target, dnwind, spacer, pam), sep = '\n')
}

guide <- sgEdit('ATGCTCcATAGATCATATTAAGG', 'TGCTCcATAGATCATATTAAGGG')
guide
