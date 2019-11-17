#' Calculate exact masses
#' 
#' @param form A character with compound formula
#' @param addct A character with adduct (if missing, neutral mass is returned)
#' @examples
#' exact_mz("H2O")
#' exact_mz("H2O", "[M+H]+")
exact_mz <- function(form = character(0), addct = NA){
  mz_neut <- getMolecule(form)$exactmass
  if(is.na(addct)){
    mz_neut
  }else{
    unlist(mass2mz(mz_neut, adduct = addct))
  }
}
