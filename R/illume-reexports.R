## The remedy generics are illumex's, shared between its data checks and
## illume's model checks: one ilm_remedies(), one ilm_apply_remedy(), one
## table and one print. illume registers its methods on them and exports
## them again, so illume::ilm_remedies() keeps working and attaching illume
## after illumex masks nothing -- the objects are the same.

#' @importFrom illumex ilm_remedies
#' @export
illumex::ilm_remedies

#' @importFrom illumex ilm_apply_remedy
#' @export
illumex::ilm_apply_remedy

#' @importFrom illumex ilm_remedy_table
#' @export
illumex::ilm_remedy_table

#' @importFrom illumex iml_remedies
#' @export
illumex::iml_remedies

#' @importFrom illumex iml_apply_remedy
#' @export
illumex::iml_apply_remedy

#' @importFrom illumex iml_remedy_table
#' @export
illumex::iml_remedy_table
