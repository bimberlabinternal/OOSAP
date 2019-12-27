#' Table mapping gene symbol to CD gene name
#'
#'
#' @format A data frame with 394 rows and 4 variables:
#' \describe{
#'   \item{GeneSymbol}{Official GeneSymbol}
#'   \item{ApprovedName}{A description}
#'   \item{PreviousSymbols}{Previous symbols.  This is the field typically used for aliasing gene symbols}
#'   \item{Synonyms}{Additional aliases}
#'   ...
#' }
#' @source \url{https://www.genenames.org/data/genegroup/#!/group/471}
"cdGenes"

#' Table mapping ENS to SYM for MMul10
#'
#'
#' @format A data frame with 35427 rows and 2 variables:
#' \describe{
#'   \item{ENSID}{ENS ID}
#'   \item{SYMB}{Gene Symbol}
#'   ...
#' }
#' @source \url{https://}
"Mmul10_conv1"

#' Table mapping ENS to SYM for MMul8
#'
#'
#' @format A data frame with 32386 rows and 2 variables:
#' \describe{
#'   \item{ENSID}{ENS ID}
#'   \item{SYMB}{Gene Symbol}
#'   ...
#' }
#' @source \url{https://}
"Mmul8_conv1"