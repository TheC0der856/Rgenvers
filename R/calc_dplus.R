#' Calculate Dplus for Allele Abundances
#'
#' This function calculates the Dplus diversity metric using the
#' \code{taxondive()} function from the vegan package.
#'
#' @param genpop A genpop object (from the adegenet package) containing an additional row "pop" containing the KBA in which the individual occurs.
#'
#' @return A numeric vector containing the Dplus values for each KBA.
#' @import vegan
#' @export
#'
#' @examples
#'   data("nancycats", package = "adegenet")
#'   df <- adegenet::genind2df(nancycats)
#'   pop <- as.vector(df$pop)
#'   gen <- adegenet::df2genind(df[2:ncol(df)], pop = pop, ncode = 3)
#'   catpop <- adegenet::genind2genpop(gen)
#'   result <- calc_dplus(catpop)
#'   print(result)

calc_dplus <- function(genpop) {
  allele_abundances <- genpop@tab
  t_allele_abundances <- t(allele_abundances)
  taxdis <- vegan::taxa2dist(t_allele_abundances)
  mod <- vegan::taxondive(allele_abundances, taxdis)
  return(mod$Dplus)
}
