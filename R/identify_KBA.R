#' Identify KBA(s) with distinct genetic diversity
#'
#' This function identifies KBA with distinct genetic diversity for different criteria
#'
#' @param genpop A genpop object (from the adegenet package) containing an additional row "pop" containing the KBA in which the individual occurs.
#' @param criterium A string specifying the criterion to apply (A1: Threatened species (A1a: Critically endangered/ Endangered species, A1b: Vulnerable species, A1c: Critically endangered/ Endangered species due only to past/current decline [Red List A only, but not A3 only], A1d: Vulnerable species due only to past/current decline [Red List A only, but not A3 only]), B1: Individual geographically restricted species)
#'
#' @return A vector containing the names of the identified KBA(s). The vector is empty if no KBA meets the criterium (character(0)).
#' @import vegan
#' @export
#'
#' @examples
#'   data(nancycats, package = "adegenet")
#'   df <- adegenet::genind2df(nancycats)
#'   pop <- as.vector(df$pop)
#'   gen <- adegenet::df2genind(df[2:ncol(df)], pop = pop, ncode = 3)
#'   catpop <- adegenet::genind2genpop(gen)
#'   result <- identify_KBA(catpop, "A1a")
#'   print(result)
#'
identify_KBA <- function(genpop, criterium) {
  thresholds <- list(
    A1a = 0.5,
    A1b = 1,
    A1c = 0.1,
    A1d = 0.2,
    B1 = 10
  )
  if (!(criterium %in% names(thresholds))) {
    stop("unknown criterium: ", criterium)
  }
  threshold <- thresholds[[criterium]]
  allele_abundances <- genpop@tab
  t_allele_abundances <- t(allele_abundances)
  taxdis <- vegan::taxa2dist(t_allele_abundances)
  mod <- vegan::taxondive(allele_abundances, taxdis)
  sumDplus <- sum(mod$Dplus)
  ratios <- mod$Dplus/sumDplus * 100
  KBA_ratios <- ratios[ratios >= threshold]
  return(names(KBA_ratios))
}
