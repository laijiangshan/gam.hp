
.onAttach <- function(libname, pkgname) {
  cite_info <- utils::citation(pkgname)
  cite_info <- "Jiangshan Lai, Jing Tang, Tingyuan Li, Aiying Zhan, Lingfeng Mao. 2024. Evaluating the relative importance of predictors in Generalized Additive Models using the gam.hp R package. Plant Diversity, 46(4): 542-546"
  packageStartupMessage("Thank you for using this package! If you use this package in your research, please cite the following references \n")
  packageStartupMessage(cite_info)
}
