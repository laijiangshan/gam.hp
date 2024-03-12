#' Plot for a \code{\link{gam.hp}} object
#'
#' @param x A \code{\link{gam.hp}} object.
#' @param  plot.perc Logical;if TRUE, the bar plot (based on ggplot2 package) of the percentage to individual effects of variables towards total explained variation, the default is FALSE to show plot with original individual effects.
 
#' @param ... unused
#' @return a ggplot object
#' @author {Jiangshan Lai} \email{lai@njfu.edu.cn}
#' @export
#' @examples
#'library(mgcv)
#'mod1 <- gam(Sepal.Length ~ s(Petal.Length) + s(Petal.Width) + Sepal.Width,data = iris)
#'plot(gam.hp(mod1))

plot.gamhp <- function(x, plot.perc = FALSE,...){
  if (!inherits(x, "gamhp")){
    stop("x should be the output of gam.hp()")
  }
if(x$type=="hierarchical.partitioning")
  {if (plot.perc){
    tips3 = data.frame(variable = rownames(x[[2]]),
                       value = as.numeric(x[[2]][, "I.perc(%)"]))
    gg = ggplot2::ggplot(tips3, ggplot2::aes(x = stats::reorder(variable,-value), y = value)) + ggplot2::geom_bar(stat = "identity") +
      ggplot2::theme_minimal() + ggplot2::labs(x = "Variables", y = "% Individual effect to Total (%I)") + ggplot2::theme(axis.text = element_text(size = 10)) + ggplot2::theme(axis.title = element_text(size = 13))+ggplot2::labs(title =paste(names(x)[1]))
  } else {
    tips2 = data.frame(variable = rownames(x[[2]]),
                       value = as.numeric(x[[2]][, "Individual"]))
    gg = ggplot2::ggplot(tips2, ggplot2::aes(x = stats::reorder(variable,-value), y = value)) + ggplot2::geom_bar(stat = "identity") +
      ggplot2::theme_minimal() + ggplot2::labs(x = "Variables", y = "Individual effect") + ggplot2::theme(axis.text = element_text(size = 10)) + ggplot2::theme(axis.title = element_text(size = 13))+ggplot2::labs(title =paste(names(x)[1]))
  }
return(gg)
}
if(x$type=="commonality.analysis")
{ 
  stop("Commonality analysis plot is unavailable in current version")
}
}