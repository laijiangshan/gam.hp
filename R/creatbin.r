#'Internal function for glmm.hp() to create diagonal matrix
#' @param  col Imput number.
#' @param  binmatrix Imput empty matrix.

#' @return a matrix
#' @return \item{a matix}{A diagonal matrix}
#' @keywords internal
creatbin <-
function(col, binmatrix) {
row<-1
val<-col
	while (val!=0){              
	if (odd(val)) {
		binmatrix[row,col]=1 
	}
	val<-as.integer(val/2)
	row<-row+1
}
##Return matrix
return(binmatrix)
}


#'Internal function for glmm.hp() to determine whether the odd number
#' @param  val Imput number.
#' @return a logical value
#' @return \item{Logical value}{TRUE or FALSE}
#' @keywords internal
odd <- function (val) 
{
    if (val%%2 == 0) 
     return(FALSE)
	 else
    return(TRUE)
}