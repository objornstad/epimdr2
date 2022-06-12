#' A function to caculate asymptotic growth, senitivy and elasticity for age-structured populations
#' @param L the Leslie matrix
#' @return A list consisting of the following components: 
#' \item{lambda}{the dominant eigen value of the Leslie matrix.}
#' \item{right.eigenvector}{the dominant right eigen vector of the Leslie matrix, proportional to the stable age-distribution.}
#' \item{left.eigenvector}{the dominant left eigen vector of the Leslie matrix representing the age-specific reproductive values.}
#' \item{elasticity}{the elasticities.}
#' \item{sensitivity}{the sensitivities.}
#' @examples
#' fa<-c(0, 0.5, 1.2)
#' sa<-c(0.8, 0.8, 0)
#' L<-matrix(0, nrow=3, ncol=3)
#' #inserting fa vector in first row
#' L[1,]<-fa
#' #inserting sa in the subdiagonal:
#' L[row(L)==col(L)+1] <-sa[1:2]
#' leslie(L)
#' @references Caswell, H. 2001.  Matrix Population Models: Construction, Analysis, and Interpretation. 2nd edn Sinauer Associates Inc., Sunderland, MA, 
#' @export
leslie = function(L){
   Ex = eigen(L)    #Eigendecompostition of matrix
   w = Re(Ex$vectors[, 1]) #right eigenvector
   lambda = Re(Ex$values[1])  #dominant eigenvalue
   v = Re(eigen(t(L))$vectors[, 1])    #left eigenvector 
   sens = (v %*% t(w))/sum(v * w)   #sensitivities
   elast = L  * sens/lambda      #elasticities
   #list of results
   res = list(lambda = lambda,    
   right.eigenvector = w,    
   left.eigenvector = v,   
   elasticity = elast,                       
   sensitivity = sens            
   )
   return(res)
}