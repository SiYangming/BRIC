#' @describeIn QUBIC Performs a QUalitative BIClustering.
#'
#' @usage qubic(i, R = FALSE, F = FALSE, d = FALSE, f = 0.85, k = 13, c = 0.90, o = 5000)
#' 
#' @export
qubic <- function(i, R = FALSE, F = FALSE, d = FALSE, f = 0.85, k = 13, c = 0.90, o = 5000) {
  vec <- c("./qubic", "-i", i)
  if(R) vec <- c(vec, "-R")
  if(F) vec <- c(vec, "-F") 
  if(d) vec <- c(vec, "-d") 
  vec <- c(vec, "-f", as.character(f))
  vec <- c(vec, "-k", as.character(k))
  vec <- c(vec, "-c", as.character(c))
  vec <- c(vec, "-o", as.character(o))
  
  return(.main(vec))
}
