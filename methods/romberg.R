source('C:/projects/NumericalAnalysis/helpers/draw.R')
source('C:/projects/NumericalAnalysis/methods/newtonCotes.R')

romberg <- function(f, a, b, n) {
  print("Aplicando metodo de Romberg")
  
  h <- (b-a)/n
  
  kArray <- 1:n
  hArray <- (b-a)/(2^(kArray-1))
  
  R <- matrix(data = NA, nrow = n, ncol = n, byrow = FALSE, dimnames = NULL)
  
  for(j in kArray) {
    if (j == 1) 
      for(k in kArray) 
        if (k == 1) {
          R[1, 1] <- trapezoid(f, a, b, 1, FALSE)[[1]]
        } else {
          i <- seq(1, 2^(k-2), by = 1)
          R[k, 1] <- 1/2 * (R[k-1, 1] + hArray[k-1] * sum(f(a + (2*i-1)*hArray[k])))  
        }
    
    if (j > 1) 
      for (k in kArray) {
        value <- R[k, j-1] + (R[k, j-1] - R[k-1, j-1])/(4^(j - 1) - 1)
        
        if (length(value) > 0)
          R[k, j] <- value
      }
  }
  
  errors <- R[n, ][2:n] - R[n, ][1:(n-1)]
  
  draw.integrate(f, a, b)
  
  print(R)
  c(R[n,n], min(errors))
}
