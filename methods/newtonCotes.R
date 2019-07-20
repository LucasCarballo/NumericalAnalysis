install.packages("Deriv")
require(Deriv)

trapezoid <- function(f, a, b) {
  h <- (b-a)
  
  t <- h/2 * (f(a) + f(b))
  
  f2dx <- Deriv(f, 'x', nderiv = 2)

  error <- 1/12 * (h^3) * fminbnd(f2dx, a, b, maximum = TRUE)[[2]]
  
  c(t, error)
}

trapezoidCompound <- function(f, a, b, n) {
  h <- (b-a)/n
  
  xi <- seq(a+h, b-h, by = h)
  
  t <- h/2 * (f(a) + f(b) + 2*sum(f(xi)))
  
  f2dx <- Deriv(f, 'x', nderiv = 2)
  
  error <- (b-a) * (h^2)/12 * fminbnd(f2dx, a, b, maximum = TRUE)[[2]]
  
  c(t, error)
}

simpson13compound <- function(f, a, b, n) {
  if (n%%2 != 0) stop("En la regla de Simpson, n tiene que ser par")
  
  h <- (b-a)/n
  
  i1 <- seq(1, n-1, by = 2) # impares
  i2 <- seq(2, n-2, by = 2) # pares
  
  s <- h/3 * ( f(a) + f(b) + 4*sum( f(a+i1*h) ) + 2*sum( f(a+i2*h) ) )
  
  f4dx <- Deriv(f, 'x', nderiv = 4)
  
  error <- (b-a) * (h^4)/180 * fminbnd(f4dx, a, b, maximum = TRUE)[[2]]
  
  c(s, error)
}

simpson38 <- function(f, a, b) {
  h <- (b-a)/3
  
  s <- 3*h/8 * (f(a) + 3*f((2*a+b)/3) + 3*f((a+2*b)/3) + f(b))
  
  f4dx <- Deriv(f, 'x', nderiv = 4)
  
  error <- -3/80 * h^5 * fminbnd(f4dx, a, b, maximum = TRUE)[[2]]
  
  c(s, error)
}

simpson38compound <- function(f, a, b, n) {
  if (n%%3 != 0) stop("En la regla de Simpson de 3/8 Compuesta, n tiene ser multiplo de 3")
  
  h <- (b-a)/n
  
  i1 <- seq(1, n-2, by = 3)
  i2 <- seq(2, n-1, by = 3)
  i3 <- seq(3, n-3, by = 3)
  
  s <- 3*h/8 * (f(a) + 3 * sum(f(a+i1*h)) + 3 * sum(f(a+i2*h)) + 2 * sum(f(a+i3*h)) + f(b))
  
  f4dx <- Deriv(f, 'x', nderiv = 4)
  
  error <- n/80 * h^5 * fminbnd(f4dx, a, b, maximum = TRUE)[[2]]
  
  plot(f)
  plot(f4dx)
  c(s, error)
}

romberg <- function(f, a, b, n) {
  h <- (b-a)/n
  
  kArray <- seq(1, n, by = 1)
  hArray <- (b-a)/(2^(kArray-1))
  
  R <- matrix(data = NA, nrow = n, ncol = n, byrow = FALSE, dimnames = NULL)
  
  for(j in kArray) {
    if (j == 1) 
      for(k in kArray) 
        if (k == 1) {
          R[1, 1] <- trapezoid(f, a, b)
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
  
  R
}

f = function(x) exp(x^2)


