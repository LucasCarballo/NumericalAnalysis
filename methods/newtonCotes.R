require(Deriv)
require(pracma)

source('C:/projects/NumericalAnalysis/helpers/draw.R')
source('C:/projects/NumericalAnalysis/methods/lagrange.R')

trapezoid <- function(f, a, b, n, plot = TRUE) {
  print("Aplicando Regla del Trapecio")
  
  h <- (b-a)/n
  
  if((b-a) == n) {
    sumfxi <- 0
  } else {
    xi <- seq(a+h, b-h, by = h)  
    sumfxi <- sum(f(xi))
  }
  
  t <- h/2 * (f(a) + f(b) + 2*sumfxi)
  
  f2dx <- Deriv(f, 'x', nderiv = 2)
  
  error <- (b-a) * (h^2)/12 * fminbnd(f2dx, a, b, maximum = TRUE)[[2]]
  
  if (plot) draw.trapezoid(f, a, b, n)

  c(t, error)
}

simpson13 <- function(f, a, b, n) {
  if (n%%2 != 0) stop("En la regla de Simpson, n tiene que ser par")
  
  print("Aplicando Regla de Simpson de 1/3")
  
  h <- (b-a)/n
  
  i1 <- seq(1, n-1, by = 2) # impares
  
  if (n == 2) {
    sumfi2 <- 0
  } else {
    i2 <- seq(2, n-2, by = 2) # pares  
    sumfi2 <- sum( f(a+i2*h) )
  }
  
  s <- h/3 * ( f(a) + f(b) + 4*sum( f(a+i1*h) ) + 2*sumfi2 ) 
  
  draw.simpson(f, a, b, n)
  
  f4dx <- Deriv(f, 'x', nderiv = 4)
  
  error <- (b-a) * (h^4)/180 * fminbnd(f4dx, a, b, maximum = TRUE)[[2]]
  
  c(s, error)
}

simpson38 <- function(f, a, b, n) {
  if (n%%3 != 0) stop("En la regla de Simpson de 3/8 Compuesta, n tiene ser multiplo de 3")
  
  print("Aplicando Regla de Simpson de 3/8")
  
  h <- (b-a)/n
  
  i1 <- seq(1, n-2, by = 3)
  i2 <- seq(2, n-1, by = 3)
  
  if (n == 3) {
    sumfi3 <- 0
  } else {
    i3 <- seq(3, n-3, by = 3)
    sumfi3 <- sum( f(a+i3*h) )
  }
  
  s <- 3*h/8 * (f(a) + 3 * sum(f(a+i1*h)) + 3 * sum(f(a+i2*h)) + 2 * sumfi3 + f(b))
  
  f4dx <- Deriv(f, 'x', nderiv = 4)
  
  error <- n/80 * h^5 * fminbnd(f4dx, a, b, maximum = TRUE)[[2]]
  
  draw.simpson(f, a, b, n)
  
  c(s, error)
}

simpson <- function(f, a, b, n) {
  if (n%%2 == 0) {
    print(simpson13(f, a, b, n))
  }
  
  if (n%%3 == 0) {
    print(simpson38(f, a, b, n))
  } else if (n%%2 != 0 & n%%3 != 0) {
    stop("Por favor, utilice 'n' par para aplicar el metodo de Simpson de 1/3 o bien 'n' multiplo de 3 para aplicar el metodo de Simpson de 3/8")
  }
}
