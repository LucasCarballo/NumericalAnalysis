require('orthopolynom')

weights <- function(p, x, n) {
  (2*(1-x^2))/(n^2*p(x)^2)
}

legendreCoeficients <- function(n) {
  legendrePolynomials <- legendre.polynomials(n)
  
  zeros <- solve(legendrePolynomials[[n+1]])
  
  P <- as.function(legendrePolynomials[[n]])
  weights <- weights(P, zeros, n)
  
  matrix(c(zeros, weights), nrow = n, ncol = 2, byrow = FALSE)
}

gaussianQuadrature <- function(f, a, b, n) {
  lCoeficients <- legendreCoeficients(n)

  f2ndx <- Deriv(f, 'x', nderiv = (2*n))
  
  error <- (((b-a)^(2*n+1) * fact(n)^4)/((2*n+1)*((fact(2*n)^3)))) * fminbnd(f2ndx, a, b, maximum = TRUE)[[2]]  
  
  g <- (b-a)/2 * sum(lCoeficients[1:n, 2] * f(((b-a) * lCoeficients[1:n, 1] + b + a)/2))
  
  c(g, error)
}
