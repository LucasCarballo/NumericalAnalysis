library('Deriv')
library('pracma')
library('orthopolynom')

source('C:/projects/NumericalAnalysis/methods/newtonCotes.R')
source('C:/projects/NumericalAnalysis/methods/romberg.R')
source('C:/projects/NumericalAnalysis/methods/gaussianQuadrature.R')

source('C:/projects/NumericalAnalysis/methods/euler.R')
source('C:/projects/NumericalAnalysis/methods/taylor.R')
source('C:/projects/NumericalAnalysis/methods/rungeKutta.R')

integrate <- function(f, a, b, n, method) {
  switch(method,
         trapezoid = trapezoid(f, a, b, n),
         simpson = simpson(f, a, b, n),
         romberg = romberg(f, a, b, n),
         gauss = gaussianQuadrature(f, a, b, n))
}

solve <- function(f, a, b, n, y0, m, method) {
  switch(method,
         euler = euler(f, a, b, n, y0),
         taylor = taylor(f, a, b, n, y0, m),
         rungeKutta = rungeKutta(f, a, b, y0, n))
}


f <- function(x) exp(x^2)

g <- function(x) sin(x^2)

func <- function(t, y) 0.7 * y - t^2 + 1# > old.par <- par(mfrow=c(1, 2))
# > plot(faithful, main="Faithful eruptions")
# > plot(large.islands, main="Islands", ylab="Area")
# > par(old.par)