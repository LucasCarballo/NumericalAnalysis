euler <- function(f, a, b, n, y0) {
  h <- (b-a) / n
  
  y <- rep(NA, n)
  y[1] <- y0
  
  t <- seq(a, b, by = h)
  
  for(i in 1:n) {
    y[i+1] <- y[i] + h * f(t[i], y[i])
  }
  
  plot(t, y, type="b", lty = 10, lwd = 5, col = "blue", main = "Euler")
  
  result <- matrix(c(t, y), ncol = 2, byrow = FALSE)
  colnames(result) <- c('t', 'y')
  
  result
}
