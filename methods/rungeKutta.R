rungeKutta <- function(f, a, b, y0, n) {
  h <- (b-a)/n
  
  t <- seq(a, b, by = h)
  y <- rep(NA, times=(n+1))
  
  y[1] = y0
  
  for (i in 2:(n+1)) {
    k1 <- h/2 * f(t[i-1], y[i-1])
    k2 <- h/2 * f(t[i-1] + h/2, y[i-1] + k1)
    k3 <- h/2 * f(t[i-1] + h/2, y[i-1] + k2)
    k4 <- h/2 * f(t[i-1] + h, y[i-1] + 2 * k3)
    
    y[i] <- y[i-1] + 1/3 * (k1 + 2 * k2 + 2 * k3 + k4)
  }
  
  plot(t, y, type="b")
  
  dat = cbind(t,y)
  dat
}