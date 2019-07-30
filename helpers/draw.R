draw.trapezoid <- function(f, a, b, n = 1) {
  h <- (b-a) / n
  x <- seq(a, b, by=h)
  y <- f(x)
  
  plot(1, type="n", xlab="x", ylab="y", xlim=c(a-0.5, b+0.5), ylim=c(0, trunc(f(b))+1))
  
  for(i in seq(1, n, 1)) {
    polygon(c(x[i:(i+1)], rev(x[i:(i+1)])), c(0, 0, rev(y[i:(i+1)])), col = '#d0fffe', border = 'blue')  
  }
  
  curve(f, a, b, n = 101, add = TRUE, type = "l", col = "green") 
}
  
draw.simpson <- function(f, a, b, n = 1) {
  h <- (b-a)/n
  
  x = seq(a, b, by = h)
  xPolygon = seq(a, b, by = h/10)
  
  y = f(x)
  yPolygon = f(xPolygon)
  
  lExpression <- lagrange.poly(x, y)[[1]]
  
  L <- function(x) eval(lExpression)
  
  Ly <- L(x)
  LyPolygon <- L(xPolygon)
  
  plot(1, type = "n", xlab = "x", ylab = "y", xlim = c(a-0.5, b+0.5), ylim = c(0, trunc(f(b))+1))  
  
  polygon(c(min(x), xPolygon, max(x)), c(0, LyPolygon, 0), col = "#1b98e0") 
  
  for(i in seq(1, n, 1)) {
    segments(x[i], 0, x[i], Ly[i])
  }
  
  curve(f, a, b, n = 101, add = TRUE, type = "l", col = "green", lwd = "2")
  curve(L, a, b, n = 101, add = TRUE, type = "l", col = "red", lwd = "2")
}

draw.integrate <- function(f, a, b) {
  n <- 101
  x <- seq(a, b, by = (b-a)/n)
  
  plot(1, type = "n", xlab = "x", ylab = "y", xlim = c(a-0.5, b+0.5), ylim = c(0, trunc(f(b))+1))  
  
  polygon(c(min(x), x, max(x)), c(0, f(x), 0), col = "#1b98e0")
  
  curve(f, a, b, n = 101, add = TRUE, type = "l", col = "red", lwd = "2")
}
