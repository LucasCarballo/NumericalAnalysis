euler <- function(f, a, b, n, y0) {
  h <- (b-a) / n
  
  y <- rep(NA, n)
  y[1] <- y0
  
  x <- seq(a, b, by = h)
  
  for(i in 1:n) {
    y[i+1] <- y[i] + h * f(x[i], y[i])
  }
  
  print(x)
  print(y)
  
  plot(x, y, type="b")
  #curve(f, a, b, n = 101, add = TRUE, type = "l", col = "green") 
  
  result <- matrix(c(x, y), ncol = 2, byrow = FALSE)
  colnames(result) <- c('x', 'y')
  
  result
}

func <- function(t, y) 0.7 * y - t^2 + 1
