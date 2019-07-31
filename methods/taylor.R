require(Deriv)

taylor <- function(f, a, b, n, y0, m) {
  h <- (b-a) / n
  
  y <- rep(NA, n)
  y[1] <- y0
  
  t <- seq(a, b, by = h)

  fdx <- list(as.expression(body(f)))
  
  for (i in 2:m) {
    fdt <- Deriv(fdx[[i-1]], "t")
    fdy <- Deriv(fdx[[i-1]], "y")
    
    fdx <- c(fdx, paste(as.expression(fdt), " + ((", as.expression(fdy), ") * (", as.expression(body(f)), "))"))
  }
  
  for(i in 1:n) {
    fkdx <- rep(NA, m)
    seqH <- rep(NA, m)
    for (j in 1:m) {
      fdxFunction <- function(t, y) eval(parse(text=fdx[[j]]))
      
      fkdx[j] <- fdxFunction(t[i], y[i])
      
      seqH[j] <- (h^(j-1))/factorial(j)
    }
    
    y[i+1] <- y[i] + h * sum(seqH * fkdx)
  }
  
  plot(t, y, type="b", lty = 1, lwd = 2, col = "green", main = "Taylor")
  
  result <- matrix(c(t, y), ncol = 2, byrow = FALSE)
  colnames(result) <- c('t', 'y')
  
  result
}
