lagrange.poly <- function(x, y) {
  
  l <- list()
  k <- 1
  
  for (i in x) {
    num <- 1
    denom <- 1
    
    p <- x[! x %in% i]
    
    for (j in p) {
      num <- paste(num, "*(x-", as.character(j), ")")
      denom <- paste(denom, "*(", as.character(i), "-", as.character(j), ")")
    }
    
    l[k] <- paste("(", num, ")", "/", "(", denom, ")")
    k <- k + 1
  }
  
  eq <- 0
  
  for (i in 1:length(y)) {
    eq <- paste(eq, '+', as.character(y[i]), "*", l[[i]])
  }
  
  parse(text=eq)
}