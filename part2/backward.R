data = read.csv("data_r.csv")

a = matrix(c(0.54, 0.49, 0.46, 0.51),nrow = 2,ncol = 2)
b = matrix(c(0.16, 0.25, 0.26, 0.28, 0.58, 0.47),nrow = 2,ncol = 3)

backward = function(V, A, B){
  T = length(V)
  m = nrow(A)
  beta = matrix(1, T, m)
  
  for(t in (T-1):1){
    tmp = as.matrix(beta[t+1, ] * B[, V[t+1]])
    beta[t, ] = t(A %*% tmp)
  }
  return(beta)
}

backward(data$Visible,a,b)
