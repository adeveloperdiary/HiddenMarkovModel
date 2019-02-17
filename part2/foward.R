

data = read.csv("data_r.csv")

a = matrix(c(0.54, 0.49, 0.46, 0.51),nrow = 2,ncol = 2)
b = matrix(c(0.16, 0.25, 0.26, 0.28, 0.58, 0.47),nrow = 2,ncol = 3)
initial_distribution = c(1/2, 1/2)

forward = function(v, a, b, initial_distribution){
  
  T = length(v)
  m = nrow(a)
  alpha = matrix(0, T, m)
  
  alpha[1, ] = initial_distribution*b[, v[1]]
  
  for(t in 2:T){
    tmp = alpha[t-1, ] %*% a
    alpha[t, ] = tmp * b[, v[t]]
  }
  return(alpha)
}

forward(data$Visible,a,b,initial_distribution)


