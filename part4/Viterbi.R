forward = function(v, a, b, initial_distribution){
  
  T = length(v)
  M = nrow(a)
  alpha = matrix(0, T, M)
  
  alpha[1, ] = initial_distribution*b[, v[1]]
  
  for(t in 2:T){
    tmp = alpha[t-1, ] %*% a
    alpha[t, ] = tmp * b[, v[t]]
  }
  return(alpha)
}

backward = function(v, a, b){
  T = length(v)
  M = nrow(a)
  beta = matrix(1, T, M)
  
  for(t in (T-1):1){
    tmp = as.matrix(beta[t+1, ] * b[, v[t+1]])
    beta[t, ] = t(a %*% tmp)
  }
  return(beta)
}


BaumWelch = function(v, a, b, initial_distribution, n.iter = 100){

  for(i in 1:n.iter){
    T = length(v)
    M = nrow(a)
    K=ncol(b)
    alpha = forward(v, a, b, initial_distribution)
    beta = backward(v, a, b)
    xi = array(0, dim=c(M, M, T-1))
    
    for(t in 1:T-1){
      denominator = ((alpha[t,] %*% a) * b[,v[t+1]]) %*% matrix(beta[t+1,]) 
      for(s in 1:M){
        numerator = alpha[t,s] * a[s,] * b[,v[t+1]] * beta[t+1,]
        xi[s,,t]=numerator/as.vector(denominator)
      }
    }
    
    
    xi.all.t = rowSums(xi, dims = 2)
    a = xi.all.t/rowSums(xi.all.t)
    
    gamma = apply(xi, c(1, 3), sum)  
    gamma = cbind(gamma, colSums(xi[, , T-1]))
    for(l in 1:K){
      b[, l] = rowSums(gamma[, which(v==l)])
    }
    b = b/rowSums(b)
    
  }
  return(list(a = a, b = b, initial_distribution = initial_distribution))
}


Viterbi=function(v,a,b,initial_distribution) {
  
  T = length(v)
  M = nrow(a)
  prev = matrix(0, T-1, M)
  omega = matrix(0, M, T)
  
  omega[, 1] = log(initial_distribution * b[, v[1]])
  for(t in 2:T){
    for(s in 1:M) {
      probs = omega[, t - 1] + log(a[, s]) + log(b[s, v[t]])
      prev[t - 1, s] = which.max(probs)
      omega[s, t] = max(probs)
    }
  }
  
  S = rep(0, T)
  last_state=which.max(omega[,ncol(omega)])
  S[1]=last_state
  
  j=2
  for(i in (T-1):1){
    S[j]=prev[i,last_state] 
    last_state=prev[i,last_state] 
    j=j+1
  }
  
  S[which(S==1)]='A'
  S[which(S==2)]='B'
  
  S=rev(S)
  
  return(S)
  
}

data = read.csv("data_r.csv")

M=2; K=3
A = matrix(1, M, M)
A = A/rowSums(A)
B = matrix(1:6, M, K)
B = B/rowSums(B)
initial_distribution = c(1/2, 1/2)

myout = BaumWelch(data$Visible, A, B, initial_distribution, n.iter = 100)
myout.hidden=Viterbi(data$Visible,myout$a,myout$b,initial_distribution)

library(HMM)
hmm =initHMM(c("A", "B"), c(1, 2, 3), 
              startProbs = initial_distribution,
              transProbs = A, emissionProbs = B)

true.out = baumWelch(hmm, data$Visible, maxIterations=100, pseudoCount=0)

true.viterbi = viterbi(true.out$hmm, data$Visible)
sum(true.viterbi != myout.hidden)

