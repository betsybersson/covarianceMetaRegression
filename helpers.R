##########################################################
## Reflect value to positive side of A
##########################################################
reflect = function(X,A){
  if (X<A) {
    out = A + (A-X)
  } else {
    out = X 
  }
  return(out)
}
##########################################################
## Squared Error Loss
# A: Estimator
# B: Truth
##########################################################
loss_sq = function(A,B){
  mat_trace(crossprod(A-B))
}
##########################################################
## Identity
# N: length of diagonal
##########################################################
eye = function(N){
  diag(rep(1,N))
}
##########################################################
## rwish
## copied directly from hoff amen
# input S0: centrality parameter
# input nu: degrees of freedom parameter
#
# output one random sample from wishart with mean S0*nu
##########################################################
rwish <-function(S0,nu=dim(S0)[1]+2){
  # sample from a Wishart distribution 
  # with expected value nu*S0 
  sS0<-(chol(S0))
  Z<-matrix(rnorm(nu*dim(S0)[1]),nu,dim(S0)[1])%*% sS0
  t(Z)%*%Z
}
##########################################################
## Random Multivariate or Matrix-variate Normal Sample
# M: nxp mean of normal density
# U: nxn row covariance matrix
# V: pxp column covariance matrix
##########################################################
rmatnorm = function(M,U = eye(nrow(M)),V = eye(ncol(M))) {
  # function for density of Y_{ n\times p} \sim N(M, V kron U)
  N = nrow(M)
  P = ncol(M)
  
  E = matrix(rnorm(N*P),ncol=P)
  
  out = M + t(chol(U)) %*% E %*% chol(V)
  
  return(out)
}
rmvnorm = function(M,V) {
  # function for density of Y_{px1} \sim N_p(M, V) (ind rows)
  P = length(M)
  
  E = rnorm(P)
  
  out = c(M) + t(chol(V)) %*% E  
  
  return(c(out))
}
##########################################################
## Trace of matrix
# X: square matrix
##########################################################
mat_trace = function(X){
  sum(diag(X))
}

##########################################################
### Un-vectorize matrix-variate data 
# input: x: n x p1p2 matrix
# for- p1: row dimension of matrix-variate data
#      p2: column dimension of matrix-variate data
#      n: number of samples
# output: X: p1 x p2 x n array 
##########################################################
vec.inv.array = function(x,p1,p2){
  
  N = nrow(x)
  XX = array(NA,dim=c(p1,p2,N))
  
  for ( nn in 1:N ){
    XX[,,nn] = matrix(x[nn,],ncol=p2,nrow=p1)
  }
  
  return(XX)
  
}
###########################
## rbind.list - rbind a list of vectors into a single matrix 
###########################
rbind.list = function(L){
  n.col = length(L[[1]])
  n.row = length(L)
  matrix(unlist(L),ncol = n.col, nrow = n.row,
         byrow = T)
}