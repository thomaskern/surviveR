hazard.ratio = function(n,c) log(n)/log(c)
z = function(n) qnorm(1-n)
num.req.events = function(alpha,beta,r,theta){
  denominator = r / (1+r)^2 * log(theta)^2
  sum(z(alpha/2), z(beta))^2 / denominator
}
div = function(x) x^-1

var.lik = function(lambda,eta,gamma,R,t){
  a = lambda / (lambda + eta)
  b = prod(lambda,
           eta,
           exp(-(lambda+eta)*t),
           1-exp(lambda+eta+gamma)*R)
  c = prod(1-exp(-gamma*R),
           lambda+eta,
           lambda+eta-gamma)
  lambda^2*div(a + b/c)
}

schoenfeld = function(alpha, beta, st.n, st.c, n.1,n.2){
  r = n.2/n.1 
  d = num.req.events(alpha,beta,r,hazard.ratio(st.n,st.c))
  d*(1+r)/(r*(2-st.n-st.c))
}

lachin = function(n.c,n.e,lambda.c,lambda.e,eta.c,eta.e,gamma,alpha,beta,s,t,R){
  n = n.c+n.e
  q.c = n.c/n
  q.e = n.e/n
  lambda.mean = q.c*lambda.c + q.e*lambda.e

  total = function(l.1,l.2) sqrt(var.lik(l.1,eta.e,gamma,R,t)*div(q.e) +
                               var.lik(l.2,eta.c,gamma,R,t)*div(q.c))

  a = z(alpha) * total(lambda.mean,lambda.mean)
  b = z(beta) * total(lambda.e,lambda.c)
  c = (lambda.c-lambda.e)^2
  (a + b)/c
}
