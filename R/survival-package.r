hazard.ratio = function(n,c) log(n)/log(c)
hazard.diff = function(n,c) log(n)-log(c)
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


freedman = function(s.c,s.n,n.c,loss,alpha,beta,s){
  stopifnot(s %in% 1:2)
  r = (1-n.c)/n.c
  m = function(l) 1 - l
  d = prod((r*hazard.ratio(s.n,s.c) + 1)^2,
         (z(1-alpha/s) + z(1-beta))^2)/
    (r*(hazard.ratio(s.n,s.c)-1))

  a = d*(1+r)
  b = (r*m(s.n) + m(s.c)) * m(loss)
  a/b
}

rubinstein = function(eta.c,eta.e,alpha,beta,lambda.c,lambda.e,s,p.e,p.c,R,t){
  stopifnot(s %in% 1:2)
  e = function(p,l,eta){
    f = function(x) exp(-s*x)
    s = (l+eta)
    prod(l/s,
         1-(f(t-R)-f(t))/(s*R))
  }
  a = (z(alpha/s) + z(beta))/hazard.diff(lambda.c,lambda.e)
  prod(a^2,
       1/e(p.c,lambda.c,eta.c)+1/e(p.e,lambda.e,eta.e))
}
