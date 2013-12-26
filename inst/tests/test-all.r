alpha_ <- 0.05
beta_ <- 0.2
St.n <- 1-0.25 # krankheitsfreie Überlebensrate zum Zeitpunkt t in der Novumgruppe
St.c <- 1-0.3 # krankheitsfreie Überlebensrate zum Zeitpunkt t in der Kontrollgruppe
n.c <- 1
n.e <- 1
N = 100
eta.e = 0.025
eta.c = 0.025
lambda.e=0.2
lambda.c=0.3
R=3
g=-0.5
S = 8
t_ = 5

test_that("correct log is used",
          expect_equal(hazard.ratio(8,3),log(8)/log(3)))

test_that("schoenfeld works",{
          expect_equal(schoenfeld(alpha_,beta_,St.n,St.c,n.c,n.e),2470.494,tolerance=0.002)
          })

test_that("lachin works",{
          expect_equal(lachin(n.c,n.e,lambda.c,lambda.e,eta.c,eta.e,g,alpha_,beta_,S,t_,R),130.9775,tolerance=0.002)
          expect_equal(lachin(n.c,n.e,lambda.c,lambda.e,eta.e,eta.e,g,alpha_,beta_,S,t_,R),130.9775,tolerance=0.002)
          expect_equal(lachin(n.c,n.e,lambda.c,lambda.e,0,0,g,alpha_,beta_,S,t_,R),125.1571,tolerance=0.002)
          })


