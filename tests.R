xx <- cbind(g=gl(2, 5), data.frame(matrix(rnorm(50), nrow=10)))
tt <- with(xx, split(xx, g))
tt
d1 <- tt[[1]][,-1]
d2 <- tt[[2]][,-1]
aa <- hotelling(d1, d2, type="nonexact")
aa1 <- hotelling(d1, d2, type="BS")
aa <- matrix(1:30, ncol=3)
idx <- combn(3:10, 2)


Y <- data.table(Y)
idx
a

theta1 <- function(Y) {
  foreach (i=1:NCOL(idx)) %dopar% {# angle of every two vectors
      thTmp <- Y[idx[, i], ]
      y1 <- thTmp[1, ]
      y2 <- thTmp[2, ]
      acos(sum(y1 * y2) / (sqrt(sum(y1 * y1) * sum(y2 * y2))))
    }
}
theta2 <- function(Y) {
  foreach (i=1:NCOL(idx)) %dopar% {# angle of every two vectors
      thTmp <- Y[idx[, i], ]
      y1 <- thTmp[1, ]
      y2 <- thTmp[2, ]
      Q <- Qs[idx[, i]]
      acos(sum(y1 * y2) / (sqrt(Q[1] * Q[2])))
    }
}

microbenchmark(theta1, theta2, times=500)
rr <- 11:12
microbenchmark(sqrt(rr[1])*sqrt(rr[2]), sqrt(rr[1]*rr[2]), times=500)

# Minhajuddin, A.T.M, Harris, I.R and Schucany,W.R.(2004) Simulating multivariate distributions with specific correlations. Journal of statistical computation and simulation. Vol. 74, No.8, Aug 2004, pp. 599-607
rmvgamma <- function(n, p, shape, rate=1/scale, rho=0, scale=1) {
  stopifnot(length(rho)==1,length(shape)==1,length(rate)==1,
            length(scale)==1,length(n)==1,length(p)==1,
            rate>0, shape>0, rho>=0, rho<=1, scale>0, n>=0, p>0)
  n <- round(n)
  p <- ceiling(p)
  theta <- rate*rho/(1-rho)
  k <- rnbinom(n, shape, rate/(rate+theta))
  matrix(rgamma(p*n, shape+k, rate+theta), n)
}


x <- rmvgamma(100000, 2, 1.9, 4, .5)
str(x)
hist(x,pr=T,br=100)
curve(dgamma(x,1.9,4), 0, 5, add=T)
ks.test(x,pgamma,shape=1.9,rate=4)
cor(x)
qqplot(x[,1],x[,2]);abline(0,1)

cor(matrix(rgamma(30000, 2), ncol=3))
