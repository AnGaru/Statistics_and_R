## make up some ‘true’ data
carrier<-rep(c(0,1), c(100,200))
null.y<-rnorm(300)
alt.y<-rnorm(300, mean=carrier/2)

#mean: t-test
t.test(null.y~carrier, var.equal=TRUE)
t.test(alt.y~carrier, var.equal=TRUE)

#mean: permutation test
null.diff<-mean(null.y[carrier==1])-mean(null.y[carrier==0])
alt.diff<-mean(alt.y[carrier==1])-mean(alt.y[carrier==0])
one.test <- function(x,y) {
  xstar<-sample(x)
  mean(y[xstar==1])-mean(y[xstar==0])
}
many.truenull <- replicate(1000, one.test(carrier, null.y))
many.falsenull <- replicate(1000, one.test(carrier, alt.y))
hist(many.truenull)
abline(v=null.diff, lwd=2, col="purple")
mean(abs(many.truenull) > abs(null.diff))
hist(many.falsenull)
abline(v=alt.diff, lwd=2, col="purple")
mean(abs(many.falsenull) > abs(alt.diff))
