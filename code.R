
dat<-read.csv("femaleMiceWeights.csv")
cnames<-colnames(dat)
twt<-dat[12,2]
elf<-dat$Bodyweight[11]
nMice<-length(dat$Diet)
meanHF<-mean(dat[13:24,2])
set.seed(1)
nrandMouse<-dat$Bodyweight[sample(13:24,size=1)]

library(dplyr)
controls<-filter(dat,Diet=="chow")
controls<-select(controls,Bodyweight)
unlist(controls)

#pipeline
chowdudes<-filter(dat,Diet=="chow")%>%select(Bodyweight)%>%unlist
mean(chowdudes)

library(downloader)
url="https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/msleep_ggplot2.csv"
filename <- basename(url)
download(url,filename)
dat<-read.csv("msleep_ggplot2.csv")
datPrim<-filter(dat,order=="Primates")
nPrimates<-nrow(datPrim)
class(datPrim)
class(select(datPrim,sleep_total))
meanSLtot<-select(datPrim,sleep_total)%>%unlist%>%mean

#Exploratory Data Analysis#
library(UsingR)
library(rafalib)
data(father.son,package="UsingR") ##available from CRAN
x <- father.son$fheight
round(sample(x,20),1)
hist(x,breaks=seq(floor(min(x)),ceiling(max(x))),main="Height histogram",xlab="Height in inches")
xs<-seq(floor(min(x)),ceiling(max(x)),0.1)
#compute empirical cumulative distribution
plot(xs,ecdf(x)(xs),type="l",xlab="Height in inches",ylab="F(x)")
#quantile-quantile plots
me<-mean(x)
esty<-sd(x)
m70<-mean(x>70)
pm70<-pnorm(70,me,esty)
ps<-seq(0.01,0.99,0.01)
qs<-quantile(x,ps)
normalqs<-qnorm(ps,me,esty)
plot(normalqs,qs,xlab="Normal percentiles",ylab="Height percentiles")
abline(0,1) ##identity line
#alternatively:
qqnorm(x)
qqline(x)
#multiple quantile-quantile plots
load("skew.RData")
par(mfrow=c(3,3)) #multifigure plot, 3x3 grid of plots
for (i in 1:9) {
  qqnorm(dat[,i])
}
##par(mfrow=c(1,1)) #shows only a plot
#boxplots
library(UsingR)
library(rafalib)
dist(exec.pay)
qqnorm(exec.pay)
boxplot(exec.pay,ylab="10,000s of dollars",ylim=c(0,400))
#scatterplots
data(father.son,package="UsingR")
x=father.son$fheight
y=father.son$sheight
plot(x,y, xlab="Father's height in inches", 
     ylab="Son's height in inches", 
     main=paste("correlation =",signif(cor(x,y),2)))
#stratification (boxplot comparison)
groups <- split(y,round(x)) 
boxplot(groups)
print(mean(y[ round(x) == 72]))
#bivariate normal distrib (Y conditioned on X)
groups <- split(y,round(x)) 
mypar(2,2)
for(i in c(5,8,11,14)){
  qqnorm(groups[[i]],main=paste0("X=",names(groups)[i]," strata"),
         ylim=range(y),xlim=c(-2.5,2.5))
  qqline(groups[[i]])
}



head(InsectSprays)
median(exec.pay)

set.seed(1)
dat <- read.csv("femaleMiceWeights.csv")
View(dat)
library(dplyr)
control <- filter(dat,Diet=="chow") %>% select(Bodyweight) %>% unlist
treatment <- filter(dat,Diet=="hf") %>% select(Bodyweight) %>% unlist
print( mean(treatment) )
print( mean(control) )
obsdiff <- mean(treatment) - mean(control)
print(obsdiff)

#means of random samples, they are random variables
library(downloader)
dir <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/"
filename <- "femaleControlsPopulation.csv"
url <- paste0(dir, filename)
##check if file exists and if it does not, download it:
if (!file.exists(filename)) download(url,destfile=filename)
##we don't usually have access to the data of the entire population
population <- read.csv("femaleControlsPopulation.csv")%>%unlist
control <- sample(population,12)
mean(control)
control <- sample(population,12)
mean(control)
control <- sample(population,12)
mean(control)
##each time, it returns a different, albeit similar, mean

#null hypothesis - null distribution
##12 control mice
control <- sample(population,12)
##another 12 control mice that we act as if they were not
treatment <- sample(population,12)
print(mean(treatment) - mean(control))
##let's repeat this process 10,000 times
n <- 10000
null <- vector("numeric",n)
for (i in 1:n) {
  control <- sample(population,12)
  treatment <- sample(population,12)
  null[i] <- mean(treatment) - mean(control)
}
##how many are larger than obsdiff? (percentage)
mean(abs(null) >= obsdiff)
##very few

#distributions
data(father.son,package="UsingR")
x <- father.son$fheight
round(sample(x,10),1)
#empirical cumulative distribution function
smallest <- floor( min(x) )
largest <- ceiling( max(x) )
values <- seq(smallest, largest,len=300)
heightecdf <- ecdf(x)
plot(values, heightecdf(values), type="l",
     xlab="a (Height in inches)",ylab="Pr(x <= a)")

#histograms
hist(x)
bins <- seq(smallest, largest)
hist(x,breaks=bins,xlab="Height (in inches)",main="Adult men heights")

#probability distributions
##Montecarlo simulation
n <- 100
library(rafalib)
nullplot(-5,5,1,30, xlab="Observed differences (grams)", ylab="Frequency")
totals <- vector("numeric",11)
for (i in 1:n) {
  control <- sample(population,12)
  treatment <- sample(population,12)
  nulldiff <- mean(treatment) - mean(control)
  j <- pmax(pmin(round(nulldiff)+6,11),1)
  totals[j] <- totals[j]+1
  text(j-6,totals[j],pch=15,round(nulldiff,1))
  ##if(i < 15) Sys.sleep(1) ##You can add this line to see values appear slowly
}
hist(null, freq=TRUE)
abline(v=obsdiff, col="red", lwd=2)

#normal distribution
1 - pnorm(obsdiff,mean(null),sd(null)) #normal approx.
#actual formula not shown

#ex1
library(downloader) 
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/femaleControlsPopulation.csv"
filename <- basename(url)
download(url, destfile=filename)
x <- unlist( read.csv(filename) )
mean(x)
set.seed(1)
abs(mean(sample(x,5))-mean(x))
set.seed(5)
abs(mean(sample(x,5))-mean(x))

#ex2
library(downloader) 
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/femaleControlsPopulation.csv"
filename <- basename(url)
download(url, destfile=filename)
x <- unlist( read.csv(filename) )
set.seed(1)
n<-1000
avgs<-vector("numeric",n)
for (i in 1:n){
  avgs[i]<-mean(sample(x,5))
}
av<-mean(x)
mean(abs(avgs-av)>1)
#repeat the above, with n<-10000
#repeat the above, with sample(x,50)

#ex3
install.packages("gapminder")
library(gapminder)
data(gapminder)
head(gapminder)
x<-subset(gapminder,year=="1952")$lifeExp
hist(x)
mean(x<=40)
mean(x<=60)-mean(x<=40)
qs=seq(from=min(x), to=max(x), length=20)
prop=function(q){mean(x<=q)}
props=sapply(qs,prop)
plot(qs,props)
#alternatively:
plot(qs, sapply(qs, function(q){mean(x<=q)}))

#ex Normal Distrib
x<-unlist(read.csv("femaleControlsPopulation.csv"))
# make averages5
set.seed(1)
n <- 1000
averages5 <- vector("numeric",n)
for(i in 1:n){
  X <- sample(x,5)
  averages5[i] <- mean(X)
}
# make averages50
set.seed(1)
n <- 1000
averages50 <- vector("numeric",n)
for(i in 1:n){
  X <- sample(x,50)
  averages50[i] <- mean(X)
}
#histograms
hist(average5)
hist(average50)
mean(averages50<=25)-mean(averages50<=23)
#repeat assuming that it's normally distributed (use pnorm)
pnorm(25,23.9,0.43)-pnorm(23,23.9,0.43)



#Populations and Samples#
dat<-read.csv("mice_pheno.csv")
library(dplyr)
controlPopulation <- filter(dat,Sex == "F" & Diet == "chow") %>% 
  select(Bodyweight) %>% unlist
length(controlPopulation)
hfPopulation <- filter(dat,Sex == "F" & Diet == "hf") %>%  
  select(Bodyweight) %>% unlist
length(hfPopulation)
#sample estimates
X<-sample(controlPopulation,25)
mean(X)
Y<-sample(hfPopulation,25)
mean(Y)
#t-test
library(dplyr)
dat<-read.csv("femaleMiceWeights.csv")
controls<-filter(dat,Diet=="chow")%>%select(Bodyweight)%>%unlist
treatment<-filter(dat,Diet=="hf")%>%select(Bodyweight)%>%unlist
N<-length(treatment) #in this case, elngth(controls)=length(treatment)
obs<-mean(treatment)-mean(controls)
se<-sqrt(var(treatment)/N + var(controls)/N)
tstat<-obs/se
2*(1-pnorm(tstat))
popul<-unlist(dat)
n<-10000
nulls<-vector("numeric",n)
for(i in 1:n){
  ctrl<-sample(popul,N)
  trtm<-sample(popul,N)
  serr<-sqrt(var(trtm)/N + var(ctrl)/N)
  nulls[i]<-(mean(trtm)-mean(ctrl))/serr
}
library(rafalib)
mypar()
qqnorm(nulls)
abline(0,1)
qqline(nulls)
ttest<-t.test(treatment,control) #does a t-test comparing a group to the control group


#ex1
x<-read.csv("mice_pheno.csv")
x<-na.omit(x)
library(dplyr)
#create a vector d with the bodyweights of all males on the chow diet
d<-subset(x,Diet=="chow")%>%subset(Sex=="M")
d<-d$Bodyweight
md<-mean(d)
library(rafalib)
popsd(d)
set.seed(1)
X<-sample(d,25)
mX<-mean(X)
e<-subset(x,Diet=="hf")%>%subset(Sex=="M")
e<-e$Bodyweight
me<-mean(e)
popsd(e)
set.seed(1)
Y<-sample(e,25)
mY<-mean(Y)
abs((me-md)-(mY-mX))
f<-subset(x,Diet=="chow")%>%subset(Sex=="F")
f<-f$Bodyweight
g<-subset(x,Diet=="hf")%>%subset(Sex=="F")
g<-g$Bodyweight
mf<-mean(f)
mg<-mean(g)
set.seed(2)
Z<-sample(f,25)
mZ<-mean(Z)
set.seed(2)
W<-sample(g,25)
mW<-mean(W)
abs((mg-mf)-(mW-mZ))

#ex2
pnorm(1,0,1)-pnorm(-1,0,1)
pnorm(2,0,1)-pnorm(-2,0,1)
pnorm(3,0,1)-pnorm(-3,0,1)
library(dplyr)
library(rafalib)
dat<-na.omit(read.csv("mice_pheno.csv"))
y<-filter(dat, Sex=="M" & Diet=="chow")%>%select(Bodyweight)%>%unlist
my<-mean(y)
vy<-popsd(y)
##non funzionano:
pnorm(vy,my,vy)-pnorm(-vy,my,vy)
pnorm(2*vy,my,vy)-pnorm(-2*vy,my,vy)
pnorm(3*vy,my,vy)-pnorm(-3*vy,my,vy)
##FUNZIONANO:
mean(y-my<=vy)-mean(y-my<=-vy)
mean(y-my<=2*vy)-mean(y-my<=-2*vy)
mean(y-my<=3*vy)-mean(y-my<=-3*vy)
z<-(y-my)/vy
head(z)
qqnorm(z)
abline(0,1)
mypar(2,2)
y <- filter(dat, Sex=="M" & Diet=="chow") %>% select(Bodyweight) %>% unlist
z <- ( y - mean(y) ) / popsd(y)
qqnorm(z);abline(0,1)
y <- filter(dat, Sex=="F" & Diet=="chow") %>% select(Bodyweight) %>% unlist
z <- ( y - mean(y) ) / popsd(y)
qqnorm(z);abline(0,1)
y <- filter(dat, Sex=="M" & Diet=="hf") %>% select(Bodyweight) %>% unlist
z <- ( y - mean(y) ) / popsd(y)
qqnorm(z);abline(0,1)
y <- filter(dat, Sex=="F" & Diet=="hf") %>% select(Bodyweight) %>% unlist
z <- ( y - mean(y) ) / popsd(y)
qqnorm(z);abline(0,1)
y <- filter(dat, Sex=="M" & Diet=="chow") %>% select(Bodyweight) %>% unlist
set.seed(1)
avgs <- replicate(10000, mean( sample(y, 25)))
mypar(1,2)
hist(avgs)
qqnorm(avgs)
qqline(avgs)
mean(avgs)
popsd(avgs)

#ex3
library(dplyr)
dat<-read.csv("femaleMiceWeights.csv")
n<-100
x<-sample(1:6, n, replace=TRUE)
p<-1/6
z<-(mean(x==6)-p)/(p+(1-p)/n)
m<-10000
set.seed(1)
rolls<-replicate(m, sample(1:6, n, replace=TRUE))#not quite what it asked
rolls<-replicate(m, sample())
zs<-vector("numeric",m)
for(i in 1:m){
  zs[i]<-(mean(rolls[,i]==6)-p)/sqrt(p*(1-p)/n)
}
mean(abs(zs)>2)#answ
qqnorm(zs)
abline(0,1)
#q2, actual answer
ps <- c(0.5,0.5,0.01,0.01)
ns <- c(5,30,30,100)
library(rafalib)
mypar(4,2)
for(i in 1:4){
  p <- ps[i]
  sides <- 1/p
  n <- ns[i]
  zs <- replicate(10000,{
    x <- sample(1:sides,n,replace=TRUE)
    (mean(x==1) - p) / sqrt(p*(1-p)/n)
  }) 
  hist(zs,nclass=7)
  qqnorm(zs)
  abline(0,1)
}
X <- filter(dat, Diet=="chow") %>% select(Bodyweight) %>% unlist
Y <- filter(dat, Diet=="hf") %>% select(Bodyweight) %>% unlist
mX<-mean(X)
#mean(Z~sqrt(12)*(mean(X)-avg(population))/(std(population)))=0
sadX<-sd(X)
2*(1-pnorm(2,mean=0,sd=sadX/sqrt(12)))#q7
mY<-mean(Y)
sadY<-sd(Y)
SE<-sqrt((var(X)+var(Y))/12)
tstat<-(mY-mX)/SE
2*(1-pnorm(2.055174,mean=0,sd=1)) #p-value of N(0,1)
t.test(Y,X)

#Inference#
#confidence intervals
chowPopulation<-read.csv("femaleControlsPopulation.csv")
chowPopulation<-unlist(chowPopulation)
mu_chow<-mean(chowPopulation)
print(mu_chow)
N<-30
chow<-sample(chowPopulation, N)
print(mean(chow))
se<-sd(chow)/sqrt(N)
print(se)
pnorm(2)-pnorm(-2)
(mean(chow)-mu_chow)/se
Q<-qnorm(1-0.05/2)
print(Q)
-Q<(mean(chow)-mu_chow)/se<Q #from this...
interval<-c(mean(chow)-Q*se, mean(chow)+Q*se)#to this.
interval[1]<mu_chow & interval[2]>mu_chow
library(rafalib)
B<-250
mypar()
plot(mu_chow+c(-7,7), c(1,1), type="n", xlab="weight",
     ylab="interval", ylim=c(1,B))
abline(v=mu_chow)
for (i in 1:B) {
  chow<-sample(chowPopulation, N)
  se<-sd(chow)/sqrt(N)
  interval<-c(mean(chow)-Q*se, mean(chow)+Q*se)
  covered<-
    mu_chow<=interval[2]&mu_chow>=interval[1]
  color<-ifelse(covered,1,2)
  lines(interval, c(i,i), col=color)
} #plots B confidence intervals, noting when
  # the mean of the entire population is not
  # inside a confInt
#repeat for N=5
#the Central Limit Theorem won't guarantee
# any convergence for N=5
#repeat for Q=qt(1-0.05/2, df=4)
Q<-qt(1-0.05/2, df=4)
#when you can't use the CLT, you use the t distribution

#rem.: %>% comes from dplyr

#statistical power
library(dplyr)
dat <- read.csv("mice_pheno.csv") #Previously downloaded 

controlPopulation <- filter(dat,Sex == "F" & Diet == "chow") %>%  
  select(Bodyweight) %>% unlist

hfPopulation <- filter(dat,Sex == "F" & Diet == "hf") %>%  
  select(Bodyweight) %>% unlist

mu_hf <- mean(hfPopulation)
mu_control <- mean(controlPopulation)
print(mu_hf - mu_control) #!=0, the null hyp.
# doesn't hold
print((mu_hf - mu_control)/mu_control * 100) # percent increase

set.seed(1)
N <- 5
hf <- sample(hfPopulation,N)
control <- sample(controlPopulation,N)
t.test(hf,control)$p.value #larger than alpha
#power is rejecting the null when the null is false
N <- 12
alpha <- 0.05
B<-2000
reject <- function(N, alpha=0.05){
  hf <- sample(hfPopulation,N) 
  control <- sample(controlPopulation,N)
  pval <- t.test(hf,control)$p.value
  pval < alpha
}
reject(12)
rejections <- replicate(B,reject(N))
# do note that replicate performs a Monte Carlo simulation
mean(rejections) #small power
Ns <- seq(5, 50, 5)
power<-sapply(Ns, function(N){
  rejections<-replicate(B,reject(N))
  mean(rejections)
})
plot(Ns,power, type="b")
#now, try for different alphas
N <- 30
alphas <- c(0.1,0.05,0.01,0.001,0.0001)
power <- sapply(alphas,function(alpha){
  rejections <- replicate(B,reject(N,alpha=alpha))
  mean(rejections)
})
plot(alphas, power, xlab="alpha", type="b", log="x")
calculatePvalue <- function(N) {
  hf <- sample(hfPopulation,N) 
  control <- sample(controlPopulation,N)
  t.test(hf,control)$p.value
}
Ns <- seq(10,200,by=10)
Ns_rep <- rep(Ns, each=10)
pvalues <- sapply(Ns_rep, calculatePvalue)
plot(Ns_rep, pvalues, log="y", xlab="sample size",
     ylab="p-values")
abline(h=c(.01, .05), col="red", lwd=2)
#a large sample size helps confirm that the null
# hypothesis doesn't hold, but the corresponding
# p-values will be very small
#a better statistic to report is the effect size
# with a confidence interval or some statistic which
# gives the reader a sense of the change in a
# meaningful scale, like reporting the effect size
# as a percentage by dividing the difference and
# the confidence interval by the control population
# mean
N <- 12
hf <- sample(hfPopulation, N)
control <- sample(controlPopulation, N)
diff <- mean(hf) - mean(control)
diff / mean(control) * 100
t.test(hf, control)$conf.int / mean(control) * 100
#Cohen's 'd', the difference between the groups
# divided by the pooled std of the two groups
sd_pool <- sqrt(((N-1)*var(hf) + (N-1)*var(control))/(2*N - 2))
diff / sd_pool

#ex1
library(downloader)
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/babies.txt"
filename <- basename(url)
download(url, destfile=filename)
babies <- read.table("babies.txt", header=TRUE)
bwt.nonsmoke <- filter(babies, smoke==0) %>% select(bwt) %>% unlist 
bwt.smoke <- filter(babies, smoke==1) %>% select(bwt) %>% unlist
library(rafalib)
mean(bwt.nonsmoke)-mean(bwt.smoke)
popsd(bwt.nonsmoke)
popsd(bwt.smoke)
set.seed(1)
dat.ns<-sample(bwt.nonsmoke, 25)
dat.s<-sample(bwt.smoke, 25)
tval<-(mean(dat.ns)-mean(dat.s))/sqrt((var(dat.ns)+var(dat.s))/25)
pval <- 1-(pnorm(abs(tval))-pnorm(-abs(tval)))
#also:
pval==2*pnorm(-abs(tval)) #more or less

#ex2
babies <- read.table("babies.txt", header=TRUE)
bwt.nonsmoke <- filter(babies, smoke==0) %>% select(bwt) %>% unlist 
bwt.smoke <- filter(babies, smoke==1) %>% select(bwt) %>% unlist
library(rafalib)
mean(bwt.nonsmoke)-mean(bwt.smoke)
popsd(bwt.nonsmoke)
popsd(bwt.smoke)
#q1
N<-25
set.seed(1)
dat.ns<-sample(bwt.nonsmoke, N)
dat.s<-sample(bwt.smoke, N)
Q<-qt(1-0.01/2, df=2*N-2)
seNS<-mean(dat.ns)/sqrt(N)
seSM<-mean(dat.s)/sqrt(N)
intervalNS<-c(mean(dat.ns)-Q*seNS, mean(dat.ns)+Q*seNS)
intervalSM<-c(mean(dat.s)-Q*seSM, mean(dat.s)+Q*seSM)
#q3
N<-5
set.seed(1)
dat.ns<-sample(bwt.nonsmoke, N)
dat.s<-sample(bwt.smoke, N)
t.test(dat.ns, dat.s)

#ex3
babies <- read.table("babies.txt", header=TRUE)
bwt.nonsmoke <- filter(babies, smoke==0) %>% select(bwt) %>% unlist 
bwt.smoke <- filter(babies, smoke==1) %>% select(bwt) %>% unlist
library(rafalib)
mean(bwt.nonsmoke)-mean(bwt.smoke)
popsd(bwt.nonsmoke)
popsd(bwt.smoke)
#q1
set.seed(1)
N<-5
dat.s<-sample(bwt.smoke, N)
dat.ns<-sample(bwt.nonsmoke, N)
t.test(dat.s, dat.ns)
#q2
set.seed(1)
reject <- function(N, alpha=0.05){
  dat.s<-sample(bwt.smoke, N)
  dat.ns<-sample(bwt.nonsmoke, N)
  pval <- t.test(dat.s, dat.ns)$p.value
  pval < alpha
}
B<-10000
rejections <- replicate(B,reject(N))
mean(rejections)
#q3: try q2 for N in (30, 60, 90, 120)
# which value of N returns a power of about 0.8?

#q4: as q3, with alpha=0.01

#Inference II#
#Monte Carlo
library(dplyr)
dat <- read.csv("mice_pheno.csv")
controlPopulation <- filter(dat,Sex == "F" & Diet == "chow") %>%  
  select(Bodyweight) %>% unlist
ttestgenerator <- function(n) {
  #note that here we have a false "high fat" group where we actually
  #sample from the nonsmokers. this is because we are modeling the *null*
  cases <- sample(controlPopulation,n)
  controls <- sample(controlPopulation,n)
  tstat <- (mean(cases)-mean(controls)) / 
    sqrt( var(cases)/n + var(controls)/n ) 
  return(tstat)
}
ttests <- replicate(1000, ttestgenerator(10)) #MC sim.
hist(ttests)
qqnorm(ttests)
abline(0,1)
ttests <- replicate(1000, ttestgenerator(3)) #MC sim.
hist(ttests)
qqnorm(ttests) # a worse sim.
abline(0,1)
ps <- (seq(0,999)+0.5)/1000
qqplot(qt(ps,df=2*3-2),ttests,xlim=c(-6,6),ylim=c(-6,6))
abline(0,1)
qqnorm(controlPopulation)
qqline(controlPopulation)
#if you don't have access to the data of the entire
# population, you can simulate that data
# (in this case, we'll assume that the entire
#   population is normally distributed)
controls<- rnorm(5000, mean=24, sd=3.5) #parametric MC
ttestgenerator <- function(n, mean=24, sd=3.5) {
  cases <- rnorm(n,mean,sd)
  controls <- rnorm(n,mean,sd)
  tstat <- (mean(cases)-mean(controls)) / 
    sqrt( var(cases)/n + var(controls)/n ) 
  return(tstat)
}
ttests <- replicate(1000, ttestgenerator(3))
qqnorm(ttests) #how the simulated control pop. fares
               # compared to the normal distrib.
abline(0,1)
ps <- (seq(0,999)+0.5)/1000 #now, with t distrib.
                    # and 4 degrees of freedom
qqplot(qt(ps,df=2*3-2),ttests,xlim=c(-6,6),ylim=c(-6,6))
abline(0,1)

#ex1
#q1:
set.seed(1)
sasa<-rnorm(5)
msa<-mean(sasa)
sdsa<-sd(sasa)
t<-sqrt(5)*msa/sdsa
#q2
set.seed(1)
N<-30
B<-1000
sasas<-replicate(B, rnorm(N))
msas<-sapply(1:B, function(x) mean(sasas[,x]))
sdsas<-sapply(1:B, function(x) sd(sasas[,x]))
ts<-sqrt(N)*msas
ts<-ts/sdsas
mean(ts>2)
#q3
C=100; ps = seq(1/(C+1), 1-1/(C+1),len=C)
qqplot(qt(ps,df=2*N-2),ts,xlim=c(-6,6),ylim=c(-6,6))
abline(0,1)
Ns<-seq(5,30,5)#ripeti q2 con questi N
sasasas<-sapply(Ns, function(x) replicate(B, rnorm(x)))
#dice che ciascuna approx è molto vicina alla distrib.
# normale, ma l'approx. peggiora secondo il codice
# che ho usato
#codice della risposta:
library(rafalib)
mypar(3,2)
Ns<-seq(5,30,5)
B <- 1000
mypar(3,2)
LIM <- c(-4.5,4.5)
for(N in Ns){
  ts <- replicate(B, {
    X <- rnorm(N)
    #sqrt(N)*mean(X)/sd(X)
    median(X)
  })
  ps <- seq(1/(B+1),1-1/(B+1),len=B)
  #qnorm(ps, sd=1/sqrt(N))
  qqplot(qt(ps,df=N-1),ts,main=N,
         xlab="Theoretical",ylab="Observed",
         xlim=LIM, ylim=LIM)
  abline(0,1)
} 
#infatti, doveva essere df=N-1!
#q4:
# come q3, ma df=2*N-2
#q5:

#codice della risposta:
set.seed(1)
N <- 15
B <- 10000
tstats <- replicate(B,{
  X <- sample(c(-1,1), N, replace=TRUE)
  sqrt(N)*mean(X)/sd(X)
  #median(X)
})
ps=seq(1/(B+1), 1-1/(B+1), len=B) 
qqplot(qt(ps,N-1), tstats, xlim=range(tstats))
abline(0,1)
#The population data is not normal thus the theory does not apply.
#We check with a Monte Carlo simulation. The qqplot shows a large tail. 
#Note that there is a small but positive chance that all the X are the same.
##In this case the denominator is 0 and the t-statistics is not defined
#q6:
# come q5, ma N=1000
#With N=1000, CLT kicks in and the t-statistic is approximated with normal 0,1
##Furthermore, t-distribution with df=999 and normal are practically the same.
#q7:
# how does the median of a sample taken from normally
# distributed population with mean 0 and standard
# deviation 1 behave?
#codice della risposta:
set.seed(1)
Ns <- seq(5,45,5)
library(rafalib)
mypar(3,3)
for(N in Ns){
  medians <- replicate(10000, median ( rnorm(N) ) )
  title <- paste("N=",N,", avg=",round( mean(medians), 2) , ", sd*sqrt(N)=", round( sd(medians)*sqrt(N),2) )
  qqnorm(medians, main = title )
  qqline(medians)
}
##there is an asymptotic result that says SD is sqrt(N*4*dnorm(0)^2)

#permutation tests
mypar()
dat=read.csv("femaleMiceWeights.csv")

library(dplyr)

control <- filter(dat,Diet=="chow") %>% select(Bodyweight) %>% unlist
treatment <- filter(dat,Diet=="hf") %>% select(Bodyweight) %>% unlist
obsdiff <- mean(treatment)-mean(control)
N <- 12
#shuffle test AND control and sample from there
avgdiff <- replicate(1000, {
  all <- sample(c(control,treatment))
  newcontrols <- all[1:N]
  newtreatments <- all[(N+1):(2*N)]
  return(mean(newtreatments) - mean(newcontrols))
})
hist(avgdiff)
abline(v=obsdiff, col="red", lwd=2)
#the proportion of permutations with larger difference
(sum(abs(avgdiff) > abs(obsdiff)) + 1) / (length(avgdiff) + 1)

N <- 5
control <- sample(control,N)
treatment <- sample(treatment,N)
obsdiff <- mean(treatment)- mean(control)
avgdiff <- replicate(1000, {
  all <- sample(c(control,treatment))
  newcontrols <- all[1:N]
  newtreatments <- all[(N+1):(2*N)]
  return(mean(newtreatments) - mean(newcontrols))
})
hist(avgdiff)
abline(v=obsdiff, col="red", lwd=2)
#there is no theoretical guarantee that the null
# distribution estimated from permutations
# approximates the actual null distribution
#also, permutations tests still have assumptions:
# samples are assumed to be independent and
# “exchangeable”. If there is hidden structure
# in your data, then permutation tests can result
# in estimated null distributions that underestimate
# the size of tails because the permutations may
# destroy the existing structure in the original data.

#ex2
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/babies.txt"
filename <- basename(url)
download(url, destfile=filename)
babies <- read.table("babies.txt", header=TRUE)
bwt.nonsmoke <- filter(babies, smoke==0) %>% select(bwt) %>% unlist 
bwt.smoke <- filter(babies, smoke==1) %>% select(bwt) %>% unlist
#q1:
N<-10
#for only one sample
set.seed(1)
nonsmokers <- sample(bwt.nonsmoke , N)
smokers <- sample(bwt.smoke , N)
obs <- mean(smokers) - mean(nonsmokers)
dat <- c(smokers,nonsmokers)
#for a null distrib.
#shuffle <- sample( dat )
#smokerssing <- shuffle[1:N]
#nonsmokerssing <- shuffle[(N+1):(2*N)]
#obsdiff<-mean(smokerssing)-mean(nonsmokerssing)
set.seed(1)
avgdiff <- replicate(1000, {
  shuffle <- sample(dat)
  smokersstar <- shuffle[1:N]
  nonsmokersstar <- shuffle[(N+1):(2*N)]
  return(mean(smokersstar) - mean(nonsmokersstar))
})
permPVal<-(sum(abs(avgdiff) > abs(obs)) + 1) / (length(avgdiff) + 1)
permPVal
hist(avgdiff)
abline(v=obs, col="red", lwd=2)
#mdiff<-mean(avgdiff)
#sddiff<-sd(avgdiff)
#tstat<-sqrt(N)*obsdiff/sqrt(var(smokerssing)+var(nonsmokerssing))
#q2: con la mediana al posto della media
N<-10
set.seed(1)
#one sample
nonsmokers <- sample(bwt.nonsmoke , N)
smokers <- sample(bwt.smoke , N)
obs <- median(smokers) - median(nonsmokers)
#null distrib.
dat <- c(smokers,nonsmokers)
#shuffle <- sample( dat )
#smokerssing <- shuffle[1:N]
#nonsmokerssing <- shuffle[(N+1):(2*N)]
#obsdiff<-mean(smokerssing)-mean(nonsmokerssing)
set.seed(1)
avgdiff <- replicate(1000, {
  shuffle <- sample(dat)
  smokersstar <- shuffle[1:N]
  nonsmokersstar <- shuffle[(N+1):(2*N)]
  return(median(smokersstar) - median(nonsmokersstar))
})
permPVal<-(sum(abs(avgdiff) > abs(obs)) + 1) / (length(avgdiff) + 1)
permPVal
hist(avgdiff)
abline(v=obs, col="red", lwd=2)

#association tests
##Fisher's test
tab <- matrix(c(3,1,1,3),2,2)
rownames(tab)<-c("Poured Before","Poured After")
colnames(tab)<-c("Guessed before","Guessed after")
tab
fisher.test(tab,alternative="greater")
##chi square test
disease=factor(c(rep(0,180),rep(1,20),rep(0,40),rep(1,10)),
               labels=c("control","cases"))
genotype=factor(c(rep("AA/Aa",200),rep("aa",50)),
                levels=c("AA/Aa","aa"))
dat <- data.frame(disease, genotype)
dat <- dat[sample(nrow(dat)),] #shuffle them up
head(dat)
table(genotype)
table(disease)
tab <- table(genotype,disease)
tab
odds_ratio<-(tab[2,2]/tab[2,1]) / (tab[1,2]/tab[1,1])
odds_ratio
p_val=mean(disease=="cases")
p_val
###expected table for chi square test
expected <- rbind(c(1-p_val,p_val)*sum(genotype=="AA/Aa"),
                  c(1-p_val,p_val)*sum(genotype=="aa"))
dimnames(expected)<-dimnames(tab)
expected
chisq.test(tab)$p.value
###more data makes the chi square test more accurate
tab<-tab*10
chisq.test(tab)$p.value
###using the theory of generalized linear models to
### estimate a confidence interval for the odds ratio
fit <- glm(disease~genotype,family="binomial",data=dat)
coeftab<- summary(fit)$coef
coeftab
###The second row of the table shown above gives you
### the estimate and SE of the log odds ratio. Mathematical
### theory tells us the this estimate is approximately normally
### distributed. We can therefore form a confidence interval
###and then exponentiate to provide a confidence interval
### for the OR.
ci <- coeftab[2,1] + c(-2,2)*coeftab[2,2]
exp(ci)

#ex3
d = read.csv("assoctest.csv")
#q1:
tab<-table(d)
odds_ratio<-(tab[2,2]/tab[2,1]) / (tab[1,2]/tab[1,1])
odds_ratio
chisq.test(tab)$statistic
#q2:
fisher.test(tab)$p.value

#Exploratory Data Analysis#
#scatterplot
library(UsingR)
data("father.son")
x=father.son$fheight
y=father.son$sheight
plot(x,y,xlab="Father's height in inches",ylab="Son's height in inches",main=paste("correlation =",signif(cor(x,y),2)))
##stratification: what might be the height of a son chosen at random?
## and what if that son had a father 72 inches tall?
groups <- split(y,round(x)) 
boxplot(groups)
print(mean(y[ round(x) == 72]))
##conditioning and bivariate normal distrib.
library(rafalib)
mypar(2,2)
for(i in c(5,8,11,14)){
  qqnorm(groups[[i]],main=paste0("X=",names(groups)[i]," strata"),
         ylim=range(y),xlim=c(-2.5,2.5))
  qqline(groups[[i]])
} #qqline(groups[[i]]) in this case is called regression line
##plotting the avg.s of strata of son heights
x=(x-mean(x))/sd(x)
y=(y-mean(y))/sd(y)
means=tapply(y,round(x*4)/4,mean)
fatherheights=as.numeric(names(means))
mypar(1,1)
plot(fatherheights,means,ylab="average of strata of son heights",ylim=range(fatherheights))
abline(0,cor(x,y))
##why correlation by itself can deceive
a=rnorm(100);a[1]=10
b=rnorm(100);b[1]=11
plot(a,b,main=paste("correlation =",signif(cor(a,b),2)))

#ex1
data(nym.2002, package="UsingR")
#q1:
library(dplyr)
males<-filter(nym.2002, gender=="Male")
females<-filter(nym.2002, gender=="Female")
mage=males$age
mime=males$time
cor(mage, mime)
#q2:
fage=females$age
fime=females$time
cor(fage, fime)
#q3:
groups <- split(mime,round(mage)) 
boxplot(groups)
groups <- split(fime,round(fage)) 
boxplot(groups)
print(mean(y[ round(x) == 72]))
mime=(mime-mean(mime))/sd(mime)
mage=(mage-mean(mage))/sd(mage)
fime=(fime-mean(fime))/sd(fime)
fage=(fage-mean(fage))/sd(fage)
plot(mage, mime, ylab="average of strata of time of completion for males", ylim=range(mime))
plot(fage, fime, ylab="average of strata of time of completion for females", ylim=range(fime))
#codice della risposta:
library(rafalib)
mypar(2,2)
plot(females$age, females$time)
plot(males$age, males$time)
group <- floor(females$age/5) * 5 #a better approach
boxplot(females$time~group)
group <- floor(males$age/5) * 5
boxplot(males$time~group)

#symmetry of log ratios
x <- 2^(rnorm(100))
y <- 2^(rnorm(100)) 
ratios <- x / y 
mypar(1,2)
hist(ratios)
logratios <- log2(ratios)
hist(logratios)
#even if this example is too direct, most of the
# time, it's better to compute the ratios of the
# logarithms of the compared values rather than
# the values themselves

#ex2
data(nym.2002, package="UsingR")
time = sort(nym.2002$time)
#q1:
min(time)/median(time)
#q2:
max(time)/median(time)
#post
plot(time/median(time), ylim=c(1/4,4)) #linear plot
abline(h=c(1/2,1,2))
plot(log2(time/median(time)),ylim=c(-2,2)) #log plot
abline(h=-1:1) #these lines are equally spaced
# in the log plot

#plots to avoid
pie(browsers,main="Browser Usage (August 2013)")
#Pie charts are a very bad way of displaying information.
# The eye is good at judging linear measures and bad at
# judging relative areas. A bar chart or dot chart is a
# preferable way of displaying this type of data.
barplot(browsers, main="Browser Usage (August 2013)", ylim=c(0,55))
abline(h=1:5 * 10)
barplot(browsers, add=TRUE)
#A barplot is a more appropriate way to plot the
# previous chart, because you can compare more
# easily the percentages
#Do avoid 3D barplots, as they are more confusing
#Do not ever use donut charts, as the one way
# one could compare the quantities, the angles,
# are removed from this chart (compared to the
# pie chart)

#Barplots are useful for displaying percentages.
#Not so much to compare two groups.
library("downloader")
filename <- "fig1.RData"
url <- "https://github.com/kbroman/Talk_Graphs/raw/master/R/fig1.RData"
if (!file.exists(filename)) download(url,filename)
load(filename)
library(rafalib)
mypar()
dat <- list(Treatment=x,Control=y)
barplot(dat, xlab="Group",ylab="Response")
abline(h=1:5 * 10)
barplot(dat, add=TRUE)
#In this case it's more appropriate to create boxplots,
# as they show how the data is distributed inside
# each group
boxplot(dat,xlab="Group",ylab="Response",cex=0)
stripchart(dat,vertical=TRUE,method="jitter",pch=16,add=TRUE,col=1)
#This problem is compounded when the data has outliers
# or very large tails
library(downloader)
url <- "https://github.com/kbroman/Talk_Graphs/raw/master/R/fig3.RData"
filename <- "fig3.RData"
if (!file.exists(filename)) download(url, filename)
load(filename)
library(rafalib)
mypar()
dat <- list(Treatment=x,Control=y)
barplot(dat, xlab="Group",ylab="Response")
abline(h=1:5 * 10)
barplot(dat, add=TRUE)
#The boxplots will show the quantiles and the outliers
mypar(1,2)
boxplot(dat,xlab="Group",ylab="Response",cex=0)
stripchart(dat,vertical=TRUE,method="jitter",pch=16,add=TRUE,col=1)
boxplot(dat,xlab="Group",ylab="Response",log="y",cex=0)
stripchart(dat,vertical=TRUE,method="jitter",pch=16,add=TRUE,col=1)

#Just showing the regression line isn't enough
url <- "https://github.com/kbroman/Talk_Graphs/raw/master/R/fig4.RData"
filename <- "fig4.RData"
if (!file.exists(filename)) download(url, filename)
load(filename)
mypar(1,2)
plot(x,y,lwd=2,type="n")
fit <- lm(y~x)
abline(fit$coef,lwd=2)
b <- round(fit$coef,4)
text(78, 200, paste("y =", b[1], "+", b[2], "x"), adj=c(0,0.5))
rho <- round(cor(x,y),4) 
text(78, 187,expression(paste(rho," = 0.8567")),adj=c(0,0.5))
plot(x,y,lwd=2)
fit <- lm(y~x)
abline(fit$coef,lwd=2)
#The graph on the right is more complete, as it shows
# the data coords. the regr. is based on
#For large amounts of points, one should use
# the hexbin function from the package hexbin

#High correlation doesn't imply replication.
# For ex., a high corr. plot might reveal a lower
# corr. if scaled to a log-log plot.
#A way to counter this is to study the Minus-Average
# plot between two groups, comparing the average of
# the samples on the log scale and the difference
# between two samples on the log scale.

#Using paired barplots to compare before-after groups
# is also a bad idea. A better approach would be to
# make a scatterplot. An alternative one is to
# make a scatterplot of the before sample and
# the before-after change
set.seed(12201970)
before <- runif(6, 5, 8)
after <- rnorm(6, before*1.05, 2)
li <- range(c(before, after))
ymx <- max(abs(after-before))

mypar(1,2)
plot(before, after, xlab="Before", ylab="After",
     ylim=li, xlim=li)
abline(0,1, lty=2, col=1)

plot(before, after-before, xlab="Before", ylim=c(-ymx, ymx),
     ylab="Change (After - Before)", lwd=2)
abline(h=0, lty=2, col=1)
#You can also use line plots and boxplots
z <- rep(c(0,1), rep(6,2))
mypar(1,2)
plot(z, c(before, after),
     xaxt="n", ylab="Response",
     xlab="", xlim=c(-0.5, 1.5))
axis(side=1, at=c(0,1), c("Before","After"))
segments(rep(0,6), before, rep(1,6), after, col=1)     

boxplot(before,after,names=c("Before","After"),ylab="Response")

#Please, don't use 3D graphs
##First read data
library(downloader)
filename <- "fig8dat.csv"
url <- "https://github.com/kbroman/Talk_Graphs/raw/master/R/fig8dat.csv"
if (!file.exists(filename)) download(url, filename)
x <- read.table(filename, sep=",", header=TRUE)

##Now make alternative plot
mypar()
plot(x[,1],x[,2],xlab="log Dose",ylab="Proportion survived",ylim=c(0,1),
     type="l",lwd=2,col=1)
lines(x[,1],x[,3],lwd=2,col=2)
lines(x[,1],x[,4],lwd=2,col=3)
legend(1,0.4,c("Drug A","Drug B","Drug C"),lwd=2, col=1:3)
#Humans actually can't see very well in 3D

#Use less digits
heights <- cbind(rnorm(8,73,3),rnorm(8,73,3),rnorm(8,80,3),
                 rnorm(8,78,3),rnorm(8,78,3))
colnames(heights)<-c("SG","PG","C","PF","SF")
rownames(heights)<- paste("team",1:8)
heights
round(heights,1)

#Robust Summaries#
set.seed(1)
x=c(rnorm(100,0,1)) ##real distribution
x[23] <- 100 ##mistake made in 23th measurement
boxplot(x)
cat("The average is",mean(x),"and the SD is",sd(x))
#better metrics in this case would be the median
# and the median average deviation (MAD)
median(x)
mad(x)
#Spearman correlation (reprise)
set.seed(1)
x=c(rnorm(100,0,1)) ##real distribution
x[23] <- 100 ##mistake made in 23th measurement
y=c(rnorm(100,0,1)) ##real distribution
y[23] <- 84 ##similar mistake made in 23th measurement
library(rafalib)
mypar()
plot(x,y,main=paste0("correlation=",round(cor(x,y),3)),pch=21,bg=1,xlim=c(-3,100),ylim=c(-3,100))
abline(0,1)
#x and y should not be correlated at all, but they
# look correlated
#The Spearman correlation follows the general idea
# of median and MAD, that of using quantiles.
#The idea is simple: we convert each dataset to
# ranks and then compute correlation:
mypar(1,2)
plot(x,y,main=paste0("correlation=",round(cor(x,y),3)),pch=21,bg=1,xlim=c(-3,100),ylim=c(-3,100))
#spearmen corr.
plot(rank(x),rank(y),main=paste0("correlation=",round(cor(x,y,method="spearman"),3)),pch=21,bg=1,xlim=c(-3,100),ylim=c(-3,100))
abline(0,1)
#rem.: the Pearson correlation is just cor(x,y)

#ex1:
mypar()
data(ChickWeight)
head(ChickWeight)
plot( ChickWeight$Time, ChickWeight$weight, col=ChickWeight$Diet)
chick = reshape(ChickWeight, idvar=c("Chick","Diet"), timevar="Time",
                direction="wide")
head(chick)
chick = na.omit(chick)
#q1:
mean(chick$weight.4)/median(chick$weight.4)
c4<-c(chick$weight.4, 3000)
mean(c4)/mean(chick$weight.4)
#q2:
median(c4)/median(chick$weight.4)
#q3:
sd(c4)/sd(chick$weight.4)
#q4:
mad(c4)/mad(chick$weight.4)
#q5:
plot(chick$weight.4, chick$weight.21)
c21<-c(chick$weight.21,3000)
drugcorr<-cor(c4, c21)
actcorr<-cor(chick$weight.4, chick$weight.21)
drugcorr/actcorr

#Wilcoxon Rank Sum Test#
#Like the mean and the standard dev., the t-test is
# susceptible to outliers.
#The Wilcoxon test is an alternative to the t-test
# that is resilient to outliers.
set.seed(779) ##779 picked for illustration purposes
N=25
x<- rnorm(N,0,1)
y<- rnorm(N,0,1)
x[1] <- 5
x[2] <- 7
cat("t-test pval:",t.test(x,y)$p.value)
cat("Wilcox test pval:",wilcox.test(x,y)$p.value)
#The basic idea behind the z-score is to:
#   1) combine all the data
#   2) turn the values into ranks
#   3) separate them back into their groups, and
#   4) compute the sum or average rank and
#         perform a test.
library(rafalib)
mypar(1,2)

stripchart(list(x,y),vertical=TRUE,ylim=c(-7,7),ylab="Observations",pch=21,bg=1)
abline(h=0)

xrank<-rank(c(x,y))[seq(along=x)]
yrank<-rank(c(x,y))[-seq(along=y)]

stripchart(list(xrank,yrank),vertical=TRUE,ylab="Ranks",pch=21,bg=1,cex=1.25)

ws <- sapply(x,function(z) rank(c(z,y))[1]-1)
text( rep(1.05,length(ws)), xrank, ws, cex=0.8)
W<-sum(ws)
n1<-length(x);n2<-length(y)
Z <- (mean(ws)-n2/2)/ sqrt(n2*(n1+n2+1)/12/n1)
print(Z)
#The Wilcoxon test is based on the z-score

#ex2
library(dplyr)
data(ChickWeight)
head(ChickWeight)
mypar()
plot( ChickWeight$Time, ChickWeight$weight, col=ChickWeight$Diet)
chick = reshape(ChickWeight, idvar=c("Chick","Diet"), timevar="Time",
                direction="wide")
head(chick)
chick=na.omit(chick)
#q1:
cc<-as.data.frame(chick)
x<-filter(cc, Diet==1) #dplyr
x<-x$weight.4
y<-filter(cc, Diet==4)
y<-y$weight.4
t.test(x,y)
wilcox.test(x,y)
x[17]<-200
t.test(x,y)$p.value
#q2:
wilcox.test(x,y)$p.value
#q3:
x<-filter(cc, Diet==1) #dplyr
x<-x$weight.4
library(rafalib)
mypar(1,3)
boxplot(x,y)
boxplot(x,y+10)
boxplot(x,y+100)
t10<-t.test(x,y+10)$statistic
t100<-t.test(x,y+100)$statistic
t10-t100
#q4:
w10<-wilcox.test(x,y+10)$statistic
w100<-wilcox.test(x,y+100)$statistic
w10-w100 #serve solo a vedere che non cambia
#Il test di Wilcoxon è più debole del t-test
# nel caso in cui i campioni siano di piccole dimensioni,
# anche quando la differenza tra i due è elevata,
# proprio perché l'agire sulle mediane fa sì
# che il valore p possa essere molto piccolo
wilcox.test(c(1,2,3), c(4,5,6))$p.value
#q6:
wilcox.test(c(1,2,3), c(400,500,600))$p.value
