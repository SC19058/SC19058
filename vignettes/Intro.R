## ----eval=FALSE----------------------------------------------------------
#  d_stat <- function(size,group,index,distance){
#    total <- size * group
#    ind <- index;ind <- matrix(ind,byrow = TRUE,nrow = group)
#    add.w <- numeric(group)
#    for (i in 1:group) {
#      add.w[i] <- size/2 * mean(distance[ind[i,],ind[i,]])
#    }
#    w <- sum(add.w)
#    s <- numeric(1)
#    for(i in 1:(group-1)){
#      for(j in (i+1):group){
#        s = s + size^2/total *
#          ( 2*mean(distance[ind[i,],ind[j,]]) -
#              mean(distance[ind[j,],ind[j,]]) -
#              mean(distance[ind[i,],ind[i,]]) )
#      }
#    }
#    s/(group - 1)/(w/(total - group))
#  }

## ------------------------------------------------------------------------
equal_test <- function(dat,size,group,dim,R = 199){
  total.size <- size*group
  distance <- as.matrix(dist(dat,diag = TRUE,upper = TRUE,
                             method = "minkowski",p = 1))
  d_stat_storage <- numeric(R+1)
  index <- 1:(size*group)
  d_stat_storage[1] <- d_stat(size,group,index,distance)
  for(k in 1:R){
    index <- sample(1:total.size,total.size,replace = FALSE)
    d_stat_storage[k+1] <- d_stat(size,group,index,distance)
  }
  mean(d_stat_storage > d_stat_storage[1])
}

## ------------------------------------------------------------------------
library(SC19058)
pdim <- 10;delta <- 0.3;size = 30;group = 4
dat <- c(dat <- c(rt(n = size*pdim,df = 4,ncp = delta),
                  rt(n = size*pdim*(group - 1),df = 4)))
dat <- matrix(dat,byrow = TRUE,nrow = size*group)
equal_test(dat,size,group,dim)

## ----echo = FALSE--------------------------------------------------------
data(simulation)
library(ggplot2)

## ----echo=FALSE----------------------------------------------------------
delta <- seq(0,0.6,by = 0.05)
dat.power <- data.frame("delta" = delta,"DISCO"=powers1[1,],
                        "Pillai" = powers1[2,],
                        "Wilks" = powers1[3,])
ggplot(data = dat.power,aes(x = delta,y = empirical_powers)) + 
  lims(y=c(0,1)) +
  geom_line(linetype = 3,aes(y = DISCO),size = 0.8,color = "red") +
  geom_point(aes(y = DISCO,color = "DISCO"),size = 2) + 
  geom_line(linetype = 2,aes(y = Wilks),size = 0.8,color = "blue")  +
  geom_point(aes(y = Wilks,color = "Wilks"),size = 2) +
  geom_line(linetype = 4,aes(y = Pillai),size = 0.8,color = "green") + 
  geom_point(aes(y = Pillai,color = "Pillai"),size = 2) 

## ----echo=FALSE----------------------------------------------------------
pdim <- c(seq(10,110,by = 10))
dat.power <- data.frame("dimension" = pdim,"DISCO"=powers2[1,],
                        "Pillai" = powers2[2,],
                        "Wilks" = powers2[3,])
ggplot(data = dat.power,aes(x = dimension,y = empirical_powers)) + 
  lims(y=c(0,1)) +
  geom_line(linetype = 3,aes(y = DISCO),size = 0.8,color = "red") +
  geom_point(aes(y = DISCO,color = "DISCO"),size =2) + 
  geom_line(linetype = 2,aes(y = Wilks),size = 0.8,color = "blue")  +
  geom_point(aes(y = Wilks,color = "Wilks"),size = 2) +
  geom_line(linetype = 4,aes(y = Pillai),size = 0.8,color = "green")  +
  geom_point(aes(y = Pillai,color = "Pillai"),size = 2) 


## ----echo = FALSE--------------------------------------------------------
sigma <- c(seq(0,1.2,by = 0.1))
dat.power <- data.frame("sigma"=sigma,"DISCO"=powers3[1,],
                        "Pillai" = powers3[2,],
                        "Wilks" = powers3[3,])
ggplot(data = dat.power,aes(x = sigma,y = empirical_powers)) + 
  lims(y=c(0,1)) +
  geom_line(linetype = 3,aes(y = DISCO),size = 0.8,color = "red") +
  geom_point(aes(y = DISCO,color = "DISCO"),size = 2) + 
  geom_line(linetype = 2,aes(y = Wilks),size = 0.8,color = "blue")  +
  geom_point(aes(y = Wilks,color = "Wilks"),size = 2) +
  geom_line(linetype = 4,aes(y = Pillai),size = 0.8,color = "green")  +
  geom_point(aes(y = Pillai,color = "Pillai"),size = 2) 

## ----echo=FALSE----------------------------------------------------------
pdim <- c(seq(10,110,by = 10))
dat.power <- data.frame("dimension" = pdim,"DISCO"=powers4[1,],
                        "Pillai" = powers4[2,],
                        "Wilks" = powers4[3,])
ggplot(data = dat.power,aes(x = dimension,y = empirical_powers)) + 
  lims(y=c(0,1)) +
  geom_line(linetype = 3,aes(y = DISCO),size = 0.5,color = "red") +
  geom_point(aes(y = DISCO,color = "DISCO"),size = 2) + 
    geom_line(linetype = 2,aes(y = Wilks),size = 0.8,color = "blue")  +
  geom_point(aes(y = Wilks,color = "Wilks"),size = 2) +
  geom_line(linetype = 4,aes(y = Pillai),size = 0.8,color = "green")  +
  geom_point(aes(y = Pillai,color = "Pillai"),size = 2) 

## ------------------------------------------------------------------------
set.seed(5)
x = rnorm(10000)
hist(x)

## ------------------------------------------------------------------------
x <- c(0.10,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.20,0.21,0.23)
y <- c(42.0,43.5,45.0,45.5,45,47.5,49.0,53.0,50.0,55.0,55.0,60.0)
plot(x,y)

## ------------------------------------------------------------------------
lm.a <- lm(y ~ x+1)
plot(lm.a)

## ------------------------------------------------------------------------
# try another packages to plot
library(ggplot2)
# Stacked barchart
(pie <- ggplot(mtcars, aes(x = factor(1), fill = factor(cyl))) +
geom_bar(width = 1))
# Pie chart
pie + coord_polar(theta = "y")
# The bullseye chart
pie + coord_polar()

## ------------------------------------------------------------------------
library(survival)
library(knitr)
data(veteran)
kable(head(veteran),format = "markdown")

## ------------------------------------------------------------------------
kable(head(iris))

## ------------------------------------------------------------------------
kable(head(mtcars), caption = "Data of Cars")

## ------------------------------------------------------------------------
#根据定义生成服从Rayleigh分布的理论随机变量X1
set.seed(2019)
sigma = 2
n1 = 1e4
Y = rnorm(n1,0,sigma)
Z = rnorm(n1,0,sigma)
X1 = sqrt(Y^2+Z^2) #X1服从参数sigma为2的Rayleigh分布

## ------------------------------------------------------------------------
#采用逆变换法生成服从Rayleigh分布对的随机变量X2
U = runif(n1,0,1)
X2 = sigma*sqrt(-2*log(1-U))

## ------------------------------------------------------------------------
hist(X1,prob = TRUE,col = "pink",main = "Rayleigh distribution",xlab = "X1(pink) and X2(blue)")
hist(X2,prob = TRUE,add = TRUE, col = rgb(0,0.5,0.5,0.4))

## ------------------------------------------------------------------------
n2 <- 1e3
p1 = 0.75
y1 <- rnorm(n2, mean = 0, sd = 1)
y2 <- rnorm(n2, mean = 3, sd = 1)
r <- sample(0:1, size = 1000, replace = TRUE, prob = c(1-p1 , p1))
y <- y1*r + y2*(1-r)
hist(y, seq(-8,8,.05), prob = TRUE)
lines(density(y))

## ------------------------------------------------------------------------
n3 <- 1e5 
y1 <- rnorm(n3, mean = 0, sd = 1)
y2 <- rnorm(n3, mean = 3, sd = 1)
for (i in 1:9) {
  p1 = 0.1*i
  r <- sample(0:1, size = n3, replace = TRUE, prob = c(1-p1 , p1))
  y <- y1*r + y2*(1-r)
  hist(y, seq(-8,8,.05), prob = TRUE)
}

## ------------------------------------------------------------------------
# 写函数 利用barlett分解生成服从warshart分布的随机变量
generate_wishart <- function(n,d,Sigma){
  #首先保证输入的参数满足题设条件
  if((n<=d+1)||(d<=0)){
    stop("n,d不符合要求！")
  }
  #生成d维零方阵
  T1 <- matrix(0,nrow=d,ncol=d)
  #赋值使其成为符合要求的下三角随机矩阵
  for(i in 1:d){
    T1[i,i] <- sqrt(rchisq(1,n-i+1))
  }
  for(i in 2:d){
    for(j in 1:i-1){
      T1[i,j] <- rnorm(1,mean = 0, sd = 1)
    }
  }
  A <- T1 %*% t(T1)
  #对协方差阵进行分解
  L <- chol(Sigma)
  W <- L %*% A %*% t(L)
  W
}

## ------------------------------------------------------------------------

Sigma1 = matrix(c(6,2,3,4,2,10,5,1,3,5,9,2,4,1,2,7),nrow = 4)
Sigma1

## ------------------------------------------------------------------------
W2 <- generate_wishart(6,4,Sigma1)
W2

## ------------------------------------------------------------------------
set.seed(50)
n <- 1e5 #蒙特卡罗方法在区间中生成变量个数
estimate.1 <- mean(sin(runif(n,0,pi/3))) * (pi/3)
cat(paste0("蒙特卡罗估计值为",estimate.1,"\n","积分精确值为",0.5))

## ------------------------------------------------------------------------
set.seed(100)
n <- 1e4 #蒙特卡罗积分法在[0,1]区间的变量个数
m <- 1e3 #重复试验次数
f <- function(x) exp(-x)/(1+x^2) #被积函数
#设计函数生成对偶变量估计积分值
MC.anti <- function( n = 1e4, antithetic = TRUE){
  u <- runif(n/2,0,1)
  if(!antithetic) 
    #输入值为false,生成独立同分布的另一组变量
    v <- runif(n/2,0,1)  
  else
    #输入值为true,采用对偶变量法生成另一组变量
    v <- 1-u
  u <- c(u,v)
  #返回积分估计值
  return(mean(f(u)))
}
MC1 <- MC2 <- numeric(m)
#多次重复试验
for(i in 1:m){
   MC1[i] <- MC.anti(n = 1e4, antithetic = TRUE)
   MC2[i] <- MC.anti(n = 1e4, antithetic = FALSE)
}
cat(" 对偶变量法估计值为",mean(MC1),"\n","单一蒙特卡罗法估计值为",mean(MC2),"\n","对偶变量法标准差为",sd(MC1),"\n","单一蒙特卡罗法标准差为",sd(MC2),"\n","方差缩减比为",(var(MC2)-var(MC1))/var(MC2))

## ------------------------------------------------------------------------
#找到子区间的分位数
quantile <- numeric(6)
for(i in 1:5){
  quantile[i+1] <- -log(1-0.2*i*(1-exp(-1)))
}
quantile

## ------------------------------------------------------------------------
set.seed(50)
n <- 1e4 #生成变量个数
m <- 1e4 #重复试验次数
#h为原被积函数除以重要函数
h <- function(x) (1-exp(-1))/(5*(1+x^2))
#MC.strat储存多次重复试验的估计结果
MC.strat <- numeric(m)
for(j in 1:m){
  x <- matrix(0,nrow = n,ncol = 5)
  for(i in 1:5){
    y <- runif(n/5,0,1)
    #采用逆变换法在子区间上生成满足重要函数密度的抽样
    x[,i] <- -log(y*(exp(-1)-1)/5+exp(-quantile[i]))
  }
  #将各子区间的积分估计值相加得到[0,1]区间上积分值
  MC.strat[j] <- sum(colMeans(h(x)))
}
#与来自于例5.13的原估计值和估计标准误差进行比较
sd <- 0.0970314
theta <- 0.5257801
cat(" 原方法估计值为",theta,"\n","分层重要抽样方法的积分估计值为",mean(MC.strat),"\n","原方法估计标准差为",sd^2,"\n","分层重要抽样方法的标准差估计值为",sd(MC.strat),"\n","方差缩减比为",(sd^2-var(MC.strat))/sd^2)

## ------------------------------------------------------------------------
n <- 20  # the sample size
m <- 1000 # the times of replications
alpha <- .05
UCL <- replicate(m,expr = {
  x <- rchisq(n, df = 2)
  (n-1) * var(x) / qchisq(alpha, df = n-1)
})
cat("the confidence interval covers the variance is", mean(UCL > 4))

## ------------------------------------------------------------------------
n <- 20 # the sample size
alpha <- .05
x <- rchisq(n, df = 2)
# if.in.confidence is a vector containing TRUE or FALSE, indicate if the confidence interval covers the mean 
if.in.confidence.interval <- replicate(m,expr = {
  x <- rchisq(n, df = 2)
  abs(mean(x)-2) < sd(x) * qt(alpha/2, df = n-1,lower.tail = FALSE)/sqrt(n)
})
cat("the confidence interval covers the mean is", mean(if.in.confidence.interval))

## ------------------------------------------------------------------------
n <- 20 ## the sample size
q <- c(0.025, 0.05, 0.95, 0.975) # the selected quantiles

## ------------------------------------------------------------------------
skewness <- replicate(10*m,expr = {
  x <- rnorm(n, mean = 0,sd = 1)
  xbar <- mean(x)
  m3 <- mean((x-xbar)^3)
  m2 <- mean((x-xbar)^2)
  m3 / m2^1.5
})

## ------------------------------------------------------------------------
var <- matrix ( q*(1-q)/n/pnorm(q,mean = 0,sd = sqrt( 6*(n-2)/(n+1)/(n+3) )), nrow = 1)
sd <- sqrt(var)
colnames(sd) <- q
rownames(sd) <- 'standard error of the estimates'
knitr::kable(sd)

## ------------------------------------------------------------------------
A <- matrix(c(quantile(skewness,q),qnorm(q,mean = 0,sd = sqrt(6/n))),nrow = 2,byrow = TRUE)
colnames(A) <- q
rownames(A) <- c('estimated quantiles','quantiles of large sample')
knitr::kable(A)

## ------------------------------------------------------------------------
B <- matrix(c(quantile(skewness,q),qnorm(q,mean = 0,sd = sqrt( 6*(n-2)/(n+1)/(n+3) ))),nrow = 2,byrow = TRUE)
colnames(B) <- q
rownames(B) <- c('estimated quantiles','quantiles with accurate variance')
knitr::kable(B)

## ------------------------------------------------------------------------
skewness <- function(x){
  #compute the sample skewness
  xbar <- mean(x)
  m3 <- mean((x-xbar)^3)
  m2 <- mean((x-xbar)^2)
  return( m3 / m2^1.5 )
}

## ------------------------------------------------------------------------
set.seed(100)
alpha <- 0.05 # the significance value
a = seq(1,100,by = 1) # parameter in beta distributions
N <- length(a)
powers1 <- numeric(N)
n <- 100  # numbers of variable for each simulation
m <- 1000 # replication times
# critical value for the skewness test 
cv <- qnorm(1-alpha/2,0,sqrt(6*(n-2)/(n+1)/(n+3)))

for(i in 1:N){
  sktest <- numeric(m)
  for(k in 1:m){
    x <- rbeta(n,a[i],a[i])
    sktest[k] <- as.integer(abs(skewness(x)) >= cv)
  }
  powers1[i] <- mean(sktest)
}

## ------------------------------------------------------------------------
library(ggplot2)
ggplot(data = data.frame(powers1), aes(x = a, y = powers1)) + 
         geom_point() +
         labs(x = "parameter of Beta(alpha,alpha)",y = "powers") +
         geom_smooth()

## ------------------------------------------------------------------------
set.seed(100)
alpha <- 0.05 # the significance value
mu = seq(1,20,by = 0.1) # parameter in t distributions
N <- length(mu)
powers2 <- numeric(N)
n <- 100  # numbers of variable for each simulation
m <- 1000 # replication times
# critical value for the skewness test 
cv <- qnorm(1-alpha/2,0,sqrt(6*(n-2)/(n+1)/(n+3)))

for(i in 1:N){
  sktest <- numeric(m)
  for(k in 1:m){
    x <- rt(n,mu[i])
    sktest[k] <- as.integer(abs(skewness(x)) >= cv)
  }
  powers2[i] <- mean(sktest)
}

## ------------------------------------------------------------------------
ggplot(data = data.frame(powers2), aes(x = mu, y = powers2)) + 
         geom_point() +
         labs(x = "parameter of t(mu)",y = "powers") +
         geom_smooth()

## ------------------------------------------------------------------------
m <- 50 # numbers of variable for each simulation
n <- 1e4 # replication times
mu0 <- 1
t1e <- matrix(nrow = 1,ncol = 3)
colnames(t1e) <- c("chisq","uniform","exponential")
test1 <- test2 <- test3 <-numeric(n)
for(i in 1:n){
  x <- rchisq(m,df = 1)
  test1[i] <- t.test(x,mu = mu0)$p.value
  y <- runif(m,0,2)
  test2[i] <- t.test(y,mu = mu0)$p.value
  z <- rexp(m, rate = 1)
  test3[i] <- t.test(z,mu = mu0)$p.value
}
t1e[1,1] <- mean(test1 <= 0.05)
t1e[1,2] <- mean(test2 <= 0.05) 
t1e[1,3] <- mean(test3 <= 0.05)
knitr::kable(t1e)

## ----echo=FALSE----------------------------------------------------------
data <- matrix(nrow = 3,ncol = 3)
data[1,1:3] <- c("a","b","a + b")
data[2,1:3] <- c("c","d","c + d")
data[3,1:3] <- c("a + c","b + d","10000")
colnames(data) <- c("Test2 reject","Test2 accpet","row sum")
rownames(data) <- c("Test1 reject","Test1 accpet","column sum")
knitr::kable(data)

## ------------------------------------------------------------------------
library(bootstrap)
pairs(scor)
cor(scor)

## ------------------------------------------------------------------------
library(ggplot2)
library(GGally)
ggpairs(scor)

## ------------------------------------------------------------------------
boot.std.grades.cor <- function(replicates = 200, sample.size = dim(scor)[1],i,j){
  ## `replicates` is the number of replication
  ## `sample.size` is the size of generated sample in each replications
  ## `i`,`j` represents the indice of subjects in dataframe "scor"\
  storage.cor <- numeric(replicates) #storage for replicates
  for(r in 1:replicates){
    index <- sample(1:sample.size, size = sample.size, replace = TRUE)
    storage.cor[r] <- cor(scor[index,i],scor[index,j])
  }
  return(list(std = sd(storage.cor), mean = mean(storage.cor), simulation = storage.cor))
}

## ------------------------------------------------------------------------
ro12 <- boot.std.grades.cor(i = 1,j = 2)
ro34 <- boot.std.grades.cor(i = 3,j = 4)
ro35 <- boot.std.grades.cor(i = 3,j = 5)
ro45 <- boot.std.grades.cor(i = 4,j = 5)

## ------------------------------------------------------------------------
hist(ro12$simulation,probability = TRUE, xlab = "cor(mec,vec)",main = "cor(mec,vec) from replications")
hist(ro34$simulation,probability = TRUE, xlab = "cor(alg,ana)",main = "cor(alg,ana) from replications")
hist(ro35$simulation,probability = TRUE, xlab = "cor(alg,sta)",main = "cor(alg,sta) from replications")
hist(ro45$simulation,probability = TRUE, xlab = "cor(ana,sta)",main = "cor(ana,sta) from replications")

## ------------------------------------------------------------------------
standard.error <- c(ro12$std,ro34$std,ro35$std,ro45$std)
mean <- c(ro12$mean,ro34$mean,ro35$mean,ro45$mean)
real <- c(cor(scor$mec,scor$vec),cor(scor$alg,scor$ana),cor(scor$alg,scor$sta),cor(scor$ana,scor$sta))
standard.error.matrix <- data.frame("subjects"=c("(mec,vec)","(alg,ana)","(alg,sta)","(ana,sta)"),
                                    "cor from sample" = real,
                                "estimated cor from bootstrap" = mean,
                                  "standard error" = standard.error
                                    )
knitr::kable(standard.error.matrix)

## ------------------------------------------------------------------------
skewness <- function(x){
  #compute the sample skewness
  xbar <- mean(x)
  m3 <- mean((x-xbar)^3)
  m2 <- mean((x-xbar)^2)
  return( m3 / m2^1.5 )
}

## ------------------------------------------------------------------------
library(boot);set.seed(555)
mu <- 0  # the skewness of normal populations
n <- 20  # the sample size in bootstrap
m <- 1e3 # replicate time in Monte Carlo experiments
boot.skew <- function(x,i) skewness(x[i])
ci.norm <- ci.basic <- ci.perc <- matrix(NA,m,2)
for(i in 1:m){
  X <- rnorm(n,mean = 0,sd = 1)
  de <- boot(data = X, statistic = boot.skew, R = 1e3)
  ci <- boot.ci(de,type = c("norm","basic","perc"))
  ci.norm[i,] <- ci$norm[2:3]
  ci.basic[i,] <- ci$basic[4:5]
  ci.perc[i,] <- ci$percent[4:5]
}

## ------------------------------------------------------------------------
cover.prob <- c(mean(ci.norm[,1]<= mu & ci.norm[,2]>= mu),
                mean(ci.basic[,1]<= mu & ci.basic[,2]>= mu),
                mean(ci.perc[,1]<= mu & ci.perc[,2]>= mu))
left.omit <- c(mean(ci.norm[,1]>= mu),
               mean(ci.basic[,1]>= mu),
               mean(ci.perc[,1]>= mu))
right.omit <- c(mean(ci.norm[,2]<= mu),
               mean(ci.basic[,2]<= mu),
               mean(ci.perc[,2]<= mu))
cover.norm <- matrix(data = c(cover.prob,left.omit,right.omit),nrow = 3,byrow = TRUE,)
rownames(cover.norm) <- c("cover probability","miss on the left","miss on the right")
colnames(cover.norm) <- c("standard normal bootstrap confidence interval","basic bootstrap confidence interval","percentile bootstrap confidence interval")
knitr::kable(t(cover.norm))

## ------------------------------------------------------------------------
library(boot);set.seed(555)
mu <- sqrt(8/5)  # the skewness of chi-squared distribution
n <- 20  # the sample size in bootstrap
m <- 1e3 # replicate time in Monte Carlo experiments
boot.skew <- function(x,i) skewness(x[i])
ci.norm <- ci.basic <- ci.perc <- matrix(NA,m,2)
for(i in 1:m){
  X <- rchisq(n,df = 5)
  de <- boot(data = X, statistic = boot.skew, R = 1e3)
  ci <- boot.ci(de,type = c("norm","basic","perc"))
  ci.norm[i,] <- ci$norm[2:3]
  ci.basic[i,] <- ci$basic[4:5]
  ci.perc[i,] <- ci$percent[4:5]
}

## ------------------------------------------------------------------------
cover.prob <- c(mean(ci.norm[,1]<= mu & ci.norm[,2]>= mu),
                mean(ci.basic[,1]<= mu & ci.basic[,2]>= mu),
                mean(ci.perc[,1]<= mu & ci.perc[,2]>= mu))
left.omit <- c(mean(ci.norm[,1]>= mu),
               mean(ci.basic[,1]>= mu),
               mean(ci.perc[,1]>= mu))
right.omit <- c(mean(ci.norm[,2]<= mu),
               mean(ci.basic[,2]<= mu),
               mean(ci.perc[,2]<= mu))
cover.chisq <- matrix(data = c(cover.prob,left.omit,right.omit),nrow = 3,byrow = TRUE,)
rownames(cover.chisq) <- c("cover probability","miss on the left","miss on the right")
colnames(cover.chisq) <- c("standard normal bootstrap confidence interval","basic bootstrap confidence interval","percentile bootstrap confidence interval")
knitr::kable(t(cover.chisq))

## ------------------------------------------------------------------------
cov.mle <- function(a){
  a <- as.matrix(a)
  k <- ncol(a)
  for(i in 1:k){
    a[,i] <- a[,i] - colMeans(a)[i]
  }
  t(a) %*% a / nrow(a)
}

## ------------------------------------------------------------------------
library(bootstrap)
n <- nrow(scor)
eigen.value <- eigen(cov.mle(scor))$value
lamda.hat <- max(eigen.value)/sum(eigen.value)
lamda.jack <- numeric(n)
# resampling using jackknife method
for (i in 1:n) {
  scor.jacknife <- scor[-i,]
  eigen.value <- eigen(cov.mle(scor.jacknife))$value
  lamda.jack[i] <- max(eigen.value)/sum(eigen.value)
}
bias <- (n-1)*(mean(lamda.jack)-lamda.hat)
se <- sqrt((n-1) * mean((lamda.jack - mean(lamda.jack))^2))
knitr::kable(data.frame("estimate of bias"=bias,"estimate of standard error"=se))

## ------------------------------------------------------------------------
library(DAAG);attach(ironslag)
L1 <- lm(magnetic ~ chemical)
L2 <- lm(magnetic ~ chemical + I(chemical^2))
L3 <- lm(magnetic ~ chemical + I(chemical^2) + I(chemical^3))
L4 <- lm(log(magnetic) ~ chemical)

## ----echo=FALSE----------------------------------------------------------
a <- seq(10,40, .1)  #sequence for plotting fits

plot(chemical, magnetic, main="Linear", pch=16)
yhat1 <- L1$coef[1] + L1$coef[2] * a
lines(a, yhat1, lwd=2)

plot(chemical, magnetic, main="Quadratic", pch=16)
yhat2 <- L2$coef[1] + L2$coef[2] * a + L2$coef[3] * a^2
lines(a, yhat2, lwd=2)

plot(chemical, magnetic, main="Cubic", pch=16)
yhat3 <- L3$coef[1] + L3$coef[2] * a + L3$coef[3] * a^2 + L3$coef[4] * a^3
lines(a, yhat3, lwd=2)

plot(chemical, magnetic, main="Exponential", pch=16)
logyhat4 <- L4$coef[1] + L4$coef[2] * a
yhat4 <- exp(logyhat4)
lines(a, yhat4, lwd=2)

## ------------------------------------------------------------------------
n <- length(magnetic)
e1 <- e2 <- e3 <- e4 <- numeric(n)
#for n-fold cross validation
#fit models on leave-one-out samples
for(k in 1:n){
  y <- magnetic[-k]
  x <- chemical[-k]
  #linear
  J1 <- lm(y~x)
  yhat1 <- J1$coef[1] + J1$coef[2]*chemical[k]
  e1[k] <- magnetic[k] - yhat1
  #quadratic
  J2 <- lm(y~x + I(x^2))
  yhat2 <- J2$coef[1] + J2$coef[2]*chemical[k] + J2$coef[3] * chemical[k]^2
  e2[k] <- magnetic[k] - yhat2
  #cubic
  J3 <- lm(y~x + I(x^2) + I(x^3))
  yhat3 <- J3$coef[1] + J3$coef[2]*chemical[k] + 
           J3$coef[3] * chemical[k]^2 + J3$coef[4] * chemical[k]^3
  e3[k] <- magnetic[k] - yhat3
  #exponential
  J4 <- lm(log(y)~log(x))
  logyhat4 <- J4$coef[1] + J4$coef[2] * log(chemical[k])
  yhat4 <- exp(logyhat4)
  e4[k] <- magnetic[k] - yhat4
}
bias <- data.frame("linear" = mean(e1^2),
                   "quadratic" = mean(e2^2),
                   "cubic" = mean(e3^2),
                   "exponential" = mean(e4^2))
knitr::kable(bias)

## ------------------------------------------------------------------------
R.squared <- data.frame("linear" = summary(L1)$adj.r.sq,
               "quadratic" = summary(L2)$adj.r.sq,
               "cubic" = summary(L3)$adj.r.sq,
               "exponential" = summary(L4)$adj.r.sq)
knitr::kable(R.squared)

## ------------------------------------------------------------------------
count_extreme_points <- function(sample1,sample2){
  x <- sample1;y <- sample2
  X <- x - mean(x)
  Y <- y - mean(y)
  outx <- sum(X > max(Y)) + sum(X < min(Y))
  outy <- sum(Y > max(X)) + sum(Y < min(X))
  max(outx,outy)
} 

## ------------------------------------------------------------------------
alpha <- 0.05    #significance level
n1 <- 30;n2 <-50 #two different sample size
mu1 <- mu2 <- 0;sigma1 <- sigma2 <- 1
m <- 1000   # the times of Monte Carlo experiments

p_value <- replicate(m,expr={
  x1 <- rnorm(n1,mu1,sigma1)
  x2 <- rnorm(n2,mu2,sigma2)
  ts <- numeric(199+1)
  ts[1] <- count_extreme_points(x1,x2)
  for(i in 1:199){
    ind <- sample(1:(n1+n2),size = n1,replace = FALSE)
    x.perm <- c(x1,x2)[ind]; y.perm <- c(x1,x2)[-ind]
    ts[i+1] <- count_extreme_points(x.perm,y.perm)
  }
  mean(ts >= ts[1])
})
#estimate the type1 error
print(mean(p_value < alpha))

## ------------------------------------------------------------------------
library(MASS);library(boot);library(Ball)
dCov <-function(x, y) {
  x <-as.matrix(x);y <-as.matrix(y)
  n <-nrow(x); m <-nrow(y)
  if(n!=m||n<2) stop("Sample sizes must agree")
  if(!(all(is.finite(c(x, y)))))
    stop("Data contains missing or infinite values")
  Akl <-function(x) {
    d <-as.matrix(dist(x))
    m <-rowMeans(d); 
    M <-mean(d)
    a <-sweep(d, 1, m);
    b <-sweep(a, 2, m)
    return(b+M)
  }
  A <-Akl(x); 
  B <-Akl(y)
  sqrt(mean(A*B))
}

ndCov2 <-function(z, ix, dims) {#dims contains dimensions of x and y
  p <- dims[1]
  q <- dims[2]
  d <- p+q
  x <- z[ , 1:p]#leave x as is
  y <- z[ix,-(1:p)]#permute rows of y
  return(nrow(z)* dCov(x, y)^2)
  }

## ------------------------------------------------------------------------
n <- c(30,50,80,100,120,150,200)
alpha <- 0.1 #significance level
m <- 200 # Monte Carlo times
powers1 <- data.frame(matrix(0,nrow = 2,ncol = length(n)))
for(i in 1:length(n)){
  p.cor <- numeric(m);p.ball <- numeric(m)
  for(j in 1:m){
    mu = c(0,0); sigma = matrix(c(1,0,0,1),nrow = 2)
    x1 <- mvrnorm(n[i],mu,sigma); e <- mvrnorm(n[i],mu,sigma)
    y1 <- x1/4 + e
    z <- matrix(c(x1,y1),nrow = n[i])
    # permutation: resampling without replacement
    boot.obj <- boot(data = z, statistic = ndCov2, R = 99,
                 sim = "permutation", dims =c(2, 2))
    tb <- c(boot.obj$t0, boot.obj$t)
    p.cor[j] <- mean(tb >= tb[1])
    p.ball[j] <- as.numeric(bcov.test(z[,1:2],z[,3:4],R = 99,seed = j)$p.value)
  }
  powers.cor <- mean(p.cor < alpha)
  powers.ball <- mean(p.ball < alpha)
  powers1[,i] <- c(powers.cor,powers.ball)
}

## ------------------------------------------------------------------------
rownames(powers1) <- c("distance correlation test","ball covariance test")
name.str <- c()
for(i in 1:length(n)){
  name.str <- c(name.str,paste("n is",n[i]))
}
colnames(powers1) <- name.str
knitr::kable(powers1)

## ------------------------------------------------------------------------
powers2 <- data.frame(matrix(0,nrow = 2,ncol = length(n)))
for(i in 1:length(n)){
  p.cor <- numeric(m);p.ball <- numeric(m)
  for(j in 1:m){
    mu = c(0,0); sigma = matrix(c(1,0,0,1),nrow = 2)
    x1 <- mvrnorm(n[i],mu,sigma); e <- mvrnorm(n[i],mu,sigma)
    y1 <- x1/4 * e
    z <- matrix(c(x1,y1),nrow = n[i])
    # permutation: resampling without replacement
    boot.obj <- boot(data = z, statistic = ndCov2, R = 99,
                 sim = "permutation", dims =c(2, 2))
    tb <- c(boot.obj$t0, boot.obj$t)
    p.cor[j] <- mean(tb >= tb[1])
    p.ball[j] <- as.numeric(bcov.test(z[,1:2],z[,3:4],R = 99,seed = j)$p.value)
  }
  powers.cor <- mean(p.cor < alpha)
  powers.ball <- mean(p.ball < alpha)
  powers2[,i] <- c(powers.cor,powers.ball)
}

## ------------------------------------------------------------------------
rownames(powers2) <- c("distance correlation test","ball covariance test")
name.str <- c()
for(i in 1:length(n)){
  name.str <- c(name.str,paste("n is",n[i]))
}
colnames(powers2) <- name.str
knitr::kable(powers2)

## ------------------------------------------------------------------------
par(mfrow = c(1,2))
plot(n,as.vector(powers1[1,]),xlim = c(10,210),ylim = c(0,1),xlab = "sample size",ylab = "power",main = "power comparision in model 1",type = "l",col = "red")
lines(n,as.vector(powers1[2,]),col = "blue")
legend("bottomright",lty = c(1,1),col=c("red","blue"),legend = c("distance correlation","ball covariance"))
plot(n,as.vector(powers2[1,]),xlim = c(10,210),ylim = c(0,1),xlab = "sample size",ylab = "power",main = "power comparision in model 2",type = "l",col = "red")
lines(n,as.vector(powers2[2,]),col = "blue")
legend("bottomright",lty = c(1,1),col=c("red","blue"),legend = c("distance correlation","ball covariance"))

## ------------------------------------------------------------------------
r <- function(x,y){
  exp(abs(x)-abs(y))
}

## ------------------------------------------------------------------------
rw.Metropolis <- function(x0,sigma,N) {
  #x0: the initial point
  #sigma: the standard deviation in the normal distribution
  #N: the length of the chain
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  accept <- 1 #count reject
  for(i in 2:N){
    y <- rnorm(1,x[i-1],sigma)
    if(u[i] <= r(x[i-1],y)){
       x[i] <- y;accept <- accept + 1
    }else{
      x[i] <- x[i-1]
    }
  }
  return(list(x = x,accept = accept))
}

## ------------------------------------------------------------------------
N <- 5000
sigma <- c(0.05,0.5,2,4,8,16)
x0 <- 20
for(i in 1:length(sigma)){
  assign(paste0("rw",i),rw.Metropolis(x0,sigma[i],N))  
}

## ------------------------------------------------------------------------
# quantile for laplace distribution(0.025,0.975)
quantile <- c(log(0.05),-log(0.05))
index <- 1:N
for(i in 1:length(sigma)){
  plot(index,get(paste0("rw",i))$x,type = "l",ylab = "x",xlab = "",
       main = bquote(sigma == .(sigma[i])))
  abline(h = quantile,col ="red",lty = 3)
}

## ------------------------------------------------------------------------
a <- matrix(0,nrow = 1,ncol = length(sigma))
rownames(a) <- "accpetance rate"
col_name <- character(length(sigma))
for(i in 1:length(sigma)){
  col_name[i] <- paste0("sigma = ",sigma[i])
  a[1,i] <- get(paste0("rw",i))$accept / N
}
colnames(a) <- col_name
knitr::kable(a)

## ------------------------------------------------------------------------
a <- c(0.1,0.2,0.3,0.4,0.45)
# quantile for laplace distribution(0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95)
q <- c(log(1-2*rev(a)),0,-log(1-2*a))
a1 <- c(rev(0.5-a),0.5,0.5 + a)
rw <- cbind(rw1$x,rw2$x,rw3$x,rw4$x,rw5$x,rw6$x)
mc <- rw[501:N,]
# quantiles of generated chain
Qrw <- apply(mc,2,function(x) quantile(x,a1))
b <- round(cbind(q,Qrw),3)
colnames(b) <- c("quantile",col_name)
knitr::kable(b)

## ------------------------------------------------------------------------
# q: quantiles of Laplace
# Qrw[,4]: quantiles of the chain with sigma is 4
  par(mfrow=c(1,2))
    qqplot(q, Qrw[,4], main="",
        xlab="standard Laplace Quantiles", ylab="Sample Quantiles")
    abline(0,1,col='blue',lwd=2)
    hist(rw4$x[501:N], breaks="scott", main="", xlab="", freq=FALSE)
    lines(q, 0.5*exp(-abs(q)),col = "red")

## ------------------------------------------------------------------------
x <- 3.14
# the formula below is not true in computer arithmetic
exp(log(x)) == log(exp(x))
# the all.equal test if two objects are nearly equal
all.equal(exp(log(x)),log(exp(x)))

## ------------------------------------------------------------------------
# the integration
integ = function(k,u){
  return((1+u^2/(k))^(-(k+1)/2))
}
#calculate c_k
c_k <- function(k,a){
  sqrt(a^2*k/(k+1-a^2))
}
# the equation in log 
equation <- function(k,a){
   (log(k) - log(k-1))/2 + 
   (2*lgamma(k/2) - lgamma((k+1)/2)- lgamma((k-1)/2)) +
   log(integrate(integ,0,c_k(k-1,a),k = k-1)$value) -
   log(integrate(integ,0, c_k(k,a),k = k)$value)
}
root.equation <- sapply(c(4:25,100,500,1000),function(k){
                           uniroot(equation,interval = c(1,2),k=k)$root})

## ------------------------------------------------------------------------
equation_11_4 <- function(k,a){
  pt(c_k(k-1,a),df = k-1) - pt(c_k(k,a),df = k)
}
root.curve <- sapply(c(4:25,100,500,1000),function(k){uniroot(equation_11_4,interval = c(1,2),k=k)$root})

## ------------------------------------------------------------------------
result <- cbind(root.equation,root.curve)
colnames(result) <- c("root in 11.5","root in 11.4")
rownames(result) <- as.character(c(4:25,100,500,1000))
knitr::kable(result)

## ------------------------------------------------------------------------
set.seed(125)
n_a. <- 28
n_b. <- 24
n_ab <- 70
n_oo <- 41
p0 <- runif(1,0,1)
q0 <- runif(1,0,1-p0)

## ------------------------------------------------------------------------
likelihood_e <- function(prob,p0,q0){
  r0 <- 1-p0-q0 
  p <- prob[1]; q <- prob[2] ; r <- 1-p-q
  - n_a. * (2*log(p)*(p0^2/(p0^2+2*p0*r0)) + log(2*p*r)*(2*p0*r0/(p0^2+2*p0*r0))) -
    n_b. * (2*log(q)*(q0^2/(q0^2+2*q0*r0)) + log(2*q*r)*(2*q0*r0/(q0^2+2*q0*r0))) -
    n_ab * log(2*p*q) - 2*n_oo * log(r^2) 
}

## ------------------------------------------------------------------------
iter <- 0;E1 <- 0;E2 <- 1
while(iter < 200 && abs(E1-E2)> 1e-6){
  output <- optim(par = c(0.1,0.1),likelihood_e,p0 = p0,q0 = q0)
  E1 <- E2;E2 <- output$value
  p0 <- output$par[1]
  q0 <- output$par[2]
  iter <- iter + 1
}
estimate <- data.frame(p0,q0,iter)
colnames(estimate) <- c("p","q","iteration times")
knitr::kable(estimate)

## ------------------------------------------------------------------------
formulas <- list(
  mpg ~ disp,
  mpg ~ I(1 / disp),
  mpg ~ disp +wt,
  mpg ~ I(1 / disp) + wt
)

## ------------------------------------------------------------------------
models.loop <- list()
for(i in 1:length(formulas)){
  models.loop <- c(models.loop,list(lm(formulas[[i]],data = mtcars)))
}

## ------------------------------------------------------------------------
models <- lapply(formulas,lm,data = mtcars)
models

## ------------------------------------------------------------------------
bootstraps <- lapply(1:10, function(i){
  rows <- sample(1:nrow(mtcars), rep =TRUE)
  mtcars[rows, ]
})

## ------------------------------------------------------------------------
models2.loop <- list()
for(i in 1:10){
  models2.loop <- c(models2.loop,list(lm(mpg ~ disp,data = bootstraps[[i]])))
}

## ------------------------------------------------------------------------
models2 <- lapply(bootstraps,lm,formula = mpg~disp)
models2

## ------------------------------------------------------------------------
rsq <- function(mod) summary.lm(mod)$r.squared
# The r.squared from four different models
rsq1 <- matrix(c(unlist(sapply(models.loop,rsq)),unlist(sapply(models,rsq))),
               byrow = TRUE, nrow = 2,dimnames = c(list(c("loop","apply")),
                                                   list(c("mpg ~ disp",
                                                          "mpg ~ I(1 / disp)",
                                                          "mpg ~ disp +wt",
                                                          "mpg ~ I(1 / disp) + wt"))))
knitr::kable(rsq1)
# The r.squared from bootstrap replicates
rsq2<-matrix(c(round(unlist(sapply(models2.loop,rsq)),3),round(unlist(sapply(models2,rsq)),3)),
               byrow = TRUE, nrow = 2,dimnames = c(list(c("loop","apply")),
                                                   list(1:10)))
knitr::kable(rsq2)

## ------------------------------------------------------------------------
trials <- replicate(
  100,
  t.test(rpois(10,10),rpois(7,10)),
  simplify = FALSE
)

## ------------------------------------------------------------------------
round(sapply(1:100,function(i){trials[[i]]$p.value}),3)

## ------------------------------------------------------------------------
round(sapply(trials,"[[","p.value"),3)

## ------------------------------------------------------------------------
sapply

## ------------------------------------------------------------------------
library(parallelsugar)
mcsapply <- function (X, FUN, ..., mc.cores,simplify = TRUE, USE.NAMES = TRUE) {
  FUN <- match.fun(FUN)
  answer <- mclapply(X = X, FUN = FUN, ...,mc.cores = mc.cores)
  if (USE.NAMES && is.character(X) && is.null(names(answer))) 
    names(answer) <- X
  if (!isFALSE(simplify) && length(answer)) 
    simplify2array(answer, higher = (simplify == "array"))
  else answer
}

## ------------------------------------------------------------------------
system.time(sapply(1:1e3, function(i) rnorm(5e5)))

## ------------------------------------------------------------------------
system.time(mcsapply(1:1e3, function(i) rnorm(5e5), mc.cores = 4))

## ------------------------------------------------------------------------
rw.Metropolis <- function(x0,sigma,N) {
  #x0: the initial point
  #sigma: the standard deviation in the normal distribution
  #N: the length of the chain
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  accept <- 1 #count reject
  for(i in 2:N){
    y <- rnorm(1,x[i-1],sigma)
    if(u[i] < exp(abs(x[i-1])-abs(y))){
       x[i] <- y;accept <- accept + 1
    }else{
      x[i] <- x[i-1]
    }
  }
  return(list(x = x,accept = accept))
}

## ------------------------------------------------------------------------
library(Rcpp)
cppFunction('List rw_Metropolis_c(double x0, double sigma, int N) {
  NumericVector x(N);
  as<DoubleVector>(x)[0] = x0;
  NumericVector u(N);
  u = as<DoubleVector>(runif(N));
  List out(2);
  int accept = 1;
  for(int i=1;i<N;i++){
    double y = as<double>(rnorm(1,x[i-1],sigma));
    if(u[i] <= exp(abs(x[i-1])-abs(y))){
        x[i] = y;accept = accept + 1;
    }
    else{
        x[i] = x[i-1];
    }
  }  
  out[0] = x;
  out[1] = accept;
  return(out);
}')

## ------------------------------------------------------------------------
N = 5000; sigma = c(0.05,0.5,2,8);x0 = 10;
for(i in 1:length(sigma)){
  assign(paste0("chain",i),rw.Metropolis(x0,sigma[i],N))  
  assign(paste0("chain_c",i),rw_Metropolis_c(x0,sigma[i],N))  
}
for(i in 1:length(sigma)){
  par(mfrow = c(1,2))
  plot(get(paste0("chain",i))$x,type = "l", ylab = "from R",
       main = bquote(sigma == .(sigma[i])))
  plot(get(paste0("chain_c",i))[[1]],type = "l", ylab = "from Rcpp",
       main = bquote(sigma == .(sigma[i])))
}


## ------------------------------------------------------------------------
a <- data.frame(0)
for(i in 1:length(sigma)){
  a <- cbind(a,round(c(get(paste0("chain",i))$accept/N
                       ,get(paste0("chain_c",i))[[2]]/N),2))
  colnames(a)[i+1] <- paste0("sigma = ",sigma[i])
}
a <- a[2:5]
rownames(a) <- c("from R","from Rcpp")
knitr::kable(a)

## ------------------------------------------------------------------------
for(i in 1:length(sigma)){
  qqplot(get(paste0("chain",i))$x,
         get(paste0("chain_c",i))[[1]],
         xlab = "from R",ylab = "from Rcpp",
         main = bquote(sigma == .(sigma[i])))
  f <- function(x) x
  curve(f, col = 'red',add = TRUE)
}


## ------------------------------------------------------------------------
library(microbenchmark)
b <- data.frame(0)
ts1 <- microbenchmark(chain = rw.Metropolis(x0,sigma[1],N),
                     chain_c = rw_Metropolis_c(x0,sigma[1],N))
ts2 <- microbenchmark(chain = rw.Metropolis(x0,sigma[2],N),
                     chain_c = rw_Metropolis_c(x0,sigma[2],N))
ts3 <- microbenchmark(chain = rw.Metropolis(x0,sigma[3],N),
                     chain_c = rw_Metropolis_c(x0,sigma[3],N))
ts4 <- microbenchmark(chain = rw.Metropolis(x0,sigma[4],N),
                     chain_c = rw_Metropolis_c(x0,sigma[4],N))
for(i in 1:length(sigma)){
  b <- cbind(b,summary(get(paste0("ts",i)))$median)
  colnames(b)[i+1] <- paste0("sigma = ",sigma[i])
}
b <- b[2:5]
rownames(b) <- c("from R","from Rcpp")
knitr::kable(b)

