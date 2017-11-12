# Computational-Program-Codes-Developed-in-R-and-SAS-for-GOLLTL-G-paper###

#######shape of density Weibull#######

OLLTLWd<-function(x,alpha,teta,a,b){
 G<-rep(0,0)
 for(i in 1:length(x)){
    G[i]<-(2*alpha*teta*dweibull(x[i], shape=a, scale = (1/b))
    *(1-pweibull(x[i], shape=a, scale = (1/b)))*((1-(1-pweibull(x[i], 
    shape=a, scale = (1/b)))^2)^((alpha*teta)-1))*(1-((1-(1-pweibull(x[i], shape=a, 
    scale = (1/b)))^2))^teta)^(alpha-1))/
((((1-(1-pweibull(x[i], shape=a, scale = (1/b)))^2)^(alpha*teta))
+((1-((1-(1-pweibull(x[i], shape=a, scale = (1/b)))^2))^teta)^(alpha)))^2)
 }
return(G)
}

#####P(X>3)#######

OLLTLWF<-function(x,alpha,teta,a,b){
 G<-rep(0,0)
 F<-rep(0,0)
    G<-function(x){
(2*alpha*teta*dweibull(x, shape=a, scale = (1/b))
*(1-pweibull(x, shape=a, scale = (1/b)))
*((1-(1-pweibull(x, shape=a, scale = (1/b)))^2)^((alpha*teta)-1))
*(1-((1-(1-pweibull(x, shape=a, scale = (1/b)))^2))^teta)^(alpha-1))/
((((1-(1-pweibull(x, shape=a, scale = (1/b)))^2)^(alpha*teta))
+((1-((1-(1-pweibull(x, shape=a, scale = (1/b)))^2))^teta)^(alpha)))^2)
}
F<-integrate(G,lower=-Inf,upper=x)$value
return(F)
}

#######Heavy tail of Weibull#######

##a=1,b=4
1-OLLTLWF(3,.2,5,1,4)

1-pweibull(3,shape=1, scale = (1/4))
##a=1, b=2
1-OLLTLWF(3,.3,7.5,1,2)

1-pweibull(3,shape=1, scale = (1/2))
##a=1.5, b=1
1-OLLTLWF(3,.4,3,1.5,1)

1-pweibull(3,shape=1.5, scale = (1/1))
##
##a=1, b=1
1-OLLTLWF(3,.4,3,1,1)

1-pweibull(3,shape=1, scale = (1/1))
##

##########################################################

OLLTLNd<-function(x,alpha,teta,mu,sigma){
 G<-rep(0,0)
 for(i in 1:length(x)){
    G[i]<-(2*alpha*teta*dnorm(x[i],mu,sigma)*(1-pnorm(x[i],mu,sigma))
    *((1-(1-pnorm(x[i],mu,sigma))^2)^((alpha*teta)-1))
    *(1-((1-(1-pnorm(x[i],mu,sigma))^2))^teta)^(alpha-1))/
((((1-(1-pnorm(x[i],mu,sigma))^2)^(alpha*teta))
+((1-((1-(1-pnorm(x[i],mu,sigma))^2))^teta)^(alpha)))^2)
 }
return(G)
}
OLLTLNF<-function(x,alpha,teta,mu,sigma){
 G<-rep(0,0)
 F<-rep(0,0)
    G<-function(x){
(2*alpha*teta*dnorm(x,mu,sigma)*(1-pnorm(x,mu,sigma))
*((1-(1-pnorm(x,mu,sigma))^2)^((alpha*teta)-1))
*(1-((1-(1-pnorm(x,mu,sigma))^2))^teta)^(alpha-1))/
((((1-(1-pnorm(x,mu,sigma))^2)^(alpha*teta))
+((1-((1-(1-pnorm(x,mu,sigma))^2))^teta)^(alpha)))^2)
}
F<-integrate(G,lower=-Inf,upper=x)$value
return(F)
}

#######Right tail of normal#######

##mu=0,sigma=1
1-OLLTLNF(3,.2,5,0,1)
1-pnorm(3,0,1)
######Left tail of normal
OLLTLNF(-3,.2,5,0,1)
pnorm(-3,0,1)
#########Right tail of normal
##mu=0,sigma=2
1-OLLTLNF(3,.35,3,0,2)
1-pnorm(3,0,2)
######Left tail of normal
OLLTLNF(-3,.35,3,0,2)
pnorm(-3,0,2)

#######asymptotic distribution of estimator#######

rOLLTLW<-function(n,alpha,teta,a,b){
G<-rep(0,0)
u<-runif(n,0,1)
G<-qweibull((1-(1-((u^(1/(alpha*teta)))/((((u^(1/alpha))+
((1-u)^(1/alpha))))^(1/teta))))^(1/2))
,shape=b, scale = (1/a))
return(G)
}
#####
OLLTLWd<-function(x,alpha,teta,a,b){
 G<-rep(0,0)
 for(i in 1:length(x)){
G[i]<-(2*alpha*teta*dweibull(x[i],shape=b, scale = (1/a))*(1-pweibull(x[i],shape=b,
scale = (1/a)))*((1-(1-pweibull(x[i],shape=b, 
scale = (1/a)))^2)^((alpha*teta)-1))*(1-((1-(1-pweibull(x[i],shape=b,
scale = (1/a)))^2))^teta)^(alpha-1))/((((1-(1-pweibull(x[i],shape=b,
scale = (1/a)))^2)^(alpha*teta))
+((1-((1-(1-pweibull(x[i],shape=b, scale = (1/a)))^2))^teta)^(alpha)))^2)
 }
return(G)
}
######
OLLTLWlikelihood<-function(par){
alpha<-par[1]
teta<-par[2]
a<-par[3]
b<-par[4]       
    G<- -sum(log(OLLTLWd(x,alpha,teta,a,b)))
return(G)
}
#####
M<-1
MLestimate<-matrix(c(0),ncol=4,nrow=1000)
###real value
alpha<-0.5
teta<-0.5
a<-2
b<-2
size<-1500
for (i in 1:1000){
estimate<-matrix(c(0),ncol=4,nrow=10000)
count1<-0
repeat {
    if (length(estimate[,1][estimate[,1]!=0])==M) break
    count1 <- count1 + 1
    #print(count)
    #print(length(estimate[,1][estimate[,1]!=0]))
    x<-rOLLTLW(size,alpha,teta,a,b)
OLLTLWest <- try(optim(c(alpha,teta,a,b), OLLTLWlikelihood,method="L-BFGS-B",
hessian=TRUE,lower = c(.01,.01,.01,.01), upper =c(Inf,Inf,Inf,Inf)), silent=TRUE)
if ('try-error' %in% class(OLLTLWest)) next
else{ 
estimate[count1,]<-OLLTLWest$par
}
}
MLestimate[i,1]<-estimate[,1][estimate[,1]!=0]
MLestimate[i,2]<-estimate[,2][estimate[,2]!=0]
MLestimate[i,3]<-estimate[,3][estimate[,3]!=0]
MLestimate[i,4]<-estimate[,4][estimate[,4]!=0]
print(i)
}
par(mfrow=c(2,2))
qqnorm(MLestimate[1:1000,1]); qqline(MLestimate[1:1000,1], col = 2)
qqnorm(MLestimate[1:1000,2]); qqline(MLestimate[1:1000,2], col = 2)
qqnorm(MLestimate[1:1000,3]); qqline(MLestimate[1:1000,3], col = 2)
qqnorm(MLestimate[1:1000,4]); qqline(MLestimate[1:1000,4], col = 2)


#######Graphical Simulation Codes#######

library(AdequacyModel)
cdfgolltlweibull=function(par,x)
{
alpha=par[1]
theta=par[2]
a=par[3]
b=par[4]
G=pweibull(x,a,b)
g=dweibull(x,a,b)

f=(1-(1-G)^2)^(alpha*theta)/((1-(1-G)^2)^(alpha*theta)
+(1-(1-(1-G)^2)^(theta))^alpha)

return(f)
}

pdfgolltlweibull=function(par,x)
{
alpha=par[1]
theta=par[2]
a=par[3]
b=par[4]
G=pweibull(x,a,b)
g=dweibull(x,a,b)

f=2*alpha*theta*(1-G)*(1-(1-G)^2)^(alpha*theta-1)
*(1-(1-(1-G)^2)^(theta))^(alpha-1)/(((1-(1-G)^2)^(alpha*theta)
+(1-(1-(1-G)^2)^(theta))^alpha))^2


return(f)
}

qq=function(alpha,theta,a,b,p) {

qweibull(1-(1-p^(1/(alpha*theta))
/(p^(1/alpha)+(1-p)^(1/alpha))^(1/theta))^(1/2),mu,sigma)

}

for( k in 1:191)
{
for( j in 1:1000)
{

for(i in 1:(45+k*5)) {

p=runif(1,0,1)

alpha1=0.5
theta1=0.5
a1=2
b1=1

datt[i]=qq(alpha1,theta1,a1,b1,p)
}

fgolltl=goodness.fit(pdf=pdfnolln, cdf=cdfnolln,
starts = c(alpha1,theta1,a1,b1), data = datt,
method="N", domain=c(0,Inf))

alpha[j]=fgolltl$mle[1]
theta[j]=fgolltl$mle[2]
a[j]=fgolltl$mle[3]
b[j]=fgolltl$mle[4]

salpha[j]=fgolltl$Erro[1]
stheta[j]=fgolltl$Erro[2]
sa[j]=fgolltl$Erro[3]
sb[j]=fgolltl$Erro[4]

lalpha[j]=alpha[j]-1.96*salpha[j]
ualpha[j]=alpha[j]+1.96*salpha[j]
ltheta[j]=theta[j]-1.96*stheta[j]
utheta[j]=theta[j]+1.96*stheta[j]
la[j]=a[j]-1.96*sa[j]
ua[j]=a[j]+1.96*sa[j]

lb[j]=b[j]-1.96*sb[j]
ub[j]=b[j]+1.96*sb[j]

}

ALlambda[k]=(3.92/NROW((salpha)))*sum((salpha))
ALbeta[k]=(3.92/NROW((stheta)))*sum((stheta))
ALa[k]=(3.92/NROW((sa)))*sum((sa))
ALb[k]=(3.92/NROW((sb)))*sum((sb))

cplambda[k]=mean((lalpha) < alpha1 & (ualpha) > alpha1);
cpbeta[k]=mean((ltheta) < theta1 & (utheta) > theta1);
cpa[k]=mean((la) < a1 & (ua) >a1);
cpb[k]=mean((lb) < b1 & (ub) >b1);

biaslambda[k]=mean(alpha-alpha1)
biasbeta[k]=mean(theta-theta1)
biasa[k]=mean(a-a1)
biasb[k]=mean(b-b1)


mselambda[k]=sum((alpha-alpha1)^2)/1000
msebeta[k]=sum((theta-theta1)^2)/1000
msea[k]=sum((a-a1)^2)/1000
mseb[k]=sum((b-b1)^2)/1000

}

#######Codes of Regression Models developed in SAS####### 

#######data paper martinez 2013#######

rm(list=ls(all=TRUE))
library(survival)

#######Loading data set #######

data<- read.csv("https://goo.gl/Zj78Zb",sep=";")
attach(dados)

#######Defining dummis variables#######

x2<-ifelse(x1=='alone',0,1)

#logit function
ll<-function(x) exp(x)/(exp(x)+1);  ll2<-function(p) log(p/(1-p))


#######cdf OLLTL-G#######

GW<-function(x,alpha,theta,mu,sigma){
  ((1-(1-G(x,mu,sigma))^2)^(alpha*theta))/((1-(1-G(x,mu,sigma))^2)^(alpha*theta) 
  +(1-(1-(1-G(x,mu,sigma))^2)^theta)^alpha)}

#######pdf OLLTL-G#######

gw<-function(x,alpha,theta,mu,sigma){
  (2*alpha*theta*g(x,mu,sigma)*(1-G(x,mu,sigma))
  *(1-(1-G(x,mu,sigma))^2)^(alpha*theta-1)
  *(1-(1-(1-G(x,mu,sigma))^2)^theta)^(alpha-1))/
    (((1-(1-G(x,mu,sigma))^2)^(alpha*theta) 
    + (1-(1-(1-G(x,mu,sigma))^2)^theta)^alpha  )^2)}


#######Cure rate family cdf#######

GWp<-function(x,alpha,theta,mu,sigma,p){
  (1-p)*(  ((1-(1-G(x,mu,sigma))^2)^(alpha*theta))
  /((1-(1-G(x,mu,sigma))^2)^(alpha*theta) 
  +(1-(1-(1-G(x,mu,sigma))^2)^theta)^alpha))}

#######Cure rate family pdf#######

gwp<-function(x,alpha,theta,mu,sigma,p){
  (1-p)*(  (2*alpha*theta*g(x,mu,sigma)
  *(1-G(x,mu,sigma))
  *(1-(1-G(x,mu,sigma))^2)^(alpha*theta-1)
  *(1-(1-(1-G(x,mu,sigma))^2)^theta)^(alpha-1))
  /(((1-(1-G(x,mu,sigma))^2)^(alpha*theta) 
             + (1-(1-(1-G(x,mu,sigma))^2)^theta)^alpha  )^2))}


#######Log Weibull as basis distribution#######

g<-function(x,mu,sigma){exp(x)*dweibull(exp(x),shape=1/sigma,scale=exp(mu))}
G<-function(x,mu,sigma){pweibull(exp(x),shape =1/sigma,scale=exp(mu))}


km <- survfit(Surv(log(time),censur)~x2)
plot(km,lwd=2,conf.int=F,ylab='Survival',xlab='Time',col=c(1))

#######MLEs obtained using SAS (LOLLTL-W location model)#######

curve(1-GWp(x,0.728,    0.274, 3.036,
exp(-1.4455),ll(-0.1115)),-2,4,add=T,col='royalblue3',lwd=4)

curve(1-GWp(x,0.728,    0.274, 3.036+0.296,
exp(-1.4455),ll(-0.1115+0.291)),-2,4,add=T,col='orangered3',lwd=4)

legend("bottomleft",c('Kaplan-Meier','x1=0',"x1=1"),
lwd=c(4,4),col=c(1,'royalblue3','orangered3'),
bty="n",seg.len=2,cex=1.5)

plot(km,lwd=2,conf.int=F,ylab='Survival',xlab='Time',col=c(1))

#######MLEs obtained using SAS (LOLLTL-W heteroscedastic model)#######

curve(1-GWp(x,0.297,    1.537, 2.279,
exp(-0.664),ll(-0.0986)),-2,4,add=T,col='royalblue1',lwd=4)

curve(1-GWp(x,0.297,    1.537, 2.279
+0.803,exp(-0.664-1.102),ll(-0.0986+0.263)),-2,4,add=T,col='orangered1',lwd=4)

legend("bottomleft",c('Kaplan-Meier','x1=0',"x1=1"),
lwd=c(4,4),col=c(1,'royalblue1','orangered1'),
bty="n",seg.len=2,cex=1.5)
