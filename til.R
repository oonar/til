#Different coefficient for studying dependence
#ormrai 2022-07-21
#til.R begins

#the necessary libraries (install if not already existing!)
library(acepack)
library(infotheo)
library(minerva)
library(energy)

#function for generating data
simxy1<-function(j,n,sigma){
  x<-rnorm(n,0,1) #normal distribution
  #x<-runif(n,-3,3) #uniform distribution
  #x<-rexp(n,0.5) #exponential distribution
  #x<-rpois(n,3) #poisson distribution
  #Linear dependence
  if(j==1){
    y<-x+rnorm(n,0,sigma)
  }
  #Logarithmic dependence
  if(j==2){
    y<-5*log(abs(5+x))+rnorm(n,0,sigma)
  }
  #Cubic dependence
  if(j==3){
    y<-0.3*x^3+rnorm(n,0,sigma)
  }
  #Quadratic dependence
  if(j==4){
    y<-0.7*x^2+rnorm(n,0,sigma)
  }
  #Sinusoidal dependence
  if(j==5){
    y<-1.3*sin(3*x)+rnorm(n,0,sigma)
  }
  #Rational function
  if(j==6){
    for(i in 1:n){
      y[i]<-min(max(1/x[i],-3),3)+rnorm(1,0,sigma)
    }
  }
  #Cross-shaped dependence
  if(j==7){
    x<-c(rnorm(floor(n/2),0,sigma/3),rnorm(n-floor(n/2),0,1))
    y<-c(rnorm(floor(n/2),0,1),rnorm(n-floor(n/2),0,sigma/3))
  }
  #Circular dependence
  if(j==8){
    k<-rnorm(n,0,1)
    h<-rnorm(n,1,sigma/7)
    x<-h*cos(k)
    y<-h*sin(k)
  }
  #Checkerboard dependence
  if(j==9){
    x<-rnorm(n,0,1)
    y<-rnorm(n,0,1)
    for(u in 1:n){
      while(((floor(0.7*x[u])-floor(0.7*y[u]))%%2)==1){
        x[u]<-rnorm(1,0,1)
        y[u]<-rnorm(1,0,1)
      }
    }
    y<-y+rnorm(n,0,sigma/2)
  }
  dxy<-cbind(x,y)
  return(dxy)
}

#studying generality
j<-1 #type of dependence
k<-1000 #number of simulations
n<-30 #number of observations in each simulation
sigma<-1 #noise parameter
par(pty='s',mfrow=c(1,1))
plot(simxy1(j,n,sigma),xaxt='n',yaxt='n',xlab='',ylab='')
df<-as.data.frame(matrix(NA,k,6))
for(i in 1:k){
  dxy<-simxy1(j,n,sigma)
  x<-dxy[,1]
  y<-dxy[,2]
  df[i,1]<-abs(cor(x,y))
  df[i,2]<-abs(cor(x,y,method='spearman'))
  fxy<-ace(x,y)
  df[i,3]<-cor(fxy$tx,fxy$ty)
  df[i,4]<-dcor(x,y,index=1)
  disc<-discretize(data.frame(x,y))
  mi<-mutinformation(disc$x,disc$y)
  df[i,5]<-sqrt(1-exp(-2*mi))
  df[i,6]<-mine(x,y)$MIC
}
colnames(df)<-c('|r|','|r_s|','max','dist','r_1','MIC')
#average values of coefficients for fixed parameters
print(colMeans(df))

#function for generating data from independent variables
simxy0<-function(n){
  #Independent normally distributed variables
  x<-rnorm(n,0,1)
  y<-rnorm(n,0,1)
  dxy<-cbind(x,y)
  return(dxy)
}
k<-1000
n<-30
alpha<-0.01 #level of significance
df<-as.data.frame(matrix(NA,k,6))
for(i in 1:k){
  dxy<-simxy0(n)
  x<-dxy[,1]
  y<-dxy[,2]
  df[i,1]<-abs(cor(x,y))
  df[i,2]<-abs(cor(x,y,method='spearman'))
  fxy<-ace(x,y)
  df[i,3]<-cor(fxy$tx,fxy$ty)
  df[i,4]<-dcor(x,y,index=1)
  disc<-discretize(data.frame(x,y))
  mi<-mutinformation(disc$x,disc$y)
  df[i,5]<-sqrt(1-exp(-2*mi))
  df[i,6]<-mine(x,y)$MIC
}
critVs<-c() #critical values
for(i in 1:6){
  critVs[i]<-as.numeric(quantile(df[,i],probs=1-alpha))
}

#studying power
k<-1000
n<-30
j<-3
sig1<-seq(0,10,by=1)
df<-as.data.frame(matrix(NA,length(sig1),6))
for(s in 1:length(sig1)){
  dfs<-as.data.frame(matrix(NA,k,6))
  sigma<-sig1[s]
  for(i in 1:k){
    dxy<-simxy1(j,n,sigma)
    x<-dxy[,1]
    y<-dxy[,2]
    dfs[i,1]<-abs(cor(x,y))
    dfs[i,2]<-abs(cor(x,y,method='spearman'))
    fxy<-ace(x,y)
    dfs[i,3]<-cor(fxy$tx,fxy$ty)
    dfs[i,4]<-dcor(x,y,index=1)
    disc<-discretize(data.frame(x,y))
    mi<-mutinformation(disc$x,disc$y)
    dfs[i,5]<-sqrt(1-exp(-2*mi))
    dfs[i,6]<-mine(x,y)$MIC
  }
  for(i in 1:6){
    df[s,i]<-sum(dfs[,i]>=critVs[i])/k
  }
}
par(pty='m',mfrow=c(1,1))
ltypes<-c(1,2,3,4,5,6)
ptypes<-c(1,4,3,5,6,8)
plot(c(0:(length(sig1)-1)),df[,1],type='l',lty=ltypes[1],
     xlim=c(0,(length(sig1)-1)),ylim=c(0,1),
     ylab='Power',xlab=expression(sigma),
     lwd=1,cex.lab=1.3,cex.axis=1.1,
     main='Power for the cubic dependence')
points(c(0:(length(sig1)-1)),df[,1],pch=ptypes[1])
for(i in 2:6){
  par(new=TRUE)
  plot(c(0:(length(sig1)-1)),df[,i],type='l',
       xlim=c(0,(length(sig1)-1)),ylim=c(0,1),lty=ltypes[i],
       axes=FALSE,xlab='',ylab='')
  points(c(0:(length(sig1)-1)),df[,i],pch=ptypes[i])
} 
legend('topright',legend=c('|r|','|r_s|','rho_max',
                           'rho_dist','r1','MIC'),
       pch=ptypes,lty=ltypes,
       lwd=c(1,1,1,1,1),cex=1.3)

#coefficient of determination from the data
coefOfDet<-function(j,dxy){
  x<-dxy[,1]
  y<-dxy[,2]
  if(j==1){
    k<-cor(x,y)^2
  }
  if(j==2){
    k<-cor(y,5*log(abs(5+x)))^2
  }
  if(j==3){
    k<-cor(y,0.3*x^3)^2
  }
  if(j==4){
    k<-cor(y,0.7*x^2)^2
  }
  if(j==5){
    k<-cor(y,1.3*sin(3*x))^2
  }
  if(j==6){
    fx<-c()
    for(i in 1:length(x)){
      fx[i]<-min(max(1/x[i],-3),3)
    }
    k<-cor(y,fx)^2
  }
  return(k)
}

#quantity (MIC/r1/MCC/DC) from data
qFromdxy<-function(q,dxy){
  x<-dxy[,1]
  y<-dxy[,2]
  if(q==1){
    fxy<-ace(x,y)
    mcc<-cor(fxy$tx,fxy$ty)
    return(mcc)
  }
  if(q==2){
    return(dcor(x,y,index=1))
  }
  if(q==3){
    disc<-discretize(data.frame(x,y))
    mi<-mutinformation(disc$x,disc$y)
    r1<-sqrt(1-exp(-2*mi))
    return(r1)
  }
  if(q==4){
    return(mine(x,y)$MIC)
  }
}

#studying equitability
k<-1000
n<-30
q<-1 #1.MCC,2.DC,3.r1,4.MIC
df<-as.data.frame(matrix(NA,k,12))
for(i in 1:k){
  sigma<-abs(rnorm(1,0,3))
  for(j in 1:6){
    if(j==6){sigma<-1.5*sigma}
    dxy<-simxy1(j,n,sigma)
    df[i,2*j-1]<-coefOfDet(j,dxy)
    df[i,2*j]<-qFromdxy(q,dxy)
  }
}
par(pty='s',mfrow=c(1,1))
if(q<3){q1<-expression(rho)}
if(q==3){q1<-'r1'}
if(q==4){q1<-'MIC'}
plot(1,type='n',xlim=c(0,1),ylim=c(0,1),pch=3,
     ylab=q1,xlab='1-R^2',cex.lab=1.7,cex.axis=1.3)
col1<-c('black','darkblue','blue1','steelblue1','lightsteelblue1','lightgray')
for(i in 1:k){
  for(j in 1:6){
    points(1-df[i,2*j-1],df[i,2*j],pch=3,col=col1[j])
  }
}

#til.R ends
