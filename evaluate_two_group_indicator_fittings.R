
library(HiddenMarkov)

mu0_gibbs <- read.table("mu0_gibbs.txt")
sigma20_gibbs <- read.table("sigma20_gibbs.txt")
P_gibbs <- read.table("P_gibbs.txt")

pdf("mu0sigma20P.pdf")

par(mfrow=c(3,2))
for(k in 1:3){
    plot(mu0_gibbs[,k],type="l",main=paste("mu0",k))
    plot(sigma20_gibbs[,k],type="l",main=paste("sigma20",k))
    acf(mu0_gibbs[,k])
    acf(sigma20_gibbs[,k])
    hist(mu0_gibbs[,k],40,col="grey")
    hist(sigma20_gibbs[,k],40,col="grey")
}


t <- 1
for(k in 1:3){
for(j in 1:3){
if(j != k){
par(mfrow=c(2,2))
    plot(P_gibbs[,t],type="l",main=paste("P",k,j))
    acf(P_gibbs[,t])
    hist(P_gibbs[,t])
    t <- t + 1
}
}
}
# check detailed balance
Ag <- array(0,c(dim(P_gibbs)[1],3,3))
for(i in 1:dim(P_gibbs)[1]){
    tt <- 1
    for(j in 1:3){
    for(k in 1:3){
    if(k != j){
    Ag[i,j,k] <- P_gibbs[i,tt]
    Ag[i,j,j] <- Ag[i,j,j] - P_gibbs[i,tt]
    tt <- tt + 1
    }
    }
    Ag[i,j,j] <- Ag[i,j,j] + 1
    }
}
piA <- array(0,c(dim(P_gibbs)[1],3,3))
Api <- array(0,c(dim(P_gibbs)[1],3,3))
for(tt in 1:dim(P_gibbs)[1]){
    piA[tt,,]=diag(compdelta(Ag[tt,,]))%*%Ag[tt,,]
    Api[tt,,]=t(Ag[tt,,])%*%diag(compdelta(Ag[tt,,]))
}
par(mfrow=c(3,3));
for(ii in 1:3){
for(jj in 1:3){
    plot(piA[,ii,jj],Api[,ii,jj],main=paste(ii,"<~~>",jj),xlab=paste("Pi",ii,"P(",ii,jj,")"),ylab=paste("P(",jj,ii,")Pi",jj));
    lines(c(min(piA[,ii,jj],Api[,ii,jj]),max(piA[,ii,jj],Api[,ii,jj])),c(min(piA[,ii,jj],Api[,ii,jj]),max(piA[,ii,jj],Api[,ii,jj])),col="red");
}
}
for(ii in 1:2){
for(jj in 2:3){
if(jj!=ii){
   x=hist(piA[,ii,jj],col="grey",xlab=paste("Pi",ii,"P(",ii,jj,")"),main=paste(ii,"to",jj),20);
   xx=hist(Api[,ii,jj],col="grey",xlab=paste("Pi",jj,"P(",jj,ii,")"),main=paste(jj,"to",ii),20);
   hist(piA[,ii,jj]-Api[,ii,jj],col="grey",xlab=paste("Pi",ii,"P(",ii,jj,")-Pi",jj,"P(",jj,ii,")"),main="difference",20,xlim=c(-(range(x$mids)[2]-range(x$mids)[1])/2,(range(x$mids)[2]-range(x$mids)[1])/2));
   abline(v=0,col="red",lwd=3);
}
}
}
dev.off()


indicator_gibbs <- read.table("state_indicators.txt")

T <- dim(indicator_gibbs)[2]/4
I <- array(0,c(dim(indicator_gibbs)[1],T))
for(n in 1:(dim(indicator_gibbs)[1])){
    for(t in 1:T){
        if(indicator_gibbs[n,4 * (t-1) + 2]==3) {I[n,t] <- 1}
        if(indicator_gibbs[n,4 * (t-1) + 2]==2){
           if(indicator_gibbs[n,4 * (t-1) + 3]==1 & indicator_gibbs[n,4 * (t-1) + 4]==2){I[n,t]<- 2}
           if(indicator_gibbs[n,4 * (t-1) + 3]==1 & indicator_gibbs[n,4 * (t-1) + 4]==3){I[n,t]<- 3}
           if(indicator_gibbs[n,4 * (t-1) + 3]==2 & indicator_gibbs[n,4 * (t-1) + 4]==3){I[n,t]<- 4}
        }
    }
}
prob <- array(0,c(T,4))
select <- rep(0,T)
for(t in 1:T){
    prob[t,] <- hist(I[,t],breaks=c(0.5,1.5,2.5,3.5,4.5))$intensities
    if(which(prob[t,]>0.5)==1){select[t] <- 3}
    if(which(prob[t,]>0.5)>1) {select[t] <- 2}
}
write.table(select,"model_select_hier.txt",row.names=FALSE,col.names=FALSE)
pdf("state_indicators.pdf")
par(mfrow=c(4,4))
for(t in 1:T){
    hist(I[,t],breaks=c(0.5,1.5,2.5,3.5,4.5),main=paste("Trace", t),xlab = "state indicator")
}
dev.off()


log_post <- read.table("logposterior.txt")[,1]
log_post <- as.numeric(log_post)
burnin <- 1
pdf("log_post_burn1.pdf")
plot(burnin:length(log_post),log_post[burnin:length(log_post)],type = "l")
dev.off()
burnin <- 2000
pdf("log_post_burn2000.pdf")
plot(burnin:length(log_post),log_post[burnin:length(log_post)],type = "l")
dev.off()


# read data Htraces

load("Htraces.RData")
thin <- 10
burnin <- 2000
num_gibbs <- 2000

for(i in 1:T){

pdf(paste("Trace",i,"fittings.pdf"))

for(num in 1: (num_gibbs/100)){

    current <- (num-1) * thin * 100 + burnin + 1
    name2 <- paste("2_S_paras",current,".txt",sep="")
    S2 <- read.table(name2)
    name3 <- paste("3_S_paras",current,".txt",sep="")
    S3 <- read.table(name3)
    which_two <- S2[,2:3]
    idx2 <- S2[,1] + 1
    idx3 <- S3[,1] + 1

    P <- read.table(paste("global_paras",current,".txt",sep=""))[3:5,]
    P <- array(as.numeric(as.matrix(P)),c(3,3))


        if(sum(i==idx2)==0){
           temp <- which(i==idx3)
           mu <- as.numeric(S3[temp,2:4])
           sig2 <- as.numeric(S3[temp,5:7])
           pii <- as.numeric(S3[temp,8:10])
           P_temp <- P
           pii[which(pii<0)] <- 0
        }
        if(sum(i==idx3)==0){
           temp <- which(i==idx2)
           mu <- as.numeric(S2[temp,4:5])
           sig2 <- as.numeric(S2[temp,6:7])
           pii <- as.numeric(S2[temp,8:9])
           pii[which(pii<0)] <- 0
           P_temp <- array(0,c(2,2))
           low <- which_two[temp,1]
           high <- which_two[temp,2]
           P_temp[1,1] <- P[low,low]/(P[low,low]+P[low,high])
           P_temp[1,2] <- P[low,high]/(P[low,low]+P[low,high])
           P_temp[2,1] <- P[high,low]/(P[high,low]+P[high,high])
           P_temp[2,2] <- P[high,high]/(P[high,low]+P[high,high])
        } 
        z <- list(x=Htraces[[i]]$F, Pi=P_temp, delta=pii, distn = "norm", pm =list(mean=mu,sd=sqrt(sig2)) ,pn = NULL, discrete = NULL, nonstat = TRUE)
        class(z) <- "dthmm"
        y <- Viterbi(z)
        if(length(mu)==3)
           p <- c(sum(y==1),sum(y==2),sum(y==3))/length(y)
        if(length(mu)==2)
           p <- c(sum(y==1),sum(y==2))/length(y)
#       par(mfrow=c(2,1))
#       plot(Htraces[[i]]$T,Htraces[[i]]$D,type='l',ylab="Donor Acceptor",xlab="time",main="Donor(green) and Acceptor(red)",col="green",ylim=c(min(Htraces[[i]]$D,Htraces[[i]]$A)-20,max(Htraces[[i]]$D,Htraces[[i]]$A)+20))
#       lines(Htraces[[i]]$T,Htraces[[i]]$A,type='l',xlab="time",col="red")
#       plot(Htraces[[i]]$T,Htraces[[i]]$F,type='l',xlab="time",ylab="fret",main=paste("Fret trace ",i),col="blue")
#       lines(Htraces[[i]]$T,mu[y],col="black")
       par(mfrow=c(1,1))
       histind <- hist(Htraces[[i]]$F,60,freq= FALSE,main = paste("Histogram of Trajectory" ,i),xlab ="FRET",col=240)
       if(length(mu)==3){
       func <- function(x) dnorm(x,mu[1],sqrt(sig2[1]))*p[1] + dnorm(x,mu[2],sqrt(sig2[2]))*p[2] + dnorm(x,mu[3],sqrt(sig2[3]))*p[3]
       lines(histind$mids, func(histind$mids),col="red",lwd=2)
       lines(histind$mids, p[1]*dnorm(histind$mids,mu[1],sqrt(sig2[1])),col="yellow")
       lines(histind$mids, p[2]*dnorm(histind$mids,mu[2],sqrt(sig2[2])),col="yellow")
       lines(histind$mids, p[3]*dnorm(histind$mids,mu[3],sqrt(sig2[3])),col="yellow")
       }
       if(length(mu)==2){
       func <- function(x) dnorm(x,mu[1],sqrt(sig2[1]))*p[1] + dnorm(x,mu[2],sqrt(sig2[2]))*p[2]
       lines(histind$mids, func(histind$mids),col="red",lwd=2)
       lines(histind$mids, p[1]*dnorm(histind$mids,mu[1],sqrt(sig2[1])),col="yellow")
       lines(histind$mids, p[2]*dnorm(histind$mids,mu[2],sqrt(sig2[2])),col="yellow")
       }
}
for(num in 1: (num_gibbs/100)){

    current <- (num-1) * thin * 100 + burnin + 1
    name2 <- paste("2_S_paras",current,".txt",sep="")
    S2 <- read.table(name2)
    name3 <- paste("3_S_paras",current,".txt",sep="")
    S3 <- read.table(name3)
    which_two <- S2[,2:3]
    idx2 <- S2[,1] + 1
    idx3 <- S3[,1] + 1

    P <- read.table(paste("global_paras",current,".txt",sep=""))[3:5,]
    P <- array(as.numeric(as.matrix(P)),c(3,3))


        if(sum(i==idx2)==0){
           temp <- which(i==idx3)
           mu <- as.numeric(S3[temp,2:4])
           sig2 <- as.numeric(S3[temp,5:7])
           pii <- as.numeric(S3[temp,8:10])
           P_temp <- P
           pii[which(pii<0)] <- 0
        }
        if(sum(i==idx3)==0){
           temp <- which(i==idx2)
           mu <- as.numeric(S2[temp,4:5])
           sig2 <- as.numeric(S2[temp,6:7])
           pii <- as.numeric(S2[temp,8:9])
           pii[which(pii<0)] <- 0
           P_temp <- array(0,c(2,2))
           low <- which_two[temp,1]
           high <- which_two[temp,2]
           P_temp[1,1] <- P[low,low]/(P[low,low]+P[low,high])
           P_temp[1,2] <- P[low,high]/(P[low,low]+P[low,high])
           P_temp[2,1] <- P[high,low]/(P[high,low]+P[high,high])
           P_temp[2,2] <- P[high,high]/(P[high,low]+P[high,high])
        } 
        z <- list(x=Htraces[[i]]$F, Pi=P_temp, delta=pii, distn = "norm", pm =list(mean=mu,sd=sqrt(sig2)) ,pn = NULL, discrete = NULL, nonstat = TRUE)
        class(z) <- "dthmm"
        y <- Viterbi(z)
        if(length(mu)==3)
           p <- c(sum(y==1),sum(y==2),sum(y==3))/length(y)
        if(length(mu)==2)
           p <- c(sum(y==1),sum(y==2))/length(y)
       par(mfrow=c(2,1))
       plot(Htraces[[i]]$T,Htraces[[i]]$D,type='l',ylab="Donor Acceptor",xlab="time",main="Donor(green) and Acceptor(red)",col="green",ylim=c(min(Htraces[[i]]$D,Htraces[[i]]$A)-20,max(Htraces[[i]]$D,Htraces[[i]]$A)+20))
       lines(Htraces[[i]]$T,Htraces[[i]]$A,type='l',xlab="time",col="red")
       plot(Htraces[[i]]$T,Htraces[[i]]$F,type='l',xlab="time",ylab="fret",main=paste("Fret trace ",i),col="blue")
       lines(Htraces[[i]]$T,mu[y],col="black")
}
dev.off()
}







