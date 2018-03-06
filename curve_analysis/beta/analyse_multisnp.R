PC_0=read.table("file:///home/vcabeli/stage/curve_analysis/beta/result_1/0_beta/resultat/prs_low_ld_5.res",sep=",")
PC_1=read.table("file:///home/vcabeli/stage/curve_analysis/beta/result_1/1_beta/resultat/prs_low_ld_5.res",sep=",")
PC_2=read.table("file:///home/vcabeli/stage/curve_analysis/beta/result_1/2_beta/resultat/prs_low_ld_5.res",sep=",")
PC_5=read.table("file:///home/vcabeli/stage/curve_analysis/beta/result_1/5_beta/resultat/prs_low_ld_5.res",sep=",")

SNP_10000=c(PC_0[21],PC_1[21],PC_2[21],PC_5[21])
vec_10000=as.data.frame(SNP_10000)
colnames(vec_10000)=c('PC0','PC1','PC2','PC5')
POINT=lapply(vec_10000,mean)

SNP_20000=c(PC_0[41],PC_1[41],PC_2[41],PC_5[41])
vec_20000=as.data.frame(SNP_20000)
colnames(vec_20000)=c('PC0','PC1','PC2','PC5')
POINT_2=lapply(vec_20000,mean)

SNP_5000=c(PC_0[11],PC_1[11],PC_2[11],PC_5[11])
vec_5000=as.data.frame(SNP_5000)
colnames(vec_5000)=c('PC0','PC1','PC2','PC5')
POINT_3=lapply(vec_5000,mean)

plot(seq(1, length(PC_0)*500, by=500),lapply(PC_0,mean),type="l",xlab="# SNPs",ylab="AUC", main="AUC en fonction du # SNPs (SNPs classé par beta)")
lines(seq(1, length(PC_0)*500, by=500),lapply(PC_1,mean),col="red")
lines(seq(1, length(PC_0)*500, by=500),lapply(PC_2,mean),col="green")
lines(seq(1, length(PC_0)*500, by=500),lapply(PC_5,mean),col="blue")
legend("topright",legend=c("0 PCs","1 PCs","2 PCs","5 PCs"),col=c("black",'red','green','blue'),lty=1,lwd = 2)


plot(c(1,500, 1000),lapply(PC_0[1:3],mean),type="l",xlab="# SNPs",ylab="AUC", main="AUC en fonction du # SNPs (SNPs classé par beta)\n zoom sur les 1000 premiers SNPs")
lines(c(1,500, 1000),lapply(PC_1[1:3],mean),col="red")
lines(c(1,500, 1000),lapply(PC_2[1:3],mean),col="green")
lines(c(1,500, 1000),lapply(PC_5[1:3],mean),col="blue")
legend("topleft",legend=c("0 PCs","1 PCs","2 PCs","5 PCs"),col=c("black",'red','green','blue'),lty=1,lwd = 2)

x=c(0,1,2,5)
plot(x,POINT, type="l", ylim=c(0.5,0.58),xaxt="n",xlab="# PCs",ylab="AUC", main="AUC en fonction du # PCs suivant un nombre de SNPs")
axis(1, at=x)
lines(x,POINT_2,col="red")
lines(x,POINT_3,col="green")
legend("topright",legend=c("5000 SNPs","10000 SNPs","20000 SNPs"),col=c("green",'black','red'),lty=1,lwd = 2)





