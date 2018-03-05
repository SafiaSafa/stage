PC_0=read.table("file:///home/vcabeli/stage/curve_analysis/result_1/0/resultat/prs_low_ld_5.res",sep=",")
PC_2=read.table("file:///home/vcabeli/stage/curve_analysis/result_1/2/resultat/prs_low_ld_5.res",sep=",")
PC_10=read.table("file:///home/vcabeli/stage/curve_analysis/result_1/10/resultat/prs_low_ld_5.res",sep=",")
PC_20=read.table("file:///home/vcabeli/stage/curve_analysis/result_1/20/resultat/prs_low_ld_5.res",sep=",")
PC_30=read.table("file:///home/vcabeli/stage/curve_analysis/result_1/30/resultat/prs_low_ld_5.res",sep=",")
PC_40=read.table("file:///home/vcabeli/stage/curve_analysis/result_1/40/resultat/prs_low_ld_5.res",sep=",")
PC_50=read.table("file:///home/vcabeli/stage/curve_analysis/result_1/50/resultat/prs_low_ld_5.res",sep=",")
PC_100=read.table("file:///home/vcabeli/stage/curve_analysis/result_1/100/resultat/prs_low_ld_5.res",sep=",")

SNP_10000=c(PC_0[10],PC_2[10],PC_10[10],PC_20[10],PC_30[10],PC_40[10],PC_50[10],PC_100[10])
vec_10000=as.data.frame(SNP_10000)
colnames(vec_10000)=c('PC0','PC2','PC10','PC20','PC30','PC40','PC50','PC100')
POINT=lapply(vec_10000,mean)

SNP_20000=c(PC_0[20],PC_2[20],PC_10[20],PC_20[20],PC_30[20],PC_40[20],PC_50[20],PC_100[20])
vec_20000=as.data.frame(SNP_20000)
colnames(vec_20000)=c('PC0','PC2','PC10','PC20','PC30','PC40','PC50','PC100')
POINT_2=lapply(vec_20000,mean)

SNP_5000=c(PC_0[5],PC_2[5],PC_10[5],PC_20[5],PC_30[5],PC_40[5],PC_50[5],PC_100[5])
vec_5000=as.data.frame(SNP_5000)
colnames(vec_5000)=c('PC0','PC2','PC10','PC20','PC30','PC40','PC50','PC100')
POINT_3=lapply(vec_5000,mean)

plot(seq(1, length(PC_0)*1000, by=1000),lapply(PC_0,mean),type="l",xlab="# SNPs",ylab="AUC", main="AUC en fonction du # SNPs (SNPs classé par P-value)")
lines(seq(1, length(PC_0)*1000, by=1000),lapply(PC_2,mean),col="red")
lines(seq(1, length(PC_0)*1000, by=1000),lapply(PC_10,mean),col="green")
lines(seq(1, length(PC_0)*1000, by=1000),lapply(PC_20,mean),col="blue")
lines(seq(1, length(PC_0)*1000, by=1000),lapply(PC_30,mean),col="brown")
lines(seq(1, length(PC_0)*1000, by=1000),lapply(PC_40,mean),col="pink")
lines(seq(1, length(PC_0)*1000, by=1000),lapply(PC_50,mean),col="yellow")
lines(seq(1, length(PC_0)*1000, by=1000),lapply(PC_100,mean),col="grey")
legend("bottomright",legend=c("0 PCs","2 PCs","10 PCs","20 PCs","30 PCs","40 PCs","50 PCS","100 PCs"),col=c("black",'red','green','blue',"brown","pink","yellow","grey"),lty=1,lwd = 2)


plot(c(1,1000),lapply(PC_0[1:2],mean),type="l",xlab="# SNPs",ylim=c(0.5,0.55),ylab="AUC", main="AUC en fonction du # SNPs (SNPs classé par P-value)\n zoom sur les 1000 premiers SNPs")
lines(c(1,1000),lapply(PC_2[1:2],mean),col="red")
lines(c(1,1000),lapply(PC_10[1:2],mean),col="green")
lines(c(1,1000),lapply(PC_20[1:2],mean),col="blue")
lines(c(1,1000),lapply(PC_30[1:2],mean),col="brown")
lines(c(1,1000),lapply(PC_40[1:2],mean),col="pink")
lines(c(1,1000),lapply(PC_50[1:2],mean),col="yellow")
lines(c(1,1000),lapply(PC_100[1:2],mean),col="grey")
legend("topleft",legend=c("0 PCs","2 PCs","10 PCs","20 PCs","30 PCs","40 PCs","50 PCS","100 PCs"),col=c("black",'red','green','blue',"brown","pink","yellow","grey"),lty=1,lwd = 2)

x=c(0,2,10,20,30,40,50,100)
plot(x,POINT, type="l", ylim=c(0.52,0.58),xaxt="n",xlab="# PCs",ylab="AUC", main="AUC en fonction du # PCs suivant un nombre de SNPs")
axis(1, at=x)
lines(x,POINT_2,col="red")
lines(x,POINT_3,col="green")
legend("topright",legend=c("5000 SNPs","10000 SNPs","20000 SNPs"),col=c("green",'black','red'),lty=1,lwd = 2)





