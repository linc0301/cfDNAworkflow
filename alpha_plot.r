#alpha_plot.r#

file <- c("BH01","IH01","IH02","CH01")
for (i in file)
{ 
iref <- print(paste0("fragment_",i))
ref <- read.table(paste0(iref,".csv"),header=T,as.is=T)
ref <- ref[1:8,]

setwd(paste0("/users/ludwig/fgm883/fgm883/cftissue/treated/",i,"_processed"))

kmer2 <- read.table("tmp_kmer2.fragpatterns.txt",header=T,as.is=T)
kmer2 <- kmer2[1:8,]
kmer3 <- read.table("tmp_kmer3.fragpatterns.txt",header=T,as.is=T)
kmer3 <- kmer3[1:8,]
kmer4 <- read.table("tmp_kmer4.fragpatterns.txt",header=T,as.is=T)
kmer4 <- kmer4[1:8,]
kmer5 <- read.table("tmp_kmer5.fragpatterns.txt",header=T,as.is=T)
kmer5 <- kmer5[1:8,]
kmer6 <- read.table("tmp_kmer6.fragpatterns.txt",header=T,as.is=T)
kmer6 <- kmer6[1:8,]
kmer7 <- read.table("tmp_kmer7.fragpatterns.txt",header=T,as.is=T)
kmer7 <- kmer7[1:8,]
kmer8 <- read.table("tmp_kmer8.fragpatterns.txt",header=T,as.is=T)
kmer8 <- kmer8[1:8,]

png("plot.png",width=800,height=1000,type="cairo")
par(mfrow=c(4,2),c(3, 2, 2, 1))

cols = c("black",rainbow(9)[1:7])
matplot(as.numeric(sub("X","",sub("X.","-",names(ref[,c(-2,-1)]),fixed=T),fixed=T)),cbind(t(ref[ref$Type=="ThreePrime",c(-2,-1)])[,1],t(kmer2[kmer2$Type=="ThreePrime",c(-2,-1)])[,1],t(kmer3[kmer3$Type=="ThreePrime",c(-2,-1)])[,1],t(kmer4[kmer4$Type=="ThreePrime",c(-2,-1)])[,1],t(kmer5[kmer5$Type=="ThreePrime",c(-2,-1)])[,1],t(kmer6[kmer6$Type=="ThreePrime",c(-2,-1)])[,1],t(kmer7[kmer7$Type=="ThreePrime",c(-2,-1)])[,1],t(kmer8[kmer8$Type=="ThreePrime",c(-2,-1)])[,1]),type="b",xlab="Position",ylab="Base composition",lty=1,main="Fwd read: A bases",ylim=c(0.0,0.6),xlim=c(-5,12),pch=19,col=cols)
abline(v=1)
abline(v=8)
legend("topright",c("original",sapply(2:8,FUN=function(x)sprintf("%dmer",x))),fill=cols)

matplot(as.numeric(sub("X","",sub("X.","-",names(ref[,c(-2,-1)]),fixed=T),fixed=T)),cbind(t(ref[ref$Type!="ThreePrime",c(-2,-1)])[,1],t(kmer2[kmer2$Type!="ThreePrime",c(-2,-1)])[,1],t(kmer3[kmer3$Type!="ThreePrime",c(-2,-1)])[,1],t(kmer4[kmer4$Type!="ThreePrime",c(-2,-1)])[,1],t(kmer5[kmer5$Type!="ThreePrime",c(-2,-1)])[,1],t(kmer6[kmer6$Type!="ThreePrime",c(-2,-1)])[,1],t(kmer7[kmer7$Type!="ThreePrime",c(-2,-1)])[,1],t(kmer8[kmer8$Type!="ThreePrime",c(-2,-1)])[,1]),type="b",xlab="Position",ylab="Base composition",lty=1,main="Rev read: A bases",ylim=c(0.0,0.6),xlim=c(-5,12),pch=19,col=cols)
abline(v=1)
abline(v=8)
legend("topright",c("original",sapply(2:8,FUN=function(x)sprintf("%dmer",x))),fill=cols)

matplot(as.numeric(sub("X","",sub("X.","-",names(ref[,c(-2,-1)]),fixed=T),fixed=T)),cbind(t(ref[ref$Type=="ThreePrime",c(-2,-1)])[,2],t(kmer2[kmer2$Type=="ThreePrime",c(-2,-1)])[,2],t(kmer3[kmer3$Type=="ThreePrime",c(-2,-1)])[,2],t(kmer4[kmer4$Type=="ThreePrime",c(-2,-1)])[,2],t(kmer5[kmer5$Type=="ThreePrime",c(-2,-1)])[,2],t(kmer6[kmer6$Type=="ThreePrime",c(-2,-1)])[,2],t(kmer7[kmer7$Type=="ThreePrime",c(-2,-1)])[,2],t(kmer8[kmer8$Type=="ThreePrime",c(-2,-1)])[,2]),type="b",xlab="Position",ylab="Base composition",lty=1,main="Fwd read: C bases",ylim=c(0.0,0.6),xlim=c(-5,12),pch=19,col=cols)
abline(v=1)
abline(v=8)
legend("topright",c("original",sapply(2:8,FUN=function(x)sprintf("%dmer",x))),fill=cols)

matplot(as.numeric(sub("X","",sub("X.","-",names(ref[,c(-2,-1)]),fixed=T),fixed=T)),cbind(t(ref[ref$Type!="ThreePrime",c(-2,-1)])[,2],t(kmer2[kmer2$Type!="ThreePrime",c(-2,-1)])[,2],t(kmer3[kmer3$Type!="ThreePrime",c(-2,-1)])[,2],t(kmer4[kmer4$Type!="ThreePrime",c(-2,-1)])[,2],t(kmer5[kmer5$Type!="ThreePrime",c(-2,-1)])[,2],t(kmer6[kmer6$Type!="ThreePrime",c(-2,-1)])[,2],t(kmer7[kmer7$Type!="ThreePrime",c(-2,-1)])[,2],t(kmer8[kmer8$Type!="ThreePrime",c(-2,-1)])[,2]),type="b",xlab="Position",ylab="Base composition",lty=1,main="Rev read: C bases",ylim=c(0.0,0.6),xlim=c(-5,12),pch=19,col=cols)
abline(v=1)
abline(v=8)
legend("topright",c("original",sapply(2:8,FUN=function(x)sprintf("%dmer",x))),fill=cols)

matplot(as.numeric(sub("X","",sub("X.","-",names(ref[,c(-2,-1)]),fixed=T),fixed=T)),cbind(t(ref[ref$Type=="ThreePrime",c(-2,-1)])[,3],t(kmer2[kmer2$Type=="ThreePrime",c(-2,-1)])[,3],t(kmer3[kmer3$Type=="ThreePrime",c(-2,-1)])[,3],t(kmer4[kmer4$Type=="ThreePrime",c(-2,-1)])[,3],t(kmer5[kmer5$Type=="ThreePrime",c(-2,-1)])[,3],t(kmer6[kmer6$Type=="ThreePrime",c(-2,-1)])[,3],t(kmer7[kmer7$Type=="ThreePrime",c(-2,-1)])[,3],t(kmer8[kmer8$Type=="ThreePrime",c(-2,-1)])[,3]),type="b",xlab="Position",ylab="Base composition",lty=1,main="Fwd read: G bases",ylim=c(0.0,0.6),xlim=c(-5,12),pch=19,col=cols)
abline(v=1)
abline(v=8)
legend("topright",c("original",sapply(2:8,FUN=function(x)sprintf("%dmer",x))),fill=cols)

matplot(as.numeric(sub("X","",sub("X.","-",names(ref[,c(-2,-1)]),fixed=T),fixed=T)),cbind(t(ref[ref$Type!="ThreePrime",c(-2,-1)])[,3],t(kmer2[kmer2$Type!="ThreePrime",c(-2,-1)])[,3],t(kmer3[kmer3$Type!="ThreePrime",c(-2,-1)])[,3],t(kmer4[kmer4$Type!="ThreePrime",c(-2,-1)])[,3],t(kmer5[kmer5$Type!="ThreePrime",c(-2,-1)])[,3],t(kmer6[kmer6$Type!="ThreePrime",c(-2,-1)])[,3],t(kmer7[kmer7$Type!="ThreePrime",c(-2,-1)])[,3],t(kmer8[kmer8$Type!="ThreePrime",c(-2,-1)])[,3]),type="b",xlab="Position",ylab="Base composition",lty=1,main="Rev read: G bases",ylim=c(0.0,0.6),xlim=c(-5,12),pch=19,col=cols)
abline(v=1)
abline(v=8)
legend("topright",c("original",sapply(2:8,FUN=function(x)sprintf("%dmer",x))),fill=cols)

matplot(as.numeric(sub("X","",sub("X.","-",names(ref[,c(-2,-1)]),fixed=T),fixed=T)),cbind(t(ref[ref$Type=="ThreePrime",c(-2,-1)])[,4],t(kmer2[kmer2$Type=="ThreePrime",c(-2,-1)])[,4],t(kmer3[kmer3$Type=="ThreePrime",c(-2,-1)])[,4],t(kmer4[kmer4$Type=="ThreePrime",c(-2,-1)])[,4],t(kmer5[kmer5$Type=="ThreePrime",c(-2,-1)])[,4],t(kmer6[kmer6$Type=="ThreePrime",c(-2,-1)])[,4],t(kmer7[kmer7$Type=="ThreePrime",c(-2,-1)])[,4],t(kmer8[kmer8$Type=="ThreePrime",c(-2,-1)])[,4]),type="b",xlab="Position",ylab="Base composition",lty=1,main="Fwd read: T bases",ylim=c(0.0,0.6),xlim=c(-5,12),pch=19,col=cols)
abline(v=1)
abline(v=8)
legend("topright",c("original",sapply(2:8,FUN=function(x)sprintf("%dmer",x))),fill=cols)

matplot(as.numeric(sub("X","",sub("X.","-",names(ref[,c(-2,-1)]),fixed=T),fixed=T)),cbind(t(ref[ref$Type!="ThreePrime",c(-2,-1)])[,4],t(kmer2[kmer2$Type!="ThreePrime",c(-2,-1)])[,4],t(kmer3[kmer3$Type!="ThreePrime",c(-2,-1)])[,4],t(kmer4[kmer4$Type!="ThreePrime",c(-2,-1)])[,4],t(kmer5[kmer5$Type!="ThreePrime",c(-2,-1)])[,4],t(kmer6[kmer6$Type!="ThreePrime",c(-2,-1)])[,4],t(kmer7[kmer7$Type!="ThreePrime",c(-2,-1)])[,4],t(kmer8[kmer8$Type!="ThreePrime",c(-2,-1)])[,4]),type="b",xlab="Position",ylab="Base composition",lty=1,main="Rev read: T bases",ylim=c(0.0,0.6),xlim=c(-5,12),pch=19,col=cols)
abline(v=1)
abline(v=8)
legend("topright",c("original",sapply(2:8,FUN=function(x)sprintf("%dmer",x))),fill=cols)

dev.off()

}
