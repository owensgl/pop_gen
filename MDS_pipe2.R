arg <- commandArgs(trailingOnly = TRUE)
tmpdata <- read.table(arg[1], header=F)
result <- t.test(tmpdata[which(tmpdata$V2=='P1'),1], tmpdata[which(tmpdata$V2=='P2'),1])
P1ave <- mean(tmpdata[which(tmpdata$V2=='P1'),1])
P2ave <- mean(tmpdata[which(tmpdata$V2=='P2'),1])
final <- NULL
if(result$p.value < 0.05){
  if (P1ave > P2ave){
    final <- c("P2")
  }else{
    final <- c("P1")
  }
}else{
    final <- c("H")
}
out <- c(P1ave, P2ave, result$p.value, final)
write.table(out, file=arg[2], quote=F, eol=":",sep="\t",row.names=F,col.names=F)
