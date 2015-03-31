arg <- commandArgs(trailingOnly = TRUE)
data <- read.table(arg[1],header=F)
mds1 <- cmdscale(as.dist(data), k=2, add=T)
distances <- dist(mds1$points, diag=T) + mds1$ac
write.table(as.matrix(distances), file=arg[2], quote=F, sep="\t",row.names=F,col.names=F)