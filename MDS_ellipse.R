library(car)
library(sp)
arg <- commandArgs(trailingOnly = TRUE)

data <- read.table(arg[1],header=F)
mds1 <- cmdscale(as.dist(data), k=2, add=T)
labels <- read.table(arg[2],header=T)
results <- cbind(mds1$points, labels)
#Calculate all confidence ellipses
P1.ellipse.99 <- dataEllipse(results[which(results[,3] == "P1"),1], results[which(results[,3] == "P1"),2], levels=.99, draw=F)
P1.ellipse.95 <- dataEllipse(results[which(results[,3] == "P1"),1], results[which(results[,3] == "P1"),2], levels=.95,draw=F)
P1.ellipse.90 <- dataEllipse(results[which(results[,3] == "P1"),1], results[which(results[,3] == "P1"),2], levels=.90,draw=F)
P1.ellipse.85 <- dataEllipse(results[which(results[,3] == "P1"),1], results[which(results[,3] == "P1"),2], levels=.85,draw=F)
P1.ellipse.80 <- dataEllipse(results[which(results[,3] == "P1"),1], results[which(results[,3] == "P1"),2], levels=.80,draw=F)
P1.ellipse.75 <- dataEllipse(results[which(results[,3] == "P1"),1], results[which(results[,3] == "P1"),2], levels=.75,draw=F)



#Calculate all confidence ellipses
P2.ellipse.99 <- dataEllipse(results[which(results[,3] == "P2"),1], results[which(results[,3] == "P2"),2], levels=.99,draw=F)
P2.ellipse.95 <- dataEllipse(results[which(results[,3] == "P2"),1], results[which(results[,3] == "P2"),2], levels=.95,draw=F)
P2.ellipse.90 <- dataEllipse(results[which(results[,3] == "P2"),1], results[which(results[,3] == "P2"),2], levels=.90,draw=F)
P2.ellipse.85 <- dataEllipse(results[which(results[,3] == "P2"),1], results[which(results[,3] == "P2"),2], levels=.85,draw=F)
P2.ellipse.80 <- dataEllipse(results[which(results[,3] == "P2"),1], results[which(results[,3] == "P2"),2], levels=.80,draw=F)
P2.ellipse.75 <- dataEllipse(results[which(results[,3] == "P2"),1], results[which(results[,3] == "P2"),2], levels=.75,draw=F)

#calculate whether points are in ellipses for hybrids
P1.hyb <- as.matrix(point.in.polygon(results[which(results[,3] == "H"),1],results[which(results[,3] == "P1"),2], P1.ellipse.99[,1],P1.ellipse.99[,2]))
P1.hyb <- cbind(point.in.polygon(results[which(results[,3] == "H"),1],results[which(results[,3] == "P1"),2], P1.ellipse.95[,1],P1.ellipse.95[,2]), P1.hyb)
P1.hyb <- cbind(point.in.polygon(results[which(results[,3] == "H"),1],results[which(results[,3] == "P1"),2], P1.ellipse.90[,1],P1.ellipse.90[,2]), P1.hyb)
P1.hyb <- cbind(point.in.polygon(results[which(results[,3] == "H"),1],results[which(results[,3] == "P1"),2], P1.ellipse.85[,1],P1.ellipse.85[,2]), P1.hyb)
P1.hyb <- cbind(point.in.polygon(results[which(results[,3] == "H"),1],results[which(results[,3] == "P1"),2], P1.ellipse.80[,1],P1.ellipse.80[,2]), P1.hyb)
P1.hyb <- cbind(point.in.polygon(results[which(results[,3] == "H"),1],results[which(results[,3] == "P1"),2], P1.ellipse.75[,1],P1.ellipse.75[,2]), P1.hyb)

P2.hyb <- as.matrix(point.in.polygon(results[which(results[,3] == "H"),1],results[which(results[,3] == "P2"),2], P2.ellipse.99[,1],P2.ellipse.99[,2]))
P2.hyb <- cbind(point.in.polygon(results[which(results[,3] == "H"),1],results[which(results[,3] == "P2"),2], P2.ellipse.95[,1],P2.ellipse.95[,2]), P2.hyb)
P2.hyb <- cbind(point.in.polygon(results[which(results[,3] == "H"),1],results[which(results[,3] == "P2"),2], P2.ellipse.90[,1],P2.ellipse.90[,2]), P2.hyb)
P2.hyb <- cbind(point.in.polygon(results[which(results[,3] == "H"),1],results[which(results[,3] == "P2"),2], P2.ellipse.85[,1],P2.ellipse.85[,2]), P2.hyb)
P2.hyb <- cbind(point.in.polygon(results[which(results[,3] == "H"),1],results[which(results[,3] == "P2"),2], P2.ellipse.80[,1],P2.ellipse.80[,2]), P2.hyb)
P2.hyb <- cbind(point.in.polygon(results[which(results[,3] == "H"),1],results[which(results[,3] == "P2"),2], P2.ellipse.75[,1],P2.ellipse.75[,2]), P2.hyb)


overlap <- as.matrix(as.logical(sum(point.in.polygon(P1.ellipse.99[,1],P1.ellipse.99[,2], P2.ellipse.99[,1],P2.ellipse.99[,2]))))
overlap <- cbind(as.logical(sum(point.in.polygon(P1.ellipse.95[,1],P1.ellipse.95[,2], P2.ellipse.95[,1],P2.ellipse.95[,2]))),overlap)
overlap <- cbind(as.logical(sum(point.in.polygon(P1.ellipse.90[,1],P1.ellipse.90[,2], P2.ellipse.90[,1],P2.ellipse.90[,2]))),overlap)
overlap <- cbind(as.logical(sum(point.in.polygon(P1.ellipse.85[,1],P1.ellipse.85[,2], P2.ellipse.85[,1],P2.ellipse.85[,2]))),overlap)
overlap <- cbind(as.logical(sum(point.in.polygon(P1.ellipse.80[,1],P1.ellipse.80[,2], P2.ellipse.80[,1],P2.ellipse.80[,2]))),overlap)
overlap <- cbind(as.logical(sum(point.in.polygon(P1.ellipse.75[,1],P1.ellipse.75[,2], P2.ellipse.75[,1],P2.ellipse.75[,2]))),overlap)

write.table(P1.hyb, file=arg[3], quote=F, eol="\n",sep="\t",row.names=F,col.names=F)
write.table(P2.hyb, file=arg[4], quote=F, eol="\n",sep="\t",row.names=F,col.names=F)
write.table(overlap, file=arg[5], quote=F, eol="\n",sep="\t",row.names=F,col.names=F)