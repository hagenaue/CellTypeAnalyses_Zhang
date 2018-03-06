#These code documents contain all of the code documenting the cell type analyses that we applied to Zhang's purified cell RNAseq data to validate our methodology
#Megan Hagenauer and Alek Pankonin
#October 27, 2016

#**************************************

#determining how tightly we can track in silico ratios:

AstrocyteMixtures<-matrix(0, length(names(table(ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag))), 30)
row.names(AstrocyteMixtures)<-names(table(ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag))
AstrocyteMixtureRatios<-c(rep(15, 5), rep(20, 5), rep(25, 5), rep(30,5), rep(35,5), rep(40,5))

for(i in c(15, 20, 25, 30, 35, 40)){
  
  temp3<-matrix(0, length(names(table(ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag))), 5)  
  
  for(j in 1:5){
    temp<-cbind(ZscoreZhang_Expression_CellType_NoPrimaryOverlap[, sample(c(15:16), i, replace=T)],
                ZscoreZhang_Expression_CellType_NoPrimaryOverlap[, sample(c(17:28), (100-i), replace=T)])
    temp2<-apply(temp, 1, mean)
    temp3[,j]<-tapply(temp2, ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag, mean)
  }
  
  AstrocyteMixtures[,c(i:(i+4))-14]<-temp3
}

head(AstrocyteMixtures)

plot(AstrocyteMixtures[1,]~AstrocyteMixtureRatios, pch=18)
points(AstrocyteMixtures[2,]~AstrocyteMixtureRatios, col=2, pch=18)
points(AstrocyteMixtures[3,]~AstrocyteMixtureRatios, col=3, pch=18)
points(AstrocyteMixtures[4,]~AstrocyteMixtureRatios, col=4, pch=18)
points(AstrocyteMixtures[5,]~AstrocyteMixtureRatios, col=5, pch=18)

png("AstrocyteMixtureRatiosVsIndex_NoZhang.png")
plot(apply(AstrocyteMixtures[c(1:4),], 2, mean)~AstrocyteMixtureRatios, pch=18, ylab="AstrocyteIndex")
BestFitLine<-lm(apply(AstrocyteMixtures[c(1:4),], 2, mean)~AstrocyteMixtureRatios)
abline(BestFitLine)
abline(a=(BestFitLine$coefficients[1]+summary.lm(BestFitLine)[[6]]), b=BestFitLine$coefficients[2], col="grey")
abline(a=(BestFitLine$coefficients[1]-summary.lm(BestFitLine)[[6]]), b=BestFitLine$coefficients[2], col="grey")

abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*30), b=0, col="red")
abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*27.5), b=0, col="green")
abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*25), b=0, col="red")
dev.off()

#prettier version:

pdf("AstrocyteMixtureRatiosVsIndex_NoZhang.pdf", height=5, width=4)
plot(apply(AstrocyteMixtures[c(1:4),], 2, mean)~AstrocyteMixtureRatios, pch=18, ylab="Astrocyte Index", xlab="Astrocyte Mixture (% Astrocyte)")
BestFitLine<-lm(apply(AstrocyteMixtures[c(1:4),], 2, mean)~AstrocyteMixtureRatios)
abline(BestFitLine, lwd=2)
abline(a=(BestFitLine$coefficients[1]+summary.lm(BestFitLine)[[6]]), b=BestFitLine$coefficients[2], col="grey", lwd=2)
abline(a=(BestFitLine$coefficients[1]-summary.lm(BestFitLine)[[6]]), b=BestFitLine$coefficients[2], col="grey", lwd=2)

abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*30), b=0, col="red", lwd=2)
abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*27.5), b=0, col="green", lwd=2)
abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*25), b=0, col="red", lwd=2)
dev.off()

#Version 2: Using whole brain as the background instead of a random sampling of cells
AstrocyteMixtures<-matrix(0, length(names(table(ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag))), 30)
row.names(AstrocyteMixtures)<-names(table(ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag))
AstrocyteMixtureRatios<-c(rep(15, 5), rep(20, 5), rep(25, 5), rep(30,5), rep(35,5), rep(40,5))

for(i in c(15, 20, 25, 30, 35, 40)){
  
  temp3<-matrix(0, length(names(table(ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag))), 5)  
  
  for(j in 1:5){
    temp<-cbind(ZscoreZhang_Expression_CellType_NoPrimaryOverlap[, sample(c(15:16), i, replace=T)],
                ZscoreZhang_Expression_CellType_NoPrimaryOverlap[, sample(c(29:31), (100-i), replace=T)])
    temp2<-apply(temp, 1, mean)
    temp3[,j]<-tapply(temp2, ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag, mean)
  }
  
  AstrocyteMixtures[,c(i:(i+4))-14]<-temp3
}

head(AstrocyteMixtures)

plot(AstrocyteMixtures[1,]~AstrocyteMixtureRatios, pch=18)
points(AstrocyteMixtures[2,]~AstrocyteMixtureRatios, col=2, pch=18)
points(AstrocyteMixtures[3,]~AstrocyteMixtureRatios, col=3, pch=18)
points(AstrocyteMixtures[4,]~AstrocyteMixtureRatios, col=4, pch=18)
points(AstrocyteMixtures[5,]~AstrocyteMixtureRatios, col=5, pch=18)

png("AstrocyteMixtureRatiosVsIndex_NoZhang_WB.png")
plot(apply(AstrocyteMixtures[c(1:4),], 2, mean)~AstrocyteMixtureRatios, pch=18, ylab="AstrocyteIndex")
BestFitLine<-lm(apply(AstrocyteMixtures[c(1:4),], 2, mean)~AstrocyteMixtureRatios)
abline(BestFitLine)
abline(a=(BestFitLine$coefficients[1]+summary.lm(BestFitLine)[[6]]), b=BestFitLine$coefficients[2], col="grey")
abline(a=(BestFitLine$coefficients[1]-summary.lm(BestFitLine)[[6]]), b=BestFitLine$coefficients[2], col="grey")

abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*30), b=0, col="red")
abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*27.5), b=0, col="green")
abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*25), b=0, col="red")
dev.off()



OligodendrocyteMixtures<-matrix(0, length(names(table(ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag))), 30)
row.names(OligodendrocyteMixtures)<-names(table(ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag))
OligodendrocyteMixtureRatios<-c(rep(15, 5), rep(20, 5), rep(25, 5), rep(30,5), rep(35,5), rep(40,5))

for(i in c(15, 20, 25, 30, 35, 40)){
  
  temp3<-matrix(0, length(names(table(ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag))), 5)  
  
  for(j in 1:5){
    temp<-cbind(ZscoreZhang_Expression_CellType_NoPrimaryOverlap[, sample(c(23:24), i, replace=T)],
                ZscoreZhang_Expression_CellType_NoPrimaryOverlap[, sample(c(15:22, 25:28), (100-i), replace=T)])
    temp2<-apply(temp, 1, mean)
    temp3[,j]<-tapply(temp2, ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag, mean)
  }
  
  OligodendrocyteMixtures[,c(i:(i+4))-14]<-temp3
}

head(OligodendrocyteMixtures)

png("OligodendrocyteMixtureRatiosVsIndex_NoZhang.png")
plot(apply(OligodendrocyteMixtures[c(29:33),], 2, mean)~OligodendrocyteMixtureRatios, pch=18, ylab="OligodendrocyteIndex")
BestFitLine<-lm(apply(OligodendrocyteMixtures[c(29:33),], 2, mean)~OligodendrocyteMixtureRatios)
abline(BestFitLine)
abline(a=(BestFitLine$coefficients[1]+summary.lm(BestFitLine)[[6]]), b=BestFitLine$coefficients[2], col="grey")
abline(a=(BestFitLine$coefficients[1]-summary.lm(BestFitLine)[[6]]), b=BestFitLine$coefficients[2], col="grey")

abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*30), b=0, col="red")
abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*27.5), b=0, col="green")
abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*25), b=0, col="red")
dev.off()

#prettier version:

pdf("OligodendrocyteMixtureRatiosVsIndex_NoZhang.pdf", height=5, width=4)
plot(apply(OligodendrocyteMixtures[c(29:33),], 2, mean)~OligodendrocyteMixtureRatios, pch=18, ylab="Oligodendrocyte Index", xlab="Oligodendrocyte Mixture (% Oligodendrocyte)")
BestFitLine<-lm(apply(OligodendrocyteMixtures[c(29:33),], 2, mean)~OligodendrocyteMixtureRatios)
abline(BestFitLine, lwd=2)
abline(a=(BestFitLine$coefficients[1]+summary.lm(BestFitLine)[[6]]), b=BestFitLine$coefficients[2], col="grey", lwd=2)
abline(a=(BestFitLine$coefficients[1]-summary.lm(BestFitLine)[[6]]), b=BestFitLine$coefficients[2], col="grey", lwd=2)

abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*30), b=0, col="red", lwd=2)
abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*27.5), b=0, col="green", lwd=2)
abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*25), b=0, col="red", lwd=2)
dev.off()


OligodendrocyteMixtures<-matrix(0, length(names(table(ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag))), 30)
row.names(OligodendrocyteMixtures)<-names(table(ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag))
OligodendrocyteMixtureRatios<-c(rep(15, 5), rep(20, 5), rep(25, 5), rep(30,5), rep(35,5), rep(40,5))

for(i in c(15, 20, 25, 30, 35, 40)){
  
  temp3<-matrix(0, length(names(table(ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag))), 5)  
  
  for(j in 1:5){
    temp<-cbind(ZscoreZhang_Expression_CellType_NoPrimaryOverlap[, sample(c(23:24), i, replace=T)],
                ZscoreZhang_Expression_CellType_NoPrimaryOverlap[, sample(c(29:31), (100-i), replace=T)])
    temp2<-apply(temp, 1, mean)
    temp3[,j]<-tapply(temp2, ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag, mean)
  }
  
  OligodendrocyteMixtures[,c(i:(i+4))-14]<-temp3
}

head(OligodendrocyteMixtures)

png("OligodendrocyteMixtureRatiosVsIndex_NoZhang_WB.png")
plot(apply(OligodendrocyteMixtures[c(29:33),], 2, mean)~OligodendrocyteMixtureRatios, pch=18, ylab="OligodendrocyteIndex")
BestFitLine<-lm(apply(OligodendrocyteMixtures[c(29:33),], 2, mean)~OligodendrocyteMixtureRatios)
abline(BestFitLine)
abline(a=(BestFitLine$coefficients[1]+summary.lm(BestFitLine)[[6]]), b=BestFitLine$coefficients[2], col="grey")
abline(a=(BestFitLine$coefficients[1]-summary.lm(BestFitLine)[[6]]), b=BestFitLine$coefficients[2], col="grey")

abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*30), b=0, col="red")
abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*27.5), b=0, col="green")
abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*25), b=0, col="red")
dev.off()

#I'm going to stop doing the whole brain mixtures - i think there isn't a reasonable amount of variability there, since there were only 3 WB samples.  I general, that is a definite problem with this entire analysis.


MicrogliaMixtures<-matrix(0, length(names(table(ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag))), 30)
row.names(MicrogliaMixtures)<-names(table(ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag))
MicrogliaMixtureRatios<-c(rep(5, 5), rep(10, 5), rep(15, 5), rep(20,5), rep(25,5), rep(30,5))

for(i in c(5, 10, 15, 20, 25, 30)){
  
  temp3<-matrix(0, length(names(table(ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag))), 5)  
  
  for(j in 1:5){
    temp<-cbind(ZscoreZhang_Expression_CellType_NoPrimaryOverlap[, sample(c(25:26), i, replace=T)],
                ZscoreZhang_Expression_CellType_NoPrimaryOverlap[, sample(c(15:24, 27:28), (100-i), replace=T)])
    temp2<-apply(temp, 1, mean)
    temp3[,j]<-tapply(temp2, ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag, mean)
  }
  
  MicrogliaMixtures[,c(i:(i+4))-4]<-temp3
}

head(MicrogliaMixtures)

plot(MicrogliaMixtures[1,]~MicrogliaMixtureRatios, pch=18)
points(MicrogliaMixtures[2,]~MicrogliaMixtureRatios, col=2, pch=18)
points(MicrogliaMixtures[3,]~MicrogliaMixtureRatios, col=3, pch=18)
points(MicrogliaMixtures[4,]~MicrogliaMixtureRatios, col=4, pch=18)
points(MicrogliaMixtures[5,]~MicrogliaMixtureRatios, col=5, pch=18)

png("MicrogliaMixtureRatiosVsIndex_NoZhang.png")
plot(apply(MicrogliaMixtures[c(10:11),], 2, mean)~MicrogliaMixtureRatios, pch=18, ylab="MicrogliaIndex")
BestFitLine<-lm(apply(MicrogliaMixtures[c(10:11),], 2, mean)~MicrogliaMixtureRatios)
abline(BestFitLine)
abline(a=(BestFitLine$coefficients[1]+summary.lm(BestFitLine)[[6]]), b=BestFitLine$coefficients[2], col="grey")
abline(a=(BestFitLine$coefficients[1]-summary.lm(BestFitLine)[[6]]), b=BestFitLine$coefficients[2], col="grey")

abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*15), b=0, col="red")
abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*17.5), b=0, col="green")
abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*20), b=0, col="red")
dev.off()

#Prettier version:

pdf("MicrogliaMixtureRatiosVsIndex_NoZhang.pdf", height=5, width=4)
plot(apply(MicrogliaMixtures[c(10:11),], 2, mean)~MicrogliaMixtureRatios, pch=18, ylab="Microglia Index", xlab="Microglia Mixture (% Microglia)")
BestFitLine<-lm(apply(MicrogliaMixtures[c(10:11),], 2, mean)~MicrogliaMixtureRatios)
abline(BestFitLine, lwd=2)
abline(a=(BestFitLine$coefficients[1]+summary.lm(BestFitLine)[[6]]), b=BestFitLine$coefficients[2], col="grey", lwd=2)
abline(a=(BestFitLine$coefficients[1]-summary.lm(BestFitLine)[[6]]), b=BestFitLine$coefficients[2], col="grey", lwd=2)

abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*15), b=0, col="red", lwd=2)
abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*17.5), b=0, col="green", lwd=2)
abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*20), b=0, col="red", lwd=2)
dev.off()


EndothelialMixtures<-matrix(0, length(names(table(ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag))), 30)
row.names(EndothelialMixtures)<-names(table(ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag))
EndothelialMixtureRatios<-c(rep(1, 5), rep(6, 5), rep(11, 5), rep(16,5), rep(21,5), rep(26,5))

for(i in c(1, 6, 11, 16, 21, 26)){
  
  temp3<-matrix(0, length(names(table(ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag))), 5)  
  
  for(j in 1:5){
    temp<-cbind(ZscoreZhang_Expression_CellType_NoPrimaryOverlap[, sample(c(27:28), i, replace=T)],
                ZscoreZhang_Expression_CellType_NoPrimaryOverlap[, sample(c(15:26), (100-i), replace=T)])
    temp2<-apply(temp, 1, mean)
    temp3[,j]<-tapply(temp2, ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag, mean)
  }
  
  EndothelialMixtures[,c(i:(i+4))]<-temp3
}

head(EndothelialMixtures)

plot(EndothelialMixtures[1,]~EndothelialMixtureRatios, pch=18)
points(EndothelialMixtures[2,]~EndothelialMixtureRatios, col=2, pch=18)
points(EndothelialMixtures[3,]~EndothelialMixtureRatios, col=3, pch=18)
points(EndothelialMixtures[4,]~EndothelialMixtureRatios, col=4, pch=18)
points(EndothelialMixtures[5,]~EndothelialMixtureRatios, col=5, pch=18)

png("EndothelialMixtureRatiosVsIndex_NoZhang.png")
plot(apply(EndothelialMixtures[c(6:8),], 2, mean)~EndothelialMixtureRatios, pch=18, ylab="EndothelialIndex")
BestFitLine<-lm(apply(EndothelialMixtures[c(6:8),], 2, mean)~EndothelialMixtureRatios)
abline(BestFitLine)
abline(a=(BestFitLine$coefficients[1]+summary.lm(BestFitLine)[[6]]), b=BestFitLine$coefficients[2], col="grey")
abline(a=(BestFitLine$coefficients[1]-summary.lm(BestFitLine)[[6]]), b=BestFitLine$coefficients[2], col="grey")

abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*5), b=0, col="red")
abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*7.5), b=0, col="green")
abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*10), b=0, col="red")
dev.off()

#prettier version:

pdf("EndothelialMixtureRatiosVsIndex_NoZhang.pdf", height=5, width=4)
plot(apply(EndothelialMixtures[c(6:8),], 2, mean)~EndothelialMixtureRatios, pch=18, ylab="Endothelial Index", xlab="Endothelial Mixture (% Endothelial)")
BestFitLine<-lm(apply(EndothelialMixtures[c(6:8),], 2, mean)~EndothelialMixtureRatios)
abline(BestFitLine, lwd=2)
abline(a=(BestFitLine$coefficients[1]+summary.lm(BestFitLine)[[6]]), b=BestFitLine$coefficients[2], col="grey", lwd=2)
abline(a=(BestFitLine$coefficients[1]-summary.lm(BestFitLine)[[6]]), b=BestFitLine$coefficients[2], col="grey", lwd=2)

abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*5), b=0, col="red", lwd=2)
abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*7.5), b=0, col="green", lwd=2)
abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*10), b=0, col="red", lwd=2)
dev.off()



NeuronMixtures<-matrix(0, length(names(table(ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag))), 30)
row.names(NeuronMixtures)<-names(table(ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag))
NeuronMixtureRatios<-c(rep(25, 5), rep(30, 5), rep(35, 5), rep(40,5), rep(45,5), rep(50,5))

for(i in c(25, 30, 35, 40, 45, 50)){
  
  temp3<-matrix(0, length(names(table(ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag))), 5)  
  
  for(j in 1:5){
    temp<-cbind(ZscoreZhang_Expression_CellType_NoPrimaryOverlap[, sample(c(17:18), i, replace=T)],
                ZscoreZhang_Expression_CellType_NoPrimaryOverlap[, sample(c(15:16, 19:28), (100-i), replace=T)])
    temp2<-apply(temp, 1, mean)
    temp3[,j]<-tapply(temp2, ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag, mean)
  }
  
  NeuronMixtures[,c(i:(i+4))-24]<-temp3
}

head(NeuronMixtures)

plot(NeuronMixtures[1,]~NeuronMixtureRatios, pch=18)
points(NeuronMixtures[2,]~NeuronMixtureRatios, col=2, pch=18)
points(NeuronMixtures[3,]~NeuronMixtureRatios, col=3, pch=18)
points(NeuronMixtures[4,]~NeuronMixtureRatios, col=4, pch=18)
points(NeuronMixtures[5,]~NeuronMixtureRatios, col=5, pch=18)

png("NeuronMixtureRatiosVsIndex_NoZhang.png")
plot(apply(NeuronMixtures[c(16:17),], 2, mean)~NeuronMixtureRatios, pch=18, ylab="NeuronIndex")
BestFitLine<-lm(apply(NeuronMixtures[c(16:17),], 2, mean)~NeuronMixtureRatios)
abline(BestFitLine)
abline(a=(BestFitLine$coefficients[1]+summary.lm(BestFitLine)[[6]]), b=BestFitLine$coefficients[2], col="grey")
abline(a=(BestFitLine$coefficients[1]-summary.lm(BestFitLine)[[6]]), b=BestFitLine$coefficients[2], col="grey")

abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*35), b=0, col="red")
abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*37.5), b=0, col="green")
abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*40), b=0, col="red")
dev.off()

#Prettier version:

pdf("NeuronMixtureRatiosVsIndex_NoZhang.pdf", height=5, width=4)
plot(apply(NeuronMixtures[c(16:17),], 2, mean)~NeuronMixtureRatios, pch=18, ylab="Neuron Index", xlab="Neuron Mixture (% Neuron)")
BestFitLine<-lm(apply(NeuronMixtures[c(16:17),], 2, mean)~NeuronMixtureRatios)
abline(BestFitLine, lwd=2)
abline(a=(BestFitLine$coefficients[1]+summary.lm(BestFitLine)[[6]]), b=BestFitLine$coefficients[2], col="grey", lwd=2)
abline(a=(BestFitLine$coefficients[1]-summary.lm(BestFitLine)[[6]]), b=BestFitLine$coefficients[2], col="grey", lwd=2)

abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*35), b=0, col="red", lwd=2)
abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*37.5), b=0, col="green", lwd=2)
abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*40), b=0, col="red", lwd=2)
dev.off()


ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag<-matrix(0, nrow=length(table(ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag)), ncol=(length(ZscoreZhang_Expression_CellType_NoPrimaryOverlap[1,])-14))

row.names(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag)<-names(table(ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag))

colnames(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag)<-colnames(ZscoreZhang_Expression_CellType_NoPrimaryOverlap)[-c(1:14)]

for(i in c(15:length(ZscoreZhang_Expression_CellType_NoPrimaryOverlap[1,]))){
  ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag[,i-14]<-tapply(ZscoreZhang_Expression_CellType_NoPrimaryOverlap[,i], ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag, mean)
}


################################
