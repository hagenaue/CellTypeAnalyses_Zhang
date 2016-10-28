#These code documents contain all of the code documenting the cell type analyses that we applied to Zhang's purified cell RNAseq data to validate our methodology
#Megan Hagenauer and Alek Pankonin
#October 27, 2016

#**************************************

JoinedZhang<-read.csv("Zhang_FullDataset.csv", row.names=1, header=T)


#Entering in cell identities:
colnames(JoinedZhang_AsNumMatrix_Log2)
ExperimentalCellType<-c("Astrocyte", "Astrocyte", "Neuron", "Neuron", "OPC", "OPC", "NFO", "NFO", "Oligodendrocyte", "Oligodendrocyte", "Microglia", "Microglia", "Endothelial", "Endothelial", "WholeBrain", "WholeBrain", "WholeBrain")

table(ExperimentalCellType)

cbind(ExperimentalCellType, colnames(JoinedZhang_AsNumMatrix_Log2))



#Taking a look at the Zhang data:
is.numeric(JoinedZhang)
is.numeric(as.matrix(JoinedZhang))
#That works.

GeneNamesForJoinedZhang<-row.names(JoinedZhang)

JoinedZhang_AsNumMatrix<-as.matrix(JoinedZhang)
is.numeric(JoinedZhang_AsNumMatrix)

#Let's take a peek at their distribution:

#I wonder if the principal components of variation in the data correlate with those variables?  

#I should probably log transform the reads first before messing around with PCA:
#Let's check out the average distribution of the signal first:
temp<-apply(JoinedZhang_AsNumMatrix, 1, mean)
png("Histogram_AverageReadsPerProbe.png")
hist(temp, breaks=1000)
dev.off()
max(temp)
#14317.77
median(temp)
[1] 1.996951

temp<-apply(JoinedZhang_AsNumMatrix, 1, max)
png("Histogram_MaxReadsPerProbe.png")
hist(temp, breaks=1000)
dev.off()
#Wow, that is super skewed with a huge number of probes with an average, max, or median of 0. Yeah, let's log transform that data.

#Oh wait - log transformation in RNAseq data is awkward, because there are 0's, which convert into -Inf. It seems that folks typically log transform shifted data instead (data+1)
JoinedZhang_AsNumMatrix_Log2<-log2((JoinedZhang_AsNumMatrix+1))
row.names(JoinedZhang_AsNumMatrix_Log2)<-GeneNamesForJoinedZhang
write.csv(JoinedZhang_AsNumMatrix_Log2, "JoinedZhang_AsNumMatrix_Log2.csv")

temp<-apply(JoinedZhang_AsNumMatrix_Log2, 1, mean)
png("Histogram_AverageReadsPerProbeLog2.png")
hist(temp, breaks=1000)
dev.off()

temp<-apply(JoinedZhang_AsNumMatrix_Log2, 1, max)
png("Histogram_MaxReadsPerProbeLog2.png")
hist(temp, breaks=1000)
dev.off()

#That's still super skewed, but not as bad as before. I'm guessing that it is difficult to get reads >0 for many transcripts in single cell type data. That means that any analysis based on sample-sample correlations is going to be largely driven by a few highly-expressed transcripts.
#Apparently there are other transformation methods out there that work better for RNAseq - VST is implemented by DESeq. VOOM applies a tranformation to the read counts that supposedly then makes the RNAseq data compatible with analyses included in the limma package. https://seqqc.wordpress.com/2015/02/16/should-you-transform-rna-seq-data-log-vst-voom/


#Time to recycle some code:
library(gdata)
library(fields)
library(stats)
library(car)
library(affy)
library(preprocessCore)
library(multtest)

#The fact that this data (unlike our microarray data) isn't quantile normalized may make the results of PCA a little funky. Let's start with sample-sample correlations:

#Visualize the sample-sample correlations using a heatmap:

temp<-cor(JoinedZhang_AsNumMatrix_Log2)
colnames(temp)<-colnames(JoinedZhang_AsNumMatrix_Log2)
row.names(temp)<-colnames(JoinedZhang_AsNumMatrix_Log2)

png("09 Sample Sample Correlations Heatmap.png")
heatmap(temp, main="Visualizing correlations between entire samples", xlab="Red=Less correlated, Light yellow=Highly correlated")
dev.off()
#Note that the heatmap can be tailored to focus on a certain level of correlation by using the command zlim=c(lower limit, upper limit)

#Visualize the sample-sample correlations using a boxplot:
png("09 Boxplot Sample Sample Correlations.png", width=2000, height=300)
boxplot(data.frame(cor(JoinedZhang_AsNumMatrix_Log2)), cex=0.25, las=3, par(cex.axis=0.75), main="Boxplot of sample-sample correlations", xlab="Subject", ylab="Sample-Sample Correlations")
Median10thQuantile<-median(apply((cor(JoinedZhang_AsNumMatrix_Log2)), 1, quantile, 0.1))
MedianQuantile<-median(apply((cor(JoinedZhang_AsNumMatrix_Log2)), 1, quantile, 0.5))
abline(a=Median10thQuantile, b=0, col=2)
abline(a=MedianQuantile, b=0, col=3)
mtext(paste("Median Sample-Sample Correlation=", round(MedianQuantile, digits=3), sep=" ")) 
dev.off()

#It might be worthwile to perhaps try running things with all the genes with extremely low reads thrown out:

temp<-apply(JoinedZhang_AsNumMatrix_Log2, 1, max)
sum(temp<2)
[1] 9787
#That's a lot of genes with 2 reads (log2=1) or less to begin with.  
sum(temp<3)
[1] 11889
#And almost half of them don't meet a standard cut off of 8 reads. Ouch.
sum(temp<3)/length(temp)
#[1] 0.5292939
#...or 53%. Dang. That's worse than the single-cell data. Does that mean that the read depth was shallow?



#Run principal components analysis (PCA) to determine which major gradients of sample-sample correlations exist in the data (i.e. who is similar to whom):
pcaNormFiltered<-prcomp(t(JoinedZhang_AsNumMatrix_Log2))
tmp<-pcaNormFiltered$x[,1:4]
write.table(tmp, "PCA_1_4.txt", sep="\t")

PC1<-pcaNormFiltered$x[,1]
PC2<-pcaNormFiltered$x[,2]

PC3<-pcaNormFiltered$x[,3]
PC4<-pcaNormFiltered$x[,4]

#Output a scree plot for the PCA:
png("09 PCA Scree Plot1.png")
plot(summary(pcaNormFiltered)$importance[2,]~(c(1:ncol(JoinedZhang_AsNumMatrix_Log2))), main="Variance Explained by Each Principal Component", xlab="PC #", ylab="Proportion of Variance Explained", col=2)
dev.off()

png("09 PCA Scree Plot2.png")
plot(summary(pcaNormFiltered)$importance[3,]~(c(1:ncol(JoinedZhang_AsNumMatrix_Log2))), main="Variance Explained by Each Principal Component", xlab="PC #", ylab="Cumulative Proportion of Variance Explained", col=3)
dev.off()

#Output a scatterplot illustrating the relationship between Principal components 1 & 2 (PC1 & PC2):
png("09 PC1 vs PC2.png")
plot(PC1~PC2, main="Principal Components Analysis of Normalized Filtered Data", col=as.factor(ExperimentalCellType))
dev.off()

#Output a scatterplot illustrating the relationship between Principal components 3 & 4 (PC3 & PC4):
png("09 PC3 vs PC4.png")
plot(PC3~PC4, main="Principal Components Analysis of Normalized Filtered Data", col=as.factor(ExperimentalCellType))
dev.off()

#Nice clustering.

###############################################################
#I'm going to z-score my log-transformed data.
#I wonder where the NAs are coming from.
#Ah -after messing around, I believe either NaNs or Inf are produced whenever a gene doesn't have any variability associated with it
#So let's remove the completely invariable data first:

JoinedZhang_StDev<-apply(JoinedZhang_AsNumMatrix_Log2, 1, sd) 
sum(JoinedZhang_StDev==0)
#[1] 5314
#Dang....

JoinedZhang_AsNumMatrix_Log2_NoSD0<-JoinedZhang_AsNumMatrix_Log2[JoinedZhang_StDev>0,]
temp<-GeneNamesForJoinedZhang[-c(22086:22088)]
GeneNamesForJoinedZhang_NoSD0<-temp[JoinedZhang_StDev>0]


ZscoreZhang<-t(scale(t(JoinedZhang_AsNumMatrix_Log2_NoSD0), center=T, scale=T))#Zscores the data 
write.csv(ZscoreZhang, "ZscoreZhang.csv")
sum(is.na(ZscoreZhang))
#looks good now

