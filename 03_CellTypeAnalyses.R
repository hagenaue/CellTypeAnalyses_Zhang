#These code documents contain all of the code documenting the cell type analyses that we applied to Zhang's purified cell RNAseq data to validate our methodology
#Megan Hagenauer and Alek Pankonin
#October 27, 2016

#**************************************

#My newer version of the cell type analysis: stolen from the ABA code

CellTypeSpecificGenes_Master3<-read.csv("CellTypeSpecificGenes_Master3.csv", header=T)

colnames(CellTypeSpecificGenes_Master3)

# [1] "Umbrella.Cell.Type"    "Specific.Cell.Type"    "Brain.Region"          "Gene.Symbol..Human."  
# [5] "Gene.Symbol..Mouse."   "Species"               "Age"                   "Statistical.Criterion"
# [9] "Specificity"           "Comparison"            "Platform"              "Citation"             
# [13] "Tag"                   "CellType_Primary" 

table(CellTypeSpecificGenes_Master3$CellType_Primary)

#Note: I had to reverse the order of the next 3 lines to properly remove NAs - I should change that before attempting to release the code in a way that other people will use.

colnames(CellTypeSpecificGenes_Master3)[4]<-"GeneSymbol_Human"
colnames(CellTypeSpecificGenes_Master3)[5]<-"GeneSymbol_Mouse"

sum(is.na(CellTypeSpecificGenes_Master3$GeneSymbol_Human))
#[1] 364
sum(is.na(CellTypeSpecificGenes_Master3$GeneSymbol_Mouse))
#[1] 7

#This line depends on the species being used:
CellTypeSpecificGenes_Master3NoNA<-CellTypeSpecificGenes_Master3[is.na(CellTypeSpecificGenes_Master3$GeneSymbol_Mouse)==F,]

CellTypeSpecificGenes_Master3NoNA[,4]<-as.character(CellTypeSpecificGenes_Master3NoNA[,4])
CellTypeSpecificGenes_Master3NoNA[,5]<-as.character(CellTypeSpecificGenes_Master3NoNA[,5])


str(CellTypeSpecificGenes_Master3NoNA)

#Note to self - if I try to make a generalizable version of this code, I will need to change it so that it references the rodent gene symbols when  using rodent data and the human gene symbols when using human data.

#Note: some data will need to be averaged by gene symbol at this point.  This typically is not true of processed RNASeq data
sum(unique(row.names(ZscoreZhang))==F)
#[1] 0

#New version (after averaging by gene symbol if necessary):
temp<-data.frame(row.names(ZscoreZhang), ZscoreZhang, stringsAsFactors=F)

#This code depends on the species of the original dataset
#colnames(temp)[1]<-"GeneSymbol_Human"
colnames(temp)[1]<-"GeneSymbol_Mouse"

sum(is.na(temp[,1]))
#[1] 0

#This code also depends on the species in the original dataset:
#If Human: sum(temp[,1] %in% CellTypeSpecificGenes_Master3[,4])
sum(temp[,1] %in% CellTypeSpecificGenes_Master3[,5])
#[1] 2513

#If human: sum(CellTypeSpecificGenes_Master3[,4]  %in%  temp[,1])
sum(CellTypeSpecificGenes_Master3[,5]  %in%  temp[,1])
# [1] 2914

#Note: NAs were causing a serious problem with this join function.  Fixed now. :)
library(plyr)
#If human: ZscoreZhang_Expression_CellType<-join(CellTypeSpecificGenes_Master3, temp, by="GeneSymbol_Human", type="inner")
ZscoreZhang_Expression_CellType<-join(CellTypeSpecificGenes_Master3, temp, by="GeneSymbol_Mouse", type="inner")
dim(ZscoreZhang_Expression_CellType)
# [1] 2914   31
#It is making all possible combinations - some of the cell type specific genes are found in more than one index.

write.csv(ZscoreZhang_Expression_CellType, "ZscoreZhang_Expression_CellType.csv")

###############################################

AVE_Expression_CellType_Primary_bySample<-matrix(NA, nrow=length(names(table(ZscoreZhang_Expression_CellType$CellType_Primary))), ncol=(ncol(ZscoreZhang_Expression_CellType)-14))

row.names(AVE_Expression_CellType_Primary_bySample)<-names(table(ZscoreZhang_Expression_CellType$CellType_Primary))
colnames(AVE_Expression_CellType_Primary_bySample)<-colnames(temp)[-1]


for(i in c(15:ncol(ZscoreZhang_Expression_CellType))){
  AVE_Expression_CellType_Primary_bySample[,(i-14)]<-tapply(ZscoreZhang_Expression_CellType[,i], ZscoreZhang_Expression_CellType$CellType_Primary, function(y) mean(y, na.rm=T))
}

head(AVE_Expression_CellType_Primary_bySample)

png("CorrMatrixCellTypeVsCellType_HeatMap.png")
heatmap(cor(t(AVE_Expression_CellType_Primary_bySample[,is.na(AVE_Expression_CellType_Primary_bySample[1,])==F])), margins = c(15, 15), cex.lab=0.5)
dev.off()

CorrelationMatrixCellTypeVsCellType<-cor(t(AVE_Expression_CellType_Primary_bySample[,is.na(AVE_Expression_CellType_Primary_bySample[1,])==F]))

write.csv(CorrelationMatrixCellTypeVsCellType, "CorrelationMatrixCellTypeVsCellType.csv")

#let's see some examples:

png("AstrocyteByEndothelial.png")
plot(AVE_Expression_CellType_Primary_bySample[1,]~AVE_Expression_CellType_Primary_bySample[2,])
dev.off()
#I think the sample size per cell type maybe too small to properly run these analyses 

AVE_Expression_CellType_Tag_bySample<-matrix(NA, nrow=length(names(table(ZscoreZhang_Expression_CellType$Tag))), ncol=ncol(ZscoreZhang_Expression_CellType)-14)
row.names(AVE_Expression_CellType_Tag_bySample)<-names(table(ZscoreZhang_Expression_CellType$Tag))
colnames(AVE_Expression_CellType_Tag_bySample)<-colnames(temp)[-1]

for(i in c(15:ncol(ZscoreZhang_Expression_CellType))){
  AVE_Expression_CellType_Tag_bySample[,(i-14)]<-tapply(ZscoreZhang_Expression_CellType[,i], ZscoreZhang_Expression_CellType$Tag, mean)
}

head(AVE_Expression_CellType_Tag_bySample)

png("CorrMatrixCellIndexVsCellIndex_HeatMap.png", width=1000, height=1000)
heatmap(cor(t(AVE_Expression_CellType_Tag_bySample[,is.na(AVE_Expression_CellType_Tag_bySample[1,])==F])), cex.lab=0.3, margins = c(20, 20))
dev.off()

CorrelationMatrixCellIndexVsCellIndex<-cor(t(AVE_Expression_CellType_Tag_bySample[,is.na(AVE_Expression_CellType_Tag_bySample[1,])==F]))

write.csv(CorrelationMatrixCellIndexVsCellIndex, "CorrelationMatrixCellIndexVsCellIndex.csv")

################################################

#Alright, so part of the trouble here is that there hasn't been any removal of overlapping probes yet, and we aren't averaging by tag. Let's go ahead and do that.


#Making a storage matrix to store information about overlap between primary indices:
CellTypeSpecificGenes_Master3_Overlap<-matrix(0, length(table(ZscoreZhang_Expression_CellType$CellType_Primary)), length(table(ZscoreZhang_Expression_CellType$CellType_Primary)) )

colnames(CellTypeSpecificGenes_Master3_Overlap)<-names(table(ZscoreZhang_Expression_CellType$CellType_Primary))
row.names(CellTypeSpecificGenes_Master3_Overlap)<-names(table(ZscoreZhang_Expression_CellType$CellType_Primary))

#Quantifying overlap between primary cell type indices:
for(i in 1: length(table(ZscoreZhang_Expression_CellType$CellType_Primary))){
  for(j in 1: length(table(ZscoreZhang_Expression_CellType$CellType_Primary))){
    
    CellTypeSpecificGenes_Master3_Overlap[i,j]<-sum(ZscoreZhang_Expression_CellType[ZscoreZhang_Expression_CellType$CellType_Primary==names(table(ZscoreZhang_Expression_CellType$CellType_Primary)[i]), 4]%in%ZscoreZhang_Expression_CellType[ZscoreZhang_Expression_CellType$CellType_Primary==names(table(ZscoreZhang_Expression_CellType$CellType_Primary)[j]), 4])/length(ZscoreZhang_Expression_CellType[ZscoreZhang_Expression_CellType$CellType_Primary==names(table(ZscoreZhang_Expression_CellType$CellType_Primary)[i]), 4])
    
  }
}

write.csv(CellTypeSpecificGenes_Master3_Overlap, "CellTypeSpecificGenes_Master3_Overlap.csv")



#What happens if we eliminate overlap between primary categories and then make master indices:

dim(ZscoreZhang_Expression_CellType)
# [1] 2914   31

#Making an empty first row for the storage matrix:
ZscoreZhang_Expression_CellType_NoPrimaryOverlap<-matrix(0, 1, (length(ZscoreZhang_Expression_CellType[1,])))
colnames(ZscoreZhang_Expression_CellType_NoPrimaryOverlap)<-colnames(ZscoreZhang_Expression_CellType)

for(i in 1: length(table(ZscoreZhang_Expression_CellType$CellType_Primary))){
  
  #Choosing all data for a particular primary cell type:
  TempCurrentIndexAllInfo<-ZscoreZhang_Expression_CellType[ZscoreZhang_Expression_CellType$CellType_Primary==names(table(ZscoreZhang_Expression_CellType$CellType_Primary)[i]), ] 
  
  #All of the gene symbols within the current primary cell type:
  TempCurrentIndex<-ZscoreZhang_Expression_CellType[ZscoreZhang_Expression_CellType$CellType_Primary==names(table(ZscoreZhang_Expression_CellType$CellType_Primary)[i]), 4] 
  
  #All of the gene symbols within all other primary cell types:
  TempAllOtherIndices<-ZscoreZhang_Expression_CellType[ZscoreZhang_Expression_CellType$CellType_Primary%in%names(table(ZscoreZhang_Expression_CellType$CellType_Primary)[-i]), 4]
  
  #Grabs only rows of data with gene symbols not found in other primary cell type indices:
  ZscoreZhang_Expression_CellType_NoPrimaryOverlap<-rbind(ZscoreZhang_Expression_CellType_NoPrimaryOverlap, TempCurrentIndexAllInfo[(TempCurrentIndex%in%TempAllOtherIndices)==F,])
  
}

dim(ZscoreZhang_Expression_CellType_NoPrimaryOverlap)
# [1] 2435   31

#removing that one dummy row:
ZscoreZhang_Expression_CellType_NoPrimaryOverlap<-ZscoreZhang_Expression_CellType_NoPrimaryOverlap[-1,]

dim(ZscoreZhang_Expression_CellType_NoPrimaryOverlap)
# [1] 2434   31

write.csv(ZscoreZhang_Expression_CellType_NoPrimaryOverlap, "ZscoreZhang_Expression_CellType_NoPrimaryOverlap.csv")


CellTypeSpecificGenes_Master3_PrimaryTable_NoPrimaryOverlap<-table(ZscoreZhang_Expression_CellType_NoPrimaryOverlap$CellType_Primary)

write.csv(CellTypeSpecificGenes_Master3_PrimaryTable_NoPrimaryOverlap, "CellTypeSpecificGenes_Master3_PrimaryTable_NoPrimaryOverlap.csv")

ZscoreZhang_Expression_CellType_NoPrimaryOverlap_Mean<-matrix(0, nrow=length(table(ZscoreZhang_Expression_CellType_NoPrimaryOverlap$CellType_Primary)), ncol=(length(ZscoreZhang_Expression_CellType_NoPrimaryOverlap[1,])-14))

row.names(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_Mean)<-names(table(ZscoreZhang_Expression_CellType_NoPrimaryOverlap$CellType_Primary))

colnames(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_Mean)<-colnames(ZscoreZhang_Expression_CellType_NoPrimaryOverlap)[-c(1:14)]

#Old version of code:
# for(i in c(15:length(ZscoreZhang_Expression_CellType_NoPrimaryOverlap[1,]))){
# ZscoreZhang_Expression_CellType_NoPrimaryOverlap_Mean[,i-14]<-tapply(ZscoreZhang_Expression_CellType_NoPrimaryOverlap[,i], ZscoreZhang_Expression_CellType_NoPrimaryOverlap$CellType_Primary, mean)
# }


#I went back and changed this so that it averaged by tag first, then by primary cell category

temp<-data.frame(ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag, ZscoreZhang_Expression_CellType_NoPrimaryOverlap$CellType_Primary) 
dim(temp)
# [1] 2434    2

CellTypePrimaryVsTag<-unique(temp)
dim(CellTypePrimaryVsTag)
# [1] 38  2

colnames(CellTypePrimaryVsTag)<-c("Tag","CellType_Primary")
head(CellTypePrimaryVsTag)
#                                    Tag CellType_Primary
# 1     Astrocyte_All_Zhang_PNAS_2015        Astrocyte
# 21     Astrocyte_All_Cahoy_JNeuro_2008        Astrocyte
# 72     Astrocyte_All_Zhang_JNeuro_2014        Astrocyte
# 102      Astrocyte_All_Doyle_Cell_2008        Astrocyte
# 118  Astrocyte_All_Zeisel_Science_2015        Astrocyte
# 309 Endothelial_All_Zhang_PNAS_2015      Endothelial

ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag<-matrix(0, nrow=length(table(ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag)), ncol=(length(ZscoreZhang_Expression_CellType_NoPrimaryOverlap[1,])-14))

row.names(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag)<-names(table(ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag))

colnames(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag)<-colnames(ZscoreZhang_Expression_CellType_NoPrimaryOverlap)[-c(1:14)]

for(i in c(15:length(ZscoreZhang_Expression_CellType_NoPrimaryOverlap[1,]))){
  ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag[,i-14]<-tapply(ZscoreZhang_Expression_CellType_NoPrimaryOverlap[,i], ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag, mean)
}

head(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag)

write.csv(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag, "ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag.csv")


#Making histograms for each cell type tag:

for(i in 1:nrow(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag)){
  png(paste("Histogram_", row.names(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag)[i], ".png", sep=""))
  hist(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag[i,], main=row.names(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag)[i], breaks=16, col=i)
  dev.off()
}

#Plotting each cell type tag vs. Actual cell type:
for(i in 1:nrow(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag)){
  png(paste("CellTypeVs_", row.names(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag)[i], ".png", sep=""), height=500, width=1000)
  plot(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag[i,]~as.factor(ExperimentalCellType), main=row.names(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag)[i], ylab=row.names(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag)[i], xlab="Actual Cell Type", col=i)
  dev.off()
}


                                                            
CorrelationCoefficients_CellTypeIndexVsCellType<-matrix(0, 38, 7)
colnames(CorrelationCoefficients_CellTypeIndexVsCellType)<-names(table(ExperimentalCellType[ExperimentalCellType!="WholeBrain"]))
row.names(CorrelationCoefficients_CellTypeIndexVsCellType)<-row.names(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag)
for(i in c(1:38)){
for(j in c(1:7)){
temp<-summary.lm(lm(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag[i,ExperimentalCellType!="WholeBrain"]~(ExperimentalCellType[ExperimentalCellType!="WholeBrain"]==names(table(ExperimentalCellType[ExperimentalCellType!="WholeBrain"]))[j])))
CorrelationCoefficients_CellTypeIndexVsCellType[i,j]<-temp$r.squared
}
}

write.csv(CorrelationCoefficients_CellTypeIndexVsCellType, "CorrelationCoefficients_CellTypeIndexVsCellType.csv")



#Plotting all the cell type tags for each cell type:
ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag_byActualCellType<-matrix(0, nrow=nrow(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag), ncol=length(unique(ExperimentalCellType)))
row.names(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag_byActualCellType)<-row.names(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag)
colnames(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag_byActualCellType)<-names(table(ExperimentalCellType))

for(i in 1: nrow(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag)){
  ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag_byActualCellType[i,]<-tapply(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag[i,], ExperimentalCellType, mean)
}

write.csv(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag_byActualCellType, "ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag_byActualCellType.csv")

ColorsForCellTypes<-as.character(CellTypePrimaryVsTag[,2])
ColorsForCellTypes[ColorsForCellTypes=="Astrocyte"]<-"lavender"
ColorsForCellTypes[ColorsForCellTypes=="Endothelial"]<-"orange"
ColorsForCellTypes[ColorsForCellTypes=="Microglia"]<-"darkolivegreen3"
ColorsForCellTypes[ColorsForCellTypes=="Mural"]<-"yellow"
ColorsForCellTypes[ColorsForCellTypes=="Neuron_All"]<-"darkviolet"
ColorsForCellTypes[ColorsForCellTypes=="Neuron_Interneuron"]<-"blue"
ColorsForCellTypes[ColorsForCellTypes=="Neuron_Projection"]<-"red"
ColorsForCellTypes[ColorsForCellTypes=="Oligodendrocyte"]<-"lightpink"
ColorsForCellTypes[ColorsForCellTypes=="Oligodendrocyte_Immature"]<-"ivory2"
ColorsForCellTypes[ColorsForCellTypes=="RBC"]<-"indianred4"

#There is a mysterious legend that keeps popping up... that is just the *original* margin parameters. Um. Yeah. Weird, right?

for(i in 1:length(unique(ExperimentalCellType))){
  png(paste("CellTypeTagVs_", colnames(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag_byActualCellType)[i], ".png", sep=""), height=800, width=800)
  barplot(height=ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag_byActualCellType[,i], width=0.3, space = NULL, main=colnames(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag_byActualCellType)[i], ylab="Average Cell Type Index (Mean Z-score for All Genes)", xlab=NULL, names.arg=row.names(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag_byActualCellType), las=2, par(mar=c(25,6,4,2)+0.1), cex.lab=1.3, cex.main=1.3, font.axis=2, font.lab=2, font.main=2,col=ColorsForCellTypes, args.legend = c(x="topright", legend=""))
  dev.off()
}


################################

png("Heatmap_CellType_NoPrimaryOverlap_MeanTag.png", height=1000, width=1000)
heatmap(cor(t(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag[,is.na(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag[1,])==F])), cex.lab=0.3, margins = c(20, 20))
dev.off()

CellType_NoPrimaryOverlap_MeanTag_CorrMatrix<-cor(t(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag[,is.na(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag[1,])==F]))

head(CellType_NoPrimaryOverlap_MeanTag_CorrMatrix)

write.csv(CellType_NoPrimaryOverlap_MeanTag_CorrMatrix, "CellType_NoPrimaryOverlap_MeanTag_CorrMatrix.csv")
###########

temp2<-tapply(ZscoreZhang_Expression_CellType_NoPrimaryOverlap[,15], ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag, mean)
Tag<-names(temp2)

CellTypePrimaryVsTag2<-join(as.data.frame(Tag), as.data.frame(CellTypePrimaryVsTag), by="Tag")


for(i in c(1:length(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag[1,]))){
  ZscoreZhang_Expression_CellType_NoPrimaryOverlap_Mean[,i]<-tapply(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag[,i], CellTypePrimaryVsTag2[,2], mean)
}

head(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_Mean)


write.csv(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_Mean, "ZscoreZhang_Expression_CellType_NoPrimaryOverlap_Mean.csv")

#Making histograms for each primary cell type:

for(i in 1:nrow(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_Mean)){
  png(paste("Histogram_", row.names(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_Mean)[i], ".png", sep=""))
  hist(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_Mean[i,], main=row.names(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_Mean)[i], breaks=16, col=i)
  dev.off()
}

is.numeric(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_Mean)

png("Heatmap_CorMatrixPrimaryCellsNoOverlap.png", height=1000, width=1000)
heatmap(cor(t(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_Mean[,is.na(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_Mean[1,])==F])), cex.lab=0.3, margins = c(20, 20))
dev.off()

CellTypeIndexNoNA3_NormBest_NoPrimaryOverlap_CorrMatrix<-cor(t(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_Mean[,is.na(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_Mean[1,])==F]))

head(CellTypeIndexNoNA3_NormBest_NoPrimaryOverlap_CorrMatrix)

write.csv(CellTypeIndexNoNA3_NormBest_NoPrimaryOverlap_CorrMatrix, "CellTypeIndexNoNA3_NormBest_NoPrimaryOverlap_CorrMatrix.csv")

