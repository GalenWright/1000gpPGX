library(ggplot2)
library(reshape2)
library(plyr)
library(cowplot)

## PERFORM PLOTTING of PCA RESULTS ##

# Read in PCA results generated with smartpca
pca <- read.table("../pgx_analyses/pca/pgx1000GP.VEP.samples.pruned.evec")
# Drop column that is used for case control status
pca <- subset(pca, select = 1:11)
colnames(pca) <- c("sample","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")

# Read in 1000GP sample information
samples <- read.table("../pgx_analyses/freq/integrated_call_samples_v3.20130502.ALL.panel",header = T)

# Merge PCA dataframe with sample data frame based on sample ID
pcaAnnotated <- merge(pca, samples, by.x="sample", by.y="sample")

# Make data frame with unique populations and their superpops
popDetail <- unique(samples[,2:3])
popDetail <- popDetail[with(popDetail, order(popDetail$pop)), ]
populations <- popDetail$pop

# Make COLOUR PALETTES
# Add a column with pre-selected colours for populations
popDetail$popColours <- c("#3A8B8B","#7CFC00","#CD9B1D","#7FFFD4","#FF8C00","#6495ED","#008B8B","#FF0000",
                          "#008B45","#8B4500","#8B7355","#CDCD00","#98FB98","#8B4513","#FFD700","#00FFFF",
                          "#98F5FF","#228B22","#548B54","#8B3A62","#CD6090","#CDC673","#FF69B4","#FFEC8B",
                          "#CD6600","#2E8B57")

# Plot PCA data and colour by superpopulation
superpop.plot.colours <- c("#009E73","#CC79A7", "#56B4E9","#E69F00","#000000")
# Make superpop palette for x axis labels
superpopXcolours <- c(rep("#009E73",7), rep("#CC79A7",4), rep("#56B4E9",5), rep("#E69F00",5), rep("#000000",5))
# Make two colour (colour blind friendly) palette for two colour plots
two.colour.palette <- c("#CC79A7", "#0072B2")

# Calculate variance explained by PCs
pgxPCAeigenvalues <- scan("../pgx_analyses/pca/pgx1000GP.VEP.pruned.eval")
# Percentage variance explained by PC1
(pgxPCAeigenvalues[1])/sum(pgxPCAeigenvalues)*100
# Percentage variance explained by PC2
(pgxPCAeigenvalues[2])/sum(pgxPCAeigenvalues)*100
# Percentage variance explained by PC3
(pgxPCAeigenvalues[3])/sum(pgxPCAeigenvalues)*100
# Percentage variance explained by PC4
(pgxPCAeigenvalues[4])/sum(pgxPCAeigenvalues)*100

# Reorder data frame by super population
popDetail <- popDetail[with(popDetail, order(popDetail$super_pop)), ]
row.names(popDetail) <- NULL
# Order factors for plotting with ggplot
pcaAnnotated$pop <- factor(pcaAnnotated$pop, levels=popDetail$pop[order(popDetail$super_pop)])

# Plot first four principal components, coloured by population
ggplot(data=pcaAnnotated, aes(x=PC1, y=PC2)) +
  geom_point(alpha=0.8, size=2.5, aes(colour = pop)) +
  theme_bw() + xlab("PC 1") + ylab("PC 2") +
  guides(col=guide_legend(ncol = 2)) +
  theme(legend.key=element_blank(), legend.title=element_blank()) +
  scale_colour_manual(values=popDetail$popColours)
ggsave("../pgx_analyses/plots/pgxPopulationsPC1and2.tiff", 
       height=5.36, width=7.58, units='in', dpi=120) 

ggplot(data=pcaAnnotated, aes(x=PC3, y=PC4)) +
  geom_point(alpha=0.8, size=2.5, aes(colour = pop)) +
  theme_bw() + xlab("PC 3") + ylab("PC 4") +
  guides(col=guide_legend(ncol = 2)) +
  theme(legend.key=element_blank(), legend.title=element_blank()) +
  scale_colour_manual(values=popDetail$popColours)
ggsave("../pgx_analyses/plots/pgxPopulationsPC3and4.tiff", 
      height=5.36, width=7.58, units='in', dpi=120)

ggplot(data=pcaAnnotated, aes(x=PC1, y=PC2)) +
  geom_point(alpha=0.8,size=2.5, aes(colour = super_pop)) +
  theme_bw() + xlab("PC 1") + ylab("PC 2") +
  theme(legend.key=element_blank(), legend.title=element_blank()) +
  scale_colour_manual(values=superpop.plot.colours)
ggsave("../pgx_analyses/plots/pgxSuperPopulationsPC1and2.tiff", 
       height=5.36, width=7.58, units='in', dpi=120)

ggplot(data=pcaAnnotated, aes(x=PC3, y=PC4)) +
  geom_point(alpha=0.8, size=2.5, aes(colour = super_pop)) +
  theme_bw() + xlab("PC 3") + ylab("PC 4") +
  theme(legend.key=element_blank(), legend.title=element_blank()) +
  scale_colour_manual(values=superpop.plot.colours)
ggsave("../pgx_analyses/plots/pgxSuperPopulationsPC3and4.tiff", 
       height=5.36, width=7.58, units='in', dpi=120)

# Analyse the functional annotations of variants
functionalSO <- read.table("../pgx_analyses/consequences/SO.terms.MAF.summary", header=T)
# Drop 0 count terms
functionalSO <- functionalSO[functionalSO$countAll!=0,]
# Tidy up terms for plotting
functionalSO$terms <- gsub("_", " ", functionalSO$terms)
functionalSO$terms <- gsub(" variant", "", functionalSO$terms)
functionalSO$terms <- gsub("non coding", "NC", functionalSO$terms)
functionalSO$terms <- Hmisc::capitalize(functionalSO$terms)
# Reorder based on count
functionalSO <- functionalSO[with(functionalSO, order(functionalSO$countAll)), ]
functionalSO$terms <- factor(functionalSO$terms, level=functionalSO[order(functionalSO$countAll),"terms"]) 
rownames(functionalSO) <- NULL

# Plot the total counts for each of the functional classes
functionalSOplot <- ggplot(functionalSO, aes(x = terms, y = countAll) ) +
  geom_point(shape=19) +
  theme_bw() + xlab("Consequence") +  ylab("Count") +
  annotate("text", size=3, fontface="bold", y=functionalSO$countAll+350, 
           x=1:nrow(functionalSO), label=functionalSO$countAll) +
  theme(legend.position="none", axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5, size=12, face=c("bold")))
functionalSOplot
ggsave("../pgx_analyses/plots/pgxFunctionalClassesCounts.tiff", 
       height=5.36*0.75, width=7.58, units='in', dpi=120)

# Make Frequency dataframe
functionalSOFreq <- NULL
functionalSOFreq$terms <- functionalSO$terms
functionalSOFreq$freqAll <- functionalSO$countAll/sum(functionalSO$countAll)
functionalSOFreq$freqSingletons <- functionalSO$countSingletons/sum(functionalSO$countSingletons)
functionalSOFreq$freqnonSing.0.01MAF <- functionalSO$countnonSing.0.01MAF/sum(functionalSO$countnonSing.0.01MAF)
functionalSOFreq$freq0.01.0.05MAF <- functionalSO$count0.01.0.05MAF/sum(functionalSO$count0.01.0.05MAF)
functionalSOFreq$freqGreater0.05 <- functionalSO$countGreater0.05/sum(functionalSO$countGreater0.05)
functionalSOFreq <- as.data.frame(functionalSOFreq)
# Add spaces in names for better plotting
colnames(functionalSOFreq) <- c("terms", "freqAll", "Singletons",
                                "Singletons - MAF 0.01", "MAF 0.01 - 0.05","MAF greater than 0.05")

# Test for differences between the different counts in the different frequency classes
# Make a dataframe to put significance test results in a set a count for indexing
functionalSignificanceDF <- data.frame(ncol=2, nrow=length(unique(functionalSO$term)))
count <- 1 
for (functional in unique(functionalSO$term)){
  sumTerms <- colSums(functionalSO[,3:6])
  functionalCount <- functionalSO[,3:6][functionalSO$term==functional,]
  remainderCount <- sumTerms-functionalCount
  functionalMatrix <- as.matrix(rbind(functionalCount, remainderCount))
  # increase workspace to prevent fisher test from crashing
  fisherFunctional <- fisher.test(functionalMatrix,workspace=100000000)
  functionalSignificanceDF[count,1] <- functional
  functionalSignificanceDF[count,2] <- fisherFunctional$p.value
  count <- count +1
}

# Name columns
colnames(functionalSignificanceDF) <- c("terms","fisherPvalue")
# Perform Bonferonni adjustmet
functionalSignificanceDF$fisherPvalueBonferroni <- p.adjust(functionalSignificanceDF$fisherPvalue, method="bonferroni", n=length(functionalSignificanceDF$fisherPvalue))
# See which classes differ significantly after multiple testing correction
functionalSignificanceDF[functionalSignificanceDF$fisherPvalueBonferroni<0.05,]

# Append p-values to frequency dataframe
functionalSOFreq <- merge(functionalSOFreq, functionalSignificanceDF, by="terms")

####### Order dataframe and make a index for plotting
# Add index for later plotting
functionalSignificanceDF$index <- 1:nrow(functionalSignificanceDF)
significantLabels <- functionalSignificanceDF[functionalSignificanceDF$fisherPvalueBonferroni<=0.05,]
significantLabels$fisherPvalueBonferroni <- format(significantLabels$fisherPvalueBonferroni,digits=3)

# Melt and reorder for plotting
mfunctionalSOFreq <- melt(functionalSOFreq, id=c("terms","fisherPvalueBonferroni","fisherPvalue" ))
# Remove combined allele frequency 
mfunctionalSOFreq <- mfunctionalSOFreq[mfunctionalSOFreq$variable!="freqAll",]
# Plot consequences that occur at a frequency of at least 1% in total
functionalFreqPlot <- ggplot(mfunctionalSOFreq, aes(x = terms, y = value) ) +
  geom_point(aes(shape=variable, col=variable), size=4) +
  theme_bw() + xlab("Consequence") +  ylab("Proportion") +
  scale_shape_manual(values=c(15:18)) +
  coord_cartesian(ylim = c(0.00, 0.45)) +
  scale_y_continuous(limits = c(0.00, 0.45)) +
  annotate("text", size=3, fontface="bold", y=0.12, angle=90,
           x=significantLabels$index, label=significantLabels$fisherPvalueBonferroni) +
  theme(legend.title=element_blank(), legend.position="bottom",legend.key=element_blank(),
        axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5, size=12, face=c("bold")))
functionalFreqPlot
ggsave("../pgx_analyses/plots/pgxFunctionalClassesByFrequency.tiff", 
             height=5.36*1.25, width=7.58, units='in', dpi=120)

plot_grid(functionalSOplot, functionalFreqPlot, 
                                  labels = c("A", "B"), align = "v", nrow=2)
ggsave("../pgx_analyses/plots/pgxFunctionalClassesGrid.tiff", 
             height=5.36*2, width=7.58, units='in', dpi=120)

# Load Global allele frequency information
setwd("../pgx_analyses/freq")
myFiles <- list.files(path = ".", pattern="*.Ind.frq")

for (file in myFiles){
  # if the merged dataset doesn't exist, create it
  if (!exists("dataset")){
    #Still need to remove first population data
    dataset <- read.table(file, header=F, sep="\t", skip=1)
    dataset <- dataset[,1:2]
  }
  # if the merged dataset does exist, append to it
  if (exists("dataset")){
    temp_dataset <-read.table(file, header=F, sep="\t", skip=1)
    dataset<-cbind(dataset, temp_dataset[,6])
    rm(temp_dataset)
  }
}
setwd("../../scripts/")
popColnames <-  c("CHROM", "POS")

for (column in (length(dataset)-2)) {
  for (i in 1:length(populations)){
    popFreq <- paste0(populations[i], "_Alt_FREQ")
    popColnames <- append(popColnames, popFreq) 
  }
}
colnames(dataset) <- popColnames

# Upload info and rsID file
INFOandrsID <- read.table("../pgx_analyses/freq/pgx1000GP.VEP.biall.INFO.rsID", header=T)
# Combine annotated with dataset, removing duplicated first two columns of dataset
datasetAnnotated <- cbind(INFOandrsID, dataset[,-(1:2)])

# If mean allele frequency is >0.5 flip alleles in those columns
for (row in 1:nrow(datasetAnnotated)) {
    if (datasetAnnotated$AF[row]>0.5) {
      datasetAnnotated[row,grep('*FREQ',names(datasetAnnotated))] <- 
        lapply(datasetAnnotated[row,grep('*FREQ',names(datasetAnnotated))], function(x) 1-x)
         }
}

### GET NUMBER of POLYMORPHIC ALLELES IN EACH POPULATION
datasetFreq <- datasetAnnotated[,grep('*FREQ',names(datasetAnnotated))]
polymorphic <- colSums(datasetFreq != 0)
polymorphic

### GET NUMBER of SINGLETON ALLELES IN EACH POPULATION

# Create unique ID column for matching
datasetAnnotated$UniqueID <- paste(datasetAnnotated$CHROM, datasetAnnotated$POS, datasetAnnotated$ALT, sep=":")

# Add a column indicating whether the variant is a singleton or not
datasetAnnotated$Singleton <- with(datasetAnnotated, ifelse(AC==1, "Yes","No"))
# Add a column indicating whether the variant has a global MAF>=0.01
datasetAnnotated$Maf0.01 <- with(datasetAnnotated, ifelse(AF>=0.1, "Yes","No"))

# Make a tall dataset for plotting
# Remove annotated columns and unique ID one
mDatasetFreq <- melt(datasetAnnotated[,-c(1:8,35)], id.vars=c("Singleton","Maf0.01"))
# Add populations and super populations
mDatasetFreq$Pop <- rep(sort(popDetail$pop),each=nrow(datasetAnnotated))
mDatasetFreq$Super <- rep(popDetail$super_pop[order(popDetail$pop)],each=nrow(datasetAnnotated))
# Order levels of populations based on superpopulations
mDatasetFreq$Pop <- factor(mDatasetFreq$Pop,levels=popDetail$pop[order(popDetail$super_pop)])

# Add CADD scores
cadd <- read.table("../data/pgx1000GP.VEP.cadd.annotations.tsv")
colnames(cadd) <- c("CHROM", "POS", "REF", "ALT", "RawScore", "PHRED")
cadd$UniqueID <- paste(cadd$CHROM, cadd$POS, cadd$ALT, sep=":")
cadd <- merge(cadd, datasetAnnotated, by="UniqueID")
cadd.delterious <- cadd[cadd$PHRED>=20,]
summary(cadd$PHRED[cadd$Singleton=="Yes"])
summary(cadd$PHRED[cadd$Singleton=="No"])
summary(cadd$PHRED[cadd$AF>=0.005])
summary(cadd$PHRED[cadd$AF<0.005])
cadd$Rare <- cadd$AF<=0.005
t.test(cadd$PHRED ~ cadd$Rare) ### So different that it goes to min R p value 2.2e-16

## ANALYSE NUMBER of SINGLETONS/POPULATION and PLOT TOTAL POLYMORPHIC/POPULATION
datasetFreqSing <- datasetAnnotated[datasetAnnotated$Singleton=="Yes",c(grep('*_Alt_FREQ',names(datasetAnnotated)))]
# Count number of non-zero occurences in each
popSingles <- apply(datasetFreqSing, 2, function(x) sum(x!=0))
polymorphicTable <- cbind(polymorphic, popSingles)
rownames(polymorphicTable) <- populations
colnames(polymorphicTable) <- c("Polymorphic", "Singletons")
polymorphicTable <- as.data.frame(polymorphicTable)
polymorphicTable$NonSingleton <- polymorphicTable$Polymorphic - polymorphicTable$Singletons
polymorphicTable$Populations <- populations
polymorphicTable$Polymorphic <- NULL
polymorphicTableM <- melt(polymorphicTable, id.vars = "Populations")
# Reorder for plotting
polymorphicTableM$Populations <- factor(polymorphicTableM$Populations,levels=popDetail$pop[order(popDetail$super_pop)])

polymorphicTableM$variable <- gsub("NonSingleton", "Non-singletons", polymorphicTableM$variable)
# Plot singletons vs non-singletons per population using a stacked bar graph
ggplot(polymorphicTableM, aes(x=Populations, y=value, fill=variable))  + 
  geom_bar(stat="identity") + xlab("Population") + ylab("Count") + theme_bw() +
  theme(legend.title=element_blank(),axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=12,face="bold", 
   colour=superpopXcolours)) + scale_fill_manual(values=two.colour.palette)
ggsave("../pgx_analyses/plots/pgxPolymorphicPerPopulation.tiff", 
       height=5.36, width=7.58, units='in', dpi=120) 

# Calculate number of singletons / individual / population
popSinglesCount <- read.table("../pgx_analyses/freq/phase3.singletons.per.sample", header =T)
popSinglesCount <- popSinglesCount[,-c(3,4)]

popSinglesSummary <- ddply(popSinglesCount, "pop", summarise,
                           mean = mean(singletons), sd = sd(singletons),
                           sem = sd(singletons)/sqrt(length(singletons)))

popSinglesSummary$pop <- factor(popSinglesSummary$pop,levels=popDetail$pop[order(popDetail$super_pop)])

ggplot(popSinglesSummary, aes(x=pop, y=mean))  + 
  geom_bar(stat="identity", fill="#009E73") + xlab("Population") + 
  ylab("Count") + theme_bw() +
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem)) +
  theme(legend.title=element_blank(), 
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=12,
                                 face="bold", colour=superpopXcolours))
ggsave("../pgx_analyses/plots/pgxAverageSingletonsPerPopulation.tiff", 
       height=5.36, width=7.58, units='in', dpi=120)
  
# Mean total number of polymorphisms across all populations
mean(polymorphicTable$Singletons+polymorphicTable$NonSingleton)

## ANALYSE THE CLINICALLY RELEVANT VARIANTS
rsIDGene <- read.table("../data/pgxClinicalLevel1.txt", header=F, sep="\t")
colnames(rsIDGene) <- c("ID", "Gene")
clinDataset <- merge(datasetAnnotated, rsIDGene, by="ID")

# Select all columns with ref allele frequency
RefAlleleFreqClin <- clinDataset[,c(grep('*_Alt_FREQ',names(clinDataset)))]
RefAlleleFreqClin$ID <- clinDataset$ID
RefAlleleFreqClin$Gene <- clinDataset$Gene
RefAlleleFreqClin <- RefAlleleFreqClin[order(RefAlleleFreqClin$Gene),]

# Plot allele frequencies
# Make tall data frame for plotting in ggplot
mRefAlleleFreqClin <- melt(RefAlleleFreqClin, id=c("ID","Gene"))
mRefAlleleFreqClin$Pop <- rep(sort(popDetail$pop),each=nrow(RefAlleleFreqClin))
mRefAlleleFreqClin$Super <- rep(popDetail$super_pop[order(popDetail$pop)],each=nrow(RefAlleleFreqClin))
mRefAlleleFreqClin$GeneRs <- paste(mRefAlleleFreqClin$Gene, mRefAlleleFreqClin$ID, sep=" ")

# Plot allele frequency stripchart
# Allele frequencies of pharmacogenomic variants with high levels of clinical evidence
clinFreqPlot <- ggplot(mRefAlleleFreqClin, aes(x = GeneRs, y = value) ) +
  geom_point( aes(colour = Super), shape=19) +
  theme_bw() + xlab("Variant") + 
  scale_colour_manual(values=superpop.plot.colours) +
  ylab("Frequency") +
  theme(legend.title=element_blank(), 
        axis.text.x=element_text(angle = 90, hjust = 1,
                                 vjust = 0.5, size=12, face="bold"))
clinFreqPlot
ggsave("../pgx_analyses/plots/pgxClinicalEvidencePopulation.tiff", 
       height=5.36, width=7.58, units='in', dpi=120) 

# Run Fst analyses and plot in ggplot
weightFst <- read.table("../pgx_analyses/fst/weighted.fst.estimates.txt")
colnames(weightFst) <- c("pop1", "pop2", "Fst")
weightFst$FstPercent <- weightFst$Fst*100
FstData <- matrix(weightFst$FstPercent,nrow=26,ncol=26, byrow = TRUE)
colnames(FstData) <- populations
rownames(FstData) <- populations
FstData <- as.data.frame(FstData)

# Order population Fst by superpopulation
FstData <- FstData[order(popDetail$super_pop),]
FstData <- FstData[,order(popDetail$super_pop)]
FstMatrix <-data.matrix(FstData)

# Calculate summary of FST x 100 
mean(FstMatrix,na.rm =T)
max(FstMatrix,na.rm =T)
min(FstMatrix,na.rm =T)
# Calculate Mean FST x 100 AFR 
mean(FstMatrix[1:7,1:7],na.rm =T)
# Calculate Mean FST x 100 AMR
mean(FstMatrix[8:11,8:11],na.rm =T)
# Calculate Mean FST x 100  ASN
mean(FstMatrix[12:16,12:16],na.rm =T)
# Calculate Mean FST x 100 EUR
mean(FstMatrix[17:21,17:21],na.rm =T)
# Calculate Mean FST x 100 SAN
mean(FstMatrix[22:26,22:26],na.rm =T)

# Analyse highly differentiated polymorphisms
# Upload the variant annotation details
pgxVariants <- read.table("../pgx_analyses/consequences/pgx1000GP.VEP.CSQ.annotation.INFO", header=T, fill=T)

# Upload FST analyses
AFR.EAS.Fst <- read.table("../pgx_analyses/fst/AFR_EAS.super.fst.txt.weir.fst", header=T)
AFR.EUR.Fst <- read.table("../pgx_analyses/fst/AFR_EUR.super.fst.txt.weir.fst", header=T)
AFR.SAS.Fst <- read.table("../pgx_analyses/fst/AFR_SAS.super.fst.txt.weir.fst", header=T)
EUR.EAS.Fst <- read.table("../pgx_analyses/fst/EUR_EAS.super.fst.txt.weir.fst", header=T)
EUR.SAS.Fst <- read.table("../pgx_analyses/fst/EUR_SAS.super.fst.txt.weir.fst", header=T)
SAS.EAS.Fst <- read.table("../pgx_analyses/fst/SAS_EAS.super.fst.txt.weir.fst", header=T)
AMR.AFR.Fst <- read.table("../pgx_analyses/fst/AMR_AFR.super.fst.txt.weir.fst", header=T)
AMR.EAS.Fst <- read.table("../pgx_analyses/fst/AMR_EAS.super.fst.txt.weir.fst", header=T)
AMR.EUR.Fst <- read.table("../pgx_analyses/fst/AMR_EUR.super.fst.txt.weir.fst", header=T)
AMR.SAS.Fst <- read.table("../pgx_analyses/fst/AMR_SAS.super.fst.txt.weir.fst", header=T)

# Add Fst values from each analysis to main data frame
pgxVariants$AFR.EAS.Fst <- AFR.EAS.Fst$WEIR_AND_COCKERHAM_FST
pgxVariants$AFR.EUR.Fst <- AFR.EUR.Fst$WEIR_AND_COCKERHAM_FST
pgxVariants$AFR.SAS.Fst <- AFR.SAS.Fst$WEIR_AND_COCKERHAM_FST
pgxVariants$EUR.EAS.Fst <- EUR.EAS.Fst$WEIR_AND_COCKERHAM_FST
pgxVariants$EUR.SAS.Fst <- EUR.SAS.Fst$WEIR_AND_COCKERHAM_FST
pgxVariants$SAS.EAS.Fst <- SAS.EAS.Fst$WEIR_AND_COCKERHAM_FST
pgxVariants$AMR.AFR.Fst <- AMR.AFR.Fst$WEIR_AND_COCKERHAM_FST
pgxVariants$AMR.EAS.Fst <- AMR.EAS.Fst$WEIR_AND_COCKERHAM_FST
pgxVariants$AMR.EUR.Fst <- AMR.EUR.Fst$WEIR_AND_COCKERHAM_FST
pgxVariants$AMR.SAS.Fst <- AMR.SAS.Fst$WEIR_AND_COCKERHAM_FST

# Replace NAs with 0
pgxVariants[is.na(pgxVariants)] <- 0

# Select columns where at least one of the Fst comparisons 
# is highly differentiated (i.e. >0.5)
pgxVariantsHighFst <- pgxVariants[pgxVariants$AFR.EAS.Fst>0.5 |
                                    pgxVariants$AFR.EUR.Fst>0.5 |
                                    pgxVariants$AFR.SAS.Fst>0.5 |
                                    pgxVariants$EUR.EAS.Fst>0.5 |
                                    pgxVariants$EUR.SAS.Fst>0.5 |
                                    pgxVariants$SAS.EAS.Fst>0.5 |
                                    pgxVariants$AMR.AFR.Fst>0.5 |
                                    pgxVariants$AMR.EAS.Fst>0.5 |
                                    pgxVariants$AMR.EUR.Fst>0.5 |
                                    pgxVariants$AMR.SAS.Fst>0.5, ]

# Calculate average Fst for each variant
pgxVariantsHighFst$Fst.mean <- rowMeans(pgxVariantsHighFst[,grep('*Fst',names(pgxVariantsHighFst))])
row.names(pgxVariantsHighFst) <- NULL
# Remove function aliases from those that have them
pgxVariantsHighFst$FUNCTION <- gsub("&.*","",pgxVariantsHighFst$FUNCTION)

# Select rows with the maximum mean FST per gene
pgxVariantsHighFstMax <- ddply(pgxVariantsHighFst, .(GENE), function(x) x[which.max(x$Fst.mean),])
# Remove multiallelic
pgxVariantsHighFstMax <- pgxVariantsHighFstMax[-grep(",", pgxVariantsHighFstMax$ALT),]
# Add CADD scores
pgxVariantsHighFstMax <- merge(pgxVariantsHighFstMax, cadd[,c(7,15)])
pgxVariantsHighFstMax <- pgxVariantsHighFstMax[,-c(2:7)]
pgxVariantsHighFstMax <- pgxVariantsHighFstMax[order(pgxVariantsHighFstMax$GENE),]
write.csv(pgxVariantsHighFstMax, "../pgx_analyses/fst/highFstMaxGenes.csv", row.names = FALSE)

# Add frequencies for these polymorphisms for plotting
highFstDataset <- merge(pgxVariantsHighFstMax, datasetAnnotated, by="ID")

# Select all columns with ref allele frequency
RefhighFstDataset <- highFstDataset[,c(grep('*_Alt_FREQ',names(highFstDataset)))]
RefhighFstDataset$ID <- highFstDataset$ID
RefhighFstDataset$Gene <- highFstDataset$GENE
RefhighFstDataset$Function <- as.factor(highFstDataset$FUNCTION)
RefhighFstDataset$Fst.mean <- highFstDataset$Fst.mean
RefhighFstDataset$GeneRs <- paste(RefhighFstDataset$Gene, RefhighFstDataset$ID, sep=" ")
RefhighFstDataset <- RefhighFstDataset[order(RefhighFstDataset$Gene),]

# Plot allele frequencies
# Make tall data frame for plotting in ggplot
mRefhighFstDataset <- melt(RefhighFstDataset, id=c("ID","Gene","Function","Fst.mean","GeneRs"))
mRefhighFstDataset$Pop <- rep(sort(popDetail$pop),each=nrow(RefhighFstDataset))
mRefhighFstDataset$Super <- rep(popDetail$super_pop[order(popDetail$pop)],each=nrow(highFstDataset))

# Reorder factors for later plotting in ggplot from max mean FST
mRefhighFstDataset$GeneRs <- factor(mRefhighFstDataset$GeneRs,levels=RefhighFstDataset[order(RefhighFstDataset$Fst.mean),"GeneRs"])

# Plot allele frequency stripchart
# Allele frequencies of highly differentiated pharmacogenomic variants
ggplot(mRefhighFstDataset, aes(x = GeneRs, y = value) ) +
  geom_point( aes(colour = Super), shape= 19) +
  theme_bw() + xlab("Variant") + 
  scale_colour_manual(values=superpop.plot.colours) +
  ylab("Frequency") + 
  theme(legend.title=element_blank(), 
        axis.text.x=element_text(angle = 90, hjust = 1,
                                 vjust = 0.5, size=12, face="bold"))
ggsave("../pgx_analyses/plots/pgxHighDiffPGxPlot.tiff", 
      height=5.36, width=7.58, units='in', dpi=120) 

# Make FST heatmap with ggplot2
FstDataTall <- data.frame(Pop1 = rep(colnames(FstData), each = nrow(FstData)),
                        Pop2 = rownames(FstData),
                        FST = unlist(FstData))

# Order Factors for ggplot to group according to population
FstDataTall$Pop1 <- factor(FstDataTall$Pop1,levels=popDetail$pop[order(popDetail$super_pop)])
FstDataTall$Pop2 <- factor(FstDataTall$Pop2,levels=popDetail$pop[order(popDetail$super_pop)])

#Matrix of Weir and Cockerham weighted Fst estimates (x100)
ggplot(FstDataTall, aes(x = Pop2, y = Pop1)) + geom_tile(aes(fill = FST),
colour = "white") + scale_fill_gradient(low = "white", high = "steelblue") +
  xlab("Population") + ylab("Population") + 
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, 
      size=12, face="bold",colour=superpopXcolours),
      axis.text.y = element_text(size=12, face="bold", colour=superpopXcolours))
ggsave("../pgx_analyses/plots/pgxPopulationFstMatrix.tiff", 
      height=6.7, width=9.95, units='in', dpi=120) 

# ANALYSE LoF VARIANTS
# Read in the data from vcftools --012 for all individuals
# "Genotypes are represented as 0, 1 and 2, 
# where the numbers represent that number of non-reference alleles"
# Missing genotypes are represented by -1
snpScoresLof <- read.table('../pgx_analyses/lof/LofPerSample.012')
# Remove the first column which is the index of sample
snpScoresLof <- snpScoresLof[,-1]
snpScoresLof$totalLof <- apply(snpScoresLof,1,sum)

# Read in the individuals sheet which represents 
# This file represents the individuals who were analysed
individualsLof <- read.table('../pgx_analyses/lof/LofPerSample.012.indv')

# Make a data frame with individuals and total LoF they carry
allLof <- cbind(individualsLof, snpScoresLof$totalLof)
colnames(allLof) <- c("sample", "totalLof")

# ADD POPULATION DATA to the DATAFRAME
allLofPop <- merge(samples, allLof, by="sample")

# MAKE PLOTS in GGPLOT
# Order for plotting in ggplot
allLofPop$pop <- factor(allLofPop$pop,levels=popDetail$pop[order(popDetail$super_pop)])

# Number of individuals carrying at least one LoF variant
lofSuper <- table(allLofPop$totalLof, allLofPop$super_pop)
prop.table(lofSuper,2)

# Test for differences in means using Kruskal-Wallis rank sum test
kruskal.test(totalLof ~ super_pop, data = allLofPop) 

# Plot number of LoF variants per gene and their frequency
lofConsequences <- read.table("../pgx_analyses/lof/LoFConsequencesIDacGENE")
colnames(lofConsequences) <- c("UniqueID", "AC", "EAS_AF", "AMR_AF", "AFR_AF", "EUR_AF", "SAS_AF","Gene")

lofACperGene <- ddply(lofConsequences,c("Gene"),summarise, N=sum(AC), LOF=length(AC), 
                      EAS_combined_AF=sum(EAS_AF), AMR_combined_AF=sum(AMR_AF), AFR_combined_AF=sum(AFR_AF),
                      EUR_combined_AF=sum(EUR_AF), SAS_combined_AF=sum(SAS_AF))
lofACperGene$GeneVariants <- with(lofACperGene, ifelse(lofACperGene$LOF>1 ,paste(lofACperGene$Gene, " (", lofACperGene$LOF, " variants)", sep=""),
                                  paste(lofACperGene$Gene, " (", lofACperGene$LOF, " variant)", sep="")))

# Reorder factors for plotting with ggplot
lofACperGene$GeneVariants <- factor(lofACperGene$GeneVariants, levels=lofACperGene[order(lofACperGene$N),"GeneVariants"])

# Bubble plot of AC/no. of LoF variants in ggplot
lofACplot <- ggplot(lofACperGene, aes(x = GeneVariants, y = N) ) +
  geom_point(shape=19, aes(size=LOF)) +
  theme_bw() + xlab("Gene") + ylab("Combined allele count") + 
  theme(legend.position="none", axis.text.x=element_text(angle=90,
        hjust = 1, vjust = 0.5, size=12, face=c("bold")))
lofACplot
ggsave("../pgx_analyses/plots/pgxLoFPerGeneCount.tiff", 
      height=5.36, width=7.58, units='in', dpi=120) 

# Allele frequencies of LoF variants in different super populations
# Melt and reorder for plotting
mLofACperGene <- melt(lofACperGene[,c(2,4:9)], id.vars=c("GeneVariants","N"))
mLofACperGene$variable <- gsub("_combined_AF","", mLofACperGene$variable)
#Plot combined allele frequency
lofAFplot <- ggplot(mLofACperGene, aes(x = GeneVariants, y = value)) +
  geom_point(shape=19, aes(colour=variable)) +
  theme_bw() + xlab("Gene") + ylab("Combined allele frequency") + 
  theme(legend.title=element_blank(), legend.position="bottom",
        axis.text.x=element_text(angle=90,
        hjust = 1, vjust = 0.5, size=12, face=c("bold")))
lofAFplot
ggsave("../pgx_analyses/plots/pgxLoFPerCombinedAF.tiff", 
              height=5.36, width=7.58, units='in', dpi=120) 
# SLC19A1 in AMR is driven mainly by PEL individuals (rs200964006)

plot_grid(lofACplot, lofAFplot, 
          labels = c("A", "B"), align = "v", nrow=2)
ggsave("../pgx_analyses/plots/pgxLoFPerGene.tiff", 
       height=5.36*2, width=7.58*2, units='in', dpi=120)

# ANALYSE CLINICAL VARIANTS
snpScoresClin <- read.table('../pgx_analyses/clinical/ClinPerSample.012')
# Remove the first column which is the index of sample
snpScoresClin <- snpScoresClin[,-1]

# Flip alleles for those where the non-reference allele is the major allele
# Recode as "a" and "b" first to avoid ambiguity
for (column in 1:ncol(snpScoresClin)) {
  if (colSums(snpScoresClin[column]) > nrow(snpScoresClin)) {
    snpScoresClin[,column][snpScoresClin[,column]==2] <- "a"
    snpScoresClin[,column][snpScoresClin[,column]==0] <- "b"
    snpScoresClin[,column][snpScoresClin[,column]=="a"] <- 0
    snpScoresClin[,column][snpScoresClin[,column]=="b"] <- 2
    snpScoresClin[,column] <- as.numeric(snpScoresClin[,column])
  }
}

snpScoresClin$totalClin <- apply(snpScoresClin,1,sum)
# Upload individuals
individualsClin <- read.table('../pgx_analyses/clinical/ClinPerSample.012.indv')
# Make a data frame with individuals and total clinical PGx vaiants they carry
allClin <- cbind(individualsClin, snpScoresClin$totalClin)
colnames(allClin) <- c("sample", "totalClin")
allLofClinPop <- merge(allLofPop, allClin, by="sample")
summary(allLofClinPop$totalClin)
clinPopTable <- table(allLofClinPop$totalClin,allLofClinPop$pop)
clinPopTable
write.csv(clinPopTable,"../pgx_analyses/clinical/clinPopTable.csv")

# PLOT CLINICAL VARIANTS / SAMPLE
# Violin and jitter plot: Number of clinically-relevant variants/sample/population
clinPerSamplePlot <- ggplot(allLofClinPop, aes(x=pop, y=totalClin)) + 
  geom_violin(aes(colour=super_pop)) + theme_bw() +
  theme(legend.title=element_blank(), 
        axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5, size=12, face="bold")) +
  xlab("Super population") + ylab("Count") + 
  scale_colour_manual(values=superpop.plot.colours) +
  geom_jitter(alpha=0.5, aes(color=super_pop),
               position = position_jitter(width = 0.1))
clinPerSamplePlot
ggsave("../pgx_analyses/plots/pgxClinPerSamplePerPopulation.tiff", 
       height=5.36, width=7.58, units='in', dpi=120) 

plot_grid(clinFreqPlot, clinPerSamplePlot, labels = c("A", "B"), align = "v", nrow=2)
ggsave("../pgx_analyses/plots/pgxClinicalGrid.tiff", 
       height=5.36*2, width=7.58, units='in', dpi=120)

# Compare carriers to non-carriers of LoF variants
allLofClinPop$LofCarry <- with(allLofClinPop, ifelse(totalLof==0, "No","Yes"))
allLofClinPop$ClinCarry <- with(allLofClinPop, ifelse(totalClin==0, "No","Yes"))

tableLoFcarry <- table(allLofClinPop$LofCarry, allLofClinPop$super_pop)
prop.table(tableLoFcarry,2)
table(allLofClinPop$LofCarry)

# Analyse the amount of variation per gene 
# Add rare variant data
pgxGeneSummary <- read.csv("../pgx_analyses/consequences/genes/geneSummaryVariantsMask.csv", header = TRUE)
rareMissense <- read.table("../pgx_analyses/consequences/genes/rare.missense.gene.txt", header = T)
pgxGeneSummary <- merge(pgxGeneSummary, rareMissense, by = "gene")
pgxGeneSummary$percRareMiss <- pgxGeneSummary$rare_missense/pgxGeneSummary$missense_variant

# Select top five and bottom 5 per missense rate
pgxGeneSummary <- pgxGeneSummary[order(pgxGeneSummary$missensePerCDSrank),]
mostMissenseGenes <- pgxGeneSummary[1:5,]
leastMissenseGenes <- pgxGeneSummary[(nrow(pgxGeneSummary)-4):nrow(pgxGeneSummary),]

# Select top five and bottom 5 per polymorphic
pgxGeneSummary <- pgxGeneSummary[order(pgxGeneSummary$variantsPerLengthRank),]
mostPolyGenes <- pgxGeneSummary[1:5,]
leastPolyGenes <- pgxGeneSummary[(nrow(pgxGeneSummary)-4):nrow(pgxGeneSummary),]

# Look at most polymorphic/missense
mostPolyMissense <- rbind(mostPolyGenes,mostMissenseGenes)
# Drop duplicates
mostPolyMissense <- mostPolyMissense[!duplicated(mostPolyMissense),]
# Order by missense rate
mostPolyMissense <- mostPolyMissense[order(mostPolyMissense$missensePerCDSrank),]

# Look at least polymorphic/missense
leastPolyMissense <- rbind(leastPolyGenes,leastMissenseGenes)
# Drop duplicates
leastPolyMissense <- leastPolyMissense[!duplicated(leastPolyMissense),]
# Order by missense rate
leastPolyMissense <- leastPolyMissense[order(leastPolyMissense$missensePerCDSrank),]

# Most inaccess 
pgxGeneSummary <- pgxGeneSummary[order(-pgxGeneSummary$percentageInaccess),]
mostInaccess <- pgxGeneSummary[1:5,]

# Analyse CYP genes in segmental/inaccessible loci
# Segmental
segmentalGenes <- pgxGeneSummary$gene[pgxGeneSummary$segmental_variants>0]
nonSegmentalGenes <- pgxGeneSummary$gene[pgxGeneSummary$segmental_variants==0]

CYP <- c(length((grep("^CYP", segmentalGenes, invert = F))), length((grep("^CYP", nonSegmentalGenes, invert = F))))
nonCYP <- c(length((grep("^CYP", segmentalGenes, invert = T))), length((grep("^CYP", nonSegmentalGenes, invert = T))))
seg <- cbind(CYP, nonCYP)
fisher.test(seg)

# Plot highly differentiated SNPs (5% in one population, but 0.5% in the global population)
highDiffPop <- read.table("../pgx_analyses/popdiff/population.pop.diff.snps.genotypes", header=F)
colnames(highDiffPop) <- c("Population", "ID")
highDiffGenes <- read.table("../pgx_analyses/popdiff/genes.pop.diff.snps.genotypes", header=F)
colnames(highDiffGenes) <- c("Gene", "Function", "ID", "UniqueID")
highDiffGenes$Function <- gsub("_", " ", highDiffGenes$Function)

# Write table of population differentiated SNPs
highDiff <- merge(highDiffPop, highDiffGenes, by="ID")
highDiff <- merge(highDiff, cadd[,c(1,7)], by = "UniqueID")
highDiff <- highDiff[order(highDiff$Gene),]
write.csv(highDiff, "../pgx_analyses/popdiff/highDiffSNPs.csv", row.names = FALSE, quote = FALSE)

# Make dataset for plotting
highDiffrsIDFreq <- merge(highDiff, datasetAnnotated, by="ID")
highDiffrsIDFreq$Gene <- paste("(", highDiffrsIDFreq$Gene, ":", sep ="")
highDiffrsIDFreq$ID <- paste(highDiffrsIDFreq$ID, ")", sep ="")
highDiffrsIDFreq$Label <- paste(highDiffrsIDFreq$Population, highDiffrsIDFreq$Gene, highDiffrsIDFreq$ID, sep = " ")
highDiffrsIDFreq <- highDiffrsIDFreq[,c(43,grep('*_Alt_FREQ',names(highDiffrsIDFreq)))]

# Make tall data frame for plotting in ggplot
mHighDiffrsIDFreq <- melt(highDiffrsIDFreq, id="Label")
mHighDiffrsIDFreq$Pop <- rep(sort(popDetail$pop),each=nrow(highDiffrsIDFreq))
mHighDiffrsIDFreq$Super <- rep(popDetail$super_pop[order(popDetail$pop)],each=nrow(highDiffrsIDFreq))

# Plot allele frequency stripchart
# Allele frequencies of highly differentiated pharmacogenomic variants
highDiffPlot <- ggplot(mHighDiffrsIDFreq, aes(x = Label, y = value) ) +
  geom_point( aes(colour = Super), shape=19) +
  theme_bw() + xlab("Variant") + 
  scale_colour_manual(values=superpop.plot.colours) +
  ylab("Frequency") + xlab("Population and variant") + geom_hline(aes(yintercept=0.05), linetype="dashed") +
  theme(legend.title=element_blank(), 
       axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5, size=12, face="bold"))
highDiffPlot
ggsave("../pgx_analyses/plots/highDiffPlot.tiff", height=5.36, width=7.58, units='in', dpi=120)

# COMPARE ALLELE FREQUENCIES WITH OTHER STUDIES
# Look at global allele frequencies for 1000GP and ExAC
ExAC.1000GP <- read.table("../pgx_analyses/freq_comp/freq.comp.pos.1000gp.exac.af", header=T)
ExAC.1000GP$freq_1000gp <- as.numeric(as.character(ExAC.1000GP$freq_1000GP))
ExAC.1000GP$diff <- abs(ExAC.1000GP$freq_1000GP - ExAC.1000GP$freq_exac)

# Plot results 
# Calculate r2 value
exacR2 <- round((summary(lm(ExAC.1000GP$freq_1000GP ~ ExAC.1000GP$freq_exac))$r.squared),digits = 2)
exacPlot <- ggplot(data=ExAC.1000GP, aes(x=freq_exac, y=freq_1000GP)) + geom_point(alpha=0.4, size=2.5) +
  theme_bw() + xlab("ExAC allele frequency") + ylab("1000GP allele frequency") +
  geom_smooth(method='lm') + 
  annotate("text", size=5, x = 0.1, y = 0.95, label = paste("R2 =",exacR2))

# Load HGDP 1000GP comparisons
AFR.hgdp.1000GP <- read.table("../pgx_analyses/freq_comp/africa.hgdp.1000GP.frq", header=T)
AMR.hgdp.1000GP <- read.table("../pgx_analyses/freq_comp/america.hgdp.1000GP.frq", header=T)
EAS.hgdp.1000GP <- read.table("../pgx_analyses/freq_comp/east-asia.hgdp.1000GP.frq", header=T)
EUR.hgdp.1000GP <- read.table("../pgx_analyses/freq_comp/europe.hgdp.1000GP.frq", header=T)
SAS.hgdp.1000GP <- read.table("../pgx_analyses/freq_comp/south-asia.hgdp.1000GP.frq", header=T)

# Plot and calculate R2
afrR2 <- round((summary(lm(AFR.hgdp.1000GP$freq_1000GP ~ AFR.hgdp.1000GP$freq_hgdp))$r.squared),digits = 2)
afrPlot <-ggplot(data=AFR.hgdp.1000GP, aes(x=freq_hgdp, y=freq_1000GP)) + geom_point() +
  theme_bw() + xlab("HGDP African frequency") + ylab("1000GP African frequency") +
  geom_smooth(method='lm') + 
  annotate("text", size=5, x = 0.1, y = 0.95, label = paste("R2 =",afrR2))

amrR2 <- round((summary(lm(AMR.hgdp.1000GP$freq_1000GP ~ AMR.hgdp.1000GP$freq_hgdp))$r.squared),digits = 2)
amrPlot <-ggplot(data=AMR.hgdp.1000GP, aes(x=freq_hgdp, y=freq_1000GP)) + geom_point() +
  theme_bw() + xlab("HGDP Americas frequency") + ylab("1000GP Admixed American frequency") +
  geom_smooth(method='lm') + 
  annotate("text", size=5, x = 0.1, y = 0.95, label = paste("R2 =",amrR2))

easR2 <- round((summary(lm(EAS.hgdp.1000GP$freq_1000GP ~ EAS.hgdp.1000GP$freq_hgdp))$r.squared),digits = 2)
easPlot <-ggplot(data=EAS.hgdp.1000GP, aes(x=freq_hgdp, y=freq_1000GP)) + geom_point() +
  theme_bw() + xlab("HGDP East Asian frequency") + ylab("1000GP East Asian frequency") +
  geom_smooth(method='lm') + 
  annotate("text", size=5, x = 0.1, y = 0.95, label = paste("R2 =",easR2))

eurR2 <- round((summary(lm(EUR.hgdp.1000GP$freq_1000GP ~ EUR.hgdp.1000GP$freq_hgdp))$r.squared),digits = 2)
eurPlot <-ggplot(data=EUR.hgdp.1000GP, aes(x=freq_hgdp, y=freq_1000GP)) + geom_point() +
  theme_bw() + xlab("HGDP European frequency") + ylab("1000GP European frequency") +
  geom_smooth(method='lm') + 
  annotate("text", size=5, x = 0.1, y = 0.95, label = paste("R2 =",eurR2))

sasR2 <- round((summary(lm(SAS.hgdp.1000GP$freq_1000GP ~ SAS.hgdp.1000GP$freq_hgdp))$r.squared),digits = 2)
sasPlot <-ggplot(data=SAS.hgdp.1000GP, aes(x=freq_hgdp, y=freq_1000GP)) + geom_point() +
  theme_bw() + xlab("HGDP South Asian frequency") + ylab("1000GP South Asian frequency") +
  geom_smooth(method='lm') + 
  annotate("text", size=5, x = 0.1, y = 0.95, label = paste("R2 =",sasR2))

plot_grid(exacPlot, afrPlot, amrPlot, easPlot, eurPlot, sasPlot, labels = c("A", "B","C","D","E","F"), align = "h", nrow=2)
ggsave("../pgx_analyses/plots/pgxExternalAlleleFreqComparison.tiff", 
       height=5.36*2, width=7.58*2, units='in', dpi=120)

# Evaluate the quality of 1000GP variant calls
site_eval <- read.table("../pgx_analyses/site_eval/pgx1000GP.svm.analyses.txt",header = T)
marginal_variant_genes <- site_eval[site_eval$marginal_variants>0,]
# Generate proportion of high qc variants for ordering factor for plotting
site_eval$prop_pass <- site_eval$high_qc_variants/site_eval$all_variants
site_eval$gene <- factor(site_eval$gene, levels = site_eval$gene[order(site_eval$prop_pass, decreasing = T)])
site_eval <- site_eval[,-c(2,6)]
# Make tall for plotting and plot
msite_eval <- ddply(melt(site_eval, id.vars = 'gene'), .(gene), mutate, prop = value / sum(value))
ggplot(msite_eval, aes(x = gene, y = prop,fill=variable)) +
  geom_bar(stat='identity') + theme_bw() + 
  theme(legend.title=element_blank(), axis.text.x=element_text(angle = 90, hjust = 1, 
                                                               vjust = 0.5, size=8, face="bold.italic"))
ggsave("../pgx_analyses/plots/pgx_site_eval.tiff", height=5.36, width=2*7.58, units='in', dpi=120) 

# Calculate proportion of different variants for adding to PGx summary
site_eval$all_variants <- site_eval$high_qc_variants +  site_eval$marginal_variants + site_eval$fail_variants
site_eval$percent_high_qc_variants <- site_eval$high_qc_variants/site_eval$all_variants
site_eval$percent_marginal_variants <- site_eval$marginal_variants/site_eval$all_variants
site_eval$percent_fail_variants <- site_eval$fail_variants/site_eval$all_variants
site_eval <- site_eval[,c(1,6,7,8)]
pgxGeneSummary <- merge(pgxGeneSummary, site_eval)

# Look at the flagged inaccessible pharmacogenes
outlierInaccess <- IQR(pgxGeneSummary$percentageInaccessOrSegmental)*3 + quantile(pgxGeneSummary$percentageInaccessOrSegmental)[4]
inaccessPGxgenes <- pgxGeneSummary$gene[pgxGeneSummary$percentageInaccessOrSegmental>outlierInaccess]
pgxGeneSummary$inaccessPGxgene  <- with(datasetAnnotated, ifelse(pgxGeneSummary$gene %in% inaccessPGxgenes, "Yes","No"))

summary(pgxGeneSummary$percent_marginal_variants~pgxGeneSummary$inaccessPGxgene)
summary(pgxGeneSummary$percent_fail_variants~pgxGeneSummary$inaccessPGxgene)

# Plot the site evaluation of the flagged inaccessible pharmacogenes
msite_eval_inaccess <- msite_eval[msite_eval$gene %in% inaccessPGxgenes,]
# Tidy key labels
msite_eval_inaccess$variable <- gsub("high_qc_variants", "High_QC_variants", msite_eval_inaccess$variable)
msite_eval_inaccess$variable <- gsub("_", " ", msite_eval_inaccess$variable)
msite_eval_inaccess$variable <- Hmisc::capitalize(msite_eval_inaccess$variable)

ggplot(msite_eval_inaccess, aes(x = gene, y = prop,fill=variable)) +
  geom_bar(stat='identity') + theme_bw() + xlab("Gene") + ylab("Proportion") +
  theme(legend.title=element_blank(), axis.text.x=element_text(angle = 90, hjust = 1, 
                                                               vjust = 0.5, size=12, face="bold.italic"))
ggsave("../pgx_analyses/plots/pgx_site_eval_inaccess.tiff", height=5.36, width=7.58, units='in', dpi=120) 

# Load calculated coverage data
coveragePgx <- read.table("../pgx_analyses/coverage/pgx1000gp.coverage.txt", header=T, stringsAsFactors=FALSE)
# Melt data for plotting in ggplot
mCoveragePgx <- melt(coveragePgx, id = "gene")
mCoveragePgx$value <- as.numeric(as.character(mCoveragePgx$value))
mCoveragePgx$gene <- as.factor(mCoveragePgx$gene)

# Look at populations separately
samples <- read.table("../pgx_analyses/freq/integrated_call_samples_v3.20130502.ALL.panel",header = T)
mCoveragePgxPops <- merge(mCoveragePgx, samples, by.x="variable", by.y="sample")

# Calculate mean coverage for all samples and all genes
pgxGenesMean <- mean(mCoveragePgxPops$value, na.rm = T)
pgxGenesMean
# Calculate mean coverage per gene
perGeneMeanCoverage <- ddply(mCoveragePgxPops,~gene,summarise,mean_coverage=mean(value, na.rm = T),sd_coverage=sd(value, na.rm = T), median_coverage=median(value, na.rm = T))
perGeneMeanCoverage <- perGeneMeanCoverage[order(perGeneMeanCoverage$mean),]
# Add whether the gene's mean coverage is either >50% or <50% less than global mean
perGeneMeanCoverage$meancoverage50percentLess <- with(perGeneMeanCoverage, ifelse(mean_coverage<pgxGenesMean*0.5, "Yes","No"))
perGeneMeanCoverage$meancoverage50percentMore <- with(perGeneMeanCoverage, ifelse(mean_coverage>pgxGenesMean*1.5, "Yes","No"))
perGeneMeanCoverage$coverageExtremes <- with(perGeneMeanCoverage, ifelse(meancoverage50percentLess==meancoverage50percentMore, "No","Yes"))

# Add to PGx gene summary table
pgxGeneSummary <- merge(pgxGeneSummary,perGeneMeanCoverage)
# Calculate mean coverage per superpopulation
ddply(mCoveragePgxPops,~super_pop,summarise,mean=mean(value, na.rm = T),sd=sd(value, na.rm = T), median=median(value, na.rm = T))

# Plot coverage per gene
# Reorder based on mean coverage
mCoveragePgxPops$gene <- factor(unique(mCoveragePgxPops$gene), 
                                levels=unique(mCoveragePgxPops$gene)[perGeneMeanCoverage$gene])

perGeneCov <- ggplot(data=mCoveragePgxPops, aes(x=gene, y=value)) + geom_boxplot() +
  theme_bw() + xlab("Gene") + ylab("Mean coverage") +
  theme(legend.title=element_blank(), axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5, size=8, face="bold.italic"))
perGeneCov
ggsave("../pgx_analyses/plots/pgxGeneMeanCoverage.tiff", height=5.36, width=7.58*2, units='in', dpi=120)

# Look at coverage by superpopulation for genes with clinical evidence of 1
mCoveragePgxPopsClinical <- mCoveragePgxPops[mCoveragePgxPops$gene %in% rsIDGene$Gene,]

# Plot mean coverage of clinical genes per superpopulation
ggplot(data=mCoveragePgxPopsClinical, aes(x=gene, y=value, fill=super_pop)) + 
  geom_boxplot() +
  theme_bw() + xlab("Gene") + ylab("Mean coverage") +
  scale_fill_manual(values=c("#009E73","#CC79A7", "#56B4E9","#E69F00","grey72")) +
  theme(legend.title=element_blank(), 
        axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5, size=12, face="bold.italic"))
ggsave("../pgx_analyses/plots/pgxGeneMeanClinicalPop.tiff", 
       height=5.36, width=7.58*2, units='in', dpi=120)

# Tidy PGx gene summary for saving
pgxGeneSummary <- pgxGeneSummary[order(pgxGeneSummary$gene),]
write.csv(pgxGeneSummary, "../pgx_analyses/gene_ens_details/pgxGeneSummary.csv", row.names = FALSE)

pgxGeneSummaryVariation <- pgxGeneSummary[,c(1,2,9,18,4,10:17,28,19,25,24,26)]
pgxGeneSummaryVariation[,c(14,15,17)] <- round(pgxGeneSummaryVariation[,c(14,15,17)], digits = 3)
write.csv(pgxGeneSummaryVariation, "../pgx_analyses/gene_ens_details/pgxGeneSummaryVariation.csv", row.names = FALSE)

pgxGeneSummaryQC <- pgxGeneSummary[,c(1,2,21:23,29:37)]
pgxGeneSummaryQC[,3:11] <- round(pgxGeneSummaryQC[,3:11], digits = 3)
write.csv(pgxGeneSummaryQC, "../pgx_analyses/gene_ens_details/pgxGeneSummaryQC.csv", row.names = FALSE)

