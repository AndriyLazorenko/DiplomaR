#This script visualizes cross-comparison of freqs of dinucleotides from OMIM and SNP
#It requires output of Visualize.R script, so run Visualize.R first on all data - 
#then this script, otherwise Rafter.RData files won't be available and no
#cross-comparison will occur

#Determining r or y allele to use in script

Allele <- readline("Which allele? R or Y? :")
if (Allele %in% c("r","R")){ 
    Allele <- "R"
}
if (Allele %in% c("y","Y")){
    Allele <- "Y"
}

#Loading data

path <- paste("resources/",sep = "")
#load OMIM data
pathOMIM <- paste(path,"OMIM/",Allele,sep = "")
apathOMIM <- paste(pathOMIM,"after",".rds",sep="")
bpathOMIM <- paste(pathOMIM,"before",".rds",sep="")
afterOMIM <- readRDS(apathOMIM)
beforeOMIM <- readRDS(bpathOMIM)
#load SNP data
pathSNP <- paste(path,"SNP/",Allele,sep = "")
apathSNP <- paste(pathSNP,"after",".rds",sep="")
bpathSNP <- paste(pathSNP,"before",".rds",sep="")
afterSNP <- readRDS(apathSNP)
beforeSNP <- readRDS(bpathSNP)

#Visualizing data
#cooking datasets
#after
afterSNP$FrequencyClass <- paste(afterSNP$FrequencyClass,"SNP",sep = "")
afterOMIM$FrequencyClass <- paste(afterOMIM$FrequencyClass,"OMIM",sep = "")
acompare <- merge(afterSNP,afterOMIM,all = T)
colnames(acompare)[which(names(acompare) == "Frequency")] <- "Percent_Difference_From_Expected"
#before
beforeSNP$FrequencyClass <- paste(beforeSNP$FrequencyClass,"SNP",sep = "")
beforeOMIM$FrequencyClass <- paste(beforeOMIM$FrequencyClass,"OMIM",sep = "")
bcompare <- merge(beforeSNP,beforeOMIM,all = T)
colnames(bcompare)[which(names(bcompare) == "Frequency")] <- "Percent_Difference_From_Expected"

#Plotting comparison of before and after frequencies
path <- paste("resources/comparison/",Allele,sep = "")
#plot after
img_comp <- ggplot(acompare, aes(x=Dinucleotides, y=Percent_Difference_From_Expected, fill=FrequencyClass)) +
    geom_bar(position=position_dodge(), stat="identity")
ggsave("acomp.png", path = path)
#plot before
img_comp <- ggplot(bcompare, aes(x=Dinucleotides, y=Percent_Difference_From_Expected, fill=FrequencyClass)) +
    geom_bar(position=position_dodge(), stat="identity")
ggsave("bcomp.png", path = path)

#clear all data. Disable in case of debug
rm(list=ls())
