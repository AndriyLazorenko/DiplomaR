#Visualizations and statistical tests from data received from Java code
#Just copy the relevant output file (e.g. R_Results_Database.csv) 
#into appropriate folder of resources and enjoy!

#Determining a database of results

Database <- readline("Which database? OMIM or SNP? :")
if (Database %in% c("OMIM","Omim","omim", "o", "O")){
    Database <- "OMIM"
    #Determining r or y allele to use in script
    Allele <- readline("Which allele? R or Y? :")
    if (Allele %in% c("r","R")){ 
        Allele <- "R"
        Res <- read.csv("resources/OMIM/R_Results_Database.csv", sep = ";", stringsAsFactors = F) #reading csv
    }
    if (Allele %in% c("y","Y")){
        Allele <- "Y"
        Res <- read.csv("resources/OMIM/Y_Results_Database.csv", sep = ";", stringsAsFactors = F) #reading csv
    }
} else {
    Database <- "SNP"
    #Determining r or y allele to use in script
    Allele <- readline("Which allele? R or Y? :")
    if (Allele %in% c("r","R")){ 
        Allele <- "R"
        Res <- read.csv("resources/SNP/R_Results_Database.csv", sep = ";", stringsAsFactors = F) #reading csv
    }
    if (Allele %in% c("y","Y")){
        Allele <- "Y"
        Res <- read.csv("resources/SNP/Y_Results_Database.csv", sep = ";", stringsAsFactors = F) #reading csv
    }
}

#Reading and pre-processing csv for R

Res[1:16,4] = as.integer(as.double(Res[1:16,2]) * as.double(Res[17,4])) #computing ME after SNP
Res[1:16,5] = as.integer(as.double(Res[1:16,2]) * as.double(Res[17,5])) #computing ME before SNP
Res[1:17,7] = as.double(Res[1:17,6]) / as.double(Res[17,6]) #computing freq after SNP
Res[1:17,9] = as.double(Res[1:17,8]) / as.double(Res[17,8]) #computing freq before SNP
Res[18,7] = as.double(Res[18,6]) / as.double(Res[19,2]) #computing freq after SNP, double SNP
Res[18,9] = as.double(Res[18,8]) / as.double(Res[19,2]) #computing freq before SNP, double SNP
Res[1:16,2] = as.double(Res[1:16,2]) #converting freqs to double

#Creating a tidy df 

strip <- Res[1:16,1:9] #only data included

#Only useful data transfered to new df with cool col names

Dinucleotides <- strip[,1]
Expected <- strip[,2]
#Abs <- strip[,3]
#MEAfter <- strip[,4]
#MEBefore <- strip[,5]
#NAfter <- strip[,6]
After <- strip[,7]
#NBefore <- strip[,8]
Before <- strip[,9]
tidy <- data.frame(Dinucleotides,Before,Expected,After)

#Tidying data a bit

require(tidyr)
tidy <- gather(tidy, FrequencyClass, Frequency, -Dinucleotides)
tidy$Frequency <- as.double(tidy$Frequency)

#Computing subsets - datasets for relevant plots

exp <- tidy[tidy$FrequencyClass=="Expected",-2]
naive <- exp
naive[,2] <- 1/16
before <- tidy[tidy$FrequencyClass=="Before",]
after <- tidy[tidy$FrequencyClass=="After",]
before[,3] <- (before[,3]-exp[,2])/exp[,2]*100
after[,3] <- (after[,3]-exp[,2])/exp[,2]*100
compare <- merge(before,after,all = T)
colnames(compare)[which(names(compare) == "Frequency")] <- "Percent_Difference_From_Expected"

#Serializing data for reuse by SNPvsOMIM.R

path <- paste("resources/",Database,"/",Allele, sep = "")
patha <- paste(path,"after",".rds", sep = "")
pathb <- paste(path,"before",".rds", sep = "")
saveRDS(after, file=patha)
saveRDS(before, file=pathb)

#Calculating statistical significance

atotal <- as.numeric(Res[17,4])
btotal <- as.numeric(Res[17,5])
total[1:16] <- btotal
total[17:32] <- atotal
forProp <-tidy
forProp <- forProp[forProp$FrequencyClass!="Expected",]
forProp$Occurrence <- forProp$Frequency * total
forProp$Total <- total
forProp$Expected_Frequency <- exp$Frequency
forProp$Actual_Frequency <- forProp$Frequency
forPropA <- forProp[forProp$FrequencyClass=="After",-c(2,3,7)]
forPropB <- forProp[forProp$FrequencyClass=="Before",-c(2,3,7)]

f <- function(x) {
    prop.test(x = as.numeric(x[2]),n = as.numeric(x[3]),p = as.numeric(x[4]))
}

#alternative way to work with prop test
#propsA <- Map(prop.test,x=as.numeric(forPropA$Occurrence),
#              n=as.numeric(forPropA$Total), p= as.numeric(forPropA$Expected_Frequency))
#propsB <- Map(prop.test,x=as.numeric(forPropB$Occurrence),
#              n=as.numeric(forPropB$Total), p= as.numeric(forPropB$Expected_Frequency))

propsA <- apply(forPropA, 1,f)
propsB <- apply(forPropB, 1,f)

#Visualizing statistical significance

forProp <- forProp[forProp$FrequencyClass!="Expected",-3]
#getting pval list for after vals
fapval<- function(x){
    apval <- x[[3]][[1]]
}
apval <-sapply(propsA,fapval)
apval <-as.numeric(apval)
#getting pval list for before vals
fbpval<- function(x){
    bpval <- x[[3]][[1]]
}
bpval <-sapply(propsB,fbpval)
bpval <-as.numeric(bpval)
#getting final list of p vals for all

cpval <- c(bpval,apval)

#getting lower confidence interval for after vals
faconflow <- function(x){
    aconflow <- x[[6]][[1]]
}
aconflow <-sapply(propsA,faconflow)
aconflow <-as.numeric(aconflow)
#getting lower confidence interval for before vals
fbconflow <- function(x){
    bconflow <- x[[6]][[1]]
}
bconflow <-sapply(propsB,fbconflow)
bconflow <-as.numeric(bconflow)
#getting final list of lower confidence intervals for all

cconflow <- c(bconflow,aconflow)

#getting upper confidence interval for after vals
faconfup <- function(x){
    aconfup <- x[[6]][[2]]
}
aconfup <-sapply(propsA,faconfup)
aconfup <-as.numeric(aconfup)
#getting upper confidence interval for before vals
fbconfup <- function(x){
    bconfup <- x[[6]][[2]]
}
bconfup <-sapply(propsB,fbconfup)
bconfup <-as.numeric(bconfup)
#getting final list of lower confidence intervals for all

cconfup <- c(bconfup,aconfup)

#getting it all together into one representative dataframe

forProp$Pvalues <- cpval
forProp$Lower_Confidence_Interval <-cconflow
forProp$Upper_Confidence_Interval <-cconfup
forProp$Reject_Random_Frequency <- !(forProp$Pvalues > 0.025 & forProp$Pvalues < 0.975)
require(gridExtra)

forProp[,-c(1,2,3,4,10)] <- signif(forProp[,-c(1,2,3,4,10)], digits = 4)
#optionally, script is able to produce tables as pdf, not as png.
#path <- paste("resources/",Database,"/",Allele,"stat.pdf",sep = "")
#pdf(path, height=11, width=18.5)
path <- paste("resources/",Database,"/img/",Allele,"/stat.png",sep = "")
png(path,height=768, width=1366)
grid.table(forProp)
dev.off()

#Visualizing data

require(ggplot2)
path <- paste("resources/",Database,"/img/",Allele,sep = "")
#Plotting naive frequencies
img_naive_freq <-
    ggplot(data=naive, aes(x=Dinucleotides, y=Frequency, fill=Dinucleotides)) +
    geom_bar(stat="identity")+
    guides(fill=FALSE)
ggsave("naive_freq.png", path = path)
#Plotting mean genomic frequencies
img_exp_freq <- 
    ggplot(data=exp, aes(x=Dinucleotides, y=Frequency, fill=Dinucleotides)) +
    geom_bar(stat="identity")+
    guides(fill=FALSE)
ggsave("exp_freq.png", path = path)
#Plotting comparison of before and after frequencies
img_comp <- ggplot(compare, aes(x=Dinucleotides, y=Percent_Difference_From_Expected, fill=FrequencyClass)) +
    geom_bar(position=position_dodge(), stat="identity")
ggsave("comp.png", path = path)
#Weird and hard-to-handle plot
img <- 
    ggplot(tidy, aes(x=Dinucleotides, y=Frequency, fill=FrequencyClass)) +
    geom_bar(position=position_dodge(), stat="identity")
ggsave("mess.png", path = path)

#Using seqLogo

require(seqLogo)
#Plotting expected frequencies
#setting up data
P1 <- c(sum(exp[1:4,2]),sum(exp[5:8,2]),sum(exp[9:12,2]),sum(exp[13:16,2]))
P2 <- c(sum(exp[c(1,5,9,13),2]),sum(exp[c(2,6,10,14),2]),
        sum(exp[c(3,7,11,15),2]),sum(exp[c(4,8,12,16),2]))
P3 <- c(0.25,0.25,0.25,0.25)
P4 <-c(sum(exp[1:4,2]),sum(exp[5:8,2]),sum(exp[9:12,2]),sum(exp[13:16,2]))
P5 <-c(sum(exp[c(1,5,9,13),2]),sum(exp[c(2,6,10,14),2]),
       sum(exp[c(3,7,11,15),2]),sum(exp[c(4,8,12,16),2]))
forPWM <- data.frame(P1,P2,P3,P4,P5)
#transforming data to relevant format
proportion <- function(x){
    rs <- sum(x);
    return(x / rs);
}
pwm <- apply(forPWM, 2, proportion)
#plotting and saving data
pwm <- makePWM(pwm)
png(paste(path,"/logoExp.png",sep = ""))
logo <- seqLogo(pwm, ic.scale = F)
dev.off()

#Plotting actual frequencies

#setting up data
before <- tidy[tidy$FrequencyClass=="Before",]
after <- tidy[tidy$FrequencyClass=="After",]
#creating df for pwm function
P1 <- c(sum(before[1:4,3]),sum(before[5:8,3]),sum(before[9:12,3]),sum(before[13:16,3]))
P2 <- c(sum(before[c(1,5,9,13),3]),sum(before[c(2,6,10,14),3]),
        sum(before[c(3,7,11,15),3]),sum(before[c(4,8,12,16),3]))
P3 <- c(0.25,0.25,0.25,0.25)
P4 <-c(sum(after[1:4,3]),sum(after[5:8,3]),sum(after[9:12,3]),sum(after[13:16,3]))
P5 <-c(sum(after[c(1,5,9,13),3]),sum(after[c(2,6,10,14),3]),
       sum(after[c(3,7,11,15),3]),sum(after[c(4,8,12,16),3]))
forPWM <- data.frame(P1,P2,P3,P4,P5)
#transforming data to relevant format
proportion <- function(x){
    rs <- sum(x);
    return(x / rs);
}
pwm <- apply(forPWM, 2, proportion)
#plotting and saving data
pwm <- makePWM(pwm)
png(paste(path,"/logoAct.png",sep = ""))
logo <- seqLogo(pwm, ic.scale = F)
dev.off()

#clear all data. Disable in case of debug
#rm(list=ls())
