#Determining r or y allele to use in script

Answer <- readline("Which allele? R or Y?")
if (Answer %in% c("r","R")){ 
    Res <- read.csv("resources/R_Results_Database.csv", sep = ";", stringsAsFactors = F) #reading csv
}
if (Answer %in% c("y","Y")){
    Res <- read.csv("resources/Y_Results_Database.csv", sep = ";", stringsAsFactors = F) #reading csv
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

#only useful data transfered to new df with cool col names

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

#tidying data a bit

require(tidyr)
tidy <- gather(tidy, FrequencyClass, Frequency, -Dinucleotides)
tidy$Frequency <- as.double(tidy$Frequency)

#Visualizing data

require(ggplot2)
img <- ggplot(tidy, aes(x=Dinucleotides, y=Frequency, fill=FrequencyClass)) + geom_bar(position=position_dodge(), stat="identity")
img


