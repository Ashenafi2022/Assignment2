


#..............data import from github repository...........
library(dplyr)
fang_genotypes <- read.table("C:/Data analysis/Bioinfo_spring2021/fang_et_al_genotypes.txt", header = FALSE)

  # fang_et_al>genotypes.txt
fang_genotypesXXXX <- ("C:/Data analysis/Bioinfo_spring2021/fang_et_al_genotypes.txt")
fangGeno <- read.table(fang_genotypesXXXX)
dim(fangGeno) # 2782 rows/observations and 986 columns/variables
glimpse(fangGeno)

  # snp_position.txt
snpXXXX <- ("C:/Data analysis/Bioinfo_spring2021/snp_position.txt")
snpPos <- read.table(snpXXXX)
dim(snpPos)
glimpse(fangGeno)

fang_genotypes111 <- ("https://github.com/EEOB-BioData/BCB546-Spring2021/blob/main/assignments/UNIX_Assignment/fang_et_al_genotypes.txt")
fang_genotypes222 <- read_table(fang_genotypes111)
dim(fang_genotypes222) # 851 X 1

data111 <- ("https://github.com/EEOB-BioData/BCB546-Spring2021.git/fang_et_al_genotypes.txt")

data111 <- ("https://github.com/Ashenafi2022/BCB546-Spring2021.git/fang_et_al_genotypes.txt")
data222 <- read.table(data111)
dim(data222)

snp_position1 <- "https://github.com/EEOB-BioData/BCB546-Spring2021/blob/main/assignments/UNIX_Assignment/snp_position.txt"
snp_position12 <- read.table(data111)
dim(snp_position12)

PATH <- "https://raw.githubusercontent.com/guru99-edu/R-Programming/master/travel_times.csv"
df <- read.csv(PATH)
glimpse(df)






#Something is wrong with importing these 2 datasets from github
fang_genotypes1 <- read_csv("https://github.com/EEOB-BioData/BCB546-Spring2021/blob/main/assignments/UNIX_Assignment/fang_et_al_genotypes.txt")
dim(fang_genotypes1) # 851 X 1

snp_position1 <- read_csv("https://github.com/EEOB-BioData/BCB546-Spring2021/blob/main/assignments/UNIX_Assignment/snp_position.txt")
dim(snp_position1) # 4809 X 1




#..............data import...........

# Installing
install.packages("readr")
# Loading
library("readr")
fang_genotypes <- read_tsv(file.choose())
str(fang_genotypes)
dim(fang_genotypes) # 2782 rows/observations and 986 columns/variables

snp_position <- read_tsv(file.choose())
class(snp_position)
str(snp_position) 
dim(snp_position) # 983 rows/observations and  15 columns/variables
head(snp_position) 


# Select SNP_ID, chromosome and Position from SNP_position data set
snp_position.selected <- select(snp_position, c("SNP_ID", "Chromosome", "Position"))
str(snp_position.selected) # 983 rows and 3 columns
glimpse(snp_position.selected) # similar function with str


#...............data processing.........................

# Number of rows and columns
nrow(fang_genotypes) # 2782
ncol(fang_genotypes) # 986
NROW(na.omit(fang_genotypes)) # 2782
NCOL(na.omit(fang_genotypes)) # 986


library(dplyr)
summary(select(fang_genotypes, Group))
levels(as.factor(fang_genotypes$Group)) # 15 groups: "TRIPS" "ZDIPL" "ZLUXR" "ZMHUE" "ZMMIL" "ZMMLR" 
                                        # "ZMMMR" "ZMPBA" "ZMPIL" "ZMPJA" "ZMXCH" "ZMXCP" "ZMXIL" 
                                        # "ZMXNO" "ZMXNT" "ZPERR"
summary(as.factor(fang_genotypes$Group))
# TRIPS ZDIPL ZLUXR ZMHUE ZMMIL ZMMLR ZMMMR ZMPBA ZMPIL ZMPJA ZMXCH ZMXCP ZMXIL ZMXNO ZMXNT ZPERR 
# 22    15    17    10   290  1256    27   900    41    34    75    69     6     7     4     9 
table(as.factor(fang_genotypes$Group))

# Number of chromosomes in snp_position data
summary(as.factor(snp_position$Chromosome))
  # chrosome 1 = 155, 2 =127, 3=107, 4=91, 5=122,  6=76, 7=97, 8=62, 9=60, 10=53, multiple=6, unknown=27 

# change the variable "Group" to a factor
fang_genotypes.selected$Group <- as.factor(fang_genotypes.selected$Group)
levels(fang_genotypes.selected$Group)

#Subset data 3 maize genotypes > transpose > merge with snp data
maize_genotypes <- filter(fang_genotypes, Group == 'ZMMIL' | Group == 'ZMMLR' | Group == 'ZMMMR')
dim(maize_genotypes) # 1573 rows/observations  986 columns/variables
head(maize_genotypes, 10) # top ten rows
levels(as.factor(maize_genotypes$Group)) # "ZMMIL" "ZMMLR" "ZMMMR"
summary(as.factor(maize_genotypes$Group)) # ZMMIL = 290, ZMMLR = 1256, ZMMMR = 27 

      #Transpose
library(tidyverse)
# reading: https://tibble.tidyverse.org/reference/rownames.html
maize_genotypes <- column_to_rownames(maize_genotypes, var = "Sample_ID")
maize_genotypes.tr <- t(maize_genotypes)%>%as.data.frame()%>%rownames_to_column(., var = "SNP_ID")
maize_genotypes.tr <- maize_genotypes.tr[3:nrow(maize_genotypes.tr),]

maize_snp <- merge(snp_position.selected, maize_genotypes.tr, by="SNP_ID")

maize_snp <- select(maize_snp, SNP_ID, Chromosome, Position, everything())

dim(maize_snp) # 983 528

summary(as.factor(maize_snp$Chromosome))

#missed sum
sum(maize_snp$Chromosome == "") # 0


# Subset data 3 teosinte genotypes > transpose > merge with snp data

teosinte_genotypes <- filter(fang_genotypes, Group == 'ZMPBA' | Group == 'ZMPIL' | Group == 'ZMPJA')
dim(teosinte_genotypes) # 975 rows/observations and 986 columns/variables
head(teosinte_genotypes)
levels(as.factor(teosinte_genotypes$Group)) # "ZMPBA" "ZMPIL" "ZMPJA"
summary(as.factor(teosinte_genotypes$Group)) # ZMPBA = 900, ZMPIL = 41, ZMPJA = 34

    #Transpose
teosinte_genotypes <- column_to_rownames(teosinte_genotypes, var = "Sample_ID")
teosinte_genotypes.tr <- t(teosinte_genotypes)%>%as.data.frame()%>%rownames_to_column(., var = "SNP_ID")
teosinte_genotypes.tr <- teosinte_genotypes.tr[3:nrow(teosinte_genotypes.tr),]

teosinte_snp <- merge(snp_position.selected, teosinte_genotypes.tr, by="SNP_ID")

teosinte_snp <- select(teosinte_snp, SNP_ID, Chromosome, Position, everything())

dim(teosinte_snp) # 983 528

# missed sum
sum(teosinte_snp$Chromosome == "") # 0


# Subset data by chromosome, order in increasing position
maize_chromosome1 <- subset(maize_snp, Chromosome==1)%>%arrange(Position)
teosinte_chromosome1 <- subset(teosinte_snp, Chromosome==1)%>%arrange(Position)

teosinte_chromosome2 <- subset(teosinte_snp, Chromosome==2)%>%arrange(teosinte_chromosome2$Position)

#Subset each chromosome, order in position decreasing order, replace missed values
maize_chromosome.dec1 <- subset(maize_snp, Chromosome==1)%>%arrange(desc(Position))
teosinte_chromosome.dec1 <- subset(teosinte_snp, Chromosome==1)%>%arrange(desc(Position))

# missing values
teosinte_chromosome.dec1$Position[(is.na(teosinte_chromosome.dec1$Position))] = "?"

df$y[is.na(df$y)]<-mean(df$y,na.rm=TRUE)

teosinte_chromosome.dec1%>%replace_na(teosinte_chromosome.dec1$Position =="?")


df %>% replace_na(list(x = 0, y = "unknown"))



#sapply

sapply(teosinte_snp, Position, min)
sapply(teosinte_snp$Position, min)






