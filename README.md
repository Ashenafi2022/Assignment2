


#..............data import from github repository...........
library(dplyr)
fang_genotypes1 <- ("https://raw.githubusercontent.com/EEOB-BioData/BCB546-Spring2021/master/assignments/UNIX_Assignment/fang_et_al_genotypes.txt")
fang_genotypes <- read.table(fang_genotypes1, head = TRUE, sep="\t")
dim(fang_genotypes) # 2782 observations/rows and 986 variables/columns

snp_position1 <- ("https://raw.githubusercontent.com/EEOB-BioData/BCB546-Spring2021/master/assignments/UNIX_Assignment/snp_position.txt")
snp_position <- read.table(snp_position1, head = TRUE, sep="\t")
dim(snp_position) # 983 observations or rows and 15 variables or columns

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


# Homozygous and heterozygous
# Four catageroies: homozygous, hetrozygous, missed, and NA (e.g. A/C, G/T)

# Gene1 = ZDP_0752a
ZDP_0752a <- transmute(maize_snp, nucleotides =case_when( 
              ZDP_0752a %in% c('A/A', 'C/C', 'G/G', 'T/T')~"homozygous",
              ZDP_0752a %in% c('A/T', 'C/G', 'G/C', 'T/A') ~"hetrozygous",
              ZDP_0752a %in% '?/?' ~ "Missed"))
table(ZDP_0752a)
# hetrozygous  homozygous      Missed 
#     36         793          22 

ggplot(ZDP_0752a, aes(nucleotides)) +
  geom_bar(fill="#0073C2FF") 
  
ggplot(ZDP_0752a)+ geom_point(mapping=aes(ZDP_0752a))
hist(as.numeric(ZDP_0752a$nucleotides)) + labs (title="Gene 1 - ZDP_0752a")

#Gene2 = ZDP_0793a
ZDP_0793a <- transmute(maize_snp, nucleotides =case_when( 
                  ZDP_0793a %in% c('A/A', 'C/C', 'G/G', 'T/T')~"homozygous",
                  ZDP_0793a %in% c('A/T', 'C/G', 'G/C', 'T/A') ~"hetrozygous",
                  ZDP_0752a %in% '?/?' ~ "Missed"))
table(ZDP_0793a)
# hetrozygous  homozygous      Missed 
#   46         723           7 

ggplot(ZDP_0793a, aes(nucleotides)) +
  geom_bar(fill="#0073C2FF") + labs (title="Gene 2 - ZDP_0793a")



