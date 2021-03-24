

# Author: Ashenafi Beyi
# Date: March 19, 2021
# Purpose: EEOB-BCB-546 assignment 2: data processing and vizualizing


# .....STEP 1. Import two data files from the course repository..........
#########################################################################

fang_genotypes <- read.table("https://raw.githubusercontent.com/EEOB-BioData/BCB546-Spring2021/master/assignments/UNIX_Assignment/fang_et_al_genotypes.txt",
                              header=T, sep="\t", stringsAsFactors =F)
dim(fang_genotypes) # 2782 observations/rows and 986 variables/columns

snp_position <- read.table("https://raw.githubusercontent.com/EEOB-BioData/BCB546-Spring2021/master/assignments/UNIX_Assignment/snp_position.txt",
                  header=T, sep="\t", stringsAsFactors =F)
dim(snp_position) # 983 observations or rows and 15 variables or columns
#........................................................................

# ....STEP 2. Data processing............................................
#########################################################################

# Select SNP_ID, chromosome and Position from SNP_position data set
library(dplyr)
snp_position.selected <- select(snp_position, c("SNP_ID", "Chromosome", "Position"))
str(snp_position.selected) # 983 rows and 3 columns
glimpse(snp_position.selected) # similar function with str

# Number of rows and columns in each data file
nrow(fang_genotypes) # 2782
ncol(fang_genotypes) # 986
NROW(na.omit(fang_genotypes)) # 2782
NCOL(na.omit(fang_genotypes)) # 986

# Summarize the frequencies of maize and teaosinte genotypes
  # Three approaches/function are used: levels, summary, and table

levels(as.factor(fang_genotypes$Group)) # 15 groups: "TRIPS" "ZDIPL" "ZLUXR" "ZMHUE" "ZMMIL" "ZMMLR" 
  # "ZMMMR" "ZMPBA" "ZMPIL" "ZMPJA" "ZMXCH" "ZMXCP" "ZMXIL" 
  # "ZMXNO" "ZMXNT" "ZPERR"

summary(as.factor(fang_genotypes$Group))
  # TRIPS ZDIPL ZLUXR ZMHUE ZMMIL ZMMLR ZMMMR ZMPBA ZMPIL ZMPJA ZMXCH ZMXCP ZMXIL ZMXNO ZMXNT ZPERR 
  # 22    15    17    10   290  1256    27   900    41    34    75    69     6     7     4     9 
table(as.factor(fang_genotypes$Group))

# Summarize frequencies of chromosomes in snp_position data file
summary(as.factor(snp_position$Chromosome))
  # chrosome 1 = 155, 2 =127, 3=107, 4=91, 5=122,  6=76, 7=97, 8=62, 9=60, 10=53, multiple=6, unknown=27 

                    # Data subsetting

  #......Maize genotypes..........................................

#Subset data 3 maize genotypes > transpose > merge with snp data
maize_genotypes <- filter(fang_genotypes, Group == 'ZMMIL' | Group == 'ZMMLR' | Group == 'ZMMMR')
dim(maize_genotypes) # 1573 rows/observations  986 columns/variables
  #  head(maize_genotypes, 10) # top ten rows
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

dim(maize_snp) # 983 1576

summary(as.factor(maize_snp$Chromosome)) # to varify the chromosome numbers

#missed sum
sum(maize_snp$Chromosome == "") # 0

   #......Teosinte genotypes....................................
# Subset data 3 teosinte genotypes > transpose > merge with snp data

teosinte_genotypes <- filter(fang_genotypes, Group == 'ZMPBA' | Group == 'ZMPIL' | Group == 'ZMPJA')
dim(teosinte_genotypes) # 975 rows/observations and 986 columns/variables
 # head(teosinte_genotypes)
levels(as.factor(teosinte_genotypes$Group)) # "ZMPBA" "ZMPIL" "ZMPJA"
summary(as.factor(teosinte_genotypes$Group)) # ZMPBA = 900, ZMPIL = 41, ZMPJA = 34

#Transpose
teosinte_genotypes <- column_to_rownames(teosinte_genotypes, var = "Sample_ID")
teosinte_genotypes.tr <- t(teosinte_genotypes)%>%as.data.frame()%>%rownames_to_column(., var = "SNP_ID")
teosinte_genotypes.tr <- teosinte_genotypes.tr[3:nrow(teosinte_genotypes.tr),]

teosinte_snp <- merge(snp_position.selected, teosinte_genotypes.tr, by="SNP_ID")

teosinte_snp <- select(teosinte_snp, SNP_ID, Chromosome, Position, everything())

dim(teosinte_snp) # 983 978
#..............................................................................


# ....STEP 3..........subset joined data set by 'Chromosome'...................
###############################################################################
  #..............order them by 'position'......................................
    # ..........replace missed value by "?" or "-".............................


# Subset data by chromosome, order in increasing by "position" and 
  # replace missed values in columns 4 to 978 by "?"

 # Maize genotypes
maize_chromosome.incr1 <- subset(maize_snp, Chromosome==1)%>%arrange(Position)%>%
  mutate_at(4:978, ~replace(., is.na(.), "?"))
maize_chromosome.incr2 <- subset(maize_snp, Chromosome==2)%>%arrange(Position)%>%
  mutate_at(4:978, ~replace(., is.na(.), "?"))
maize_chromosome.incr3 <- subset(maize_snp, Chromosome==3)%>%arrange(Position)%>%
  mutate_at(4:978, ~replace(., is.na(.), "?"))
maize_chromosome.incr4 <- subset(maize_snp, Chromosome==4)%>%arrange(Position)%>%
  mutate_at(4:978, ~replace(., is.na(.), "?"))
maize_chromosome.incr5 <- subset(maize_snp, Chromosome==5)%>%arrange(Position)%>%
  mutate_at(4:978, ~replace(., is.na(.), "?"))
maize_chromosome.incr6 <- subset(maize_snp, Chromosome==6)%>%arrange(Position)%>%
  mutate_at(4:978, ~replace(., is.na(.), "?"))
maize_chromosome.incr7 <- subset(maize_snp, Chromosome==7)%>%arrange(Position)%>%
  mutate_at(4:978, ~replace(., is.na(.), "?"))
maize_chromosome.incr8 <- subset(maize_snp, Chromosome==8)%>%arrange(Position)%>%
  mutate_at(4:978, ~replace(., is.na(.), "?"))
maize_chromosome.incr9 <- subset(maize_snp, Chromosome==9)%>%arrange(Position)%>%
  mutate_at(4:978, ~replace(., is.na(.), "?"))
maize_chromosome.incr10 <- subset(maize_snp, Chromosome==10)%>%arrange(Position)%>%
  mutate_at(4:978, ~replace(., is.na(.), "?"))

 # Teosinte genotypes
teosinte_chromosome.incr1 <- subset(teosinte_snp, Chromosome==1)%>%arrange(Position) %>% 
  mutate_at(4:978, ~replace(., is.na(.), "?"))
teosinte_chromosome.incr2 <- subset(teosinte_snp, Chromosome==2)%>%arrange(Position) %>% 
  mutate_at(4:978, ~replace(., is.na(.), "?"))
teosinte_chromosome.incr3 <- subset(teosinte_snp, Chromosome==3)%>%arrange(Position) %>% 
  mutate_at(4:978, ~replace(., is.na(.), "?"))
teosinte_chromosome.incr4 <- subset(teosinte_snp, Chromosome==4)%>%arrange(Position) %>% 
  mutate_at(4:978, ~replace(., is.na(.), "?"))
teosinte_chromosome.incr5 <- subset(teosinte_snp, Chromosome==5)%>%arrange(Position) %>% 
  mutate_at(4:978, ~replace(., is.na(.), "?"))
teosinte_chromosome.incr6 <- subset(teosinte_snp, Chromosome==6)%>%arrange(Position) %>% 
  mutate_at(4:978, ~replace(., is.na(.), "?"))
teosinte_chromosome.incr7 <- subset(teosinte_snp, Chromosome==7)%>%arrange(Position) %>% 
  mutate_at(4:978, ~replace(., is.na(.), "?"))
teosinte_chromosome.incr8 <- subset(teosinte_snp, Chromosome==8)%>%arrange(Position) %>% 
  mutate_at(4:978, ~replace(., is.na(.), "?"))
teosinte_chromosome.incr9 <- subset(teosinte_snp, Chromosome==9)%>%arrange(Position) %>% 
  mutate_at(4:978, ~replace(., is.na(.), "?"))
teosinte_chromosome.incr10 <- subset(teosinte_snp, Chromosome==10)%>%arrange(Position) %>% 
  mutate_at(4:978, ~replace(., is.na(.), "?"))


# Subset data by chromosome, order in increasing by "position" and 
  # replace missed values in columns 4 to 978 by "-"

# Maize genotypes
maize_chromosome.dec1 <- subset(maize_snp, Chromosome==1)%>%arrange(desc(Position))%>%
  mutate_at(4:978, ~replace(., is.na(.), "-"))
maize_chromosome.dec2 <- subset(maize_snp, Chromosome==2)%>%arrange(desc(Position))%>%
  mutate_at(4:978, ~replace(., is.na(.), "-"))
maize_chromosome.dec3 <- subset(maize_snp, Chromosome==3)%>%arrange(desc(Position))%>%
  mutate_at(4:978, ~replace(., is.na(.), "-"))
maize_chromosome.dec4 <- subset(maize_snp, Chromosome==4)%>%arrange(desc(Position))%>%
  mutate_at(4:978, ~replace(., is.na(.), "-"))
maize_chromosome.dec5 <- subset(maize_snp, Chromosome==5)%>%arrange(desc(Position))%>%
  mutate_at(4:978, ~replace(., is.na(.), "-"))
maize_chromosome.dec6 <- subset(maize_snp, Chromosome==6)%>%arrange(desc(Position))%>%
  mutate_at(4:978, ~replace(., is.na(.), "-"))
maize_chromosome.dec7 <- subset(maize_snp, Chromosome==7)%>%arrange(desc(Position))%>%
  mutate_at(4:978, ~replace(., is.na(.), "-"))
maize_chromosome.dec8 <- subset(maize_snp, Chromosome==8)%>%arrange(desc(Position))%>%
  mutate_at(4:978, ~replace(., is.na(.), "-"))
maize_chromosome.dec9 <- subset(maize_snp, Chromosome==9)%>%arrange(desc(Position))%>%
  mutate_at(4:978, ~replace(., is.na(.), "-"))
maize_chromosome.dec10 <- subset(maize_snp, Chromosome==10)%>%arrange(desc(Position))%>%
  mutate_at(4:978, ~replace(., is.na(.), "-"))

# Teosinte genotypes
teosinte_chromosome.dec1 <- subset(teosinte_snp, Chromosome==1)%>%arrange(desc(Position))%>%
  mutate_at(4:978, ~replace(., is.na(.), "-"))
teosinte_chromosome.dec2 <- subset(teosinte_snp, Chromosome==2)%>%arrange(desc(Position))%>%
  mutate_at(4:978, ~replace(., is.na(.), "-"))
teosinte_chromosome.dec3 <- subset(teosinte_snp, Chromosome==3)%>%arrange(desc(Position))%>%
  mutate_at(4:978, ~replace(., is.na(.), "-"))
teosinte_chromosome.dec4 <- subset(teosinte_snp, Chromosome==4)%>%arrange(desc(Position))%>%
  mutate_at(4:978, ~replace(., is.na(.), "-"))
teosinte_chromosome.dec5 <- subset(teosinte_snp, Chromosome==5)%>%arrange(desc(Position))%>%
  mutate_at(4:978, ~replace(., is.na(.), "-"))
teosinte_chromosome.dec6 <- subset(teosinte_snp, Chromosome==6)%>%arrange(desc(Position))%>%
  mutate_at(4:978, ~replace(., is.na(.), "-"))
teosinte_chromosome.dec7 <- subset(teosinte_snp, Chromosome==7)%>%arrange(desc(Position))%>%
  mutate_at(4:978, ~replace(., is.na(.), "-"))
teosinte_chromosome.dec8 <- subset(teosinte_snp, Chromosome==8)%>%arrange(desc(Position))%>%
  mutate_at(4:978, ~replace(., is.na(.), "-"))
teosinte_chromosome.dec9 <- subset(teosinte_snp, Chromosome==9)%>%arrange(desc(Position))%>%
  mutate_at(4:978, ~replace(., is.na(.), "-"))
teosinte_chromosome.dec10 <- subset(teosinte_snp, Chromosome==10)%>%arrange(desc(Position))%>%
  mutate_at(4:978, ~replace(., is.na(.), "-"))
#..................................................................................

write.table(teosinte_chromosome.dec10, "https://raw.githubusercontent.com/EEOB-BioData/BCB546-Spring2021/master/assignments/UNIX_Assignment/teosinte_chromosome.dec10.txt", sep="\t")

write.table(maize_chromosome.dec10, "C:/Data analysis/Bioinfo_spring2021/Assignment2/maize_chromosome.dec10.txt", sep="\t")



write.table(teosinte_chromosome.incr10, "C:/Data analysis/Bioinfo_spring2021/Assignment2/teosinte_chromosome.incr10.txt", sep="\t")



# Apply functions

# sapply(as.numeric(teosinte_snp$Position), range)
# sapply(teosinte_snp$Position, max)

# ....STEP 4................Data vizualization.....................................
###################################################################################

# 1. Bar cahrts for SNPs per chromosome in maize and teosinte genotypes

maize_snp_bar <- ggplot(data = maize_snp) + geom_point(mapping=aes(x=Chromosome, y=Position)) +
  ggtitle("SNP distribution in Chromosome - Maize genotypes")
maize_snp_bar + theme(plot.title = element_text(color = "green", size = 12, face = "bold", hjust = 0.5))
teosinte_snp_bar <- ggplot(data = teosinte_snp) + geom_point(mapping=aes(x=Chromosome, y=Position)) +
  ggtitle("SNP distribution in Chromosome - Teosinte genotypes")
teosinte_snp_bar + theme(plot.title = element_text(color = "green", size = 12, face = "bold", hjust = 0.5))

# 2. Bar charts showing Homozygous and heterozygous distribution in selected genes
# Four catageroies: homozygous, hetrozygous, missed, and NA (e.g. A/C, G/T)

# Gene1 = ZDP_0752a
ZDP_0752a <- transmute(maize_snp, nucleotides =case_when( 
  ZDP_0752a %in% c('A/A', 'C/C', 'G/G', 'T/T')~"homozygous",
  ZDP_0752a %in% c('A/T', 'C/G', 'G/C', 'T/A') ~"hetrozygous",
  ZDP_0752a %in% '?/?' ~ "Missed"))
table(ZDP_0752a)
# hetrozygous  homozygous      Missed 
#     36         793          22 

Gene1 <- ggplot(ZDP_0752a, aes(nucleotides)) + geom_bar(fill="#0073C2FF") +
            ggtitle ("Gene 1 - ZDP_0752a")
Gene1 + theme(plot.title = element_text(color = "black", size = 12, face = "bold", hjust = 0.5))


#Gene2 = ZDP_0793a
ZDP_0793a <- transmute(maize_snp, nucleotides =case_when( 
  ZDP_0793a %in% c('A/A', 'C/C', 'G/G', 'T/T')~"homozygous",
  ZDP_0793a %in% c('A/T', 'C/G', 'G/C', 'T/A') ~"hetrozygous",
  ZDP_0752a %in% '?/?' ~ "Missed"))
table(ZDP_0793a)
# hetrozygous  homozygous      Missed 
#   46         723           7 

Gene2 <- ggplot(ZDP_0793a, aes(nucleotides)) + geom_bar(fill="#0073C2FF") + 
           ggtitle("Gene 2 - ZDP_0793a")
Gene2 + theme(plot.title = element_text(color = "black", size = 12, face = "bold", hjust = 0.5))


# 3. bar charts showing chromosome counts and SNP distributions in the entire fang_genotypes file 

fang_genotypes1 <- column_to_rownames(fang_genotypes, var = "Sample_ID")
fang_genotypes.tr <- t(fang_genotypes1)%>%as.data.frame()%>%rownames_to_column(., var = "SNP_ID")
fang_genotypes.tr <- fang_genotypes.tr[3:nrow(maize_genotypes.tr),]

joined_snp_genotypes <- merge(snp_position, fang_genotypes.tr, by= "SNP_ID")

# Bar chart of Chromosome frequencies
chromosome_count <- ggplot(joined_snp_genotypes, aes(x= Chromosome) ) + 
  geom_histogram(stat= "Count", color = "blue", fill = "blue") + ggtitle ("Number of SNPs per chromosome")

chromosome_count + theme(plot.title = element_text(color = "red", size = 12, face = "bold", hjust = 0.5))
# The title is in red color, size 12, centered and bold.

# Bar chart of SNP distribution across chromosomes

chromosome_snp <- ggplot(joined_snp_genotypes, aes(x= Chromosome, y= Position))+ 
  geom_point(stat=, color = "red", alpha= 0.01)+ ggtitle ("SNPs distribution per chromosome")

chromosome_snp + theme(plot.title = element_text(color = "red", size = 12, face = "bold", hjust = 0.5))
# The title is in red color, size 12, centered and bold.
#..........................................................................................

#..............Summary of the work flows and what has been done.............................
############################################################################################

#  1. Two data sets were imported from the course repository

#  2. Data processing was conducted:
      # i. Data size and structure were explored.
      # ii. Three target variables selected in SNP_position data set.
      # iii. Target maize and teosinte genotypes were filtered, and 
         # transposed and merged with SNP_position data set.

#  3. The merged SNP position and genotypes were subset by chromosome, order by "Position", and
      # missed values replaced by "?" or "-" for maize and teosinte genotypes separately.

#  4. Data vizualization
      # i. Chromosome counts and SNP distribution in maize and teosinte genotypes files
      # ii.two genes were selected for data vizualization exercise
         # > the selected gene was catagorized into "homozygous", "hetrozygous" and "missed"
          # > Bar chart showing nucleotide distribution was plotted.
      # iii. the entire genotype data set was merged with snp position data
       #   > bar charts depcting chromosome ferquencies and snp distribution across 
             # the chromsomes were plotted
#...................................................................................

