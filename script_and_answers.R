setwd("/Users/dominique/Documents/JOBs/Artios/test")
install.packages("tidyverse")
library(tidyverse)
install.packages("readxl")
library(readxl)
install.packages("ggrepel")
library(ggrepel)

################################
### 1. Read in both datasets ###
################################
# Read in the Olivieri dataset. 
dat.oli <- read.csv("Olivieri2020_drugz_subset.csv")
# Read in the Hustedt dataset and reorganise for easier processing.
dat.hus.0 <- read_xlsx("Hustedt et al. 2019 - results - rsob190156supp2.xlsx",sheet = 1,skip = 1)
dat.hus <- bind_rows((dat.hus.0[-1,c(1,2:4)])%>%`colnames<-`(.,c("GENE","normZ","Pval","FDR"))%>%mutate(cell_line = "RPE1-hTER", treatment = "AZD6738"),
                     (dat.hus.0[-1,c(1,5:7)])%>%`colnames<-`(.,c("GENE","normZ","Pval","FDR"))%>%mutate(cell_line = "RPE1-hTER", treatment = "VE821"),
                     (dat.hus.0[-1,c(1,8:10)])%>%`colnames<-`(.,c("GENE","normZ","Pval","FDR"))%>%mutate(cell_line = "Hela", treatment = "VE821"),
                     (dat.hus.0[-1,c(1,11:13)])%>%`colnames<-`(.,c("GENE","normZ","Pval","FDR"))%>%mutate(cell_line = "HCT116", treatment = "VE821"))
dat.hus$normZ<-as.numeric(dat.hus$normZ)
dat.hus$Pval<-as.numeric(dat.hus$Pval)
dat.hus$FDR<-as.numeric(dat.hus$FDR)
################################
### 2. Explore both datasets ###
################################
# a. Is there anything unusual about the data?
#   There are quite a few larger than 1 FDRs in both datasets. Which is unusual
#   by the definition of FDR and I have not encountered this before.
#   I noticed in the Olivieri dataset, except for Cisplatin1 treatment, 
#   all other treatments have NA entries for some genes. The Hustedt 
#   dataset has no NA entry.
length(dat.oli[is.na(dat.oli)==T])
# [1] 1274
for (i in unique(dat.oli$treatment)) {
  print(i)
  print(paste0("normZ = NA ",length(dat.oli$normZ[is.na(dat.oli$normZ)==T&dat.oli$treatment==i])))
  print(paste0("FDR = NA ",length(dat.oli$FDR[is.na(dat.oli$FDR)==T&dat.oli$treatment==i])))
}
# [1] "Cisplatin1"
# [1] "normZ = NA 0"
# [1] "FDR = NA 0"
# [1] "IR"
# [1] "normZ = NA 67"
# [1] "FDR = NA 67"
# [1] "UV"
# [1] "normZ = NA 21"
# [1] "FDR = NA 21"
# [1] "Olaparib"
# [1] "normZ = NA 105"
# [1] "FDR = NA 105"
# [1] "AZD6738"
# [1] "normZ = NA 110"
# [1] "FDR = NA 110"
# [1] "Cisplatin2"
# [1] "normZ = NA 133"
# [1] "FDR = NA 133"
# [1] "Cisplatin3"
# [1] "normZ = NA 110"
# [1] "FDR = NA 110"
# [1] "Formaldehyde"
# [1] "normZ = NA 89"
# [1] "FDR = NA 89"
# [1] "KBrO3"
# [1] "normZ = NA 2"
# [1] "FDR = NA 2"
# Get rid the genes with NA across treatments in the dataset:
gene.NA <- unique(dat.oli$GENE[is.na(dat.oli$normZ)])
dat.oli <- dat.oli%>%filter(!GENE%in%gene.NA)
length(dat.hus[is.na(dat.hus)==T])
# [1] 0


# b. Are the data comparable?
# The effect sizes (here, synthetic lethality) are 
# comparable as normalisation is done across replicates within
# each independent experiment, the z scores are centralised and 
# scaled by the standard deviations, respectively. However, 
# caution should be practised when comparing datasets from
# different experiments. The batch effect is the most overt
# cause of experimental variations. It can be minimised from
# the experimental design or later be controlled during the
# analyses. 

#########################################
### 3. Answer the following question: ###
#########################################
# a. How many cell lines are there in these two datasets? 
    # There is the RPE1-hTERT cell line in the Olivieri dataset as mentioned in the background material.
    # There are 3 cell line s in the Hustedt dataset: RPE1-hTER, Hela and HCT116.
unique(dat.hus$cell_line)
    # [1] "RPE1-hTER" "Hela"      "HCT116"

# b. How many treatments?
    # Olivieri dataset has 9 treatments.
unique(dat.oli$treatment)
# [1] "Cisplatin1"   "IR"           "UV"           "Olaparib"     "AZD6738"      "Cisplatin2"   "Cisplatin3"  
# [8] "Formaldehyde" "KBrO3" 
    # Hustedt dataset has 2 treatments.
unique(dat.hus$treatment)
# [1] "AZD6738" "VE821"

#   i. What two treatments are the most similar?
    # The most similar treatments are Cisplation2 and Cisplatin3, 
    # according to the correlation between the Z scores of the treatments. 

# Create an empty dataframe to store the calculated comparisons:
comp <- data.frame(comparison = character(0), correlation = numeric(0))
# Combine the 2 datasets for easy access. Get rid of NA data and redundant entries:
dat <- bind_rows(
  x = dat.oli%>%filter(GENE%in%(intersect(unique(dat.hus$GENE),unique(dat.oli$GENE)))),
  y = dat.hus%>%filter(cell_line == "RPE1-hTER",GENE%in%(intersect(unique(dat.hus$GENE),unique(dat.oli$GENE))))%>%select(-cell_line, -Pval),
)%>%distinct(GENE, treatment, normZ, .keep_all = T)
# Calculate the correlation of Z scores between each pair of treatments:
for (i in 1:(length(unique(dat$treatment))-1)) {
  for (j in (i+1):length(unique(dat$treatment))) {
    comp[nrow(comp) + 1, ] <- c(paste0("cor(",unique(dat$treatment)[i],",",unique(dat$treatment)[j],")"),
                                       cor((dat%>%filter(treatment == (unique(treatment)[i])))$normZ,
                                           (dat%>%filter(treatment == (unique(treatment)[j])))$normZ,
                                           use = "pairwise.complete.obs",method = "pearson"))
  }
}
# Extract the most correlated pair:
comp$comparison[comp$correlation == max(comp$correlation)]
# [1] "cor(Cisplatin2,Cisplatin3)"
pdf("Correlation_between_treatments.pdf")
ggplot(data = merge.data.frame(x = dat%>%filter(treatment=="Cisplatin2")%>%select(GENE,normZ),
                               y = dat%>%filter(treatment=="Cisplatin3")%>%select(GENE,normZ),
                               by = "GENE",suffixes = c("Cisplatin2","Cisplatin3")),
       aes(x = normZCisplatin2,
           y = normZCisplatin3))+
  geom_point()+
  geom_abline(slope = 1,intercept = 0,colour = "red",linetype=3)+
  annotate(geom = "text", x = -19,y = 1,label = "r = 0.58",size = 6)+
  xlab("Cisplatin2 normalised Z score") + ylab("Cisplatin3 normalised Z score")+
  theme_bw()
dev.off()


#   c. Assuming a hit is any gene meeting the following thresholds of normZ <= -2.5 and
# FDR < 0.1, how many hits does each treatment generate? 
dat%>%filter(normZ <= -2.5, FDR < 0.1)%>%count(treatment)
# treatment   n
# 1       AZD6738  76
# 2    Cisplatin1  49
# 3    Cisplatin2  96
# 4    Cisplatin3 115
# 5  Formaldehyde  37
# 6            IR  14
# 7         KBrO3 147
# 8      Olaparib  20
# 9            UV 104
# 10        VE821 120

#   d. How would you prioritise genes?
    # I would rank the genes by their normZ scores, prioritise the more negative (more lethal) ones.
#   e. How would you visualise top genes?
    # I would present a volcano plot by plotting the effect sizes (Z scores) on the x axis,
    # agaisnt the -log10(P values) on the y axis for each gene, similar to the visualisation 
    # of differentially expressed genes. This will allow the examination of the effects of 
    # the genes in relation to the significance of their effects.
pdf("Visualisation_of_top_genes_example.pdf")
ggplot(data = dat.hus%>%filter(treatment=="AZD6738",cell_line=="RPE1-hTER")%>%
         mutate(colour = ifelse(test = normZ <= -2.5&FDR < 0.1,yes = "red",no = "black")), 
       aes(x = normZ,y = -log10(Pval)))+
  geom_point(aes(colour = colour))+
  scale_color_identity()+
  geom_text_repel(aes(label = ifelse(test = normZ <= -2.5&FDR < 0.1,yes = GENE,no = "")),
                  # max.overlaps = Inf,
                  colour = "red")+
  labs(title = "Hustedt_RPE1-hTER_AZD6738")+
  theme_bw()
dev.off()

# 4. What other questions could you explore across both datasets?
# 5. ATR inhibitor (ATRi) specific analysis:
#   a. VE821 and AZD6738 are both ATRi, are they more similar to each other compared to other non-ATRi treatments?
      # Compared to other treatments, VE821 and AZD6738 render more similar effects. 
      # The correlation of normalised Z scores between VE821 and AZD6738 are
      # the second strongest among their comparisons to other treatments.
comp%>%filter(grepl("VE821",comparison)==T|grepl("AZD6738",comparison)==T)%>%arrange(desc(correlation))
# comparison        correlation
# 1           cor(KBrO3,VE821)  0.247771913905029
# 2         cor(AZD6738,VE821)  0.218487985366889
# 3         cor(AZD6738,KBrO3)  0.181503798260185
# 4    cor(AZD6738,Cisplatin2)  0.162766588732243
# 5    cor(AZD6738,Cisplatin3)  0.159446453134058
# 6      cor(Cisplatin3,VE821)  0.146048425250444
# 7              cor(UV,VE821)  0.142751956038224
# 8      cor(Cisplatin2,VE821)  0.131584463451455
# 9      cor(Cisplatin1,VE821) 0.0611031709054065
# 10   cor(Cisplatin1,AZD6738) 0.0602614872617928
# 11 cor(AZD6738,Formaldehyde) 0.0587230428486212
# 12     cor(Olaparib,AZD6738) 0.0424888166241379
# 13       cor(Olaparib,VE821)  0.040832452248118
# 14             cor(IR,VE821) 0.0337116988510621
# 15           cor(UV,AZD6738) 0.0284127261120482
# 16           cor(IR,AZD6738) 0.0211348663295456
# 17   cor(Formaldehyde,VE821)  0.013678530018786

#   b. What are the top scoring genes in ATRi?
dat%>%filter(normZ <= -2.5, FDR < 0.1, treatment == "VE821"|treatment == "AZD6738")%>%select(GENE)
pdf("Visualisation_of_top_genes_ATRi.pdf")
ggplot(data = dat%>%filter(treatment=="AZD6738"|treatment=="VE821")%>%
         mutate(colour = ifelse(test = normZ <= -2.5&FDR < 0.1,yes = "red",no = "black")), 
       aes(x = normZ,y = -log10(FDR)))+
  facet_wrap(~treatment)+
  geom_point(aes(colour = colour))+
  scale_color_identity()+
  geom_text_repel(aes(label = ifelse(test = normZ <= -2.5&FDR < 0.1,yes = GENE,no = "")),
                  # max.overlaps = Inf,
                  colour = "red")+
  labs(title = "AZD6738")+ 
  theme_bw()+
  theme(strip.background = element_rect(fill = "white"))
dev.off()
#   c. How would you define the consensus genes?
unique(intersect((dat%>%filter(normZ <= -2.5, FDR < 0.1, treatment == "AZD6738"))$GENE,
          (dat%>%filter(normZ <= -2.5, FDR < 0.1, treatment == "VE821"))$GENE))
#       i. Do you expect any cell line specific effects?
    # Yes. The RPE1-hTER cell line have many more genes present high sythetic lethality in both ATRi
    # treatments, with common genes such as POLE4 and POLE3. The HCT116 and Hela cell lines have 
    # fewer significant hits with high lethality in ATRi treatment.
pdf("Visualisation_of_top_genes_in_different_cell_lines.pdf",width = 15,height = 10)
ggplot(data = dat.hus%>%
         mutate(colour = ifelse(test = normZ <= -2.5&FDR < 0.1,yes = "red",no = "black")), 
       aes(x = normZ,y = -log10(Pval)))+
  facet_grid(cell_line~treatment)+
  geom_point(aes(colour = colour))+
  scale_color_identity()+
  geom_text_repel(aes(label = ifelse(test = normZ <= -2.5&FDR < 0.1,yes = GENE,no = "")),
                  max.overlaps = 50,
                  colour = "red")+
  theme_bw()+
  theme(strip.background = element_rect(fill = "white"))
dev.off()
#   d. What other questions could you explore?
      # Explore the functional pathways that the consensus hits of the ATRi
      # treatments could be involved in.
      # Explore the differences in genetic backgrounds between HCT116, Hela 
      # and RPE1-hTER lines that could be responsible for the cell line 
      # specific effects reponding to ATRi treatments.







