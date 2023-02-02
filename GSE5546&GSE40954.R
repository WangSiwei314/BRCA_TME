options(stringsAsFactors = F)
library(survminer);library(survival);
library(openxlsx)
library(magrittr); library(tidyr);library(tibble);library(dplyr);
library(ggplot2);library(ggpubr);library(patchwork); 

# GSE5546 ---------------------------------------------------------------------------------------------------------
GSE5546_Anno <- read.xlsx("./GSE5546/GSE5546.xlsx",sheet = 1); 
rownames(GSE5546_Anno) <- GSE5546_Anno$GEO.accession;

GSE5546_Exp <- read.xlsx("./GSE5546/GSE5546.xlsx",sheet = 2);
GSE5546_Exp <- GSE5546_Exp[GSE5546_Exp$GENE_SYMBOL %in% Cytokines,]
GSE5546_Exp[is.na(GSE5546_Exp)] <- 0;
GSE5546_Exp <- group_by(GSE5546_Exp[,-1], GENE_SYMBOL) %>% summarise_all(., max) %>% as.data.frame(.) %>% column_to_rownames(.,"GENE_SYMBOL")

GSE5546_Anno$RS <- apply(GSE5546_Exp,2,function(x){0.359*x[1] - 0.324*x[2] - 0.367*x[3] - 0.328*x[4] - 0.343*x[5]})[rownames(GSE5546_Anno)]
GSE5546_Anno$DT <- paste(GSE5546_Anno$TP53_Status, ifelse(GSE5546_Anno$RS > median(GSE5546_Anno$RS), "High", "Low"), sep = "_") 
GSE5546_Anno$DT <- ifelse(GSE5546_Anno$DT == "mt_High", "DT", "Other")
  
ggsurvplot(survfit(Surv(Time.survival, Survival == "1") ~ TP53_Status, data = GSE5546_Anno), 
           data = GSE5546_Anno, conf.int = F, pval = T, legend.title = "HR: ")

ggsurvplot(survfit(Surv(Time.recurrence, Recurrence == "1") ~ TP53_Status, data = GSE5546_Anno), 
           data = GSE5546_Anno, conf.int = F, pval = T, legend.title = "HR: ")
ggsurvplot(survfit(Surv(Time.recurrence, Recurrence == "1") ~ DT, data = GSE5546_Anno), 
           data = GSE5546_Anno, conf.int = F, pval = T, legend.title = "HR: ")


# GSE40954 ---------------------------------------------------------------------------------------------------------
GSE40954_Anno <- read.xlsx("./GSE40954/GSE40954.xlsx",sheet = 3,rowNames = TRUE);

GSE40954_Exp1 <- read.xlsx("./GSE40954/GSE40954.xlsx",sheet = 1, startRow = 2);
GSE40954_Exp1 <- GSE40954_Exp1[GSE40954_Exp1$Gene_Symbol %in% Cytokines, -1]
GSE40954_Exp1[GSE40954_Exp1 == "null"] <- 0;
GSE40954_Exp1[,2:46] <- as.data.frame(apply(GSE40954_Exp1[,-1],2,as.numeric))
GSE40954_Exp1 <- group_by(GSE40954_Exp1, Gene_Symbol) %>% summarise_all(., mean) %>% as.data.frame(.) %>% column_to_rownames(.,"Gene_Symbol")
colnames(GSE40954_Exp1) <- gsub("-.-"," - ",colnames(GSE40954_Exp1))
GSE40954_Exp1 <- apply(GSE40954_Exp1,2,function(x){0.359*x[1] - 0.324*x[2] - 0.367*x[3] - 0.328*x[4] - 0.343*x[5]})

GSE40954_Exp2 <- read.xlsx("./GSE40954/GSE40954.xlsx",sheet = 2, startRow = 2);
GSE40954_Exp2 <- GSE40954_Exp2[GSE40954_Exp2$Gene_Symbol %in% Cytokines, -1]
GSE40954_Exp2[GSE40954_Exp2 == "null"] <- 0;
GSE40954_Exp2[,2:37] <- as.data.frame(apply(GSE40954_Exp2[,-1],2,as.numeric))
GSE40954_Exp2 <- group_by(GSE40954_Exp2, Gene_Symbol) %>% summarise_all(., mean) %>% as.data.frame(.) %>% column_to_rownames(.,"Gene_Symbol")
colnames(GSE40954_Exp2) <- gsub("-.-"," - ",colnames(GSE40954_Exp2))
GSE40954_Exp2 <- apply(GSE40954_Exp2,2,function(x){0.359*x[1]- 0.367*x[2]})

GSE40954_Anno$RS <- c(GSE40954_Exp1,GSE40954_Exp2)[rownames(GSE40954_Anno)]
GSE40954_Anno <- GSE40954_Anno[!is.na(GSE40954_Anno$RS),];
GSE40954_Anno$DT <- paste(GSE40954_Anno$TP53, ifelse(GSE40954_Anno$RS > median(GSE40954_Anno$RS), "High", "Low"), sep = "_") 
GSE40954_Anno$DT <- ifelse(GSE40954_Anno$DT == "1_High", "DT", "Other")

ggsurvplot(survfit(Surv(RFS.time, RFS.status == "1") ~ TP53, data = GSE40954_Anno), 
           data = GSE40954_Anno, conf.int = F, pval = T, legend.title = "HR: ") # 3*3.5 GSE40954_TP53_Surv

ggsurvplot(survfit(Surv(RFS.time, RFS.status == "1") ~ DT, data = GSE40954_Anno), 
           data = GSE40954_Anno, conf.int = F, pval = T, legend.title = "HR: ") # 3*3.5 GSE40954_TP53RS_Surv
















