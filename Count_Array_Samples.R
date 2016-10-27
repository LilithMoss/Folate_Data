################################################
# Read in, count, compare
# samples in all 4 arrays (Illumina, Axiom, )
################################################
#Disable Scientific Notation
options(scipen=999)
#options(scipen=0) #This is the default
#Local Paths
il_path <- "C:/Users/Lilith Moss/Documents/MATH/RESEARCH/Folate/Arrays/illumina/"
il_imp_path <- "C:/Users/Lilith Moss/Documents/MATH/RESEARCH/Folate/Arrays/illumina/imputed/"
ax_path <- "C:/Users/Lilith Moss/Documents/MATH/RESEARCH/Folate/Arrays/axiom/"
on_path <- "C:/Users/Lilith Moss/Documents/MATH/RESEARCH/Folate/Arrays/onco/"
#HPC Paths
# il_path <- "/home/pmd-01/chemenya/FOLATE/DATA/FAM/illumina/"
# il_imp_path <- "/home/pmd-01/chemenya/FOLATE/DATA/FAM/illumina/imputed/"
# ax_path <- "/home/pmd-01/chemenya/FOLATE/DATA/FAM/axiom/"
# on_path <- "/home/pmd-01/chemenya/FOLATE/DATA/FAM/onco/"

#Illumina
#Typed
list.files(il_path)
oneM <- read.table(paste0(il_path,"ccfr_1M_20091211_fwd.fam"))
oneMduo <- read.table(paste0(il_path,"ccfr_1Mduo_20091211_fwd.fam"))
oneStage2A <- read.table(paste0(il_path,"ccfr_Stage2A_20100421_fwd.fam"))
#Imputed
list.files(il_imp_path)
imp.oneMoneMduo <- read.table(paste0(il_imp_path,"ccfr_1M-1Mduo_impute_20120614.fam"))
imp.omni.cc <- read.table(paste0(il_imp_path,"ccfr_Omni_cases-and-controls_impute_20120614.fam"))
imp.omni.cases <- read.table(paste0(il_imp_path,"ccfr_Omni_cases_impute_20120614.fam"))

#Axiom
list.files(ax_path)
axiom <- read.table(paste0(ax_path,"axiom_corect_mecc_cfr_ky_pre-impute_2013-01-23.fam"))

#Onco
list.files(on_path)
onco <- read.table(paste0(on_path,"corect_oncoarray_quality-controlled_genotypes_v3_2016-03-29.CFR.fam"),header=F)
names(onco) <- c("FID","IID","PAT","MAT","SEX","PHEN")
drops <- read.table(paste0(on_path,"oncoarray_sample_qc_drops_2016-05-19.CFR.txt"),header=T)
phen <- read.csv(paste0(on_path,"phen_trunc.csv"))
#Test ID 110036015634
"110036015634" %in% phen.trunc$subjectid
oncom <- merge(onco,phen,by="IID",all.x=T)
onco2 <- oncom[,c(1,9,3:6)]
sum(drops$subject_id %in% onco2$IID)

#Biospecimen data with case/control/ineligible status
bio <- read.csv(paste0(on_path,"bio_status.csv"))
unique(bio$Status)
sum(bio$Status=="Ineligible")

#Count how many in each
oneMs <- nrow(oneM)
oneMduos <- nrow(oneMduo)
oneStage2As <- nrow(oneStage2A)
axioms <- nrow(axiom)
oncos <- nrow(oncom)
naive.totals <- cbind(oneMs,oneMduos,oneStage2As,axioms,oncos)

#Extract names only
names(oneM) <- c("FID","IID","PAT","MAT","SEX","PHEN")
names(oneMduo) <- c("FID","IID","PAT","MAT","SEX","PHEN")
names(oneStage2A) <- c("FID","IID","PAT","MAT","SEX","PHEN")
names(axiom) <- c("FID","IID","PAT","MAT","SEX","PHEN")
names(onco2) <- c("FID","IID","PAT","MAT","SEX","PHEN")

#Chip column
oneM$chip <- rep("oneM",nrow(oneM))
oneMduo$chip <- rep("oneMduo",nrow(oneMduo))
oneStage2A$chip <- rep("oneStage2A",nrow(oneStage2A))
axiom$chip <- rep("axiom",nrow(axiom))
onco2$chip <- rep("onco",nrow(onco2))

#Imputed column - Do this last
oneM$chip <- 
oneMduo$chip <- 
oneStage2A$chip <- 
axiom$imputed <- rep("Y",nrow(axiom))
onco2$imputed <- rep("Y",nrow(onco2))

#Status column
oneM_ <- merge(oneM,bio[,3:4],by="IID",all.x=T)
oneMduo_ <- merge(oneMduo,bio[,3:4],by="IID",all.x=T)
oneStage2A_ <- merge(oneStage2A,bio[,3:4],by="IID",all.x=T)
axiom_ <- merge(axiom,bio[,3:4],by="IID",all.x=T)
onco2_ <- merge(onco2,bio[,3:4],by="IID",all.x=T)

#Drop onco samples from drop list
onco2_$drop <- ifelse(onco2_$IID %in% drops$subject_id,1,0)
onco2_ <- onco2_[onco2_$drop==0,]
onco2_ <- onco2_[,1:8]

#Drop ineligibles from all the arrays
#No need to do this as none have ineligibles

#Combine and clear up case-status
all <- rbind(oneM_,oneMduo_,oneStage2A_,axiom_,onco2_)
all$Status <- as.character(all$Status)
all$Status[is.na(all$Status)] <- "NotIn"
all$nom <- ifelse(all$PHEN==1 & all$Status=="Case",1,0) #8 Misclassified controls that are actually cases
all$status.new <- ifelse(all$PHEN==2,"Case","Control")
#Replace case status if incident case
for(i in 1:nrow(all)){
  if(all[i,8]=="Case"){
    all[i,10] <- "Case"
  }
}

for(i in 1:nrow(all)){
  if(all[i,8]=="Ineligible"){
    all[i,10] <- "Ineligible"
  }
}

#Create unduplicated Data
allu <- all
count <- as.data.frame(table(allu$IID))
names(count) <- c("IID","Freq")
allum <- merge(allu,count,by="IID",all.x=T)
rep <- allum[allum$Freq>1,]

#Who is ineligible?
inel <- allum[allum$Status=="Ineligible",]
allum2 <- allum[allum$Status!="Ineligible",]
rep <- allum2[allum2$Freq>1,]

#Output numbers
table(allum2$chip,allum2$status.new)
table(allum2$chip)

#Get Unique Numbers
uniq <- allum2[allum2$Freq==1,]
rep <- rep[order(rep$IID),]

table(uniq$chip,uniq$status.new)
table(uniq$chip)

table(rep$status.new)
table(uniq$status.new)


