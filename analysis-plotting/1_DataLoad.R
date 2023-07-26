
setwd("/AWI_MPI/FRAM/CTD/Rstats")
setwd("/Users/mwietz/ownCloud - mwietz@owncloud.mpi-bremen.de/AWI_MPI/FRAM/CTD/Rstats")
load("CTD.Rdata")

#save.image("CTD.Rdata")

# Load packages and colors
library(gtools)
library(dplyr)
library(tibble)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggfortify)
library(ggrepel)
library(data.table)
library(R.matlab)
library(ShortRead)
library(vegan)
library(ampvis2)
library(mixOmics)
library(amap)
library(ape)
library(psych)
library(iNEXT)
library(olsrr)
library(fishualize)
library(cowplot)


#######################################################
  ### LOAD ASVs + TAX ###
#######################################################

# Load ASV table
ASV <- read.table(
  "seqtab.txt",
  h = T, sep = "\t",
  check.names = F)

# Load taxonomy table
TAX <- read.table(
  "tax.txt",
  h = T, sep = "\t")

# Remove likely contaminants
# same as done for RAS data
TAX <- TAX %>%
  filter(!Family %in% c(
    "Corynebacteriaceae","Bacillaceae",
    "Weeksellaceae","Enterococcaceae",
    "Streptococcaceae","Propionibacteriaceae",
    "Xanthobacteraceae","Staphylococcaceae"))

# remove EUKs, mitochondria and chloroplast
TAX <- TAX[-grep(
  'Mitochondria', TAX$Family),]
TAX <- TAX[-grep(
  'Chloroplast', TAX$Order),]
TAX <- TAX[-grep(
  'Eukaryota', TAX$Kingdom),]

# match TAX after contaminant removal
ASV <- ASV[row.names(TAX),]

# Set abundance filter: >=3 counts in >=3 samples
filter <- apply(ASV, 1, function(row) {
  sum(row >= 3) >= 3})

# Filter ASVs accordingly
ASV <- ASV[filter, ]

# match TAX after abundance filter
TAX <- TAX[row.names(ASV),]

#################################################

# Rename NAs with last known taxrank + "uc"
k <- ncol(TAX)-1
for (i in 2:k) {
if (sum(is.na(TAX[, i])) >1) {
  temp <- TAX[is.na(TAX[, i]), ]
  for (j in 1:nrow(temp)) {
    if (sum(is.na(
      temp[j, i:(k+1)])) == length(temp[j, i:(k+1)])) {
      temp[j, i] <- paste(temp[j, (i-1)], " uc", sep = "")
      temp[j, (i+1):(k+1)] <- temp[j, i]
  }
}
  TAX[is.na(TAX[, i]), ] <- temp}
  if (sum(is.na(TAX[, i]))==1) {
    temp <- TAX[is.na(TAX[, i]), ]
    if (sum(is.na(temp[i:(k+1)])) == length(temp[i:(k+1)])) {
    temp[i] <- paste(temp[(i-1)], " uc", sep="")
    temp[(i+1):(k+1)] <- temp[i]
    }
    TAX[is.na(TAX[, i]),] <- temp
  }
}
TAX[is.na(TAX[, (k+1)]), (k+1)] <- paste(
  TAX[is.na(TAX[, (k+1)]), k], " uc", sep="")

# shorten/modify names
TAX <- TAX %>%
  mutate(across(everything(),~gsub("Clade_","SAR11 Clade ", .))) %>%
  mutate(across(everything(),~gsub("_clade","", .))) %>%
  mutate(across(everything(),~gsub("Candidatus","Cand", .))) %>%
  mutate(across(everything(),~gsub("Roseobacter_NAC11-7_lineage","NAC11-7", .))) %>%
  mutate(across(everything(),~gsub("_marine_group","", .))) %>%
  mutate(across(everything(),~gsub("_terrestrial_group","", .))) %>%
  mutate(across(everything(),~gsub("_CC9902","", .))) %>%
  mutate(across(everything(),~gsub("(Marine_group_B)","", ., fixed=T))) %>%
  mutate(across(everything(),~gsub("(SAR406)","SAR406", ., fixed=T)))

# Add dummy *Species* column (needed for ampvis)
TAX$Species <- TAX$Genus


#######################################################
   ###  SAMPLE INFO / CTD + SUBSTRATE DATA
#######################################################

# load metadata
# remove failed seq-sample
# add longitudinal ranges
# Add PangaeaShort (without subsample) for merging
ENV <- read.table(
  "metadata.txt", h=T, sep = "\t", 
  stringsAsFactors=F, skipNul=T) %>%
  mutate(date=as.Date(
    date, "%Y-%m-%d")) %>%
  mutate(jday=as.numeric(
    format(date, "%j"))) %>%
  mutate(lon_short=cut(lon, c(seq(
    -6, 12, 3), Inf), labels=c(
      "6-3","3-0","0-3","3-6",
      "6-9","9-12","12-inf"))) %>%
  filter(sample_title!="2017-N5-100") %>%
  mutate(PangaeaShort=gsub("-.*","",PangaeaEvent))
  
# Factorize parameters of interest
# Add "year" column
ENV$region <- factor(ENV$region, levels=c(
  "western Fram Strait","eastern Fram Strait"))
ENV$layer <- factor(ENV$layer, levels=c(
  "surface","chl-max","bel.chl-max","lower-photic"))
ENV$lon_short <- factor(ENV$lon_short, levels=c(
  "6-3","3-0","0-3","3-6","6-9","9-12"))
ENV$year <- factor(format(
  ENV$date, "%Y"))

# Load fractionated chlorophyll
# Add PangaeaShort (without subsample) for merging
fracChl <- read.csv(
  "fracChl.txt", h=T, sep = "\t", check.names=F,
  stringsAsFactors=F, skipNul=T) %>%
  mutate(date=as.Date(
    date, "%Y-%m-%d")) %>%
  mutate(PangaeaShort = gsub("-.*","",PangaeaEvent)) %>%
  dplyr::select(-c("date","lat","lon","PangaeaEvent"))

# Load biogeo data / cell numbers
biogeo <- read.csv(
  "biogeochemistry.txt", h=T, sep="\t", check.names=F,
  stringsAsFactors=F, skipNul=T)

# Average depth of layers
avgDepth <- ENV %>%
  group_by(layer) %>%
  summarise(mean(depth))
  

#######################################################
   ###  SATELLITE ICE + CHLOROPHYLL
#######################################################

# Read Matlab into R
ice_chl <- readMat(
  "ice_chl.mat") 

# unnest list
ice_chl <- lapply(rapply(
  ice_chl, enquote, how="unlist"), eval)

# Rename
names(ice_chl) <- c(
    "date","lat","lon","chl_conc","chl_date",
    "ice_dist","ice_conc","ice_date") 

# Extract CTD sampling dates; convert from Matlab 
dates <- data.frame(ice_chl$date) %>%
  mutate(across(everything(), ~as.POSIXct((
   . - 719529)*86400, origin="1970-01-01", tz="UTC"))) %>%
  dplyr::rename_at(1, ~"date") 

##########################

## ICE 
# concentration and distance to ice edge 
# available for each sampling date + weeks before/after

Ice <- ice_chl[grepl("ice", names(ice_chl))] 

# Keep only columns until actual sampling date (col #31)
Ice <- lapply(Ice, function(x) x[, c(1:31)])

# To individual DFs
IceDist <- data.frame(Ice$ice_dist)
IceConc <- data.frame(Ice$ice_conc)

# Calculate past ice conditions
# Add dates & lat
IceConcPast <- IceConc %>%
  mutate(IceConcPast = rowMeans(.)) %>%
  dplyr::select(c("IceConcPast")) %>%
  cbind(dates, lat = ice_chl$lat) 
IceDistPast <- IceDist %>%
  mutate(IceDistPast = rowMeans(.)) %>%
  dplyr::select(c("IceDistPast")) %>%
  cbind(dates, lat = ice_chl$lat) 

# Reformat: only actual sampling date
IceConc <- IceConc %>%
  dplyr::select(c(31)) %>%
  dplyr::rename(iceConc = X31) %>%
  cbind(dates, lat = ice_chl$lat)
IceDist <- IceDist %>%
  dplyr::select(c(31)) %>%
  dplyr::rename(iceDist = X31) %>%
  cbind(dates, lat = ice_chl$lat)

##############################

## CHLOROPHYLL
# available for each sampling date + weeks before/after

Chl <- ice_chl[grepl("chl", names(ice_chl))] 

# Keep only columns until actual sampling date (col #6)
Chl <- lapply(Chl, function(x) x[, c(1:6)])

# To individual DF
ChlConc <- data.frame(Chl$chl_conc)

# Calculate past Chl concentrations
# Add dates & lat
ChlPast <- ChlConc %>%
  mutate(chlPast = rowMeans(.)) %>%
  dplyr::select(c("chlPast")) %>%
  cbind(dates, lat = ice_chl$lat) 

# Reformat: only actual sampling date
ChlConc <- ChlConc %>%
  dplyr::select(c(6)) %>%
  dplyr::rename(chlConc = X6) %>%
  cbind(dates, lat = ice_chl$lat)

###############################

# Append all ice + chl data to ENV
# convert NaN to NA
ENV <- ENV %>%
 left_join(biogeo, by=c("sample_title")) %>%
 left_join(fracChl) %>%
 left_join(IceConc, by=c("date","lat")) %>%
 left_join(IceDist, by=c("date","lat")) %>%
 left_join(IceConcPast, by=c("date","lat")) %>%
 left_join(IceDistPast, by=c("date","lat")) %>%
 left_join(ChlConc, by=c("date","lat")) %>%
 left_join(ChlPast, by=c("date","lat")) %>%
 distinct() %>%
 mutate_if(is.numeric, ~ifelse(is.nan(.), NA, .)) 


#######################################################
   ###  PROCESSING / FORMATTING
#######################################################

# sort; match ASV+ENV 
ASV <- ASV[,mixedsort(names(ASV))]
ENV <- ENV[mixedorder(ENV$seqname),]
ASV <- ASV[,c((match(ENV$seqname, 
  colnames(ASV))))]

# rename 
colnames(ASV) = ENV$sample_title

# Verify same order of ENV and ASV tables
ENV <- ENV[mixedorder(ENV$sample_title),]
ASV <- ASV[, ENV$sample_title]

#####################################

# Relative abundances
ASV.rel = as.data.frame(
  apply(ASV, 2, function(x) x / sum(x) * 100)) 

# Hellinger-transform
ASV.hel = as.data.frame(
  apply(ASV, 2, function(x) sqrt(x / sum(x))))

#####################################

# Ampvis load
ampvis <- data.frame(
  OTU = rownames(ASV),
  ASV, TAX, check.names=F)
ampvis <- amp_load(ampvis, ENV)

# Export for SI table
ASV.rel %>% mutate_if(is.numeric, round, 2) %>%
  write.table(file="../ASV_rel.txt",sep="\t", 
              row.names=T, col.names=T, quote=F)

# Export for SI table
ENV %>% mutate_if(is.numeric, round, 2) %>%
  write.table(file="../metadata.txt",sep="\t", 
              row.names=F, col.names=T, quote=F)


#######################################################
  ###  ALPHA-DIVERSITY
#######################################################

# calculate on full set (i.e. without subsetting)
# run on AWI-VM for speed
setwd("/isibhv/projects/FRAMdata/MolObs/WaterCol_CTD_1519/")

read.table(
  "seqtab.txt",
  h=T, sep="\t", check.names=F) %>%
  iNEXT(
    datatype="abundance", 
    q=c(0), conf=0.95, nboot=100) -> iNEXT

# save+export
save(iNEXT, file="iNEXT.Rdata")

###################################

# Rest done on local R

# extract from iNEXT object
load("iNEXT.Rdata")

# Extract diversity indicators
richness <- iNEXT$AsyEst[
  iNEXT$AsyEst$Diversity=="Species richness",] %>%
  arrange(Site) 
simpson <- iNEXT$AsyEst[
  iNEXT$AsyEst$Diversity=="Simpson diversity",] %>%
  arrange(Site) 

# compile
AlphaDiv <- data.frame(
  seqname = richness$Site,
  richness = richness$Observed,
  simpson = simpson$Observed) %>%
  left_join(ENV)


#######################################################
  ### DISTANCES
#######################################################

# Distance matrix -- all taxa
dist.asv <- t(ASV.hel) %>% as.data.frame
dist.asv <- vegdist(dist.asv, method="bray") 

# Distance matrix -- BACTERIA only
dist.bac <- ASV.hel %>% 
  filter(row.names(.) %in% row.names(TAX)[grep(
    "Bacteria", TAX$Kingdom)]) %>% as.data.frame()
dist.bac <- t(dist.bac) %>% as.data.frame
dist.bac <- vegdist(dist.bac, method="bray") 

# Distance matrix -- ARCHAEA only
dist.arch <- ASV.hel %>% 
  filter(row.names(.) %in% row.names(TAX)[grep(
    "Archaea", TAX$Kingdom)]) %>% as.data.frame()
dist.arch <- t(dist.arch) %>% as.data.frame
dist.arch <- vegdist(dist.arch, method="bray") 

###############################

# Distance matrix - Substrates
dist.env = ENV 
rownames(dist.env) <- NULL
dist.env <- dist.env %>%
  column_to_rownames("sample_title") %>%
  #filter(sample_title!="2018-HG4-35") %>%
  select_if(is.numeric) %>%
  dplyr::select(-c(
    "lat","lon","O2_sat","TDN","O2_conc","temp",
    "depth","iceDist","chlPast","chlConc",
    "chl", "chl_0.4-3um","chl_larger3um",
    "iceConc","IceConcPast","IceDistPast",
    #,"DHAA-C","DHAA-N","year","jday",
    "sal", "Bacterial cells")) %>% 
  drop_na() 
dist.env <- vegdist(dist.env, method="euclidean")  %>%
  as.matrix(dist.env)  %>% as.data.frame() 
  

############################################################################################
  ###  Prepare for ASV-MAG mapping
############################################################################################

# Set abundance filter
filter <- apply(ASV.rel, 1, function(row) {
  sum(row >= 0.01) >= 3})

# Filter accordingly
tax <- TAX[filter, ] %>%
  rownames_to_column("asv")

# read fasta
asv.fa <- as.data.frame(sread(readFasta(
  "uniques.fasta"))) %>%
  add_column(asv=1:nrow(.)) %>%
  mutate(asv=paste0("asv", asv)) %>%
  right_join(tax) %>%
  mutate(asv=paste0(">", asv)) %>%
  rowwise %>%
  mutate(id = paste(sort(c(
    asv, Genus)), collapse="_")) %>%
  mutate(across(c(id), ~gsub(
    "bac_|_uc|_clade", "", .))) 

# export combined data
setNames(as.character(
  asv.fa$x), asv.fa$id) %>%
  write.table(
    ., file="../ASV_seqs.fasta",sep="\n", 
    row.names=T, col.names=F, quote=F)


############################################################################################
###  Export for FRAM-eDNA Metfies et al ###
############################################################################################

ASV.rel %>% write.table(
  file="../../../collaborations/eDNA-FRAM/ASV-ctd.txt",
  sep="\t", row.names=T, col.names=T, quote=F)
TAX %>% write.table(
  file="../../../collaborations/eDNA-FRAM/TAX-ctd.txt",
  sep="\t", row.names=T, col.names=T, quote=F)
ENV %>% write.table(
  file="../../../collaborations/eDNA-FRAM/meta-ctd.txt",
  sep="\t", row.names=F, col.names=T, quote=F)


############################################################

## modified theme_classic
theme_classic2 <- function(
  base_size=12, base_family=""){
  theme_bw(base_size=base_size, base_family=base_family) %+replace%
    theme(panel.border = element_blank(),
          axis.line = element_line(colour="black"),
          panel.grid.major = element_line(),
          #panel.grid.major.x = element_blank(),
          #panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          strip.background = element_rect(colour="black", size = 0.5),
          legend.key = element_blank()
    )
}

# remove temporary datasets
rm(filter,files,temp,i,j,k,asv,tax)

# save data
save.image("CTD.Rdata")
