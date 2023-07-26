


############################################################################################
   ###  RAREFACTION AND COVERAGE  ###
############################################################################################

# calculate on full set (i.e. without subsetting)
# run on AWI-VM for speed
setwd("/isibhv/projects/FRAMdata/MolObs/WaterCol_CTD_1519/")

read.table(
  "seqtab.txt",
  h=T, sep="\t", check.names=F) %>%
  iNEXT(
    datatype="abundance", 
    q=c(0), conf=0.95, nboot=100) -> iNEXT

# save+export; rest done on local R
save(iNEXT, file="iNEXT.Rdata")

###################################

# extract form inext, reforam for mergin

load("iNEXT.Rdata")

rarefac <- do.call(rbind, iNEXT$iNextEst) %>%
  rownames_to_column("seqname") %>%
  mutate(seqname = gsub("\\.\\d{1,2}","", seqname))

rarefac.point <- rarefac[which(
  rarefac$method == "observed"),]
rarefac.line <- rarefac[which(
  rarefac$method != "observed"),]
rarefac.line$method <- factor(
  rarefac.line$method,
    c("interpolated","extrapolated"),
    c("interpolation","extrapolation"))

# Add sorting columns
rarefac <- rarefac %>%
  left_join(dplyr::select(
    ENV, seqname, region, Station, layer)) %>%
  distinct()
rarefac.line <- rarefac.line %>%
  left_join(dplyr::select(
    ENV, seqname, region, Station, layer)) %>%
  distinct()

###################################

cover <- fortify(
  iNEXT, type=2)

cover.point <- cover [which(
  cover$method == "observed"),]
cover.line <- cover [which(
  cover$method != "observed"),]
cover.line$method <- factor(
  cover.line$method,
    c("interpolated","extrapolated"),
    c("interpolation","extrapolation"))

# add sorting columns
cover <- cover %>%
  left_join(dplyr::select(
    ENV, sample_title, 
    region, Station, layer), 
     by=c("site"="sample_title")) %>%
  distinct()

cover.line <- cover.line %>%
  left_join(dplyr::select(
    ENV, sample_title, 
    region, Station, layer), 
    by=c("site"="sample_title")) %>%
  distinct()

###################################

rarefac  %>% 
#left_join(ENV, by=c(
#  "site"="sample_title",
#  "region"="region")) %>% 
#  mutate(region2=case_when(lat >0~"WSC",TRUE~"EGC")) %>%
ggplot(aes(
  x=x, y=y, colour=site)) +
geom_line(aes(
  linetype = method), lwd = 0.5, 
  data = rarefac.line) +
scale_x_continuous(limits = c(0,1e+5)) +
  scale_color_fish_d(option="Cirrhilabrus_solorensis") +
facet_grid(layer~region) +
labs(x="Sample size", 
     y="Species richness") +
theme_bw() + 
theme(
    axis.text.x=element_text(
      angle=90, hjust=1, vjust=0.5),
    legend.position="none")

cover %>%
  left_join(ENV, by=c(
    "site"="sample_title",
    "region"="region")) %>% 
  mutate(region2=case_when(lat >0~"WSC",TRUE~"EGC")) %>%
ggplot(aes(x=x, y=y, colour=site))+ 
geom_line(aes(
  linetype = method), 
  lwd = 0.5, data=cover.line) +
scale_x_continuous(
  limits = c(0,1e+5)) +
scale_y_continuous(
  breaks = seq(0.9,1,0.05), 
  limits = c(0.9,1)) +
#scale_color_manual(
 # values=scico::scico(205, palette='roma'),
 # guide="none") +
  scale_color_fish_d(option="Cirrhilabrus_solorensis") +
facet_grid(layer~region2) +
labs(x="Sample size", y="Sample coverage") +
theme_bw() + 
theme(axis.text.x=element_text(
        angle=90, hjust=1, vjust=0.5),
      legend.position="none")


###########################################################################################
   ###  ALPHA-DIV  ###
############################################################################################

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

###################################

