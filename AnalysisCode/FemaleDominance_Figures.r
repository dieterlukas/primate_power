library(viridis)
library(ggtree)
library(caper)
library(tidyr)
library(ggstance)
library(ggdist)
library(cowplot)
library(ggplot2)


# Colour schemes:

# Strict three-way classification of dominance:
  
male_dominance_color<-"#31688EFF"
  # Co dominance 2: viridis(3)[2] "#21908CFF"
co_dominance_color<-"#A6D854"
  # Female dominance 3: viridis(3)[1] "#440154FF"
female_dominance_color<-"#FC8D62"
  dominance_colors<-c("#FC8D62","#A6D854", "#31688EFF" )
  
  alt_dominance_colors<-c("steelblue","grey80","orange")


dominance_colors1<-c("#601A4A",  "#63ACBE","#EE442F")
dominance_colors2<-c("#31688EFF","#8FD744FF","#D55E00")
dominance_colors3<-c("#56B4E9","#E69F00","#CC79A7")
dominance_colors4<-c("#2271B2", "#E69F00","#359B73")
dominance_colors5<-c("#FC8D62","#A6D854", "#8DA0CB" )


# sex-specific distribution of aggression
ff_aggression<- "#FDE725FF"
mm_aggression<- "#443A83FF"
fm_aggression<- "#21908CFF"


palette1<- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

palette2<- c("#440154FF", "#443A83FF" ,"#31688EFF", "#21908CFF", "#35B779FF", "#8FD744FF", "#FDE725FF")






################################################################################
################################################################################
# Figure 1: phylogenetic distribution of sex-specific aggression and dominance

fightswondata<-as.data.frame((combined %>% group_by(corrected_species_id) %>% summarize(mean(perc_won_females)) ))
colnames(fightswondata)<-c("Species","perc_won_females")
rownames(fightswondata)<-fightswondata$Species

combined$numericalstrictdom<-1
combined[combined$strictfdom=="2",]$numericalstrictdom<-2
combined[combined$strictfdom=="3",]$numericalstrictdom<-3
dominancedata<-as.data.frame((combined %>% group_by(corrected_species_id) %>% summarize(mean(numericalstrictdom)) ))
colnames(dominancedata)<-c("Species","strictfdom")
rownames(dominancedata)<-dominancedata$Species
dominancedata$strictfdom<-round(dominancedata$strictfdom)

combined$numericalrelaxeddom<-1
combined[combined$strictfdom=="3",]$numericalrelaxeddom<-2

relaxeddominancedata<-as.data.frame((combined %>% group_by(corrected_species_id) %>% summarize(mean(numericalrelaxeddom)) ))
colnames(relaxeddominancedata)<-c("Species","relaxedfdom")
rownames(relaxeddominancedata)<-relaxeddominancedata$Species
relaxeddominancedata$relaxedfdom<-round(relaxeddominancedata$relaxedfdom)


missing<-treedata(inputtree,data=fightswondata,warnings=FALSE)
# We remove the species from the tree for which we have no data
perctree<-missing$phy
# We remove the species from the dataset which are not in the tree
speciesnames<-perctree$tip.label
percdata<-data.frame(fightswondata[speciesnames,])
colnames(percdata)<-"perc_won_females"
row.names(percdata)<-speciesnames

# Plot the strict 3-way classification of intersexual dominance on a phylogeny, including reconstructed states

strict_females_phylodata_1<-as.data.frame(combined %>% group_by(corrected_species_id) %>% summarise(mean(as.numeric(strictfdom),na.rm=T )))
strict_females_phylodata_2<-as.data.frame(strict_females_phylodata_1[,2])
row.names(strict_females_phylodata_2)<-strict_females_phylodata_1[,1]
colnames(strict_females_phylodata_2)<-"strictfdom"
strict_females_phylodata_3<-as.data.frame(strict_females_phylodata_2[is.na(strict_females_phylodata_2$strictfdom)==FALSE,])
row.names(strict_females_phylodata_3)<-row.names(strict_females_phylodata_2)[is.na(strict_females_phylodata_2$strictfdom)==FALSE]
colnames(strict_females_phylodata_3)<-"strictfdom"
missing<-treedata(inputtree,data=strict_females_phylodata_3,warnings=FALSE)
stricttree<-missing$phy
speciesnames<-stricttree$tip.label
strictdata<-data.frame(strict_females_phylodata_3[speciesnames,])
colnames(strictdata)<-"strictfdom"
row.names(strictdata)<-speciesnames

phylostrictdata<-round(strictdata,0)
phylostrictdata<-as.data.frame(phylostrictdata)
grafentree<-compute.brlen(stricttree,method="Grafen")
onetree<-compute.brlen(stricttree,1)
ancestraldominance<-ace(phylostrictdata$strictfdom,stricttree,type="discrete",model="ARD")


xvalues<-rep(2,nrow(phylostrictdata))
xvalues[phylostrictdata$strictfdom ==1]<-3
xvalues[phylostrictdata$strictfdom ==3]<-1
names(xvalues)<-rownames(phylostrictdata)
obj<-contMap(stricttree,xvalues,plot=FALSE)
obj<-setMap(obj,colors=c(dominance_colors[1],dominance_colors[2],dominance_colors[3]))






treeplot<-ggtree(perctree,ladderize=F) 

p2 <- facet_plot(treeplot, panel="Percentage fights won by females", data=fightswondata, geom=geom_point, aes(x=perc_won_females), cex=3,color=c("black"))

fightswondata_range_max<-as.data.frame((combined %>% group_by(corrected_species_id) %>% summarize(max(perc_won_females)) ))
colnames(fightswondata_range_max)<-c("Species","max_perc_won_females")
rownames(fightswondata_range_max)<-fightswondata_range_max$Species
fightswondata_range_max<-data.frame(fightswondata_range_max[speciesnames,])
fightswondata_range_max$Count<-c(1:nrow(fightswondata_range_max))

fightswondata_range_min<-as.data.frame((combined %>% group_by(corrected_species_id) %>% summarize(min(perc_won_females)) ))
colnames(fightswondata_range_min)<-c("Species","min_perc_won_females")
rownames(fightswondata_range_min)<-fightswondata_range_min$Species
fightswondata_range_min<-data.frame(fightswondata_range_min[speciesnames,])

fightswondata_range_min$Count<-c(1:nrow(fightswondata_range_min))
fightswondata_range_min$delta_long<-fightswondata_range_max$max_perc_won_females- fightswondata_range_min$min_perc_won_females
fightswondata_range_min$delta_lat<-rep(0,nrow(fightswondata_range_min))
fightswondata_range_min[fightswondata_range_min$delta_long %in% 0 ,]$delta_long<-2
fightswondata_range_min[fightswondata_range_min$min_perc_won_females %in% 100 & fightswondata_range_min$delta_long %in% 2,]$min_perc_won_females<-98


fightdistributiondata<-as.data.frame(combined %>% group_by(corrected_species_id) %>% summarize(fmaggression=mean(perc_aggression_fm,na.rm=T),ffaggression=mean(perc_aggression_ff,na.rm=T),mmaggression=mean(perc_aggression_mm,na.rm=T) ))
rownames(fightdistributiondata)<-fightswondata$Species

fightdistributiondata<-data.frame(fightdistributiondata[speciesnames,])
colnames(fightdistributiondata)<-c("Species","between the sexes","among females","among males")
row.names(fightdistributiondata)<-speciesnames

fightdistributiondata_long<-as.data.frame(pivot_longer(fightdistributiondata,cols=c("between the sexes","among females","among males")))
fightdistributiondata_long$numericspecies<-NA
for (i in 1:length(unique(fightdistributiondata_long$Species))){
  fightdistributiondata_long[ ((i-1)*3+1):((i-1)*3+3),]$numericspecies<-i
}


fightdistributiondata_long[is.na(fightdistributiondata_long$value)==T,]$value<-0


p3 <- facet_plot(treeplot, panel="Distribution of aggression", data=fightdistributiondata_long, geom=geom_barh, mapping=aes(fill = name,y = numericspecies, x = value),stat="identity", cex=1,fill=rep(c(fm_aggression,mm_aggression,ff_aggression),length(unique(fightdistributiondata_long$Species))))




dominancedata<-as.data.frame((combined %>% group_by(corrected_species_id) %>% summarize(mean(numericalstrictdom)) ))
colnames(dominancedata)<-c("Species","strictfdom")
rownames(dominancedata)<-dominancedata$Species
dominancedata$strictfdom<-round(dominancedata$strictfdom)
dominancedata<-data.frame(dominancedata[speciesnames,])
dominancedata$colours<-1
dominancedata[dominancedata$strictfdom==1,]$colours<-male_dominance_color
dominancedata[dominancedata$strictfdom==2,]$colours<-co_dominance_color
dominancedata[dominancedata$strictfdom==3,]$colours<-female_dominance_color
dominancedata$placeholder<-1

p4 <- facet_plot(treeplot, panel="Intersexual dominance", data=dominancedata, geom=geom_point, aes(x=placeholder), color=c(dominancedata$colours),pch=15,cex=2)



################################################################################


# Figure 1: three panels with distribution of data


pdf("figures/L_Figure1.pdf")
treeplot+geom_facet(panel="Distribution of aggression", data=fightdistributiondata_long, geom=geom_barh, mapping=aes(fill = name,y = numericspecies, x = value),stat="identity", cex=1,fill=rep(c(fm_aggression,mm_aggression,ff_aggression),121))+geom_facet(panel="Percentage fights won by females", data=fightswondata_range_min, geom=geom_segment, aes(x=min_perc_won_females,xend=min_perc_won_females+delta_long,y=Count,yend=Count+delta_lat), size=2,color=c("black"))+geom_facet(panel="Intersexual dominance", data=dominancedata, geom=geom_point, aes(x=placeholder), cex=2,color=c(dominancedata$colours),pch=15)+theme_tree()
dev.off()

# Figure 1: one panel with phylogeny

pdf("figures/R_Figure1.pdf")
plot(obj,lwd=5,outline=FALSE,direction="leftwards")
dev.off()









################################################################################
################################################################################
################################################################################

# Figure 2 - only raw data
# top row: a) social organization, b) fission-fusion, c) sexual size dimorphism, d) canine size dimorphism,  e) home range overlap, f) number of females
# bottom row:  g) female evictions, h) population origin (ns), i) harshness (ns), j) rainfall seasonality (ns), k) rainfall unpredictability(ns), l) female infanticide (ns), m) relative canine size (ns)

# Figure 3 - only raw data
# single row: a) mating system, b) foraging location, c) sexual receptivity, d) adult sex ratio, e) male reproductive skew, f) receptive synchrony (ns), g) relative testes mass (ns)

# Figure 4 - only raw data
# single row: a) sex bias in dispersal, b) female relatedness (ns), c) female coalitions (ns), d) adult sex ratio, e) male-male conflicts, f) number of males (ns)



################################################################################
##### Figure 2: female competition


summarizedtable<-combined %>%
  group_by(SocOrgPMK,strictfdom) %>%
  summarize(Total = n())

summarizedtable<-as.data.frame(summarizedtable)
summarizedtable[8,]<-c("P",1,0)
summarizedtable[9,]<-c("S",1,0)
colnames(summarizedtable)<-c("SocOrgPMK","StrictFemdom","Observations")
summarizedtable$Observations<-as.integer(summarizedtable$Observations)
summarizedtable[summarizedtable$SocOrgPMK =="G", ]$SocOrgPMK<-"group"
summarizedtable[summarizedtable$SocOrgPMK =="S", ]$SocOrgPMK<-"solitary"
summarizedtable[summarizedtable$SocOrgPMK =="P", ]$SocOrgPMK<-"pair"
summarizedtable[summarizedtable$StrictFemdom ==1, ]$StrictFemdom<-"3) male dominance"
summarizedtable[summarizedtable$StrictFemdom ==2, ]$StrictFemdom<-"2) co-dominance"
summarizedtable[summarizedtable$StrictFemdom ==3, ]$StrictFemdom<-"1) female dominance"
summarizedtable<-summarizedtable[order(summarizedtable$SocOrgPMK,summarizedtable$StrictFemdom),]
summarizedtable$StrictFemdom<-as.factor(summarizedtable$StrictFemdom)

plot_socialsystem <-ggplot(summarizedtable, aes(x = factor(SocOrgPMK,levels=c("group","solitary","pair")), y = Observations, fill = StrictFemdom)) + 
  geom_bar(stat = "identity",fill=c(female_dominance_color,co_dominance_color,male_dominance_color,female_dominance_color,co_dominance_color,male_dominance_color,female_dominance_color,co_dominance_color,male_dominance_color))+
  theme(       axis.text.y=element_blank(),
               axis.ticks.y=element_blank(),
               axis.title.y = element_blank(),
               axis.title.x = element_blank(),
               axis.text=element_text(size=10)
  )

# Fission fusion

summarizedtable<-combined %>%
  group_by(fissionfusion,strictfdom) %>%
  summarize(Total = n())

summarizedtable<-as.data.frame(summarizedtable)
summarizedtable<-summarizedtable[c(1:5),]
summarizedtable[6,]<-c("Yes",3,0)
colnames(summarizedtable)<-c("FissionFusion","StrictFemdom","Observations")
summarizedtable$Observations<-as.integer(summarizedtable$Observations)
summarizedtable[summarizedtable$FissionFusion =="Yes", ]$FissionFusion<-"present"
summarizedtable[summarizedtable$FissionFusion =="No", ]$FissionFusion<-"absent"
summarizedtable[summarizedtable$StrictFemdom ==1, ]$StrictFemdom<-"3) male dominance"
summarizedtable[summarizedtable$StrictFemdom ==2, ]$StrictFemdom<-"2) co-dominance"
summarizedtable[summarizedtable$StrictFemdom ==3, ]$StrictFemdom<-"1) female dominance"
summarizedtable<-summarizedtable[order(summarizedtable$FissionFusion,summarizedtable$StrictFemdom),]
summarizedtable$StrictFemdom<-as.factor(summarizedtable$StrictFemdom)

plot_fissionfusion <-ggplot(summarizedtable, aes(x = factor(FissionFusion,levels=c("absent","present")), y = Observations, fill = StrictFemdom)) + 
  geom_bar(stat = "identity",fill=c(female_dominance_color,co_dominance_color,male_dominance_color,female_dominance_color,co_dominance_color,male_dominance_color))+
  theme(       axis.text.y=element_blank(),
               axis.ticks.y=element_blank(),
               axis.title.y = element_blank(),
               axis.title.x = element_blank(),
               axis.text=element_text(size=10)
  )


# body size dimorphism
df_bodysizedimorphism <- combined[ complete.cases(combined$strictfdom,combined$SexualDimorphism_MaleWeight_over_FemaleWeight),]
df_bodysizedimorphism[df_bodysizedimorphism$strictfdom==1,]$strictfdom<-"c) Strict male dominance"
df_bodysizedimorphism[df_bodysizedimorphism$strictfdom==2,]$strictfdom<-"b) Co dominance"
df_bodysizedimorphism[df_bodysizedimorphism$strictfdom==3,]$strictfdom<-"a) Strict female dominance"

df = data.frame(
  strictfdom=df_bodysizedimorphism$strictfdom,
  SexualDimorphism_MaleWeight_over_FemaleWeight=exp(df_bodysizedimorphism$SexualDimorphism_MaleWeight_over_FemaleWeight)
)

plot_bodysizedimorphism<-ggplot(df)+aes(y=strictfdom,x=SexualDimorphism_MaleWeight_over_FemaleWeight,fill=strictfdom)+stat_halfeye(aes(thickness = after_stat(pdf*n)))+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text=element_text(size=10)
  )+
  scale_y_discrete(labels = c('Strict female dominance', 'Co dominance', 'Strict male dominance'))+
  scale_x_continuous(name="Male body size relative to female body size")+scale_fill_manual(values=c(female_dominance_color,co_dominance_color,male_dominance_color))+theme(legend.position="none")




# Canine size dimorphism
df_sexratiocanines <- combined[ complete.cases(combined$strictfdom,combined$CanineDimorphism),]
df_sexratiocanines$CanineDimorphism<-df_sexratiocanines$CanineDimorphism
df_sexratiocanines[df_sexratiocanines$strictfdom==1,]$strictfdom<-"c) Strict male dominance"
df_sexratiocanines[df_sexratiocanines$strictfdom==2,]$strictfdom<-"b) Co dominance"
df_sexratiocanines[df_sexratiocanines$strictfdom==3,]$strictfdom<-"a) Strict female dominance"

df = data.frame(
  strictfdom=df_sexratiocanines$strictfdom,
  CanineDimorphism=exp(df_sexratiocanines$CanineDimorphism)
)

plot_caninesizedimorphism<-ggplot(df)+aes(y=strictfdom,x=CanineDimorphism,fill=strictfdom)+stat_halfeye(aes(thickness = after_stat(pdf*n)))+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text=element_text(size=10)
  )+
  scale_x_continuous(name="Male canine size relative to female canine size")+scale_fill_manual(values=c(female_dominance_color,co_dominance_color,male_dominance_color))+theme(legend.position="none")


# home range overlap
df_homerange_overlap <- combined[ complete.cases(combined$strictfdom,combined$homerange_overlap),]

df_homerange_overlap[df_homerange_overlap$strictfdom==1,]$strictfdom<-"c) Strict male dominance"
df_homerange_overlap[df_homerange_overlap$strictfdom==2,]$strictfdom<-"b) Co dominance"
df_homerange_overlap[df_homerange_overlap$strictfdom==3,]$strictfdom<-"a) Strict female dominance"

df = data.frame(
  strictfdom=df_homerange_overlap$strictfdom,
  homerange_overlap=df_homerange_overlap$homerange_overlap
)

plot_homerange_overlap<-ggplot(df)+aes(y=strictfdom,x=homerange_overlap,fill=strictfdom)+stat_halfeye(aes(thickness = after_stat(pdf*n)))+
  theme(        axis.text.y=element_blank(),
                axis.ticks.y=element_blank(),
                axis.title.y = element_blank(),
                axis.title.x = element_blank(),
                axis.text=element_text(size=10)
  )+
  scale_y_discrete(labels = c('Strict female dominance', 'Co dominance', 'Strict male dominance'))+scale_fill_manual(values=c(female_dominance_color,co_dominance_color,male_dominance_color),)+theme(legend.position="none")



# number of females
df_females <- combined[ complete.cases(combined$strictfdom,combined$females),]
df_females[df_females$strictfdom==1,]$strictfdom<-"c) Strict male dominance"
df_females[df_females$strictfdom==2,]$strictfdom<-"b) Co dominance"
df_females[df_females$strictfdom==3,]$strictfdom<-"a) Strict female dominance"
df = data.frame(
  strictfdom=df_females$strictfdom,
  females=df_females$females
)

plot_females<-ggplot(df)+aes(y=strictfdom,x=females,fill=strictfdom)+stat_halfeye(aes(thickness = after_stat(pdf*n)))+
  theme(        axis.text.y=element_blank(),
                axis.ticks.y=element_blank(),
                axis.title.y = element_blank(),
                axis.title.x = element_blank(),
                axis.text=element_text(size=10)
  )+
  scale_y_discrete(labels = c('Strict female dominance', 'Co dominance', 'Strict male dominance'))+scale_fill_manual(values=c(female_dominance_color,co_dominance_color,male_dominance_color),)+theme(legend.position="none")



#  Fig 2 combined plot 
pdf("figures/R_Fig2top.pdf",width=18,height=4)
plot_grid(plot_socialsystem,plot_fissionfusion,plot_bodysizedimorphism,plot_caninesizedimorphism,plot_homerange_overlap,plot_females, rel_widths = c(4,3,3,3,3,3),nrow=1,scale=0.9)
dev.off()


# Fig 2 bottom
# Female eviction

summarizedtable<-combined %>%
  group_by(female_evictions,strictfdom) %>%
  summarize(Total = n())

summarizedtable<-as.data.frame(summarizedtable)

summarizedtable<-summarizedtable[c(1:5),]
summarizedtable[6,]<-c("No",3,0)
colnames(summarizedtable)<-c("female_evictions","StrictFemdom","Observations")
summarizedtable$Observations<-as.integer(summarizedtable$Observations)
summarizedtable[summarizedtable$female_evictions =="Yes", ]$female_evictions<-"present"
summarizedtable[summarizedtable$female_evictions =="No", ]$female_evictions<-"absent"
summarizedtable[summarizedtable$StrictFemdom ==1, ]$StrictFemdom<-"3) male dominance"
summarizedtable[summarizedtable$StrictFemdom ==2, ]$StrictFemdom<-"2) co-dominance"
summarizedtable[summarizedtable$StrictFemdom ==3, ]$StrictFemdom<-"1) female dominance"
summarizedtable<-summarizedtable[order(summarizedtable$female_evictions,summarizedtable$StrictFemdom),]
summarizedtable$StrictFemdom<-as.factor(summarizedtable$StrictFemdom)

plot_femaleevictions <-ggplot(summarizedtable, aes(x = factor(female_evictions,levels=c("absent","present")), y = Observations, fill = StrictFemdom)) + 
  geom_bar(stat = "identity",fill=c(female_dominance_color,co_dominance_color,male_dominance_color,female_dominance_color,co_dominance_color,male_dominance_color))+
  theme(       axis.text.y=element_blank(),
               axis.ticks.y=element_blank(),
               axis.title.y = element_blank(),
               axis.title.x = element_blank(),
               axis.text=element_text(size=10)
  )


# Captivity

summarizedtable<-combined %>%
  group_by(origin,strictfdom) %>%
  summarize(Total = n())

summarizedtable<-as.data.frame(summarizedtable)

colnames(summarizedtable)<-c("origin","StrictFemdom","Observations")
summarizedtable$Observations<-as.integer(summarizedtable$Observations)
summarizedtable[summarizedtable$StrictFemdom ==1, ]$StrictFemdom<-"3) male dominance"
summarizedtable[summarizedtable$StrictFemdom ==2, ]$StrictFemdom<-"2) co-dominance"
summarizedtable[summarizedtable$StrictFemdom ==3, ]$StrictFemdom<-"1) female dominance"
summarizedtable<-summarizedtable[order(summarizedtable$origin,summarizedtable$StrictFemdom),]
summarizedtable$StrictFemdom<-as.factor(summarizedtable$StrictFemdom)

plot_origin <-ggplot(summarizedtable, aes(x = factor(origin,levels=c("wild","captive")), y = Observations, fill = StrictFemdom)) + 
  geom_bar(stat = "identity",fill=c(col.alpha(female_dominance_color,0.4),col.alpha(co_dominance_color,0.4),col.alpha(male_dominance_color,0.4),col.alpha(female_dominance_color,0.4),col.alpha(co_dominance_color,0.4),col.alpha(male_dominance_color,0.4) ))+
  theme(       axis.text.y=element_blank(),
               axis.ticks.y=element_blank(),
               axis.title.y = element_blank(),
               axis.title.x = element_blank(),
               axis.text=element_text(size=10)
  )

# environmental harshness
df_envharshness <- combined[ complete.cases(combined$strictfdom,combined$env_harshness),]

df_envharshness[df_envharshness$strictfdom==1,]$strictfdom<-"c) Strict male dominance"
df_envharshness[df_envharshness$strictfdom==2,]$strictfdom<-"b) Co dominance"
df_envharshness[df_envharshness$strictfdom==3,]$strictfdom<-"a) Strict female dominance"

df = data.frame(
  strictfdom=df_envharshness$strictfdom,
  env_harshness=(df_envharshness$env_harshness)
)

plot_env_harshness<-ggplot(df)+aes(y=strictfdom,x=env_harshness,fill=strictfdom)+stat_halfeye(aes(thickness = after_stat(pdf*n)))+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text=element_text(size=10)
  )+
  scale_y_discrete(labels = c('Strict female dominance', 'Co dominance', 'Strict male dominance'))+
  scale_x_continuous(name="sexual receptivity (days)")+scale_fill_manual(values=c(col.alpha(female_dominance_color,0.4),col.alpha(co_dominance_color,0.4),col.alpha(male_dominance_color,0.4)))+theme(legend.position="none")




# rainfall seasonality
df_rainfallseasonality <- combined[ complete.cases(combined$strictfdom,combined$rainfall_annualvariation),]

df_rainfallseasonality[df_rainfallseasonality$strictfdom==1,]$strictfdom<-"c) Strict male dominance"
df_rainfallseasonality[df_rainfallseasonality$strictfdom==2,]$strictfdom<-"b) Co dominance"
df_rainfallseasonality[df_rainfallseasonality$strictfdom==3,]$strictfdom<-"a) Strict female dominance"

df = data.frame(
  strictfdom=df_rainfallseasonality$strictfdom,
  rainfall_annualvariation=(df_rainfallseasonality$rainfall_annualvariation)
)

plot_rainfall_seasonality<-ggplot(df)+aes(y=strictfdom,x=rainfall_annualvariation,fill=strictfdom)+stat_halfeye(aes(thickness = after_stat(pdf*n)))+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text=element_text(size=10)
  )+
  scale_y_discrete(labels = c('Strict female dominance', 'Co dominance', 'Strict male dominance'))+
  scale_x_continuous(name="rainfall seasonality")+scale_fill_manual(values=c(col.alpha(female_dominance_color,1),col.alpha(co_dominance_color,1),col.alpha(male_dominance_color,1)))+theme(legend.position="none")



# rainfall unpredictability
df_rainfallunpredictability <- combined[ complete.cases(combined$strictfdom,combined$rainfall_unpredictability),]

df_rainfallunpredictability[df_rainfallunpredictability$strictfdom==1,]$strictfdom<-"c) Strict male dominance"
df_rainfallunpredictability[df_rainfallunpredictability$strictfdom==2,]$strictfdom<-"b) Co dominance"
df_rainfallunpredictability[df_rainfallunpredictability$strictfdom==3,]$strictfdom<-"a) Strict female dominance"

df = data.frame(
  strictfdom=df_rainfallunpredictability$strictfdom,
  rainfall_unpredictability=(df_rainfallunpredictability$rainfall_unpredictability)
)

plot_rainfall_unpredictability<-ggplot(df)+aes(y=strictfdom,x=rainfall_unpredictability,fill=strictfdom)+stat_halfeye(aes(thickness = after_stat(pdf*n)))+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text=element_text(size=10)
  )+
  scale_y_discrete(labels = c('Strict female dominance', 'Co dominance', 'Strict male dominance'))+
  scale_x_continuous(name="rainfall unpredictability")+scale_fill_manual(values=c(col.alpha(female_dominance_color,0.4),col.alpha(co_dominance_color,0.4),col.alpha(male_dominance_color,0.4)))+theme(legend.position="none")



# Female infanticide

summarizedtable<-combined %>%
  group_by(female_infanticide,strictfdom) %>%
  summarize(Total = n())

summarizedtable<-as.data.frame(summarizedtable)
summarizedtable<-summarizedtable[c(1:6),]
colnames(summarizedtable)<-c("female_infanticide","StrictFemdom","Observations")
summarizedtable$Observations<-as.integer(summarizedtable$Observations)
summarizedtable[summarizedtable$female_infanticide =="Yes", ]$female_infanticide<-"present"
summarizedtable[summarizedtable$female_infanticide =="No", ]$female_infanticide<-"absent"
summarizedtable[summarizedtable$StrictFemdom ==1, ]$StrictFemdom<-"3) male dominance"
summarizedtable[summarizedtable$StrictFemdom ==2, ]$StrictFemdom<-"2) co-dominance"
summarizedtable[summarizedtable$StrictFemdom ==3, ]$StrictFemdom<-"1) female dominance"
summarizedtable<-summarizedtable[order(summarizedtable$female_infanticide,summarizedtable$StrictFemdom),]
summarizedtable$StrictFemdom<-as.factor(summarizedtable$StrictFemdom)

plot_femaleinfanticide <-ggplot(summarizedtable, aes(x = factor(female_infanticide,levels=c("absent","present")), y = Observations, fill = StrictFemdom)) + 
  geom_bar(stat = "identity",fill=c(col.alpha(female_dominance_color,0.5),col.alpha(co_dominance_color,0.4),col.alpha(male_dominance_color,0.4),col.alpha(female_dominance_color,0.4),col.alpha(co_dominance_color,0.4),col.alpha(male_dominance_color,0.4)))+
  theme(       axis.text.y=element_blank(),
               axis.ticks.y=element_blank(),
               axis.title.y = element_blank(),
               axis.title.x = element_blank(),
               axis.text=element_text(size=10)
  )


# relative canine size
df_caninesize <- combined[ complete.cases(combined$strictfdom,combined$relative_femalecaninesize),]
df_caninesize[df_caninesize$strictfdom==1,]$strictfdom<-"c) Strict male dominance"
df_caninesize[df_caninesize$strictfdom==2,]$strictfdom<-"b) Co dominance"
df_caninesize[df_caninesize$strictfdom==3,]$strictfdom<-"a) Strict female dominance"
df = data.frame(
  strictfdom=df_caninesize$strictfdom,
  relative_femalecaninesize=df_caninesize$relative_femalecaninesize
)

plot_caninesize<-ggplot(df)+aes(y=strictfdom,x=relative_femalecaninesize,fill=strictfdom)+stat_halfeye(aes(thickness = after_stat(pdf*n)))+
  theme(        axis.text.y=element_blank(),
                axis.ticks.y=element_blank(),
                axis.title.y = element_blank(),
                axis.title.x = element_blank(),
                axis.text=element_text(size=10)
  )+
  scale_y_discrete(labels = c('Strict female dominance', 'Co dominance', 'Strict male dominance'))+scale_fill_manual(values=c(col.alpha(female_dominance_color,0.4),col.alpha(co_dominance_color,0.4),col.alpha(male_dominance_color,0.4)))+theme(legend.position="none")




pdf("figures/R_Fig2bottom.pdf",width=18,height=4)
plot_grid(plot_femaleevictions,plot_origin,plot_env_harshness,plot_rainfall_seasonality,plot_rainfall_unpredictability,plot_femaleinfanticide,plot_caninesize, rel_widths = c(3,3,3,3,3,3,3),nrow=1,scale=0.9)
dev.off()





################################################################################
##### Figure 3: sexual selection

# mating system
mostly_females_phylodata_1<-combined
mostly_females_phylodata_1$strictfdom<-as.numeric(as.factor(combined$strictfdom))
mostly_females_phylodata_1<-as.data.frame(mostly_females_phylodata_1 %>% group_by(corrected_species_id) %>% summarise(mean(strictfdom,na.rm=T )) )
mostly_females_phylodata_2<-as.data.frame(mostly_females_phylodata_1[,2])
row.names(mostly_females_phylodata_2)<-mostly_females_phylodata_1[,1]
colnames(mostly_females_phylodata_2)<-"strictfdom"
mostly_females_phylodata_3<-as.data.frame(mostly_females_phylodata_2[is.na(mostly_females_phylodata_2$strictfdom)==FALSE,])
row.names(mostly_females_phylodata_3)<-row.names(mostly_females_phylodata_2)[is.na(mostly_females_phylodata_2$strictfdom)==FALSE]
colnames(mostly_females_phylodata_3)<-"strictfdom"

monogamy<-combined
monogamy$MatSysPMK<-as.numeric(as.factor(monogamy$MatSysPMK))
monogamy<-monogamy[is.na(monogamy$MatSysPMK)==F,]
monogamy_phylodata<-as.data.frame(monogamy %>% group_by(corrected_species_id) %>% summarise(mean(MatSysPMK,na.rm=T )) )
mostly_females_phylodata<-left_join(mostly_females_phylodata_1,monogamy_phylodata,by="corrected_species_id")
rownames(mostly_females_phylodata)<-mostly_females_phylodata$corrected_species_id
colnames(mostly_females_phylodata)<-c("Species","StrictFemdom","MatingSystem")
mostly_females_phylodata<-mostly_females_phylodata[is.na(mostly_females_phylodata$MatingSystem)==F,]
mostly_females_phylodata<-mostly_females_phylodata[is.na(mostly_females_phylodata$StrictFemdom)==F,]
mostly_females_phylodata$StrictFemdom<-round(mostly_females_phylodata$StrictFemdom,0)

summarizedtable<-mostly_females_phylodata %>%
  group_by(MatingSystem,StrictFemdom) %>%
  summarize(NumberOfSpecies = n())

summarizedtable<-combined %>%
  group_by(MatSysPMK,strictfdom) %>%
  summarize(Total = n())

summarizedtable<-as.data.frame(summarizedtable)

summarizedtable<-summarizedtable[1:11,]

summarizedtable<-rbind(summarizedtable,c("MON",1,0))
colnames(summarizedtable)<-c("MatingSystem","StrictFemdom","Observations")
summarizedtable$Observations<-as.integer(summarizedtable$Observations)
summarizedtable[summarizedtable$MatingSystem =="MON", ]$MatingSystem<-"Monogamy"
summarizedtable[summarizedtable$MatingSystem =="PAN", ]$MatingSystem<-"Polyandry"
summarizedtable[summarizedtable$MatingSystem =="POL", ]$MatingSystem<-"Polygyny"
summarizedtable[summarizedtable$MatingSystem =="PRO", ]$MatingSystem<-"Promiscuity"

summarizedtable[summarizedtable$StrictFemdom ==1, ]$StrictFemdom<-"3) male dominance"
summarizedtable[summarizedtable$StrictFemdom ==2, ]$StrictFemdom<-"2) co-dominance"
summarizedtable[summarizedtable$StrictFemdom ==3, ]$StrictFemdom<-"1) female dominance"

summarizedtable<-summarizedtable[order(summarizedtable$MatingSystem,summarizedtable$StrictFemdom),]

summarizedtable$StrictFemdom<-as.factor(summarizedtable$StrictFemdom)

plot_matingsystem <-ggplot(summarizedtable, aes(x = factor(MatingSystem,levels=c("Polygyny","Promiscuity","Polyandry","Monogamy")), y = Observations, fill = StrictFemdom)) + 
  geom_bar(stat = "identity",fill=c(female_dominance_color,co_dominance_color,male_dominance_color,female_dominance_color,co_dominance_color,male_dominance_color,female_dominance_color,co_dominance_color,male_dominance_color,female_dominance_color,co_dominance_color,male_dominance_color))+
  theme(       axis.text.y=element_blank(),
               axis.ticks.y=element_blank(),
               axis.title.y = element_blank(),
               axis.title.x = element_blank(),
               axis.text=element_text(size=10)
  )


# arboreality

summarizedtable<-combined %>%
  group_by(Strata_Wilman,strictfdom) %>%
  summarize(Total = n())

summarizedtable<-as.data.frame(summarizedtable)

summarizedtable[9,]<-c("G",3,0)
summarizedtable<-summarizedtable[1:9,]


colnames(summarizedtable)<-c("Arboreality","StrictFemdom","Observations")
summarizedtable$Observations<-as.integer(summarizedtable$Observations)
summarizedtable[summarizedtable$Arboreality =="Ar", ]$Arboreality<-"Arboreal"
summarizedtable[summarizedtable$Arboreality =="G", ]$Arboreality<-"Ground"
summarizedtable[summarizedtable$Arboreality =="S", ]$Arboreality<-"Scansorial"


summarizedtable[summarizedtable$StrictFemdom ==1, ]$StrictFemdom<-"3) male dominance"
summarizedtable[summarizedtable$StrictFemdom ==2, ]$StrictFemdom<-"2) co-dominance"
summarizedtable[summarizedtable$StrictFemdom ==3, ]$StrictFemdom<-"1) female dominance"

summarizedtable<-summarizedtable[order(summarizedtable$Arboreality,summarizedtable$StrictFemdom),]

summarizedtable$StrictFemdom<-as.factor(summarizedtable$StrictFemdom)

plot_arboreality <-ggplot(summarizedtable, aes(x = factor(Arboreality,levels=c("Ground","Scansorial","Arboreal")), y = Observations, fill = StrictFemdom)) + 
  geom_bar(stat = "identity",fill=c(female_dominance_color,co_dominance_color,male_dominance_color,female_dominance_color,co_dominance_color,male_dominance_color,female_dominance_color,co_dominance_color,male_dominance_color))+
  theme(       axis.text.y=element_blank(),
               axis.ticks.y=element_blank(),
               axis.title.y = element_blank(),
               axis.title.x = element_blank(),
               axis.text=element_text(size=10)
  )


# sexual receptivity
df_sexualreceptivity <- combined[ complete.cases(combined$strictfdom,combined$sexualreceptivity_hours),]

df_sexualreceptivity[df_sexualreceptivity$strictfdom==1,]$strictfdom<-"c) Strict male dominance"
df_sexualreceptivity[df_sexualreceptivity$strictfdom==2,]$strictfdom<-"b) Co dominance"
df_sexualreceptivity[df_sexualreceptivity$strictfdom==3,]$strictfdom<-"a) Strict female dominance"

df = data.frame(
  strictfdom=df_sexualreceptivity$strictfdom,
  sexualreceptivity_hours=exp(df_sexualreceptivity$sexualreceptivity_hours)/24
)

plot_sexualreceptivity<-ggplot(df)+aes(y=strictfdom,x=sexualreceptivity_hours,fill=strictfdom)+stat_halfeye(aes(thickness = after_stat(pdf*n)))+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text=element_text(size=10)
  )+
  scale_y_discrete(labels = c('Strict female dominance', 'Co dominance', 'Strict male dominance'))+
  scale_x_continuous(name="sexual receptivity (days)")+scale_fill_manual(values=c(female_dominance_color,co_dominance_color,male_dominance_color))+theme(legend.position="none")



# sex ratio
df_sexratio <- combined[ complete.cases(combined$strictfdom,combined$sexratio),]

df_sexratio[df_sexratio$strictfdom==1,]$strictfdom<-"c) Strict male dominance"
df_sexratio[df_sexratio$strictfdom==2,]$strictfdom<-"b) Co dominance"
df_sexratio[df_sexratio$strictfdom==3,]$strictfdom<-"a) Strict female dominance"

df = data.frame(
  strictfdom=df_sexratio$strictfdom,
  sexratio_malesperfemale=df_sexratio$sexratio
)

plot_sexratio<-ggplot(df)+aes(y=strictfdom,x=sexratio_malesperfemale,fill=strictfdom)+stat_halfeye(aes(thickness = after_stat(pdf*n)))+
  theme(        axis.text.y=element_blank(),
                axis.ticks.y=element_blank(),
                axis.title.y = element_blank(),
                axis.title.x = element_blank(),
                axis.text=element_text(size=10)
  )+
  scale_y_discrete(labels = c('Strict female dominance', 'Co dominance', 'Strict male dominance'))+scale_x_continuous(name="Sex ratio in social groups",lim=c(0.09,0.76),breaks=c(0.1, 0.33,0.5, 0.66),labels=c("1♂ / 10♀♀","1♂ / 2♀♀","1♂ / 1♀","2♂♂ / 1♀"))+scale_fill_manual(values=c(female_dominance_color,co_dominance_color,male_dominance_color),)+theme(legend.position="none")



# reproductive skew
df_skew <- combined[ complete.cases(combined$strictfdom,combined$M_skew_index),]

df_skew[df_skew$strictfdom==1,]$strictfdom<-"c) Strict male dominance"
df_skew[df_skew$strictfdom==2,]$strictfdom<-"b) Co dominance"
df_skew[df_skew$strictfdom==3,]$strictfdom<-"a) Strict female dominance"

df = data.frame(
  strictfdom=df_skew$strictfdom,
  skew=df_skew$M_skew_index
)

plot_reproductiveskew<-ggplot(df)+aes(y=strictfdom,x=skew,fill=strictfdom)+stat_halfeye(aes(thickness = after_stat(pdf*n)))+
  theme(       axis.text.y=element_blank(),
               axis.ticks.y=element_blank(),
               axis.title.y = element_blank(),
               axis.title.x = element_blank(),
               axis.text=element_text(size=10)
  )+
  scale_y_discrete(labels = c('Strict female dominance', 'Co dominance', 'Strict male dominance'))+scale_x_continuous(name="Male reproductive skew less or more than expected by chance")+scale_fill_manual(values=c(female_dominance_color,co_dominance_color,male_dominance_color),)+theme(legend.position="none")+geom_point(x=100,y=1.1,pch=24, fill=female_dominance_color, alpha=0.5,size=8, colour=female_dominance_color)+geom_point(x=100,y=1,pch=21, fill="black", alpha=0.5,size=3, colour="black") 




# reproductive synchrony
df_synchrony <- combined[ complete.cases(combined$strictfdom,combined$Synchrony),]

df_synchrony[df_synchrony$strictfdom==1,]$strictfdom<-"c) Strict male dominance"
df_synchrony[df_synchrony$strictfdom==2,]$strictfdom<-"b) Co dominance"
df_synchrony[df_synchrony$strictfdom==3,]$strictfdom<-"a) Strict female dominance"

df = data.frame(
  strictfdom=df_synchrony$strictfdom,
  Synchrony=df_synchrony$Synchrony
)

plot_synchrony<-ggplot(df)+aes(y=strictfdom,x=Synchrony,fill=strictfdom)+stat_halfeye(aes(thickness = after_stat(pdf*n)))+
  theme(       axis.text.y=element_blank(),
               axis.ticks.y=element_blank(),
               axis.title.y = element_blank(),
               axis.title.x = element_blank(),
               axis.text=element_text(size=10)
  )+
  scale_y_discrete(labels = c('Strict female dominance', 'Co dominance', 'Strict male dominance'))+scale_x_continuous(name="Reproductive synchrony")+scale_fill_manual(values=c( col.alpha(female_dominance_color,0.4),col.alpha(co_dominance_color,0.4),col.alpha(male_dominance_color,0.4) ),)+theme(legend.position="none")



# testes size
df_testessize <- combined[ complete.cases(combined$strictfdom,combined$relative_testes_mass),]

df_testessize[df_testessize$strictfdom==1,]$strictfdom<-"c) Strict male dominance"
df_testessize[df_testessize$strictfdom==2,]$strictfdom<-"b) Co dominance"
df_testessize[df_testessize$strictfdom==3,]$strictfdom<-"a) Strict female dominance"

df = data.frame(
  strictfdom=df_testessize$strictfdom,
  relative_testes_mass=df_testessize$relative_testes_mass
)

plot_testesmass<-ggplot(df)+aes(y=strictfdom,x=relative_testes_mass,fill=strictfdom)+stat_halfeye(aes(thickness = after_stat(pdf*n)))+
  theme(        axis.text.y=element_blank(),
                axis.ticks.y=element_blank(),
                axis.title.y = element_blank(),
                axis.title.x = element_blank(),
                axis.text=element_text(size=10)
  )+
  scale_y_discrete(labels = c('Strict female dominance', 'Co dominance', 'Strict male dominance'))+scale_fill_manual(values=c(col.alpha(female_dominance_color,0.4),col.alpha(co_dominance_color,0.4),col.alpha(male_dominance_color,0.4)))+theme(legend.position="none")




# combined plot
pdf("figures/R_Fig3.pdf",width=18,height=4)
plot_grid(plot_matingsystem,plot_arboreality,plot_sexualreceptivity, plot_sexratio,plot_reproductiveskew,plot_synchrony,plot_testesmass, rel_widths = c(5,3,3,3,3,3,3),nrow=1)
dev.off()



################################################################################
##### Figure 4 offspring loss


# Park
summarizedtable<-combined %>%
  group_by(Park,strictfdom) %>%
  summarize(Total = n())

summarizedtable<-as.data.frame(summarizedtable)
summarizedtable<-summarizedtable[1:5,]
summarizedtable[6,]<-c("Yes",1,0)
colnames(summarizedtable)<-c("Park","StrictFemdom","Observations")
summarizedtable$Observations<-as.integer(summarizedtable$Observations)
summarizedtable[summarizedtable$Park =="Yes", ]$Park<-"parked"
summarizedtable[summarizedtable$Park =="No", ]$Park<-"carried"
summarizedtable[summarizedtable$StrictFemdom ==1, ]$StrictFemdom<-"3) male dominance"
summarizedtable[summarizedtable$StrictFemdom ==2, ]$StrictFemdom<-"2) co-dominance"
summarizedtable[summarizedtable$StrictFemdom ==3, ]$StrictFemdom<-"1) female dominance"
summarizedtable<-summarizedtable[order(summarizedtable$Park,summarizedtable$StrictFemdom),]
summarizedtable$StrictFemdom<-as.factor(summarizedtable$StrictFemdom)

plot_park <-ggplot(summarizedtable, aes(x = factor(Park,levels=c("parked","carried")), y = Observations, fill = StrictFemdom)) + 
  geom_bar(stat = "identity",fill=c(female_dominance_color,co_dominance_color,male_dominance_color,female_dominance_color,co_dominance_color,male_dominance_color))+
  theme(       axis.text.y=element_blank(),
               axis.ticks.y=element_blank(),
               axis.title.y = element_blank(),
               axis.title.x = element_blank(),
               axis.text=element_text(size=10)
  )



# relative lactation duration
df_RelativeLactationDuration <- combined[ complete.cases(combined$strictfdom,combined$RelativeLactationDuration),]

df_RelativeLactationDuration[df_RelativeLactationDuration$strictfdom==1,]$strictfdom<-"c) Strict male dominance"
df_RelativeLactationDuration[df_RelativeLactationDuration$strictfdom==2,]$strictfdom<-"b) Co dominance"
df_RelativeLactationDuration[df_RelativeLactationDuration$strictfdom==3,]$strictfdom<-"a) Strict female dominance"

df = data.frame(
  strictfdom=df_RelativeLactationDuration$strictfdom,
  RelativeLactationDuration=df_RelativeLactationDuration$RelativeLactationDuration
)

plot_RelativeLactationDuration<-ggplot(df[1:153,])+aes(y=strictfdom,x=RelativeLactationDuration,fill=strictfdom)+stat_halfeye(aes(thickness = after_stat(pdf*n)))+
  theme(        axis.text.y=element_blank(),
                axis.ticks.y=element_blank(),
                axis.title.y = element_blank(),
                axis.title.x = element_blank(),
                axis.text=element_text(size=10)
  )+
  scale_y_discrete(labels = c('Strict female dominance', 'Co dominance', 'Strict male dominance'))+scale_fill_manual(values=c(col.alpha(female_dominance_color,1),col.alpha(co_dominance_color,1),col.alpha(male_dominance_color,1) ))+theme(legend.position="none")



# infanticide

summarizedtable<-combined %>%
  group_by(maleinfanticide,strictfdom) %>%
  summarize(Total = n())
summarizedtable<-as.data.frame(summarizedtable)
summarizedtable<-summarizedtable[1:6,]
colnames(summarizedtable)<-c("maleinfanticide","StrictFemdom","Observations")
summarizedtable$Observations<-as.integer(summarizedtable$Observations)
summarizedtable[summarizedtable$maleinfanticide =="Yes", ]$maleinfanticide<-"present"
summarizedtable[summarizedtable$maleinfanticide =="No", ]$maleinfanticide<-"absent"
summarizedtable[summarizedtable$StrictFemdom ==1, ]$StrictFemdom<-"3) male dominance"
summarizedtable[summarizedtable$StrictFemdom ==2, ]$StrictFemdom<-"2) co-dominance"
summarizedtable[summarizedtable$StrictFemdom ==3, ]$StrictFemdom<-"1) female dominance"
summarizedtable<-summarizedtable[order(summarizedtable$maleinfanticide,summarizedtable$StrictFemdom),]
summarizedtable$StrictFemdom<-as.factor(summarizedtable$StrictFemdom)

plot_maleinfanticide <-ggplot(summarizedtable, aes(x = factor(maleinfanticide,levels=c("absent","present")), y = Observations, fill = StrictFemdom)) + 
  geom_bar(stat = "identity",fill=c(col.alpha(female_dominance_color,0.4),col.alpha(co_dominance_color,0.4),col.alpha(male_dominance_color,0.4),col.alpha(female_dominance_color,0.4),col.alpha(co_dominance_color,0.4),col.alpha(male_dominance_color,0.4)))+
  theme(       axis.text.y=element_blank(),
               axis.ticks.y=element_blank(),
               axis.title.y = element_blank(),
               axis.title.x = element_blank(),
               axis.text=element_text(size=10)
  )


# allomaternal care

summarizedtable<-combined %>%
  group_by(AlloMaternalCare,strictfdom) %>%
  summarize(Total = n())
summarizedtable<-as.data.frame(summarizedtable)
summarizedtable<-summarizedtable[1:6,]
colnames(summarizedtable)<-c("AlloMaternalCare","StrictFemdom","Observations")
summarizedtable$Observations<-as.integer(summarizedtable$Observations)
summarizedtable[summarizedtable$AlloMaternalCare =="Yes", ]$AlloMaternalCare<-"present"
summarizedtable[summarizedtable$AlloMaternalCare =="No", ]$AlloMaternalCare<-"absent"
summarizedtable[summarizedtable$StrictFemdom ==1, ]$StrictFemdom<-"3) male dominance"
summarizedtable[summarizedtable$StrictFemdom ==2, ]$StrictFemdom<-"2) co-dominance"
summarizedtable[summarizedtable$StrictFemdom ==3, ]$StrictFemdom<-"1) female dominance"
summarizedtable<-summarizedtable[order(summarizedtable$AlloMaternalCare,summarizedtable$StrictFemdom),]
summarizedtable$StrictFemdom<-as.factor(summarizedtable$StrictFemdom)

plot_AlloMaternalCare <-ggplot(summarizedtable, aes(x = factor(AlloMaternalCare,levels=c("absent","present")), y = Observations, fill = StrictFemdom)) + 
  geom_bar(stat = "identity",fill=c(col.alpha(female_dominance_color,0.4),col.alpha(co_dominance_color,0.4),col.alpha(male_dominance_color,0.4),col.alpha(female_dominance_color,0.4),col.alpha(co_dominance_color,0.4),col.alpha(male_dominance_color,0.4)))+
  theme(       axis.text.y=element_blank(),
               axis.ticks.y=element_blank(),
               axis.title.y = element_blank(),
               axis.title.x = element_blank(),
               axis.text=element_text(size=10)
  )

# allomaternal care continuous

df_AlloMaternalCare_continuous <- combined[ complete.cases(combined$strictfdom,combined$AlloMaternalCare_continuous),]

df_AlloMaternalCare_continuous[df_AlloMaternalCare_continuous$strictfdom==1,]$strictfdom<-"c) Strict male dominance"
df_AlloMaternalCare_continuous[df_AlloMaternalCare_continuous$strictfdom==2,]$strictfdom<-"b) Co dominance"
df_AlloMaternalCare_continuous[df_AlloMaternalCare_continuous$strictfdom==3,]$strictfdom<-"a) Strict female dominance"

df = data.frame(
  strictfdom=df_AlloMaternalCare_continuous$strictfdom,
  AlloMaternalCare_continuous=df_AlloMaternalCare_continuous$AlloMaternalCare_continuous
)

plot_AlloMaternalCare_continuous<-ggplot(df[1:153,])+aes(y=strictfdom,x=AlloMaternalCare_continuous,fill=strictfdom)+stat_halfeye(aes(thickness = after_stat(pdf*n)))+
  theme(        axis.text.y=element_blank(),
                axis.ticks.y=element_blank(),
                axis.title.y = element_blank(),
                axis.title.x = element_blank(),
                axis.text=element_text(size=10)
  )+
  scale_y_discrete(labels = c('Strict female dominance', 'Co dominance', 'Strict male dominance'))+scale_fill_manual(values=c(col.alpha(female_dominance_color,1),col.alpha(co_dominance_color,1),col.alpha(male_dominance_color,1) ))+theme(legend.position="none")




# combined plot
pdf("figures/R_Fig4_left_new.pdf",width=14.4,height=4)
plot_grid(plot_park, plot_RelativeLactationDuration,plot_maleinfanticide,plot_AlloMaternalCare,plot_AlloMaternalCare_continuous, rel_widths = c(3,4,3,3,4),nrow=1,scale=0.9)
dev.off()








################################################################################
##### Figure 5 female sociality



# sex bias in dispersal

summarizedtable<-combined %>%
  group_by(sexbias_dispersal,strictfdom) %>%
  summarize(Total = n())

summarizedtable<-as.data.frame(summarizedtable)
summarizedtable<-summarizedtable[1:9,]
summarizedtable[9,]<-c("Female",3,0)
colnames(summarizedtable)<-c("sexbias_dispersal","StrictFemdom","Observations")
summarizedtable$Observations<-as.integer(summarizedtable$Observations)
summarizedtable[summarizedtable$sexbias_dispersal =="Female", ]$sexbias_dispersal<-"female"
summarizedtable[summarizedtable$sexbias_dispersal =="Male", ]$sexbias_dispersal<-"male"
summarizedtable[summarizedtable$sexbias_dispersal =="Both", ]$sexbias_dispersal<-"both"
summarizedtable[summarizedtable$StrictFemdom ==1, ]$StrictFemdom<-"3) male dominance"
summarizedtable[summarizedtable$StrictFemdom ==2, ]$StrictFemdom<-"2) co-dominance"
summarizedtable[summarizedtable$StrictFemdom ==3, ]$StrictFemdom<-"1) female dominance"
summarizedtable<-summarizedtable[order(summarizedtable$sexbias_dispersal,summarizedtable$StrictFemdom),]
summarizedtable$StrictFemdom<-as.factor(summarizedtable$StrictFemdom)

plot_sexbias_dispersal <-ggplot(summarizedtable, aes(x = factor(sexbias_dispersal,levels=c("female","both","male")), y = Observations, fill = StrictFemdom)) + 
  geom_bar(stat = "identity",fill=c(female_dominance_color,co_dominance_color,male_dominance_color,female_dominance_color,co_dominance_color,male_dominance_color,female_dominance_color,co_dominance_color,male_dominance_color))+
  theme(       axis.text.y=element_blank(),
               axis.ticks.y=element_blank(),
               axis.title.y = element_blank(),
               axis.title.x = element_blank(),
               axis.text=element_text(size=10)
  )



# female average relatedness
df_female_relatedness <- combined[ complete.cases(combined$strictfdom,combined$female_average_relatedness),]

df_female_relatedness[df_female_relatedness$strictfdom==1,]$strictfdom<-"c) Strict male dominance"
df_female_relatedness[df_female_relatedness$strictfdom==2,]$strictfdom<-"b) Co dominance"
df_female_relatedness[df_female_relatedness$strictfdom==3,]$strictfdom<-"a) Strict female dominance"

df = data.frame(
  strictfdom=df_female_relatedness$strictfdom,
  female_average_relatedness=df_female_relatedness$female_average_relatedness
)

plot_female_relatedness<-ggplot(df[1:47,])+aes(y=strictfdom,x=female_average_relatedness,fill=strictfdom)+stat_halfeye(aes(thickness = after_stat(pdf*n)))+
  theme(        axis.text.y=element_blank(),
                axis.ticks.y=element_blank(),
                axis.title.y = element_blank(),
                axis.title.x = element_blank(),
                axis.text=element_text(size=10)
  )+
  scale_y_discrete(labels = c('Strict female dominance', 'Co dominance', 'Strict male dominance'))+scale_fill_manual(values=c(col.alpha(female_dominance_color,0.4),col.alpha(co_dominance_color,0.4),col.alpha(male_dominance_color,0.4) ))+theme(legend.position="none")



# female coalitions

summarizedtable<-combined %>%
  group_by(jointaggression_females,strictfdom) %>%
  summarize(Total = n())
summarizedtable<-as.data.frame(summarizedtable)
summarizedtable<-summarizedtable[1:6,]
colnames(summarizedtable)<-c("jointaggression_females","StrictFemdom","Observations")
summarizedtable$Observations<-as.integer(summarizedtable$Observations)
summarizedtable[summarizedtable$jointaggression_females =="Yes", ]$jointaggression_females<-"present"
summarizedtable[summarizedtable$jointaggression_females =="No", ]$jointaggression_females<-"absent"
summarizedtable[summarizedtable$StrictFemdom ==1, ]$StrictFemdom<-"3) male dominance"
summarizedtable[summarizedtable$StrictFemdom ==2, ]$StrictFemdom<-"2) co-dominance"
summarizedtable[summarizedtable$StrictFemdom ==3, ]$StrictFemdom<-"1) female dominance"
summarizedtable<-summarizedtable[order(summarizedtable$jointaggression_females,summarizedtable$StrictFemdom),]
summarizedtable$StrictFemdom<-as.factor(summarizedtable$StrictFemdom)

plot_jointaggression_females <-ggplot(summarizedtable, aes(x = factor(jointaggression_females,levels=c("absent","present")), y = Observations, fill = StrictFemdom)) + 
  geom_bar(stat = "identity",fill=c(col.alpha(female_dominance_color,0.4),col.alpha(co_dominance_color,0.4),col.alpha(male_dominance_color,0.4),col.alpha(female_dominance_color,0.4),col.alpha(co_dominance_color,0.4),col.alpha(male_dominance_color,0.4)))+
  theme(       axis.text.y=element_blank(),
               axis.ticks.y=element_blank(),
               axis.title.y = element_blank(),
               axis.title.x = element_blank(),
               axis.text=element_text(size=10)
  )



# combined plot
pdf("figures/R_Fig5_left.pdf",width=14.4,height=4)
plot_grid(plot_sexbias_dispersal, plot_female_relatedness,plot_jointaggression_females, rel_widths = c(4,3,3),nrow=1,scale=0.9)
dev.off()







################################################################################
##### Figure 6 self organisation

# sex ratio
df_sexratio <- combined[ complete.cases(combined$strictfdom,combined$sexratio),]

df_sexratio[df_sexratio$strictfdom==1,]$strictfdom<-"c) Strict male dominance"
df_sexratio[df_sexratio$strictfdom==2,]$strictfdom<-"b) Co dominance"
df_sexratio[df_sexratio$strictfdom==3,]$strictfdom<-"a) Strict female dominance"

df = data.frame(
  strictfdom=df_sexratio$strictfdom,
  sexratio_malesperfemale=df_sexratio$sexratio
)

plot_sexratio<-ggplot(df)+aes(y=strictfdom,x=sexratio_malesperfemale,fill=strictfdom)+stat_halfeye(aes(thickness = after_stat(pdf*n)))+
  theme(        axis.text.y=element_blank(),
                axis.ticks.y=element_blank(),
                axis.title.y = element_blank(),
                axis.title.x = element_blank(),
                axis.text=element_text(size=10)
  )+
  scale_y_discrete(labels = c('Strict female dominance', 'Co dominance', 'Strict male dominance'))+scale_x_continuous(name="Sex ratio in social groups",lim=c(0.09,0.76),breaks=c(0.1, 0.33,0.5, 0.66),labels=c("1♂ / 10♀♀","1♂ / 2♀♀","1♂ / 1♀","2♂♂ / 1♀"))+scale_fill_manual(values=c(female_dominance_color,co_dominance_color,male_dominance_color),)+theme(legend.position="none")




# proportion of male-male conflicts
df_perc_aggress_mm <- combined[ complete.cases(combined$strictfdom,combined$perc_aggression_mm),]

df_perc_aggress_mm[df_perc_aggress_mm$strictfdom==1,]$strictfdom<-"c) Strict male dominance"
df_perc_aggress_mm[df_perc_aggress_mm$strictfdom==2,]$strictfdom<-"b) Co dominance"
df_perc_aggress_mm[df_perc_aggress_mm$strictfdom==3,]$strictfdom<-"a) Strict female dominance"

df = data.frame(
  strictfdom=df_perc_aggress_mm$strictfdom,
  perc_aggression_mm=df_perc_aggress_mm$perc_aggression_mm
)

plot_percaggress_mm<-ggplot(df)+aes(y=strictfdom,x=perc_aggression_mm,fill=strictfdom)+stat_halfeye(aes(thickness = after_stat(pdf*n)))+
  theme(        axis.text.y=element_blank(),
                axis.ticks.y=element_blank(),
                axis.title.y = element_blank(),
                axis.title.x = element_blank(),
                axis.text=element_text(size=10)
  )+
  scale_y_discrete(labels = c('Strict female dominance', 'Co dominance', 'Strict male dominance'))+scale_fill_manual(values=c(col.alpha(female_dominance_color,1),col.alpha(co_dominance_color,1),col.alpha(male_dominance_color,1) ))+theme(legend.position="none")





# males
df_males <- combined[ complete.cases(combined$strictfdom,combined$males),]

df_males[df_males$strictfdom==1,]$strictfdom<-"c) Strict male dominance"
df_males[df_males$strictfdom==2,]$strictfdom<-"b) Co dominance"
df_males[df_males$strictfdom==3,]$strictfdom<-"a) Strict female dominance"

df = data.frame(
  strictfdom=df_males$strictfdom,
  males=exp(df_males$males)
)

plot_males<-ggplot(df)+aes(y=strictfdom,x=males,fill=strictfdom)+stat_halfeye(aes(thickness = after_stat(pdf*n)))+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text=element_text(size=10)
  )+
  scale_y_discrete(labels = c('Strict female dominance', 'Co dominance', 'Strict male dominance'))+
  scale_x_continuous(name="number of males")+scale_fill_manual(values=c(col.alpha(female_dominance_color,0.4),col.alpha(co_dominance_color,0.4),col.alpha(male_dominance_color,0.4) ))+theme(legend.position="none")


# combined plot
pdf("figures/R_Fig5_right.pdf",width=10.8,height=4)
plot_grid(plot_sexratio,plot_percaggress_mm, plot_males, rel_widths = c(3,3,3),nrow=1,scale=0.9)
dev.off()












# 
# ### bottom: model output
# 
# # mating system
# 
# speciesaverage<-combined[is.na(combined$MatSysPMK)==F,]
# 
# dat_list_strict <- list(
#   R = as.integer(as.factor(1/as.integer(speciesaverage$strictfdom))),
#   matingsystem = as.integer(as.factor(speciesaverage$MatSysPMK)),
#   species = as.integer(as.factor(speciesaverage$corrected_species_id))
# )
# 
# m_matingsystems <- ulam(
#   alist(
#     R ~ dordlogit( phi , cutpoints ),
#     phi <-a + b[matingsystem],
#     a ~ normal( 0 , 5 ),
#     b[matingsystem] ~ dnorm(0,5),
#     cutpoints ~ dnorm( 0 , 5 )
#   ) , data=dat_list_strict , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)
# 
# 
# matingsystems <- seq(from=1,to=4,by=1)
# pdat <- data.frame(matingsystems=matingsystems)
# 
# overallphi_speciesaverages<-matrix(ncol=1, nrow=4)
# for(i in 1:length(matingsystems)){
#   overallphi_speciesaverages[i,]<-precis(m_matingsystems,depth=2)[1,1]+precis(m_matingsystems,depth=2)[i+1,1]*matingsystems[i]
# }
# 
# overallprobs_speciesaverages_mating<-pordlogit( 1:3 , overallphi_speciesaverages , precis(m_matingsystems,depth=2)[6:7,1] )
# 
# overallprobs_speciesaverages_mating[,3]<-overallprobs_speciesaverages_mating[,3]-overallprobs_speciesaverages_mating[,2]
# overallprobs_speciesaverages_mating[,2]<-overallprobs_speciesaverages_mating[,2]-overallprobs_speciesaverages_mating[,1]
# overallprobs_speciesaverages_mating <-t(overallprobs_speciesaverages_mating)
# 
# 
# overallprobs_speciesaverages_mating <-overallprobs_speciesaverages_mating[,c(1,2,4,3)]
# colnames(overallprobs_speciesaverages_mating)<-c("Polygynous","Promiscous","Polyandrous", "Monogamous")
# 
# 
# # dimorphism body size
# 
# speciesaverage<-as.data.frame(combined %>% group_by(corrected_species_id) %>% summarise(mean(SexualDimorphism_MaleWeight_over_FemaleWeight,na.rm=T),mean(as.numeric(strictfdom))))
# 
# colnames(speciesaverage)<-c("corrected_species_id","sizedimorphism","strictfdom")
# speciesaverage<-speciesaverage[is.na(speciesaverage$sizedimorphism)==F,]
# speciesaverage$strictfdom<-round(speciesaverage$strictfdom,0)
# 
# dat_list_strict_bodydimorphism <- list(
#   R = as.integer(as.factor(1/speciesaverage$strictfdom)),
#   sizedimorphism = standardize(exp(speciesaverage$sizedimorphism)),
#   species = as.integer(as.factor(speciesaverage$corrected_species_id))
# )
# 
# m_speciesaverage_sizedimorphism <- ulam(
#   alist(
#     R ~ dordlogit( phi , cutpoints ),
#     phi <-a + b*sizedimorphism ,
#     a ~ normal( 0 , 5 ),
#     b ~ dnorm(0,5),
#     cutpoints ~ dnorm( 0 , 5 )
#   ) , data= dat_list_strict_bodydimorphism , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)
# 
# 
# # standardized size dimorphism
# sizedimorphism <- seq(from=-1.5,to=3.5,by=0.5)
# pdat <- data.frame(sizedimorphism=sizedimorphism)
# 
# overallphi_speciesaverages<-matrix(ncol=1, nrow=11)
# for(i in 1:length(sizedimorphism)){
#   overallphi_speciesaverages[i,]<-precis(m_speciesaverage_sizedimorphism,depth=2)[1,1]+precis(m_speciesaverage_sizedimorphism,depth=2)[2,1]*sizedimorphism[i]
# }
# 
# overallprobs_speciesaverages_sizedimorphism<-pordlogit( 1:3 , overallphi_speciesaverages , precis(m_speciesaverage_sizedimorphism,depth=2)[3:4,1] )
# 
# 
# 
# 
# # dimorphism canine size
# 
# speciesaverage<-as.data.frame(combined %>% group_by(corrected_species_id) %>% summarise(mean(CanineDimorphism,na.rm=T),mean(as.numeric(strictfdom))))
# 
# colnames(speciesaverage)<-c("corrected_species_id","sizedimorphism","strictfdom")
# speciesaverage<-speciesaverage[is.na(speciesaverage$sizedimorphism)==F,]
# speciesaverage$strictfdom<-round(speciesaverage$strictfdom,0)
# 
# dat_list_strict_caninedimorphism <- list(
#   R = as.integer(as.factor(1/speciesaverage$strictfdom)),
#   sizedimorphism = standardize(exp(speciesaverage$sizedimorphism)),
#   species = as.integer(as.factor(speciesaverage$corrected_species_id))
# )
# 
# m_speciesaverage_caninesizedimorphism <- ulam(
#   alist(
#     R ~ dordlogit( phi , cutpoints ),
#     phi <-a + b*sizedimorphism ,
#     a ~ normal( 0 , 5 ),
#     b ~ dnorm(0,5),
#     cutpoints ~ dnorm( 0 , 5 )
#   ) , data= dat_list_strict_caninedimorphism , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)
# 
# 
# # standardized size dimorphism
# sizedimorphism_canine <- seq(from=-1.5,to=4.5,by=0.5)
# pdat <- data.frame(sizedimorphism=sizedimorphism)
# 
# overallphi_speciesaverages<-matrix(ncol=1, nrow=13)
# for(i in 1:length(sizedimorphism_canine)){
#   overallphi_speciesaverages[i,]<-precis(m_speciesaverage_caninesizedimorphism,depth=2)[1,1]+precis(m_speciesaverage_caninesizedimorphism,depth=2)[2,1]* sizedimorphism_canine[i]
# }
# 
# overallprobs_speciesaverages_caninedimorphism<-pordlogit( 1:3 , overallphi_speciesaverages , precis(m_speciesaverage_sizedimorphism,depth=2)[3:4,1] )
# 
# 
# 
# # sex ratio 
# 
# speciesaverage<-combined[is.na(combined$sexratio)==F,]
# 
# dat_list_strict_sexratio <- list(
#   R = as.integer(as.factor(1/as.integer(speciesaverage$strictfdom))),
#   sexratio = standardize(speciesaverage$sexratio),
#   species = as.integer(as.factor(speciesaverage$corrected_species_id))
# )
# 
# m_speciesaverage_sexratio <- ulam(
#   alist(
#     R ~ dordlogit( phi , cutpoints ),
#     phi <-a + b*sexratio ,
#     a ~ normal( 0 , 5 ),
#     b ~ dnorm(0,5),
#     cutpoints ~ dnorm( 0 , 5 )
#   ) , data=dat_list_strict_sexratio , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)
# 
# 
# sexratio <- seq(from=-2.5,to=2.5,by=0.5)
# pdat <- data.frame(sexratio=sexratio)
# 
# overallphi_speciesaverages<-matrix(ncol=1, nrow=11)
# for(i in 1:length(sexratio)){
#   overallphi_speciesaverages[i,]<-precis(m_speciesaverage_sexratio,depth=2)[1,1]+precis(m_speciesaverage_sexratio,depth=2)[2,1]*sexratio[i]
# }
# 
# overallprobs_speciesaverages_sexratio<-pordlogit( 1:3 , overallphi_speciesaverages , precis(m_speciesaverage_sexratio,depth=2)[3:4,1] )
# 
# 
# 
# # male reproductive skew
# 
# 
# speciesaverage<-combined[is.na(combined$M_skew_index)==F,]
# 
# dat_list_strict_skew <- list(
#   R = as.integer(as.factor(1/as.integer(speciesaverage$strictfdom))),
#   male_skew = standardize(speciesaverage$M_skew_index),
#   species = as.integer(as.factor(speciesaverage$corrected_species_id))
# )
# 
# m_speciesaverage_skew <- ulam(
#   alist(
#     R ~ dordlogit( phi , cutpoints ),
#     phi <-a + b*male_skew ,
#     a ~ normal( 0 , 5 ),
#     b ~ dnorm(0,5),
#     cutpoints ~ dnorm( 0 , 5 )
#   ) , data=dat_list_strict_skew , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)
# 
# 
# male_skew <- seq(from=-3.0,to=2.0,by=0.5)
# pdat <- data.frame(male_skew=male_skew)
# 
# overallphi_speciesaverages<-matrix(ncol=1, nrow=11)
# for(i in 1:length(male_skew)){
#   overallphi_speciesaverages[i,]<-precis(m_speciesaverage_skew,depth=2)[1,1]+precis(m_speciesaverage_skew,depth=2)[2,1]*male_skew[i]
# }
# 
# overallprobs_speciesaverages_skew<-pordlogit( 1:3 , overallphi_speciesaverages , precis(m_speciesaverage_skew,depth=2)[3:4,1] )
# 
# 
# 
# 
# 
# 
# # combined plot
# 
# pdf("R_Fig2a_bottom.pdf",width=18,height=4)
# 
# previouspar<-par()
# op <- par(oma=c(0.2,0.2,0.2,0.2), mar=c(0.3,1.0,0.3,1.0), mfrow=c(1,5))
# barplot(overallprobs_speciesaverages_mating,col=dominance_colors,axisnames = F)
# 
# plot( NULL , type="n" , xlab="Male body weight relative to female body weight" ,
#       xlim=c(min(dat_list_strict_bodydimorphism $sizedimorphism,na.rm=T),max(dat_list_strict_bodydimorphism $sizedimorphism,na.rm=T)) , ylim=c(0,1) ,  yaxp=c(0,1,4) ,bty="n",xaxt="n",ylab="",yaxt="n")
# axis(side=1,at=c(-1.347572,-0.768047,-0.1885223,0.3910024,0.9705271,1.550052,2.129577,2.709101,3.288626),labels=FALSE)
# polygon(x=c(sizedimorphism,rev(sizedimorphism)),y=c(overallprobs_speciesaverages_sizedimorphism[,2],rev(overallprobs_speciesaverages_sizedimorphism[,1])),col=co_dominance_color,border=NA)
# polygon(x=c(sizedimorphism,rev(sizedimorphism)),y=c(overallprobs_speciesaverages_sizedimorphism[,1],rep(0,11)),col=female_dominance_color,border=NA)
# polygon(x=c(sizedimorphism,rev(sizedimorphism)),y=c(rep(1,11),rev(overallprobs_speciesaverages_sizedimorphism[,2])),col=male_dominance_color,border=NA)
# 
# plot( NULL , type="n" , xlab="Male canine size relative to female canine size" ,
#       xlim=c(min(dat_list_strict_caninedimorphism $sizedimorphism,na.rm=T),max(dat_list_strict_caninedimorphism $sizedimorphism,na.rm=T)) , ylim=c(0,1) ,  yaxp=c(0,1,4) ,bty="n",xaxt="n",ylab="",yaxt="n")
# axis(side=1,at=c( (0.75-1.678664)/0.7459911, (1.00-1.678664)/0.7459911,(2-1.678664)/0.7459911,(3-1.678664)/0.7459911 ,(4-1.678664)/0.7459911,(5-1.678664)/0.7459911)  ,labels=FALSE)
# polygon(x=c(sizedimorphism_canine,rev(sizedimorphism_canine)),y=c(overallprobs_speciesaverages_caninedimorphism[,2],rev(overallprobs_speciesaverages_caninedimorphism[,1])),col=co_dominance_color,border=NA)
# polygon(x=c(sizedimorphism_canine,rev(sizedimorphism_canine)),y=c(overallprobs_speciesaverages_caninedimorphism[,1],rep(0,13)),col=female_dominance_color,border=NA)
# polygon(x=c(sizedimorphism_canine,rev(sizedimorphism_canine)),y=c(rep(1,13),rev(overallprobs_speciesaverages_caninedimorphism[,2])),col=male_dominance_color,border=NA)
# 
# plot( NULL , type="n" , xlab="Number of males relative to number of females" ,      xlim=c(min(dat_list_strict_sexratio$sexratio,na.rm=T),max(dat_list_strict_sexratio$sexratio,na.rm=T)) , ylim=c(0,1) ,  yaxp=c(0,1,4) ,bty="n",xaxt="n",ylab="",yaxt="n")
# axis(side=1,at=c(-2.166648,-0.068047,1,2.324043),labels=FALSE)
# polygon(x=c(sexratio,rev(sexratio)),y=c(overallprobs_speciesaverages_sexratio[,2],rev(overallprobs_speciesaverages_sexratio[,1])),col=co_dominance_color,border=NA)
# polygon(x=c(sexratio,rev(sexratio)),y=c(overallprobs_speciesaverages_sexratio[,1],rep(0,11)),col=female_dominance_color,border=NA)
# polygon(x=c(sexratio,rev(sexratio)),y=c(rep(1,11),rev(overallprobs_speciesaverages_sexratio[,2])),col=male_dominance_color,border=NA)
# 
# 
# plot( NULL , type="n" , xlab="Male reproductive skew less or more than expected by chance" ,  xlim=c(min(dat_list_strict_skew$male_skew,na.rm=T),max(dat_list_strict_skew$male_skew,na.rm=T)) , ylim=c(0,1) ,  yaxp=c(0,1,4) ,bty="n",xaxt="n",ylab="",yaxt="n")
# axis(side=1,at=c( (-2.5-0.1507726)/0.8886107, (-1.25-0.1507726)/0.8886107, (0-0.1507726)/0.8886107, (1.25-0.1507726)/0.8886107, (1.75-0.1507726)/0.8886107  ),labels=FALSE)
# polygon(x=c(male_skew,rev(male_skew)),y=c(overallprobs_speciesaverages_skew[,2],rev(overallprobs_speciesaverages_skew[,1])),col=co_dominance_color,border=NA)
# polygon(x=c(male_skew,rev(male_skew)),y=c(overallprobs_speciesaverages_skew[,1],rep(0,11)),col=female_dominance_color,border=NA)
# polygon(x=c(male_skew,rev(male_skew)),y=c(rep(1,11),rev(overallprobs_speciesaverages_skew[,2])),col=male_dominance_color,border=NA)
# 
# par<-previouspar
# 
# dev.off()
# 
# 
# 
# 
# 
# ### bottom: model output
# 
# # arboreality
# 
# speciesaverage<-combined[is.na(combined$Strata_Wilman)==F,]
# 
# dat_list_strict_arboreality <- list(
#   R = as.integer(as.factor(1/as.integer(speciesaverage$strictfdom))),
#   arboreality = as.integer(as.factor(speciesaverage$Strata_Wilman)),
#   species = as.integer(as.factor(speciesaverage$corrected_species_id))
# )
# 
# m_arboreality <- ulam(
#   alist(
#     R ~ dordlogit( phi , cutpoints ),
#     phi <-a + b[arboreality],
#     a ~ normal( 0 , 5 ),
#     b[arboreality] ~ dnorm(0,5),
#     cutpoints ~ dnorm( 0 , 5 )
#   ) , data=dat_list_strict_arboreality , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)
# 
# 
# arboreality <- seq(from=1,to=3,by=1)
# pdat <- data.frame(arboreality=arboreality)
# 
# overallphi_speciesaverages<-matrix(ncol=1, nrow=3)
# for(i in 1:length(arboreality)){
#   overallphi_speciesaverages[i,]<-precis(m_arboreality,depth=2)[1,1]+precis(m_arboreality,depth=2)[i+1,1]*arboreality[i]
# }
# 
# overallprobs_speciesaverages_arboreality<-pordlogit( 1:3 , overallphi_speciesaverages , precis(m_arboreality,depth=2)[5:6,1] )
# 
# overallprobs_speciesaverages_arboreality[,3]<-overallprobs_speciesaverages_arboreality[,3]-overallprobs_speciesaverages_arboreality[,2]
# overallprobs_speciesaverages_arboreality[,2]<-overallprobs_speciesaverages_arboreality[,2]-overallprobs_speciesaverages_arboreality[,1]
# overallprobs_speciesaverages_arboreality <-t(overallprobs_speciesaverages_arboreality)
# 
# 
# overallprobs_speciesaverages_arboreality <-overallprobs_speciesaverages_arboreality[,c(1,3,2)]
# 
# 
# 
# # ovulationsigns
# 
# speciesaverage<-combined[is.na(combined$ovulation_signs)==F,]
# 
# dat_list_strict_ovulation_signs <- list(
#   R = as.integer(as.factor(1/as.integer(speciesaverage$strictfdom))),
#   ovulation_signs = as.integer(as.factor(speciesaverage$ovulation_signs)),
#   species = as.integer(as.factor(speciesaverage$corrected_species_id))
# )
# 
# m_ovulation_signs <- ulam(
#   alist(
#     R ~ dordlogit( phi , cutpoints ),
#     phi <-a + b[ovulation_signs],
#     a ~ normal( 0 , 5 ),
#     b[ovulation_signs] ~ dnorm(0,5),
#     cutpoints ~ dnorm( 0 , 5 )
#   ) , data=dat_list_strict_ovulation_signs , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)
# 
# 
# ovulation_signs <- seq(from=1,to=3,by=1)
# pdat <- data.frame(ovulation_signs=ovulation_signs)
# 
# overallphi_speciesaverages<-matrix(ncol=1, nrow=3)
# for(i in 1:length(ovulation_signs)){
#   overallphi_speciesaverages[i,]<-precis(m_ovulation_signs,depth=2)[1,1]+precis(m_ovulation_signs,depth=2)[i+1,1]*ovulation_signs[i]
# }
# 
# overallprobs_speciesaverages_ovulation_signs<-pordlogit( 1:3 , overallphi_speciesaverages , precis(m_ovulation_signs,depth=2)[5:6,1] )
# 
# overallprobs_speciesaverages_ovulation_signs[,3]<-overallprobs_speciesaverages_ovulation_signs[,3]-overallprobs_speciesaverages_ovulation_signs[,2]
# overallprobs_speciesaverages_ovulation_signs[,2]<-overallprobs_speciesaverages_ovulation_signs[,2]-overallprobs_speciesaverages_ovulation_signs[,1]
# overallprobs_speciesaverages_ovulation_signs <-t(overallprobs_speciesaverages_ovulation_signs)
# 
# 
# overallprobs_speciesaverages_ovulation_signs <-overallprobs_speciesaverages_ovulation_signs[,c(1,3,2)]
# 
# 
# 
# 
# 
# # sexual receptivity
# 
# sexualreceptivity_hours_data<-combined[is.na(combined$sexualreceptivity_hours)==F,]
# sexualreceptivity_hours_data$strictfdom<-as.integer(sexualreceptivity_hours_data$strictfdom)
# 
# dat_list_strict_sexualreceptivity_hours <- list(
#   R = as.integer(as.factor(1/sexualreceptivity_hours_data$strictfdom)),
#   sexualreceptivity_hours = standardize(exp(sexualreceptivity_hours_data$sexualreceptivity_hours)),
#   species = as.integer(as.factor(sexualreceptivity_hours_data$corrected_species_id))
# )
# 
# m_speciesaverage_sexualreceptivity_hours <- ulam(
#   alist(
#     R ~ dordlogit( phi , cutpoints ),
#     phi <-a + b*sexualreceptivity_hours ,
#     a ~ normal( 0 , 5 ),
#     b ~ dnorm(0,5),
#     cutpoints ~ dnorm( 0 , 5 )
#   ) , data= dat_list_strict_sexualreceptivity_hours , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)
# 
# 
# 
# sexualreceptivity_hours <- seq(from=-1.5,to=3.5,by=0.5)
# pdat <- data.frame(sexualreceptivity_hours=sexualreceptivity_hours)
# 
# overallphi_speciesaverages<-matrix(ncol=1, nrow=11)
# for(i in 1:length(sexualreceptivity_hours)){
#   overallphi_speciesaverages[i,]<-precis(m_speciesaverage_sexualreceptivity_hours,depth=2)[1,1]+precis(m_speciesaverage_sexualreceptivity_hours,depth=2)[2,1]*sexualreceptivity_hours[i]
# }
# 
# overallprobs_speciesaverages_sexualreceptivity_hours<-pordlogit( 1:3 , overallphi_speciesaverages , precis(m_speciesaverage_sexualreceptivity_hours,depth=2)[3:4,1] )
# 
# 
# 
# 
# # testes size
# 
# testes_data<-combined[is.na(combined$relative_testes_mass)==F,]
# testes_data$strictfdom<-as.integer(testes_data$strictfdom)
# 
# dat_list_strict_relative_testes_mass<- list(
#   R = as.integer(as.factor(1/testes_data$strictfdom)),
#   relative_testes_mass = standardize((testes_data$relative_testes_mass)),
#   species = as.integer(as.factor(testes_data$corrected_species_id))
# )
# 
# m_speciesaverage_relative_testes_mass <- ulam(
#   alist(
#     R ~ dordlogit( phi , cutpoints ),
#     phi <-a + b*relative_testes_mass ,
#     a ~ normal( 0 , 5 ),
#     b ~ dnorm(0,5),
#     cutpoints ~ dnorm( 0 , 5 )
#   ) , data= dat_list_strict_relative_testes_mass , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)
# 
# 
# relative_testes_mass <- seq(from=-2.5,to=1.5,by=0.5)
# pdat <- data.frame(relative_testes_mass=relative_testes_mass)
# 
# overallphi_speciesaverages<-matrix(ncol=1, nrow=length(relative_testes_mass))
# for(i in 1:length(relative_testes_mass)){
#   overallphi_speciesaverages[i,]<-precis(m_speciesaverage_relative_testes_mass,depth=2)[1,1]+precis(m_speciesaverage_relative_testes_mass,depth=2)[2,1]* relative_testes_mass[i]
# }
# 
# overallprobs_speciesaverages_relative_testes_mass<-pordlogit( 1:3 , overallphi_speciesaverages , precis(m_speciesaverage_relative_testes_mass,depth=2)[3:4,1] )
# 
# 
# 
# # synchrony
# 
# speciesaverage<-combined[is.na(combined$Synchrony)==F,]
# 
# dat_list_strict_Synchrony <- list(
#   R = as.integer(as.factor(1/as.integer(speciesaverage$strictfdom))),
#   Synchrony = standardize(speciesaverage$Synchrony),
#   species = as.integer(as.factor(speciesaverage$corrected_species_id))
# )
# 
# m_speciesaverage_Synchrony <- ulam(
#   alist(
#     R ~ dordlogit( phi , cutpoints ),
#     phi <-a + b*Synchrony ,
#     a ~ normal( 0 , 5 ),
#     b ~ dnorm(0,5),
#     cutpoints ~ dnorm( 0 , 5 )
#   ) , data=dat_list_strict_Synchrony , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)
# 
# 
# Synchrony <- seq(from=-1,to=2.0,by=0.5)
# pdat <- data.frame(Synchrony=Synchrony)
# 
# overallphi_speciesaverages<-matrix(ncol=1, nrow=length(Synchrony))
# for(i in 1:length(Synchrony)){
#   overallphi_speciesaverages[i,]<-precis(m_speciesaverage_Synchrony,depth=2)[1,1]+precis(m_speciesaverage_Synchrony,depth=2)[2,1]*Synchrony[i]
# }
# 
# overallprobs_speciesaverages_synchrony<-pordlogit( 1:3 , overallphi_speciesaverages , precis(m_speciesaverage_Synchrony,depth=2)[3:4,1] )
# 
# 
# 
# 
# 
# # combined plot
# 
# 
# # Fig 2b bottom
# pdf("figures/R_Fig2b_bottom.pdf",width=18,height=4)
# previouspar<-par()
# op <- par(oma=c(0.2,0.2,0.2,0.2), mar=c(0.3,1.0,0.3,1.0), mfrow=c(1,5))
# barplot(overallprobs_speciesaverages_arboreality,col=dominance_colors,axisnames = F)
# 
# barplot(overallprobs_speciesaverages_ovulation_signs,col=c( col.alpha(dominance_colors[1],1), col.alpha(dominance_colors[2],1),col.alpha(dominance_colors[3],1)),axisnames = F,yaxt="n")
# 
# plot( NULL , type="n" , xlab="sexual receptivity" ,
#       xlim=c(min(dat_list_strict_sexualreceptivity_hours $sexualreceptivity_hours,na.rm=T),max(dat_list_strict_sexualreceptivity_hours $sexualreceptivity_hours,na.rm=T)) , ylim=c(0,1) ,  yaxp=c(0,1,4) ,bty="n",xaxt="n",ylab="",yaxt="n")
# polygon(x=c(sexualreceptivity_hours,rev(sexualreceptivity_hours)),y=c(overallprobs_speciesaverages_sexualreceptivity_hours[,2],rev(overallprobs_speciesaverages_sexualreceptivity_hours[,1])),col=co_dominance_color,border=NA)
# polygon(x=c(sexualreceptivity_hours,rev(sexualreceptivity_hours)),y=c(overallprobs_speciesaverages_sexualreceptivity_hours[,1],rep(0,length(sexualreceptivity_hours))),col=female_dominance_color,border=NA)
# polygon(x=c(sexualreceptivity_hours,rev(sexualreceptivity_hours)),y=c(rep(1,length(sexualreceptivity_hours)),rev(overallprobs_speciesaverages_sexualreceptivity_hours[,2])),col=male_dominance_color,border=NA)
# 
# plot( NULL , type="n" , xlab="testes size" ,
#       xlim=c(min(dat_list_strict_relative_testes_mass$relative_testes_mass,na.rm=T),max(dat_list_strict_relative_testes_mass$relative_testes_mass,na.rm=T)) , ylim=c(0,1) ,  yaxp=c(0,1,4) ,bty="n",xaxt="n",ylab="",yaxt="n")
# polygon(x=c(relative_testes_mass,rev(relative_testes_mass)),y=c(overallprobs_speciesaverages_relative_testes_mass[,2],rev(overallprobs_speciesaverages_relative_testes_mass[,1])),col=col.alpha(co_dominance_color,0.4),border=NA)
# polygon(x=c(relative_testes_mass,rev(relative_testes_mass)),y=c(overallprobs_speciesaverages_relative_testes_mass[,1],rep(0,length(relative_testes_mass))),col=col.alpha(female_dominance_color,0.4),border=NA)
# polygon(x=c(relative_testes_mass,rev(relative_testes_mass)),y=c(rep(1,length(relative_testes_mass)),rev(overallprobs_speciesaverages_relative_testes_mass[,2])),col=col.alpha(male_dominance_color,0.4),border=NA)
# 
# 
# plot( NULL , type="n" , xlab="synchrony" ,
#       xlim=c(min(dat_list_strict_Synchrony$Synchrony,na.rm=T),max(dat_list_strict_Synchrony$Synchrony,na.rm=T)) , ylim=c(0,1) ,  yaxp=c(0,1,4) ,bty="n",xaxt="n",ylab="",yaxt="n")
# polygon(x=c(Synchrony,rev(Synchrony)),y=c(overallprobs_speciesaverages_synchrony[,2],rev(overallprobs_speciesaverages_synchrony[,1])),col=col.alpha(co_dominance_color,0.4),border=NA)
# polygon(x=c(Synchrony,rev(Synchrony)),y=c(overallprobs_speciesaverages_synchrony[,1],rep(0,length(Synchrony))),col=col.alpha(female_dominance_color,0.4),border=NA)
# polygon(x=c(Synchrony,rev(Synchrony)),y=c(rep(1,length(Synchrony)),rev(overallprobs_speciesaverages_synchrony[,2])),col=col.alpha(male_dominance_color,0.4),border=NA)
# 
# par<-previouspar
# dev.off()
# 
# 
# 
# ### bottom: model output
# 
# # SocOrgPMK
# 
# speciesaverage<-combined[is.na(combined$SocOrgPMK)==F,]
# 
# dat_list_strict_SocOrgPMK <- list(
#   R = as.integer(as.factor(1/as.integer(speciesaverage$strictfdom))),
#   SocOrgPMK = as.integer(as.factor(speciesaverage$SocOrgPMK)),
#   species = as.integer(as.factor(speciesaverage$corrected_species_id))
# )
# 
# m_SocOrgPMK <- ulam(
#   alist(
#     R ~ dordlogit( phi , cutpoints ),
#     phi <-a + b[SocOrgPMK],
#     a ~ normal( 0 , 5 ),
#     b[SocOrgPMK] ~ dnorm(0,5),
#     cutpoints ~ dnorm( 0 , 5 )
#   ) , data=dat_list_strict_SocOrgPMK , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)
# 
# 
# SocOrgPMK <- seq(from=1,to=3,by=1)
# pdat <- data.frame(SocOrgPMK=SocOrgPMK)
# 
# overallphi_speciesaverages<-matrix(ncol=1, nrow=3)
# for(i in 1:length(SocOrgPMK)){
#   overallphi_speciesaverages[i,]<-precis(m_SocOrgPMK,depth=2)[1,1]+precis(m_SocOrgPMK,depth=2)[i+1,1]*SocOrgPMK[i]
# }
# 
# overallprobs_speciesaverages_SocOrgPMK<-pordlogit( 1:3 , overallphi_speciesaverages , precis(m_SocOrgPMK,depth=2)[5:6,1] )
# 
# overallprobs_speciesaverages_SocOrgPMK[,3]<-overallprobs_speciesaverages_SocOrgPMK[,3]-overallprobs_speciesaverages_SocOrgPMK[,2]
# overallprobs_speciesaverages_SocOrgPMK[,2]<-overallprobs_speciesaverages_SocOrgPMK[,2]-overallprobs_speciesaverages_SocOrgPMK[,1]
# overallprobs_speciesaverages_SocOrgPMK <-t(overallprobs_speciesaverages_SocOrgPMK)
# 
# 
# overallprobs_speciesaverages_SocOrgPMK <-overallprobs_speciesaverages_SocOrgPMK[,c(1,3,2)]
# 
# 
# 
# # joint aggression
# 
# speciesaverage<-combined[is.na(combined$jointaggression_females)==F,]
# 
# dat_list_jointaggression_females <- list(
#   R = as.integer(as.factor(1/as.integer(speciesaverage$strictfdom))),
#   jointaggression_females = as.integer(as.factor(speciesaverage$jointaggression_females)),
#   species = as.integer(as.factor(speciesaverage$corrected_species_id))
# )
# 
# m_jointaggression_females <- ulam(
#   alist(
#     R ~ dordlogit( phi , cutpoints ),
#     phi <-a + b[jointaggression_females],
#     a ~ normal( 0 , 5 ),
#     b[jointaggression_females] ~ dnorm(0,5),
#     cutpoints ~ dnorm( 0 , 5 )
#   ) , data=dat_list_jointaggression_females , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)
# 
# 
# jointaggression_females <- seq(from=1,to=2,by=1)
# pdat <- data.frame(jointaggression_females=jointaggression_females)
# 
# overallphi_speciesaverages<-matrix(ncol=1, nrow=2)
# for(i in 1:length(jointaggression_females)){
#   overallphi_speciesaverages[i,]<-precis(m_jointaggression_females,depth=2)[1,1]+precis(m_jointaggression_females,depth=2)[i+1,1]*jointaggression_females[i]
# }
# 
# overallprobs_speciesaverages_jointaggression_females<-pordlogit( 1:3 , overallphi_speciesaverages , precis(m_jointaggression_females,depth=2)[4:5,1] )
# 
# overallprobs_speciesaverages_jointaggression_females[,3]<-overallprobs_speciesaverages_jointaggression_females[,3]-overallprobs_speciesaverages_jointaggression_females[,2]
# overallprobs_speciesaverages_jointaggression_females[,2]<-overallprobs_speciesaverages_jointaggression_females[,2]-overallprobs_speciesaverages_jointaggression_females[,1]
# overallprobs_speciesaverages_jointaggression_females <-t(overallprobs_speciesaverages_jointaggression_females)
# 
# 
# 
# 
# 
# 
# 
# 
# # sex bias dispersal
# 
# speciesaverage<-combined[is.na(combined$sexbias_dispersal)==F,]
# 
# dat_list_sexbias_dispersal <- list(
#   R = as.integer(as.factor(1/as.integer(speciesaverage$strictfdom))),
#   sexbias_dispersal = as.integer(as.factor(speciesaverage$sexbias_dispersal)),
#   species = as.integer(as.factor(speciesaverage$corrected_species_id))
# )
# 
# m_jointaggression_sexbiasdispersal <- ulam(
#   alist(
#     R ~ dordlogit( phi , cutpoints ),
#     phi <-a + b[sexbias_dispersal],
#     a ~ normal( 0 , 5 ),
#     b[sexbias_dispersal] ~ dnorm(0,5),
#     cutpoints ~ dnorm( 0 , 5 )
#   ) , data=dat_list_sexbias_dispersal , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)
# 
# 
# sexbias_dispersal <- seq(from=1,to=3,by=1)
# pdat <- data.frame(sexbias_dispersal=sexbias_dispersal)
# 
# overallphi_speciesaverages<-matrix(ncol=1, nrow=3)
# for(i in 1:length(sexbias_dispersal)){
#   overallphi_speciesaverages[i,]<-precis(m_jointaggression_sexbiasdispersal,depth=2)[1,1]+precis(m_jointaggression_sexbiasdispersal,depth=2)[i+1,1]*sexbias_dispersal[i]
# }
# 
# overallprobs_speciesaverages_sexbias_dispersal<-pordlogit( 1:3 , overallphi_speciesaverages , precis(m_jointaggression_sexbiasdispersal,depth=2)[5:6,1] )
# 
# overallprobs_speciesaverages_sexbias_dispersal[,3]<-overallprobs_speciesaverages_sexbias_dispersal[,3]-overallprobs_speciesaverages_sexbias_dispersal[,2]
# overallprobs_speciesaverages_sexbias_dispersal[,2]<-overallprobs_speciesaverages_sexbias_dispersal[,2]-overallprobs_speciesaverages_sexbias_dispersal[,1]
# overallprobs_speciesaverages_sexbias_dispersal <-t(overallprobs_speciesaverages_sexbias_dispersal)
# 
# overallprobs_speciesaverages_sexbias_dispersal<-overallprobs_speciesaverages_sexbias_dispersal[,c(3,1,2)]
# 
# 
# 
# # female average relatedness
# 
# speciesaverage<-combined[is.na(combined$female_average_relatedness)==F,]
# 
# dat_list_strict_female_relatedness <- list(
#   R = as.integer(as.factor(1/as.integer(speciesaverage$strictfdom))),
#   female_average_relatedness = standardize(speciesaverage$female_average_relatedness),
#   species = as.integer(as.factor(speciesaverage$corrected_species_id))
# )
# 
# m_speciesaverage_female_relatedness <- ulam(
#   alist(
#     R ~ dordlogit( phi , cutpoints ),
#     phi <-a + b*female_average_relatedness ,
#     a ~ normal( 0 , 5 ),
#     b ~ dnorm(0,5),
#     cutpoints ~ dnorm( 0 , 5 )
#   ) , data=dat_list_strict_female_relatedness , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)
# 
# 
# female_average_relatedness <- seq(from=min(dat_list_strict_female_relatedness$female_average_relatedness),to=max((dat_list_strict_female_relatedness$female_average_relatedness)),length.out=9)
# pdat <- data.frame(female_average_relatedness=female_average_relatedness)
# 
# overallphi_speciesaverages<-matrix(ncol=1, nrow=length(female_average_relatedness))
# for(i in 1:length(female_average_relatedness)){
#   overallphi_speciesaverages[i,]<-precis(m_speciesaverage_female_relatedness,depth=2)[1,1]+precis(m_speciesaverage_female_relatedness,depth=2)[2,1]*female_average_relatedness[i]
# }
# 
# overallprobs_speciesaverages_female_relatedness<-pordlogit( 1:3 , overallphi_speciesaverages , precis(m_speciesaverage_female_relatedness,depth=2)[3:4,1] )
# 
# 
# 
# 
# 
# 
# # combined plot
# pdf("figures/R_Fig2c_bottom.pdf",width=14.4,height=4)
# previouspar<-par()
# op <- par(oma=c(0.2,0.2,0.2,0.2), mar=c(0.3,1.0,0.3,1.0), mfrow=c(1,4))
# barplot(overallprobs_speciesaverages_SocOrgPMK,col=dominance_colors,axisnames = F)
# 
# barplot(overallprobs_speciesaverages_jointaggression_females,col=c(col.alpha(dominance_colors[1],alpha=0.4),col.alpha(dominance_colors[2],alpha=0.4),col.alpha(dominance_colors[3],alpha=0.4)),axisnames = F,yaxt="n")
# 
# barplot(overallprobs_speciesaverages_sexbias_dispersal,col=dominance_colors,axisnames = F,yaxt="n")
# 
# plot( NULL , type="n" , xlab="female_average_relatedness" ,
#       xlim=c(min(dat_list_strict_female_relatedness$female_average_relatedness,na.rm=T),max(dat_list_strict_female_relatedness$female_average_relatedness,na.rm=T)) , ylim=c(0,1) ,  yaxp=c(0,1,4) ,bty="n",xaxt="n",ylab="",yaxt="n")
# polygon(x=c(female_average_relatedness,rev(female_average_relatedness)),y=c(overallprobs_speciesaverages_female_relatedness[,2],rev(overallprobs_speciesaverages_female_relatedness[,1])),col=col.alpha(co_dominance_color,0.4),border=NA)
# polygon(x=c(female_average_relatedness,rev(female_average_relatedness)),y=c(overallprobs_speciesaverages_female_relatedness[,1],rep(0,length(female_average_relatedness))),col=col.alpha(female_dominance_color,0.4),border=NA)
# polygon(x=c(female_average_relatedness,rev(female_average_relatedness)),y=c(rep(1,length(female_average_relatedness)),rev(overallprobs_speciesaverages_female_relatedness[,2])),col=col.alpha(male_dominance_color,0.4),border=NA)
# par<-previouspar
# 
# dev.off()
# 
# 
# 
# ### bottom: model output
# 
# # sex ratio 
# 
# speciesaverage<-combined[is.na(combined$sexratio)==F,]
# 
# dat_list_strict_sexratio <- list(
#   R = as.integer(as.factor(1/as.integer(speciesaverage$strictfdom))),
#   sexratio = standardize(speciesaverage$sexratio),
#   species = as.integer(as.factor(speciesaverage$corrected_species_id))
# )
# 
# m_speciesaverage_sexratio <- ulam(
#   alist(
#     R ~ dordlogit( phi , cutpoints ),
#     phi <-a + b*sexratio ,
#     a ~ normal( 0 , 5 ),
#     b ~ dnorm(0,5),
#     cutpoints ~ dnorm( 0 , 5 )
#   ) , data=dat_list_strict_sexratio , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)
# 
# 
# sexratio <- seq(from=-2.5,to=2.5,by=0.5)
# pdat <- data.frame(sexratio=sexratio)
# 
# overallphi_speciesaverages<-matrix(ncol=1, nrow=11)
# for(i in 1:length(sexratio)){
#   overallphi_speciesaverages[i,]<-precis(m_speciesaverage_sexratio,depth=2)[1,1]+precis(m_speciesaverage_sexratio,depth=2)[2,1]*sexratio[i]
# }
# 
# overallprobs_speciesaverages_sexratio<-pordlogit( 1:3 , overallphi_speciesaverages , precis(m_speciesaverage_sexratio,depth=2)[3:4,1] )
# 
# 
# 
# 
# # percent aggression male-male
# 
# speciesaverage<-combined[is.na(combined$perc_aggression_mm)==F,]
# 
# dat_list_strict_perc_aggression_mm <- list(
#   R = as.integer(as.factor(1/as.integer(speciesaverage$strictfdom))),
#   perc_aggression_mm = standardize(speciesaverage$perc_aggression_mm),
#   species = as.integer(as.factor(speciesaverage$corrected_species_id))
# )
# 
# m_speciesaverage_perc_aggression_mm <- ulam(
#   alist(
#     R ~ dordlogit( phi , cutpoints ),
#     phi <-a + b*perc_aggression_mm ,
#     a ~ normal( 0 , 5 ),
#     b ~ dnorm(0,5),
#     cutpoints ~ dnorm( 0 , 5 )
#   ) , data=dat_list_strict_perc_aggression_mm , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)
# 
# 
# perc_aggression_mm <- seq(from=min(dat_list_strict_perc_aggression_mm$perc_aggression_mm),to=max((dat_list_strict_perc_aggression_mm$perc_aggression_mm)),length.out=9)
# pdat <- data.frame(perc_aggression_mm=perc_aggression_mm)
# 
# overallphi_speciesaverages<-matrix(ncol=1, nrow=length(perc_aggression_mm))
# for(i in 1:length(perc_aggression_mm)){
#   overallphi_speciesaverages[i,]<-precis(m_speciesaverage_perc_aggression_mm,depth=2)[1,1]+precis(m_speciesaverage_perc_aggression_mm,depth=2)[2,1]*perc_aggression_mm[i]
# }
# 
# overallprobs_speciesaverages_perc_aggression_mm<-pordlogit( 1:3 , overallphi_speciesaverages , precis(m_speciesaverage_perc_aggression_mm,depth=2)[3:4,1] )
# 
# 
# 
# # number males
# 
# males_data<-combined[is.na(combined$males)==F,]
# males_data$strictfdom<-as.integer(males_data$strictfdom)
# 
# dat_list_strict_males<- list(
#   R = as.integer(as.factor(1/males_data$strictfdom)),
#   males = standardize((males_data$males)),
#   species = as.integer(as.factor(males_data$corrected_species_id))
# )
# 
# m_speciesaverage_relative_males <- ulam(
#   alist(
#     R ~ dordlogit( phi , cutpoints ),
#     phi <-a + b*males ,
#     a ~ normal( 0 , 5 ),
#     b ~ dnorm(0,5),
#     cutpoints ~ dnorm( 0 , 5 )
#   ) , data= dat_list_strict_males , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)
# 
# 
# males <- seq(from=min(dat_list_strict_males$males),to=max(dat_list_strict_males$males),length.out=9)
# pdat <- data.frame(males=males)
# 
# overallphi_speciesaverages<-matrix(ncol=1, nrow=length(males))
# for(i in 1:length(males)){
#   overallphi_speciesaverages[i,]<-precis(m_speciesaverage_relative_males,depth=2)[1,1]+precis(m_speciesaverage_relative_males,depth=2)[2,1]* males[i]
# }
# 
# overallprobs_speciesaverages_relative_males<-pordlogit( 1:3 , overallphi_speciesaverages , precis(m_speciesaverage_relative_males,depth=2)[3:4,1] )
# 
# 
# 
# # combined plot
# pdf("figures/R_Fig2d_bottom.pdf",width=10.8,height=4)
# previouspar<-par()
# op <- par(oma=c(0.2,0.2,0.2,0.2), mar=c(0.3,1.0,0.3,1.0), mfrow=c(1,3))
# 
# plot( NULL , type="n" , xlab="Number of males relative to number of females" ,      xlim=c(min(dat_list_strict_sexratio$sexratio,na.rm=T),max(dat_list_strict_sexratio$sexratio,na.rm=T)) , ylim=c(0,1) ,  yaxp=c(0,1,4) ,bty="n",xaxt="n",ylab="",yaxt="n")
# axis(side=1,at=c(-2.166648,-0.068047,1,2.324043),labels=FALSE)
# polygon(x=c(sexratio,rev(sexratio)),y=c(overallprobs_speciesaverages_sexratio[,2],rev(overallprobs_speciesaverages_sexratio[,1])),col=co_dominance_color,border=NA)
# polygon(x=c(sexratio,rev(sexratio)),y=c(overallprobs_speciesaverages_sexratio[,1],rep(0,11)),col=female_dominance_color,border=NA)
# polygon(x=c(sexratio,rev(sexratio)),y=c(rep(1,11),rev(overallprobs_speciesaverages_sexratio[,2])),col=male_dominance_color,border=NA)
# 
# 
# plot( NULL , type="n" , xlab="perc_aggression_mm" ,
#       xlim=c(min(dat_list_strict_perc_aggression_mm$perc_aggression_mm,na.rm=T),max(dat_list_strict_perc_aggression_mm$perc_aggression_mm,na.rm=T)) , ylim=c(0,1) ,  yaxp=c(0,1,4) ,bty="n",xaxt="n",ylab="",yaxt="n")
# polygon(x=c(perc_aggression_mm,rev(perc_aggression_mm)),y=c(overallprobs_speciesaverages_perc_aggression_mm[,2],rev(overallprobs_speciesaverages_perc_aggression_mm[,1])),col=col.alpha(co_dominance_color,1),border=NA)
# polygon(x=c(perc_aggression_mm,rev(perc_aggression_mm)),y=c(overallprobs_speciesaverages_perc_aggression_mm[,1],rep(0,length(perc_aggression_mm))),col=col.alpha(female_dominance_color,1),border=NA)
# polygon(x=c(perc_aggression_mm,rev(perc_aggression_mm)),y=c(rep(1,length(perc_aggression_mm)),rev(overallprobs_speciesaverages_perc_aggression_mm[,2])),col=col.alpha(male_dominance_color,1),border=NA)
# par<-previouspar
# 
# plot( NULL , type="n" , xlab="males" ,
#       xlim=c(min(dat_list_strict_males $males,na.rm=T),max(dat_list_strict_males $males,na.rm=T)) , ylim=c(0,1) ,  yaxp=c(0,1,4) ,bty="n",xaxt="n",ylab="",yaxt="n")
# polygon(x=c(males,rev(males)),y=c(overallprobs_speciesaverages_relative_males[,2],rev(overallprobs_speciesaverages_relative_males[,1])),col=co_dominance_color,border=NA)
# polygon(x=c(males,rev(males)),y=c(overallprobs_speciesaverages_relative_males[,1],rep(0,length(males))),col=female_dominance_color,border=NA)
# polygon(x=c(males,rev(males)),y=c(rep(1,length(males)),rev(overallprobs_speciesaverages_relative_males[,2])),col=male_dominance_color,border=NA)
# 
# 
# dev.off()
# 
# 
# ### bottom: model output
# 
# # origin
# 
# speciesaverage<-combined[is.na(combined$origin)==F,]
# speciesaverage<-speciesaverage[speciesaverage$origin !="provisioned",]
# 
# dat_list_strict_origin <- list(
#   R = (as.factor(1/as.integer(speciesaverage$strictfdom))),
#   origin = as.integer(as.factor(speciesaverage$origin)),
#   species = as.integer(as.factor(speciesaverage$corrected_species_id))
# )
# 
# m_origin <- ulam(
#   alist(
#     R ~ dordlogit( phi , cutpoints ),
#     phi <-a + b[origin],
#     a ~ normal( 0 , 5 ),
#     b[origin] ~ dnorm(0,5),
#     cutpoints ~ dnorm( 0 , 5 )
#   ) , data=dat_list_strict_origin , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)
# 
# 
# origin <- seq(from=1,to=2,by=1)
# pdat <- data.frame(origin=origin)
# 
# overallphi_speciesaverages<-matrix(ncol=1, nrow=2)
# for(i in 1:length(origin)){
#   overallphi_speciesaverages[i,]<-precis(m_origin,depth=2)[1,1]+precis(m_origin,depth=2)[i+1,1]*origin[i]
# }
# 
# overallprobs_speciesaverages_origin<-pordlogit( 1:3 , overallphi_speciesaverages , precis(m_origin,depth=2)[4:5,1] )
# 
# overallprobs_speciesaverages_origin[,3]<-overallprobs_speciesaverages_origin[,3]-overallprobs_speciesaverages_origin[,2]
# overallprobs_speciesaverages_origin[,2]<-overallprobs_speciesaverages_origin[,2]-overallprobs_speciesaverages_origin[,1]
# overallprobs_speciesaverages_origin <-t(overallprobs_speciesaverages_origin)
# 
# 
# 
# # environmental harshness
# 
# env_harshness_data<-combined[is.na(combined$env_harshness)==F,]
# env_harshness_data$strictfdom<-as.integer(env_harshness_data$strictfdom)
# 
# dat_list_strict_env_harshness <- list(
#   R = as.integer(as.factor(1/env_harshness_data$strictfdom)),
#   env_harshness = standardize((env_harshness_data$env_harshness)),
#   species = as.integer(as.factor(env_harshness_data$corrected_species_id))
# )
# 
# m_speciesaverage_env_harshness <- ulam(
#   alist(
#     R ~ dordlogit( phi , cutpoints ),
#     phi <-a + b*env_harshness ,
#     a ~ normal( 0 , 5 ),
#     b ~ dnorm(0,5),
#     cutpoints ~ dnorm( 0 , 5 )
#   ) , data= dat_list_strict_env_harshness , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)
# 
# 
# 
# env_harshness <- seq(from=min(dat_list_strict_env_harshness$env_harshness),to=max(dat_list_strict_env_harshness$env_harshness),length.out=9)
# pdat <- data.frame(env_harshness=env_harshness)
# 
# overallphi_speciesaverages<-matrix(ncol=1, nrow=9)
# for(i in 1:length(env_harshness)){
#   overallphi_speciesaverages[i,]<-precis(m_speciesaverage_env_harshness,depth=2)[1,1]+precis(m_speciesaverage_env_harshness,depth=2)[2,1]*env_harshness[i]
# }
# 
# overallprobs_speciesaverages_env_harshness<-pordlogit( 1:3 , overallphi_speciesaverages , precis(m_speciesaverage_env_harshness,depth=2)[3:4,1] )
# 
# 
# # rainfall variation
# 
# rainfallvariation_data<-combined[is.na(combined$rainfall_annualvariation)==F,]
# rainfallvariation_data$strictfdom<-as.integer(rainfallvariation_data$strictfdom)
# 
# dat_list_strict_rainfallvariation <- list(
#   R = as.integer(as.factor(1/rainfallvariation_data$strictfdom)),
#   rainfall_annualvariation = standardize((rainfallvariation_data$rainfall_annualvariation)),
#   species = as.integer(as.factor(rainfallvariation_data$corrected_species_id))
# )
# 
# m_speciesaverage_rainfallvariation <- ulam(
#   alist(
#     R ~ dordlogit( phi , cutpoints ),
#     phi <-a + b*rainfall_annualvariation ,
#     a ~ normal( 0 , 5 ),
#     b ~ dnorm(0,5),
#     cutpoints ~ dnorm( 0 , 5 )
#   ) , data= dat_list_strict_rainfallvariation , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)
# 
# 
# 
# rainfall_annualvariation <- seq(from=min(dat_list_strict_rainfallvariation$rainfall_annualvariation),to=max(dat_list_strict_rainfallvariation$rainfall_annualvariation),length.out=9)
# pdat <- data.frame(rainfall_annualvariation=rainfall_annualvariation)
# 
# overallphi_speciesaverages<-matrix(ncol=1, nrow=9)
# for(i in 1:length(rainfall_annualvariation)){
#   overallphi_speciesaverages[i,]<-precis(m_speciesaverage_rainfallvariation,depth=2)[1,1]+precis(m_speciesaverage_rainfallvariation,depth=2)[2,1]*rainfall_annualvariation[i]
# }
# 
# overallprobs_speciesaverages_rainfallvariation<-pordlogit( 1:3 , overallphi_speciesaverages , precis(m_speciesaverage_rainfallvariation,depth=2)[3:4,1] )
# 
# 
# 
# # rainfall predictability
# 
# rainfallpredictability_data<-combined[is.na(combined$rainfall_unpredictability)==F,]
# rainfallpredictability_data$strictfdom<-as.integer(rainfallpredictability_data$strictfdom)
# 
# dat_list_strict_rainfallpredictability <- list(
#   R = as.integer(as.factor(1/rainfallpredictability_data$strictfdom)),
#   rainfall_unpredictability = standardize((rainfallpredictability_data$rainfall_unpredictability)),
#   species = as.integer(as.factor(rainfallpredictability_data$corrected_species_id))
# )
# 
# m_speciesaverage_rainfallpredictability <- ulam(
#   alist(
#     R ~ dordlogit( phi , cutpoints ),
#     phi <-a + b*rainfall_unpredictability ,
#     a ~ normal( 0 , 5 ),
#     b ~ dnorm(0,5),
#     cutpoints ~ dnorm( 0 , 5 )
#   ) , data= dat_list_strict_rainfallpredictability , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)
# 
# rainfall_unpredictability <- seq(from=min(dat_list_strict_rainfallpredictability$rainfall_unpredictability),to=max(dat_list_strict_rainfallpredictability$rainfall_unpredictability),length.out=9)
# pdat <- data.frame(rainfall_unpredictability=rainfall_unpredictability)
# 
# overallphi_speciesaverages<-matrix(ncol=1, nrow=9)
# for(i in 1:length(rainfall_unpredictability)){
#   overallphi_speciesaverages[i,]<-precis(m_speciesaverage_rainfallpredictability,depth=2)[1,1]+precis(m_speciesaverage_rainfallpredictability,depth=2)[2,1]*rainfall_unpredictability[i]
# }
# 
# overallprobs_speciesaverages_rainfallpredictability<-pordlogit( 1:3 , overallphi_speciesaverages , precis(m_speciesaverage_rainfallpredictability,depth=2)[3:4,1] )
# 
# 
# 
# # combined plot
# pdf("figures/R_Fig2e_bottom.pdf",width=14.4,height=4)
# previouspar<-par()
# op <- par(oma=c(0.2,0.2,0.2,0.2), mar=c(0.3,1.0,0.3,1.0), mfrow=c(1,4))
# barplot(overallprobs_speciesaverages_origin,col=c(col.alpha(dominance_colors[1],0.4),col.alpha(dominance_colors[2],0.4),col.alpha(dominance_colors[3],0.4) ),axisnames = F)
# 
# plot( NULL , type="n" , xlab="environmental harshness" ,
#       xlim=c(min(dat_list_strict_env_harshness$env_harshness,na.rm=T),max(dat_list_strict_env_harshness$env_harshness,na.rm=T)) , ylim=c(0,1) ,  yaxp=c(0,1,4) ,bty="n",xaxt="n",ylab="",yaxt="n")
# polygon(x=c(env_harshness,rev(env_harshness)),y=c(overallprobs_speciesaverages_env_harshness[,2],rev(overallprobs_speciesaverages_env_harshness[,1])),col=co_dominance_color,border=NA)
# polygon(x=c(env_harshness,rev(env_harshness)),y=c(overallprobs_speciesaverages_env_harshness[,1],rep(0,length(env_harshness))),col=female_dominance_color,border=NA)
# polygon(x=c(env_harshness,rev(env_harshness)),y=c(rep(1,length(env_harshness)),rev(overallprobs_speciesaverages_env_harshness[,2])),col=male_dominance_color,border=NA)
# 
# 
# plot( NULL , type="n" , xlab="rainfall seasonality" ,
#       xlim=c(min(dat_list_strict_rainfallvariation$rainfall_annualvariation,na.rm=T),max(dat_list_strict_rainfallvariation$rainfall_annualvariation,na.rm=T)) , ylim=c(0,1) ,  yaxp=c(0,1,4) ,bty="n",xaxt="n",ylab="",yaxt="n")
# polygon(x=c(rainfall_annualvariation,rev(rainfall_annualvariation)),y=c(overallprobs_speciesaverages_rainfallvariation[,2],rev(overallprobs_speciesaverages_rainfallvariation[,1])),col=col.alpha(co_dominance_color,0.4),border=NA)
# polygon(x=c(rainfall_annualvariation,rev(rainfall_annualvariation)),y=c(overallprobs_speciesaverages_rainfallvariation[,1],rep(0,length(rainfall_annualvariation))),col=col.alpha(female_dominance_color,0.4),border=NA)
# polygon(x=c(rainfall_annualvariation,rev(rainfall_annualvariation)),y=c(rep(1,length(rainfall_annualvariation)),rev(overallprobs_speciesaverages_rainfallvariation[,2])),col=col.alpha(male_dominance_color,0.4),border=NA)
# 
# 
# plot( NULL , type="n" , xlab="rainfall unpredictability" ,    xlim=c(min(dat_list_strict_rainfallpredictability$rainfall_unpredictability,na.rm=T),max(dat_list_strict_rainfallpredictability$rainfall_unpredictability,na.rm=T)) , ylim=c(0,1) ,  yaxp=c(0,1,4) ,bty="n",xaxt="n",ylab="",yaxt="n")
# polygon(x=c(rainfall_unpredictability,rev(rainfall_unpredictability)),y=c(overallprobs_speciesaverages_rainfallpredictability[,2],rev(overallprobs_speciesaverages_rainfallpredictability[,1])),col=co_dominance_color,border=NA)
# polygon(x=c(rainfall_unpredictability,rev(rainfall_unpredictability)),y=c(overallprobs_speciesaverages_rainfallpredictability[,1],rep(0,length(rainfall_unpredictability))),col=female_dominance_color,border=NA)
# polygon(x=c(rainfall_unpredictability,rev(rainfall_unpredictability)),y=c(rep(1,length(rainfall_unpredictability)),rev(overallprobs_speciesaverages_rainfallpredictability[,2])),col=male_dominance_color,border=NA)
# 
# 
# par<-previouspar
# 
# dev.off()
# 
# ### bottom: model output
# 
# # female evictions
# 
# speciesaverage<-combined[is.na(combined$female_evictions)==F,]
# 
# dat_list_strict_female_evictions <- list(
#   R = as.integer(as.factor(1/as.integer(speciesaverage$strictfdom))),
#   female_evictions = ifelse(speciesaverage$female_evictions=="Yes",1,2),
#   species = as.integer(as.factor(speciesaverage$corrected_species_id))
# )
# 
# m_female_evictions <- ulam(
#   alist(
#     R ~ dordlogit( phi , cutpoints ),
#     phi <-a + b[female_evictions],
#     a ~ normal( 0 , 5 ),
#     b[female_evictions] ~ dnorm(0,5),
#     cutpoints ~ dnorm( 0 , 5 )
#   ) , data=dat_list_strict_female_evictions , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)
# 
# 
# female_evictions <- seq(from=1,to=2,by=1)
# pdat <- data.frame(female_evictions=female_evictions)
# 
# overallphi_speciesaverages<-matrix(ncol=1, nrow=2)
# for(i in 1:length(female_evictions)){
#   overallphi_speciesaverages[i,]<-precis(m_female_evictions,depth=2)[1,1]+precis(m_female_evictions,depth=2)[i+1,1]*female_evictions[i]
# }
# 
# overallprobs_speciesaverages_female_evictions<-pordlogit( 1:3 , overallphi_speciesaverages , precis(m_female_evictions,depth=2)[4:5,1] )
# 
# overallprobs_speciesaverages_female_evictions[,3]<-overallprobs_speciesaverages_female_evictions[,3]-overallprobs_speciesaverages_female_evictions[,2]
# overallprobs_speciesaverages_female_evictions[,2]<-overallprobs_speciesaverages_female_evictions[,2]-overallprobs_speciesaverages_female_evictions[,1]
# overallprobs_speciesaverages_female_evictions <-t(overallprobs_speciesaverages_female_evictions)
# 
# 
# # female infanticide
# 
# speciesaverage<-combined[is.na(combined$female_infanticide)==F,]
# 
# dat_list_strict_female_infanticide <- list(
#   R = as.integer(as.factor(1/as.integer(speciesaverage$strictfdom))),
#   female_infanticide = ifelse(speciesaverage$female_infanticide=="Yes",1,2),
#   species = as.integer(as.factor(speciesaverage$corrected_species_id))
# )
# 
# m_female_infanticide <- ulam(
#   alist(
#     R ~ dordlogit( phi , cutpoints ),
#     phi <-a + b[female_infanticide],
#     a ~ normal( 0 , 5 ),
#     b[female_infanticide] ~ dnorm(0,5),
#     cutpoints ~ dnorm( 0 , 5 )
#   ) , data=dat_list_strict_female_infanticide , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)
# 
# 
# female_infanticide <- seq(from=1,to=2,by=1)
# pdat <- data.frame(female_infanticide=female_infanticide)
# 
# overallphi_speciesaverages<-matrix(ncol=1, nrow=2)
# for(i in 1:length(female_infanticide)){
#   overallphi_speciesaverages[i,]<-precis(m_female_infanticide,depth=2)[1,1]+precis(m_female_infanticide,depth=2)[i+1,1]*female_infanticide[i]
# }
# 
# overallprobs_speciesaverages_female_infanticide<-pordlogit( 1:3 , overallphi_speciesaverages , precis(m_female_infanticide,depth=2)[4:5,1] )
# 
# overallprobs_speciesaverages_female_infanticide[,3]<-overallprobs_speciesaverages_female_infanticide[,3]-overallprobs_speciesaverages_female_infanticide[,2]
# overallprobs_speciesaverages_female_infanticide[,2]<-overallprobs_speciesaverages_female_infanticide[,2]-overallprobs_speciesaverages_female_infanticide[,1]
# overallprobs_speciesaverages_female_infanticide <-t(overallprobs_speciesaverages_female_infanticide)
# 
# 
# 
# # seasonal breeding
# 
# r_seasonality_value_data<-combined[is.na(combined$r_seasonality_value)==F,]
# r_seasonality_value_data$strictfdom<-as.integer(r_seasonality_value_data$strictfdom)
# 
# dat_list_strict_r_seasonality_value<- list(
#   R = as.integer(as.factor(1/r_seasonality_value_data$strictfdom)),
#   r_seasonality_value = standardize((r_seasonality_value_data$r_seasonality_value)),
#   species = as.integer(as.factor(r_seasonality_value_data$corrected_species_id))
# )
# 
# m_speciesaverage_r_seasonality_value <- ulam(
#   alist(
#     R ~ dordlogit( phi , cutpoints ),
#     phi <-a + b*r_seasonality_value ,
#     a ~ normal( 0 , 5 ),
#     b ~ dnorm(0,5),
#     cutpoints ~ dnorm( 0 , 5 )
#   ) , data= dat_list_strict_r_seasonality_value, chains=4 , cores=4 ,cmdstan=T, messages=FALSE)
# 
# 
# r_seasonality_value <- seq(from=min(dat_list_strict_r_seasonality_value$r_seasonality_value),to=max(dat_list_strict_r_seasonality_value$r_seasonality_value),length.out=9)
# pdat <- data.frame(r_seasonality_value=r_seasonality_value)
# 
# overallphi_speciesaverages<-matrix(ncol=1, nrow=length(r_seasonality_value))
# for(i in 1:length(r_seasonality_value)){
#   overallphi_speciesaverages[i,]<-precis(m_speciesaverage_r_seasonality_value,depth=2)[1,1]+precis(m_speciesaverage_r_seasonality_value,depth=2)[2,1]* r_seasonality_value[i]
# }
# 
# overallprobs_speciesaverages_r_seasonality_value<-pordlogit( 1:3 , overallphi_speciesaverages , precis(m_speciesaverage_r_seasonality_value,depth=2)[3:4,1] )
# 
# 
# 
# # home range overlap
# 
# speciesaverage<-combined[is.na(combined$homerange_overlap)==F,]
# 
# dat_list_strict_homerange_overlap <- list(
#   R = as.integer(as.factor(1/as.integer(speciesaverage$strictfdom))),
#   homerange_overlap = standardize(speciesaverage$homerange_overlap),
#   species = as.integer(as.factor(speciesaverage$corrected_species_id))
# )
# 
# m_speciesaverage_homerange_overlap <- ulam(
#   alist(
#     R ~ dordlogit( phi , cutpoints ),
#     phi <-a + b*homerange_overlap ,
#     a ~ normal( 0 , 5 ),
#     b ~ dnorm(0,5),
#     cutpoints ~ dnorm( 0 , 5 )
#   ) , data=dat_list_strict_homerange_overlap , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)
# 
# 
# homerange_overlap <- seq(from=min(dat_list_strict_homerange_overlap$homerange_overlap),to=max(dat_list_strict_homerange_overlap$homerange_overlap),length=10)
# pdat <- data.frame(homerange_overlap=homerange_overlap)
# 
# overallphi_speciesaverages<-matrix(ncol=1, nrow=length(homerange_overlap))
# for(i in 1:length(homerange_overlap)){
#   overallphi_speciesaverages[i,]<-precis(m_speciesaverage_homerange_overlap,depth=2)[1,1]+precis(m_speciesaverage_homerange_overlap,depth=2)[2,1]*homerange_overlap[i]
# }
# 
# overallprobs_speciesaverages_homerange_overlap<-pordlogit( 1:3 , overallphi_speciesaverages , precis(m_speciesaverage_homerange_overlap,depth=2)[3:4,1] )
# 
# 
# 
# # females
# 
# speciesaverage<-combined[is.na(combined$females)==F,]
# 
# dat_list_strict_females <- list(
#   R = as.integer(as.factor(1/as.integer(speciesaverage$strictfdom))),
#   females = standardize(speciesaverage$females),
#   species = as.integer(as.factor(speciesaverage$corrected_species_id))
# )
# 
# m_speciesaverage_females <- ulam(
#   alist(
#     R ~ dordlogit( phi , cutpoints ),
#     phi <-a + b*females ,
#     a ~ normal( 0 , 5 ),
#     b ~ dnorm(0,5),
#     cutpoints ~ dnorm( 0 , 5 )
#   ) , data=dat_list_strict_females , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)
# 
# 
# females <- seq(from=min(dat_list_strict_females$females),to=max(dat_list_strict_females$females),length=10)
# pdat <- data.frame(females=females)
# 
# overallphi_speciesaverages<-matrix(ncol=1, nrow=length(females))
# for(i in 1:length(females)){
#   overallphi_speciesaverages[i,]<-precis(m_speciesaverage_females,depth=2)[1,1]+precis(m_speciesaverage_females,depth=2)[2,1]*females[i]
# }
# 
# overallprobs_speciesaverages_females<-pordlogit( 1:3 , overallphi_speciesaverages , precis(m_speciesaverage_females,depth=2)[3:4,1] )
# 
# 
# # canine size
# 
# speciesaverage<-combined[is.na(combined$relative_femalecaninesize)==F,]
# 
# dat_list_strict_caninesize <- list(
#   R = as.integer(as.factor(1/as.integer(speciesaverage$strictfdom))),
#   relative_femalecaninesize = standardize(speciesaverage$relative_femalecaninesize),
#   species = as.integer(as.factor(speciesaverage$corrected_species_id))
# )
# 
# m_speciesaverage_caninesize <- ulam(
#   alist(
#     R ~ dordlogit( phi , cutpoints ),
#     phi <-a + b*relative_femalecaninesize ,
#     a ~ normal( 0 , 5 ),
#     b ~ dnorm(0,5),
#     cutpoints ~ dnorm( 0 , 5 )
#   ) , data=dat_list_strict_caninesize , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)
# 
# 
# relative_femalecaninesize <- seq(from=min(dat_list_strict_caninesize$relative_femalecaninesize),to=max(dat_list_strict_caninesize$relative_femalecaninesize),by=0.5)
# pdat <- data.frame(relative_femalecaninesize=relative_femalecaninesize)
# 
# overallphi_speciesaverages<-matrix(ncol=1, nrow=length(relative_femalecaninesize))
# for(i in 1:length(relative_femalecaninesize)){
#   overallphi_speciesaverages[i,]<-precis(m_speciesaverage_caninesize,depth=2)[1,1]+precis(m_speciesaverage_caninesize,depth=2)[2,1]*relative_femalecaninesize[i]
# }
# 
# overallprobs_speciesaverages_caninesize<-pordlogit( 1:3 , overallphi_speciesaverages , precis(m_speciesaverage_caninesize,depth=2)[3:4,1] )
# 
# 
# 
# pdf("figures/R_Fig2f_bottom.pdf",width=21.6,height=4)
# previouspar<-par()
# op <- par(oma=c(0.2,0.2,0.2,0.2), mar=c(0.3,1.0,0.3,1.0), mfrow=c(1,6))
# barplot(overallprobs_speciesaverages_female_evictions,col=c(col.alpha(dominance_colors[1],1),col.alpha(dominance_colors[2],1),col.alpha(dominance_colors[3],1) ),axisnames = F)
# 
# barplot(overallprobs_speciesaverages_female_infanticide,col=c(col.alpha(dominance_colors[1],0.4),col.alpha(dominance_colors[2],0.4),col.alpha(dominance_colors[3],0.4) ),axisnames = F)
# 
# plot( NULL , type="n" , xlab="seasonal breeding" ,
#       xlim=c(min(dat_list_strict_r_seasonality_value$r_seasonality_value,na.rm=T),max(dat_list_strict_r_seasonality_value$r_seasonality_value,na.rm=T)) , ylim=c(0,1) ,  yaxp=c(0,1,4) ,bty="n",xaxt="n",ylab="",yaxt="n")
# polygon(x=c(r_seasonality_value,rev(r_seasonality_value)),y=c(overallprobs_speciesaverages_r_seasonality_value[,2],rev(overallprobs_speciesaverages_r_seasonality_value[,1])),col=co_dominance_color,border=NA)
# polygon(x=c(r_seasonality_value,rev(r_seasonality_value)),y=c(overallprobs_speciesaverages_r_seasonality_value[,1],rep(0,length(r_seasonality_value))),col=female_dominance_color,border=NA)
# polygon(x=c(r_seasonality_value,rev(r_seasonality_value)),y=c(rep(1,length(r_seasonality_value)),rev(overallprobs_speciesaverages_r_seasonality_value[,2])),col=male_dominance_color,border=NA)
# 
# 
# plot( NULL , type="n" , xlab="home range overlap" ,
#       xlim=c(min(dat_list_strict_homerange_overlap$homerange_overlap,na.rm=T),max(dat_list_strict_homerange_overlap$homerange_overlap,na.rm=T)) , ylim=c(0,1) ,  yaxp=c(0,1,4) ,bty="n",xaxt="n",ylab="",yaxt="n")
# polygon(x=c(homerange_overlap,rev(homerange_overlap)),y=c(overallprobs_speciesaverages_homerange_overlap[,2],rev(overallprobs_speciesaverages_homerange_overlap[,1])),col=co_dominance_color,border=NA)
# polygon(x=c(homerange_overlap,rev(homerange_overlap)),y=c(overallprobs_speciesaverages_homerange_overlap[,1],rep(0,length(homerange_overlap))),col=female_dominance_color,border=NA)
# polygon(x=c(homerange_overlap,rev(homerange_overlap)),y=c(rep(1,length(homerange_overlap)),rev(overallprobs_speciesaverages_homerange_overlap[,2])),col=male_dominance_color,border=NA)
# 
# 
# plot( NULL , type="n" , xlab="number of females" ,    xlim=c(min(dat_list_strict_females$females,na.rm=T),max(dat_list_strict_females$females,na.rm=T)) , ylim=c(0,1) ,  yaxp=c(0,1,4) ,bty="n",xaxt="n",ylab="",yaxt="n")
# polygon(x=c(females,rev(females)),y=c(overallprobs_speciesaverages_females[,2],rev(overallprobs_speciesaverages_females[,1])),col=co_dominance_color,border=NA)
# polygon(x=c(females,rev(females)),y=c(overallprobs_speciesaverages_females[,1],rep(0,length(females))),col=female_dominance_color,border=NA)
# polygon(x=c(females,rev(females)),y=c(rep(1,length(females)),rev(overallprobs_speciesaverages_females[,2])),col=male_dominance_color,border=NA)
# 
# 
# plot( NULL , type="n" , xlab="relative canine size" ,    xlim=c(min(dat_list_strict_caninesize$relative_femalecaninesize,na.rm=T),max(dat_list_strict_caninesize$relative_femalecaninesize,na.rm=T)) , ylim=c(0,1) ,  yaxp=c(0,1,4) ,bty="n",xaxt="n",ylab="",yaxt="n")
# polygon(x=c(relative_femalecaninesize,rev(relative_femalecaninesize)),y=c(overallprobs_speciesaverages_caninesize[,2],rev(overallprobs_speciesaverages_caninesize[,1])),col=co_dominance_color,border=NA)
# polygon(x=c(relative_femalecaninesize,rev(relative_femalecaninesize)),y=c(overallprobs_speciesaverages_caninesize[,1],rep(0,length(relative_femalecaninesize))),col=female_dominance_color,border=NA)
# polygon(x=c(relative_femalecaninesize,rev(relative_femalecaninesize)),y=c(rep(1,length(relative_femalecaninesize)),rev(overallprobs_speciesaverages_caninesize[,2])),col=male_dominance_color,border=NA)
# 
# 
# 
# par<-previouspar
# 
# dev.off()
