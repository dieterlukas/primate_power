library(viridis)
library(ggtree)
library(caper)
library(tidyr)
library(ggstance)
library(ggdist)
library(cowplot)

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
fightswondata_range_min[fightswondata_range_min$delta_long %in% "0.2" ,]$delta_long<-2
fightswondata_range_min[fightswondata_range_min$delta_long %in% "0.41" ,]$delta_long<-2
fightswondata_range_min[fightswondata_range_min$min_perc_won_females %in% 100 & fightswondata_range_min$delta_long %in% 2,]$min_perc_won_females<-98


fightdistributiondata<-as.data.frame(combined %>% group_by(corrected_species_id) %>% summarize(fmaggression=mean(perc_aggression_fm,na.rm=T),ffaggression=mean(perc_aggression_ff,na.rm=T),mmaggression=mean(perc_aggression_mm,na.rm=T) ))
rownames(fightdistributiondata)<-fightswondata$Species

fightdistributiondata<-data.frame(fightdistributiondata[speciesnames,])
colnames(fightdistributiondata)<-c("Species","between the sexes","among females","among males")
row.names(fightdistributiondata)<-speciesnames

fightdistributiondata_long<-as.data.frame(pivot_longer(fightdistributiondata,cols=c("between the sexes","among females","among males")))
fightdistributiondata_long$numericspecies<-NA
for (i in 1:117){
  fightdistributiondata_long[ ((i-1)*3+1):((i-1)*3+3),]$numericspecies<-i
}


fightdistributiondata_long[is.na(fightdistributiondata_long$value)==T,]$value<-0


p3 <- facet_plot(treeplot, panel="Distribution of aggression", data=fightdistributiondata_long, geom=geom_barh, mapping=aes(fill = name,y = numericspecies, x = value),stat="identity", cex=1,fill=rep(c(fm_aggression,mm_aggression,ff_aggression),117))




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

p3 <- facet_plot(treeplot, panel="Intersexual dominance", data=dominancedata, geom=geom_point, aes(x=placeholder), color=c(dominancedata$colours),pch=15,cex=2)



################################################################################


# Figure 1

treeplot+geom_facet(panel="Distribution of aggression", data=fightdistributiondata_long, geom=geom_barh, mapping=aes(fill = name,y = numericspecies, x = value),stat="identity", cex=1,fill=rep(c(fm_aggression,mm_aggression,ff_aggression),117))+geom_facet(panel="Percentage fights won by females", data=fightswondata_range_min, geom=geom_segment, aes(x=min_perc_won_females,xend=min_perc_won_females+delta_long,y=Count,yend=Count+delta_lat), size=2,color=c("black"))+geom_facet(panel="Intersexual dominance", data=dominancedata, geom=geom_point, aes(x=placeholder), cex=2,color=c(dominancedata$colours),pch=15)+theme_tree()

plot(obj,lwd=5,outline=FALSE,direction="leftwards")





################################################################################
################################################################################
################################################################################

# Figure 2
# 2a) Mating system, dimorphism body size, dimorphism canine size, sex ratio, male reproductive skew

# 2b) Arboreality, sexual receptivity, concealed ovulation, testes size, reproductive synchrony

# 2c) Social system, (average female kinship), female coalitions, female philopatry

# 2d) number of males, proportion of male-male conflicts, 

# 2e) * Captivity, * environmental harshness, rainfall seasonality, rainfall unpredictability, *  seasonal breeding, * home range overlap, number of females, * female eviction, female infanticide, female canine size



# 2a) Mating system * , dimorphism body size *, dimorphism canine size *, sex ratio *, male reproductive skew *
# 2b) Arboreality *, sexual receptivity * , concealed ovulation NO, testes size NO, reproductive synchrony NO
# 2c) Social system *, female coalitions *, female philopatry *, number of males *, proportion of male-male conflicts NO
# 2d) Captivity NO, environmental harshness *, seasonal breeding *, home range overlap *, female eviction *



################################################################################
##### Figure 2 a: dominance and male mating

### top: raw data

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
summarizedtable[summarizedtable$MatingSystem =="POL", ]$MatingSystem<-"Pzogyny"
summarizedtable[summarizedtable$MatingSystem =="PRO", ]$MatingSystem<-"Promiscuity"

summarizedtable[summarizedtable$StrictFemdom ==1, ]$StrictFemdom<-"3) male dominance"
summarizedtable[summarizedtable$StrictFemdom ==2, ]$StrictFemdom<-"2) co-dominance"
summarizedtable[summarizedtable$StrictFemdom ==3, ]$StrictFemdom<-"1) female dominance"

summarizedtable<-summarizedtable[order(summarizedtable$MatingSystem,summarizedtable$StrictFemdom),]

summarizedtable$StrictFemdom<-as.factor(summarizedtable$StrictFemdom)

plot_matingsystem <-ggplot(summarizedtable, aes(x = MatingSystem, y = Observations, fill = StrictFemdom)) + 
  geom_bar(stat = "identity",fill=c(female_dominance_color,co_dominance_color,male_dominance_color,female_dominance_color,co_dominance_color,male_dominance_color,female_dominance_color,co_dominance_color,male_dominance_color,female_dominance_color,co_dominance_color,male_dominance_color))+
  theme(       axis.text.y=element_blank(),
               axis.ticks.y=element_blank(),
               axis.title.y = element_blank(),
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
        axis.text=element_text(size=10)
  )+
  scale_x_continuous(name="Male canine size relative to female canine size")+scale_fill_manual(values=c(female_dominance_color,co_dominance_color,male_dominance_color))+theme(legend.position="none")



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
               axis.text=element_text(size=10)
  )+
  scale_y_discrete(labels = c('Strict female dominance', 'Co dominance', 'Strict male dominance'))+scale_x_continuous(name="Male reproductive skew less or more than expected by chance")+scale_fill_manual(values=c(female_dominance_color,co_dominance_color,male_dominance_color),)+theme(legend.position="none")+geom_point(x=100,y=1.1,pch=24, fill=female_dominance_color, alpha=0.5,size=8, colour=female_dominance_color)+geom_point(x=100,y=1,pch=21, fill="black", alpha=0.5,size=3, colour="black") 


# combined plot
plot_grid(plot_matingsystem ,plot_bodysizedimorphism, plot_caninesizedimorphism, plot_sexratio,plot_reproductiveskew, rel_widths = c(5,3,3,3,3),nrow=1)



### bottom: model output

# mating system

speciesaverage<-combined[is.na(combined$MatSysPMK)==F,]

dat_list_strict <- list(
  R = as.integer(as.factor(1/as.integer(speciesaverage$strictfdom))),
  matingsystem = as.integer(as.factor(speciesaverage$MatSysPMK)),
  species = as.integer(as.factor(speciesaverage$corrected_species_id))
)

m_matingsystems <- ulam(
  alist(
    R ~ dordlogit( phi , cutpoints ),
    phi <-a + b[matingsystem],
    a ~ normal( 0 , 5 ),
    b[matingsystem] ~ dnorm(0,5),
    cutpoints ~ dnorm( 0 , 5 )
  ) , data=dat_list_strict , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)


matingsystems <- seq(from=1,to=4,by=1)
pdat <- data.frame(matingsystems=matingsystems)

overallphi_speciesaverages<-matrix(ncol=1, nrow=4)
for(i in 1:length(matingsystems)){
  overallphi_speciesaverages[i,]<-precis(m_matingsystems,depth=2)[1,1]+precis(m_matingsystems,depth=2)[i+1,1]*matingsystems[i]
}

overallprobs_speciesaverages_mating<-pordlogit( 1:3 , overallphi_speciesaverages , precis(m_matingsystems,depth=2)[6:7,1] )

overallprobs_speciesaverages_mating[,3]<-overallprobs_speciesaverages_mating[,3]-overallprobs_speciesaverages_mating[,2]
overallprobs_speciesaverages_mating[,2]<-overallprobs_speciesaverages_mating[,2]-overallprobs_speciesaverages_mating[,1]
overallprobs_speciesaverages_mating <-t(overallprobs_speciesaverages_mating)


overallprobs_speciesaverages_mating <-overallprobs_speciesaverages_mating[,c(1,2,4,3)]
colnames(overallprobs_speciesaverages_mating)<-c("Polygynous","Promiscous","Polyandrous", "Monogamous")


# dimorphism body size

speciesaverage<-as.data.frame(combined %>% group_by(corrected_species_id) %>% summarise(mean(SexualDimorphism_MaleWeight_over_FemaleWeight,na.rm=T),mean(as.numeric(strictfdom))))

colnames(speciesaverage)<-c("corrected_species_id","sizedimorphism","strictfdom")
speciesaverage<-speciesaverage[is.na(speciesaverage$sizedimorphism)==F,]
speciesaverage$strictfdom<-round(speciesaverage$strictfdom,0)

dat_list_strict_bodydimorphism <- list(
  R = as.integer(as.factor(1/speciesaverage$strictfdom)),
  sizedimorphism = standardize(exp(speciesaverage$sizedimorphism)),
  species = as.integer(as.factor(speciesaverage$corrected_species_id))
)

m_speciesaverage_sizedimorphism <- ulam(
  alist(
    R ~ dordlogit( phi , cutpoints ),
    phi <-a + b*sizedimorphism ,
    a ~ normal( 0 , 5 ),
    b ~ dnorm(0,5),
    cutpoints ~ dnorm( 0 , 5 )
  ) , data= dat_list_strict_bodydimorphism , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)


# standardized size dimorphism
sizedimorphism <- seq(from=-1.5,to=3.5,by=0.5)
pdat <- data.frame(sizedimorphism=sizedimorphism)

overallphi_speciesaverages<-matrix(ncol=1, nrow=11)
for(i in 1:length(sizedimorphism)){
  overallphi_speciesaverages[i,]<-precis(m_speciesaverage_sizedimorphism,depth=2)[1,1]+precis(m_speciesaverage_sizedimorphism,depth=2)[2,1]*sizedimorphism[i]
}

overallprobs_speciesaverages_sizedimorphism<-pordlogit( 1:3 , overallphi_speciesaverages , precis(m_speciesaverage_sizedimorphism,depth=2)[3:4,1] )




# dimorphism canine size

speciesaverage<-as.data.frame(combined %>% group_by(corrected_species_id) %>% summarise(mean(CanineDimorphism,na.rm=T),mean(as.numeric(strictfdom))))

colnames(speciesaverage)<-c("corrected_species_id","sizedimorphism","strictfdom")
speciesaverage<-speciesaverage[is.na(speciesaverage$sizedimorphism)==F,]
speciesaverage$strictfdom<-round(speciesaverage$strictfdom,0)

dat_list_strict_caninedimorphism <- list(
  R = as.integer(as.factor(1/speciesaverage$strictfdom)),
  sizedimorphism = standardize(exp(speciesaverage$sizedimorphism)),
  species = as.integer(as.factor(speciesaverage$corrected_species_id))
)

m_speciesaverage_caninesizedimorphism <- ulam(
  alist(
    R ~ dordlogit( phi , cutpoints ),
    phi <-a + b*sizedimorphism ,
    a ~ normal( 0 , 5 ),
    b ~ dnorm(0,5),
    cutpoints ~ dnorm( 0 , 5 )
  ) , data= dat_list_strict_caninedimorphism , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)


# standardized size dimorphism
sizedimorphism_canine <- seq(from=-1.5,to=4.5,by=0.5)
pdat <- data.frame(sizedimorphism=sizedimorphism)

overallphi_speciesaverages<-matrix(ncol=1, nrow=13)
for(i in 1:length(sizedimorphism_canine)){
  overallphi_speciesaverages[i,]<-precis(m_speciesaverage_caninesizedimorphism,depth=2)[1,1]+precis(m_speciesaverage_caninesizedimorphism,depth=2)[2,1]* sizedimorphism_canine[i]
}

overallprobs_speciesaverages_caninedimorphism<-pordlogit( 1:3 , overallphi_speciesaverages , precis(m_speciesaverage_sizedimorphism,depth=2)[3:4,1] )



# sex ratio 

speciesaverage<-combined[is.na(combined$sexratio)==F,]

dat_list_strict_sexratio <- list(
  R = as.integer(as.factor(1/as.integer(speciesaverage$strictfdom))),
  sexratio = standardize(speciesaverage$sexratio),
  species = as.integer(as.factor(speciesaverage$corrected_species_id))
)

m_speciesaverage_sexratio <- ulam(
  alist(
    R ~ dordlogit( phi , cutpoints ),
    phi <-a + b*sexratio ,
    a ~ normal( 0 , 5 ),
    b ~ dnorm(0,5),
    cutpoints ~ dnorm( 0 , 5 )
  ) , data=dat_list_strict_sexratio , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)


sexratio <- seq(from=-2.5,to=2.5,by=0.5)
pdat <- data.frame(sexratio=sexratio)

overallphi_speciesaverages<-matrix(ncol=1, nrow=11)
for(i in 1:length(sexratio)){
  overallphi_speciesaverages[i,]<-precis(m_speciesaverage_sexratio,depth=2)[1,1]+precis(m_speciesaverage_sexratio,depth=2)[2,1]*sexratio[i]
}

overallprobs_speciesaverages_sexratio<-pordlogit( 1:3 , overallphi_speciesaverages , precis(m_speciesaverage_sexratio,depth=2)[3:4,1] )



# male reproductive skew


speciesaverage<-combined[is.na(combined$M_skew_index)==F,]

dat_list_strict_skew <- list(
  R = as.integer(as.factor(1/as.integer(speciesaverage$strictfdom))),
  male_skew = standardize(speciesaverage$M_skew_index),
  species = as.integer(as.factor(speciesaverage$corrected_species_id))
)

m_speciesaverage_skew <- ulam(
  alist(
    R ~ dordlogit( phi , cutpoints ),
    phi <-a + b*male_skew ,
    a ~ normal( 0 , 5 ),
    b ~ dnorm(0,5),
    cutpoints ~ dnorm( 0 , 5 )
  ) , data=dat_list_strict_skew , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)


male_skew <- seq(from=-3.0,to=2.0,by=0.5)
pdat <- data.frame(male_skew=male_skew)

overallphi_speciesaverages<-matrix(ncol=1, nrow=11)
for(i in 1:length(male_skew)){
  overallphi_speciesaverages[i,]<-precis(m_speciesaverage_skew,depth=2)[1,1]+precis(m_speciesaverage_skew,depth=2)[2,1]*male_skew[i]
}

overallprobs_speciesaverages_skew<-pordlogit( 1:3 , overallphi_speciesaverages , precis(m_speciesaverage_skew,depth=2)[3:4,1] )






# combined plot

previouspar<-par()
op <- par(oma=c(1.5,1.5,1.5,0.2), mar=c(1.5,1.5,1.5,1.5), mfrow=c(1,5))
barplot(overallprobs_speciesaverages_mating,col=dominance_colors,axisnames = F)

plot( NULL , type="n" , xlab="Male body weight relative to female body weight" ,
      xlim=c(min(dat_list_strict_bodydimorphism $sizedimorphism,na.rm=T),max(dat_list_strict_bodydimorphism $sizedimorphism,na.rm=T)) , ylim=c(0,1) ,  yaxp=c(0,1,4) ,bty="n",xaxt="n",ylab="",yaxt="n")
axis(side=1,at=c(-1.347572,-0.768047,-0.1885223,0.3910024,0.9705271,1.550052,2.129577,2.709101,3.288626),labels=FALSE)
polygon(x=c(sizedimorphism,rev(sizedimorphism)),y=c(overallprobs_speciesaverages_sizedimorphism[,2],rev(overallprobs_speciesaverages_sizedimorphism[,1])),col=co_dominance_color,border=NA)
polygon(x=c(sizedimorphism,rev(sizedimorphism)),y=c(overallprobs_speciesaverages_sizedimorphism[,1],rep(0,11)),col=female_dominance_color,border=NA)
polygon(x=c(sizedimorphism,rev(sizedimorphism)),y=c(rep(1,11),rev(overallprobs_speciesaverages_sizedimorphism[,2])),col=male_dominance_color,border=NA)

plot( NULL , type="n" , xlab="Male canine size relative to female canine size" ,
      xlim=c(min(dat_list_strict_caninedimorphism $sizedimorphism,na.rm=T),max(dat_list_strict_caninedimorphism $sizedimorphism,na.rm=T)) , ylim=c(0,1) ,  yaxp=c(0,1,4) ,bty="n",xaxt="n",ylab="",yaxt="n")
axis(side=1,at=c( (0.75-1.678664)/0.7459911, (1.00-1.678664)/0.7459911,(2-1.678664)/0.7459911,(3-1.678664)/0.7459911 ,(4-1.678664)/0.7459911,(5-1.678664)/0.7459911)  ,labels=FALSE)
polygon(x=c(sizedimorphism_canine,rev(sizedimorphism_canine)),y=c(overallprobs_speciesaverages_caninedimorphism[,2],rev(overallprobs_speciesaverages_caninedimorphism[,1])),col=co_dominance_color,border=NA)
polygon(x=c(sizedimorphism_canine,rev(sizedimorphism_canine)),y=c(overallprobs_speciesaverages_caninedimorphism[,1],rep(0,13)),col=female_dominance_color,border=NA)
polygon(x=c(sizedimorphism_canine,rev(sizedimorphism_canine)),y=c(rep(1,13),rev(overallprobs_speciesaverages_caninedimorphism[,2])),col=male_dominance_color,border=NA)

plot( NULL , type="n" , xlab="Number of males relative to number of females" ,      xlim=c(min(dat_list_strict_sexratio$sexratio,na.rm=T),max(dat_list_strict_sexratio$sexratio,na.rm=T)) , ylim=c(0,1) ,  yaxp=c(0,1,4) ,bty="n",xaxt="n",ylab="",yaxt="n")
axis(side=1,at=c(-2.166648,-0.068047,1,2.324043),labels=FALSE)
polygon(x=c(sexratio,rev(sexratio)),y=c(overallprobs_speciesaverages_sexratio[,2],rev(overallprobs_speciesaverages_sexratio[,1])),col=co_dominance_color,border=NA)
polygon(x=c(sexratio,rev(sexratio)),y=c(overallprobs_speciesaverages_sexratio[,1],rep(0,11)),col=female_dominance_color,border=NA)
polygon(x=c(sexratio,rev(sexratio)),y=c(rep(1,11),rev(overallprobs_speciesaverages_sexratio[,2])),col=male_dominance_color,border=NA)


plot( NULL , type="n" , xlab="Male reproductive skew less or more than expected by chance" ,  xlim=c(min(dat_list_strict_skew$male_skew,na.rm=T),max(dat_list_strict_skew$male_skew,na.rm=T)) , ylim=c(0,1) ,  yaxp=c(0,1,4) ,bty="n",xaxt="n",ylab="",yaxt="n")
axis(side=1,at=c( (-2.5-0.1507726)/0.8886107, (-1.25-0.1507726)/0.8886107, (0-0.1507726)/0.8886107, (1.25-0.1507726)/0.8886107, (1.75-0.1507726)/0.8886107  ),labels=FALSE)
polygon(x=c(male_skew,rev(male_skew)),y=c(overallprobs_speciesaverages_skew[,2],rev(overallprobs_speciesaverages_skew[,1])),col=co_dominance_color,border=NA)
polygon(x=c(male_skew,rev(male_skew)),y=c(overallprobs_speciesaverages_skew[,1],rep(0,11)),col=female_dominance_color,border=NA)
polygon(x=c(male_skew,rev(male_skew)),y=c(rep(1,11),rev(overallprobs_speciesaverages_skew[,2])),col=male_dominance_color,border=NA)

par<-previouspar





################################################################################
##### Figure 2 b: dominance and female mating
# Arboreality, sexual receptivity, concealed ovulation, testes size, reproductive synchrony
### top: raw data

# arboreality

summarizedtable<-combined %>%
  group_by(Strata_Wilman,strictfdom) %>%
  summarize(Total = n())

summarizedtable<-as.data.frame(summarizedtable)

summarizedtable[9,]<-c("G",3,0)


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

plot_arboreality <-ggplot(summarizedtable, aes(x = factor(Arboreality,levels=c("Arboreal","Scansorial","Ground")), y = Observations, fill = StrictFemdom)) + 
  geom_bar(stat = "identity",fill=c(female_dominance_color,co_dominance_color,male_dominance_color,female_dominance_color,co_dominance_color,male_dominance_color,female_dominance_color,co_dominance_color,male_dominance_color))+
  theme(       axis.text.y=element_blank(),
               axis.ticks.y=element_blank(),
               axis.title.y = element_blank(),
               axis.title.x = element_blank(),
               axis.text=element_text(size=10)
  )


# concealed evolution

summarizedtable<-combined %>%
  group_by(ovulation_signs,strictfdom) %>%
  summarize(Total = n())

summarizedtable<-as.data.frame(summarizedtable)

summarizedtable<-summarizedtable[1:9,]


colnames(summarizedtable)<-c("Ovulationsigns","StrictFemdom","Observations")
summarizedtable$Observations<-as.integer(summarizedtable$Observations)


summarizedtable[summarizedtable$StrictFemdom ==1, ]$StrictFemdom<-"3) male dominance"
summarizedtable[summarizedtable$StrictFemdom ==2, ]$StrictFemdom<-"2) co-dominance"
summarizedtable[summarizedtable$StrictFemdom ==3, ]$StrictFemdom<-"1) female dominance"

summarizedtable<-summarizedtable[order(summarizedtable$Ovulationsigns,summarizedtable$StrictFemdom),]

summarizedtable$StrictFemdom<-as.factor(summarizedtable$StrictFemdom)

plot_ovulationsigns <-ggplot(summarizedtable, aes(x = factor(Ovulationsigns,levels=c("Absent","Present","Exaggerated")), y = Observations, fill = StrictFemdom)) + 
  geom_bar(stat = "identity",fill=c(col.alpha(female_dominance_color,0.4),col.alpha(co_dominance_color,0.4),col.alpha(male_dominance_color,0.4),col.alpha(female_dominance_color,0.4),col.alpha(co_dominance_color,0.4),col.alpha(male_dominance_color,0.4),col.alpha(female_dominance_color,0.4),col.alpha(co_dominance_color,0.4),col.alpha(male_dominance_color,0.4)))+
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
        axis.text=element_text(size=10)
  )+
  scale_y_discrete(labels = c('Strict female dominance', 'Co dominance', 'Strict male dominance'))+
  scale_x_continuous(name="sexual receptivity (days)")+scale_fill_manual(values=c(female_dominance_color,co_dominance_color,male_dominance_color))+theme(legend.position="none")




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
                axis.text=element_text(size=10)
  )+
  scale_y_discrete(labels = c('Strict female dominance', 'Co dominance', 'Strict male dominance'))+scale_fill_manual(values=c(col.alpha(female_dominance_color,0.4),col.alpha(co_dominance_color,0.4),col.alpha(male_dominance_color,0.4)))+theme(legend.position="none")








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
               axis.text=element_text(size=10)
  )+
  scale_y_discrete(labels = c('Strict female dominance', 'Co dominance', 'Strict male dominance'))+scale_x_continuous(name="Reproductive synchrony")+scale_fill_manual(values=c( col.alpha(female_dominance_color,0.4),col.alpha(co_dominance_color,0.4),col.alpha(male_dominance_color,0.4) ),)+theme(legend.position="none")


#  Fig 2b top combined plot 
plot_grid(plot_arboreality ,plot_ovulationsigns, plot_sexualreceptivity, plot_testesmass,plot_synchrony, rel_widths = c(5,3,3,3,3),nrow=1,scale=0.85)



### bottom: model output

# arboreality

speciesaverage<-combined[is.na(combined$Strata_Wilman)==F,]

dat_list_strict_arboreality <- list(
  R = as.integer(as.factor(1/as.integer(speciesaverage$strictfdom))),
  arboreality = as.integer(as.factor(speciesaverage$Strata_Wilman)),
  species = as.integer(as.factor(speciesaverage$corrected_species_id))
)

m_arboreality <- ulam(
  alist(
    R ~ dordlogit( phi , cutpoints ),
    phi <-a + b[arboreality],
    a ~ normal( 0 , 5 ),
    b[arboreality] ~ dnorm(0,5),
    cutpoints ~ dnorm( 0 , 5 )
  ) , data=dat_list_strict_arboreality , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)


arboreality <- seq(from=1,to=3,by=1)
pdat <- data.frame(arboreality=arboreality)

overallphi_speciesaverages<-matrix(ncol=1, nrow=3)
for(i in 1:length(arboreality)){
  overallphi_speciesaverages[i,]<-precis(m_arboreality,depth=2)[1,1]+precis(m_arboreality,depth=2)[i+1,1]*arboreality[i]
}

overallprobs_speciesaverages_arboreality<-pordlogit( 1:3 , overallphi_speciesaverages , precis(m_arboreality,depth=2)[5:6,1] )

overallprobs_speciesaverages_arboreality[,3]<-overallprobs_speciesaverages_arboreality[,3]-overallprobs_speciesaverages_arboreality[,2]
overallprobs_speciesaverages_arboreality[,2]<-overallprobs_speciesaverages_arboreality[,2]-overallprobs_speciesaverages_arboreality[,1]
overallprobs_speciesaverages_arboreality <-t(overallprobs_speciesaverages_arboreality)


overallprobs_speciesaverages_arboreality <-overallprobs_speciesaverages_arboreality[,c(1,3,2)]



# ovulationsigns

speciesaverage<-combined[is.na(combined$ovulation_signs)==F,]

dat_list_strict_ovulation_signs <- list(
  R = as.integer(as.factor(1/as.integer(speciesaverage$strictfdom))),
  ovulation_signs = as.integer(as.factor(speciesaverage$ovulation_signs)),
  species = as.integer(as.factor(speciesaverage$corrected_species_id))
)

m_ovulation_signs <- ulam(
  alist(
    R ~ dordlogit( phi , cutpoints ),
    phi <-a + b[ovulation_signs],
    a ~ normal( 0 , 5 ),
    b[ovulation_signs] ~ dnorm(0,5),
    cutpoints ~ dnorm( 0 , 5 )
  ) , data=dat_list_strict_ovulation_signs , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)


ovulation_signs <- seq(from=1,to=3,by=1)
pdat <- data.frame(ovulation_signs=ovulation_signs)

overallphi_speciesaverages<-matrix(ncol=1, nrow=3)
for(i in 1:length(ovulation_signs)){
  overallphi_speciesaverages[i,]<-precis(m_ovulation_signs,depth=2)[1,1]+precis(m_ovulation_signs,depth=2)[i+1,1]*ovulation_signs[i]
}

overallprobs_speciesaverages_ovulation_signs<-pordlogit( 1:3 , overallphi_speciesaverages , precis(m_ovulation_signs,depth=2)[5:6,1] )

overallprobs_speciesaverages_ovulation_signs[,3]<-overallprobs_speciesaverages_ovulation_signs[,3]-overallprobs_speciesaverages_ovulation_signs[,2]
overallprobs_speciesaverages_ovulation_signs[,2]<-overallprobs_speciesaverages_ovulation_signs[,2]-overallprobs_speciesaverages_ovulation_signs[,1]
overallprobs_speciesaverages_ovulation_signs <-t(overallprobs_speciesaverages_ovulation_signs)


overallprobs_speciesaverages_ovulation_signs <-overallprobs_speciesaverages_ovulation_signs[,c(1,3,2)]





# sexual receptivity

sexualreceptivity_hours_data<-combined[is.na(combined$sexualreceptivity_hours)==F,]
sexualreceptivity_hours_data$strictfdom<-as.integer(sexualreceptivity_hours_data$strictfdom)

dat_list_strict_sexualreceptivity_hours <- list(
  R = as.integer(as.factor(1/sexualreceptivity_hours_data$strictfdom)),
  sexualreceptivity_hours = standardize(exp(sexualreceptivity_hours_data$sexualreceptivity_hours)),
  species = as.integer(as.factor(sexualreceptivity_hours_data$corrected_species_id))
)

m_speciesaverage_sexualreceptivity_hours <- ulam(
  alist(
    R ~ dordlogit( phi , cutpoints ),
    phi <-a + b*sexualreceptivity_hours ,
    a ~ normal( 0 , 5 ),
    b ~ dnorm(0,5),
    cutpoints ~ dnorm( 0 , 5 )
  ) , data= dat_list_strict_sexualreceptivity_hours , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)



sexualreceptivity_hours <- seq(from=-1.5,to=3.5,by=0.5)
pdat <- data.frame(sexualreceptivity_hours=sexualreceptivity_hours)

overallphi_speciesaverages<-matrix(ncol=1, nrow=11)
for(i in 1:length(sexualreceptivity_hours)){
  overallphi_speciesaverages[i,]<-precis(m_speciesaverage_sexualreceptivity_hours,depth=2)[1,1]+precis(m_speciesaverage_sexualreceptivity_hours,depth=2)[2,1]*sexualreceptivity_hours[i]
}

overallprobs_speciesaverages_sexualreceptivity_hours<-pordlogit( 1:3 , overallphi_speciesaverages , precis(m_speciesaverage_sexualreceptivity_hours,depth=2)[3:4,1] )




# testes size

testes_data<-combined[is.na(combined$relative_testes_mass)==F,]
testes_data$strictfdom<-as.integer(testes_data$strictfdom)

dat_list_strict_relative_testes_mass<- list(
  R = as.integer(as.factor(1/testes_data$strictfdom)),
  relative_testes_mass = standardize((testes_data$relative_testes_mass)),
  species = as.integer(as.factor(testes_data$corrected_species_id))
)

m_speciesaverage_relative_testes_mass <- ulam(
  alist(
    R ~ dordlogit( phi , cutpoints ),
    phi <-a + b*relative_testes_mass ,
    a ~ normal( 0 , 5 ),
    b ~ dnorm(0,5),
    cutpoints ~ dnorm( 0 , 5 )
  ) , data= dat_list_strict_relative_testes_mass , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)


relative_testes_mass <- seq(from=-2.5,to=1.5,by=0.5)
pdat <- data.frame(relative_testes_mass=relative_testes_mass)

overallphi_speciesaverages<-matrix(ncol=1, nrow=length(relative_testes_mass))
for(i in 1:length(relative_testes_mass)){
  overallphi_speciesaverages[i,]<-precis(m_speciesaverage_relative_testes_mass,depth=2)[1,1]+precis(m_speciesaverage_relative_testes_mass,depth=2)[2,1]* relative_testes_mass[i]
}

overallprobs_speciesaverages_relative_testes_mass<-pordlogit( 1:3 , overallphi_speciesaverages , precis(m_speciesaverage_relative_testes_mass,depth=2)[3:4,1] )



# synchrony

speciesaverage<-combined[is.na(combined$Synchrony)==F,]

dat_list_strict_Synchrony <- list(
  R = as.integer(as.factor(1/as.integer(speciesaverage$strictfdom))),
  Synchrony = standardize(speciesaverage$Synchrony),
  species = as.integer(as.factor(speciesaverage$corrected_species_id))
)

m_speciesaverage_Synchrony <- ulam(
  alist(
    R ~ dordlogit( phi , cutpoints ),
    phi <-a + b*Synchrony ,
    a ~ normal( 0 , 5 ),
    b ~ dnorm(0,5),
    cutpoints ~ dnorm( 0 , 5 )
  ) , data=dat_list_strict_Synchrony , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)


Synchrony <- seq(from=-1,to=2.0,by=0.5)
pdat <- data.frame(Synchrony=Synchrony)

overallphi_speciesaverages<-matrix(ncol=1, nrow=length(Synchrony))
for(i in 1:length(Synchrony)){
  overallphi_speciesaverages[i,]<-precis(m_speciesaverage_Synchrony,depth=2)[1,1]+precis(m_speciesaverage_Synchrony,depth=2)[2,1]*Synchrony[i]
}

overallprobs_speciesaverages_synchrony<-pordlogit( 1:3 , overallphi_speciesaverages , precis(m_speciesaverage_Synchrony,depth=2)[3:4,1] )





# combined plot

plot_grid(plot_arboreality ,plot_ovulationsigns, plot_sexualreceptivity, plot_testesmass,plot_synchrony, rel_widths = c(5,3,3,3,3),nrow=1,scale=0.85)

# Fig 2b bottom
previouspar<-par()
op <- par(oma=c(1.5,1.5,1.5,0.2), mar=c(1.5,1.5,1.5,1.5), mfrow=c(1,5))
barplot(overallprobs_speciesaverages_arboreality,col=dominance_colors,axisnames = F)

barplot(overallprobs_speciesaverages_ovulation_signs,col=c( col.alpha(dominance_colors[1],0.4), col.alpha(dominance_colors[2],0.4),col.alpha(dominance_colors[3],0.4)),axisnames = F,yaxt="n")

plot( NULL , type="n" , xlab="sexual receptivity" ,
      xlim=c(min(dat_list_strict_sexualreceptivity_hours $sexualreceptivity_hours,na.rm=T),max(dat_list_strict_sexualreceptivity_hours $sexualreceptivity_hours,na.rm=T)) , ylim=c(0,1) ,  yaxp=c(0,1,4) ,bty="n",xaxt="n",ylab="",yaxt="n")
polygon(x=c(sexualreceptivity_hours,rev(sexualreceptivity_hours)),y=c(overallprobs_speciesaverages_sexualreceptivity_hours[,2],rev(overallprobs_speciesaverages_sexualreceptivity_hours[,1])),col=co_dominance_color,border=NA)
polygon(x=c(sexualreceptivity_hours,rev(sexualreceptivity_hours)),y=c(overallprobs_speciesaverages_sexualreceptivity_hours[,1],rep(0,length(sexualreceptivity_hours))),col=female_dominance_color,border=NA)
polygon(x=c(sexualreceptivity_hours,rev(sexualreceptivity_hours)),y=c(rep(1,length(sexualreceptivity_hours)),rev(overallprobs_speciesaverages_sexualreceptivity_hours[,2])),col=male_dominance_color,border=NA)

plot( NULL , type="n" , xlab="testes size" ,
      xlim=c(min(dat_list_strict_relative_testes_mass$relative_testes_mass,na.rm=T),max(dat_list_strict_relative_testes_mass$relative_testes_mass,na.rm=T)) , ylim=c(0,1) ,  yaxp=c(0,1,4) ,bty="n",xaxt="n",ylab="",yaxt="n")
polygon(x=c(relative_testes_mass,rev(relative_testes_mass)),y=c(overallprobs_speciesaverages_relative_testes_mass[,2],rev(overallprobs_speciesaverages_relative_testes_mass[,1])),col=col.alpha(co_dominance_color,0.4),border=NA)
polygon(x=c(relative_testes_mass,rev(relative_testes_mass)),y=c(overallprobs_speciesaverages_relative_testes_mass[,1],rep(0,length(relative_testes_mass))),col=col.alpha(female_dominance_color,0.4),border=NA)
polygon(x=c(relative_testes_mass,rev(relative_testes_mass)),y=c(rep(1,length(relative_testes_mass)),rev(overallprobs_speciesaverages_relative_testes_mass[,2])),col=col.alpha(male_dominance_color,0.4),border=NA)


plot( NULL , type="n" , xlab="synchrony" ,
      xlim=c(min(dat_list_strict_Synchrony$Synchrony,na.rm=T),max(dat_list_strict_Synchrony$Synchrony,na.rm=T)) , ylim=c(0,1) ,  yaxp=c(0,1,4) ,bty="n",xaxt="n",ylab="",yaxt="n")
polygon(x=c(Synchrony,rev(Synchrony)),y=c(overallprobs_speciesaverages_synchrony[,2],rev(overallprobs_speciesaverages_synchrony[,1])),col=col.alpha(co_dominance_color,0.4),border=NA)
polygon(x=c(Synchrony,rev(Synchrony)),y=c(overallprobs_speciesaverages_synchrony[,1],rep(0,length(Synchrony))),col=col.alpha(female_dominance_color,0.4),border=NA)
polygon(x=c(Synchrony,rev(Synchrony)),y=c(rep(1,length(Synchrony)),rev(overallprobs_speciesaverages_synchrony[,2])),col=col.alpha(male_dominance_color,0.4),border=NA)

par<-previouspar



################################################################################
##### Figure 2 c: dominance and female mating
# Social system, female coalitions, female philopatry, number of males, proportion of male-male conflicts 
### top: raw data

# social system

summarizedtable<-combined %>%
  group_by(SocOrgPMK,strictfdom) %>%
  summarize(Total = n())

summarizedtable<-as.data.frame(summarizedtable)

summarizedtable[8,]<-c("P",1,0)
summarizedtable[9,]<-c("S",1,0)


colnames(summarizedtable)<-c("SocOrgPMK","StrictFemdom","Observations")
summarizedtable$Observations<-as.integer(summarizedtable$Observations)
summarizedtable[summarizedtable$SocOrgPMK =="G", ]$SocOrgPMK<-"Group"
summarizedtable[summarizedtable$SocOrgPMK =="S", ]$SocOrgPMK<-"Solitary"
summarizedtable[summarizedtable$SocOrgPMK =="P", ]$SocOrgPMK<-"Pair"


summarizedtable[summarizedtable$StrictFemdom ==1, ]$StrictFemdom<-"3) male dominance"
summarizedtable[summarizedtable$StrictFemdom ==2, ]$StrictFemdom<-"2) co-dominance"
summarizedtable[summarizedtable$StrictFemdom ==3, ]$StrictFemdom<-"1) female dominance"

summarizedtable<-summarizedtable[order(summarizedtable$SocOrgPMK,summarizedtable$StrictFemdom),]

summarizedtable$StrictFemdom<-as.factor(summarizedtable$StrictFemdom)

plot_socialsystem <-ggplot(summarizedtable, aes(x = factor(SocOrgPMK,levels=c("Group","Solitary","Pair")), y = Observations, fill = StrictFemdom)) + 
  geom_bar(stat = "identity",fill=c(female_dominance_color,co_dominance_color,male_dominance_color,female_dominance_color,co_dominance_color,male_dominance_color,female_dominance_color,co_dominance_color,male_dominance_color))+
  theme(       axis.text.y=element_blank(),
               axis.ticks.y=element_blank(),
               axis.title.y = element_blank(),
               axis.title.x = element_blank(),
               axis.text=element_text(size=10)
  )


# female coalitions

summarizedtable<-combined %>%
  group_by(jointaggression_females,strictfdom) %>%
  summarize(Total = n())

summarizedtable<-as.data.frame(summarizedtable)

summarizedtable<-summarizedtable[1:6,]

colnames(summarizedtable)<-c("jointaggression_females","StrictFemdom","Observations")
summarizedtable$Observations<-as.integer(summarizedtable$Observations)


summarizedtable[summarizedtable$StrictFemdom ==1, ]$StrictFemdom<-"3) male dominance"
summarizedtable[summarizedtable$StrictFemdom ==2, ]$StrictFemdom<-"2) co-dominance"
summarizedtable[summarizedtable$StrictFemdom ==3, ]$StrictFemdom<-"1) female dominance"

summarizedtable<-summarizedtable[order(summarizedtable$jointaggression_females,summarizedtable$StrictFemdom),]

summarizedtable$StrictFemdom<-as.factor(summarizedtable$StrictFemdom)

plot_jointaggression_females <-ggplot(summarizedtable, aes(x = factor(jointaggression_females,levels=c("Yes","No")), y = Observations, fill = StrictFemdom)) + 
  geom_bar(stat = "identity",fill=c(female_dominance_color,co_dominance_color,male_dominance_color,female_dominance_color,co_dominance_color,male_dominance_color))+
  theme(       axis.text.y=element_blank(),
               axis.ticks.y=element_blank(),
               axis.title.y = element_blank(),
               axis.title.x = element_blank(),
               axis.text=element_text(size=10)
  )




# sex bias in dispersal

summarizedtable<-combined %>%
  group_by(sexbias_dispersal,strictfdom) %>%
  summarize(Total = n())

summarizedtable<-as.data.frame(summarizedtable)

summarizedtable<-summarizedtable[1:8,]
summarizedtable[9,]<-c("Female",3,0)

colnames(summarizedtable)<-c("sexbias_dispersal","StrictFemdom","Observations")
summarizedtable$Observations<-as.integer(summarizedtable$Observations)


summarizedtable[summarizedtable$StrictFemdom ==1, ]$StrictFemdom<-"3) male dominance"
summarizedtable[summarizedtable$StrictFemdom ==2, ]$StrictFemdom<-"2) co-dominance"
summarizedtable[summarizedtable$StrictFemdom ==3, ]$StrictFemdom<-"1) female dominance"

summarizedtable<-summarizedtable[order(summarizedtable$sexbias_dispersal,summarizedtable$StrictFemdom),]

summarizedtable$StrictFemdom<-as.factor(summarizedtable$StrictFemdom)

plot_sexbias_dispersal <-ggplot(summarizedtable, aes(x = factor(sexbias_dispersal,levels=c("Male","Both","Female")), y = Observations, fill = StrictFemdom)) + 
  geom_bar(stat = "identity",fill=c(female_dominance_color,co_dominance_color,male_dominance_color,female_dominance_color,co_dominance_color,male_dominance_color,female_dominance_color,co_dominance_color,male_dominance_color))+
  theme(       axis.text.y=element_blank(),
               axis.ticks.y=element_blank(),
               axis.title.y = element_blank(),
               axis.title.x = element_blank(),
               axis.text=element_text(size=10)
  )


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
  scale_x_continuous(name="number of males")+scale_fill_manual(values=c(female_dominance_color,co_dominance_color,male_dominance_color))+theme(legend.position="none")




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
  scale_y_discrete(labels = c('Strict female dominance', 'Co dominance', 'Strict male dominance'))+scale_fill_manual(values=c(col.alpha(female_dominance_color,0.4),col.alpha(co_dominance_color,0.4),col.alpha(male_dominance_color,0.4) ),)+theme(legend.position="none")




# combined plot
plot_grid(plot_socialsystem ,plot_jointaggression_females, plot_sexbias_dispersal, plot_males,plot_percaggress_mm, rel_widths = c(3,3,3,3,3),nrow=1,scale=0.85)



### bottom: model output

# SocOrgPMK

speciesaverage<-combined[is.na(combined$SocOrgPMK)==F,]

dat_list_strict_SocOrgPMK <- list(
  R = as.integer(as.factor(1/as.integer(speciesaverage$strictfdom))),
  SocOrgPMK = as.integer(as.factor(speciesaverage$SocOrgPMK)),
  species = as.integer(as.factor(speciesaverage$corrected_species_id))
)

m_SocOrgPMK <- ulam(
  alist(
    R ~ dordlogit( phi , cutpoints ),
    phi <-a + b[SocOrgPMK],
    a ~ normal( 0 , 5 ),
    b[SocOrgPMK] ~ dnorm(0,5),
    cutpoints ~ dnorm( 0 , 5 )
  ) , data=dat_list_strict_SocOrgPMK , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)


SocOrgPMK <- seq(from=1,to=3,by=1)
pdat <- data.frame(SocOrgPMK=SocOrgPMK)

overallphi_speciesaverages<-matrix(ncol=1, nrow=3)
for(i in 1:length(SocOrgPMK)){
  overallphi_speciesaverages[i,]<-precis(m_SocOrgPMK,depth=2)[1,1]+precis(m_SocOrgPMK,depth=2)[i+1,1]*SocOrgPMK[i]
}

overallprobs_speciesaverages_SocOrgPMK<-pordlogit( 1:3 , overallphi_speciesaverages , precis(m_SocOrgPMK,depth=2)[5:6,1] )

overallprobs_speciesaverages_SocOrgPMK[,3]<-overallprobs_speciesaverages_SocOrgPMK[,3]-overallprobs_speciesaverages_SocOrgPMK[,2]
overallprobs_speciesaverages_SocOrgPMK[,2]<-overallprobs_speciesaverages_SocOrgPMK[,2]-overallprobs_speciesaverages_SocOrgPMK[,1]
overallprobs_speciesaverages_SocOrgPMK <-t(overallprobs_speciesaverages_SocOrgPMK)


overallprobs_speciesaverages_SocOrgPMK <-overallprobs_speciesaverages_SocOrgPMK[,c(1,3,2)]



# joint aggression

speciesaverage<-combined[is.na(combined$jointaggression_females)==F,]

dat_list_jointaggression_females <- list(
  R = as.integer(as.factor(1/as.integer(speciesaverage$strictfdom))),
  jointaggression_females = as.integer(as.factor(speciesaverage$jointaggression_females)),
  species = as.integer(as.factor(speciesaverage$corrected_species_id))
)

m_jointaggression_females <- ulam(
  alist(
    R ~ dordlogit( phi , cutpoints ),
    phi <-a + b[jointaggression_females],
    a ~ normal( 0 , 5 ),
    b[jointaggression_females] ~ dnorm(0,5),
    cutpoints ~ dnorm( 0 , 5 )
  ) , data=dat_list_jointaggression_females , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)


jointaggression_females <- seq(from=1,to=2,by=1)
pdat <- data.frame(jointaggression_females=jointaggression_females)

overallphi_speciesaverages<-matrix(ncol=1, nrow=2)
for(i in 1:length(jointaggression_females)){
  overallphi_speciesaverages[i,]<-precis(m_jointaggression_females,depth=2)[1,1]+precis(m_jointaggression_females,depth=2)[i+1,1]*jointaggression_females[i]
}

overallprobs_speciesaverages_jointaggression_females<-pordlogit( 1:3 , overallphi_speciesaverages , precis(m_jointaggression_females,depth=2)[4:5,1] )

overallprobs_speciesaverages_jointaggression_females[,3]<-overallprobs_speciesaverages_jointaggression_females[,3]-overallprobs_speciesaverages_jointaggression_females[,2]
overallprobs_speciesaverages_jointaggression_females[,2]<-overallprobs_speciesaverages_jointaggression_females[,2]-overallprobs_speciesaverages_jointaggression_females[,1]
overallprobs_speciesaverages_jointaggression_females <-t(overallprobs_speciesaverages_jointaggression_females)








# sex bias dispersal

speciesaverage<-combined[is.na(combined$sexbias_dispersal)==F,]

dat_list_sexbias_dispersal <- list(
  R = as.integer(as.factor(1/as.integer(speciesaverage$strictfdom))),
  sexbias_dispersal = as.integer(as.factor(speciesaverage$sexbias_dispersal)),
  species = as.integer(as.factor(speciesaverage$corrected_species_id))
)

m_jointaggression_sexbiasdispersal <- ulam(
  alist(
    R ~ dordlogit( phi , cutpoints ),
    phi <-a + b[sexbias_dispersal],
    a ~ normal( 0 , 5 ),
    b[sexbias_dispersal] ~ dnorm(0,5),
    cutpoints ~ dnorm( 0 , 5 )
  ) , data=dat_list_sexbias_dispersal , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)


sexbias_dispersal <- seq(from=1,to=3,by=1)
pdat <- data.frame(sexbias_dispersal=sexbias_dispersal)

overallphi_speciesaverages<-matrix(ncol=1, nrow=3)
for(i in 1:length(sexbias_dispersal)){
  overallphi_speciesaverages[i,]<-precis(m_jointaggression_sexbiasdispersal,depth=2)[1,1]+precis(m_jointaggression_sexbiasdispersal,depth=2)[i+1,1]*sexbias_dispersal[i]
}

overallprobs_speciesaverages_sexbias_dispersal<-pordlogit( 1:3 , overallphi_speciesaverages , precis(m_jointaggression_sexbiasdispersal,depth=2)[5:6,1] )

overallprobs_speciesaverages_sexbias_dispersal[,3]<-overallprobs_speciesaverages_sexbias_dispersal[,3]-overallprobs_speciesaverages_sexbias_dispersal[,2]
overallprobs_speciesaverages_sexbias_dispersal[,2]<-overallprobs_speciesaverages_sexbias_dispersal[,2]-overallprobs_speciesaverages_sexbias_dispersal[,1]
overallprobs_speciesaverages_sexbias_dispersal <-t(overallprobs_speciesaverages_sexbias_dispersal)

overallprobs_speciesaverages_sexbias_dispersal<-overallprobs_speciesaverages_sexbias_dispersal[,c(3,1,2)]


# number males

males_data<-combined[is.na(combined$males)==F,]
males_data$strictfdom<-as.integer(males_data$strictfdom)

dat_list_strict_males<- list(
  R = as.integer(as.factor(1/males_data$strictfdom)),
  males = standardize((males_data$males)),
  species = as.integer(as.factor(males_data$corrected_species_id))
)

m_speciesaverage_relative_males <- ulam(
  alist(
    R ~ dordlogit( phi , cutpoints ),
    phi <-a + b*males ,
    a ~ normal( 0 , 5 ),
    b ~ dnorm(0,5),
    cutpoints ~ dnorm( 0 , 5 )
  ) , data= dat_list_strict_males , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)


males <- seq(from=min(dat_list_strict_males$males),to=max(dat_list_strict_males$males),length.out=9)
pdat <- data.frame(males=males)

overallphi_speciesaverages<-matrix(ncol=1, nrow=length(males))
for(i in 1:length(males)){
  overallphi_speciesaverages[i,]<-precis(m_speciesaverage_relative_males,depth=2)[1,1]+precis(m_speciesaverage_relative_males,depth=2)[2,1]* males[i]
}

overallprobs_speciesaverages_relative_males<-pordlogit( 1:3 , overallphi_speciesaverages , precis(m_speciesaverage_relative_males,depth=2)[3:4,1] )



# percent aggression male-male

speciesaverage<-combined[is.na(combined$perc_aggression_mm)==F,]

dat_list_strict_perc_aggression_mm <- list(
  R = as.integer(as.factor(1/as.integer(speciesaverage$strictfdom))),
  perc_aggression_mm = standardize(speciesaverage$perc_aggression_mm),
  species = as.integer(as.factor(speciesaverage$corrected_species_id))
)

m_speciesaverage_perc_aggression_mm <- ulam(
  alist(
    R ~ dordlogit( phi , cutpoints ),
    phi <-a + b*perc_aggression_mm ,
    a ~ normal( 0 , 5 ),
    b ~ dnorm(0,5),
    cutpoints ~ dnorm( 0 , 5 )
  ) , data=dat_list_strict_perc_aggression_mm , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)


perc_aggression_mm <- seq(from=min(dat_list_strict_perc_aggression_mm$perc_aggression_mm),to=max((dat_list_strict_perc_aggression_mm$perc_aggression_mm)),length.out=9)
pdat <- data.frame(perc_aggression_mm=perc_aggression_mm)

overallphi_speciesaverages<-matrix(ncol=1, nrow=length(perc_aggression_mm))
for(i in 1:length(perc_aggression_mm)){
  overallphi_speciesaverages[i,]<-precis(m_speciesaverage_perc_aggression_mm,depth=2)[1,1]+precis(m_speciesaverage_perc_aggression_mm,depth=2)[2,1]*perc_aggression_mm[i]
}

overallprobs_speciesaverages_perc_aggression_mm<-pordlogit( 1:3 , overallphi_speciesaverages , precis(m_speciesaverage_perc_aggression_mm,depth=2)[3:4,1] )





# combined plot

previouspar<-par()
op <- par(oma=c(1.5,1.5,1.5,0.2), mar=c(1.5,1.5,1.5,1.5), mfrow=c(1,5))
barplot(overallprobs_speciesaverages_SocOrgPMK,col=dominance_colors,axisnames = F)


barplot(overallprobs_speciesaverages_jointaggression_females,col=c(col.alpha(dominance_colors[1],alpha=0.4),col.alpha(dominance_colors[2],alpha=0.4),col.alpha(dominance_colors[3],alpha=0.4)),axisnames = F,yaxt="n")

barplot(overallprobs_speciesaverages_sexbias_dispersal,col=dominance_colors,axisnames = F,yaxt="n")




plot( NULL , type="n" , xlab="males" ,
      xlim=c(min(dat_list_strict_males $males,na.rm=T),max(dat_list_strict_males $males,na.rm=T)) , ylim=c(0,1) ,  yaxp=c(0,1,4) ,bty="n",xaxt="n",ylab="",yaxt="n")

polygon(x=c(males,rev(males)),y=c(overallprobs_speciesaverages_relative_males[,2],rev(overallprobs_speciesaverages_relative_males[,1])),col=co_dominance_color,border=NA)
polygon(x=c(males,rev(males)),y=c(overallprobs_speciesaverages_relative_males[,1],rep(0,length(males))),col=female_dominance_color,border=NA)
polygon(x=c(males,rev(males)),y=c(rep(1,length(males)),rev(overallprobs_speciesaverages_relative_males[,2])),col=male_dominance_color,border=NA)




plot( NULL , type="n" , xlab="perc_aggression_mm" ,
      xlim=c(min(dat_list_strict_perc_aggression_mm$perc_aggression_mm,na.rm=T),max(dat_list_strict_perc_aggression_mm$perc_aggression_mm,na.rm=T)) , ylim=c(0,1) ,  yaxp=c(0,1,4) ,bty="n",xaxt="n",ylab="",yaxt="n")

polygon(x=c(perc_aggression_mm,rev(perc_aggression_mm)),y=c(overallprobs_speciesaverages_perc_aggression_mm[,2],rev(overallprobs_speciesaverages_perc_aggression_mm[,1])),col=col.alpha(co_dominance_color,0.4),border=NA)
polygon(x=c(perc_aggression_mm,rev(perc_aggression_mm)),y=c(overallprobs_speciesaverages_perc_aggression_mm[,1],rep(0,length(perc_aggression_mm))),col=col.alpha(female_dominance_color,0.4),border=NA)
polygon(x=c(perc_aggression_mm,rev(perc_aggression_mm)),y=c(rep(1,length(perc_aggression_mm)),rev(overallprobs_speciesaverages_perc_aggression_mm[,2])),col=col.alpha(male_dominance_color,0.4),border=NA)




par<-previouspar







################################################################################
##### Figure 2 d: dominance and female competition
# Captivity, environmental harshness, seasonal breeding, home range overlap, female eviction
### top: raw data

# Captivity

summarizedtable<-combined %>%
  group_by(origin,strictfdom) %>%
  summarize(Total = n())

summarizedtable<-as.data.frame(summarizedtable)

summarizedtable<-summarizedtable[c(1,2,3,5,6,7),]


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


# Female eviction

summarizedtable<-combined %>%
  group_by(female_evictions,strictfdom) %>%
  summarize(Total = n())

summarizedtable<-as.data.frame(summarizedtable)

summarizedtable<-summarizedtable[c(1:5),]
summarizedtable[6,]<-c("No",3,0)

colnames(summarizedtable)<-c("female_evictions","StrictFemdom","Observations")
summarizedtable$Observations<-as.integer(summarizedtable$Observations)



summarizedtable[summarizedtable$StrictFemdom ==1, ]$StrictFemdom<-"3) male dominance"
summarizedtable[summarizedtable$StrictFemdom ==2, ]$StrictFemdom<-"2) co-dominance"
summarizedtable[summarizedtable$StrictFemdom ==3, ]$StrictFemdom<-"1) female dominance"

summarizedtable<-summarizedtable[order(summarizedtable$female_evictions,summarizedtable$StrictFemdom),]

summarizedtable$StrictFemdom<-as.factor(summarizedtable$StrictFemdom)

plot_femaleevictions <-ggplot(summarizedtable, aes(x = factor(female_evictions,levels=c("Yes","No")), y = Observations, fill = StrictFemdom)) + 
  geom_bar(stat = "identity",fill=c(female_dominance_color,co_dominance_color,male_dominance_color,female_dominance_color,co_dominance_color,male_dominance_color))+
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
        axis.text=element_text(size=10)
  )+
  scale_y_discrete(labels = c('Strict female dominance', 'Co dominance', 'Strict male dominance'))+
  scale_x_continuous(name="sexual receptivity (days)")+scale_fill_manual(values=c(female_dominance_color,co_dominance_color,male_dominance_color))+theme(legend.position="none")


# seasonal breeding
df_seasonalbreeding <- combined[ complete.cases(combined$strictfdom,combined$r_seasonality_value),]

df_seasonalbreeding[df_seasonalbreeding$strictfdom==1,]$strictfdom<-"c) Strict male dominance"
df_seasonalbreeding[df_seasonalbreeding$strictfdom==2,]$strictfdom<-"b) Co dominance"
df_seasonalbreeding[df_seasonalbreeding$strictfdom==3,]$strictfdom<-"a) Strict female dominance"

df = data.frame(
  strictfdom=df_seasonalbreeding$strictfdom,
  r_seasonality_value=(df_seasonalbreeding$r_seasonality_value)
)

plot_seasonalbreeding<-ggplot(df)+aes(y=strictfdom,x=r_seasonality_value,fill=strictfdom)+stat_halfeye(aes(thickness = after_stat(pdf*n)))+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y = element_blank(),
        axis.text=element_text(size=10)
  )+
  scale_y_discrete(labels = c('Strict female dominance', 'Co dominance', 'Strict male dominance'))+
  scale_x_continuous(name="seasonality in breeding")+scale_fill_manual(values=c(female_dominance_color,co_dominance_color,male_dominance_color))+theme(legend.position="none")




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
                axis.text=element_text(size=10)
  )+
  scale_y_discrete(labels = c('Strict female dominance', 'Co dominance', 'Strict male dominance'))+scale_fill_manual(values=c(female_dominance_color,co_dominance_color,male_dominance_color),)+theme(legend.position="none")




# combined plot
plot_grid(plot_origin ,plot_femaleevictions, plot_env_harshness, plot_seasonalbreeding,plot_homerange_overlap, rel_widths = c(3,3,3,3,3),nrow=1,scale=0.85)



### bottom: model output

# origin

speciesaverage<-combined[is.na(combined$origin)==F,]
speciesaverage<-speciesaverage[speciesaverage$origin !="provisioned",]

dat_list_strict_origin <- list(
  R = (as.factor(1/as.integer(speciesaverage$strictfdom))),
  origin = as.integer(as.factor(speciesaverage$origin)),
  species = as.integer(as.factor(speciesaverage$corrected_species_id))
)

m_origin <- ulam(
  alist(
    R ~ dordlogit( phi , cutpoints ),
    phi <-a + b[origin],
    a ~ normal( 0 , 5 ),
    b[origin] ~ dnorm(0,5),
    cutpoints ~ dnorm( 0 , 5 )
  ) , data=dat_list_strict_origin , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)


origin <- seq(from=1,to=2,by=1)
pdat <- data.frame(origin=origin)

overallphi_speciesaverages<-matrix(ncol=1, nrow=2)
for(i in 1:length(origin)){
  overallphi_speciesaverages[i,]<-precis(m_origin,depth=2)[1,1]+precis(m_origin,depth=2)[i+1,1]*origin[i]
}

overallprobs_speciesaverages_origin<-pordlogit( 1:3 , overallphi_speciesaverages , precis(m_origin,depth=2)[4:5,1] )

overallprobs_speciesaverages_origin[,3]<-overallprobs_speciesaverages_origin[,3]-overallprobs_speciesaverages_origin[,2]
overallprobs_speciesaverages_origin[,2]<-overallprobs_speciesaverages_origin[,2]-overallprobs_speciesaverages_origin[,1]
overallprobs_speciesaverages_origin <-t(overallprobs_speciesaverages_origin)




# female evictions

speciesaverage<-combined[is.na(combined$female_evictions)==F,]

dat_list_strict_female_evictions <- list(
  R = as.integer(as.factor(1/as.integer(speciesaverage$strictfdom))),
  female_evictions = as.integer(as.factor(speciesaverage$female_evictions)),
  species = as.integer(as.factor(speciesaverage$corrected_species_id))
)

m_female_evictions <- ulam(
  alist(
    R ~ dordlogit( phi , cutpoints ),
    phi <-a + b[female_evictions],
    a ~ normal( 0 , 5 ),
    b[female_evictions] ~ dnorm(0,5),
    cutpoints ~ dnorm( 0 , 5 )
  ) , data=dat_list_strict_female_evictions , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)


female_evictions <- seq(from=1,to=2,by=1)
pdat <- data.frame(female_evictions=female_evictions)

overallphi_speciesaverages<-matrix(ncol=1, nrow=2)
for(i in 1:length(female_evictions)){
  overallphi_speciesaverages[i,]<-precis(m_female_evictions,depth=2)[1,1]+precis(m_female_evictions,depth=2)[i+1,1]*female_evictions[i]
}

overallprobs_speciesaverages_female_evictions<-pordlogit( 1:3 , overallphi_speciesaverages , precis(m_female_evictions,depth=2)[4:5,1] )

overallprobs_speciesaverages_female_evictions[,3]<-overallprobs_speciesaverages_female_evictions[,3]-overallprobs_speciesaverages_female_evictions[,2]
overallprobs_speciesaverages_female_evictions[,2]<-overallprobs_speciesaverages_female_evictions[,2]-overallprobs_speciesaverages_female_evictions[,1]
overallprobs_speciesaverages_female_evictions <-t(overallprobs_speciesaverages_female_evictions)







# environmental harshness

env_harshness_data<-combined[is.na(combined$env_harshness)==F,]
env_harshness_data$strictfdom<-as.integer(env_harshness_data$strictfdom)

dat_list_strict_env_harshness <- list(
  R = as.integer(as.factor(1/env_harshness_data$strictfdom)),
  env_harshness = standardize((env_harshness_data$env_harshness)),
  species = as.integer(as.factor(env_harshness_data$corrected_species_id))
)

m_speciesaverage_env_harshness <- ulam(
  alist(
    R ~ dordlogit( phi , cutpoints ),
    phi <-a + b*env_harshness ,
    a ~ normal( 0 , 5 ),
    b ~ dnorm(0,5),
    cutpoints ~ dnorm( 0 , 5 )
  ) , data= dat_list_strict_env_harshness , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)



env_harshness <- seq(from=min(dat_list_strict_env_harshness$env_harshness),to=max(dat_list_strict_env_harshness$env_harshness),length.out=9)
pdat <- data.frame(env_harshness=env_harshness)

overallphi_speciesaverages<-matrix(ncol=1, nrow=9)
for(i in 1:length(env_harshness)){
  overallphi_speciesaverages[i,]<-precis(m_speciesaverage_env_harshness,depth=2)[1,1]+precis(m_speciesaverage_env_harshness,depth=2)[2,1]*env_harshness[i]
}

overallprobs_speciesaverages_env_harshness<-pordlogit( 1:3 , overallphi_speciesaverages , precis(m_speciesaverage_env_harshness,depth=2)[3:4,1] )




# seasonal breeding

r_seasonality_value_data<-combined[is.na(combined$r_seasonality_value)==F,]
r_seasonality_value_data$strictfdom<-as.integer(r_seasonality_value_data$strictfdom)

dat_list_strict_r_seasonality_value<- list(
  R = as.integer(as.factor(1/r_seasonality_value_data$strictfdom)),
  r_seasonality_value = standardize((r_seasonality_value_data$r_seasonality_value)),
  species = as.integer(as.factor(r_seasonality_value_data$corrected_species_id))
)

m_speciesaverage_r_seasonality_value <- ulam(
  alist(
    R ~ dordlogit( phi , cutpoints ),
    phi <-a + b*r_seasonality_value ,
    a ~ normal( 0 , 5 ),
    b ~ dnorm(0,5),
    cutpoints ~ dnorm( 0 , 5 )
  ) , data= dat_list_strict_r_seasonality_value, chains=4 , cores=4 ,cmdstan=T, messages=FALSE)


r_seasonality_value <- seq(from=min(dat_list_strict_r_seasonality_value$r_seasonality_value),to=max(dat_list_strict_r_seasonality_value$r_seasonality_value),length.out=9)
pdat <- data.frame(r_seasonality_value=r_seasonality_value)

overallphi_speciesaverages<-matrix(ncol=1, nrow=length(r_seasonality_value))
for(i in 1:length(r_seasonality_value)){
  overallphi_speciesaverages[i,]<-precis(m_speciesaverage_r_seasonality_value,depth=2)[1,1]+precis(m_speciesaverage_r_seasonality_value,depth=2)[2,1]* r_seasonality_value[i]
}

overallprobs_speciesaverages_r_seasonality_value<-pordlogit( 1:3 , overallphi_speciesaverages , precis(m_speciesaverage_r_seasonality_value,depth=2)[3:4,1] )



# home range overlap

speciesaverage<-combined[is.na(combined$homerange_overlap)==F,]

dat_list_strict_homerange_overlap <- list(
  R = as.integer(as.factor(1/as.integer(speciesaverage$strictfdom))),
  homerange_overlap = standardize(speciesaverage$homerange_overlap),
  species = as.integer(as.factor(speciesaverage$corrected_species_id))
)

m_speciesaverage_homerange_overlap <- ulam(
  alist(
    R ~ dordlogit( phi , cutpoints ),
    phi <-a + b*homerange_overlap ,
    a ~ normal( 0 , 5 ),
    b ~ dnorm(0,5),
    cutpoints ~ dnorm( 0 , 5 )
  ) , data=dat_list_strict_homerange_overlap , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)


homerange_overlap <- seq(from=-1,to=2.0,by=0.5)
pdat <- data.frame(homerange_overlap=homerange_overlap)

overallphi_speciesaverages<-matrix(ncol=1, nrow=length(homerange_overlap))
for(i in 1:length(homerange_overlap)){
  overallphi_speciesaverages[i,]<-precis(m_speciesaverage_homerange_overlap,depth=2)[1,1]+precis(m_speciesaverage_homerange_overlap,depth=2)[2,1]*homerange_overlap[i]
}

overallprobs_speciesaverages_homerange_overlap<-pordlogit( 1:3 , overallphi_speciesaverages , precis(m_speciesaverage_homerange_overlap,depth=2)[3:4,1] )





# combined plot
previouspar<-par()
op <- par(oma=c(1.5,1.5,1.5,0.2), mar=c(1.5,1.5,1.5,1.5), mfrow=c(1,5))
barplot(overallprobs_speciesaverages_origin,col=c(col.alpha(dominance_colors[1],0.4),col.alpha(dominance_colors[2],0.4),col.alpha(dominance_colors[3],0.4) ),axisnames = F)

barplot(overallprobs_speciesaverages_female_evictions,col=dominance_colors,axisnames = F,yaxt="n")

plot( NULL , type="n" , xlab="environmental harshness" ,
      xlim=c(min(dat_list_strict_env_harshness$env_harshness,na.rm=T),max(dat_list_strict_env_harshness$env_harshness,na.rm=T)) , ylim=c(0,1) ,  yaxp=c(0,1,4) ,bty="n",xaxt="n",ylab="",yaxt="n")
polygon(x=c(env_harshness,rev(env_harshness)),y=c(overallprobs_speciesaverages_env_harshness[,2],rev(overallprobs_speciesaverages_env_harshness[,1])),col=co_dominance_color,border=NA)
polygon(x=c(env_harshness,rev(env_harshness)),y=c(overallprobs_speciesaverages_env_harshness[,1],rep(0,length(env_harshness))),col=female_dominance_color,border=NA)
polygon(x=c(env_harshness,rev(env_harshness)),y=c(rep(1,length(env_harshness)),rev(overallprobs_speciesaverages_env_harshness[,2])),col=male_dominance_color,border=NA)


plot( NULL , type="n" , xlab="seasonal breeding" ,
      xlim=c(min(dat_list_strict_r_seasonality_value$r_seasonality_value,na.rm=T),max(dat_list_strict_r_seasonality_value$r_seasonality_value,na.rm=T)) , ylim=c(0,1) ,  yaxp=c(0,1,4) ,bty="n",xaxt="n",ylab="",yaxt="n")
polygon(x=c(r_seasonality_value,rev(r_seasonality_value)),y=c(overallprobs_speciesaverages_r_seasonality_value[,2],rev(overallprobs_speciesaverages_r_seasonality_value[,1])),col=co_dominance_color,border=NA)
polygon(x=c(r_seasonality_value,rev(r_seasonality_value)),y=c(overallprobs_speciesaverages_r_seasonality_value[,1],rep(0,length(r_seasonality_value))),col=female_dominance_color,border=NA)
polygon(x=c(r_seasonality_value,rev(r_seasonality_value)),y=c(rep(1,length(r_seasonality_value)),rev(overallprobs_speciesaverages_r_seasonality_value[,2])),col=male_dominance_color,border=NA)


plot( NULL , type="n" , xlab="home range overlap" ,    xlim=c(min(dat_list_strict_homerange_overlap$homerange_overlap,na.rm=T),max(dat_list_strict_homerange_overlap$homerange_overlap,na.rm=T)) , ylim=c(0,1) ,  yaxp=c(0,1,4) ,bty="n",xaxt="n",ylab="",yaxt="n")
polygon(x=c(homerange_overlap,rev(homerange_overlap)),y=c(overallprobs_speciesaverages_homerange_overlap[,2],rev(overallprobs_speciesaverages_homerange_overlap[,1])),col=co_dominance_color,border=NA)
polygon(x=c(homerange_overlap,rev(homerange_overlap)),y=c(overallprobs_speciesaverages_homerange_overlap[,1],rep(0,length(homerange_overlap))),col=female_dominance_color,border=NA)
polygon(x=c(homerange_overlap,rev(homerange_overlap)),y=c(rep(1,length(homerange_overlap)),rev(overallprobs_speciesaverages_homerange_overlap[,2])),col=male_dominance_color,border=NA)


par<-previouspar


































################################################################################

################################################################################
# Previous code
################################################################################

################################################################################

# Figure 2: dominance and mating system

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

summarizedtable<-as.data.frame(summarizedtable)

summarizedtable<-rbind(summarizedtable,c(1,1,0))
summarizedtable<-rbind(summarizedtable,c(2,1,0))

summarizedtable[summarizedtable$MatingSystem ==1, ]$MatingSystem<-"Monogamy"
summarizedtable[summarizedtable$MatingSystem ==2, ]$MatingSystem<-"Polyandry"
summarizedtable[summarizedtable$MatingSystem ==3, ]$MatingSystem<-"Polygyny"
summarizedtable[summarizedtable$MatingSystem ==4, ]$MatingSystem<-"Promiscuity"

summarizedtable[summarizedtable$StrictFemdom ==1, ]$StrictFemdom<-"c) male dominance"
summarizedtable[summarizedtable$StrictFemdom ==2, ]$StrictFemdom<-"b) co-dominance"
summarizedtable[summarizedtable$StrictFemdom ==3, ]$StrictFemdom<-"a) female dominance"

ggplot (as.data.frame(summarizedtable), aes(x=MatingSystem, y=NumberOfSpecies, fill=StrictFemdom)) + geom_bar (stat="identity", position ="dodge")+scale_fill_manual(values=dominance_colors)+ labs(
  x = "Mating System", y = "Number of Species", fill="Intersexual dominance")+theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+theme(legend.position = "none")


polyandryplot<-ggplot (as.data.frame(summarizedtable[summarizedtable$MatingSystem %in% "Polyandry",]), aes(x=NumberOfSpecies, y=StrictFemdom, fill=StrictFemdom)) + geom_barh (stat="identity")+labs(
  x = "Number of species", y = "", fill=c("strict male dominance","no sex bias in dominance","strict female dominance"),title="Polyandry")+theme(axis.text=element_text(size=12),axis.title=element_text(size=14),plot.title=element_text(size=20,face="bold",hjust=0.5),legend.position = "none")+scale_fill_manual(values=(dominance_colors))+xlim(0,35)

monogamyplot<-ggplot (as.data.frame(summarizedtable[summarizedtable$MatingSystem %in% "Monogamy",]), aes(x=NumberOfSpecies, y=StrictFemdom, fill=StrictFemdom)) + geom_barh (stat="identity")+labs(
  x = "Number of species", y = "", fill=c("strict male dominance","no sex bias in dominance","strict female dominance"),title="Monogamy")+theme(axis.text.x=element_text(size=12),axis.title.x=element_text(size=14),axis.title.y=element_blank(),axis.text.y=element_blank(),plot.title=element_text(size=20,face="bold",hjust=0.5),legend.position = "none")+scale_fill_manual(values=(dominance_colors))+xlim(0,35)

promiscuityplot<-ggplot (as.data.frame(summarizedtable[summarizedtable$MatingSystem %in% "Promiscuity",]), aes(x=NumberOfSpecies, y=StrictFemdom, fill=StrictFemdom)) + geom_barh (stat="identity")+labs(
  x = "Number of species", y = "", fill=c("strict male dominance","no sex bias in dominance","strict female dominance"),title="Promiscuity")+theme(axis.text.x=element_text(size=12),axis.title.x=element_text(size=14),axis.title.y=element_blank(),axis.text.y=element_blank(),plot.title=element_text(size=20,face="bold",hjust=0.5),legend.position = "none")+scale_fill_manual(values=(dominance_colors))+xlim(0,35)


polygynyplot<-ggplot(as.data.frame(summarizedtable[summarizedtable$MatingSystem %in% "Polygyny",]), aes(x=NumberOfSpecies, y=StrictFemdom, fill=StrictFemdom)) + geom_barh (stat="identity")+labs(
  x = "Number of species", y = "", fill=c("strict male dominance","no sex bias in dominance","strict female dominance"),title="Polygyny")+theme(axis.text.x=element_text(size=12),axis.title.x=element_text(size=14),axis.title.y=element_blank(),axis.text.y=element_blank(),plot.title=element_text(size=20,face="bold",hjust=0.5),legend.position = "none")+scale_fill_manual(values=(dominance_colors))+xlim(0,35)

plot_grid(polyandryplot, monogamyplot, promiscuityplot,polygynyplot,rel_widths = c(5,3,3,3),nrow=1)


summarizedtable[summarizedtable$MatingSystem=="Polygyny",]$MatingSystem<-"Pzogyny"

ggplot(summarizedtable, aes(x = MatingSystem, y = NumberOfSpecies, fill = StrictFemdom)) + 
  geom_bar(stat = "identity",fill=c(co_dominance_color,female_dominance_color,co_dominance_color,female_dominance_color,male_dominance_color,co_dominance_color,female_dominance_color,male_dominance_color,co_dominance_color,female_dominance_color,co_dominance_color,female_dominance_color))+
  theme(       axis.text.y=element_blank(),
               axis.ticks.y=element_blank(),
               axis.title.y = element_blank(),
               axis.title.x = element_blank(),
               axis.text=element_text(size=10)
  )




################################################################################################
################################################################################################

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
                axis.text=element_text(size=10)
  )+
  scale_y_discrete(labels = c('Strict female dominance', 'Co dominance', 'Strict male dominance'))+scale_x_continuous(name="Sex ratio in social groups",lim=c(0.09,0.76),breaks=c(0.1, 0.33,0.5, 0.66),labels=c("1♂ / 10♀♀","1♂ / 2♀♀","1♂ / 1♀","2♂♂ / 1♀"))+scale_fill_manual(values=c(female_dominance_color,co_dominance_color,male_dominance_color),)+theme(legend.position="none")



# testes sizes

df_testessize <- combined[ complete.cases(combined$strictfdom,combined$relative_testes_mass),]

df_testessize[df_testessize$strictfdom==1,]$strictfdom<-"c) Strict male dominance"
df_testessize[df_testessize$strictfdom==2,]$strictfdom<-"b) Co dominance"
df_testessize[df_testessize$strictfdom==3,]$strictfdom<-"a) Strict female dominance"

df = data.frame(
  strictfdom=df_testessize$strictfdom,
  testes_size=df_testessize$relative_testes_mass
)

plot_testessize<-ggplot(df)+aes(y=strictfdom,x=testes_size,fill=strictfdom)+stat_halfeye(aes(thickness = after_stat(pdf*n)))+
  theme(        axis.text.y=element_blank(),
                axis.ticks.y=element_blank(),
                axis.title.y = element_blank(),
                axis.text=element_text(size=10)
  )+
  scale_y_discrete(labels = c('Strict female dominance', 'Co dominance', 'Strict male dominance'))+scale_x_continuous(name="Testes sizes relative to body size",lim=c(-2,2))+scale_fill_manual(values=c(female_dominance_color,co_dominance_color,male_dominance_color),)+theme(legend.position="none")



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
                axis.text=element_text(size=10)
  )+
  scale_y_discrete(labels = c('Strict female dominance', 'Co dominance', 'Strict male dominance'))+scale_x_continuous(name="Male reproductive skew less or more than expected by chance")+scale_fill_manual(values=c(female_dominance_color,co_dominance_color,male_dominance_color),)+theme(legend.position="none")+geom_point(x=100,y=1.1,pch=24, fill=female_dominance_color, alpha=0.5,size=8, colour=female_dominance_color)+geom_point(x=100,y=1,pch=21, fill="black", alpha=0.5,size=3, colour="black") 



################################################################################
# Figure 2 a: dominance and male mating system



df_bodysizedimorphism <- combined[ complete.cases(combined$strictfdom,combined$SexualDimorphism_MaleWeight_over_FemaleWeight),]

df_bodysizedimorphism[df_bodysizedimorphism$strictfdom==1,]$strictfdom<-"c) Strict male dominance"
df_bodysizedimorphism[df_bodysizedimorphism$strictfdom==2,]$strictfdom<-"b) Co dominance"
df_bodysizedimorphism[df_bodysizedimorphism$strictfdom==3,]$strictfdom<-"a) Strict female dominance"

df = data.frame(
  strictfdom=df_bodysizedimorphism$strictfdom,
  SexualDimorphism_MaleWeight_over_FemaleWeight=exp(df_bodysizedimorphism$SexualDimorphism_MaleWeight_over_FemaleWeight)
)

plot_bodysizedimorphism<-ggplot(df)+aes(y=strictfdom,x=SexualDimorphism_MaleWeight_over_FemaleWeight,fill=strictfdom)+stat_halfeye(aes(thickness = after_stat(pdf*n)))+
  theme(        axis.ticks.y=element_blank(),
                axis.title.y = element_blank(),
                axis.text=element_text(size=10)
  )+
  scale_y_discrete(labels = c('Strict female dominance', 'Co dominance', 'Strict male dominance'))+
  scale_x_continuous(name="Male body size relative to female body size")+scale_fill_manual(values=c(female_dominance_color,co_dominance_color,male_dominance_color))+theme(legend.position="none")




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
        axis.text=element_text(size=10)
  )+
  scale_x_continuous(name="Male canine size relative to female canine size")+scale_fill_manual(values=c(female_dominance_color,co_dominance_color,male_dominance_color))+theme(legend.position="none")



#### Make a combined plot:

plot_grid(plot_bodysizedimorphism, plot_caninesizedimorphism, plot_testessize, plot_sexratio,plot_reproductiveskew, rel_widths = c(5,3,3,3,3),nrow=1)







#################################################
#################################################
# Model output figures:


# Sexual size dimorphism
speciesaverage<-as.data.frame(combined %>% group_by(corrected_species_id) %>% summarise(mean(SexualDimorphism_MaleWeight_over_FemaleWeight,na.rm=T),mean(as.numeric(strictfdom))))

colnames(speciesaverage)<-c("corrected_species_id","sizedimorphism","strictfdom")
speciesaverage<-speciesaverage[is.na(speciesaverage$sizedimorphism)==F,]
speciesaverage$strictfdom<-round(speciesaverage$strictfdom,0)

dat_list_strict <- list(
  R = as.integer(as.factor(1/speciesaverage$strictfdom)),
  sizedimorphism = standardize(exp(speciesaverage$sizedimorphism)),
  species = as.integer(as.factor(speciesaverage$corrected_species_id))
)

m_speciesaverage_sizedimorphism <- ulam(
  alist(
    R ~ dordlogit( phi , cutpoints ),
    phi <-a + b*sizedimorphism ,
    a ~ normal( 0 , 5 ),
    b ~ dnorm(0,5),
    cutpoints ~ dnorm( 0 , 5 )
  ) , data=dat_list_strict , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)

precis(m_speciesaverage_sizedimorphism,depth=2)


# standardized size dimorphism
sizedimorphism <- seq(from=-1.5,to=3.5,by=0.5)
pdat <- data.frame(sizedimorphism=sizedimorphism)

overallphi_speciesaverages<-matrix(ncol=1, nrow=11)
for(i in 1:length(sizedimorphism)){
  overallphi_speciesaverages[i,]<-precis(m_speciesaverage_sizedimorphism,depth=2)[1,1]+precis(m_speciesaverage_sizedimorphism,depth=2)[2,1]*sizedimorphism[i]
}

overallprobs_speciesaverages<-pordlogit( 1:3 , overallphi_speciesaverages , precis(m_speciesaverage_sizedimorphism,depth=2)[3:4,1] )

#(0.8-mean(exp(speciesaverage$sizedimorphism)))/sd(exp(speciesaverage$sizedimorphism))


plot( NULL , type="n" , xlab="Male body weight relative to female body weight" ,
      xlim=c(min(dat_list_strict$sizedimorphism,na.rm=T),max(dat_list_strict$sizedimorphism,na.rm=T)) , ylim=c(0,1) ,  yaxp=c(0,1,4) ,bty="n",xaxt="n",ylab="")
mtext("Cumulative probability for",side=2,at=0.5,line=3)
mtext("the three dominance systems",side=2,at=0.5,line=2)
axis(side=1,at=c(-1.347572,-0.768047,-0.1885223,0.3910024,0.9705271,1.550052,2.129577,2.709101,3.288626),labels=FALSE)

mtext("0.8",side=1,at=-1.347572,line=1)
mtext("1.0",side=1,at=-0.768047,line=1)
mtext("1.2",side=1,at=-0.1885223,line=1)
mtext("1.4",side=1,at=0.3910024,line=1)
mtext("1.6",side=1,at=0.9705271,line=1)
mtext("1.8",side=1,at=1.550052,line=1)
mtext("2.0",side=1,at=2.129577,line=1)
mtext("2.2",side=1,at=2.709101,line=1)
mtext("2.4",side=1,at=3.288626,line=1)


polygon(x=c(sizedimorphism,rev(sizedimorphism)),y=c(overallprobs_speciesaverages[,2],rev(overallprobs_speciesaverages[,1])),col=co_dominance_color,border=NA)

polygon(x=c(sizedimorphism,rev(sizedimorphism)),y=c(overallprobs_speciesaverages[,1],rep(0,11)),col=female_dominance_color,border=NA)

polygon(x=c(sizedimorphism,rev(sizedimorphism)),y=c(rep(1,11),rev(overallprobs_speciesaverages[,2])),col=male_dominance_color,border=NA)

mtext("Male dominance",side=1,line=-5,at=1.75)
mtext("Co dominance",side=1,line=-14,at=0)
mtext("Female dominance",side=1,line=-28,at=-0.9)


hist(exp(combined$CanineDimorphism),breaks=40)

barplot(as.matrix(table(combined$strictfdom)/sum(table(combined$strictfdom))),col=dominance_colors)





# Canine size dimorphism
speciesaverage<-as.data.frame(combined %>% group_by(corrected_species_id) %>% summarise(mean(CanineDimorphism,na.rm=T),mean(as.numeric(strictfdom))))

colnames(speciesaverage)<-c("corrected_species_id","sizedimorphism","strictfdom")
speciesaverage<-speciesaverage[is.na(speciesaverage$sizedimorphism)==F,]
speciesaverage$strictfdom<-round(speciesaverage$strictfdom,0)

dat_list_strict <- list(
  R = as.integer(as.factor(speciesaverage$strictfdom)),
  sizedimorphism = standardize(exp(speciesaverage$sizedimorphism)),
  species = as.integer(as.factor(speciesaverage$corrected_species_id))
)

m_speciesaverage_sizedimorphism <- ulam(
  alist(
    R ~ dordlogit( phi , cutpoints ),
    phi <-a + b*sizedimorphism ,
    a ~ normal( 0 , 5 ),
    b ~ dnorm(0,5),
    cutpoints ~ dnorm( 0 , 5 )
  ) , data=dat_list_strict , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)

precis(m_speciesaverage_sizedimorphism,depth=2)


# standardized size dimorphism
sizedimorphism <- seq(from=-1.5,to=4.5,by=0.5)
pdat <- data.frame(sizedimorphism=sizedimorphism)

overallphi_speciesaverages<-matrix(ncol=1, nrow=13)
for(i in 1:length(sizedimorphism)){
  overallphi_speciesaverages[i,]<-precis(m_speciesaverage_sizedimorphism,depth=2)[1,1]+precis(m_speciesaverage_sizedimorphism,depth=2)[2,1]*sizedimorphism[i]
}

overallprobs_speciesaverages<-pordlogit( 1:3 , overallphi_speciesaverages , precis(m_speciesaverage_sizedimorphism,depth=2)[3:4,1] )

(0.9-mean(exp(speciesaverage$sizedimorphism)))/sd(exp(speciesaverage$sizedimorphism))


plot( NULL , type="n" , xlab="Male canine size relative to female canine size" ,
      xlim=c(min(dat_list_strict$sizedimorphism,na.rm=T),max(dat_list_strict$sizedimorphism,na.rm=T)) , ylim=c(0,1) ,  yaxp=c(0,1,4) ,bty="n",xaxt="n",ylab="")
mtext("Cumulative probability for",side=2,at=0.5,line=3)
mtext("the three dominance systems",side=2,at=0.5,line=2)
axis(side=1,at=c((0.8-mean(exp(speciesaverage$sizedimorphism)))/sd(exp(speciesaverage$sizedimorphism)),(1.0-mean(exp(speciesaverage$sizedimorphism)))/sd(exp(speciesaverage$sizedimorphism)),(1.5-mean(exp(speciesaverage$sizedimorphism)))/sd(exp(speciesaverage$sizedimorphism)),(2.0-mean(exp(speciesaverage$sizedimorphism)))/sd(exp(speciesaverage$sizedimorphism)),(2.5-mean(exp(speciesaverage$sizedimorphism)))/sd(exp(speciesaverage$sizedimorphism)),(3.0-mean(exp(speciesaverage$sizedimorphism)))/sd(exp(speciesaverage$sizedimorphism)),(3.5-mean(exp(speciesaverage$sizedimorphism)))/sd(exp(speciesaverage$sizedimorphism)),(4.0-mean(exp(speciesaverage$sizedimorphism)))/sd(exp(speciesaverage$sizedimorphism)),(4.5-mean(exp(speciesaverage$sizedimorphism)))/sd(exp(speciesaverage$sizedimorphism)),(5.0-mean(exp(speciesaverage$sizedimorphism)))/sd(exp(speciesaverage$sizedimorphism))),labels=FALSE)

mtext("0.8",side=1,at=(0.8-mean(exp(speciesaverage$sizedimorphism)))/sd(exp(speciesaverage$sizedimorphism)),line=1)
mtext("1.0",side=1,at=(1.0-mean(exp(speciesaverage$sizedimorphism)))/sd(exp(speciesaverage$sizedimorphism)),line=1)
mtext("1.5",side=1,at=(1.5-mean(exp(speciesaverage$sizedimorphism)))/sd(exp(speciesaverage$sizedimorphism)),line=1)
mtext("2.0",side=1,at=(2.0-mean(exp(speciesaverage$sizedimorphism)))/sd(exp(speciesaverage$sizedimorphism)),line=1)
mtext("2.5",side=1,at=(2.5-mean(exp(speciesaverage$sizedimorphism)))/sd(exp(speciesaverage$sizedimorphism)),line=1)
mtext("3.0",side=1,at=(3.0-mean(exp(speciesaverage$sizedimorphism)))/sd(exp(speciesaverage$sizedimorphism)),line=1)
mtext("3.5",side=1,at=(3.5-mean(exp(speciesaverage$sizedimorphism)))/sd(exp(speciesaverage$sizedimorphism)),line=1)
mtext("4.0",side=1,at=(4.0-mean(exp(speciesaverage$sizedimorphism)))/sd(exp(speciesaverage$sizedimorphism)),line=1)
mtext("4.5",side=1,at=(4.5-mean(exp(speciesaverage$sizedimorphism)))/sd(exp(speciesaverage$sizedimorphism)),line=1)
mtext("5.0",side=1,at=(5.0-mean(exp(speciesaverage$sizedimorphism)))/sd(exp(speciesaverage$sizedimorphism)),line=1)

polygon(x=c(sizedimorphism,rev(sizedimorphism)),y=c(overallprobs_speciesaverages[,2],rev(overallprobs_speciesaverages[,1])),col=co_dominance_color,border=NA)

polygon(x=c(sizedimorphism,rev(sizedimorphism)),y=c(overallprobs_speciesaverages[,1],rep(0,13)),col=male_dominance_color,border=NA)

polygon(x=c(sizedimorphism,rev(sizedimorphism)),y=c(rep(1,13),rev(overallprobs_speciesaverages[,2])),col=female_dominance_color,border=NA)

mtext("Male dominance",side=1,line=-5,at=1.75)
mtext("Co dominance",side=1,line=-14,at=0)
mtext("Female dominance",side=1,line=-28,at=-0.9)


hist(exp(combined$SexualDimorphism_MaleWeight_over_FemaleWeight),breaks=40)

barplot(as.matrix(table(combined$strictfdom)/sum(table(combined$strictfdom))),col=dominance_colors)




# sex ratio

speciesaverage<-combined[is.na(combined$sexratio)==F,]

dat_list_strict <- list(
  R = as.integer(as.factor(speciesaverage$strictfdom)),
  sexratio = standardize(speciesaverage$sexratio),
  species = as.integer(as.factor(speciesaverage$corrected_species_id))
)

m_speciesaverage_testessize <- ulam(
  alist(
    R ~ dordlogit( phi , cutpoints ),
    phi <-a + b*sexratio ,
    a ~ normal( 0 , 5 ),
    b ~ dnorm(0,5),
    cutpoints ~ dnorm( 0 , 5 )
  ) , data=dat_list_strict , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)

precis(m_speciesaverage_testessize,depth=2)


# standardized sexratio
sexratio <- seq(from=-2.5,to=2.5,by=0.5)
pdat <- data.frame(sexratio=sexratio)

overallphi_speciesaverages<-matrix(ncol=1, nrow=11)
for(i in 1:length(sexratio)){
  overallphi_speciesaverages[i,]<-precis(m_speciesaverage_testessize,depth=2)[1,1]+precis(m_speciesaverage_testessize,depth=2)[2,1]*sexratio[i]
}

overallprobs_speciesaverages<-pordlogit( 1:3 , overallphi_speciesaverages , precis(m_speciesaverage_testessize,depth=2)[3:4,1] )

#(0.8-mean(exp(speciesaverage$sizedimorphism)))/sd(exp(speciesaverage$sizedimorphism))


plot( NULL , type="n" , xlab="Number of males relative to number of females" ,
      xlim=c(min(dat_list_strict$sexratio,na.rm=T),max(dat_list_strict$sexratio,na.rm=T)) , ylim=c(0,1) ,  yaxp=c(0,1,4) ,bty="n",xaxt="n",ylab="")
mtext("Cumulative probability for",side=2,at=0.5,line=3)
mtext("the three dominance systems",side=2,at=0.5,line=2)
axis(side=1,at=c(-2.166648,-0.068047,1,2.324043),labels=FALSE)

mtext("1 ♂ / 10 ♀♀",side=1,at=-2.166648,line=1)
mtext("1 ♂ / 2 ♀♀",side=1,at=-0.068047,line=1)
mtext("1 ♂ / 1 ♀",side=1,at=1,line=1)
mtext("2 ♂♂ / 1 ♀",side=1,at=2.324043,line=1)




polygon(x=c(sexratio,rev(sexratio)),y=c(overallprobs_speciesaverages[,2],rev(overallprobs_speciesaverages[,1])),col=co_dominance_color,border=NA)

polygon(x=c(sexratio,rev(sexratio)),y=c(overallprobs_speciesaverages[,1],rep(0,11)),col=male_dominance_color,border=NA)

polygon(x=c(sexratio,rev(sexratio)),y=c(rep(1,11),rev(overallprobs_speciesaverages[,2])),col=female_dominance_color,border=NA)

mtext("Male dominance",side=1,line=-5,at=1.75)
mtext("Co dominance",side=1,line=-14,at=0)
mtext("Female dominance",side=1,line=-28,at=-0.9)


hist(exp(combined$CanineDimorphism),breaks=40)

barplot(as.matrix(table(combined$strictfdom)/sum(table(combined$strictfdom))),col=dominance_colors)


# male reproductive skew

speciesaverage<-combined[is.na(combined$M_skew_index)==F,]

dat_list_strict <- list(
  R = as.integer(as.factor(speciesaverage$strictfdom)),
  male_skew = standardize(speciesaverage$M_skew_index),
  species = as.integer(as.factor(speciesaverage$corrected_species_id))
)

m_speciesaverage_testessize <- ulam(
  alist(
    R ~ dordlogit( phi , cutpoints ),
    phi <-a + b*male_skew ,
    a ~ normal( 0 , 5 ),
    b ~ dnorm(0,5),
    cutpoints ~ dnorm( 0 , 5 )
  ) , data=dat_list_strict , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)

precis(m_speciesaverage_testessize,depth=2)


# standardized testes size
male_skew <- seq(from=-3.0,to=2.0,by=0.5)
pdat <- data.frame(male_skew=male_skew)

overallphi_speciesaverages<-matrix(ncol=1, nrow=11)
for(i in 1:length(male_skew)){
  overallphi_speciesaverages[i,]<-precis(m_speciesaverage_testessize,depth=2)[1,1]+precis(m_speciesaverage_testessize,depth=2)[2,1]*male_skew[i]
}

overallprobs_speciesaverages<-pordlogit( 1:3 , overallphi_speciesaverages , precis(m_speciesaverage_testessize,depth=2)[3:4,1] )

#(0.8-mean(exp(speciesaverage$sizedimorphism)))/sd(exp(speciesaverage$sizedimorphism))


plot( NULL , type="n" , xlab="Male reproductive less or more than expected by chance" ,
      xlim=c(min(dat_list_strict$male_skew,na.rm=T),max(dat_list_strict$male_skew,na.rm=T)) , ylim=c(0,1) ,  yaxp=c(0,1,4) ,bty="n",xaxt="n",ylab="")
mtext("Cumulative probability for",side=2,at=0.5,line=3)
mtext("the three dominance systems",side=2,at=0.5,line=2)
axis(side=1,at=c(-3,-1.75,-0.5,0.75,2),labels=FALSE)

mtext("-2",side=1,at=-3,line=1)
mtext("-1",side=1,at=-1.75,line=1)
mtext("0",side=1,at=-0.5,line=1)
mtext("1",side=1,at=0.75,line=1)
mtext("2",side=1,at=2,line=1)



polygon(x=c(male_skew,rev(male_skew)),y=c(overallprobs_speciesaverages[,2],rev(overallprobs_speciesaverages[,1])),col=co_dominance_color,border=NA)

polygon(x=c(male_skew,rev(male_skew)),y=c(overallprobs_speciesaverages[,1],rep(0,11)),col=male_dominance_color,border=NA)

polygon(x=c(male_skew,rev(male_skew)),y=c(rep(1,11),rev(overallprobs_speciesaverages[,2])),col=female_dominance_color,border=NA)

mtext("Male dominance",side=1,line=-5,at=1.75)
mtext("Co dominance",side=1,line=-14,at=0)
mtext("Female dominance",side=1,line=-28,at=-0.9)


hist(exp(combined$CanineDimorphism),breaks=40)

barplot(as.matrix(table(combined$strictfdom)/sum(table(combined$strictfdom))),col=dominance_colors)









# Relative testes size

speciesaverage<-combined[is.na(combined$relative_testes_mass)==F,]

dat_list_strict <- list(
  R = as.integer(as.factor(speciesaverage$strictfdom)),
  testessize = standardize(speciesaverage$relative_testes_mass),
  species = as.integer(as.factor(speciesaverage$corrected_species_id))
)

m_speciesaverage_testessize <- ulam(
  alist(
    R ~ dordlogit( phi , cutpoints ),
    phi <-a + b*testessize ,
    a ~ normal( 0 , 5 ),
    b ~ dnorm(0,5),
    cutpoints ~ dnorm( 0 , 5 )
  ) , data=dat_list_strict , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)

precis(m_speciesaverage_testessize,depth=2)


# standardized testes size
testessize <- seq(from=-2.5,to=1.5,by=0.5)
pdat <- data.frame(testessize=testessize)

overallphi_speciesaverages<-matrix(ncol=1, nrow=9)
for(i in 1:length(testessize)){
  overallphi_speciesaverages[i,]<-precis(m_speciesaverage_testessize,depth=2)[1,1]+precis(m_speciesaverage_testessize,depth=2)[2,1]*testessize[i]
}

overallprobs_speciesaverages<-pordlogit( 1:3 , overallphi_speciesaverages , precis(m_speciesaverage_testessize,depth=2)[3:4,1] )

#(0.8-mean(exp(speciesaverage$sizedimorphism)))/sd(exp(speciesaverage$sizedimorphism))


plot( NULL , type="n" , xlab="Residual testis size relative to body size" ,
      xlim=c(min(dat_list_strict$testessize,na.rm=T),max(dat_list_strict$testessize,na.rm=T)) , ylim=c(0,1) ,  yaxp=c(0,1,4) ,bty="n",xaxt="n",ylab="")
mtext("Cumulative probability for",side=2,at=0.5,line=3)
mtext("the three dominance systems",side=2,at=0.5,line=2)
axis(side=1,at=c(-2.5,-1.25,0,1.25),labels=FALSE)

mtext("-2",side=1,at=-2.5,line=1)
mtext("-1",side=1,at=-1.25,line=1)
mtext("0",side=1,at=0,line=1)
mtext("1",side=1,at=1.25,line=1)



polygon(x=c(testessize,rev(testessize)),y=c(overallprobs_speciesaverages[,2],rev(overallprobs_speciesaverages[,1])),col=co_dominance_color,border=NA)

polygon(x=c(testessize,rev(testessize)),y=c(overallprobs_speciesaverages[,1],rep(0,9)),col=male_dominance_color,border=NA)

polygon(x=c(testessize,rev(testessize)),y=c(rep(1,9),rev(overallprobs_speciesaverages[,2])),col=female_dominance_color,border=NA)

mtext("Male dominance",side=1,line=-5,at=1.75)
mtext("Co dominance",side=1,line=-14,at=0)
mtext("Female dominance",side=1,line=-28,at=-0.9)


hist(exp(combined$CanineDimorphism),breaks=40)

barplot(as.matrix(table(combined$strictfdom)/sum(table(combined$strictfdom))),col=dominance_colors)





# sexual receptivity

speciesaverage<-combined[is.na(combined$sexualreceptivity_hours)==F,]

dat_list_strict <- list(
  R = as.integer(as.factor(speciesaverage$strictfdom)),
  sexualreceptivity_hours = standardize(exp(speciesaverage$sexualreceptivity_hours)),
  species = as.integer(as.factor(speciesaverage$corrected_species_id))
)

m_speciesaverage_testessize <- ulam(
  alist(
    R ~ dordlogit( phi , cutpoints ),
    phi <-a + b*sexualreceptivity_hours ,
    a ~ normal( 0 , 5 ),
    b ~ dnorm(0,5),
    cutpoints ~ dnorm( 0 , 5 )
  ) , data=dat_list_strict , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)

precis(m_speciesaverage_testessize,depth=2)


sexualreceptivity_hours <- seq(from=-1.5,to=3.5,by=0.5)
pdat <- data.frame(sexualreceptivity_hours=sexualreceptivity_hours)

overallphi_speciesaverages<-matrix(ncol=1, nrow=11)
for(i in 1:length(sexualreceptivity_hours)){
  overallphi_speciesaverages[i,]<-precis(m_speciesaverage_testessize,depth=2)[1,1]+precis(m_speciesaverage_testessize,depth=2)[2,1]*sexualreceptivity_hours[i]
}

overallprobs_speciesaverages<-pordlogit( 1:3 , overallphi_speciesaverages , precis(m_speciesaverage_testessize,depth=2)[3:4,1] )

#(0.8-mean(exp(speciesaverage$sizedimorphism)))/sd(exp(speciesaverage$sizedimorphism))


plot( NULL , type="n" , xlab="Sexual receptivity duration" ,
      xlim=c(min(dat_list_strict$sexualreceptivity_hours,na.rm=T),max(dat_list_strict$sexualreceptivity_hours,na.rm=T)) , ylim=c(0,1) ,  yaxp=c(0,1,4) ,bty="n",xaxt="n",ylab="")
mtext("Cumulative probability for",side=2,at=0.5,line=3)
mtext("the three dominance systems",side=2,at=0.5,line=2)
axis(side=1,at=c(-1.347572,-0.768047,-0.1885223,0.3910024,0.9705271,1.550052,2.129577,2.709101,3.288626),labels=FALSE)

mtext("0.8",side=1,at=-1.347572,line=1)
mtext("1.0",side=1,at=-0.768047,line=1)
mtext("1.2",side=1,at=-0.1885223,line=1)
mtext("1.4",side=1,at=0.3910024,line=1)
mtext("1.6",side=1,at=0.9705271,line=1)
mtext("1.8",side=1,at=1.550052,line=1)
mtext("2.0",side=1,at=2.129577,line=1)
mtext("2.2",side=1,at=2.709101,line=1)
mtext("2.4",side=1,at=3.288626,line=1)


polygon(x=c(sexualreceptivity_hours,rev(sexualreceptivity_hours)),y=c(overallprobs_speciesaverages[,2],rev(overallprobs_speciesaverages[,1])),col=co_dominance_color,border=NA)

polygon(x=c(sexualreceptivity_hours,rev(sexualreceptivity_hours)),y=c(overallprobs_speciesaverages[,1],rep(0,11)),col=male_dominance_color,border=NA)

polygon(x=c(sexualreceptivity_hours,rev(sexualreceptivity_hours)),y=c(rep(1,11),rev(overallprobs_speciesaverages[,2])),col=female_dominance_color,border=NA)

mtext("Male dominance",side=1,line=-5,at=1.75)
mtext("Co dominance",side=1,line=-14,at=0)
mtext("Female dominance",side=1,line=-28,at=-0.9)


hist(exp(combined$CanineDimorphism),breaks=40)

barplot(as.matrix(table(combined$strictfdom)/sum(table(combined$strictfdom))),col=dominance_colors)







###############################################################
# Plots for model output with binary predictor traits: two barplots at 0 and 1
###############################################################

## Mating systems

speciesaverage<-combined[is.na(combined$MatSysPMK)==F,]

dat_list_strict <- list(
  R = as.integer(as.factor(speciesaverage$strictfdom)),
  matingsystem = as.integer(as.factor(speciesaverage$MatSysPMK)),
  species = as.integer(as.factor(speciesaverage$corrected_species_id))
)

m_matingsystems <- ulam(
  alist(
    R ~ dordlogit( phi , cutpoints ),
    phi <-a + b[matingsystem],
    a ~ normal( 0 , 5 ),
    b[matingsystem] ~ dnorm(0,5),
    cutpoints ~ dnorm( 0 , 5 )
  ) , data=dat_list_strict , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)

precis(m_matingsystems,depth=2)


# standardized testes size
matingsystems <- seq(from=1,to=4,by=1)
pdat <- data.frame(matingsystems=matingsystems)

overallphi_speciesaverages<-matrix(ncol=1, nrow=4)
for(i in 1:length(matingsystems)){
  overallphi_speciesaverages[i,]<-precis(m_matingsystems,depth=2)[1,1]+precis(m_matingsystems,depth=2)[i+1,1]*matingsystems[i]
}

overallprobs_speciesaverages<-pordlogit( 1:3 , overallphi_speciesaverages , precis(m_matingsystems,depth=2)[6:7,1] )

overallprobs_speciesaverages[,3]<-overallprobs_speciesaverages[,3]-overallprobs_speciesaverages[,2]
overallprobs_speciesaverages[,2]<-overallprobs_speciesaverages[,2]-overallprobs_speciesaverages[,1]
overallprobs_speciesaverages<-t(overallprobs_speciesaverages)


overallprobs_speciesaverages<-overallprobs_speciesaverages[,c(1,2,4,3)]
colnames(overallprobs_speciesaverages)<-c("Monogamous","Polyandrous","Promiscous","Polygynous")
barplot(overallprobs_speciesaverages,col=rev(dominance_colors))

barplot(table(combined$MatSysPMK)[c(1,2,4,3)])



## Arboreality

speciesaverage<-combined[is.na(combined$Strata_Wilman)==F,]

dat_list_strict <- list(
  R = as.integer(as.factor(speciesaverage$strictfdom)),
  strata = as.integer(as.factor(speciesaverage$Strata_Wilman)),
  species = as.integer(as.factor(speciesaverage$corrected_species_id))
)

m_strata <- ulam(
  alist(
    R ~ dordlogit( phi , cutpoints ),
    phi <-a + b[strata],
    a ~ normal( 0 , 5 ),
    b[strata] ~ dnorm(0,5),
    cutpoints ~ dnorm( 0 , 5 )
  ) , data=dat_list_strict , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)

precis(m_strata,depth=2)


# standardized testes size
strata <- seq(from=1,to=3,by=1)
pdat <- data.frame(strata=strata)

overallphi_speciesaverages<-matrix(ncol=1, nrow=3)
for(i in 1:length(strata)){
  overallphi_speciesaverages[i,]<-precis(m_strata,depth=2)[1,1]+precis(m_strata,depth=2)[i+1,1]*matingsystems[i]
}

overallprobs_speciesaverages<-pordlogit( 1:3 , overallphi_speciesaverages , precis(m_strata,depth=2)[5:6,1] )

overallprobs_speciesaverages[,3]<-overallprobs_speciesaverages[,3]-overallprobs_speciesaverages[,2]
overallprobs_speciesaverages[,2]<-overallprobs_speciesaverages[,2]-overallprobs_speciesaverages[,1]
overallprobs_speciesaverages<-t(overallprobs_speciesaverages)


overallprobs_speciesaverages<-overallprobs_speciesaverages[,c(1,3,2)]
colnames(overallprobs_speciesaverages)<-c("Arboreal","Scansorial","Ground")
barplot(overallprobs_speciesaverages,col=rev(dominance_colors))












###############################################################
# Line plots to show all the results
results<-read.csv("intersexualdominance_combinedresults.csv")
continuous_sexualselection<-c("SexualDimorphism_MaleWeight_over_FemaleWeight", "CanineDimorphism","sexratio","M_skew_index","relative_testes_mass","sexualreceptivity_hours","Synchrony")
continuous_femalesociality<-c("female_average_relatedness")
continuous_selforganisation<-c("sexratio","perc_aggression_mm","males","females")
continuous_femalecompetition<-c("env_harshness","rainfall_unpredictability","rainfall_annualvariation","NDVI_Seasonality","r_seasonality_value","homerange_overlap","females","relative_femalecaninesize")
discrete_sexualselection<-c("MatSysPMK","Strata_Wilman","ovulation_signs")
discrete_selforganisation<-c("between_groupconflict")
discrete_femalesociality<-c("SocOrgPMK","jointaggression_females","sexbias_dispersal")
discrete_femalecompetition<-c("origin","female_evictions","female_infanticide")

sexualselection_results<-results[rowSums(sapply(c(continuous_sexualselection,discrete_sexualselection), grepl, results$continuouspredictor, ignore.case=TRUE))==1,]
sexualselection_results<-sexualselection_results[sexualselection_results$outcome %in% "perc_won" & sexualselection_results$phylogeny %in% "Yes",]
sexualselection_results$median<-sapply(1 : nrow(sexualselection_results), function(x) median(c(sexualselection_results[x,]$estimate.upper,sexualselection_results[x,]$estimate.lower)) )

femalesociality_results<-results[rowSums(sapply(c(continuous_femalesociality,discrete_femalesociality), grepl, results$continuouspredictor, ignore.case=TRUE))==1,]
femalesociality_results<-femalesociality_results[femalesociality_results$outcome %in% "perc_won" & femalesociality_results$phylogeny %in% "Yes",]
femalesociality_results$median<-sapply(1 : nrow(femalesociality_results), function(x) median(c(femalesociality_results[x,]$estimate.upper,femalesociality_results[x,]$estimate.lower)) )

selforganisation_results<-results[rowSums(sapply(c(continuous_selforganisation,discrete_selforganisation), grepl, results$continuouspredictor, ignore.case=TRUE))==1,]
selforganisation_results<-selforganisation_results[selforganisation_results$outcome %in% "perc_won" & selforganisation_results$phylogeny %in% "Yes",]
selforganisation_results$median<-sapply(1 : nrow(selforganisation_results), function(x) median(c(selforganisation_results[x,]$estimate.upper,selforganisation_results[x,]$estimate.lower)) )

femalecompetition_results<-results[rowSums(sapply(c(continuous_femalecompetition,discrete_femalecompetition), grepl, results$continuouspredictor, ignore.case=TRUE))==1,]
femalecompetition_results<-femalecompetition_results[femalecompetition_results$outcome %in% "perc_won" & femalecompetition_results$phylogeny %in% "Yes",]
femalecompetition_results$median<-sapply(1 : nrow(femalecompetition_results), function(x) median(c(femalecompetition_results[x,]$estimate.upper,femalecompetition_results[x,]$estimate.lower)) )

sexualselection_results$colour<-ifelse(sexualselection_results$association_present %in% "confident","black","darkgrey")

plot( NULL , type="n" , xlab="relationship" ,
      xlim=c(-1,1) , ylim=c(-1,1) ,bty="n",xaxt="n",ylab="")
rect(xleft=-2,ybottom=-0.25,ytop=0.25,xright=2,col="grey95",border=NA)
for (i in 1:nrow(sexualselection_results)){
 abline(a=0,b=sexualselection_results[i,]$median/2,col=sexualselection_results[i,]$colour)
  }


###############################################################
# Forest plot to show all the results


# Forest plot
results<-read.csv("intersexualdominance_combinedresults.csv")
summarized_results<-results %>% group_by(continuouspredictor) %>% summarise(mean(estimate.lower),mean(estimate.upper))
summarized_results<-as.data.frame(summarized_results)
colnames(summarized_results)<-c("contrast","lower","upper")

summarized_results$index<-c(1:nrow(summarized_results))
summarized_results$median<-sapply(1 : nrow(summarized_results), function(x) median(c(summarized_results[x,]$upper,summarized_results[x,]$lower)) )


ggplot(data=summarized_results, aes(y=index, x=median, xmin=lower, xmax=upper)) +
  geom_point() +
  geom_errorbarh(height=.1) +
  scale_y_continuous(name = "", breaks=1:nrow(summarized_results), labels=summarized_results$contrast)



results$beta<-sapply(1 : nrow(results), function(x) median(c(results[x,]$estimate.upper,results[x,]$estimate.lower)) )

results$se<-sapply(1 : nrow(results), function(x) results[x,]$beta-results[x,]$estimate.lower )

sresults<-results[results$phylogeny %in% "Yes",]

sresults<-select(sresults,-X)
colnames(sresults)[1]<-"name"
sresults<-sresults[,c(2,1,3,4,5,6,7,8,9)]
colnames(sresults)[c(1,2)]<-c("name","outcome")

df_linear <-
  ggforestplot::df_linear_associations %>%
  dplyr::arrange(name) %>%
  dplyr::filter(dplyr::row_number() <= 30)

# Forestplot
forestplot(
  df = sresults[1:20,],
  estimate = beta,
  logodds = FALSE,
  colour = outcome,
  title = "Association with intersexual dominance",
  xlab = "correlation coefficient"
)







###############################################################
###############################################################



# Sex ratio dimorphism
speciesaverage<-as.data.frame(combined %>% group_by(corrected_species_id) %>% summarise(mean(CanineDimorphism,na.rm=T),mean(as.numeric(strictfdom))))

colnames(speciesaverage)<-c("corrected_species_id","sizedimorphism","strictfdom")
speciesaverage<-speciesaverage[is.na(speciesaverage$sizedimorphism)==F,]
speciesaverage$strictfdom<-round(speciesaverage$strictfdom,0)

dat_list_strict <- list(
  R = as.integer(as.factor(speciesaverage$strictfdom)),
  sizedimorphism = standardize(exp(speciesaverage$sizedimorphism)),
  species = as.integer(as.factor(speciesaverage$corrected_species_id))
)

m_speciesaverage_sizedimorphism <- ulam(
  alist(
    R ~ dordlogit( phi , cutpoints ),
    phi <-a + b*sizedimorphism ,
    a ~ normal( 0 , 5 ),
    b ~ dnorm(0,5),
    cutpoints ~ dnorm( 0 , 5 )
  ) , data=dat_list_strict , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)

precis(m_speciesaverage_sizedimorphism,depth=2)


# standardized size dimorphism
sizedimorphism <- seq(from=-1.5,to=4.5,by=0.5)
pdat <- data.frame(sizedimorphism=sizedimorphism)

overallphi_speciesaverages<-matrix(ncol=1, nrow=13)
for(i in 1:length(sizedimorphism)){
  overallphi_speciesaverages[i,]<-precis(m_speciesaverage_sizedimorphism,depth=2)[1,1]+precis(m_speciesaverage_sizedimorphism,depth=2)[2,1]*sizedimorphism[i]
}

overallprobs_speciesaverages<-pordlogit( 1:3 , overallphi_speciesaverages , precis(m_speciesaverage_sizedimorphism,depth=2)[3:4,1] )

(0.9-mean(exp(speciesaverage$sizedimorphism)))/sd(exp(speciesaverage$sizedimorphism))


plot( NULL , type="n" , xlab="Male canine size relative to female canine size" ,
      xlim=c(min(dat_list_strict$sizedimorphism,na.rm=T),max(dat_list_strict$sizedimorphism,na.rm=T)) , ylim=c(0,1) ,  yaxp=c(0,1,4) ,bty="n",xaxt="n",ylab="")
mtext("Cumulative probability for",side=2,at=0.5,line=3)
mtext("the three dominance systems",side=2,at=0.5,line=2)
axis(side=1,at=c((0.8-mean(exp(speciesaverage$sizedimorphism)))/sd(exp(speciesaverage$sizedimorphism)),(1.0-mean(exp(speciesaverage$sizedimorphism)))/sd(exp(speciesaverage$sizedimorphism)),(1.5-mean(exp(speciesaverage$sizedimorphism)))/sd(exp(speciesaverage$sizedimorphism)),(2.0-mean(exp(speciesaverage$sizedimorphism)))/sd(exp(speciesaverage$sizedimorphism)),(2.5-mean(exp(speciesaverage$sizedimorphism)))/sd(exp(speciesaverage$sizedimorphism)),(3.0-mean(exp(speciesaverage$sizedimorphism)))/sd(exp(speciesaverage$sizedimorphism)),(3.5-mean(exp(speciesaverage$sizedimorphism)))/sd(exp(speciesaverage$sizedimorphism)),(4.0-mean(exp(speciesaverage$sizedimorphism)))/sd(exp(speciesaverage$sizedimorphism)),(4.5-mean(exp(speciesaverage$sizedimorphism)))/sd(exp(speciesaverage$sizedimorphism)),(5.0-mean(exp(speciesaverage$sizedimorphism)))/sd(exp(speciesaverage$sizedimorphism))),labels=FALSE)

mtext("0.8",side=1,at=(0.8-mean(exp(speciesaverage$sizedimorphism)))/sd(exp(speciesaverage$sizedimorphism)),line=1)
mtext("1.0",side=1,at=(1.0-mean(exp(speciesaverage$sizedimorphism)))/sd(exp(speciesaverage$sizedimorphism)),line=1)
mtext("1.5",side=1,at=(1.5-mean(exp(speciesaverage$sizedimorphism)))/sd(exp(speciesaverage$sizedimorphism)),line=1)
mtext("2.0",side=1,at=(2.0-mean(exp(speciesaverage$sizedimorphism)))/sd(exp(speciesaverage$sizedimorphism)),line=1)
mtext("2.5",side=1,at=(2.5-mean(exp(speciesaverage$sizedimorphism)))/sd(exp(speciesaverage$sizedimorphism)),line=1)
mtext("3.0",side=1,at=(3.0-mean(exp(speciesaverage$sizedimorphism)))/sd(exp(speciesaverage$sizedimorphism)),line=1)
mtext("3.5",side=1,at=(3.5-mean(exp(speciesaverage$sizedimorphism)))/sd(exp(speciesaverage$sizedimorphism)),line=1)
mtext("4.0",side=1,at=(4.0-mean(exp(speciesaverage$sizedimorphism)))/sd(exp(speciesaverage$sizedimorphism)),line=1)
mtext("4.5",side=1,at=(4.5-mean(exp(speciesaverage$sizedimorphism)))/sd(exp(speciesaverage$sizedimorphism)),line=1)
mtext("5.0",side=1,at=(5.0-mean(exp(speciesaverage$sizedimorphism)))/sd(exp(speciesaverage$sizedimorphism)),line=1)

polygon(x=c(sizedimorphism,rev(sizedimorphism)),y=c(overallprobs_speciesaverages[,2],rev(overallprobs_speciesaverages[,1])),col=co_dominance_color,border=NA)

polygon(x=c(sizedimorphism,rev(sizedimorphism)),y=c(overallprobs_speciesaverages[,1],rep(0,13)),col=male_dominance_color,border=NA)

polygon(x=c(sizedimorphism,rev(sizedimorphism)),y=c(rep(1,13),rev(overallprobs_speciesaverages[,2])),col=female_dominance_color,border=NA)

mtext("Male dominance",side=1,line=-5,at=1.75)
mtext("Co dominance",side=1,line=-14,at=0)
mtext("Female dominance",side=1,line=-28,at=-0.9)


hist(exp(combined$SexualDimorphism_MaleWeight_over_FemaleWeight),breaks=40)

barplot(as.matrix(table(combined$strictfdom)/sum(table(combined$strictfdom))),col=dominance_colors)










################################################################################

#figure 1 (supplement) variation in quantitative intersexual dominance, echo=FALSE}
previouspar<-par()
op <- par(oma=c(2.5,2.5,1.5,0.2), mar=c(2.5,2.5,2.5,3.5), mfrow=c(1,1))

hist(combined$perc_won_females,breaks=50,xlab="Percentage of aggressive interactions won by females",ylab="Number of observations",main=NULL,cex.lab=1.5,xaxt="n", col = "grey75", border = "black",plot.new=TRUE)

Axis(side=1, labels=FALSE, at=c(0,10,20,30,40,50,60,70,80,90,100))
mtext("0", side=1,line=1.0,at=0,cex=1,las=1)
mtext("10", side=1,line=1.0,at=10,cex=1,las=1)
mtext("20", side=1,line=1.0,at=20,cex=1,las=1)
mtext("30", side=1,line=1.0,at=30,cex=1,las=1)
mtext("40", side=1,line=1.0,at=40,cex=1,las=1)
mtext("50", side=1,line=1.0,at=50,cex=1,las=1)
mtext("60", side=1,line=1.0,at=60,cex=1,las=1)
mtext("70", side=1,line=1.0,at=70,cex=1,las=1)
mtext("80", side=1,line=1.0,at=80,cex=1,las=1)
mtext("90", side=1,line=1.0,at=90,cex=1,las=1)
mtext("100", side=1,line=1.0,at=100,cex=1,las=1)

mtext("Number of observations",side=2,line=3.0,at=-4,cex=1.6,las=3,adj=0)

mtext("Proportion of intersexual fights won by females",side=1,line=3.0,at=12,cex=1.6,las=1,adj=0)
abline(v=33.3, col="darkred")
abline(v=66.6, col="darkred")

mtext("Males dominant over females",side=1,line=-12.0,at=4,cex=0.9,las=1,adj=0)

mtext("Co-dominance",side=1,line=-12.0,at=41,cex=0.9,las=1,adj=0)

mtext("Females dominant over males",side=1,line=-12.0,at=68,cex=0.9,las=1,adj=0)

par<-previouspar


##############################################################################
# Plot the percentage of fights won by females across the phylogeny
originalpar<-par()

op <- par(oma=c(0,2.5,0,0), mar=c(2.5,4.5,4.2,2.2), mfrow=c(1,1))

dotTree(perctree,percdata,legend=TRUE)

mtext("percent intersexual fights won by females",side=1,line=-2.5,at=15,cex=2,las=1,adj=0,col="darkblue")

par<-originalpar


##############################################################################
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

dominancecolours<-colnames(ancestraldominance$lik.anc)[max.col(ancestraldominance$lik.anc, ties.method = "random")]

dominancecolours[dominancecolours=="1"]<-viridis(3)[3]
dominancecolours[dominancecolours=="2"]<-viridis(3)[2]
dominancecolours[dominancecolours=="3"]<-viridis(3)[1]


ancestraldominancealt<-fastAnc(onetree, x=phylostrictdata$strictfdom, vars=FALSE, CI=FALSE)

dominancecolours<-round(ancestraldominancealt,0)+5

plot(stricttree, show.tip.label=TRUE,font=1,cex=0.75,no.margin=TRUE,adj=1,edge.color=dominancecolours,edge.width=5)

legend(x="topleft",legend=c("Male dominance","Co dominance","Female dominance"),lwd=3,col=c(viridis(3)[3],viridis(3)[2],viridis(3)[1]),bty="n")



##############################################################################
# Plot the relaxed 2-way classification of intersexual dominance on a phylogeny, including reconstructed states
mostly_females_phylodata_1<-combined
mostly_females_phylodata_1$mostlyfdom<-as.numeric(as.factor(combined$mostlyfdom))
mostly_females_phylodata_1<-as.data.frame(mostly_females_phylodata_1 %>% group_by(corrected_species_id) %>% summarise(mean(mostlyfdom,na.rm=T )) )
mostly_females_phylodata_2<-as.data.frame(mostly_females_phylodata_1[,2])
row.names(mostly_females_phylodata_2)<-mostly_females_phylodata_1[,1]
colnames(mostly_females_phylodata_2)<-"mostlyfdom"
mostly_females_phylodata_3<-as.data.frame(mostly_females_phylodata_2[is.na(mostly_females_phylodata_2$mostlyfdom)==FALSE,])
row.names(mostly_females_phylodata_3)<-row.names(mostly_females_phylodata_2)[is.na(mostly_females_phylodata_2$mostlyfdom)==FALSE]
colnames(mostly_females_phylodata_3)<-"mostlyfdom"
missing<-treedata(inputtree,data=mostly_females_phylodata_3,warnings=FALSE)
mostlytree<-missing$phy
speciesnames<-mostlytree$tip.label
mostlydata<-data.frame(mostly_females_phylodata_3[speciesnames,])
colnames(mostlydata)<-"mostlyfdom"
row.names(mostlydata)<-speciesnames
phylomostlydata<-round(mostlydata,0)
phylomostlydata<-as.data.frame(phylomostlydata)
grafentree<-compute.brlen(mostlytree,method="Grafen")
onetree<-compute.brlen(mostlytree,1)

ancestraldominance<-ace(phylomostlydata$mostlyfdom,mostlytree,type="discrete",model="ARD")

dominancecolours<-colnames(ancestraldominance$lik.anc)[max.col(ancestraldominance$lik.anc, ties.method = "random")]

dominancecolours[dominancecolours==1]<-"grey50"
dominancecolours[dominancecolours==2]<-"maroon3"

phylomostlydata$mostlyfdom<-as.factor(phylomostlydata$mostlyfdom)
cols<-setNames(c("grey70","maroon1"),levels(phylomostlydata$mostlyfdom))
plotTree(mostlytree,type="fan",fsize=0.7,ftype="i",lwd=1)
nodelabels(node=1:mostlytree$Nnode+Ntip(mostlytree),
           pie=ancestraldominance$lik.anc,piecol=cols,cex=0.2)
tiplabels(pie=to.matrix(phylomostlydata[mostlytree$tip.label,1],
                        levels(phylomostlydata$mostlyfdom)),piecol=cols,cex=0.2)

legend(x="bottom",legend=c("Male dominance","Female dominance"),lwd=3,col=c("maroon1","grey50"),bty="n")


########################################################################################

# Combined plot

# Plot the phylogeny
library(ggtree)
library(caper)

safe_colorblind_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

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


dominancedata<-data.frame(dominancedata[speciesnames,])
dominancedata$colours<-1
dominancedata[dominancedata$strictfdom==1,]$colours<-6
dominancedata[dominancedata$strictfdom==2,]$colours<-4
dominancedata[dominancedata$strictfdom==3,]$colours<-17


relaxeddominancedata<-data.frame(relaxeddominancedata[speciesnames,])
relaxeddominancedata$colours<-1
relaxeddominancedata[relaxeddominancedata$relaxedfdom==1,]$colours<-6
relaxeddominancedata[relaxeddominancedata$relaxedfdom==2,]$colours<-17

# Split species at the level of 'superfamilies' (https://en.wikipedia.org/wiki/Primate)

getMRCA(perctree,c("Macaca_mulatta","Rhinopithecus_bieti")) # Cercopithecoidea 177
sum(clade.matrix(perctree)$clade.matrix[195,]) # 36 species
getMRCA(perctree,c("Pan_paniscus","Gorilla_gorilla")) # Hominoidea 174
sum(clade.matrix(perctree)$clade.matrix[192,]) # 4 species
getMRCA(perctree,c("Ateles_belzebuth","Cacajao_melanocephalus")) # Ceboidea 140
sum(clade.matrix(perctree)$clade.matrix[149,]) # 43 species
getMRCA(perctree,c("Lepilemur_leucopus","Daubentonia_madagascariensis")) # Lemuroidea 109
sum(clade.matrix(perctree)$clade.matrix[118,]) #31 species
# Plus Loris_lydekkerianus Lorisoidea 108


#tree2 <- groupClade(perctree, c(109,140,174,177))
tree2 <- groupClade(perctree, c(118,149,192,195))
levels(attributes(tree2)$group)<-c(0:5)
attributes(tree2)$group[1]<-5
treeplot<-ggtree(tree2, aes(color=group)) + 
  scale_color_manual(values=safe_colorblind_palette[1:7])


p2 <- facet_plot(treeplot, panel="Percentage fights won by females", data=fightswondata, geom=geom_point, aes(x=perc_won_females), cex=3,color=c(rep(safe_colorblind_palette[6],1),rep(safe_colorblind_palette[2],31),rep(safe_colorblind_palette[4],4),rep(safe_colorblind_palette[5],36),rep(safe_colorblind_palette[3],43)))

p2+theme_tree2()

p3 <- facet_plot(treeplot, panel="Intersexual dominance", data=dominancedata[,c(1,3)], geom=geom_point, aes(x=placeholder), cex=1.5,color=c(rep(safe_colorblind_palette[6],1),rep(safe_colorblind_palette[2],31),rep(safe_colorblind_palette[3],34),rep(safe_colorblind_palette[4],4),rep(safe_colorblind_palette[5],36)),pch=c(dominancedata$strictfdom))

p3+theme_tree2()

p4 <- facet_plot(treeplot, panel="Intersexual dominance", data=relaxeddominancedata[,c(1,3)], geom=geom_point, aes(x=placeholder), cex=1.5,color=c(rep(safe_colorblind_palette[6],1),rep(safe_colorblind_palette[2],31),rep(safe_colorblind_palette[3],34),rep(safe_colorblind_palette[4],4),rep(safe_colorblind_palette[5],36)),pch=c(relaxeddominancedata$relaxedfdom))

p4+theme_tree2()


d=fortify(perctree)
dd = subset(d, isTip)
dominancedata<-dominancedata[dd$label[order(dd$y, decreasing=FALSE)],]

p5 <- facet_plot(treeplot, panel="Intersexual dominance", data=dominancedata[,c(1,2)], geom=geom_point, aes(x=strictfdom), cex=1.5,color=c(rep(safe_colorblind_palette[6],1),rep(safe_colorblind_palette[2],31),rep(safe_colorblind_palette[3],34),rep(safe_colorblind_palette[4],4),rep(safe_colorblind_palette[5],36)),pch=c(dominancedata$colours))

p5+theme_tree2()


d=fortify(perctree)
dd = subset(d, isTip)
relaxeddominancedata<-relaxeddominancedata[dd$label[order(dd$y, decreasing=FALSE)],]

p6 <- facet_plot(treeplot, panel="Intersexual dominance", data=relaxeddominancedata[,c(1,2)], geom=geom_point, aes(x=relaxedfdom), cex=1.5,color=c(rep(safe_colorblind_palette[6],1),rep(safe_colorblind_palette[2],31),rep(safe_colorblind_palette[3],34),rep(safe_colorblind_palette[4],4),rep(safe_colorblind_palette[5],36)),pch=c(relaxeddominancedata$colours))

p6+theme_tree2()




######################################################################################################



perc_won_females_phylodata_1<-as.data.frame(combined %>% group_by(corrected_species_id) %>% summarise(mean(perc_won_females,na.rm=T)))

# We turn this into a data frame with just a single column containing the values
perc_won_females_phylodata_2<-as.data.frame(perc_won_females_phylodata_1[,2])
row.names(perc_won_females_phylodata_2)<-perc_won_females_phylodata_1[,1]
colnames(perc_won_females_phylodata_2)<-"perc_won_females"

# We remove species where no quantitative data is available
perc_won_females_phylodata_3<-as.data.frame(perc_won_females_phylodata_2[is.na(perc_won_females_phylodata_2$perc_won_females)==FALSE,])
row.names(perc_won_females_phylodata_3)<-row.names(perc_won_females_phylodata_2)[is.na(perc_won_females_phylodata_2$perc_won_females)==FALSE]
colnames(perc_won_females_phylodata_3)<-"perc_won_females"

# We match the species that in the data frame to the species in the phylogenetic tree
missing<-treedata(inputtree,data=perc_won_females_phylodata_3,warnings=FALSE)
# We remove the species from the tree for which we have no data
perctree<-missing$phy
# We remove the species from the dataset which are not in the tree
speciesnames<-perctree$tip.label
percdata<-as.data.frame(perc_won_females_phylodata_3[rownames(perc_won_females_phylodata_3) == speciesnames,])
percdata<-data.frame(perc_won_females_phylodata_3[speciesnames,])
colnames(percdata)<-"perc_won_females"
row.names(percdata)<-speciesnames

percdata$species<-rownames(percdata)
strictdata$species<-rownames(strictdata)
mostlydata$species<-rownames(mostlydata)

combined_data<-left_join(percdata,strictdata,by="species")
combined_data<-left_join(combined_data,mostlydata,by="species")
combined_data$strictfdom<-round(combined_data$strictfdom,0)
combined_data[41,]$mostlyfdom<-1
combined_data$mostlyfdom<-round(combined_data$mostlyfdom,0)
combined_data$mostlyfdom<-case_when(combined_data$mostlyfdom==1 ~ 2, combined_data$mostlyfdom==2 ~ 1)
rownames(combined_data)<-combined_data$species

missing<-treedata(inputtree,data=combined_data,warnings=FALSE)
# We remove the species from the tree for which we have no data
combinedtree<-missing$phy


p <- ggtree(stricttree)
p1 <- p %<+% combined_data$strictfdom + geom_tippoint(aes(color=location))
p2 <- facet_plot(p1, panel="dot", data=d2, geom=geom_point, 
                 aes(x=val), color='firebrick') + theme_tree2()



p <- ggtree(stricttree, branch.length = "none") + 
  geom_tiplab() + theme(legend.position='none')

plot.phylo(stricttree)

plot1<-plotTree(stricttree)
plot2<-
  
  library(cowplot)
plot_grid(plot1, plot2, labels = c('A', 'B'),rel_widths = c(5,3))


tr <- rtree(10)
dd = data.frame(id=tr$tip.label, value=abs(rnorm(10)))
p <- ggtree(tr)
facet_plot(p, 'Trait', data = dd, geom=geom_point, mapping=aes(x=value))

rownames(combined_data)<-combined_data$species
dotTree(stricttree,combined_data$perc_won_females,length=10,ftype="i")

strictfdom_data<-as.matrix(combined_data$strictfdom)
rownames(strictfdom_data)<-combined_data$species
colnames(strictfdom_data)<-"strictfdom"
plotTree.barplot(combinedtree,strictfdom_data)


combined_data<-combined_data[,c(1,3,4)]
combined_data$perc_won_females<-combined_data$perc_won_females/100

colors<-setNames(c("black","grey50","grey10"),c("1","2","3"))

perc_won_females_data<-as.matrix(combined_data$perc_won_females)
rownames(perc_won_females_data)<-combined_data$species
colnames(perc_won_females_data)<-"perc_won_females"

dotTree(combinedtree,perc_won_females_data,length=10,col="blye")
dotTree(combinedtree,strictfdom_data,type="discrete",length=10,colors=colors,add=T)



##############################################################################

combined_data
combinedtree

par(mfrow=c(1,3))
plotTree(combinedtree,mar=c(5.1,1.1,2.1,0.1))
par(mar=c(5.1,0.1,2.1,1.1))
boxplot(combined_data$perc_won_females~factor(rownames(combined_data),levels=combinedtree$tip.label),horizontal=TRUE,
        axes=FALSE,xlim=c(1,Ntip(combinedtree)))
axis(1)
title(xlab="log(body size)")



##############################################################################

# Relationship between intersexual dominance and the mating system

# Percentage fights won
distribution_strictfdom_matingsystem<-table(combined$MatSysPMK,combined$strictfdom)
colnames(distribution_strictfdom_matingsystem)<-c("male_dominance","co_dominance","female_dominance")
distribution_strictfdom_matingsystem[1,]<-distribution_strictfdom_matingsystem[1,]/sum(distribution_strictfdom_matingsystem[1,])
distribution_strictfdom_matingsystem[2,]<-distribution_strictfdom_matingsystem[2,]/sum(distribution_strictfdom_matingsystem[2,])
distribution_strictfdom_matingsystem[3,]<-distribution_strictfdom_matingsystem[3,]/sum(distribution_strictfdom_matingsystem[3,])
distribution_strictfdom_matingsystem[4,]<-distribution_strictfdom_matingsystem[4,]/sum(distribution_strictfdom_matingsystem[4,])

boxplot(combined$perc_won_females~combined$MatSysPMK)


library(ggdist)

df_matingsystem <- combined[ complete.cases(combined$perc_won_females,combined$MatSysPMK),]

df = data.frame(
  group=df_matingsystem$MatSysPMK,
  Percentage_IntersexualFights_WonByFemales=df_matingsystem$perc_won_females
)

ggplot(df)+aes(y=group,x=Percentage_IntersexualFights_WonByFemales)+stat_halfeye()



##############################################################################

# Relationship between testes size and strict dominance classification


df_testessize <- combined[ complete.cases(combined$strictfdom,combined$testes_mass),]

df_testessize[df_testessize$strictfdom==1,]$strictfdom<-"male_dominance"
df_testessize[df_testessize$strictfdom==2,]$strictfdom<-"co_dominance"
df_testessize[df_testessize$strictfdom==3,]$strictfdom<-"female_dominance"

df = data.frame(
  strictfdom=df_testessize$strictfdom,
  testes_size=df_testessize$testes_mass
)

ggplot(df)+aes(y=strictfdom,x=testes_size)+stat_halfeye()





##############################################################################

# Relationship between sex ratio and strict dominance classification

df_sexratio <- combined[ complete.cases(combined$strictfdom,combined$sexratio),]

df_sexratio[df_sexratio$strictfdom==1,]$strictfdom<-"male_dominance"
df_sexratio[df_sexratio$strictfdom==2,]$strictfdom<-"co_dominance"
df_sexratio[df_sexratio$strictfdom==3,]$strictfdom<-"female_dominance"

df = data.frame(
  strictfdom=df_sexratio$strictfdom,
  sexratio_malesperfemale=df_sexratio$sexratio
)

ggplot(df)+aes(y=strictfdom,x=sexratio_malesperfemale)+stat_halfeye()



##############################################################################

# Relationship between sexual dimorphism and strict dominance classification

library(ggdist)
library(cowplot)

df_sexratio <- combined[ complete.cases(combined$strictfdom,combined$SexualDimorphism_MaleWeight_over_FemaleWeight),]

df_sexratio[df_sexratio$strictfdom==1,]$strictfdom<-"c) Strict male dominance"
df_sexratio[df_sexratio$strictfdom==2,]$strictfdom<-"b) Co dominance"
df_sexratio[df_sexratio$strictfdom==3,]$strictfdom<-"a) Strict female dominance"

df = data.frame(
  strictfdom=df_sexratio$strictfdom,
  SexualDimorphism_MaleWeight_over_FemaleWeight=df_sexratio$SexualDimorphism_MaleWeight_over_FemaleWeight
)

plot1<-ggplot(df)+aes(y=strictfdom,x=SexualDimorphism_MaleWeight_over_FemaleWeight)+stat_halfeye()+
  theme(        axis.ticks.y=element_blank(),
                axis.title.y = element_blank(),
                axis.text=element_text(size=14)
  )+
  scale_y_discrete(labels = c('Strict female dominance', 'Co dominance', 'Strict male dominance'))+
  scale_x_continuous(name="Body size dimorphism")


df_sexratiocanines <- combined[ complete.cases(combined$strictfdom,combined$CanineDimorphism),]
df_sexratiocanines$CanineDimorphism<-df_sexratiocanines$CanineDimorphism
df_sexratiocanines[df_sexratiocanines$strictfdom==1,]$strictfdom<-"c) Strict male dominance"
df_sexratiocanines[df_sexratiocanines$strictfdom==2,]$strictfdom<-"b) Co dominance"
df_sexratiocanines[df_sexratiocanines$strictfdom==3,]$strictfdom<-"a) Strict female dominance"

df = data.frame(
  strictfdom=df_sexratiocanines$strictfdom,
  CanineDimorphism=df_sexratiocanines$CanineDimorphism
)

plot2<-ggplot(df)+aes(y=strictfdom,x=CanineDimorphism)+stat_halfeye()+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y = element_blank(),
        axis.text=element_text(size=14)
  )+
  scale_x_continuous(name="Canine size dimorphism")


plot_grid(plot1, plot2, labels = c('A', 'B'),rel_widths = c(5,3))


##############################################################################

# Relationship between synchrony and strict dominance classification

df_sexratio <- combined[ complete.cases(combined$strictfdom,combined$Synchrony),]

df_sexratio[df_sexratio$strictfdom==1,]$strictfdom<-"male_dominance"
df_sexratio[df_sexratio$strictfdom==2,]$strictfdom<-"co_dominance"
df_sexratio[df_sexratio$strictfdom==3,]$strictfdom<-"female_dominance"

df = data.frame(
  strictfdom=df_sexratio$strictfdom,
  Synchrony=df_sexratio$Synchrony
)

plot1<-ggplot(df)+aes(y=strictfdom,x=Synchrony)+stat_halfeye()


df_sexratiocanines <- combined[ complete.cases(combined$strictfdom,combined$r_seasonality_value),]

df_sexratiocanines[df_sexratiocanines$strictfdom==1,]$strictfdom<-"male_dominance"
df_sexratiocanines[df_sexratiocanines$strictfdom==2,]$strictfdom<-"co_dominance"
df_sexratiocanines[df_sexratiocanines$strictfdom==3,]$strictfdom<-"female_dominance"

df = data.frame(
  strictfdom=df_sexratiocanines$strictfdom,
  Seasonality=df_sexratiocanines$r_seasonality_value
)

plot2<-ggplot(df)+aes(y=strictfdom,x=Seasonality)+stat_halfeye()


plot_grid(plot1, plot2, labels = c('A', 'B'))




# Figures for Elise 01 June 2022


# 1. A graph illustrating the relationship between testes size and intersexual dominance (strict classification) selecting only polygynous and promiscuous species (i.e. excluding monogamous and polyandrous ones).

polpro<-combined[combined$MatSysPMK %in% c("POL","PRO"),]

boxplot(polpro$testes_mass/polpro$body_mass~polpro$strictfdom)



df_testessize <- polpro[ complete.cases(polpro$strictfdom,polpro$testes_mass),]
df_testessize$reltestessize<-(df_testessize$testes_mass/df_testessize$body_mass)*100

df_testessize[df_testessize$strictfdom==1,]$strictfdom<-"c) male dominance"
df_testessize[df_testessize$strictfdom==2,]$strictfdom<-"b) co-dominance"
df_testessize[df_testessize$strictfdom==3,]$strictfdom<-"a) female dominance"

df = data.frame(
  strictfdom=df_testessize$strictfdom,
  testes_size=df_testessize$reltestessize
)

bxp<-ggplot(df)+aes(y=strictfdom,x=testes_size,fill=factor(strictfdom))+stat_halfeye()


bxp <- bxp + labs(
  x = "Testes Size (percent body weight)", y = "Intersexual dominance")+theme(axis.text=element_text(size=12),
                                                                              axis.title=element_text(size=14,face="bold"))+theme(legend.position = "none")+scale_fill_viridis(discrete = T)




# 2. The plot illustrating the relationship between SSD and %age fights won by females (the one you showed this morning, but perhaps with a regression line and larger dots?)

plot(combined$perc_won_females~combined$SexualDimorphism_MaleWeight_over_FemaleWeight,bty="n",ann=F,xlim=c(0.75,2.5),cex=2,xaxt="n")
axis(side=1, at=c(0.75,1.0,1.25,1.5,1.75,2.0,2.25,2.5), labels=c("0.75","1.00","1.25","1.50","1.75","2.00","2.25","2.50"))
mtext("Male size relative to female size",side=1,line=3.5,at=1.2,cex=1.5,las=1,adj=0)
abline(v=1,lty=2,col="grey")
mtext("females larger",side=1,line=2.0,at=0.75,cex=0.8,las=1,adj=0)
mtext("males larger",side=1,line=2.0,at=2.25,cex=0.8,las=1,adj=0)
mtext("Percent fights won by females",side=2,line=2.5,at=0,cex=1.5,las=3,adj=0)


# 2a. Variation of sexual dimorphism plot that also includes environmental harshness

plot(combined[combined$SexualDimorphism_MaleWeight_over_FemaleWeight>1.11,]$perc_won_females~combined[combined$SexualDimorphism_MaleWeight_over_FemaleWeight>1.11,]$SexualDimorphism_MaleWeight_over_FemaleWeight,bty="n",ann=F,cex=combined[combined$SexualDimorphism_MaleWeight_over_FemaleWeight>1.11,]$env_harshness+2,xaxt="n",xlim=c(0.75,2.5))
points(combined[combined$SexualDimorphism_MaleWeight_over_FemaleWeight<1.11,]$perc_won_females~combined[combined$SexualDimorphism_MaleWeight_over_FemaleWeight<1.11,]$SexualDimorphism_MaleWeight_over_FemaleWeight,bty="n",ann=F,xlim=c(0.75,2.5),cex=combined[combined$SexualDimorphism_MaleWeight_over_FemaleWeight<1.11,]$env_harshness+2,col="red")
axis(side=1, at=c(0.75,1.0,1.25,1.5,1.75,2.0,2.25,2.5), labels=c("0.75","1.00","1.25","1.50","1.75","2.00","2.25","2.50"))
mtext("Male size relative to female size",side=1,line=3.5,at=1.2,cex=1.5,las=1,adj=0)
abline(v=1.1,lty=2,col="grey")
mtext("females larger",side=1,line=2.0,at=0.75,cex=0.8,las=1,adj=0)
mtext("males larger",side=1,line=2.0,at=2.25,cex=0.8,las=1,adj=0)
mtext("Percent fights won by females",side=2,line=2.5,at=0,cex=1.5,las=3,adj=0)
legend(x="topright",legend=c("larger circle indicates harsher environment"),pch=c("0"),lwd=0,col=c("black"),bty="n")




# 3. A plot illustrating the relationship between intersexual dominance (strict classification) and mating systems -- perhaps barplots indicating the proportion of species in different categories of intersexual dominance for each mating system?

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

summarizedtable<-as.data.frame(summarizedtable)

summarizedtable<-rbind(summarizedtable,c(1,1,0))
summarizedtable<-rbind(summarizedtable,c(2,1,0))

summarizedtable[summarizedtable$MatingSystem ==1, ]$MatingSystem<-"Monogamy"
summarizedtable[summarizedtable$MatingSystem ==2, ]$MatingSystem<-"Polyandry"
summarizedtable[summarizedtable$MatingSystem ==3, ]$MatingSystem<-"Polygyny"
summarizedtable[summarizedtable$MatingSystem ==4, ]$MatingSystem<-"Promiscuity"

summarizedtable[summarizedtable$StrictFemdom ==1, ]$StrictFemdom<-"c) male dominance"
summarizedtable[summarizedtable$StrictFemdom ==2, ]$StrictFemdom<-"b) co-dominance"
summarizedtable[summarizedtable$StrictFemdom ==3, ]$StrictFemdom<-"a) female dominance"

ggplot (as.data.frame(summarizedtable), aes(x=MatingSystem, y=NumberOfSpecies, fill=StrictFemdom)) + geom_bar (stat="identity", position ="dodge")+scale_fill_viridis(discrete = T)+theme_bw()+ labs(
  x = "Mating System", y = "Number of Species", fill="Intersexual dominance")+theme(axis.text=element_text(size=12),
                                                                                    axis.title=element_text(size=14,face="bold"))+theme(legend.position = "top")



# 4. Is there any way to illustrate the phylogenic reconstruction exploring the links between monogamy and intersexual dom, and between SSD and intersexual dom? 
