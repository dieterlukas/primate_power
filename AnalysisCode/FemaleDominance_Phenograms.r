library(phytools)
library(OUwie)


# How does the percentage of fights won by females change relative to the social mating system?

percwon_females_phylodata_1<-combined
percwon_females_phylodata_1<-as.data.frame(percwon_females_phylodata_1 %>% group_by(corrected_species_id) %>% summarise(mean(perc_won_females,na.rm=T )) )
percwon_females_phylodata_2<-as.data.frame(percwon_females_phylodata_1[,2])
row.names(percwon_females_phylodata_2)<-percwon_females_phylodata_1[,1]
colnames(percwon_females_phylodata_2)<-"perc_won_females"
percwon_females_phylodata_3<-as.data.frame(percwon_females_phylodata_2[is.na(percwon_females_phylodata_2$perc_won_females)==FALSE,])
row.names(percwon_females_phylodata_3)<-row.names(percwon_females_phylodata_2)[is.na(percwon_females_phylodata_2$perc_won_females)==FALSE]
colnames(percwon_females_phylodata_3)<-"perc_won_females"
percwon_females_phylodata_3$corrected_species_id<-rownames(percwon_females_phylodata_3)

matingsystem<-select(specieslevelpredictors,MatSysPMK)
matingsystem<-as.data.frame(matingsystem)
matingsystem<-as.data.frame(matingsystem[is.na(matingsystem$MatSysPMK)==F,])
colnames(matingsystem)<-"MatSysPMK"
rownames(matingsystem)<-rownames(specieslevelpredictors[is.na(specieslevelpredictors$MatSysPMK)==F,])
matingsystem$corrected_species_id<-rownames(matingsystem)

percwon_females_phylodata<-inner_join(percwon_females_phylodata_3,matingsystem,by="corrected_species_id")
rownames(percwon_females_phylodata)<-percwon_females_phylodata$corrected_species_id


#  http://www.phytools.org/static.help/add.species.to.genus.html

# We match the species in the dataset to the species in the phylogenetic tree
missing<-treedata(inputtree,data=percwon_females_phylodata,warnings=FALSE)
# We remove species with no data from the tree
mtree<-missing$phy
speciesnames<-mtree$tip.label

mdata<-percwon_females_phylodata[percwon_females_phylodata$corrected_species_id %in% speciesnames,]


matingsystem_phylo<-mdata$MatSysPMK
names(matingsystem_phylo)<-mdata$corrected_species_id

matingsystemtree <- make.simmap(mtree, matingsystem_phylo, model="ARD")



perc_wondata<-mdata$perc_won_females
names(perc_wondata)<-mdata$corrected_species_id

phenogram(matingsystemtree,perc_wondata,ylab="percentage fights won by females")
legend(x="topleft",legend=c("Monogamy","Polyandry","Polygyny","Promiscuity"),lwd=3,col=c("black","#DF536B","#61D04F","#2297E6"),bty="n")

#   1 MON        2 PAN            3 POL              4 PRO
#  "black" "#DF536B"/pink "#61D04F"/green "#2297E6"/blue




# How does the percentage of fights won by females change relative to whether there is strong sexual dimorphism or not?

dimorphism<-combined
dimorphism<-dimorphism[is.na(dimorphism$SexualDimorphism_MaleWeight_over_FemaleWeight)==F,]
dimorphism[dimorphism$SexualDimorphism_MaleWeight_over_FemaleWeight<1.101,]$SexualDimorphism_MaleWeight_over_FemaleWeight<- 0
dimorphism[dimorphism$SexualDimorphism_MaleWeight_over_FemaleWeight>1.1,]$SexualDimorphism_MaleWeight_over_FemaleWeight<- 1
dimorphism_phylodata<-as.data.frame(dimorphism %>% group_by(corrected_species_id) %>% summarise(mean(SexualDimorphism_MaleWeight_over_FemaleWeight,na.rm=T )) )

percwon_females_phylodata_1<-combined
percwon_females_phylodata_1<-as.data.frame(percwon_females_phylodata_1 %>% group_by(corrected_species_id) %>% summarise(mean(perc_won_females,na.rm=T )) )
colnames(percwon_females_phylodata_1)<-c("corrected_species_id","perc_won_females")
percwon_females_phylodata_3<-as.data.frame(percwon_females_phylodata_1[is.na(percwon_females_phylodata_1$perc_won_females)==FALSE,])
row.names(percwon_females_phylodata_3)<-percwon_females_phylodata_3$corrected_species_id


percwon_females_phylodata<-inner_join(percwon_females_phylodata_3,dimorphism_phylodata,by="corrected_species_id")
rownames(percwon_females_phylodata)<-percwon_females_phylodata$corrected_species_id


#  http://www.phytools.org/static.help/add.species.to.genus.html

# We match the species in the dataset to the species in the phylogenetic tree
missing<-treedata(inputtree,data=percwon_females_phylodata,warnings=FALSE)
# We remove species with no data from the tree
mtree<-missing$phy
speciesnames<-mtree$tip.label

mdata<-percwon_females_phylodata[percwon_females_phylodata$corrected_species_id %in% speciesnames,]
colnames(mdata)<-c("corrected_species_id","perc_won_females","dimorphism")

dimorphism_phylo<-mdata$dimorphism
names(dimorphism_phylo)<-mdata$corrected_species_id

dimorphismtree <- make.simmap(mtree, dimorphism_phylo, model="ARD")



perc_wondata<-mdata$perc_won_females
names(perc_wondata)<-mdata$corrected_species_id

phenogram(dimorphismtree,perc_wondata,ylab="percentage fights won by females")
legend(x="topleft",legend=c("No body size dimorphism","Males more than 10% larger"),lwd=3,col=c("black","#DF536B"),bty="n")

#   not dimorphic        dimorphic
#  "black" "#DF536B"/pink



# How does the percentage of fights won by females change relative to whether species are arboreal or not?

arboreal<-select(combined,corrected_species_id,Strata_Wilman)
arboreal<-arboreal[is.na(arboreal$Strata_Wilman)==F,]
arboreal$arboreality<-0
arboreal[arboreal$Strata_Wilman=="Ar",]$arboreality<-1
arboreal_phylo<-as.data.frame(arboreal %>% group_by(corrected_species_id) %>% summarise(mean(arboreality,na.rm=T )) )
colnames(arboreal_phylo)<-c("corrected_species_id","arboreality")

percwon_females_phylodata_1<-combined
percwon_females_phylodata_1<-as.data.frame(percwon_females_phylodata_1 %>% group_by(corrected_species_id) %>% summarise(mean(perc_won_females,na.rm=T )) )
colnames(percwon_females_phylodata_1)<-c("corrected_species_id","perc_won_females")
percwon_females_phylodata_3<-as.data.frame(percwon_females_phylodata_1[is.na(percwon_females_phylodata_1$perc_won_females)==FALSE,])
row.names(percwon_females_phylodata_3)<-percwon_females_phylodata_3$corrected_species_id

percwon_females_phylodata<-inner_join(percwon_females_phylodata_3,arboreal_phylo,by="corrected_species_id")
rownames(percwon_females_phylodata)<-percwon_females_phylodata$corrected_species_id

#  http://www.phytools.org/static.help/add.species.to.genus.html

# We match the species in the dataset to the species in the phylogenetic tree
missing<-treedata(inputtree,data=percwon_females_phylodata,warnings=FALSE)
# We remove species with no data from the tree
mtree<-missing$phy
speciesnames<-mtree$tip.label

mdata<-percwon_females_phylodata[percwon_females_phylodata$corrected_species_id %in% speciesnames,]
colnames(mdata)<-c("corrected_species_id","perc_won_females","arboreality")

arbo_phylo<-mdata$arboreality
names(arbo_phylo)<-mdata$corrected_species_id

arbotree <- make.simmap(mtree, arbo_phylo, model="ARD")

perc_wondata<-mdata$perc_won_females
names(perc_wondata)<-mdata$corrected_species_id

pdf("Phenogram_PercFightsWon_Arboreality.pdf")
phenogram(arbotree,perc_wondata,ylab="percentage fights won by females")
legend(x="topleft",legend=c("Arboreal","Ground"),lwd=3,col=c("black","#DF536B"),bty="n")
dev.off()




# How does the group sex ratio change in relation to the dominance system?

dominance<-select(combined,corrected_species_id,strictfdom)
dominance<-dominance[is.na(dominance$strictfdom)==F,]
dominance$strictfdom<-as.numeric(dominance$strictfdom)
dominance_phylo<-as.data.frame(dominance %>% group_by(corrected_species_id) %>% summarise(mean(strictfdom,na.rm=T )) )
colnames(dominance_phylo)<-c("corrected_species_id","dominance")
dominance_phylo$dominance<-round(dominance_phylo$dominance,0)

sexratio_females_phylodata_1<-combined
sexratio_females_phylodata_1<-as.data.frame(sexratio_females_phylodata_1 %>% group_by(corrected_species_id) %>% summarise(mean(sexratio,na.rm=T )) )
colnames(sexratio_females_phylodata_1)<-c("corrected_species_id","sexratio")
sexratio_females_phylodata_3<-as.data.frame(sexratio_females_phylodata_1[is.na(sexratio_females_phylodata_1$sexratio)==FALSE,])
row.names(sexratio_females_phylodata_3)<-sexratio_females_phylodata_3$corrected_species_id

sexratio_females_phylodata<-inner_join(sexratio_females_phylodata_3,dominance_phylo,by="corrected_species_id")
rownames(sexratio_females_phylodata)<-sexratio_females_phylodata$corrected_species_id

#  http://www.phytools.org/static.help/add.species.to.genus.html

# We match the species in the dataset to the species in the phylogenetic tree
missing<-treedata(inputtree,data=sexratio_females_phylodata,warnings=FALSE)
# We remove species with no data from the tree
mtree<-missing$phy
speciesnames<-mtree$tip.label

mdata<-sexratio_females_phylodata[sexratio_females_phylodata$corrected_species_id %in% speciesnames,]
colnames(mdata)<-c("corrected_species_id","sexratio","dominance")

dominance_phylo<-mdata$dominance
names(dominance_phylo)<-mdata$corrected_species_id

dominancetree <- make.simmap(mtree, dominance_phylo, model="ARD")

sexratiodata<-mdata$sexratio
names(sexratiodata)<-mdata$corrected_species_id

pdf("Phenogram_SexRatio_Dominance.pdf")
phenogram(dominancetree,sexratiodata,ylab="sex ratio")
legend(x="topleft",legend=c("Male dominance","Co dominance","Female dominance"),lwd=3,col=c("black","#DF536B","#61D04F"),bty="n")
dev.off()


# How does the number of females change in relation to the dominance system?

dominance<-select(combined,corrected_species_id,strictfdom)
dominance<-dominance[is.na(dominance$strictfdom)==F,]
dominance$strictfdom<-as.numeric(dominance$strictfdom)
dominance_phylo<-as.data.frame(dominance %>% group_by(corrected_species_id) %>% summarise(mean(strictfdom,na.rm=T )) )
colnames(dominance_phylo)<-c("corrected_species_id","dominance")
dominance_phylo$dominance<-round(dominance_phylo$dominance,0)

females_phylodata_1<-combined
females_phylodata_1<-as.data.frame(females_phylodata_1 %>% group_by(corrected_species_id) %>% summarise(mean(females,na.rm=T )) )
colnames(females_phylodata_1)<-c("corrected_species_id","females")
females_phylodata_3<-as.data.frame(females_phylodata_1[is.na(females_phylodata_1$females)==FALSE,])
row.names(females_phylodata_3)<-females_phylodata_3$corrected_species_id

females_phylodata<-inner_join(females_phylodata_3,dominance_phylo,by="corrected_species_id")
rownames(females_phylodata)<-females_phylodata$corrected_species_id

#  http://www.phytools.org/static.help/add.species.to.genus.html

# We match the species in the dataset to the species in the phylogenetic tree
missing<-treedata(inputtree,data=females_phylodata,warnings=FALSE)
# We remove species with no data from the tree
mtree<-missing$phy
speciesnames<-mtree$tip.label

mdata<-females_phylodata[females_phylodata$corrected_species_id %in% speciesnames,]
colnames(mdata)<-c("corrected_species_id","females","dominance")

dominance_phylo<-mdata$dominance
names(dominance_phylo)<-mdata$corrected_species_id

dominancetree <- make.simmap(mtree, dominance_phylo, model="ARD")

femalesdata<-mdata$females
names(femalesdata)<-mdata$corrected_species_id

pdf("Phenogram_Females_Dominance.pdf")
phenogram(dominancetree,femalesdata,ylab="number of females")
legend(x="topleft",legend=c("Male dominance","Co dominance","Female dominance"),lwd=3,col=c("black","#DF536B","#61D04F"),bty="n")
dev.off()


# How does the sexual dimorphism change in relation to the dominance system?

dominance<-select(combined,corrected_species_id,strictfdom)
dominance<-dominance[is.na(dominance$strictfdom)==F,]
dominance$strictfdom<-as.numeric(dominance$strictfdom)
dominance_phylo<-as.data.frame(dominance %>% group_by(corrected_species_id) %>% summarise(mean(strictfdom,na.rm=T )) )
colnames(dominance_phylo)<-c("corrected_species_id","dominance")
dominance_phylo$dominance<-round(dominance_phylo$dominance,0)

dimorphism_phylodata_1<-combined
dimorphism_phylodata_1<-as.data.frame(dimorphism_phylodata_1 %>% group_by(corrected_species_id) %>% summarise(mean(SexualDimorphism_MaleWeight_over_FemaleWeight,na.rm=T )) )
colnames(dimorphism_phylodata_1)<-c("corrected_species_id","SexualDimorphism_MaleWeight_over_FemaleWeight")
dimorphism_phylodata_3<-as.data.frame(dimorphism_phylodata_1[is.na(dimorphism_phylodata_1$SexualDimorphism_MaleWeight_over_FemaleWeight)==FALSE,])
row.names(dimorphism_phylodata_3)<-dimorphism_phylodata_3$corrected_species_id

dimorphism_phylodata<-inner_join(dimorphism_phylodata_3,dominance_phylo,by="corrected_species_id")
rownames(dimorphism_phylodata)<-dimorphism_phylodata$corrected_species_id

#  http://www.phytools.org/static.help/add.species.to.genus.html

# We match the species in the dataset to the species in the phylogenetic tree
missing<-treedata(inputtree,data=dimorphism_phylodata,warnings=FALSE)
# We remove species with no data from the tree
mtree<-missing$phy
speciesnames<-mtree$tip.label

mdata<-dimorphism_phylodata[dimorphism_phylodata$corrected_species_id %in% speciesnames,]
colnames(mdata)<-c("corrected_species_id","SexualDimorphism_MaleWeight_over_FemaleWeight","dominance")

dominance_phylo<-mdata$dominance
names(dominance_phylo)<-mdata$corrected_species_id

dominancetree <- make.simmap(mtree, dominance_phylo, model="ARD")

dimorphismdata<-mdata$SexualDimorphism_MaleWeight_over_FemaleWeight
names(dimorphismdata)<-mdata$corrected_species_id

pdf("Phenogram_Dimorphism_Dominance.pdf")
phenogram(dominancetree,dimorphismdata,ylab="body size dimorphism")
legend(x="topleft",legend=c("Male dominance","Co dominance","Female dominance"),lwd=3,col=c("black","#DF536B","#61D04F"),bty="n")
dev.off()


# How does the sexual receptivity change in relation to the dominance system?

dominance<-select(combined,corrected_species_id,strictfdom)
dominance<-dominance[is.na(dominance$strictfdom)==F,]
dominance$strictfdom<-as.numeric(dominance$strictfdom)
dominance_phylo<-as.data.frame(dominance %>% group_by(corrected_species_id) %>% summarise(mean(strictfdom,na.rm=T )) )
colnames(dominance_phylo)<-c("corrected_species_id","dominance")
dominance_phylo$dominance<-round(dominance_phylo$dominance,0)

receptivity_phylodata_1<-combined
receptivity_phylodata_1<-as.data.frame(receptivity_phylodata_1 %>% group_by(corrected_species_id) %>% summarise(mean(sexualreceptivity_hours,na.rm=T )) )
colnames(receptivity_phylodata_1)<-c("corrected_species_id","sexualreceptivity_hours")
receptivity_phylodata_3<-as.data.frame(receptivity_phylodata_1[is.na(receptivity_phylodata_1$sexualreceptivity_hours)==FALSE,])
row.names(receptivity_phylodata_3)<-receptivity_phylodata_3$corrected_species_id

receptivity_phylodata<-inner_join(receptivity_phylodata_3,dominance_phylo,by="corrected_species_id")
rownames(receptivity_phylodata)<-receptivity_phylodata$corrected_species_id

#  http://www.phytools.org/static.help/add.species.to.genus.html

# We match the species in the dataset to the species in the phylogenetic tree
missing<-treedata(inputtree,data=receptivity_phylodata,warnings=FALSE)
# We remove species with no data from the tree
mtree<-missing$phy
speciesnames<-mtree$tip.label

mdata<-receptivity_phylodata[receptivity_phylodata$corrected_species_id %in% speciesnames,]
colnames(mdata)<-c("corrected_species_id","sexualreceptivity_hours","dominance")

dominance_phylo<-mdata$dominance
names(dominance_phylo)<-mdata$corrected_species_id

dominancetree <- make.simmap(mtree, dominance_phylo, model="ARD")

receptivitydata<-mdata$sexualreceptivity_hours
names(receptivitydata)<-mdata$corrected_species_id

pdf("Phenogram_Receptivity_Dominance.pdf")
phenogram(dominancetree,receptivitydata,ylab="sexual receptivity (hours)")
legend(x="topleft",legend=c("Male dominance","Co dominance","Female dominance"),lwd=3,col=c("black","#DF536B","#61D04F"),bty="n")
dev.off()


# How does the receptive synchrony in relation to the dominance system?

dominance<-select(combined,corrected_species_id,strictfdom)
dominance<-dominance[is.na(dominance$strictfdom)==F,]
dominance$strictfdom<-as.numeric(dominance$strictfdom)
dominance_phylo<-as.data.frame(dominance %>% group_by(corrected_species_id) %>% summarise(mean(strictfdom,na.rm=T )) )
colnames(dominance_phylo)<-c("corrected_species_id","dominance")
dominance_phylo$dominance<-round(dominance_phylo$dominance,0)

receptivity_phylodata_1<-combined
receptivity_phylodata_1<-as.data.frame(receptivity_phylodata_1 %>% group_by(corrected_species_id) %>% summarise(mean(receptive_synchrony,na.rm=T )) )
colnames(receptivity_phylodata_1)<-c("corrected_species_id","receptive_synchrony")
receptivity_phylodata_3<-as.data.frame(receptivity_phylodata_1[is.na(receptivity_phylodata_1$receptive_synchrony)==FALSE,])
row.names(receptivity_phylodata_3)<-receptivity_phylodata_3$corrected_species_id

receptivity_phylodata<-inner_join(receptivity_phylodata_3,dominance_phylo,by="corrected_species_id")
rownames(receptivity_phylodata)<-receptivity_phylodata$corrected_species_id

#  http://www.phytools.org/static.help/add.species.to.genus.html

# We match the species in the dataset to the species in the phylogenetic tree
missing<-treedata(inputtree,data=receptivity_phylodata,warnings=FALSE)
# We remove species with no data from the tree
mtree<-missing$phy
speciesnames<-mtree$tip.label

mdata<-receptivity_phylodata[receptivity_phylodata$corrected_species_id %in% speciesnames,]
colnames(mdata)<-c("corrected_species_id","receptive_synchrony","dominance")

dominance_phylo<-mdata$dominance
names(dominance_phylo)<-mdata$corrected_species_id

dominancetree <- make.simmap(mtree, dominance_phylo, model="ARD")

receptivitydata<-mdata$receptive_synchrony
names(receptivitydata)<-mdata$corrected_species_id

pdf("Phenogram_ReceptiveSynchrony_Dominance.pdf")
phenogram(dominancetree,receptivitydata,ylab="receptive synchrony")
legend(x="topleft",legend=c("Male dominance","Co dominance","Female dominance"),lwd=3,col=c("black","#DF536B","#61D04F"),bty="n")
dev.off()




# How does the seasonality in relation to the dominance system?

dominance<-select(combined,corrected_species_id,strictfdom)
dominance<-dominance[is.na(dominance$strictfdom)==F,]
dominance$strictfdom<-as.numeric(dominance$strictfdom)
dominance_phylo<-as.data.frame(dominance %>% group_by(corrected_species_id) %>% summarise(mean(strictfdom,na.rm=T )) )
colnames(dominance_phylo)<-c("corrected_species_id","dominance")
dominance_phylo$dominance<-round(dominance_phylo$dominance,0)

seasonality_phylodata_1<-combined
seasonality_phylodata_1<-as.data.frame(seasonality_phylodata_1 %>% group_by(corrected_species_id) %>% summarise(mean(r_seasonality_value,na.rm=T )) )
colnames(seasonality_phylodata_1)<-c("corrected_species_id","r_seasonality_value")
seasonality_phylodata_3<-as.data.frame(seasonality_phylodata_1[is.na(seasonality_phylodata_1$r_seasonality_value)==FALSE,])
row.names(seasonality_phylodata_3)<-seasonality_phylodata_3$corrected_species_id

seasonality_phylodata<-inner_join(seasonality_phylodata_3,dominance_phylo,by="corrected_species_id")
rownames(seasonality_phylodata)<-seasonality_phylodata$corrected_species_id

#  http://www.phytools.org/static.help/add.species.to.genus.html

# We match the species in the dataset to the species in the phylogenetic tree
missing<-treedata(inputtree,data=seasonality_phylodata,warnings=FALSE)
# We remove species with no data from the tree
mtree<-missing$phy
speciesnames<-mtree$tip.label

mdata<-seasonality_phylodata[seasonality_phylodata$corrected_species_id %in% speciesnames,]
colnames(mdata)<-c("corrected_species_id","r_seasonality_value","dominance")

dominance_phylo<-mdata$dominance
names(dominance_phylo)<-mdata$corrected_species_id

dominancetree <- make.simmap(mtree, dominance_phylo, model="ARD")

seasonalitydata<-mdata$r_seasonality_value
names(seasonalitydata)<-mdata$corrected_species_id

pdf("Phenogram_ReceptiveSynchrony_Dominance.pdf")
phenogram(dominancetree,seasonalitydata,ylab="environmental seasonality")
legend(x="topleft",legend=c("Male dominance","Co dominance","Female dominance"),lwd=3,col=c("black","#DF536B","#61D04F"),bty="n")
dev.off()



# How does the canine dimorphism in relation to the dominance system?

dominance<-select(combined,corrected_species_id,strictfdom)
dominance<-dominance[is.na(dominance$strictfdom)==F,]
dominance$strictfdom<-as.numeric(dominance$strictfdom)
dominance_phylo<-as.data.frame(dominance %>% group_by(corrected_species_id) %>% summarise(mean(strictfdom,na.rm=T )) )
colnames(dominance_phylo)<-c("corrected_species_id","dominance")
dominance_phylo$dominance<-round(dominance_phylo$dominance,0)

canines_phylodata_1<-combined
canines_phylodata_1<-as.data.frame(canines_phylodata_1 %>% group_by(corrected_species_id) %>% summarise(mean(CanineDimorphism,na.rm=T )) )
colnames(canines_phylodata_1)<-c("corrected_species_id","CanineDimorphism")
canines_phylodata_3<-as.data.frame(canines_phylodata_1[is.na(canines_phylodata_1$CanineDimorphism)==FALSE,])
row.names(canines_phylodata_3)<-canines_phylodata_3$corrected_species_id

canine_phylodata<-inner_join(canines_phylodata_3,dominance_phylo,by="corrected_species_id")
rownames(canine_phylodata)<-canine_phylodata$corrected_species_id

#  http://www.phytools.org/static.help/add.species.to.genus.html

# We match the species in the dataset to the species in the phylogenetic tree
missing<-treedata(inputtree,data=canine_phylodata,warnings=FALSE)
# We remove species with no data from the tree
mtree<-missing$phy
speciesnames<-mtree$tip.label

mdata<-canine_phylodata[canine_phylodata$corrected_species_id %in% speciesnames,]
colnames(mdata)<-c("corrected_species_id","CanineDimorphism","dominance")

dominance_phylo<-mdata$dominance
names(dominance_phylo)<-mdata$corrected_species_id

dominancetree <- make.simmap(mtree, dominance_phylo, model="ARD")

caninedata<-mdata$CanineDimorphism
names(caninedata)<-mdata$corrected_species_id

pdf("Phenogram_CanineDimorphism_Dominance.pdf")
phenogram(dominancetree,caninedata,ylab="canine size dimorphism")
legend(x="topleft",legend=c("Male dominance","Co dominance","Female dominance"),lwd=3,col=c("black","#DF536B","#61D04F"),bty="n")
dev.off()


# How does the body size change in relation to the dominance system?

dominance<-select(combined,corrected_species_id,strictfdom)
dominance<-dominance[is.na(dominance$strictfdom)==F,]
dominance$strictfdom<-as.numeric(dominance$strictfdom)
dominance_phylo<-as.data.frame(dominance %>% group_by(corrected_species_id) %>% summarise(mean(strictfdom,na.rm=T )) )
colnames(dominance_phylo)<-c("corrected_species_id","dominance")
dominance_phylo$dominance<-round(dominance_phylo$dominance,0)

body_phylodata_1<-combined
body_phylodata_1<-as.data.frame(body_phylodata_1 %>% group_by(corrected_species_id) %>% summarise(mean(body_mass,na.rm=T )) )
colnames(body_phylodata_1)<-c("corrected_species_id","body_mass")
body_phylodata_3<-as.data.frame(body_phylodata_1[is.na(body_phylodata_1$body_mass)==FALSE,])
row.names(body_phylodata_3)<-body_phylodata_3$corrected_species_id

body_phylodata<-inner_join(body_phylodata_3,dominance_phylo,by="corrected_species_id")
rownames(body_phylodata)<-body_phylodata$corrected_species_id

#  http://www.phytools.org/static.help/add.species.to.genus.html

# We match the species in the dataset to the species in the phylogenetic tree
missing<-treedata(inputtree,data=body_phylodata,warnings=FALSE)
# We remove species with no data from the tree
mtree<-missing$phy
speciesnames<-mtree$tip.label

mdata<-body_phylodata[body_phylodata$corrected_species_id %in% speciesnames,]
colnames(mdata)<-c("corrected_species_id","body_mass","dominance")

dominance_phylo<-mdata$dominance
names(dominance_phylo)<-mdata$corrected_species_id

dominancetree <- make.simmap(mtree, dominance_phylo, model="ARD")

bodydata<-mdata$body_mass
names(bodydata)<-mdata$corrected_species_id

pdf("Phenogram_BodySize_Dominance.pdf")
phenogram(dominancetree,bodydata,ylab="body size")
legend(x="topleft",legend=c("Male dominance","Co dominance","Female dominance"),lwd=3,col=c("black","#DF536B","#61D04F"),bty="n")
dev.off()



# How does the testes size in relation to the dominance system?

dominance<-select(combined,corrected_species_id,strictfdom)
dominance<-dominance[is.na(dominance$strictfdom)==F,]
dominance$strictfdom<-as.numeric(dominance$strictfdom)
dominance_phylo<-as.data.frame(dominance %>% group_by(corrected_species_id) %>% summarise(mean(strictfdom,na.rm=T )) )
colnames(dominance_phylo)<-c("corrected_species_id","dominance")
dominance_phylo$dominance<-round(dominance_phylo$dominance,0)

body_phylodata_1<-combined
body_phylodata_1$reltestesmass<-body_phylodata_1$testes_mass/body_phylodata_1$body_mass
body_phylodata_1<-as.data.frame(body_phylodata_1 %>% group_by(corrected_species_id) %>% summarise(mean(reltestesmass,na.rm=T )) )
colnames(body_phylodata_1)<-c("corrected_species_id","reltestesmass")
body_phylodata_3<-as.data.frame(body_phylodata_1[is.na(body_phylodata_1$reltestesmass)==FALSE,])
row.names(body_phylodata_3)<-body_phylodata_3$corrected_species_id

body_phylodata<-inner_join(body_phylodata_3,dominance_phylo,by="corrected_species_id")
rownames(body_phylodata)<-body_phylodata$corrected_species_id

#  http://www.phytools.org/static.help/add.species.to.genus.html

# We match the species in the dataset to the species in the phylogenetic tree
missing<-treedata(inputtree,data=body_phylodata,warnings=FALSE)
# We remove species with no data from the tree
mtree<-missing$phy
speciesnames<-mtree$tip.label

mdata<-body_phylodata[body_phylodata$corrected_species_id %in% speciesnames,]
colnames(mdata)<-c("corrected_species_id","reltestesmass","dominance")

dominance_phylo<-mdata$dominance
names(dominance_phylo)<-mdata$corrected_species_id

dominancetree <- make.simmap(mtree, dominance_phylo, model="ARD")

bodydata<-mdata$reltestesmass
names(bodydata)<-mdata$corrected_species_id

pdf("Phenogram_TestesSize_Dominance.pdf")
phenogram(dominancetree,bodydata,ylab="relative testes size dimorphism")
legend(x="topleft",legend=c("Male dominance","Co dominance","Female dominance"),lwd=3,col=c("black","#DF536B","#61D04F"),bty="n")
dev.off()