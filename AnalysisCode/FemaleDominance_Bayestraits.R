
# Phylogenetic reconstruction of the evolution of intersexual dominance using Bayestraits

# Analyses using the strict three way classification of dominance

############################################################

# Multistate model investigating the transitions between male-, co-, and female- dominance

strict_females_phylodata_1<-combined
strict_females_phylodata_1$strictfdom<-as.numeric(as.factor(combined$strictfdom))
strict_females_phylodata_1<-as.data.frame(strict_females_phylodata_1 %>% group_by(corrected_species_id) %>% summarise(mean(strictfdom,na.rm=T )) )
colnames(strict_females_phylodata_1)<-c("Species","Strictfdom")
rownames(strict_females_phylodata_1)<-strict_females_phylodata_1$Species

missing<-treedata(inputtree,data=strict_females_phylodata_1,warnings=FALSE)
# We remove species with no data from the tree
mtree<-missing$phy
speciesnames<-mtree$tip.label

mdata<-strict_females_phylodata_1[strict_females_phylodata_1$Species %in% speciesnames,]
mdata$Strictfdom<-round(mdata$Strictfdom)

bayestraitsdiscretedata<-matrix(NA,nrow=nrow(mdata),ncol=2)
colnames(bayestraitsdiscretedata)<-c("Species","Strictfdom")
rownames(bayestraitsdiscretedata)<-mtree$tip.label
bayestraitsdiscretedata[,1]<-mtree$tip.label
bayestraitsdiscretedata[,2]<- as.numeric(as.factor(mdata$Strictfdom))
bayestraitsdiscretedata<-as.data.frame(bayestraitsdiscretedata)

command_vec1 <- c("1", "1","MLTries = 100") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
results_1 <- bayestraits(bayestraitsdiscretedata, mtree, command_vec1)
log_1 <- results_1$Log$results

colnames(log_1)<-c("Tree","Likelihood","Transition Male to Co","Transition Male to Female","Transition Co to Male","Transition Co to Female","Transition Female to Male","Transition Female to Co","Likelihood Root is Male","Likelihood Root is Co","Likelihood Root is Female")

log_1$Count_MaleDominance<-table(bayestraitsdiscretedata$Strictfdom)[1]
log_1$Count_CoDominance<-table(bayestraitsdiscretedata$Strictfdom)[2]
log_1$Count_FemaleDominance<-table(bayestraitsdiscretedata$Strictfdom)[3]

write.csv(log_1,file="results/Bayestraits_Transitions_StrictIntersexDominance.csv")

# Evolutionary transition always pass through the intermediate state
# Roughly equal chances of all other transitions (to female dominance lowest)
# Distribution is 32 male dominance, 66 egalitarian, 19 female dominance species





############################################################


# Testing whether there is dependent evolution of monogamy and female dominance

mostly_females_phylodata_1<-combined
mostly_females_phylodata_1$female_dominance<-as.numeric(as.factor(combined$strictfdom))
mostly_females_phylodata_1<-as.data.frame(mostly_females_phylodata_1 %>% group_by(corrected_species_id) %>% summarise(mean(female_dominance,na.rm=T )) )
colnames(mostly_females_phylodata_1)<-c("corrected_species_id", "femaledominance")
mostly_females_phylodata_1$femaledominance<-as.numeric(mostly_females_phylodata_1$femaledominance>2.4)

monogamy<-combined
monogamy$MatSysPMK<-as.numeric(as.factor(monogamy$MatSysPMK))
monogamy<-monogamy[is.na(monogamy$MatSysPMK)==F,]
monogamy[monogamy$MatSysPMK>1,]$MatSysPMK<-0
monogamy_phylodata<-as.data.frame(monogamy %>% group_by(corrected_species_id) %>% summarise(mean(MatSysPMK,na.rm=T )) )
mostly_females_phylodata<-left_join(mostly_females_phylodata_1,monogamy_phylodata,by="corrected_species_id")
rownames(mostly_females_phylodata)<-mostly_females_phylodata$corrected_species_id
colnames(mostly_females_phylodata)<-c("Species","FemaleDominance","Monogamy")
mostly_females_phylodata<-mostly_females_phylodata[is.na(mostly_females_phylodata$Monogamy)==F,]
mostly_females_phylodata<-mostly_females_phylodata[is.na(mostly_females_phylodata$FemaleDominance)==F,]


# We match the species in the dataset to the species in the phylogenetic tree
missing<-treedata(inputtree,data=mostly_females_phylodata,warnings=FALSE)
# We remove species with no data from the tree
mtree<-missing$phy
speciesnames<-mtree$tip.label

mdata<-mostly_females_phylodata[speciesnames,]

bayestraitsdiscretedata<-matrix(NA,nrow=nrow(mdata),ncol=3)
colnames(bayestraitsdiscretedata)<-c("Species","Monogamy","FemaleDominance")
rownames(bayestraitsdiscretedata)<-mtree$tip.label
bayestraitsdiscretedata[,1]<-mtree$tip.label
bayestraitsdiscretedata[,2]<-as.numeric(as.factor(mdata$Monogamy))-1
bayestraitsdiscretedata[,3]<- as.numeric(as.factor(mdata$FemaleDominance))-1
bayestraitsdiscretedata<-as.data.frame(bayestraitsdiscretedata)

command_vec1 <- c("2", "1","mltries = 100") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
results_1 <- bayestraits(bayestraitsdiscretedata, mtree, command_vec1)
log_1 <- results_1$Log$results


command_vec2 <- c("3", "1","mltries = 100") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
results_2 <- bayestraits(bayestraitsdiscretedata, mtree, command_vec2)
log_2 <- results_2$Log$results

# 1 = 0,0 not monogamous, not female dominance - 2 = 0,1 not monogamous, female dominance   - 3 = 1,0 monogamous, not female dominance   4 = 1,1 monogamous, female dominance
# Transition to female dominance 1-2 and 3-4: more likely when monogamy present (15x more)
# Transition to monogamy 1-3 and 2-4: more likely when females not dominant (5x more)
# might suggest that monogamy comes first

colnames(log_2)<-c("Tree","Likelihood","TransitionToFemaleDominance_AbsenceMonogamy","TransitionToMonogamy_AbsenceFemaleDominance","LossFemaleDominance_AbsenceMonogamy","TransitionToMonogamy_PresenceFemaleDominance","LossMonogamy_AbsenceFemaleDominance","TransitionToFemaleDominance_PresenceMonogamy","LossMonogamy_PresenceFemaleDominance","LossFemaleDominance_PresenceMonogamy","LikelihoodRootNotMonogamyNotFemaleDominance","LikelihoodRootNotMonogamyFemaleDominance","LikelihoodRootMonogamyNotFemaleDominance","LikelihoodRootMonogamyFemaleDominance")


log_2$Count_MaleAndCoDominance<-table(bayestraitsdiscretedata$FemaleDominance)[1]
log_2$Count_FemaleDominance<-table(bayestraitsdiscretedata$FemaleDominance)[2]
log_2$NotMonogamous<-table(bayestraitsdiscretedata$Monogamy)[1]
log_2$Monogamous<-table(bayestraitsdiscretedata$Monogamy)[2]

log_2$Positive_IndependentModelBetter<-2*(log_2$Likelihood-log_1$Lh)

write.csv(log_2,file="results/Bayestraits_Dependent_StrictFemaleDominanceMonogamy.csv")



# Testing whether there is dependent evolution of monogamy and male dominance

mostly_females_phylodata_1<-combined
mostly_females_phylodata_1$female_dominance<-as.numeric(as.factor(combined$strictfdom))
mostly_females_phylodata_1<-as.data.frame(mostly_females_phylodata_1 %>% group_by(corrected_species_id) %>% summarise(mean(female_dominance,na.rm=T )) )
colnames(mostly_females_phylodata_1)<-c("corrected_species_id", "femaledominance")
mostly_females_phylodata_1$femaledominance<-as.numeric(mostly_females_phylodata_1$femaledominance<1.6)

monogamy<-combined
monogamy$MatSysPMK<-as.numeric(as.factor(monogamy$MatSysPMK))
monogamy<-monogamy[is.na(monogamy$MatSysPMK)==F,]
monogamy[monogamy$MatSysPMK>1,]$MatSysPMK<-0
monogamy_phylodata<-as.data.frame(monogamy %>% group_by(corrected_species_id) %>% summarise(mean(MatSysPMK,na.rm=T )) )
mostly_females_phylodata<-left_join(mostly_females_phylodata_1,monogamy_phylodata,by="corrected_species_id")
rownames(mostly_females_phylodata)<-mostly_females_phylodata$corrected_species_id
colnames(mostly_females_phylodata)<-c("Species","FemaleDominance","Monogamy")
mostly_females_phylodata<-mostly_females_phylodata[is.na(mostly_females_phylodata$Monogamy)==F,]
mostly_females_phylodata<-mostly_females_phylodata[is.na(mostly_females_phylodata$FemaleDominance)==F,]


# We match the species in the dataset to the species in the phylogenetic tree
missing<-treedata(inputtree,data=mostly_females_phylodata,warnings=FALSE)
# We remove species with no data from the tree
mtree<-missing$phy
speciesnames<-mtree$tip.label

mdata<-mostly_females_phylodata[speciesnames,]

bayestraitsdiscretedata<-matrix(NA,nrow=nrow(mdata),ncol=3)
colnames(bayestraitsdiscretedata)<-c("Species","Monogamy","MaleDominance")
rownames(bayestraitsdiscretedata)<-mtree$tip.label
bayestraitsdiscretedata[,1]<-mtree$tip.label
bayestraitsdiscretedata[,2]<-as.numeric(as.factor(mdata$Monogamy))-1
bayestraitsdiscretedata[,3]<- as.numeric(as.factor(mdata$FemaleDominance))-1
bayestraitsdiscretedata<-as.data.frame(bayestraitsdiscretedata)

command_vec1 <- c("2", "1","mltries = 100") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
results_1 <- bayestraits(bayestraitsdiscretedata, mtree, command_vec1)
log_1 <- results_1$Log$results


command_vec2 <- c("3", "1","mltries = 100") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
results_2 <- bayestraits(bayestraitsdiscretedata, mtree, command_vec2)
log_2 <- results_2$Log$results

# 1 = 0,0 not monogamous, not male dominance - 2 = 0,1 not monogamous, male dominance   - 3 = 1,0 monogamous, not male dominance   4 = 1,1 monogamous, male dominance
# Transition to male dominance 1-2 and 3-4: only when not monogamous
# Transition to monogamy 1-3 and 2-4: only when not male dominance

colnames(log_2)<-c("Tree","Likelihood","TransitionToMaleDominance_AbsenceMonogamy","TransitionToMonogamy_AbsenceMaleDominance","LossMaleDominance_AbsenceMonogamy","TransitionToMonogamy_PresenceMaleDominance","LossMonogamy_AbsenceMaleDominance","TransitionToMaleDominance_PresenceMonogamy","LossMonogamy_PresenceMaleDominance","LossMaleDominance_PresenceMonogamy","LikelihoodRootNotMonogamyNotMaleDominance","LikelihoodRootNotMonogamyMaleDominance","LikelihoodRootMonogamyNotMaleDominance","LikelihoodRootMonogamyFemaleDominance")

log_2$Count_FemaleAndCoDominance<-table(bayestraitsdiscretedata$FemaleDominance)[1]
log_2$Count_MaleDominance<-table(bayestraitsdiscretedata$FemaleDominance)[2]
log_2$NotMonogamous<-table(bayestraitsdiscretedata$Monogamy)[1]
log_2$Monogamous<-table(bayestraitsdiscretedata$Monogamy)[2]

log_2$Positive_IndependentModelBetter<-2*(log_2$Likelihood-log_1$Lh)

write.csv(log_2,file="results/Bayestraits_Dependent_StrictMaleDominanceMonogamy.csv")



# Testing whether there is dependent evolution of polygyny and male dominance

mostly_females_phylodata_1<-combined
mostly_females_phylodata_1$female_dominance<-as.numeric(as.factor(combined$strictfdom))
mostly_females_phylodata_1<-as.data.frame(mostly_females_phylodata_1 %>% group_by(corrected_species_id) %>% summarise(mean(female_dominance,na.rm=T )) )
colnames(mostly_females_phylodata_1)<-c("corrected_species_id", "femaledominance")
mostly_females_phylodata_1$femaledominance<-as.numeric(mostly_females_phylodata_1$femaledominance<1.6)

polygyny<-combined
polygyny$MatSysPMK<-as.numeric(combined$MatSysPMK=="POL")
polygyny<-polygyny[is.na(polygyny$MatSysPMK)==F,]
polygyny$MatSysPMK<-polygyny$MatSysPMK
polygyny_phylodata<-as.data.frame(polygyny %>% group_by(corrected_species_id) %>% summarise(mean(MatSysPMK,na.rm=T )) )
mostly_females_phylodata<-left_join(mostly_females_phylodata_1,polygyny_phylodata,by="corrected_species_id")
rownames(mostly_females_phylodata)<-mostly_females_phylodata$corrected_species_id
colnames(mostly_females_phylodata)<-c("Species","MaleDominance","Polygyny")
mostly_females_phylodata<-mostly_females_phylodata[is.na(mostly_females_phylodata$Polygyny)==F,]
mostly_females_phylodata<-mostly_females_phylodata[is.na(mostly_females_phylodata$MaleDominance)==F,]


# We match the species in the dataset to the species in the phylogenetic tree
missing<-treedata(inputtree,data=mostly_females_phylodata,warnings=FALSE)
# We remove species with no data from the tree
mtree<-missing$phy
speciesnames<-mtree$tip.label

mdata<-mostly_females_phylodata[speciesnames,]

bayestraitsdiscretedata<-matrix(NA,nrow=nrow(mdata),ncol=3)
colnames(bayestraitsdiscretedata)<-c("Species","Polygyny","MaleDominance")
rownames(bayestraitsdiscretedata)<-mtree$tip.label
bayestraitsdiscretedata[,1]<-mtree$tip.label
bayestraitsdiscretedata[,2]<-as.numeric(as.factor(mdata$Polygyny))-1
bayestraitsdiscretedata[,3]<- as.numeric(as.factor(mdata$Male))-1
bayestraitsdiscretedata<-as.data.frame(bayestraitsdiscretedata)

command_vec1 <- c("2", "1","mltries = 100") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
results_1 <- bayestraits(bayestraitsdiscretedata, mtree, command_vec1)
log_1 <- results_1$Log$results


command_vec2 <- c("3", "1","mltries = 100") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
results_2 <- bayestraits(bayestraitsdiscretedata, mtree, command_vec2)
log_2 <- results_2$Log$results

# 1 = 0,0 not polygynous, not male dominance - 2 = 0,1 not polygynous, male dominance   - 3 = 1,0 polygynous, not male dominance   4 = 1,1 polygynous, male dominance
# Transition to male dominance 1-2 and 3-4: only when not polygynous
# Transition to polygyny 1-3 and 2-4: only when not male dominance

log_2$Count_FemaleAndCoDominance<-table(bayestraitsdiscretedata$MaleDominance)[1]
log_2$Count_MaleDominance<-table(bayestraitsdiscretedata$MaleDominance)[2]
log_2$NotMonogamous<-table(bayestraitsdiscretedata$Polygyny)[1]
log_2$Monogamous<-table(bayestraitsdiscretedata$Polygyny)[2]

colnames(log_2)<-c("Tree","Likelihood","TransitionToMaleDominance_AbsencePolygyny","TransitionToPolygyny_AbsenceMaleDominance","LossMaleDominance_AbsencePolygyny","TransitionToPolygyny_PresenceMaleDominance","LossPolygyny_AbsenceMaleDominance","TransitionToMaleDominance_PresencePolygyny","LossPolygyny_PresenceMaleDominance","LossMaleDominance_PresencePolygyny","LikelihoodRootNotPolygynyNotMaleDominance","LikelihoodRootNotPolygynyMaleDominance","LikelihoodRootPolygynyNotMaleDominance","LikelihoodRootPolygynyMaleDominance","CountNotMaleDominant","CountMaleDominant","CountNotPolygynous","CountPolygynous")

log_2$Positive_IndependentModelBetter<-2*(log_2$Likelihood-log_1$Lh)

write.csv(log_2,file="results/Bayestraits_Dependent_StrictMaleDominancePolygyny.csv")




# Testing whether there is dependent evolution of single-male groups and male dominance

mostly_females_phylodata_1<-combined
mostly_females_phylodata_1$female_dominance<-as.numeric(as.factor(combined$strictfdom))
mostly_females_phylodata_1<-as.data.frame(mostly_females_phylodata_1 %>% group_by(corrected_species_id) %>% summarise(mean(female_dominance,na.rm=T )) )
colnames(mostly_females_phylodata_1)<-c("corrected_species_id", "femaledominance")
mostly_females_phylodata_1$femaledominance<-as.numeric(mostly_females_phylodata_1$femaledominance<1.6)

unimale<-combined
unimale$unimale<-as.numeric(exp(unimale$males)<1.5)
unimale$unimale<-round(unimale$unimale)
unimale<-unimale[is.na(unimale$unimale)==F,]
unimale[exp(unimale$females)<2,]$unimale<-0
unimale_phylodata<-as.data.frame(unimale %>% group_by(corrected_species_id) %>% summarise(mean(unimale,na.rm=T )) )
mostly_females_phylodata<-left_join(mostly_females_phylodata_1,unimale_phylodata,by="corrected_species_id")
rownames(mostly_females_phylodata)<-mostly_females_phylodata$corrected_species_id
colnames(mostly_females_phylodata)<-c("Species","MaleDominance","UniMale")
mostly_females_phylodata<-mostly_females_phylodata[is.na(mostly_females_phylodata$UniMale)==F,]
mostly_females_phylodata<-mostly_females_phylodata[is.na(mostly_females_phylodata$MaleDominance)==F,]
mostly_females_phylodata$UniMale<-round(mostly_females_phylodata$UniMale)

# We match the species in the dataset to the species in the phylogenetic tree
missing<-treedata(inputtree,data=mostly_females_phylodata,warnings=FALSE)
# We remove species with no data from the tree
mtree<-missing$phy
speciesnames<-mtree$tip.label

mdata<-mostly_females_phylodata[speciesnames,]

bayestraitsdiscretedata<-matrix(NA,nrow=nrow(mdata),ncol=3)
colnames(bayestraitsdiscretedata)<-c("Species","UniMale","MaleDominance")
rownames(bayestraitsdiscretedata)<-mtree$tip.label
bayestraitsdiscretedata[,1]<-mtree$tip.label
bayestraitsdiscretedata[,2]<-as.numeric(as.factor(mdata$UniMale))-1
bayestraitsdiscretedata[,3]<- as.numeric(as.factor(mdata$Male))-1
bayestraitsdiscretedata<-as.data.frame(bayestraitsdiscretedata)

command_vec1 <- c("2", "1","mltries = 100") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
results_1 <- bayestraits(bayestraitsdiscretedata, mtree, command_vec1)
log_1 <- results_1$Log$results


command_vec2 <- c("3", "1","mltries = 100") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
results_2 <- bayestraits(bayestraitsdiscretedata, mtree, command_vec2)
log_2 <- results_2$Log$results

# 1 = 0,0 not polygynous, not male dominance - 2 = 0,1 not polygynous, male dominance   - 3 = 1,0 polygynous, not male dominance   4 = 1,1 polygynous, male dominance
# Transition to male dominance 1-2 and 3-4: only when not polygynous
# Transition to polygyny 1-3 and 2-4: only when not male dominance

log_2$Count_FemaleAndCoDominance<-table(bayestraitsdiscretedata$MaleDominance)[1]
log_2$Count_MaleDominance<-table(bayestraitsdiscretedata$MaleDominance)[2]
log_2$NotMonogamous<-table(bayestraitsdiscretedata$UniMale)[1]
log_2$Monogamous<-table(bayestraitsdiscretedata$UniMale)[2]

colnames(log_2)<-c("Tree","Likelihood","TransitionToMaleDominance_AbsenceUnimale","TransitionToUnimale_AbsenceMaleDominance","LossMaleDominance_AbsenceUnimale","TransitionToUnimale_PresenceMaleDominance","LossUnimale_AbsenceMaleDominance","TransitionToMaleDominance_PresenceUnimale","LossUnimale_PresenceMaleDominance","LossMaleDominance_PresenceUnimale","LikelihoodRootNotUnimaleNotMaleDominance","LikelihoodRootNotUnimaleMaleDominance","LikelihoodRootUnimaleNotMaleDominance","LikelihoodRootUnimaleMaleDominance","CountNotMaleDominant","CountMaleDominant","CountNotUnimale","CountUnimale")

log_2$Positive_IndependentModelBetter<-2*(log_2$Likelihood-log_1$Lh)

write.csv(log_2,file="results/Bayestraits_Dependent_StrictMaleDominanceUnimale.csv")



# Testing whether there is dependent evolution of dimorphism and female dominance

mostly_females_phylodata_1<-combined
mostly_females_phylodata_1$female_dominance<-as.numeric(as.factor(combined$strictfdom))
mostly_females_phylodata_1<-as.data.frame(mostly_females_phylodata_1 %>% group_by(corrected_species_id) %>% summarise(mean(female_dominance,na.rm=T )) )
colnames(mostly_females_phylodata_1)<-c("corrected_species_id", "femaledominance")
mostly_females_phylodata_1$femaledominance<-as.numeric(mostly_females_phylodata_1$femaledominance>2.4)


dimorphism<-combined
dimorphism<-dimorphism[is.na(dimorphism$SexualDimorphism_MaleWeight_over_FemaleWeight)==F,]
dimorphism[exp(dimorphism$SexualDimorphism_MaleWeight_over_FemaleWeight)<1.201,]$SexualDimorphism_MaleWeight_over_FemaleWeight<- 0
dimorphism[exp(dimorphism$SexualDimorphism_MaleWeight_over_FemaleWeight)>1.2,]$SexualDimorphism_MaleWeight_over_FemaleWeight<- 1
dimorphism_phylodata<-as.data.frame(dimorphism %>% group_by(corrected_species_id) %>% summarise(mean(SexualDimorphism_MaleWeight_over_FemaleWeight,na.rm=T )) )

mostly_females_phylodata<-left_join(mostly_females_phylodata_1,dimorphism_phylodata,by="corrected_species_id")
rownames(mostly_females_phylodata)<-mostly_females_phylodata$corrected_species_id
colnames(mostly_females_phylodata)<-c("Species","FemaleDominance","Dimorphic")
mostly_females_phylodata<-mostly_females_phylodata[is.na(mostly_females_phylodata$Dimorphic)==F,]
mostly_females_phylodata<-mostly_females_phylodata[is.na(mostly_females_phylodata$FemaleDominance)==F,]


# We match the species in the dataset to the species in the phylogenetic tree
missing<-treedata(inputtree,data=mostly_females_phylodata,warnings=FALSE)
# We remove species with no data from the tree
mtree<-missing$phy
speciesnames<-mtree$tip.label

mdata<-mostly_females_phylodata[speciesnames,]
mdata$FemaleDominance<-round(mdata$FemaleDominance,0)

bayestraitsdiscretedata<-matrix(NA,nrow=nrow(mdata),ncol=3)
colnames(bayestraitsdiscretedata)<-c("Species","Dimorphic","FemaleDominance")
rownames(bayestraitsdiscretedata)<-mdata$Species
bayestraitsdiscretedata[,1]<-mdata$Species
bayestraitsdiscretedata[,2]<-as.numeric(as.factor(mdata$Dimorphic))-1
bayestraitsdiscretedata[,3]<- as.numeric(as.factor(mdata$FemaleDominance))-1
bayestraitsdiscretedata<-as.data.frame(bayestraitsdiscretedata)

command_vec1 <- c("2", "1") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
results_1 <- bayestraits(bayestraitsdiscretedata, mtree, command_vec1)
log_1 <- results_1$Log$results


command_vec2 <- c("3", "1","mltries 100") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
results_2 <- bayestraits(bayestraitsdiscretedata, mtree, command_vec2)
log_2 <- results_2$Log$results


# 1 = 0,0 not dimorphic, not female dominance - 2 = 0,1 not dimorphic, female dominance   - 3 = 1,0 dimorphic, not female dominance   4 = 1,1 dimorphic, female dominance


log_2$Count_FemaleAndCoDominance<-table(bayestraitsdiscretedata$FemaleDominance)[1]
log_2$Count_FemaleDominance<-table(bayestraitsdiscretedata$FemaleDominance)[2]
log_2$NotMonogamous<-table(bayestraitsdiscretedata$Dimorphic)[1]
log_2$Monogamous<-table(bayestraitsdiscretedata$Dimorphic)[2]

colnames(log_2)<-c("Tree","Likelihood","TransitionToFemaleDominance_Absencedimorphism","TransitionTodimorphism_AbsenceFemaleDominance","LossFemaleDominance_Absencedimorphism","TransitionTodimorphism_PresenceFemaleDominance","Lossdimorphism_AbsenceFemaleDominance","TransitionToFemaleDominance_Presencedimorphism","Lossdimorphism_PresenceFemaleDominance","LossFemaleDominance_Presencedimorphism","LikelihoodRootNotdimorphismNotFemaleDominance","LikelihoodRootNotdimorphismFemaleDominance","LikelihoodRootdimorphismNotFemaleDominance","LikelihoodRootdimorphismFemaleDominance","CountNotFemaleDominant","CountFemaleDominant","CountNotdimorphic","CountDimorphic")

log_2$Positive_IndependentModelBetter<-2*(log_2$Likelihood-log_1$Lh)

write.csv(log_2,file="results/Bayestraits_Dependent_StrictFemaleDominanceDimorphism.csv")




# Testing whether there is dependent evolution of dimorphism and male dominance

mostly_females_phylodata_1<-combined
mostly_females_phylodata_1$female_dominance<-as.numeric(as.factor(combined$strictfdom))
mostly_females_phylodata_1<-as.data.frame(mostly_females_phylodata_1 %>% group_by(corrected_species_id) %>% summarise(mean(female_dominance,na.rm=T )) )
colnames(mostly_females_phylodata_1)<-c("corrected_species_id", "femaledominance")
mostly_females_phylodata_1$femaledominance<-as.numeric(mostly_females_phylodata_1$femaledominance<1.6)


dimorphism<-combined
dimorphism<-dimorphism[is.na(dimorphism$SexualDimorphism_MaleWeight_over_FemaleWeight)==F,]
dimorphism[exp(dimorphism$SexualDimorphism_MaleWeight_over_FemaleWeight)<1.201,]$SexualDimorphism_MaleWeight_over_FemaleWeight<- 0
dimorphism[exp(dimorphism$SexualDimorphism_MaleWeight_over_FemaleWeight)>1.2,]$SexualDimorphism_MaleWeight_over_FemaleWeight<- 1
dimorphism_phylodata<-as.data.frame(dimorphism %>% group_by(corrected_species_id) %>% summarise(mean(SexualDimorphism_MaleWeight_over_FemaleWeight,na.rm=T )) )

mostly_females_phylodata<-left_join(mostly_females_phylodata_1,dimorphism_phylodata,by="corrected_species_id")
rownames(mostly_females_phylodata)<-mostly_females_phylodata$corrected_species_id
colnames(mostly_females_phylodata)<-c("Species","MaleDominance","Dimorphic")
mostly_females_phylodata<-mostly_females_phylodata[is.na(mostly_females_phylodata$Dimorphic)==F,]
mostly_females_phylodata<-mostly_females_phylodata[is.na(mostly_females_phylodata$MaleDominance)==F,]


# We match the species in the dataset to the species in the phylogenetic tree
missing<-treedata(inputtree,data=mostly_females_phylodata,warnings=FALSE)
# We remove species with no data from the tree
mtree<-missing$phy
speciesnames<-mtree$tip.label

mdata<-mostly_females_phylodata[speciesnames,]
mdata$MaleDominance<-round(mdata$MaleDominance,0)

bayestraitsdiscretedata<-matrix(NA,nrow=nrow(mdata),ncol=3)
colnames(bayestraitsdiscretedata)<-c("Species","Dimorphic","MaleDominance")
rownames(bayestraitsdiscretedata)<-mdata$Species
bayestraitsdiscretedata[,1]<-mdata$Species
bayestraitsdiscretedata[,2]<-as.numeric(as.factor(mdata$Dimorphic))-1
bayestraitsdiscretedata[,3]<- as.numeric(as.factor(mdata$MaleDominance))-1
bayestraitsdiscretedata<-as.data.frame(bayestraitsdiscretedata)

command_vec1 <- c("2", "1", "mltries 100") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
results_1 <- bayestraits(bayestraitsdiscretedata, mtree, command_vec1)
log_1 <- results_1$Log$results


command_vec2 <- c("3", "1","mltries 100") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
results_2 <- bayestraits(bayestraitsdiscretedata, mtree, command_vec2)
log_2 <- results_2$Log$results


# 1 = 0,0 not dimorphic, not male dominance - 2 = 0,1 not dimorphic, male dominance   - 3 = 1,0 dimorphic, not male dominance   4 = 1,1 dimorphic, male dominance


log_2$Count_FemaleAndCoDominance<-table(bayestraitsdiscretedata$MaleDominance)[1]
log_2$Count_MaleDominance<-table(bayestraitsdiscretedata$MaleDominance)[2]
log_2$NotMonogamous<-table(bayestraitsdiscretedata$Dimorphic)[1]
log_2$Monogamous<-table(bayestraitsdiscretedata$Dimorphic)[2]

colnames(log_2)<-c("Tree","Likelihood","TransitionToMaleDominance_Absencedimorphism","TransitionTodimorphism_AbsenceMaleDominance","LossMaleDominance_Absencedimorphism","TransitionTodimorphism_PresenceMaleDominance","Lossdimorphism_AbsenceMaleDominance","TransitionToMaleDominance_Presencedimorphism","Lossdimorphism_PresenceMaleDominance","LossMaleDominance_Presencedimorphism","LikelihoodRootNotdimorphismNotMaleDominance","LikelihoodRootNotdimorphismMaleDominance","LikelihoodRootdimorphismNotMaleDominance","LikelihoodRootdimorphismMaleDominance","CountNotMaleDominant","CountMaleDominant","CountNotdimorphic","CountDimorphic")

log_2$Positive_IndependentModelBetter<-2*(log_2$Likelihood-log_1$Lh)

write.csv(log_2,file="results/Bayestraits_Dependent_StrictMaleDominanceDimorphism.csv")



# Testing whether there is dependent evolution of low dimorphism and male dominance

mostly_females_phylodata_1<-combined
mostly_females_phylodata_1$female_dominance<-as.numeric(as.factor(combined$strictfdom))
mostly_females_phylodata_1<-as.data.frame(mostly_females_phylodata_1 %>% group_by(corrected_species_id) %>% summarise(mean(female_dominance,na.rm=T )) )
colnames(mostly_females_phylodata_1)<-c("corrected_species_id", "femaledominance")
mostly_females_phylodata_1$femaledominance<-as.numeric(mostly_females_phylodata_1$femaledominance<1.6)


dimorphism<-combined
dimorphism<-dimorphism[is.na(dimorphism$SexualDimorphism_MaleWeight_over_FemaleWeight)==F,]
dimorphism[exp(dimorphism$SexualDimorphism_MaleWeight_over_FemaleWeight)<1.101,]$SexualDimorphism_MaleWeight_over_FemaleWeight<- 0
dimorphism[exp(dimorphism$SexualDimorphism_MaleWeight_over_FemaleWeight)>1.1,]$SexualDimorphism_MaleWeight_over_FemaleWeight<- 1
dimorphism_phylodata<-as.data.frame(dimorphism %>% group_by(corrected_species_id) %>% summarise(mean(SexualDimorphism_MaleWeight_over_FemaleWeight,na.rm=T )) )

mostly_females_phylodata<-left_join(mostly_females_phylodata_1,dimorphism_phylodata,by="corrected_species_id")
rownames(mostly_females_phylodata)<-mostly_females_phylodata$corrected_species_id
colnames(mostly_females_phylodata)<-c("Species","MaleDominance","Dimorphic")
mostly_females_phylodata<-mostly_females_phylodata[is.na(mostly_females_phylodata$Dimorphic)==F,]
mostly_females_phylodata<-mostly_females_phylodata[is.na(mostly_females_phylodata$MaleDominance)==F,]


# We match the species in the dataset to the species in the phylogenetic tree
missing<-treedata(inputtree,data=mostly_females_phylodata,warnings=FALSE)
# We remove species with no data from the tree
mtree<-missing$phy
speciesnames<-mtree$tip.label

mdata<-mostly_females_phylodata[speciesnames,]
mdata$MaleDominance<-round(mdata$MaleDominance,0)

bayestraitsdiscretedata<-matrix(NA,nrow=nrow(mdata),ncol=3)
colnames(bayestraitsdiscretedata)<-c("Species","Dimorphic","MaleDominance")
rownames(bayestraitsdiscretedata)<-mdata$Species
bayestraitsdiscretedata[,1]<-mdata$Species
bayestraitsdiscretedata[,2]<-as.numeric(as.factor(mdata$Dimorphic))-1
bayestraitsdiscretedata[,3]<- as.numeric(as.factor(mdata$MaleDominance))-1
bayestraitsdiscretedata<-as.data.frame(bayestraitsdiscretedata)

command_vec1 <- c("2", "1", "mltries 100") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
results_1 <- bayestraits(bayestraitsdiscretedata, mtree, command_vec1)
log_1 <- results_1$Log$results


command_vec2 <- c("3", "1","mltries 100") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
results_2 <- bayestraits(bayestraitsdiscretedata, mtree, command_vec2)
log_2 <- results_2$Log$results


# 1 = 0,0 not dimorphic, not male dominance - 2 = 0,1 not dimorphic, male dominance   - 3 = 1,0 dimorphic, not male dominance   4 = 1,1 dimorphic, male dominance


log_2$Count_FemaleAndCoDominance<-table(bayestraitsdiscretedata$MaleDominance)[1]
log_2$Count_MaleDominance<-table(bayestraitsdiscretedata$MaleDominance)[2]
log_2$NotMonogamous<-table(bayestraitsdiscretedata$Dimorphic)[1]
log_2$Monogamous<-table(bayestraitsdiscretedata$Dimorphic)[2]

colnames(log_2)<-c("Tree","Likelihood","TransitionToMaleDominance_Absencedimorphism","TransitionTodimorphism_AbsenceMaleDominance","LossMaleDominance_Absencedimorphism","TransitionTodimorphism_PresenceMaleDominance","Lossdimorphism_AbsenceMaleDominance","TransitionToMaleDominance_Presencedimorphism","Lossdimorphism_PresenceMaleDominance","LossMaleDominance_Presencedimorphism","LikelihoodRootNotdimorphismNotMaleDominance","LikelihoodRootNotdimorphismMaleDominance","LikelihoodRootdimorphismNotMaleDominance","LikelihoodRootdimorphismMaleDominance","CountNotMaleDominant","CountMaleDominant","CountNotdimorphic","CountDimorphic")

log_2$Positive_IndependentModelBetter<-2*(log_2$Likelihood-log_1$Lh)

write.csv(log_2,file="results/Bayestraits_Dependent_StrictMaleDominanceLowDimorphism.csv")






##########################################################################################
##########################################################################################
##########################################################################################


# Analyses using the relaxed two way classification of dominance

############################################################

# Multistate model investigating the transitions between male- and female- dominance

mostly_females_phylodata_1<-combined
mostly_females_phylodata_1$mostlyfdom<-as.numeric(as.factor(combined$mostlyfdom))
mostly_females_phylodata_1<-as.data.frame(mostly_females_phylodata_1 %>% group_by(corrected_species_id) %>% summarise(mean(mostlyfdom,na.rm=T )) )
colnames(mostly_females_phylodata_1)<-c("Species","mostlyfdom")
rownames(mostly_females_phylodata_1)<-mostly_females_phylodata_1$Species
mostly_females_phylodata_1<-mostly_females_phylodata_1[is.na(mostly_females_phylodata_1$mostlyfdom)==F,]

missing<-treedata(inputtree,data=mostly_females_phylodata_1,warnings=FALSE)
# We remove species with no data from the tree
mtree<-missing$phy
speciesnames<-mtree$tip.label

mdata<-mostly_females_phylodata_1[mostly_females_phylodata_1$Species %in% speciesnames,]
mdata$mostlyfdom<-round(mdata$mostlyfdom)

bayestraitsdiscretedata<-matrix(NA,nrow=nrow(mdata),ncol=2)
colnames(bayestraitsdiscretedata)<-c("Species","mostlyfdom")
rownames(bayestraitsdiscretedata)<-mtree$tip.label
bayestraitsdiscretedata[,1]<-mtree$tip.label
bayestraitsdiscretedata[,2]<- as.numeric(as.factor(mdata$mostlyfdom))
bayestraitsdiscretedata<-as.data.frame(bayestraitsdiscretedata)

command_vec1 <- c("1", "1","MLTries = 100") #option 1 = 1 discrete independent; option 2 = 1 maximum likelihood
results_1 <- bayestraits(bayestraitsdiscretedata, mtree, command_vec1)
log_1 <- results_1$Log$results

colnames(log_1)<-c("Tree","Likelihood","Transition Male to Female","Transition Female to Male","Likelihood Root is Male","Likelihood Root is Female")

log_1$Count_MaleDominance<-table(bayestraitsdiscretedata$mostlyfdom)[1]
log_1$Count_FemaleDominance<-table(bayestraitsdiscretedata$mostlyfdom)[2]

write.csv(log_1,file="results/Bayestraits_Transitions_mostlyIntersexDominance.csv")





############################################################


# Testing whether there is dependent evolution of monogamy and female dominance

mostly_females_phylodata_1<-combined
mostly_females_phylodata_1$female_dominance<-as.numeric(as.factor(combined$mostlyfdom))
mostly_females_phylodata_1<-as.data.frame(mostly_females_phylodata_1 %>% group_by(corrected_species_id) %>% summarise(mean(female_dominance,na.rm=T )) )
colnames(mostly_females_phylodata_1)<-c("corrected_species_id", "femaledominance")
mostly_females_phylodata_1$femaledominance<-as.numeric(mostly_females_phylodata_1$femaledominance>1.49)

monogamy<-combined
monogamy$MatSysPMK<-as.numeric(as.factor(monogamy$MatSysPMK))
monogamy<-monogamy[is.na(monogamy$MatSysPMK)==F,]
monogamy[monogamy$MatSysPMK>1,]$MatSysPMK<-0
monogamy_phylodata<-as.data.frame(monogamy %>% group_by(corrected_species_id) %>% summarise(mean(MatSysPMK,na.rm=T )) )
mostly_females_phylodata<-left_join(mostly_females_phylodata_1,monogamy_phylodata,by="corrected_species_id")
rownames(mostly_females_phylodata)<-mostly_females_phylodata$corrected_species_id
colnames(mostly_females_phylodata)<-c("Species","FemaleDominance","Monogamy")
mostly_females_phylodata<-mostly_females_phylodata[is.na(mostly_females_phylodata$Monogamy)==F,]
mostly_females_phylodata<-mostly_females_phylodata[is.na(mostly_females_phylodata$FemaleDominance)==F,]


# We match the species in the dataset to the species in the phylogenetic tree
missing<-treedata(inputtree,data=mostly_females_phylodata,warnings=FALSE)
# We remove species with no data from the tree
mtree<-missing$phy
speciesnames<-mtree$tip.label

mdata<-mostly_females_phylodata[speciesnames,]

bayestraitsdiscretedata<-matrix(NA,nrow=nrow(mdata),ncol=3)
colnames(bayestraitsdiscretedata)<-c("Species","Monogamy","FemaleDominance")
rownames(bayestraitsdiscretedata)<-mtree$tip.label
bayestraitsdiscretedata[,1]<-mtree$tip.label
bayestraitsdiscretedata[,2]<-as.numeric(as.factor(mdata$Monogamy))-1
bayestraitsdiscretedata[,3]<- as.numeric(as.factor(mdata$FemaleDominance))-1
bayestraitsdiscretedata<-as.data.frame(bayestraitsdiscretedata)

command_vec1 <- c("2", "1","mltries = 100") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
results_1 <- bayestraits(bayestraitsdiscretedata, mtree, command_vec1)
log_1 <- results_1$Log$results


command_vec2 <- c("3", "1","mltries = 100") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
results_2 <- bayestraits(bayestraitsdiscretedata, mtree, command_vec2)
log_2 <- results_2$Log$results

# 1 = 0,0 not monogamous, not female dominance - 2 = 0,1 not monogamous, female dominance   - 3 = 1,0 monogamous, not female dominance   4 = 1,1 monogamous, female dominance
# Transition to female dominance 1-2 and 3-4: more likely when monogamy present (15x more)
# Transition to monogamy 1-3 and 2-4: more likely when females not dominant (5x more)
# might suggest that monogamy comes first

colnames(log_2)<-c("Tree","Likelihood","TransitionToFemaleDominance_AbsenceMonogamy","TransitionToMonogamy_AbsenceFemaleDominance","LossFemaleDominance_AbsenceMonogamy","TransitionToMonogamy_PresenceFemaleDominance","LossMonogamy_AbsenceFemaleDominance","TransitionToFemaleDominance_PresenceMonogamy","LossMonogamy_PresenceFemaleDominance","LossFemaleDominance_PresenceMonogamy","LikelihoodRootNotMonogamyNotFemaleDominance","LikelihoodRootNotMonogamyFemaleDominance","LikelihoodRootMonogamyNotFemaleDominance","LikelihoodRootMonogamyFemaleDominance")


log_2$SpeciesCount_NotMonogamous_NotFemaleDominant<-as.data.frame(bayestraitsdiscretedata %>% group_by(Monogamy,FemaleDominance) %>% summarise(n()))[1,3]
log_2$SpeciesCount_NotMonogamous_FemaleDominant<-as.data.frame(bayestraitsdiscretedata %>% group_by(Monogamy,FemaleDominance) %>% summarise(n()))[2,3]
log_2$SpeciesCount_Monogamous_NotFemaleDominant<-as.data.frame(bayestraitsdiscretedata %>% group_by(Monogamy,FemaleDominance) %>% summarise(n()))[3,3]
log_2$SpeciesCount_Monogamous_FemaleDominant<-as.data.frame(bayestraitsdiscretedata %>% group_by(Monogamy,FemaleDominance) %>% summarise(n()))[4,3]


log_2$Positive_IndependentModelBetter<-2*(log_2$Likelihood-log_1$Lh)

write.csv(log_2,file="results/Bayestraits_Dependent_RelaxedFemaleDominanceMonogamy.csv")




# Testing whether there is dependent evolution of polygyny and female dominance

mostly_females_phylodata_1<-combined
mostly_females_phylodata_1$female_dominance<-as.numeric(as.factor(combined$mostlyfdom))
mostly_females_phylodata_1<-as.data.frame(mostly_females_phylodata_1 %>% group_by(corrected_species_id) %>% summarise(mean(female_dominance,na.rm=T )) )
colnames(mostly_females_phylodata_1)<-c("corrected_species_id", "femaledominance")
mostly_females_phylodata_1$femaledominance<-as.numeric(mostly_females_phylodata_1$femaledominance>1.49)

polygyny<-combined
polygyny$MatSysPMK<-as.numeric(combined$MatSysPMK=="POL")
polygyny<-polygyny[is.na(polygyny$MatSysPMK)==F,]
polygyny$MatSysPMK<-polygyny$MatSysPMK
polygyny_phylodata<-as.data.frame(polygyny %>% group_by(corrected_species_id) %>% summarise(mean(MatSysPMK,na.rm=T )) )
mostly_females_phylodata<-left_join(mostly_females_phylodata_1,polygyny_phylodata,by="corrected_species_id")
rownames(mostly_females_phylodata)<-mostly_females_phylodata$corrected_species_id
colnames(mostly_females_phylodata)<-c("Species","FemaleDominance","Polygyny")
mostly_females_phylodata<-mostly_females_phylodata[is.na(mostly_females_phylodata$Polygyny)==F,]
mostly_females_phylodata<-mostly_females_phylodata[is.na(mostly_females_phylodata$FemaleDominance)==F,]


# We match the species in the dataset to the species in the phylogenetic tree
missing<-treedata(inputtree,data=mostly_females_phylodata,warnings=FALSE)
# We remove species with no data from the tree
mtree<-missing$phy
speciesnames<-mtree$tip.label

mdata<-mostly_females_phylodata[speciesnames,]

bayestraitsdiscretedata<-matrix(NA,nrow=nrow(mdata),ncol=3)
colnames(bayestraitsdiscretedata)<-c("Species","Polygyny","FemaleDominance")
rownames(bayestraitsdiscretedata)<-mtree$tip.label
bayestraitsdiscretedata[,1]<-mtree$tip.label
bayestraitsdiscretedata[,2]<-as.numeric(as.factor(mdata$Polygyny))-1
bayestraitsdiscretedata[,3]<- as.numeric(as.factor(mdata$FemaleDominance))-1
bayestraitsdiscretedata<-as.data.frame(bayestraitsdiscretedata)

command_vec1 <- c("2", "1","mltries = 100") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
results_1 <- bayestraits(bayestraitsdiscretedata, mtree, command_vec1)
log_1 <- results_1$Log$results


command_vec2 <- c("3", "1","mltries = 100") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
results_2 <- bayestraits(bayestraitsdiscretedata, mtree, command_vec2)
log_2 <- results_2$Log$results

# 1 = 0,0 not polygynous, not female dominance - 2 = 0,1 not polygynous, female dominance   - 3 = 1,0 polygynous, not female dominance   4 = 1,1 polygynous, female dominance
# Transition to female dominance 1-2 and 3-4: only when not polygynous
# Transition to polygyny 1-3 and 2-4: only when not female dominance

colnames(log_2)<-c("Tree","Likelihood","TransitionToFemaleDominance_AbsencePolygyny","TransitionToPolygyny_AbsenceFemaleDominance","LossFemaleDominance_AbsencePolygyny","TransitionToPolygyny_PresenceFemaleDominance","LossPolygyny_AbsenceFemaleDominance","TransitionToFemaleDominance_PresencePolygyny","LossPolygyny_PresenceFemaleDominance","LossFemaleDominance_PresencePolygyny","LikelihoodRootNotPolygynyNotFemaleDominance","LikelihoodRootNotPolygynyFemaleDominance","LikelihoodRootPolygynyNotFemaleDominance","LikelihoodRootPolygynyFemaleDominance")


log_2$SpeciesCount_NotPolygynous_NotFemaleDominant<-as.data.frame(bayestraitsdiscretedata %>% group_by(Polygyny,FemaleDominance) %>% summarise(n()))[1,3]
log_2$SpeciesCount_NotPolygynous_FemaleDominant<-as.data.frame(bayestraitsdiscretedata %>% group_by(Polygyny,FemaleDominance) %>% summarise(n()))[2,3]
log_2$SpeciesCount_Polygynous_NotFemaleDominant<-as.data.frame(bayestraitsdiscretedata %>% group_by(Polygyny,FemaleDominance) %>% summarise(n()))[3,3]
log_2$SpeciesCount_Polygynous_FemaleDominant<-as.data.frame(bayestraitsdiscretedata %>% group_by(Polygyny,FemaleDominance) %>% summarise(n()))[4,3]


log_2$Positive_IndependentModelBetter<-2*(log_2$Likelihood-log_1$Lh)

write.csv(log_2,file="results/Bayestraits_Dependent_RelaxedFemaleDominancePolygyny.csv")




# Testing whether there is dependent evolution of single-male groups and female dominance

mostly_females_phylodata_1<-combined
mostly_females_phylodata_1$female_dominance<-as.numeric(as.factor(combined$mostlyfdom))
mostly_females_phylodata_1<-as.data.frame(mostly_females_phylodata_1 %>% group_by(corrected_species_id) %>% summarise(mean(female_dominance,na.rm=T )) )
colnames(mostly_females_phylodata_1)<-c("corrected_species_id", "femaledominance")
mostly_females_phylodata_1$femaledominance<-as.numeric(mostly_females_phylodata_1$femaledominance>1.49)

unimale<-combined
unimale$unimale<-as.numeric(exp(unimale$males)<1.5)
unimale$unimale<-round(unimale$unimale)
unimale<-unimale[is.na(unimale$unimale)==F,]
unimale[exp(unimale$females)<2,]$unimale<-0
unimale_phylodata<-as.data.frame(unimale %>% group_by(corrected_species_id) %>% summarise(mean(unimale,na.rm=T )) )
mostly_females_phylodata<-left_join(mostly_females_phylodata_1,unimale_phylodata,by="corrected_species_id")
rownames(mostly_females_phylodata)<-mostly_females_phylodata$corrected_species_id
colnames(mostly_females_phylodata)<-c("Species","FemaleDominance","UniMale")
mostly_females_phylodata<-mostly_females_phylodata[is.na(mostly_females_phylodata$UniMale)==F,]
mostly_females_phylodata<-mostly_females_phylodata[is.na(mostly_females_phylodata$FemaleDominance)==F,]
mostly_females_phylodata$UniMale<-round(mostly_females_phylodata$UniMale)

# We match the species in the dataset to the species in the phylogenetic tree
missing<-treedata(inputtree,data=mostly_females_phylodata,warnings=FALSE)
# We remove species with no data from the tree
mtree<-missing$phy
speciesnames<-mtree$tip.label

mdata<-mostly_females_phylodata[speciesnames,]

bayestraitsdiscretedata<-matrix(NA,nrow=nrow(mdata),ncol=3)
colnames(bayestraitsdiscretedata)<-c("Species","UniMale","FemaleDominance")
rownames(bayestraitsdiscretedata)<-mtree$tip.label
bayestraitsdiscretedata[,1]<-mtree$tip.label
bayestraitsdiscretedata[,2]<-as.numeric(as.factor(mdata$UniMale))-1
bayestraitsdiscretedata[,3]<- as.numeric(as.factor(mdata$FemaleDominance))-1
bayestraitsdiscretedata<-as.data.frame(bayestraitsdiscretedata)

command_vec1 <- c("2", "1","mltries = 100") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
results_1 <- bayestraits(bayestraitsdiscretedata, mtree, command_vec1)
log_1 <- results_1$Log$results


command_vec2 <- c("3", "1","mltries = 100") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
results_2 <- bayestraits(bayestraitsdiscretedata, mtree, command_vec2)
log_2 <- results_2$Log$results

# 1 = 0,0 not unimale, not female dominance - 2 = 0,1 not unimale, female dominance   - 3 = 1,0 unimale, not female dominance   4 = 1,1 unimale, female dominance
# Transition to female dominance 1-2 and 3-4: only when not unimale
# Transition to unimale 1-3 and 2-4: only when not male dominance

colnames(log_2)<-c("Tree","Likelihood","TransitionToFemaleDominance_AbsenceUniMale","TransitionToUnimale_AbsenceFemaleDominance","LossFemaleDominance_AbsenceUnimale","TransitionToUnimale_PresenceFemaleDominance","LossUnimale_AbsenceFemaleDominance","TransitionToFemaleDominance_PresenceUnimale","LossUnimale_PresenceFemaleDominance","LossFemaleDominance_PresenceUnimale","LikelihoodRootNotUnimaleNotFemaleDominance","LikelihoodRootNotUnimaleFemaleDominance","LikelihoodRootUnimaleNotFemaleDominance","LikelihoodRootUnimaleFemaleDominance")


log_2$SpeciesCount_NotPolygynous_NotFemaleDominant<-as.data.frame(bayestraitsdiscretedata %>% group_by(UniMale,FemaleDominance) %>% summarise(n()))[1,3]
log_2$SpeciesCount_NotPolygynous_FemaleDominant<-as.data.frame(bayestraitsdiscretedata %>% group_by(UniMale,FemaleDominance) %>% summarise(n()))[2,3]
log_2$SpeciesCount_Polygynous_NotFemaleDominant<-as.data.frame(bayestraitsdiscretedata %>% group_by(UniMale,FemaleDominance) %>% summarise(n()))[3,3]
log_2$SpeciesCount_Polygynous_FemaleDominant<-as.data.frame(bayestraitsdiscretedata %>% group_by(UniMale,FemaleDominance) %>% summarise(n()))[4,3]


log_2$Positive_IndependentModelBetter<-2*(log_2$Likelihood-log_1$Lh)

write.csv(log_2,file="results/Bayestraits_Dependent_RelaxedFemaleDominanceUnimale.csv")



# Testing whether there is dependent evolution of dimorphism and female dominance

mostly_females_phylodata_1<-combined
mostly_females_phylodata_1$female_dominance<-as.numeric(as.factor(combined$mostlyfdom))
mostly_females_phylodata_1<-as.data.frame(mostly_females_phylodata_1 %>% group_by(corrected_species_id) %>% summarise(mean(female_dominance,na.rm=T )) )
colnames(mostly_females_phylodata_1)<-c("corrected_species_id", "femaledominance")
mostly_females_phylodata_1$femaledominance<-as.numeric(mostly_females_phylodata_1$femaledominance>1.49)


dimorphism<-combined
dimorphism<-dimorphism[is.na(dimorphism$SexualDimorphism_MaleWeight_over_FemaleWeight)==F,]
dimorphism[exp(dimorphism$SexualDimorphism_MaleWeight_over_FemaleWeight)<1.201,]$SexualDimorphism_MaleWeight_over_FemaleWeight<- 0
dimorphism[exp(dimorphism$SexualDimorphism_MaleWeight_over_FemaleWeight)>1.2,]$SexualDimorphism_MaleWeight_over_FemaleWeight<- 1
dimorphism_phylodata<-as.data.frame(dimorphism %>% group_by(corrected_species_id) %>% summarise(mean(SexualDimorphism_MaleWeight_over_FemaleWeight,na.rm=T )) )

mostly_females_phylodata<-left_join(mostly_females_phylodata_1,dimorphism_phylodata,by="corrected_species_id")
rownames(mostly_females_phylodata)<-mostly_females_phylodata$corrected_species_id
colnames(mostly_females_phylodata)<-c("Species","FemaleDominance","Dimorphic")
mostly_females_phylodata<-mostly_females_phylodata[is.na(mostly_females_phylodata$Dimorphic)==F,]
mostly_females_phylodata<-mostly_females_phylodata[is.na(mostly_females_phylodata$FemaleDominance)==F,]


# We match the species in the dataset to the species in the phylogenetic tree
missing<-treedata(inputtree,data=mostly_females_phylodata,warnings=FALSE)
# We remove species with no data from the tree
mtree<-missing$phy
speciesnames<-mtree$tip.label

mdata<-mostly_females_phylodata[speciesnames,]


bayestraitsdiscretedata<-matrix(NA,nrow=nrow(mdata),ncol=3)
colnames(bayestraitsdiscretedata)<-c("Species","Dimorphic","FemaleDominance")
rownames(bayestraitsdiscretedata)<-mdata$Species
bayestraitsdiscretedata[,1]<-mdata$Species
bayestraitsdiscretedata[,2]<-as.numeric(as.factor(mdata$Dimorphic))-1
bayestraitsdiscretedata[,3]<- as.numeric(as.factor(mdata$FemaleDominance))-1
bayestraitsdiscretedata<-as.data.frame(bayestraitsdiscretedata)

command_vec1 <- c("2", "1") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
results_1 <- bayestraits(bayestraitsdiscretedata, mtree, command_vec1)
log_1 <- results_1$Log$results


command_vec2 <- c("3", "1","mltries 100") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
results_2 <- bayestraits(bayestraitsdiscretedata, mtree, command_vec2)
log_2 <- results_2$Log$results


# 1 = 0,0 not dimorphic, not female dominance - 2 = 0,1 not dimorphic, female dominance   - 3 = 1,0 dimorphic, not female dominance   4 = 1,1 dimorphic, female dominance
colnames(log_2)<-c("Tree","Likelihood","TransitionToFemaleDominance_AbsenceDimorphism","TransitionToDimorphism_AbsenceFemaleDominance","LossFemaleDominance_AbsenceDimorphism","TransitionToDimorphism_PresenceFemaleDominance","LossDimorphism_AbsenceFemaleDominance","TransitionToFemaleDominance_PresenceDimorphism","LossDimorphism_PresenceFemaleDominance","LossFemaleDominance_PresenceDimorphism","LikelihoodRootNotDimorphismNotFemaleDominance","LikelihoodRootNotDimorphismFemaleDominance","LikelihoodRootDimorphismNotFemaleDominance","LikelihoodRootDimorphismFemaleDominance")


log_2$SpeciesCount_NotDimorphic_NotFemaleDominant<-as.data.frame(bayestraitsdiscretedata %>% group_by(Dimorphic,FemaleDominance) %>% summarise(n()))[1,3]
log_2$SpeciesCount_NotDimorphic_FemaleDominant<-as.data.frame(bayestraitsdiscretedata %>% group_by(Dimorphic,FemaleDominance) %>% summarise(n()))[2,3]
log_2$SpeciesCount_Dimorphic_NotFemaleDominant<-as.data.frame(bayestraitsdiscretedata %>% group_by(Dimorphic,FemaleDominance) %>% summarise(n()))[3,3]
log_2$SpeciesCount_Dimorphic_FemaleDominant<-as.data.frame(bayestraitsdiscretedata %>% group_by(Dimorphic,FemaleDominance) %>% summarise(n()))[4,3]


log_2$Positive_IndependentModelBetter<-2*(log_2$Likelihood-log_1$Lh)

write.csv(log_2,file="results/Bayestraits_Dependent_RelaxedFemaleDominanceDimorphism.csv")



# Testing whether there is dependent evolution of lowdimorphism and female dominance

mostly_females_phylodata_1<-combined
mostly_females_phylodata_1$female_dominance<-as.numeric(as.factor(combined$mostlyfdom))
mostly_females_phylodata_1<-as.data.frame(mostly_females_phylodata_1 %>% group_by(corrected_species_id) %>% summarise(mean(female_dominance,na.rm=T )) )
colnames(mostly_females_phylodata_1)<-c("corrected_species_id", "femaledominance")
mostly_females_phylodata_1$femaledominance<-as.numeric(mostly_females_phylodata_1$femaledominance>1.49)


dimorphism<-combined
dimorphism<-dimorphism[is.na(dimorphism$SexualDimorphism_MaleWeight_over_FemaleWeight)==F,]
dimorphism[exp(dimorphism$SexualDimorphism_MaleWeight_over_FemaleWeight)<1.101,]$SexualDimorphism_MaleWeight_over_FemaleWeight<- 0
dimorphism[exp(dimorphism$SexualDimorphism_MaleWeight_over_FemaleWeight)>1.1,]$SexualDimorphism_MaleWeight_over_FemaleWeight<- 1
dimorphism_phylodata<-as.data.frame(dimorphism %>% group_by(corrected_species_id) %>% summarise(mean(SexualDimorphism_MaleWeight_over_FemaleWeight,na.rm=T )) )

mostly_females_phylodata<-left_join(mostly_females_phylodata_1,dimorphism_phylodata,by="corrected_species_id")
rownames(mostly_females_phylodata)<-mostly_females_phylodata$corrected_species_id
colnames(mostly_females_phylodata)<-c("Species","FemaleDominance","Dimorphic")
mostly_females_phylodata<-mostly_females_phylodata[is.na(mostly_females_phylodata$Dimorphic)==F,]
mostly_females_phylodata<-mostly_females_phylodata[is.na(mostly_females_phylodata$FemaleDominance)==F,]


# We match the species in the dataset to the species in the phylogenetic tree
missing<-treedata(inputtree,data=mostly_females_phylodata,warnings=FALSE)
# We remove species with no data from the tree
mtree<-missing$phy
speciesnames<-mtree$tip.label

mdata<-mostly_females_phylodata[speciesnames,]


bayestraitsdiscretedata<-matrix(NA,nrow=nrow(mdata),ncol=3)
colnames(bayestraitsdiscretedata)<-c("Species","Dimorphic","FemaleDominance")
rownames(bayestraitsdiscretedata)<-mdata$Species
bayestraitsdiscretedata[,1]<-mdata$Species
bayestraitsdiscretedata[,2]<-as.numeric(as.factor(mdata$Dimorphic))-1
bayestraitsdiscretedata[,3]<- as.numeric(as.factor(mdata$FemaleDominance))-1
bayestraitsdiscretedata<-as.data.frame(bayestraitsdiscretedata)

command_vec1 <- c("2", "1") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
results_1 <- bayestraits(bayestraitsdiscretedata, mtree, command_vec1)
log_1 <- results_1$Log$results


command_vec2 <- c("3", "1","mltries 100") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
results_2 <- bayestraits(bayestraitsdiscretedata, mtree, command_vec2)
log_2 <- results_2$Log$results


# 1 = 0,0 not dimorphic, not female dominance - 2 = 0,1 not dimorphic, female dominance   - 3 = 1,0 dimorphic, not female dominance   4 = 1,1 dimorphic, female dominance
colnames(log_2)<-c("Tree","Likelihood","TransitionToFemaleDominance_AbsenceLowDimorphism","TransitionToLowDimorphism_AbsenceFemaleDominance","LossFemaleDominance_AbsenceLowDimorphism","TransitionToLowDimorphism_PresenceFemaleDominance","LossLowDimorphism_AbsenceFemaleDominance","TransitionToFemaleDominance_PresenceLowDimorphism","LossLowDimorphism_PresenceFemaleDominance","LossFemaleDominance_PresenceLowDimorphism","LikelihoodRootNotLowDimorphismNotFemaleDominance","LikelihoodRootNotLowDimorphismFemaleDominance","LikelihoodRootLowDimorphismNotFemaleDominance","LikelihoodRootLowDimorphismFemaleDominance")


log_2$SpeciesCount_NotDimorphic_NotFemaleDominant<-as.data.frame(bayestraitsdiscretedata %>% group_by(Dimorphic,FemaleDominance) %>% summarise(n()))[1,3]
log_2$SpeciesCount_NotDimorphic_FemaleDominant<-as.data.frame(bayestraitsdiscretedata %>% group_by(Dimorphic,FemaleDominance) %>% summarise(n()))[2,3]
log_2$SpeciesCount_Dimorphic_NotFemaleDominant<-as.data.frame(bayestraitsdiscretedata %>% group_by(Dimorphic,FemaleDominance) %>% summarise(n()))[3,3]
log_2$SpeciesCount_Dimorphic_FemaleDominant<-as.data.frame(bayestraitsdiscretedata %>% group_by(Dimorphic,FemaleDominance) %>% summarise(n()))[4,3]


log_2$Positive_IndependentModelBetter<-2*(log_2$Likelihood-log_1$Lh)

write.csv(log_2,file="results/Bayestraits_Dependent_RelaxedFemaleDominance_LowDimorphism.csv")




# sexual receptivity

# Testing whether there is dependent evolution of arboreality and female dominance

mostly_females_phylodata_1<-combined
mostly_females_phylodata_1$female_dominance<-as.numeric(as.factor(combined$mostlyfdom))
mostly_females_phylodata_1<-as.data.frame(mostly_females_phylodata_1 %>% group_by(corrected_species_id) %>% summarise(mean(female_dominance,na.rm=T )) )
colnames(mostly_females_phylodata_1)<-c("corrected_species_id", "femaledominance")
mostly_females_phylodata_1$femaledominance<-as.numeric(mostly_females_phylodata_1$femaledominance>1.49)


arboreal<-combined
arboreal<-arboreal[is.na(arboreal$Strata_Wilman)==F,]
arboreal$Strata_Wilman<-as.numeric(arboreal$Strata_Wilman=="Ar")
arboreal_phylodata<-as.data.frame(arboreal %>% group_by(corrected_species_id) %>% summarise(mean(Strata_Wilman,na.rm=T )) )

mostly_females_phylodata<-left_join(mostly_females_phylodata_1,arboreal_phylodata,by="corrected_species_id")
rownames(mostly_females_phylodata)<-mostly_females_phylodata$corrected_species_id
colnames(mostly_females_phylodata)<-c("Species","FemaleDominance","Arboreal")
mostly_females_phylodata<-mostly_females_phylodata[is.na(mostly_females_phylodata$Arboreal)==F,]
mostly_females_phylodata<-mostly_females_phylodata[is.na(mostly_females_phylodata$FemaleDominance)==F,]


# We match the species in the dataset to the species in the phylogenetic tree
missing<-treedata(inputtree,data=mostly_females_phylodata,warnings=FALSE)
# We remove species with no data from the tree
mtree<-missing$phy
speciesnames<-mtree$tip.label

mdata<-mostly_females_phylodata[speciesnames,]


bayestraitsdiscretedata<-matrix(NA,nrow=nrow(mdata),ncol=3)
colnames(bayestraitsdiscretedata)<-c("Species","Arboreal","FemaleDominance")
rownames(bayestraitsdiscretedata)<-mdata$Species
bayestraitsdiscretedata[,1]<-mdata$Species
bayestraitsdiscretedata[,2]<-as.numeric(as.factor(mdata$Arboreal))-1
bayestraitsdiscretedata[,3]<- as.numeric(as.factor(mdata$FemaleDominance))-1
bayestraitsdiscretedata<-as.data.frame(bayestraitsdiscretedata)

command_vec1 <- c("2", "1", "mltries 100") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
results_1 <- bayestraits(bayestraitsdiscretedata, mtree, command_vec1)
log_1 <- results_1$Log$results


command_vec2 <- c("3", "1","mltries 100") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
results_2 <- bayestraits(bayestraitsdiscretedata, mtree, command_vec2)
log_2 <- results_2$Log$results


# 1 = 0,0 not dimorphic, not female dominance - 2 = 0,1 not dimorphic, female dominance   - 3 = 1,0 dimorphic, not female dominance   4 = 1,1 dimorphic, female dominance
colnames(log_2)<-c("Tree","Likelihood","TransitionToFemaleDominance_AbsenceArboreal","TransitionToArboreal_AbsenceFemaleDominance","LossFemaleDominance_AbsenceArboreal","TransitionToArboreal_PresenceFemaleDominance","LossArboreal_AbsenceFemaleDominance","TransitionToFemaleDominance_PresenceArboreal","LossArboreal_PresenceFemaleDominance","LossFemaleDominance_PresenceArboreal","LikelihoodRootNotArborealNotFemaleDominance","LikelihoodRootNotArborealFemaleDominance","LikelihoodRootArborealNotFemaleDominance","LikelihoodRootLowDimorphismFemaleDominance")


log_2$SpeciesCount_NotDimorphic_NotFemaleDominant<-as.data.frame(bayestraitsdiscretedata %>% group_by(Arboreal,FemaleDominance) %>% summarise(n()))[1,3]
log_2$SpeciesCount_NotDimorphic_FemaleDominant<-as.data.frame(bayestraitsdiscretedata %>% group_by(Arboreal,FemaleDominance) %>% summarise(n()))[2,3]
log_2$SpeciesCount_Dimorphic_NotFemaleDominant<-as.data.frame(bayestraitsdiscretedata %>% group_by(Arboreal,FemaleDominance) %>% summarise(n()))[3,3]
log_2$SpeciesCount_Dimorphic_FemaleDominant<-as.data.frame(bayestraitsdiscretedata %>% group_by(Arboreal,FemaleDominance) %>% summarise(n()))[4,3]


log_2$Positive_IndependentModelBetter<-2*(log_2$Likelihood-log_1$Lh)

write.csv(log_2,file="results/Bayestraits_Dependent_RelaxedFemaleDominance_Arboreality.csv")



# Testing whether there is dependent evolution of seasonality and female dominance

mostly_females_phylodata_1<-combined
mostly_females_phylodata_1$female_dominance<-as.numeric(as.factor(combined$mostlyfdom))
mostly_females_phylodata_1<-as.data.frame(mostly_females_phylodata_1 %>% group_by(corrected_species_id) %>% summarise(mean(female_dominance,na.rm=T )) )
colnames(mostly_females_phylodata_1)<-c("corrected_species_id", "femaledominance")
mostly_females_phylodata_1$femaledominance<-as.numeric(mostly_females_phylodata_1$femaledominance>1.49)


arboreal<-combined
arboreal<-arboreal[is.na(arboreal$Strata_Wilman)==F,]
arboreal$Strata_Wilman<-as.numeric(arboreal$Strata_Wilman=="Ar")
arboreal_phylodata<-as.data.frame(arboreal %>% group_by(corrected_species_id) %>% summarise(mean(Strata_Wilman,na.rm=T )) )

mostly_females_phylodata<-left_join(mostly_females_phylodata_1,arboreal_phylodata,by="corrected_species_id")
rownames(mostly_females_phylodata)<-mostly_females_phylodata$corrected_species_id
colnames(mostly_females_phylodata)<-c("Species","FemaleDominance","Arboreal")
mostly_females_phylodata<-mostly_females_phylodata[is.na(mostly_females_phylodata$Arboreal)==F,]
mostly_females_phylodata<-mostly_females_phylodata[is.na(mostly_females_phylodata$FemaleDominance)==F,]


# We match the species in the dataset to the species in the phylogenetic tree
missing<-treedata(inputtree,data=mostly_females_phylodata,warnings=FALSE)
# We remove species with no data from the tree
mtree<-missing$phy
speciesnames<-mtree$tip.label

mdata<-mostly_females_phylodata[speciesnames,]


bayestraitsdiscretedata<-matrix(NA,nrow=nrow(mdata),ncol=3)
colnames(bayestraitsdiscretedata)<-c("Species","Arboreal","FemaleDominance")
rownames(bayestraitsdiscretedata)<-mdata$Species
bayestraitsdiscretedata[,1]<-mdata$Species
bayestraitsdiscretedata[,2]<-as.numeric(as.factor(mdata$Arboreal))-1
bayestraitsdiscretedata[,3]<- as.numeric(as.factor(mdata$FemaleDominance))-1
bayestraitsdiscretedata<-as.data.frame(bayestraitsdiscretedata)

command_vec1 <- c("2", "1", "mltries 100") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
results_1 <- bayestraits(bayestraitsdiscretedata, mtree, command_vec1)
log_1 <- results_1$Log$results


command_vec2 <- c("3", "1","mltries 100") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
results_2 <- bayestraits(bayestraitsdiscretedata, mtree, command_vec2)
log_2 <- results_2$Log$results


# 1 = 0,0 not dimorphic, not female dominance - 2 = 0,1 not dimorphic, female dominance   - 3 = 1,0 dimorphic, not female dominance   4 = 1,1 dimorphic, female dominance
colnames(log_2)<-c("Tree","Likelihood","TransitionToFemaleDominance_AbsenceArboreal","TransitionToArboreal_AbsenceFemaleDominance","LossFemaleDominance_AbsenceArboreal","TransitionToArboreal_PresenceFemaleDominance","LossArboreal_AbsenceFemaleDominance","TransitionToFemaleDominance_PresenceArboreal","LossArboreal_PresenceFemaleDominance","LossFemaleDominance_PresenceArboreal","LikelihoodRootNotArborealNotFemaleDominance","LikelihoodRootNotArborealFemaleDominance","LikelihoodRootArborealNotFemaleDominance","LikelihoodRootLowDimorphismFemaleDominance")


log_2$SpeciesCount_NotDimorphic_NotFemaleDominant<-as.data.frame(bayestraitsdiscretedata %>% group_by(Arboreal,FemaleDominance) %>% summarise(n()))[1,3]
log_2$SpeciesCount_NotDimorphic_FemaleDominant<-as.data.frame(bayestraitsdiscretedata %>% group_by(Arboreal,FemaleDominance) %>% summarise(n()))[2,3]
log_2$SpeciesCount_Dimorphic_NotFemaleDominant<-as.data.frame(bayestraitsdiscretedata %>% group_by(Arboreal,FemaleDominance) %>% summarise(n()))[3,3]
log_2$SpeciesCount_Dimorphic_FemaleDominant<-as.data.frame(bayestraitsdiscretedata %>% group_by(Arboreal,FemaleDominance) %>% summarise(n()))[4,3]


log_2$Positive_IndependentModelBetter<-2*(log_2$Likelihood-log_1$Lh)

write.csv(log_2,file="results/Bayestraits_Dependent_RelaxedFemaleDominance_Arboreality.csv")






##########################################################################################

library(phytools)
library(OUwie)

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

phenogram(matingsystemtree,perc_wondata)

#   1 MON        2 PAN            3 POL              4 PRO
#  "black" "#DF536B"/pink "#61D04F"/green "#2297E6"/blue


matingsystem_phylo[matingsystem_phylo=="MON"]<-1
matingsystem_phylo[matingsystem_phylo=="PRO"]<-2
matingsystem_phylo[matingsystem_phylo=="PAN"]<-3
matingsystem_phylo[matingsystem_phylo=="POL"]<-4
socialorganisationrecon<-rerootingMethod(mtree, matingsystem_phylo, model="ARD")
socialorganisationmtr<-mtree 
socialorganisationmtr$node.label<-vector() # set up an empty vector for the node labels


# this for loop goes through every node (i) and fills the empty vector socialorganisationmtr$node.label with the reconstructed most likely ancestral state at each node. We are assigning the state with the highest likelihood to each node which is identified by max(socialorganisationrecon$lik.anc[i,]) we then find which trait it belongs to by finding the column name which corresponds to the maximum estimate.
for (i in 1:mtree$Nnode) socialorganisationmtr$node.label[i]<-names(socialorganisationrecon$marginal.anc[i,])[socialorganisationrecon$marginal.anc[i,]==max(socialorganisationrecon$marginal.anc[i,])] 

complexitysocialorganisationresults<-list() 

# We set up a vector containing the names of the potential evolutionary models we want to run and use this in a simple loop.
# The models assume that there is no selection on complexity (BM1), no relationship between social organisation and complexity (OU1), or that the optimum for complexity or the rate of change in complexity depends on whether an ancestral species had a specific social organisation
mods<-c("BM1", "OU1", "BMS", "OUM") 

complexitydata<-mdata[,c(2,3,1)]
complexitydata$MatSysPMK<-matingsystem_phylo
complexitydata$perc_won_females<-standardize(complexitydata$perc_won_females)

for(i in 1:length(mods)){
  complexitysocialorganisationresults[[i]]<-OUwie(socialorganisationmtr, complexitydata, model=mods[i], simmap.tree=FALSE, root.station=FALSE)
}


complexitybestfit<-c(complexitysocialorganisationresults[[1]]$AICc, complexitysocialorganisationresults[[2]]$AICc, complexitysocialorganisationresults[[3]]$AICc, complexitysocialorganisationresults[[4]]$AICc)
names(complexitybestfit)<-mods # results are in the same order in which we ran the models
complexitybestfit-min(complexitybestfit)

estimate_complexity_eachorganisation<-complexitysocialorganisationresults[[4]]$theta
rownames(estimate_complexity_eachorganisation)<-unique(matingsystem_phylo)
estimate_complexity_eachorganisation[,1]*sd(mdata[,1])+mean(mdata[,1])
#        MON         PRO         PAN         POL 
#      94.810463 52.508275 56.193474     2.250075 








# Testing whether there is dependent evolution of UMMFs and male dominance

mostly_females_phylodata_1<-combined
mostly_females_phylodata_1$mostlyfdom<-as.numeric(as.factor(combined$mostlyfdom))
mostly_females_phylodata_1<-as.data.frame(mostly_females_phylodata_1 %>% group_by(corrected_species_id) %>% summarise(mean(mostlyfdom,na.rm=T )) )
mostly_females_phylodata_2<-as.data.frame(mostly_females_phylodata_1[,2])
row.names(mostly_females_phylodata_2)<-mostly_females_phylodata_1[,1]
colnames(mostly_females_phylodata_2)<-"mostlyfdom"
mostly_females_phylodata_3<-as.data.frame(mostly_females_phylodata_2[is.na(mostly_females_phylodata_2$mostlyfdom)==FALSE,])
row.names(mostly_females_phylodata_3)<-row.names(mostly_females_phylodata_2)[is.na(mostly_females_phylodata_2$mostlyfdom)==FALSE]
colnames(mostly_females_phylodata_3)<-"mostlyfdom"

UMMF<-combined
UMMF$singlemale<-as.numeric(UMMF$males==1)
UMMF$multiplefemales<-as.numeric(UMMF$females>1)
UMMF$UMMF<-as.numeric(UMMF$singlemale+UMMF$multiplefemales==2)
UMMF<-UMMF[is.na(UMMF$UMMF)==F,]
UMMF_phylodata<-as.data.frame(UMMF %>% group_by(corrected_species_id) %>% summarise(mean(UMMF,na.rm=T )) )
mostly_females_phylodata<-left_join(mostly_females_phylodata_1,UMMF_phylodata,by="corrected_species_id")
rownames(mostly_females_phylodata)<-mostly_females_phylodata$corrected_species_id
colnames(mostly_females_phylodata)<-c("Species","RelaxedFemdom","UMMF")
mostly_females_phylodata<-mostly_females_phylodata[is.na(mostly_females_phylodata$UMMF)==F,]
mostly_females_phylodata<-mostly_females_phylodata[is.na(mostly_females_phylodata$RelaxedFemdom)==F,]
mostly_females_phylodata$UMMF<-round(mostly_females_phylodata$UMMF,0)
mostly_females_phylodata$RelaxedFemdom<-round(mostly_females_phylodata$RelaxedFemdom,0)


#  http://www.phytools.org/static.help/add.species.to.genus.html

# We match the species in the dataset to the species in the phylogenetic tree
missing<-treedata(inputtree,data=mostly_females_phylodata,warnings=FALSE)
# We remove species with no data from the tree
mtree<-missing$phy
speciesnames<-mtree$tip.label

mdata<-mostly_females_phylodata[ speciesnames,]
mdata$RelaxedFemdom<-round(mdata$RelaxedFemdom,0)

bayestraitsdiscretedata<-matrix(NA,nrow=nrow(mdata),ncol=3)
colnames(bayestraitsdiscretedata)<-c("Species","UMMF","RelaxedFemdom")
rownames(bayestraitsdiscretedata)<-mtree$tip.label
bayestraitsdiscretedata[,1]<-mtree$tip.label
bayestraitsdiscretedata[,2]<-as.numeric(as.factor(mdata$UMMF))-1
bayestraitsdiscretedata[,3]<- as.numeric(as.factor(mdata$RelaxedFemdom))-1
bayestraitsdiscretedata<-as.data.frame(bayestraitsdiscretedata)

command_vec1 <- c("2", "1") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
results_1 <- bayestraits(bayestraitsdiscretedata, mtree, command_vec1)
log_1 <- results_1$Log
log_1$results$Lh

command_vec2 <- c("3", "1") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
results_2 <- bayestraits(bayestraitsdiscretedata, mtree, command_vec2)
log_2 <- results_2$Log
log_2$results$Lh

# 1 = 0,0 (OtherSocial + Female Dominance) - 2 = 0,1   - 3 = 1,0    4 = 1,1 (UMMF + Male Dominance)
# Female dominance only in other social systems, male dominance mainly in UMMF, UMMF more likely when already male dominance

# $results
#  Tree.No        Lh      q12      q13      q21      q24      q31      q34      q42 q43 Root...P.0.0.
# 1       1 -78.81987 0.000755 1.145145 0.029404 29.47959 37.89048 1.216348 99.15014   0      0.256331
#  Root...P.0.1. Root...P.1.0. Root...P.1.1.
# 1       0.24376      0.256154      0.243755


# 15 UMMF species, 78 with other mating system









# Testing whether there is dependent evolution of seasonality and female dominance

mostly_females_phylodata_1<-combined
mostly_females_phylodata_1$mostlyfdom<-as.numeric(as.factor(combined$mostlyfdom))
mostly_females_phylodata_1<-as.data.frame(mostly_females_phylodata_1 %>% group_by(corrected_species_id) %>% summarise(mean(mostlyfdom,na.rm=T )) )
mostly_females_phylodata_2<-as.data.frame(mostly_females_phylodata_1[,2])
row.names(mostly_females_phylodata_2)<-mostly_females_phylodata_1[,1]
colnames(mostly_females_phylodata_2)<-"mostlyfdom"
mostly_females_phylodata_3<-as.data.frame(mostly_females_phylodata_2[is.na(mostly_females_phylodata_2$mostlyfdom)==FALSE,])
row.names(mostly_females_phylodata_3)<-row.names(mostly_females_phylodata_2)[is.na(mostly_females_phylodata_2$mostlyfdom)==FALSE]
colnames(mostly_females_phylodata_3)<-"mostlyfdom"

seasonality<-combined
seasonality<-seasonality[is.na(seasonality$r_seasonality_value)==F,]
seasonality[seasonality$r_seasonality_value>0.5,]$r_seasonality_value<- 1
seasonality[seasonality$r_seasonality_value<0.51,]$r_seasonality_value<- 0
seasonality_phylodata<-as.data.frame(seasonality %>% group_by(corrected_species_id) %>% summarise(mean(r_seasonality_value,na.rm=T )) )

mostly_females_phylodata<-left_join(mostly_females_phylodata_1,seasonality_phylodata,by="corrected_species_id")
rownames(mostly_females_phylodata)<-mostly_females_phylodata$corrected_species_id
colnames(mostly_females_phylodata)<-c("Species","RelaxedFemdom","Seasonality")
mostly_females_phylodata<-mostly_females_phylodata[is.na(mostly_females_phylodata$Seasonality)==F,]
mostly_females_phylodata<-mostly_females_phylodata[is.na(mostly_females_phylodata$RelaxedFemdom)==F,]


#  http://www.phytools.org/static.help/add.species.to.genus.html

# We match the species in the dataset to the species in the phylogenetic tree
missing<-treedata(inputtree,data=mostly_females_phylodata,warnings=FALSE)
# We remove species with no data from the tree
mtree<-missing$phy
speciesnames<-mtree$tip.label

mdata<-mostly_females_phylodata[mostly_females_phylodata$Species %in% speciesnames,]
mdata$RelaxedFemdom<-round(mdata$RelaxedFemdom,0)

bayestraitsdiscretedata<-matrix(NA,nrow=nrow(mdata),ncol=3)
colnames(bayestraitsdiscretedata)<-c("Species","Seasonality","RelaxedFemdom")
rownames(bayestraitsdiscretedata)<-mdata$Species
bayestraitsdiscretedata[,1]<-mdata$Species
bayestraitsdiscretedata[,2]<-as.numeric(as.factor(mdata$Seasonality))-1
bayestraitsdiscretedata[,3]<- as.numeric(as.factor(mdata$RelaxedFemdom))-1
bayestraitsdiscretedata<-as.data.frame(bayestraitsdiscretedata)

command_vec1 <- c("2", "1") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
results_1 <- bayestraits(bayestraitsdiscretedata, mtree, command_vec1)
log_1 <- results_1$Log
log_1$results$Lh

command_vec2 <- c("3", "1") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
results_2 <- bayestraits(bayestraitsdiscretedata, mtree, command_vec2)
log_2 <- results_2$Log
log_2

# 1 = 0,0  - 2 = 0,1   - 3 = 1,0    4 = 1,1
# Transition to female dominance 2-1 and 4-3: female dominance only when there is seasonality
# Transition to seasonality 1-3 and 2-4: only when there is male dominance
# might suggest that seasonality comes first

# 4 species with female dominance






# Testing whether there is dependent evolution of dimorphism and female dominance

mostly_females_phylodata_1<-combined
mostly_females_phylodata_1$mostlyfdom<-as.numeric(as.factor(combined$mostlyfdom))
mostly_females_phylodata_1<-as.data.frame(mostly_females_phylodata_1 %>% group_by(corrected_species_id) %>% summarise(mean(mostlyfdom,na.rm=T )) )
mostly_females_phylodata_2<-as.data.frame(mostly_females_phylodata_1[,2])
row.names(mostly_females_phylodata_2)<-mostly_females_phylodata_1[,1]
colnames(mostly_females_phylodata_2)<-"mostlyfdom"
mostly_females_phylodata_3<-as.data.frame(mostly_females_phylodata_2[is.na(mostly_females_phylodata_2$mostlyfdom)==FALSE,])
row.names(mostly_females_phylodata_3)<-row.names(mostly_females_phylodata_2)[is.na(mostly_females_phylodata_2$mostlyfdom)==FALSE]
colnames(mostly_females_phylodata_3)<-"mostlyfdom"

dimorphism<-combined
dimorphism<-dimorphism[is.na(dimorphism$SexualDimorphism_MaleWeight_over_FemaleWeight)==F,]
dimorphism[dimorphism$SexualDimorphism_MaleWeight_over_FemaleWeight<1.101,]$SexualDimorphism_MaleWeight_over_FemaleWeight<- 0
dimorphism[dimorphism$SexualDimorphism_MaleWeight_over_FemaleWeight>1.1,]$SexualDimorphism_MaleWeight_over_FemaleWeight<- 1
dimorphism_phylodata<-as.data.frame(dimorphism %>% group_by(corrected_species_id) %>% summarise(mean(SexualDimorphism_MaleWeight_over_FemaleWeight,na.rm=T )) )

mostly_females_phylodata<-left_join(mostly_females_phylodata_1,dimorphism_phylodata,by="corrected_species_id")
rownames(mostly_females_phylodata)<-mostly_females_phylodata$corrected_species_id
colnames(mostly_females_phylodata)<-c("Species","RelaxedFemdom","Dimorphism")
mostly_females_phylodata<-mostly_females_phylodata[is.na(mostly_females_phylodata$Dimorphism)==F,]
mostly_females_phylodata<-mostly_females_phylodata[is.na(mostly_females_phylodata$RelaxedFemdom)==F,]


#  http://www.phytools.org/static.help/add.species.to.genus.html

# We match the species in the dataset to the species in the phylogenetic tree
missing<-treedata(inputtree,data=mostly_females_phylodata,warnings=FALSE)
# We remove species with no data from the tree
mtree<-missing$phy
speciesnames<-mtree$tip.label

mdata<-mostly_females_phylodata[mostly_females_phylodata$Species %in% speciesnames,]
mdata$RelaxedFemdom<-round(mdata$RelaxedFemdom,0)

bayestraitsdiscretedata<-matrix(NA,nrow=nrow(mdata),ncol=3)
colnames(bayestraitsdiscretedata)<-c("Species","Dimorphism","RelaxedFemdom")
rownames(bayestraitsdiscretedata)<-mdata$Species
bayestraitsdiscretedata[,1]<-mdata$Species
bayestraitsdiscretedata[,2]<-as.numeric(as.factor(mdata$Dimorphism))-1
bayestraitsdiscretedata[,3]<- as.numeric(as.factor(mdata$RelaxedFemdom))-1
bayestraitsdiscretedata<-as.data.frame(bayestraitsdiscretedata)

command_vec1 <- c("2", "1") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
results_1 <- bayestraits(bayestraitsdiscretedata, mtree, command_vec1)
log_1 <- results_1$Log
log_1$results$Lh

command_vec2 <- c("3", "1") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
results_2 <- bayestraits(bayestraitsdiscretedata, mtree, command_vec2)
log_2 <- results_2$Log
log_2

# 1 = 0,0 (no dimorphism, female dominance)  - 2 = 0,1   - 3 = 1,0    4 = 1,1 (dimorphism, male dominance)
# when there is dimorphism, female dominance does not evolve - male dominance more likely to evolve
# dimorphism more likely to evolve when there is male dominance, to be lost when there is female dominance
# suggests co-evolutionary pattern


# 25 species with female dominance / 50 without
# 31 species without dimorphism / 44 species with

# For the phenogram

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

phenogram(dimorphismtree,perc_wondata)

#   1 MON        2 PAN            3 POL              4 PRO
#  "black" "#DF536B"/pink "#61D04F"/green "#2297E6"/blue




# Testing whether there is dependent evolution of arboreality and female dominance

mostly_females_phylodata_1<-combined
mostly_females_phylodata_1$mostlyfdom<-as.numeric(as.factor(combined$mostlyfdom))
mostly_females_phylodata_1<-as.data.frame(mostly_females_phylodata_1 %>% group_by(corrected_species_id) %>% summarise(mean(mostlyfdom,na.rm=T )) )
colnames(mostly_females_phylodata_1)<-c("corrected_species_id","mostlyfdom")
rownames(mostly_females_phylodata_1)<-mostly_females_phylodata_1$corrected_species_id
mostly_females_phylodata_1<-mostly_females_phylodata_1[is.na(mostly_females_phylodata_1$mostlyfdom)==F,]

arboreal<-select(combined,corrected_species_id,Strata_Wilman)
arboreal<-arboreal[is.na(arboreal$Strata_Wilman)==F,]
arboreal$arboreality<-0
arboreal[arboreal$Strata_Wilman=="Ar",]$arboreality<-1
arboreal_phylo<-as.data.frame(arboreal %>% group_by(corrected_species_id) %>% summarise(mean(arboreality,na.rm=T )) )
colnames(arboreal_phylo)<-c("corrected_species_id","arboreality")

arboreal_phylodata<-inner_join(arboreal_phylo,mostly_females_phylodata_1,by="corrected_species_id")

rownames(arboreal_phylodata)<-arboreal_phylodata$corrected_species_id
colnames(arboreal_phylodata)<-c("corrected_species_id","Arboreality","RelaxedFemdom")



#  http://www.phytools.org/static.help/add.species.to.genus.html

# We match the species in the dataset to the species in the phylogenetic tree
missing<-treedata(inputtree,data=arboreal_phylodata,warnings=FALSE)
# We remove species with no data from the tree
mtree<-missing$phy
speciesnames<-mtree$tip.label

mdata<-arboreal_phylodata[arboreal_phylodata$corrected_species_id %in% speciesnames,]
mdata$RelaxedFemdom<-round(mdata$RelaxedFemdom,0)

bayestraitsdiscretedata<-matrix(NA,nrow=nrow(mdata),ncol=3)
colnames(bayestraitsdiscretedata)<-c("Species","Arboreality","RelaxedFemdom")
rownames(bayestraitsdiscretedata)<-mdata$corrected_species_id
bayestraitsdiscretedata[,1]<-mdata$corrected_species_id
bayestraitsdiscretedata[,2]<-as.numeric(as.factor(mdata$Arboreality))-1
bayestraitsdiscretedata[,3]<- as.numeric(as.factor(mdata$RelaxedFemdom))-1
bayestraitsdiscretedata<-as.data.frame(bayestraitsdiscretedata)

command_vec1 <- c("2", "1") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
results_1 <- bayestraits(bayestraitsdiscretedata, mtree, command_vec1)
log_1 <- results_1$Log
log_1$results$Lh

command_vec2 <- c("3", "1") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
results_2 <- bayestraits(bayestraitsdiscretedata, mtree, command_vec2)
log_2 <- results_2$Log
log_2

# 1 = 0,0 (no arboreality, male dominance)  - 2 = 0,1   - 3 = 1,0    4 = 1,1 (arboreality, female dominance)
# when there is dimorphism, female dominance does not evolve - male dominance more likely to evolve
# dimorphism more likely to evolve when there is male dominance, to be lost when there is female dominance
# suggests co-evolutionary pattern


# 25 species with female dominance / 50 without
# 31 species without dimorphism / 44 species with

# For the phenogram

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

phenogram(arbotree,perc_wondata)

#   1 MON        2 PAN            3 POL              4 PRO
#  "black" "#DF536B"/pink "#61D04F"/green "#2297E6"/blue










# Testing whether there is dependent evolution of dimorphism and UMMF

UMMF<-combined
UMMF$singlemale<-as.numeric(UMMF$males==1)
UMMF$multiplefemales<-as.numeric(UMMF$females>1)
UMMF$UMMF<-as.numeric(UMMF$singlemale+UMMF$multiplefemales==2)
UMMF<-UMMF[is.na(UMMF$UMMF)==F,]
UMMF_phylodata<-as.data.frame(UMMF %>% group_by(corrected_species_id) %>% summarise(mean(UMMF,na.rm=T )) )


dimorphism<-combined
dimorphism<-dimorphism[is.na(dimorphism$SexualDimorphism_MaleWeight_over_FemaleWeight)==F,]
dimorphism[dimorphism$SexualDimorphism_MaleWeight_over_FemaleWeight<1.101,]$SexualDimorphism_MaleWeight_over_FemaleWeight<- 0
dimorphism[dimorphism$SexualDimorphism_MaleWeight_over_FemaleWeight>1.1,]$SexualDimorphism_MaleWeight_over_FemaleWeight<- 1
dimorphism_phylodata<-as.data.frame(dimorphism %>% group_by(corrected_species_id) %>% summarise(mean(SexualDimorphism_MaleWeight_over_FemaleWeight,na.rm=T )) )

mostly_females_phylodata<-left_join(UMMF_phylodata,dimorphism_phylodata,by="corrected_species_id")
rownames(mostly_females_phylodata)<-mostly_females_phylodata$corrected_species_id
colnames(mostly_females_phylodata)<-c("Species","UMMF","Dimorphism")
mostly_females_phylodata<-mostly_females_phylodata[is.na(mostly_females_phylodata$Dimorphism)==F,]
mostly_females_phylodata<-mostly_females_phylodata[is.na(mostly_females_phylodata$UMMF)==F,]
mostly_females_phylodata$UMMF<-round(mostly_females_phylodata$UMMF,0)


#  http://www.phytools.org/static.help/add.species.to.genus.html

# We match the species in the dataset to the species in the phylogenetic tree
missing<-treedata(inputtree,data=mostly_females_phylodata,warnings=FALSE)
# We remove species with no data from the tree
mtree<-missing$phy
speciesnames<-mtree$tip.label

mdata<-mostly_females_phylodata[mostly_females_phylodata$Species %in% speciesnames,]


bayestraitsdiscretedata<-matrix(NA,nrow=nrow(mdata),ncol=3)
colnames(bayestraitsdiscretedata)<-c("Species","Dimorphism","UMMF")
rownames(bayestraitsdiscretedata)<-mdata$Species
bayestraitsdiscretedata[,1]<-mdata$Species
bayestraitsdiscretedata[,2]<-as.numeric(as.factor(mdata$Dimorphism))-1
bayestraitsdiscretedata[,3]<- as.numeric(as.factor(mdata$UMMF))-1
bayestraitsdiscretedata<-as.data.frame(bayestraitsdiscretedata)

command_vec1 <- c("2", "1") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
results_1 <- bayestraits(bayestraitsdiscretedata, mtree, command_vec1)
log_1 <- results_1$Log
log_1$results$Lh

command_vec2 <- c("3", "1") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
results_2 <- bayestraits(bayestraitsdiscretedata, mtree, command_vec2)
log_2 <- results_2$Log
log_2$results$Lh

# 1 = 0,0 (no dimorphism, no UMMF)  - 2 = 0,1   - 3 = 1,0    4 = 1,1 (dimorphism, UMMF)
# UMMF only evolves when there is dimorphism, dimorphism evovles before



# 25 species with female dominance / 50 without
# 31 species without dimorphism / 44 species with