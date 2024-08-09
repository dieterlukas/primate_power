
##########################################################################################
# Calculation of phylogenetic signal

# First, for the percentage of fights won by females

# To calculate the phylogenetic signal in line with Blomberg's K, we can only have a single value for each species. We therefore take the average for each species:
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

# We calculate the phylogenetic signal
phylosignal_percwon<-phylosig(tree=perctree,x=percdata$perc_won_females,method="K",test=T)



# Calculation of phylogenetic signal

# Second, for the strict classification of the dominance system

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
phylosignal_strict<-phylosig(tree=stricttree,x=strictdata$strictfdom,method="K",test=T)


# Calculation of phylogenetic signal

phylosignals<-matrix(nrow=2,ncol=3)
phylosignals<-as.data.frame(phylosignals)
colnames(phylosignals)<-c("variable","phylogenetic_signal_Blomberg_K","significance")
phylosignals$variable<-c("Percentage fights won by females","Strict three-way classification")
phylosignals$phylogenetic_signal_Blomberg_K<-c(phylosignal_percwon$K,phylosignal_strict$K)
phylosignals$significance<-c(phylosignal_percwon$P,phylosignal_strict$P)
write.csv(phylosignals,file="./results/FemaleDominance_PhylogeneticSignal.csv")

##########################################################################################
# Estimation of the phylogenetic covariance using STAN models

# We need a list of species to match it to the phylogenetic tree
spp_obs<-unique(combined$corrected_species_id)
spp_obs<-matrix(nrow=length(spp_obs))
rownames(spp_obs)<-unique(combined$corrected_species_id)

# We match the species in the dataset to the species in the phylogenetic tree
missing<-treedata(inputtree,data=spp_obs,warnings=FALSE)
# We remove species with no data from the tree
mtree<-missing$phy
# We remove species not included in the tree from the dataset
speciesnames<-mtree$tip.label
mdata<-combined[combined$corrected_species_id %in% speciesnames,]

# We calculate the pair-wise phylogenetic distances among the species in the tree; this gives a matrix where values are the total branch length needed to connect two species (diagonal is the distance of a species to itself so zero)
Dmat<-cophenetic(mtree)

# First, we estimate the covariance for the percentage of fights won by females
# We remove entries where we do not have this quantitative data
perc_won_females_phylodata<-mdata[ complete.cases(mdata$perc_won_females,mdata$corrected_species_id),]

# We create the input for the STAN model. The STAN model estimates covariance in a regression, so we need a factor in the model. To make sure that this factor does not shape the results, we simulate random data and call this factor G.
dat_list <- list(
  N_spp = nrow(perc_won_females_phylodata),
  perc_won_females = as.integer(perc_won_females_phylodata$perc_won_females)/100,
  G = rep(1,nrow(perc_won_females_phylodata))
)
# We add the phylogenetic distance matrix to the input. The model estimates similarity among species, so calculate the relative distances with the maximum distance (species only connect at the root) to 1.
dat_list$Dmat<-Dmat[ perc_won_females_phylodata$corrected_species_id,perc_won_females_phylodata$corrected_species_id ]/max(Dmat)

# We define the model: the phylogenetic similarity matrix Dmat is assumed to follow a Gaussian process, with closely related species being more similar than distantly related species, but the change with phylogenetic distance does not have to be linear
rethinking_phylogenetic_gaussian_origin <- ulam(
  alist(
    vector[N_spp]:perc_won_females ~ multi_normal( mu , SIGMA ),
    mu <- a+bG*G,
    matrix[N_spp,N_spp]: SIGMA <- cov_GPL2( Dmat, etasq, rhosq, 0.01),
    a~normal(0,1),
    bG~normal(0,0.5),
    etasq~exponential(1),
    rhosq~exponential(1)
  ),data=dat_list,chains=4,cores=4,cmdstan=T)

# We extract the posterior Bayesian distribution of the different estimates
post <- extract.samples(rethinking_phylogenetic_gaussian_origin)
# Given that we have multiple parameters that interact to shape the Gaussian process, it is easiest to plot the predicted outcome

pdf("figures/R_Figure_PhylogeneticCovariance_PercFightsWon.pdf")
plot( NULL , xlim=c(0,max(dat_list$Dmat)) , ylim=c(0,3.5) ,
      xlab="phylogenetic distance" , ylab="covariance" )
# posterior
for ( i in 1:230 )
  curve( post$etasq[i]*exp(-post$rhosq[i]*x^2 ) , add=TRUE , col=rangi2 )

valuesacross<-rep(NA,101)
for (x in 0:100) {
  valueshere<-rep(NA,500)
  for ( i in 1:500 ) {
    valueshere[i]<-post$etasq[i]*exp(-post$rhosq[i]*(x/100)^2 )
  } 
  valuesacross[x+1]<-mean(valueshere)
}
lines(x=seq(from = 0, to =1, by=0.01), y= valuesacross, lwd=5)
dev.off()



# Second, we estimate the covariance for the strict classification of dominance systems

strictfdom_females_phylodata<-mdata[ complete.cases(mdata$strictfdom,mdata$corrected_species_id),]

dat_list <- list(
  N_spp = nrow(strictfdom_females_phylodata),
  strictfdom = standardize(as.integer(as.factor(strictfdom_females_phylodata$strictfdom))),
  G = rep(1,nrow(strictfdom_females_phylodata))
)

dat_list$Dmat<-Dmat[ strictfdom_females_phylodata$corrected_species_id,strictfdom_females_phylodata$corrected_species_id ]/max(Dmat)

rethinking_phylogenetic_gaussian_origin <- ulam(
  alist(
    vector[N_spp]:strictfdom ~ multi_normal( mu , SIGMA ),
    mu <- a+bG*G,
    matrix[N_spp,N_spp]: SIGMA <- cov_GPL2( Dmat, etasq, rhosq, 0.01),
    a~normal(0,1),
    bG~normal(0,0.5),
    etasq~exponential(1),
    rhosq~exponential(1)
  ),data=dat_list,chains=4,cores=4,cmdstan=T)

post <- extract.samples(rethinking_phylogenetic_gaussian_origin)

pdf("figures/R_Figure_PhylogeneticCovariance_StrictClassificationIntersexualDominance.pdf")
plot( NULL , xlim=c(0,max(dat_list$Dmat)) , ylim=c(0,30) ,
      xlab="phylogenetic distance" , ylab="covariance" )
# posterior
for ( i in 1:230 )
  curve( post$etasq[i]*exp(-post$rhosq[i]*x^2 ) , add=TRUE , col=rangi2 )
valuesacross<-rep(NA,101)
for (x in 0:100) {
  valueshere<-rep(NA,500)
  for ( i in 1:500 ) {
    valueshere[i]<-post$etasq[i]*exp(-post$rhosq[i]*(x/100)^2 )
  } 
  valuesacross[x+1]<-mean(valueshere)
}
lines(x=seq(from = 0, to =1, by=0.01), y= valuesacross, lwd=5)
dev.off()
