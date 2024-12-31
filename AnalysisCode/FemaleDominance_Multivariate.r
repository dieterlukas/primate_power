# Multivariate analyses


# sexual receptivity + reproductive synchrony + females

dstan_continuous <- combined[ complete.cases(combined$perc_won_females,combined$receptive_synchrony,combined$sexualreceptivity_hours,combined$females,combined$corrected_species_id),]


spp_obs<-unique(dstan_continuous$corrected_species_id)
spp_obs<-matrix(nrow=length(spp_obs))
rownames(spp_obs)<-unique(dstan_continuous$corrected_species_id)

# We match the species in the dataset to the species in the phylogenetic tree
missing<-treedata(inputtree,data=spp_obs,warnings=FALSE)
# We remove species with no data from the tree
mtree<-missing$phy

# We calculate the pair-wise phylogenetic distances among the species in the tree; this gives a matrix where values are the total branch length needed to connect two species (diagonal is the distance of a species to itself so zero)
Dmat<-cophenetic(mtree)

dat_list_continuous_phylogenetic <- list(
  N_spp = nrow(dstan_continuous),
  perc_won_females = as.integer(dstan_continuous$perc_won_females),
  sexualreceptivity_hours = standardize(dstan_continuous$sexualreceptivity_hours),
  receptive_synchrony = standardize(dstan_continuous$receptive_synchrony),
  females = standardize(dstan_continuous$females),
  species = as.integer(as.factor(dstan_continuous$corrected_species_id)),
  total = rep(100,nrow(dstan_continuous))
)

dat_list_continuous_phylogenetic$Dmat<-Dmat[ unique(dstan_continuous$corrected_species_id),unique(dstan_continuous$corrected_species_id) ]/max(Dmat)

colnames(dat_list_continuous_phylogenetic$Dmat)<-as.integer(as.factor(colnames(dat_list_continuous_phylogenetic$Dmat)))
rownames(dat_list_continuous_phylogenetic$Dmat)<-as.integer(as.factor(rownames(dat_list_continuous_phylogenetic$Dmat)))

dat_list_continuous_phylogenetic$N_spp<-length(unique(dstan_continuous$corrected_species_id))


m_continuous_continuous_phylogenetic <- ulam(
  alist(
    perc_won_females ~ dbinom(total,p),
    logit(p) <- k[species] +b*sexualreceptivity_hours+c*receptive_synchrony+d*females,
    vector[N_spp]:k ~ multi_normal( 0 , SIGMA ),
    matrix[N_spp,N_spp]: SIGMA <- cov_GPL2( Dmat, etasq, rhosq, 0.01),
    ## adaptive priors
    b ~ dnorm(0,1),
    c ~ dnorm(0,1),
    d ~ dnorm(0,1),
    etasq~exponential(1),
    rhosq~exponential(1)
  ) , data=dat_list_continuous_phylogenetic , chains=4 , cores=4 , log_lik=TRUE , cmdstan=T)

precis(m_continuous_continuous_phylogenetic,depth=2)



# arboreality + females

dstan_strict <- combined[ complete.cases(combined$strictfdom,combined$females,combined$Strata_Wilman,combined$corrected_species_id),]

spp_obs<-unique(dstan_strict$corrected_species_id)
spp_obs<-matrix(nrow=length(spp_obs))
rownames(spp_obs)<-unique(dstan_strict$corrected_species_id)

# We match the species in the dataset to the species in the phylogenetic tree
missing<-treedata(inputtree,data=spp_obs,warnings=FALSE)
# We remove species with no data from the tree
mtree<-missing$phy

# We calculate the pair-wise phylogenetic distances among the species in the tree; this gives a matrix where values are the total branch length needed to connect two species (diagonal is the distance of a species to itself so zero)
Dmat<-cophenetic(mtree)
dat_list_continuous_phylogenetic <- list(
  N_spp = nrow(dstan_strict),
  R = as.factor(dstan_strict$strictfdom),
  females = standardize(dstan_strict$females),
  arboreal = as.integer(dstan_strict$Strata_Wilman=="Ar"),
  species = as.integer(as.factor(dstan_strict$corrected_species_id)),
  total = rep(100,nrow(dstan_strict))
)

dat_list_continuous_phylogenetic$Dmat<-Dmat[ unique(dstan_strict$corrected_species_id),unique(dstan_strict$corrected_species_id) ]/max(Dmat)

colnames(dat_list_continuous_phylogenetic$Dmat)<-as.integer(as.factor(colnames(dat_list_continuous_phylogenetic$Dmat)))
rownames(dat_list_continuous_phylogenetic$Dmat)<-as.integer(as.factor(rownames(dat_list_continuous_phylogenetic$Dmat)))

dat_list_continuous_phylogenetic$N_spp<-length(unique(dstan_strict$corrected_species_id))


m_strict_continuous_phylogenetic <- ulam(
  alist(
    R ~ dordlogit( phi , cutpoints ),
    phi <- k[species] +b*females+c*arboreal,
    vector[N_spp]:k ~ multi_normal( 0 , SIGMA ),
    matrix[N_spp,N_spp]: SIGMA <- cov_GPL2( Dmat, etasq, rhosq, 0.01),
    b ~ dnorm(0,1),
    c~dnorm(0,1),
    cutpoints ~ dnorm( 0 , 1.5 ),
    etasq~exponential(1),
    rhosq~exponential(1)
  ) , data=dat_list_continuous_phylogenetic , chains=4 , cores=4 ,cmdstan=T, messages=FALSE, refresh=0)


precis(m_strict_continuous_phylogenetic)


