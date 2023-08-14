# Multivariate analyses


# seasonality + group sex ratio

dstan_continuous <- combined[ complete.cases(combined$perc_won_females,combined$r_seasonality_value,combined$sexratio,combined$corrected_species_id),]

# Continuous outcome: percentage of fights won by females
dat_list_continuous <- list(
  N_spp = nrow(dstan_continuous),
  perc_won_females = as.integer(dstan_continuous$perc_won_females),
  r_seasonality_value = standardize(dstan_continuous$r_seasonality_value),
  sexratio = standardize(dstan_continuous$sexratio),
  species = as.integer(as.factor(dstan_continuous$corrected_species_id)),
  total = rep(100,nrow(dstan_continuous))
)

m_continuous_continuous <- ulam(
  alist(
    perc_won_females ~ dbinom(total,p),
    logit(p) <- a +b*r_seasonality_value+c*sexratio,
    ## adaptive priors
    a ~ dnorm( 0 , 1 ),
    b ~ dnorm(0,1),
    c ~ dnorm(0,1)
  ) , data=dat_list_continuous , chains=4 , cores=4 , log_lik=TRUE , cmdstan=T)

precis(m_continuous_continuous,depth=2)

# Seasonality and sex ratio independently influence how many fights females win


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
  r_seasonality_value = standardize(dstan_continuous$r_seasonality_value),
  sexratio = standardize(dstan_continuous$sexratio),
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
    logit(p) <- k[species] +b*sexratio+c*r_seasonality_value,
    vector[N_spp]:k ~ multi_normal( 0 , SIGMA ),
    matrix[N_spp,N_spp]: SIGMA <- cov_GPL2( Dmat, etasq, rhosq, 0.01),
    ## adaptive priors
    b ~ dnorm(0,1),
    c ~ dnorm(0,1),
    etasq~exponential(1),
    rhosq~exponential(1)
  ) , data=dat_list_continuous_phylogenetic , chains=4 , cores=4 , log_lik=TRUE , cmdstan=T)

precis(m_continuous_continuous_phylogenetic,depth=2)


dat_list_continuous_phylogenetic <- list(
  N_spp = nrow(dstan_continuous),
  perc_won_females = standardize(as.integer(dstan_continuous$perc_won_females)),
  r_seasonality_value = dstan_continuous$r_seasonality_value*100,
  sexratio = standardize(dstan_continuous$sexratio),
  species = as.integer(as.factor(dstan_continuous$corrected_species_id)),
  total = rep(100,nrow(dstan_continuous))
)

dat_list_continuous_phylogenetic$Dmat<-Dmat[ unique(dstan_continuous$corrected_species_id),unique(dstan_continuous$corrected_species_id) ]/max(Dmat)

colnames(dat_list_continuous_phylogenetic$Dmat)<-as.integer(as.factor(colnames(dat_list_continuous_phylogenetic$Dmat)))
rownames(dat_list_continuous_phylogenetic$Dmat)<-as.integer(as.factor(rownames(dat_list_continuous_phylogenetic$Dmat)))

dat_list_continuous_phylogenetic$N_spp<-length(unique(dstan_continuous$corrected_species_id))


m_seasonality_continuous_phylogenetic <- ulam(
  alist(
    r_seasonality_value ~ dbinom(total,p),
    logit(p) <- k[species] +b*sexratio+c*perc_won_females,
    vector[N_spp]:k ~ multi_normal( 0 , SIGMA ),
    matrix[N_spp,N_spp]: SIGMA <- cov_GPL2( Dmat, etasq, rhosq, 0.01),
    ## adaptive priors
    b ~ dnorm(0,1),
    c ~ dnorm(0,1),
    etasq~exponential(1),
    rhosq~exponential(1)
  ) , data=dat_list_continuous_phylogenetic , chains=4 , cores=4 , log_lik=TRUE , cmdstan=T)

precis(m_seasonality_continuous_phylogenetic,depth=2)

# only sex ratio linked to seasonality but not percent won females



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
