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






# canine dimorphism + females

dstan_strict <- combined[ complete.cases(combined$strictfdom,combined$females,combined$CanineDimorphism,combined$corrected_species_id),]

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
  dimorphism = standardize(dstan_strict$CanineDimorphism),
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
    phi <- k[species] +b*females+c*dimorphism,
    vector[N_spp]:k ~ multi_normal( 0 , SIGMA ),
    matrix[N_spp,N_spp]: SIGMA <- cov_GPL2( Dmat, etasq, rhosq, 0.01),
    b ~ dnorm(0,1),
    c~dnorm(0,1),
    cutpoints ~ dnorm( 0 , 1.5 ),
    etasq~exponential(1),
    rhosq~exponential(1)
  ) , data=dat_list_continuous_phylogenetic , chains=4 , cores=4 ,cmdstan=T, messages=FALSE, refresh=0)


precis(m_strict_continuous_phylogenetic)


# canine dimorphism + body size

dstan_strict <- combined[ complete.cases(combined$strictfdom,combined$body_mass,combined$CanineDimorphism,combined$corrected_species_id),]

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
  body_mass = standardize(dstan_strict$body_mass),
  dimorphism = standardize(dstan_strict$CanineDimorphism),
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
    phi <- k[species] +b*body_mass+c*dimorphism,
    vector[N_spp]:k ~ multi_normal( 0 , SIGMA ),
    matrix[N_spp,N_spp]: SIGMA <- cov_GPL2( Dmat, etasq, rhosq, 0.01),
    b ~ dnorm(0,1),
    c~dnorm(0,1),
    cutpoints ~ dnorm( 0 , 1.5 ),
    etasq~exponential(1),
    rhosq~exponential(1)
  ) , data=dat_list_continuous_phylogenetic , chains=4 , cores=4 ,cmdstan=T, messages=FALSE, refresh=0)


precis(m_strict_continuous_phylogenetic)





# mating system + sex bias in dispersal

dstan_strict <- combined[ complete.cases(combined$strictfdom,combined$sexbias_dispersal,combined$MatSysPMK,combined$corrected_species_id),]

dstan_strict$Polygyny<-dstan_strict$MatSysPMK=="POL"
dstan_strict[dstan_strict$sexbias_dispersal=="Both",]$sexbias_dispersal<-2
dstan_strict[dstan_strict$sexbias_dispersal=="Female",]$sexbias_dispersal<-1
dstan_strict[dstan_strict$sexbias_dispersal=="Male",]$sexbias_dispersal<-3

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
  dispersal = as.integer(dstan_strict$sexbias_dispersal),
  polygyny = as.integer(dstan_strict$Polygyny),
  species = as.integer(as.factor(dstan_strict$corrected_species_id))
)

dat_list_continuous_phylogenetic$Dmat<-Dmat[ unique(dstan_strict$corrected_species_id),unique(dstan_strict$corrected_species_id) ]/max(Dmat)

colnames(dat_list_continuous_phylogenetic$Dmat)<-as.integer(as.factor(colnames(dat_list_continuous_phylogenetic$Dmat)))
rownames(dat_list_continuous_phylogenetic$Dmat)<-as.integer(as.factor(rownames(dat_list_continuous_phylogenetic$Dmat)))

dat_list_continuous_phylogenetic$N_spp<-length(unique(dstan_strict$corrected_species_id))


m_strict_continuous_phylogenetic <- ulam(
  alist(
    R ~ dordlogit( phi , cutpoints ),
    phi <- k[species] +b*dispersal+c*polygyny,
    vector[N_spp]:k ~ multi_normal( 0 , SIGMA ),
    matrix[N_spp,N_spp]: SIGMA <- cov_GPL2( Dmat, etasq, rhosq, 0.01),
    b ~ dnorm(0,1),
    c~dnorm(0,1),
    cutpoints ~ dnorm( 0 , 1.5 ),
    etasq~exponential(1),
    rhosq~exponential(1)
  ) , data=dat_list_continuous_phylogenetic , chains=4 , cores=4 ,cmdstan=T, messages=FALSE, refresh=0)


precis(m_strict_continuous_phylogenetic)
# mean   sd  5.5% 94.5% n_eff Rhat4
# b (dispersal)      0.74 0.24  0.37  1.13   533  1.00
# c (polygyny)    -1.53 0.49 -2.31 -0.78  1704  1.00
# etasq  3.41 1.77  1.13  6.57   167  1.01
# rhosq  1.95 1.38  0.37  4.48   175  1.02


dstan_strict <- combined[ complete.cases(combined$strictfdom,combined$sexbias_dispersal,combined$MatSysPMK,combined$corrected_species_id),]

dstan_strict$Polygyny<-dstan_strict$MatSysPMK=="POL"
dstan_strict[dstan_strict$sexbias_dispersal=="Both",]$sexbias_dispersal<-1
dstan_strict[dstan_strict$sexbias_dispersal=="Female",]$sexbias_dispersal<-1
dstan_strict[dstan_strict$sexbias_dispersal=="Male",]$sexbias_dispersal<-0

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
  dispersal = as.integer(dstan_strict$sexbias_dispersal),
  polygyny = as.integer(dstan_strict$Polygyny),
  species = as.integer(as.factor(dstan_strict$corrected_species_id))
)

dat_list_continuous_phylogenetic$Dmat<-Dmat[ unique(dstan_strict$corrected_species_id),unique(dstan_strict$corrected_species_id) ]/max(Dmat)

colnames(dat_list_continuous_phylogenetic$Dmat)<-as.integer(as.factor(colnames(dat_list_continuous_phylogenetic$Dmat)))
rownames(dat_list_continuous_phylogenetic$Dmat)<-as.integer(as.factor(rownames(dat_list_continuous_phylogenetic$Dmat)))

dat_list_continuous_phylogenetic$N_spp<-length(unique(dstan_strict$corrected_species_id))


m_strict_continuous_phylogenetic <- ulam(
  alist(
    R ~ dordlogit( phi , cutpoints ),
    phi <- k[species] +b*dispersal+c*polygyny,
    vector[N_spp]:k ~ multi_normal( 0 , SIGMA ),
    matrix[N_spp,N_spp]: SIGMA <- cov_GPL2( Dmat, etasq, rhosq, 0.01),
    b ~ dnorm(0,1),
    c~dnorm(0,1),
    cutpoints ~ dnorm( 0 , 1.5 ),
    etasq~exponential(1),
    rhosq~exponential(1)
  ) , data=dat_list_continuous_phylogenetic , chains=4 , cores=4 ,cmdstan=T, messages=FALSE, refresh=0)


precis(m_strict_continuous_phylogenetic)




dstan_strict <- combined[ complete.cases(combined$strictfdom,combined$sexbias_dispersal,combined$CanineDimorphism,combined$corrected_species_id),]

dstan_strict[dstan_strict$sexbias_dispersal=="Both",]$sexbias_dispersal<-1
dstan_strict[dstan_strict$sexbias_dispersal=="Female",]$sexbias_dispersal<-1
dstan_strict[dstan_strict$sexbias_dispersal=="Male",]$sexbias_dispersal<-0

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
  dispersal = as.integer(dstan_strict$sexbias_dispersal),
  dimorphism = standardize(dstan_strict$CanineDimorphism),
  species = as.integer(as.factor(dstan_strict$corrected_species_id))
)

dat_list_continuous_phylogenetic$Dmat<-Dmat[ unique(dstan_strict$corrected_species_id),unique(dstan_strict$corrected_species_id) ]/max(Dmat)

colnames(dat_list_continuous_phylogenetic$Dmat)<-as.integer(as.factor(colnames(dat_list_continuous_phylogenetic$Dmat)))
rownames(dat_list_continuous_phylogenetic$Dmat)<-as.integer(as.factor(rownames(dat_list_continuous_phylogenetic$Dmat)))

dat_list_continuous_phylogenetic$N_spp<-length(unique(dstan_strict$corrected_species_id))


m_strict_continuous_phylogenetic <- ulam(
  alist(
    R ~ dordlogit( phi , cutpoints ),
    phi <- k[species] +b*dispersal+c*dimorphism,
    vector[N_spp]:k ~ multi_normal( 0 , SIGMA ),
    matrix[N_spp,N_spp]: SIGMA <- cov_GPL2( Dmat, etasq, rhosq, 0.01),
    b ~ dnorm(0,1),
    c~dnorm(0,1),
    cutpoints ~ dnorm( 0 , 1.5 ),
    etasq~exponential(1),
    rhosq~exponential(1)
  ) , data=dat_list_continuous_phylogenetic , chains=4 , cores=4 ,cmdstan=T, messages=FALSE, refresh=0)


precis(m_strict_continuous_phylogenetic)
