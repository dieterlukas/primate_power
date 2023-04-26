# Multivariate analyses

# Number of females + Mating System

filtereddata<-select(combined,corrected_species_id,strictfdom,females,MatSysPMK)

filtereddata<-filtereddata[complete.cases(filtereddata),]


dat_list_strict <- list(
  R = as.integer(as.factor(filtereddata$strictfdom)),
  females = standardize(filtereddata$females),
  matingsystem = as.integer(as.factor(filtereddata$MatSysPMK)),
  species = as.integer(as.factor(filtereddata$corrected_species_id))
)

m_strict_continuous <- ulam(
  alist(
    R ~ dordlogit( phi , cutpoints ),
    phi <-a + b*females + c[matingsystem],
    a ~ dnorm( 0 , 5 ),
    b ~ dnorm(0,5),
    c[matingsystem] ~ dnorm(c_bar,sigma_c),
    c_bar ~ dnorm(0,1),
    sigma_c ~  dexp(1),
    cutpoints ~ dnorm( 0 , 5 )
  ) , data=dat_list_strict , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)

precis(m_strict_continuous)

posterior<-extract.samples(m_strict_continuous)
precis(posterior$c[,2]-posterior$c[,1])
precis(posterior$c[,3]-posterior$c[,1])
precis(posterior$c[,4]-posterior$c[,1])
precis(posterior$c[,3]-posterior$c[,2])
precis(posterior$c[,4]-posterior$c[,2])
precis(posterior$c[,4]-posterior$c[,3])


polygyny<- link( m_strict_continuous, data=data.frame( matingsystem=rep(1,9), females=seq(from=-1,to=1,by=0.25)))
polygyny_mean <- apply( polygyny , 2 , mean ) # mean predicted
polygyny_ci <- apply( polygyny , 2 , PI , prob=0.97 ) # 97 percentile compatibility




pdat <- data.frame(matingsystem=rep(1,9), females=seq(from=-1,to=1,by=0.25))
overallphi_speciesaverages<-matrix(ncol=1, nrow=9)
for(i in 1:nrow(pdat)){
  overallphi_speciesaverages[i,]<-precis(m_strict_continuous,depth=2)[1,1]+precis(m_strict_continuous,depth=2)[2,1]*pdat$females[i]+precis(m_strict_continuous)[3,1]
}

overallprobs_speciesaverages<-pordlogit( 1:3 , overallphi_speciesaverages , precis(m_strict_continuous,depth=2)[9:10,1] )



# Number of males + Mating System

filtereddata<-select(combined,corrected_species_id,strictfdom,males,MatSysPMK)

filtereddata<-filtereddata[complete.cases(filtereddata),]


dat_list_strict <- list(
  R = as.integer(as.factor(filtereddata$strictfdom)),
  males = standardize(filtereddata$males),
  matingsystem = as.integer(as.factor(filtereddata$MatSysPMK)),
  species = as.integer(as.factor(filtereddata$corrected_species_id))
)

m_strict_continuous <- ulam(
  alist(
    R ~ dordlogit( phi , cutpoints ),
    phi <-a + b*males + c[matingsystem],
    a ~ dnorm( 0 , 5 ),
    b ~ dnorm(0,5),
    c[matingsystem] ~ dnorm(c_bar,sigma_c),
    c_bar ~ dnorm(0,1),
    sigma_c ~  dexp(1),
    cutpoints ~ dnorm( 0 , 5 )
  ) , data=dat_list_strict , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)

precis(m_strict_continuous)

posterior<-extract.samples(m_strict_continuous)
precis(posterior$c[,2]-posterior$c[,1])
precis(posterior$c[,3]-posterior$c[,1])
precis(posterior$c[,4]-posterior$c[,1])
precis(posterior$c[,3]-posterior$c[,2])
precis(posterior$c[,4]-posterior$c[,2])
precis(posterior$c[,4]-posterior$c[,3])


# Canine dimorphism + Mating System

filtereddata<-select(combined,corrected_species_id,strictfdom,CanineDimorphism,MatSysPMK)

filtereddata<-filtereddata[complete.cases(filtereddata),]


dat_list_strict <- list(
  R = as.integer(as.factor(filtereddata$strictfdom)),
  CanineDimorphism = standardize(filtereddata$CanineDimorphism),
  matingsystem = as.integer(as.factor(filtereddata$MatSysPMK)),
  species = as.integer(as.factor(filtereddata$corrected_species_id))
)

m_strict_continuous <- ulam(
  alist(
    R ~ dordlogit( phi , cutpoints ),
    phi <-a + b*CanineDimorphism + c[matingsystem],
    a ~ dnorm( 0 , 1 ),
    b ~ dnorm(0,1),
    c[matingsystem] ~ dnorm(c_bar,sigma_c),
    c_bar ~ dnorm(0,1),
    sigma_c ~  dexp(1),
    cutpoints ~ dnorm( 0 , 1.5 )
  ) , data=dat_list_strict , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)

precis(m_strict_continuous)

posterior<-extract.samples(m_strict_continuous)
precis(posterior$c[,2]-posterior$c[,1])
precis(posterior$c[,3]-posterior$c[,1])
precis(posterior$c[,4]-posterior$c[,1])
precis(posterior$c[,3]-posterior$c[,2])
precis(posterior$c[,4]-posterior$c[,2])
precis(posterior$c[,4]-posterior$c[,3])


# Canine dimorphism + Body size dimorphism

filtereddata<-select(combined,corrected_species_id,strictfdom,SexualDimorphism_MaleWeight_over_FemaleWeight,CanineDimorphism)

filtereddata<-filtereddata[complete.cases(filtereddata),]


dat_list_strict <- list(
  R = as.integer(as.factor(filtereddata$strictfdom)),
  bodydimorphism = standardize(filtereddata$SexualDimorphism_MaleWeight_over_FemaleWeight),
  caninedimorphism = standardize(filtereddata$CanineDimorphism),
  species = as.integer(as.factor(filtereddata$corrected_species_id))
)

m_strict_continuous <- ulam(
  alist(
    R ~ dordlogit( phi , cutpoints ),
    phi <-a + b*bodydimorphism + c*caninedimorphism,
    a ~ dnorm( 0 , 1 ),
    b ~ dnorm(0,1),
    c ~ dnorm(0,1),
    cutpoints ~ dnorm( 0 , 1.5 )
  ) , data=dat_list_strict , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)

precis(m_strict_continuous)




# Sex ratio + Mating System

filtereddata<-select(combined,corrected_species_id,strictfdom,sexratio,MatSysPMK)

filtereddata<-filtereddata[complete.cases(filtereddata),]


dat_list_strict <- list(
  R = as.integer(as.factor(filtereddata$strictfdom)),
  sexratio = standardize(filtereddata$sexratio),
  matingsystem = as.integer(as.factor(filtereddata$MatSysPMK)),
  species = as.integer(as.factor(filtereddata$corrected_species_id))
)

m_strict_continuous <- ulam(
  alist(
    R ~ dordlogit( phi , cutpoints ),
    phi <-a + b*sexratio + c[matingsystem],
    a ~ dnorm( 0 , 5 ),
    b ~ dnorm(0,5),
    c[matingsystem] ~ dnorm(c_bar,sigma_c),
    c_bar ~ dnorm(0,1),
    sigma_c ~  dexp(1),
    cutpoints ~ dnorm( 0 , 5 )
  ) , data=dat_list_strict , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)

precis(m_strict_continuous)

posterior<-extract.samples(m_strict_continuous)
precis(posterior$c[,2]-posterior$c[,1])
precis(posterior$c[,3]-posterior$c[,1])
precis(posterior$c[,4]-posterior$c[,1])
precis(posterior$c[,3]-posterior$c[,2])
precis(posterior$c[,4]-posterior$c[,2])
precis(posterior$c[,4]-posterior$c[,3])



# Sex ratio nested in species and additionally canine size dimorphism across observations

filtereddata<-select(combined,corrected_species_id,strictfdom,sexratio,CanineDimorphism)
filtereddata<-filtereddata[complete.cases(filtereddata),]

dat_list_strict <- list(
  R = as.integer(as.factor(filtereddata$strictfdom)),
  sexratio = standardize(filtereddata$sexratio),
  species = as.integer(as.factor(filtereddata$corrected_species_id)),
  caninedimorphism = standardize(filtereddata$CanineDimorphism)
)

m_withinspecies_sexratio <- ulam(
  alist(
    R ~ dordlogit( phi , cutpoints ),
    phi <-a_species[species] + b_species[species]*sexratio +c*caninedimorphism,
    a_species[species]~normal( a_bar , a_sigma ),
    b_species[species]~dnorm(b_bar,b_sigma),
    a_bar~dnorm(0,5),
    b_bar~dnorm(0,5),
    a_sigma~dexp(1),
    b_sigma~dexp(1),
    c ~ dnorm(-3,2),
    cutpoints ~ dnorm( 0 , 3 )
  ) , data=dat_list_strict , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)

precis(m_withinspecies_sexratio)


m_withinspecies_sexratio <- ulam(
  alist(
    R ~ dordlogit( phi , cutpoints ),
    phi <-a_species[species] + b_species[species]*sexratio +c*caninedimorphism+d*sexratio,
    a_species[species]~normal( a_bar , a_sigma ),
    b_species[species]~dnorm(b_bar,b_sigma),
    a_bar~dnorm(0,5),
    b_bar~dnorm(0,5),
    a_sigma~dexp(1),
    b_sigma~dexp(1),
    c ~ dnorm(-3,2),
    d ~ dnorm(0,2),
    cutpoints ~ dnorm( 0 , 3 )
  ) , data=dat_list_strict , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)

precis(m_withinspecies_sexratio)



# Dimorphism + Mating System + Arboreality

filtereddata<-select(combined,corrected_species_id,strictfdom,SexualDimorphism_MaleWeight_over_FemaleWeight,MatSysPMK,Strata_Wilman)

filtereddata<-filtereddata[complete.cases(filtereddata),]


dat_list_strict <- list(
  R = as.integer(as.factor(filtereddata$strictfdom)),
  dimorphism = standardize(filtereddata$SexualDimorphism_MaleWeight_over_FemaleWeight),
  matingsystem = as.integer(as.factor(filtereddata$MatSysPMK)),
  arboreality = as.integer(as.factor(filtereddata$Strata_Wilman)),
  species = as.integer(as.factor(filtereddata$corrected_species_id))
)

m_strict_continuous <- ulam(
  alist(
    R ~ dordlogit( phi , cutpoints ),
    phi <-a + b*dimorphism + c[matingsystem] + d[arboreality],
    a ~ dnorm( 0 , 1 ),
    b ~ dnorm(0,1),
    c[matingsystem] ~ dnorm(c_bar,sigma_c),
    d[arboreality] ~ dnorm(d_bar,sigma_d),
    c_bar ~ dnorm(0,1),
    d_bar ~ dnorm(0,1),
    sigma_c ~  dexp(1),
    sigma_d ~  dexp(1),
    cutpoints ~ dnorm( 0 , 1.5 )
  ) , data=dat_list_strict , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)

precis(m_strict_continuous,depth=2)





filtereddata<-select(combined,corrected_species_id,strictfdom,SexualDimorphism_MaleWeight_over_FemaleWeight,body_mass)

filtereddata<-filtereddata[complete.cases(filtereddata),]


dat_list_strict <- list(
  R = as.integer(as.factor(filtereddata$strictfdom)),
  dimorphism = standardize(filtereddata$SexualDimorphism_MaleWeight_over_FemaleWeight),
  bodymass = standardize(filtereddata$body_mass),
  species = as.integer(as.factor(filtereddata$corrected_species_id))
)

m_strict_continuous <- ulam(
  alist(
    R ~ dordlogit( phi , cutpoints ),
    phi <-a + b*dimorphism + c*bodymass,
    a ~ dnorm( 0 , 1 ),
    b ~ dnorm(0,1),
    c ~ dnorm(0,1),
    cutpoints ~ dnorm( 0 , 1.5 )
  ) , data=dat_list_strict , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)

precis(m_strict_continuous)





# Female canine height in the different dominance systems

filtereddata<-select(combined,corrected_species_id,strictfdom,female_canine_height,female_mass)

filtereddata<-filtereddata[complete.cases(filtereddata),]


dat_list_strict <- list(
  dominance = as.integer(as.factor(filtereddata$strictfdom)),
  canine = (filtereddata$female_canine_height),
  bodymass = (filtereddata$female_mass)
)

m_strict_continuous <- ulam(
  alist(
    canine ~ dnorm(mu,sigma),
    mu <-a+b[dominance] + c*bodymass,
    a ~ dnorm( 0 , 1 ),
    b[dominance] ~ dnorm( 0 , 1 ),
    c ~ dnorm(0,1),
    sigma ~ dexp(1)
  ) , data=dat_list_strict , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)

precis(m_strict_continuous,depth=2)

post<-extract.samples(m_strict_continuous)
contrast_femaledominant_maledominant<-post$b[,3]-post$b[,1]
contrast_femaledominant_codominant<-post$b[,3]-post$b[,2]
contrast_maledominant_codominant<-post$b[,1]-post$b[,2]

precis(list(contrast_femaledominant_maledominant,contrast_femaledominant_codominant,contrast_maledominant_codominant))


dat_list_strict_phylogenetic <- list(
  dominance = as.integer(as.factor(filtereddata$strictfdom)),
  canine = standardize(filtereddata$female_canine_height),
  bodymass = standardize(filtereddata$female_mass),
  species = as.integer(as.factor(filtereddata$corrected_species_id))
)


spp_obs<-unique(filtereddata$corrected_species_id)
spp_obs<-matrix(nrow=length(spp_obs))
rownames(spp_obs)<-unique(filtereddata$corrected_species_id)

# We match the species in the dataset to the species in the phylogenetic tree
missing<-treedata(inputtree,data=spp_obs,warnings=FALSE)
# We remove species with no data from the tree
mtree<-missing$phy
Dmat<-cophenetic(mtree)

# We calculate the pair-wise phylogenetic distances among the species in the tree; this gives a matrix where values are the total branch length needed to connect two species (diagonal is the distance of a species to itself so zero)
dat_list_strict_phylogenetic$Dmat<-Dmat[ (filtereddata$corrected_species_id),(filtereddata$corrected_species_id) ]/max(Dmat)

colnames(dat_list_strict_phylogenetic$Dmat)<-as.integer(as.factor(colnames(dat_list_strict_phylogenetic$Dmat)))
rownames(dat_list_strict_phylogenetic$Dmat)<-as.integer(as.factor(rownames(dat_list_strict_phylogenetic$Dmat)))

dat_list_strict_phylogenetic$N_spp<-length(unique(filtereddata$corrected_species_id))

m_strict_continuous_phylogenetic <- ulam(
  alist(
    canine ~ multi_normal(mu,SIGMA),
    mu <-a + b[dominance] + c*bodymass,
    a ~ dnorm( 0 , 1 ),
    c ~ dnorm( 0 , 1 ),
    b[dominance] ~ dnorm(0,1),
    matrix[N_spp,N_spp]: SIGMA <- cov_GPL2( Dmat, etasq, rhosq, 0.01),
    etasq~half_normal(1,0.25),
    rhosq~half_normal(3,0.25)
  ) , data=dat_list_strict_phylogenetic , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)

precis(m_strict_continuous_phylogenetic,depth=2)






# Male canine height in the different dominance systems

filtereddata<-select(combined,corrected_species_id,strictfdom,male_canine_height,male_mass)

filtereddata<-filtereddata[complete.cases(filtereddata),]


dat_list_strict <- list(
  dominance = as.integer(as.factor(filtereddata$strictfdom)),
  canine = standardize(filtereddata$male_canine_height),
  bodymass = standardize(filtereddata$male_mass)
)

m_strict_continuous <- ulam(
  alist(
    canine ~ dnorm(mu,sigma),
    mu <-a+b[dominance] + c*bodymass,
    a ~ dnorm( 0 , 1 ),
    b[dominance] ~ dnorm( 0 , 1 ),
    c ~ dnorm(0,1),
    sigma ~ dexp(1)
  ) , data=dat_list_strict , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)

precis(m_strict_continuous,depth=2)

post<-extract.samples(m_strict_continuous)
contrast_femaledominant_maledominant<-post$b[,3]-post$b[,1]
contrast_femaledominant_codominant<-post$b[,3]-post$b[,2]
contrast_maledominant_codominant<-post$b[,1]-post$b[,2]

precis(list(contrast_femaledominant_maledominant,contrast_femaledominant_codominant,contrast_maledominant_codominant))


dat_list_strict_phylogenetic <- list(
  dominance = as.integer(as.factor(filtereddata$strictfdom)),
  canine = standardize(filtereddata$male_canine_height),
  bodymass = standardize(filtereddata$male_mass),
  species = as.integer(as.factor(filtereddata$corrected_species_id))
)


spp_obs<-unique(filtereddata$corrected_species_id)
spp_obs<-matrix(nrow=length(spp_obs))
rownames(spp_obs)<-unique(filtereddata$corrected_species_id)

# We match the species in the dataset to the species in the phylogenetic tree
missing<-treedata(inputtree,data=spp_obs,warnings=FALSE)
# We remove species with no data from the tree
mtree<-missing$phy
Dmat<-cophenetic(mtree)

# We calculate the pair-wise phylogenetic distances among the species in the tree; this gives a matrix where values are the total branch length needed to connect two species (diagonal is the distance of a species to itself so zero)
dat_list_strict_phylogenetic$Dmat<-Dmat[ (filtereddata$corrected_species_id),(filtereddata$corrected_species_id) ]/max(Dmat)

colnames(dat_list_strict_phylogenetic$Dmat)<-as.integer(as.factor(colnames(dat_list_strict_phylogenetic$Dmat)))
rownames(dat_list_strict_phylogenetic$Dmat)<-as.integer(as.factor(rownames(dat_list_strict_phylogenetic$Dmat)))

dat_list_strict_phylogenetic$N_spp<-length(unique(filtereddata$corrected_species_id))

m_strict_continuous_phylogenetic <- ulam(
  alist(
    canine ~ multi_normal(mu,SIGMA),
    mu <-a + b[dominance] + c*bodymass,
    a ~ dnorm( 0 , 1 ),
    c ~ dnorm( 0 , 1 ),
    b[dominance] ~ dnorm(0,1),
    matrix[N_spp,N_spp]: SIGMA <- cov_GPL2( Dmat, etasq, rhosq, 0.01),
    etasq~half_normal(1,0.25),
    rhosq~half_normal(3,0.25)
  ) , data=dat_list_strict_phylogenetic , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)

precis(m_strict_continuous_phylogenetic,depth=2)




# canine dimorphism in the different dominance systems

filtereddata<-select(combined,corrected_species_id,strictfdom,CanineDimorphism,male_mass)

filtereddata<-filtereddata[complete.cases(filtereddata),]


dat_list_strict <- list(
  dominance = as.integer(as.factor(filtereddata$strictfdom)),
  canine = standardize(filtereddata$CanineDimorphism),
  bodymass = standardize(filtereddata$male_mass)
)

m_strict_continuous <- ulam(
  alist(
    canine ~ dnorm(mu,sigma),
    mu <-a+b[dominance] + c*bodymass,
    a ~ dnorm( 0 , 1 ),
    b[dominance] ~ dnorm( 0 , 1 ),
    c ~ dnorm(0,1),
    sigma ~ dexp(1)
  ) , data=dat_list_strict , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)

precis(m_strict_continuous,depth=2)

post<-extract.samples(m_strict_continuous)
contrast_femaledominant_maledominant<-post$b[,3]-post$b[,1]
contrast_femaledominant_codominant<-post$b[,3]-post$b[,2]
contrast_maledominant_codominant<-post$b[,1]-post$b[,2]

precis(list(contrast_femaledominant_maledominant,contrast_femaledominant_codominant,contrast_maledominant_codominant))


dat_list_strict_phylogenetic <- list(
  dominance = as.integer(as.factor(filtereddata$strictfdom)),
  canine = standardize(filtereddata$CanineDimorphism),
  bodymass = standardize(filtereddata$male_mass),
  species = as.integer(as.factor(filtereddata$corrected_species_id))
)


spp_obs<-unique(filtereddata$corrected_species_id)
spp_obs<-matrix(nrow=length(spp_obs))
rownames(spp_obs)<-unique(filtereddata$corrected_species_id)

# We match the species in the dataset to the species in the phylogenetic tree
missing<-treedata(inputtree,data=spp_obs,warnings=FALSE)
# We remove species with no data from the tree
mtree<-missing$phy
Dmat<-cophenetic(mtree)

# We calculate the pair-wise phylogenetic distances among the species in the tree; this gives a matrix where values are the total branch length needed to connect two species (diagonal is the distance of a species to itself so zero)
dat_list_strict_phylogenetic$Dmat<-Dmat[ (filtereddata$corrected_species_id),(filtereddata$corrected_species_id) ]/max(Dmat)

colnames(dat_list_strict_phylogenetic$Dmat)<-as.integer(as.factor(colnames(dat_list_strict_phylogenetic$Dmat)))
rownames(dat_list_strict_phylogenetic$Dmat)<-as.integer(as.factor(rownames(dat_list_strict_phylogenetic$Dmat)))

dat_list_strict_phylogenetic$N_spp<-length(unique(filtereddata$corrected_species_id))

m_strict_continuous_phylogenetic <- ulam(
  alist(
    canine ~ multi_normal(mu,SIGMA),
    mu <-a + b[dominance] + c*bodymass,
    a ~ dnorm( 0 , 1 ),
    c ~ dnorm( 0 , 1 ),
    b[dominance] ~ dnorm(0,1),
    matrix[N_spp,N_spp]: SIGMA <- cov_GPL2( Dmat, etasq, rhosq, 0.01),
    etasq~half_normal(1,0.25),
    rhosq~half_normal(3,0.25)
  ) , data=dat_list_strict_phylogenetic , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)

precis(m_strict_continuous_phylogenetic,depth=2)



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






# seasonality + mating system

dstan_continuous <- combined[ complete.cases(combined$perc_won_females,combined$r_seasonality_value,combined$MatSysPMK,combined$corrected_species_id),]

# Continuous outcome: percentage of fights won by females
dat_list_continuous <- list(
  N_spp = nrow(dstan_continuous),
  perc_won_females = as.integer(dstan_continuous$perc_won_females),
  r_seasonality_value = standardize(dstan_continuous$r_seasonality_value),
  MatSysPMK = as.integer(as.factor(dstan_continuous$MatSysPMK)),
  species = as.integer(as.factor(dstan_continuous$corrected_species_id)),
  total = rep(100,nrow(dstan_continuous))
)

m_continuous_continuous <- ulam(
  alist(
    perc_won_females ~ dbinom(total,p),
    logit(p) <- a +b*r_seasonality_value+c*MatSysPMK,
    ## adaptive priors
    a ~ dnorm( 0 , 1 ),
    b ~ dnorm(0,1),
    c ~ dnorm(0,1)
  ) , data=dat_list_continuous , chains=4 , cores=4 , log_lik=TRUE , cmdstan=T)

precis(m_continuous_continuous,depth=2)




# environmental harshness + mating system

dstan_continuous <- combined[ complete.cases(combined$perc_won_females,combined$env_harshness,combined$MatSysPMK,combined$corrected_species_id),]

# Continuous outcome: percentage of fights won by females
dat_list_continuous <- list(
  N_spp = nrow(dstan_continuous),
  perc_won_females = as.integer(dstan_continuous$perc_won_females),
  env_harshness = standardize(dstan_continuous$env_harshness),
  MatSysPMK = as.integer(as.factor(dstan_continuous$MatSysPMK)),
  species = as.integer(as.factor(dstan_continuous$corrected_species_id)),
  total = rep(100,nrow(dstan_continuous))
)

m_continuous_continuous <- ulam(
  alist(
    perc_won_females ~ dbinom(total,p),
    logit(p) <- a +b*env_harshness+c*MatSysPMK,
    ## adaptive priors
    a ~ dnorm( 0 , 1 ),
    b ~ dnorm(0,1),
    c ~ dnorm(0,1)
  ) , data=dat_list_continuous , chains=4 , cores=4 , log_lik=TRUE , cmdstan=T)

precis(m_continuous_continuous,depth=2)





# ssd + group sex ratio

dstan_continuous <- combined[ complete.cases(combined$perc_won_females,combined$SexualDimorphism_MaleWeight_over_FemaleWeight,combined$sexratio,combined$corrected_species_id),]

# Continuous outcome: percentage of fights won by females
dat_list_continuous <- list(
  N_spp = nrow(dstan_continuous),
  perc_won_females = as.integer(dstan_continuous$perc_won_females),
  SexualDimorphism_MaleWeight_over_FemaleWeight = standardize(dstan_continuous$SexualDimorphism_MaleWeight_over_FemaleWeight),
  sexratio = standardize(dstan_continuous$sexratio),
  species = as.integer(as.factor(dstan_continuous$corrected_species_id)),
  total = rep(100,nrow(dstan_continuous))
)

m_continuous_continuous <- ulam(
  alist(
    perc_won_females ~ dbinom(total,p),
    logit(p) <- a +b*SexualDimorphism_MaleWeight_over_FemaleWeight+c*sexratio,
    ## adaptive priors
    a ~ dnorm( 0 , 1 ),
    b ~ dnorm(0,1),
    c ~ dnorm(0,1)
  ) , data=dat_list_continuous , chains=4 , cores=4 , log_lik=TRUE , cmdstan=T)

precis(m_continuous_continuous,depth=2)



# philopatry + relatedness

dstan_continuous <- combined[ complete.cases(combined$perc_won_females,combined$female_average_relatedness,combined$female_dispersal,combined$corrected_species_id),]

# Continuous outcome: percentage of fights won by females
dat_list_continuous <- list(
  N_spp = nrow(dstan_continuous),
  perc_won_females = as.integer(dstan_continuous$perc_won_females),
  female_average_relatedness = standardize(dstan_continuous$female_average_relatedness),
  female_dispersal = as.integer(as.factor(dstan_continuous$female_dispersal)),
  species = as.integer(as.factor(dstan_continuous$corrected_species_id)),
  total = rep(100,nrow(dstan_continuous))
)

m_continuous_continuous <- ulam(
  alist(
    perc_won_females ~ dbinom(total,p),
    logit(p) <- a +b*female_average_relatedness+c*female_dispersal,
    ## adaptive priors
    a ~ dnorm( 0 , 1 ),
    b ~ dnorm(0,1),
    c ~ dnorm(0,1)
  ) , data=dat_list_continuous , chains=4 , cores=4 , log_lik=TRUE , cmdstan=T)

precis(m_continuous_continuous,depth=2)

plot(dstan_continuous$perc_won_females~dstan_continuous$female_average_relatedness,col=as.integer(as.factor(dstan_continuous$female_dispersal)))


plot(combined$perc_won_females~combined$female_average_relatedness,col=as.integer(as.factor(combined$cooperative_breeder)),pch=as.integer(as.factor(combined$cooperative_breeder)),cex=2)

plot(combined$perc_won_females~combined$female_average_relatedness,cex=combined$females/2)




# malemale aggression + sex ratio

dstan_continuous <- combined[ complete.cases(combined$perc_won_females,combined$perc_aggression_mm,combined$sexratio,combined$corrected_species_id),]

# Continuous outcome: percentage of fights won by females
dat_list_continuous <- list(
  N_spp = nrow(dstan_continuous),
  perc_won_females = as.integer(dstan_continuous$perc_won_females),
  perc_aggression_mm = standardize(dstan_continuous$perc_aggression_mm),
  sexratio = standardize((dstan_continuous$sexratio)),
  species = as.integer(as.factor(dstan_continuous$corrected_species_id)),
  total = rep(100,nrow(dstan_continuous))
)

m_continuous_continuous <- ulam(
  alist(
    perc_won_females ~ dbinom(total,p),
    logit(p) <- a +b*perc_aggression_mm+c*sexratio,
    ## adaptive priors
    a ~ dnorm( 0 , 1 ),
    b ~ dnorm(0,1),
    c ~ dnorm(0,1)
  ) , data=dat_list_continuous , chains=4 , cores=4 , log_lik=TRUE , cmdstan=T)

precis(m_continuous_continuous,depth=2)




# polygyny + female dispersal

dstan_continuous <- combined[ complete.cases(combined$perc_won_females,combined$MatSysPMK,combined$sexbias_dispersal,combined$corrected_species_id),]

# Continuous outcome: percentage of fights won by females
dat_list_continuous <- list(
  N_spp = nrow(dstan_continuous),
  perc_won_females = as.integer(dstan_continuous$perc_won_females),
  MatSysPMK = as.integer(as.factor(dstan_continuous$MatSysPMK)),
  sexbias_dispersal = as.integer(as.factor(dstan_continuous$sexbias_dispersal)),
  species = as.integer(as.factor(dstan_continuous$corrected_species_id)),
  total = rep(100,nrow(dstan_continuous))
)

m_continuous_continuous <- ulam(
  alist(
    perc_won_females ~ dbinom(total,p),
    logit(p) <- a +b[MatSysPMK]+c[sexbias_dispersal],
    ## adaptive priors
    a ~ dnorm( 0 , 1 ),
    b[MatSysPMK] ~ dnorm(0,1),
    c[sexbias_dispersal] ~ dnorm(0,1)
  ) , data=dat_list_continuous , chains=4 , cores=4 , log_lik=TRUE , cmdstan=T)

precis(m_continuous_continuous,depth=2)


posterior<-extract.samples(m_continuous_continuous)
precis(posterior$b[,2]-posterior$b[,1])
precis(posterior$b[,3]-posterior$b[,1])
precis(posterior$b[,4]-posterior$b[,1])
precis(posterior$b[,3]-posterior$b[,2])
precis(posterior$b[,4]-posterior$b[,2])
precis(posterior$b[,4]-posterior$b[,3])

precis(posterior$c[,2]-posterior$c[,1])
precis(posterior$c[,3]-posterior$c[,1])
precis(posterior$c[,3]-posterior$c[,2])




# environmental seasonality versus reproductive synchrony
dstan_continuous <- combined[ complete.cases(combined$perc_won_females,combined$NDVI_Seasonality,combined$r_seasonality_value,combined$corrected_species_id),]

# Continuous outcome: percentage of fights won by females
dat_list_continuous <- list(
  N_spp = nrow(dstan_continuous),
  perc_won_females = as.integer(dstan_continuous$perc_won_females),
  NDVI_Seasonality = standardize(dstan_continuous$NDVI_Seasonality),
  r_seasonality_value = standardize(dstan_continuous$r_seasonality_value),
  species = as.integer(as.factor(dstan_continuous$corrected_species_id)),
  total = rep(100,nrow(dstan_continuous))
)

m_continuous_continuous <- ulam(
  alist(
    perc_won_females ~ dbinom(total,p),
    logit(p) <- a +b*NDVI_Seasonality+c*r_seasonality_value,
    ## adaptive priors
    a ~ dnorm( 0 , 1 ),
    b ~ dnorm(0,1),
    c ~ dnorm(0,1)
  ) , data=dat_list_continuous , chains=4 , cores=4 , log_lik=TRUE , cmdstan=T)

precis(m_continuous_continuous,depth=2)
