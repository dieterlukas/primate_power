# Setting up the functions for the analyses

run_analyses_continuouspredictor <- function(continuouspredictor){
  
  
  colnames(combined)[colnames(combined) %in% continuouspredictor]<-c("continuouspredictor")
  
  results<-matrix(nrow=8,ncol=5)    
  results<-as.data.frame(results)
  colnames(results)<-c("outcome", "continuouspredictor","phylogeny", "estimate lower","estimate upper")
  results$outcome<-c(rep("perc_won", 2),rep("strict three way", 2),rep("strict female dominance",2),rep("strict male dominance",2))
  results$continuouspredictor<-rep(continuouspredictor,8)
  results$phylogeny<-c(rep(c("No","Yes"),4))
  
  dstan_continuous <- combined[ complete.cases(combined$perc_won_females,combined$continuouspredictor,combined$corrected_species_id),]
  
  dstan_strict <- combined[ complete.cases(combined$strictfdom,combined$continuouspredictor,combined$corrected_species_id),]
  
  dstan_relaxed <- combined[ complete.cases(combined$mostlyfdom,combined$continuouspredictor,combined$corrected_species_id),]  
  
  print("Finished setup")
  
  # Continuous outcome: percentage of fights won by females
  dat_list_continuous <- list(
    N_spp = nrow(dstan_continuous),
    perc_won_females = as.integer(dstan_continuous$perc_won_females),
    continuouspredictor = standardize(dstan_continuous$continuouspredictor),
    species = as.integer(as.factor(dstan_continuous$corrected_species_id)),
    total = rep(100,nrow(dstan_continuous))
  )
  
  m_continuous_continuous <- ulam(
    alist(
      perc_won_females ~ dbinom(total,p),
      logit(p) <- a +b*continuouspredictor,
      ## adaptive priors
      a ~ dnorm( 0 , 1 ),
      b ~ dnorm(0,1)
    ) , data=dat_list_continuous , chains=4 , cores=4 , log_lik=TRUE , cmdstan=T, messages=FALSE, refresh=0)
  
  results[1,4:5]<-precis(m_continuous_continuous)[2,3:4]
  
  
  print("finished m_continuous_continuous")
  
  
  
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
    perc_won_females = as.integer(dstan_continuous$perc_won_females),
    continuouspredictor = standardize(dstan_continuous$continuouspredictor),
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
      logit(p) <- k[species] +b*continuouspredictor,
      vector[N_spp]:k ~ multi_normal( 0 , SIGMA ),
      matrix[N_spp,N_spp]: SIGMA <- cov_GPL2( Dmat, etasq, rhosq, 0.01),
      ## adaptive priors
      b ~ dnorm(0,1),
      etasq~exponential(1),
      rhosq~exponential(1)
    ) , data=dat_list_continuous_phylogenetic , chains=4 , cores=4 , log_lik=TRUE , cmdstan=T, messages=FALSE, refresh=0)
  
  results[2,4:5]<-precis(m_continuous_continuous_phylogenetic)[1,3:4]
  
  print("finished m_continuous_continuous_phylogenetic")
  
  
  # Strict categorical classification of intersexual dominance
  
  dat_list_strict <- list(
    R = as.integer(as.factor(dstan_strict$strictfdom)),
    continuouspredictor = standardize(dstan_strict$continuouspredictor),
    species = as.integer(as.factor(dstan_strict$corrected_species_id))
  )
  
  m_strict_continuous <- ulam(
    alist(
      R ~ dordlogit( phi , cutpoints ),
      phi <-a + b*continuouspredictor,
      a ~ dnorm( 0 , 1 ),
      b ~ dnorm(0,1),
      cutpoints ~ dnorm( 0 , 1.5 )
    ) , data=dat_list_strict , chains=4 , cores=4 ,cmdstan=T, messages=FALSE, refresh=0)
  
  results[3,4:5]<-precis(m_strict_continuous)[2,3:4]
  
  # Strict categorical classification of intersexual dominance phylogenetic
  
  
  dat_list_strict_phylogenetic <- list(
    R = as.integer(as.factor(dstan_strict$strictfdom)),
    continuouspredictor = standardize(dstan_strict$continuouspredictor),
    species = as.integer(as.factor(dstan_strict$corrected_species_id))
  )
  
  
  spp_obs<-unique(dstan_strict$corrected_species_id)
  spp_obs<-matrix(nrow=length(spp_obs))
  rownames(spp_obs)<-unique(dstan_strict$corrected_species_id)
  
  # We match the species in the dataset to the species in the phylogenetic tree
  missing<-treedata(inputtree,data=spp_obs,warnings=FALSE)
  # We remove species with no data from the tree
  mtree<-missing$phy
  Dmat<-cophenetic(mtree)
  
  # We calculate the pair-wise phylogenetic distances among the species in the tree; this gives a matrix where values are the total branch length needed to connect two species (diagonal is the distance of a species to itself so zero)
  dat_list_strict_phylogenetic$Dmat<-Dmat[ unique(dstan_strict$corrected_species_id),unique(dstan_strict$corrected_species_id) ]/max(Dmat)
  
  colnames(dat_list_strict_phylogenetic$Dmat)<-as.integer(as.factor(colnames(dat_list_strict_phylogenetic$Dmat)))
  rownames(dat_list_strict_phylogenetic$Dmat)<-as.integer(as.factor(rownames(dat_list_strict_phylogenetic$Dmat)))
  
  dat_list_strict_phylogenetic$N_spp<-length(unique(dstan_strict$corrected_species_id))
  
  
  m_strict_continuous_phylogenetic <- ulam(
    alist(
      R ~ dordlogit( phi , cutpoints ),
      phi <-k[species] + b*continuouspredictor,
      vector[N_spp]:k ~ multi_normal( 0 , SIGMA ),
      matrix[N_spp,N_spp]: SIGMA <- cov_GPL2( Dmat, etasq, rhosq, 0.01),
      b ~ dnorm(0,1),
      cutpoints ~ dnorm( 0 , 1.5 ),
      etasq~exponential(1),
      rhosq~exponential(1)
    ) , data=dat_list_strict_phylogenetic , chains=4 , cores=4 ,cmdstan=T, messages=FALSE, refresh=0)
  
  results[4,4:5]<-precis(m_strict_continuous_phylogenetic)[1,3:4]
  
  print("finished m_strict_continuous_phylogenetic")
  
  
  
  # Female dominance classification of intersexual dominance
  dat_list_relaxed <- list(
    R = as.integer(as.factor(dstan_strict$strictfdom %in% 3))-1,
    continuouspredictor = standardize(dstan_strict$continuouspredictor),
    species = as.integer(as.factor(dstan_strict$corrected_species_id))
  )
  
  
  m_relaxed_continuouspredictor <- ulam(
    alist(
      R ~ dbinom(1,p),
      logit(p) <- a + b*continuouspredictor ,
      a ~ dnorm( 0 , 1 ),
      b ~ dnorm( 0 , 0.5 ),
      cutpoints ~ dnorm( 0 , 1.5 )
    ) , data=dat_list_relaxed , chains=4 , cores=4 , cmdstan = T, messages=FALSE )
  
  results[5,4:5]<-precis(m_relaxed_continuouspredictor)[2,3:4]
  
  print("finished m_femaledominance_continuouspredictor")
  
  # Female dominance classification of intersexual dominance phylogenetic
  
  spp_obs<-unique(dstan_strict$corrected_species_id)
  spp_obs<-matrix(nrow=length(spp_obs))
  rownames(spp_obs)<-unique(dstan_strict$corrected_species_id)
  
  # We match the species in the dataset to the species in the phylogenetic tree
  missing<-treedata(inputtree,data=spp_obs,warnings=FALSE)
  # We remove species with no data from the tree
  mtree<-missing$phy
  
  Dmat<-cophenetic(mtree)
  
  dat_list_relaxed_phylogenetic <- list(
    N_spp = nrow(dstan_strict),
    R = as.integer(as.factor(dstan_strict$strictfdom %in% 3))-1,
    continuouspredictor = standardize(dstan_strict$continuouspredictor),
    species = as.integer(as.factor(dstan_strict$corrected_species_id))
  )
  
  dat_list_relaxed_phylogenetic$Dmat<-Dmat[ unique(dstan_strict$corrected_species_id),unique(dstan_strict$corrected_species_id) ]/max(Dmat)
  
  colnames(dat_list_relaxed_phylogenetic$Dmat)<-as.integer(as.factor(colnames(dat_list_relaxed_phylogenetic$Dmat)))
  rownames(dat_list_relaxed_phylogenetic$Dmat)<-as.integer(as.factor(rownames(dat_list_relaxed_phylogenetic$Dmat)))
  
  dat_list_relaxed_phylogenetic$N_spp<-length(unique(dstan_strict$corrected_species_id))
  
  
  
  m_relaxed_continuous_phylogenetic <- ulam(
    alist(
      R ~ dbinom(1,p),
      logit(p) <-k[species] + b*continuouspredictor,
      vector[N_spp]:k ~ multi_normal( 0 , SIGMA ),
      matrix[N_spp,N_spp]: SIGMA <- cov_GPL2( Dmat, etasq, rhosq, 0.01),
      b ~ dnorm(0,1),
      cutpoints ~ dnorm( 0 , 1.5 ),
      etasq~exponential(1),
      rhosq~exponential(1)
    ) , data=dat_list_relaxed_phylogenetic , chains=4 , cores=4 ,cmdstan=T, messages=FALSE, refresh=0)
  
  results[6,4:5]<-precis(m_relaxed_continuous_phylogenetic)[1,3:4]
  
  print("finished m_femaledominance_continuous_phylogenetic")
  
  
  # Female dominance classification of intersexual dominance
  dat_list_relaxed <- list(
    R = as.integer(as.factor(dstan_strict$strictfdom %in% 1))-1,
    continuouspredictor = standardize(dstan_strict$continuouspredictor),
    species = as.integer(as.factor(dstan_strict$corrected_species_id))
  )
  
  
  m_relaxed_continuouspredictor <- ulam(
    alist(
      R ~ dbinom(1,p),
      logit(p) <- a + b*continuouspredictor ,
      a ~ dnorm( 0 , 1 ),
      b ~ dnorm( 0 , 0.5 ),
      cutpoints ~ dnorm( 0 , 1.5 )
    ) , data=dat_list_relaxed , chains=4 , cores=4 , cmdstan = T, messages=FALSE )
  
  results[7,4:5]<-precis(m_relaxed_continuouspredictor)[2,3:4]
  
  print("finished m_femaledominance_continuouspredictor")
  
  
  # Male dominance classification of intersexual dominance phylogenetic
  
  spp_obs<-unique(dstan_strict$corrected_species_id)
  spp_obs<-matrix(nrow=length(spp_obs))
  rownames(spp_obs)<-unique(dstan_strict$corrected_species_id)
  
  # We match the species in the dataset to the species in the phylogenetic tree
  missing<-treedata(inputtree,data=spp_obs,warnings=FALSE)
  # We remove species with no data from the tree
  mtree<-missing$phy
  
  Dmat<-cophenetic(mtree)
  
  dat_list_relaxed_phylogenetic <- list(
    N_spp = nrow(dstan_strict),
    R = as.integer(as.factor(dstan_strict$strictfdom %in% 1))-1,
    continuouspredictor = standardize(dstan_strict$continuouspredictor),
    species = as.integer(as.factor(dstan_strict$corrected_species_id))
  )
  
  dat_list_relaxed_phylogenetic$Dmat<-Dmat[ unique(dstan_strict$corrected_species_id),unique(dstan_strict$corrected_species_id) ]/max(Dmat)
  
  colnames(dat_list_relaxed_phylogenetic$Dmat)<-as.integer(as.factor(colnames(dat_list_relaxed_phylogenetic$Dmat)))
  rownames(dat_list_relaxed_phylogenetic$Dmat)<-as.integer(as.factor(rownames(dat_list_relaxed_phylogenetic$Dmat)))
  
  dat_list_relaxed_phylogenetic$N_spp<-length(unique(dstan_strict$corrected_species_id))
  
  
  
  m_relaxed_continuous_phylogenetic <- ulam(
    alist(
      R ~ dbinom( 1 , p ),
      logit(p) <-k[species] + b*continuouspredictor,
      vector[N_spp]:k ~ multi_normal( 0 , SIGMA ),
      matrix[N_spp,N_spp]: SIGMA <- cov_GPL2( Dmat, etasq, rhosq, 0.01),
      b ~ dnorm(0,1),
      cutpoints ~ dnorm( 0 , 1.5 ),
      etasq~exponential(1),
      rhosq~exponential(1)
    ) , data=dat_list_relaxed_phylogenetic , chains=4 , cores=4 ,cmdstan=T, messages=FALSE, refresh=0)
  
  results[8,4:5]<-precis(m_relaxed_continuous_phylogenetic)[1,3:4]
  
  print("finished m_femaledominance_continuous_phylogenetic")
  
  
  
  
  
  colnames(combined)[colnames(combined) %in% "continuouspredictor"]<-continuouspredictor
  
  return(results)
  
}






run_analyses_discretepredictor <- function(discretepredictor){
  
  colnames(combined)[colnames(combined) %in% discretepredictor]<-c("discretepredictor")
  
  
  dstan_continuous <- combined[ complete.cases(combined$perc_won_females,combined$discretepredictor,combined$corrected_species_id),]
  
  dstan_strict <- combined[ complete.cases(combined$strictfdom,combined$discretepredictor,combined$corrected_species_id),]
  
  dstan_relaxed <- combined[ complete.cases(combined$mostlyfdom,combined$discretepredictor,combined$corrected_species_id),]  
  
  
  totalvariants<-length(unique(dstan_continuous$discretepredictor))
  numberofcontrasts<-totalvariants*(totalvariants-1)/2
  
  results<-matrix(nrow=numberofcontrasts*4*2,ncol=5)    
  results<-as.data.frame(results)
  colnames(results)<-c("outcome", "discretepredictor","phylogeny", "contrast lower","contrast upper")
  results$outcome<-c(rep("perc_won", numberofcontrasts*2),rep("strict three way", numberofcontrasts*2),rep("strict female dominance",numberofcontrasts*2),rep("strict male dominance",numberofcontrasts*2))
  allcombinations<-NA
  for(n in 1:numberofcontrasts){
    allcombinations[n]<-paste(discretepredictor,combn(sort(unique(dstan_strict$discretepredictor)),2)[2,n],"minus",combn(sort(unique(dstan_strict$discretepredictor)),2)[1,n])
  }
  results$discretepredictor<-rep(allcombinations,8)
  count<-1
  
  
  
  # Continuous outcome: percentage of fights won by females
  dat_list_continuous <- list(
    perc_won_females = as.integer(dstan_continuous$perc_won_females),
    discretepredictor = as.integer(as.factor(dstan_continuous$discretepredictor)),
    species = as.integer(as.factor(dstan_continuous$corrected_species_id)),
    total = rep(100,nrow(dstan_continuous))
  )
  m_continuous_discrete <- ulam(
    alist(
      perc_won_females ~ dbinom(total,p),
      logit(p) <- a +b[discretepredictor],
      ## adaptive priors
      a ~ dnorm( 0 , 1 ),
      b[discretepredictor] ~ dnorm(0,1)
    ) , data=dat_list_continuous , chains=4 , cores=4 , log_lik=TRUE , cmdstan=T, messages=FALSE, refresh=0)
  
  posterior<-extract.samples(m_continuous_discrete)
  
  for(variants in 1:numberofcontrasts){
    
    logitthiscontrast<-inv_logit(posterior$b[,combn(c(1:totalvariants),2)[2,variants]])-inv_logit(posterior$b[,combn(c(1:totalvariants),2)[1,variants]])
    results[count,4]<-precis(as.numeric(logitthiscontrast))[1,3]
    results[count,5]<-precis(as.numeric(logitthiscontrast))[1,4]
    results[count,3]<-"No"
    count<-count+1
    
  }
  
  print("finished m_continuous_discrete")
  
  
  
  
  spp_obs<-unique(dstan_continuous$corrected_species_id)
  spp_obs<-matrix(nrow=length(spp_obs))
  rownames(spp_obs)<-unique(dstan_continuous$corrected_species_id)
  
  # We match the species in the dataset to the species in the phylogenetic tree
  missing<-treedata(inputtree,data=spp_obs,warnings=FALSE)
  # We remove species with no data from the tree
  mtree<-missing$phy
  # We remove species not included in the tree from the dataset
  speciesnames<-mtree$tip.label
  mdata<-dstan_continuous[dstan_continuous$corrected_species_id %in% speciesnames,]
  
  # We calculate the pair-wise phylogenetic distances among the species in the tree; this gives a matrix where values are the total branch length needed to connect two species (diagonal is the distance of a species to itself so zero)
  Dmat<-cophenetic(mtree)
  
  
  # Continuous outcome phylogenetic: percentage of fights won by females
  dat_list_continuous_discrete_phylogenetic <- list(
    perc_won_females = as.integer(dstan_continuous$perc_won_females),
    discretepredictor = as.integer(as.factor(dstan_continuous$discretepredictor)),
    species = as.integer(as.factor(dstan_continuous$corrected_species_id)),
    total = rep(100,nrow(dstan_continuous))
  )
  
  dat_list_continuous_discrete_phylogenetic$Dmat<-Dmat[ unique(dstan_continuous$corrected_species_id),unique(dstan_continuous$corrected_species_id) ]/max(Dmat)
  
  colnames(dat_list_continuous_discrete_phylogenetic$Dmat)<-as.integer(as.factor(colnames(dat_list_continuous_discrete_phylogenetic$Dmat)))
  rownames(dat_list_continuous_discrete_phylogenetic$Dmat)<-as.integer(as.factor(rownames(dat_list_continuous_discrete_phylogenetic$Dmat)))
  
  dat_list_continuous_discrete_phylogenetic$N_spp<-length(unique(dstan_continuous$corrected_species_id))
  
  m_continuous_discrete_phylogenetic <- ulam(
    alist(
      perc_won_females ~ dbinom(total,p),
      logit(p) <- k[species] +b[discretepredictor],
      vector[N_spp]:k ~ multi_normal( 0 , SIGMA ),
      matrix[N_spp,N_spp]: SIGMA <- cov_GPL2( Dmat, etasq, rhosq, 0.01),
      ## adaptive priors
      a ~ dnorm( 0 , 1 ),
      b[discretepredictor] ~ dnorm(0,1),
      etasq~exponential(1),
      rhosq~exponential(1)
    ) , data=dat_list_continuous_discrete_phylogenetic , chains=4 , cores=4 , log_lik=TRUE , cmdstan=T, messages=FALSE, refresh=0)
  
  posterior<-extract.samples(m_continuous_discrete_phylogenetic)
  
  for(variants in 1:numberofcontrasts){
    
    logitthiscontrast<-inv_logit(posterior$b[,combn(c(1:totalvariants),2)[2,variants]])-inv_logit(posterior$b[,combn(c(1:totalvariants),2)[1,variants]])
    results[count,4]<-precis(as.numeric(logitthiscontrast))[1,3]
    results[count,5]<-precis(as.numeric(logitthiscontrast))[1,4]
    results[count,3]<-"Yes"
    count<-count+1
    
  }
  
  print("finished m_continuous_discrete_phylogenetic")
  
  # Strict categorical classification of intersexual dominance
  
  dat_list <- list(
    R = as.integer(as.factor(dstan_strict$strictfdom)),
    origin = as.integer(as.factor(dstan_strict$discretepredictor)),
    species = as.integer(as.factor(dstan_strict$corrected_species_id))
  )
  
  m_strict_discrete <- ulam(
    alist(
      R ~ dordlogit( phi , cutpoints ),
      phi <- b[origin] ,
      b[origin] ~ dnorm( 0 , 0.5 ),
      cutpoints ~ dnorm( 0 , 1.5 )
    ) , data=dat_list , chains=4 , cores=4 ,cmdstan=T, messages=FALSE, refresh=0)
  
  posterior<-extract.samples(m_strict_discrete)
  
  for(variants in 1:numberofcontrasts){
    
    logitthiscontrast<-inv_logit(posterior$b[,combn(c(1:totalvariants),2)[2,variants]])-inv_logit(posterior$b[,combn(c(1:totalvariants),2)[1,variants]])
    results[count,4]<-precis(as.numeric(logitthiscontrast))[1,3]
    results[count,5]<-precis(as.numeric(logitthiscontrast))[1,4]
    results[count,3]<-"No"
    count<-count+1
    
  }
  
  print("finished m_strict_discrete")
  
  # Strict categorical classification of intersexual dominance phylogenetic
  spp_obs<-unique(dstan_strict$corrected_species_id)
  spp_obs<-matrix(nrow=length(spp_obs))
  rownames(spp_obs)<-unique(dstan_strict$corrected_species_id)
  
  # We match the species in the dataset to the species in the phylogenetic tree
  missing<-treedata(inputtree,data=spp_obs,warnings=FALSE)
  # We remove species with no data from the tree
  mtree<-missing$phy
  # We remove species not included in the tree from the dataset
  speciesnames<-mtree$tip.label
  mdata<-dstan_strict[dstan_strict$corrected_species_id %in% speciesnames,]
  
  # We calculate the pair-wise phylogenetic distances among the species in the tree; this gives a matrix where values are the total branch length needed to connect two species (diagonal is the distance of a species to itself so zero)
  Dmat<-cophenetic(mtree)
  
  
  dat_list_strict_phylogenetic <- list(
    R = as.integer(as.factor(dstan_strict$strictfdom)),
    origin = as.integer(as.factor(dstan_strict$discretepredictor)),
    species = as.integer(as.factor(dstan_strict$corrected_species_id))
  )
  
  dat_list_strict_phylogenetic$Dmat<-Dmat[ unique(dstan_strict$corrected_species_id),unique(dstan_strict$corrected_species_id) ]/max(Dmat)
  
  colnames(dat_list_strict_phylogenetic$Dmat)<-as.integer(as.factor(colnames(dat_list_strict_phylogenetic$Dmat)))
  rownames(dat_list_strict_phylogenetic$Dmat)<-as.integer(as.factor(rownames(dat_list_strict_phylogenetic$Dmat)))
  
  dat_list_strict_phylogenetic$N_spp<-length(unique(dstan_strict$corrected_species_id))
  
  
  
  m_strict_discrete_phylogenetic <- ulam(
    alist(
      R ~ dordlogit( phi , cutpoints ),
      phi <- k[species]+b[origin] ,
      b[origin] ~ dnorm( 0 , 0.5 ),
      cutpoints ~ dnorm( 0 , 1.5 ),
      vector[N_spp]:k ~ multi_normal( 0 , SIGMA ),
      matrix[N_spp,N_spp]: SIGMA <- cov_GPL2( Dmat, etasq, rhosq, 0.01),
      etasq~exponential(1),
      rhosq~exponential(1)
    ) , data=dat_list_strict_phylogenetic , chains=4 , cores=4 ,cmdstan=T, messages=FALSE, refresh=0)
  
  
  posterior<-extract.samples(m_strict_discrete_phylogenetic)
  
  for(variants in 1:numberofcontrasts){
    
    logitthiscontrast<-inv_logit(posterior$b[,combn(c(1:totalvariants),2)[2,variants]])-inv_logit(posterior$b[,combn(c(1:totalvariants),2)[1,variants]])
    results[count,4]<-precis(as.numeric(logitthiscontrast))[1,3]
    results[count,5]<-precis(as.numeric(logitthiscontrast))[1,4]
    results[count,3]<-"Yes"
    count<-count+1
    
  }
  
  print("finished m_strict_discrete_phylogenetic")
  
  
  # Female dominance classification of intersexual dominance
  dat_list_relaxed <- list(
    R = as.integer(as.factor(dstan_strict$strictfdom %in% 3))-1,
    origin = as.integer(as.factor(dstan_strict$discretepredictor)),
    species = as.integer(as.factor(dstan_strict$corrected_species_id))
  )
  
  
  m_relaxed_discrete <- ulam(
    alist(
      R ~ dbinom( 1 , p ),
      logit(p) <- b[origin] ,
      b[origin] ~ dnorm( 0 , 0.5 ),
      cutpoints ~ dnorm( 0 , 1.5 )
    ) , data=dat_list_relaxed , chains=4 , cores=4 , cmdstan = T, messages=FALSE )
  
  posterior<-extract.samples(m_relaxed_discrete)
  
  for(variants in 1:numberofcontrasts){
    
    logitthiscontrast<-inv_logit(posterior$b[,combn(c(1:totalvariants),2)[2,variants]])-inv_logit(posterior$b[,combn(c(1:totalvariants),2)[1,variants]])
    results[count,4]<-precis(as.numeric(logitthiscontrast))[1,3]
    results[count,5]<-precis(as.numeric(logitthiscontrast))[1,4]
    results[count,3]<-"No"
    count<-count+1
    
  }
  
  print("finished m_femaledominance_discrete")
  
  # Female dominance classification of intersexual dominance phylogenetic
  spp_obs<-unique(dstan_strict$corrected_species_id)
  spp_obs<-matrix(nrow=length(spp_obs))
  rownames(spp_obs)<-unique(dstan_strict$corrected_species_id)
  
  # We match the species in the dataset to the species in the phylogenetic tree
  missing<-treedata(inputtree,data=spp_obs,warnings=FALSE)
  # We remove species with no data from the tree
  mtree<-missing$phy
  # We remove species not included in the tree from the dataset
  speciesnames<-mtree$tip.label
  mdata<-dstan_strict[dstan_strict$corrected_species_id %in% speciesnames,]
  
  # We calculate the pair-wise phylogenetic distances among the species in the tree; this gives a matrix where values are the total branch length needed to connect two species (diagonal is the distance of a species to itself so zero)
  Dmat<-cophenetic(mtree)
  
  dat_list_relaxed_phylogenetic <- list(
    N_spp = nrow(dstan_strict),
    R = as.integer(as.factor(dstan_strict$strictfdom %in% 3))-1,
    origin = as.integer(as.factor(dstan_strict$discretepredictor)),
    species = as.integer(as.factor(dstan_strict$corrected_species_id))
  )
  
  dat_list_relaxed_phylogenetic$Dmat<-Dmat[ unique(dstan_strict$corrected_species_id),unique(dstan_strict$corrected_species_id) ]/max(Dmat)
  
  colnames(dat_list_relaxed_phylogenetic$Dmat)<-as.integer(as.factor(colnames(dat_list_relaxed_phylogenetic$Dmat)))
  rownames(dat_list_relaxed_phylogenetic$Dmat)<-as.integer(as.factor(rownames(dat_list_relaxed_phylogenetic$Dmat)))
  dat_list_relaxed_phylogenetic$N_spp<-length(unique(dstan_strict$corrected_species_id))
  
  
  m_relaxed_discrete_phylogenetic <- ulam(
    alist(
      R ~ dbinom( 1 , p ),
      logit(p) <- k[species]+b[origin] ,
      b[origin] ~ dnorm( 0 , 0.5 ),
      cutpoints ~ dnorm( 0 , 1.5 ),
      vector[N_spp]:k ~ multi_normal( 0 , SIGMA ),
      matrix[N_spp,N_spp]: SIGMA <- cov_GPL2( Dmat, etasq, rhosq, 0.01),
      etasq~exponential(1),
      rhosq~exponential(1)
    ) , data=dat_list_relaxed_phylogenetic , chains=4 , cores=4, cmdstan=T, messages=FALSE )
  
  posterior<-extract.samples(m_relaxed_discrete_phylogenetic)
  
  for(variants in 1:numberofcontrasts){
    
    logitthiscontrast<-inv_logit(posterior$b[,combn(c(1:totalvariants),2)[2,variants]])-inv_logit(posterior$b[,combn(c(1:totalvariants),2)[1,variants]])
    results[count,4]<-precis(as.numeric(logitthiscontrast))[1,3]
    results[count,5]<-precis(as.numeric(logitthiscontrast))[1,4]
    results[count,3]<-"Yes"
    count<-count+1
      }
  
  print("finished m_femaledominance_discrete_phylogenetic")
  
  
  # Male dominance classification of intersexual dominance
  dat_list_relaxed <- list(
    R = as.integer(as.factor(dstan_strict$strictfdom %in% 1))-1,
    origin = as.integer(as.factor(dstan_strict$discretepredictor)),
    species = as.integer(as.factor(dstan_strict$corrected_species_id))
  )
  
  
  m_relaxed_discrete <- ulam(
    alist(
      R ~ dbinom(1 , p ),
      logit(p) <- b[origin] ,
      b[origin] ~ dnorm( 0 , 0.5 ),
      cutpoints ~ dnorm( 0 , 1.5 )
    ) , data=dat_list_relaxed , chains=4 , cores=4 , cmdstan = T, messages=FALSE )
  
  posterior<-extract.samples(m_relaxed_discrete)
  
  for(variants in 1:numberofcontrasts){
    
    logitthiscontrast<-inv_logit(posterior$b[,combn(c(1:totalvariants),2)[2,variants]])-inv_logit(posterior$b[,combn(c(1:totalvariants),2)[1,variants]])
    results[count,4]<-precis(as.numeric(logitthiscontrast))[1,3]
    results[count,5]<-precis(as.numeric(logitthiscontrast))[1,4]
    results[count,3]<-"No"
    count<-count+1
    
  }
  
  print("finished m_maledominance_discrete")
  
  
  # Male dominance classification of intersexual dominance phylogenetic
  spp_obs<-unique(dstan_strict$corrected_species_id)
  spp_obs<-matrix(nrow=length(spp_obs))
  rownames(spp_obs)<-unique(dstan_strict$corrected_species_id)
  
  # We match the species in the dataset to the species in the phylogenetic tree
  missing<-treedata(inputtree,data=spp_obs,warnings=FALSE)
  # We remove species with no data from the tree
  mtree<-missing$phy
  # We remove species not included in the tree from the dataset
  speciesnames<-mtree$tip.label
  mdata<-dstan_strict[dstan_strict$corrected_species_id %in% speciesnames,]
  
  # We calculate the pair-wise phylogenetic distances among the species in the tree; this gives a matrix where values are the total branch length needed to connect two species (diagonal is the distance of a species to itself so zero)
  Dmat<-cophenetic(mtree)
  
  dat_list_relaxed_phylogenetic <- list(
    N_spp = nrow(dstan_strict),
    R = as.integer(as.factor(dstan_strict$strictfdom %in% 1))-1,
    origin = as.integer(as.factor(dstan_strict$discretepredictor)),
    species = as.integer(as.factor(dstan_strict$corrected_species_id))
  )
  
  dat_list_relaxed_phylogenetic$Dmat<-Dmat[ unique(dstan_strict$corrected_species_id),unique(dstan_strict$corrected_species_id) ]/max(Dmat)
  
  colnames(dat_list_relaxed_phylogenetic$Dmat)<-as.integer(as.factor(colnames(dat_list_relaxed_phylogenetic$Dmat)))
  rownames(dat_list_relaxed_phylogenetic$Dmat)<-as.integer(as.factor(rownames(dat_list_relaxed_phylogenetic$Dmat)))
  dat_list_relaxed_phylogenetic$N_spp<-length(unique(dstan_strict$corrected_species_id))
  
  
  m_relaxed_discrete_phylogenetic <- ulam(
    alist(
      R ~ dbinom( 1 , p ),
      logit(p) <- k[species]+b[origin] ,
      b[origin] ~ dnorm( 0 , 0.5 ),
      cutpoints ~ dnorm( 0 , 1.5 ),
      vector[N_spp]:k ~ multi_normal( 0 , SIGMA ),
      matrix[N_spp,N_spp]: SIGMA <- cov_GPL2( Dmat, etasq, rhosq, 0.01),
      etasq~exponential(1),
      rhosq~exponential(1)
    ) , data=dat_list_relaxed_phylogenetic , chains=4 , cores=4, cmdstan=T, messages=FALSE )
  
  posterior<-extract.samples(m_relaxed_discrete_phylogenetic)
  
  for(variants in 1:numberofcontrasts){
    
    logitthiscontrast<-inv_logit(posterior$b[,combn(c(1:totalvariants),2)[2,variants]])-inv_logit(posterior$b[,combn(c(1:totalvariants),2)[1,variants]])
    results[count,4]<-precis(as.numeric(logitthiscontrast))[1,3]
    results[count,5]<-precis(as.numeric(logitthiscontrast))[1,4]
    results[count,3]<-"Yes"
    count<-count+1
  }
  
  print("finished m_maledominance_discrete_phylogenetic")
  
  
  
  
  
  
  colnames(combined)[colnames(combined) %in% "discretepredictor"]<-c(discretepredictor)
  
  
  
  return(results)
  
}
