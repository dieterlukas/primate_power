# Setting up the functions for the analyses

run_analyses_continuouspredictor_nested <- function(continuouspredictor){
  
  
  colnames(combined)[colnames(combined) %in% continuouspredictor]<-c("continuouspredictor")
  colnames(combined)[colnames(combined) %in% nestingvariable]<-c("nestingvariable")
  
  results<-matrix(nrow=3*2*2,ncol=5)    
  results<-as.data.frame(results)
  colnames(results)<-c("outcome", "nestingvariable","phylogeny", "contrast lower","contrast upper")
  results$outcome<-c(rep("perc_won", 4),rep("strict", 4),rep("relaxed",4))
  results$nestingvariable<-rep(c(paste("no",nestingvariable,continuouspredictor),paste(nestingvariable, continuouspredictor)),6)
  results$phylogeny<-rep(c("no","no","yes","yes"),3)
  
  
  dstan_continuous <- combined[ complete.cases(combined$perc_won_females,combined$continuouspredictor,combined$nestingvariable,combined$corrected_species_id),]
  
  dstan_strict <- combined[ complete.cases(combined$strictfdom,combined$continuouspredictor,combined$nestingvariable,combined$corrected_species_id),]
  
  dstan_relaxed <- combined[ complete.cases(combined$mostlyfdom,combined$continuouspredictor,combined$nestingvariable,combined$corrected_species_id),]  
  
  
  # Continuous outcome: percentage of fights won by females
  dat_list_continuous <- list(
    N_spp = nrow(dstan_continuous),
    perc_won_females = as.integer(dstan_continuous$perc_won_females),
    continuouspredictor = standardize(dstan_continuous$continuouspredictor),
    nestingvariable = as.integer(as.factor(dstan_continuous$nestingvariable)),
    species = as.integer(as.factor(dstan_continuous$corrected_species_id)),
    total = rep(100,nrow(dstan_continuous))
  )
  
  m_continuous_continuous <- ulam(
    alist(
      perc_won_females ~ dbinom(total,p),
      logit(p) <- a[nestingvariable] +b[nestingvariable]*continuouspredictor,
      ## adaptive priors
      a[nestingvariable] ~ dnorm( 0 , 1 ),
      b[nestingvariable] ~ dnorm(0,1)
    ) , data=dat_list_continuous , chains=4 , cores=4 , log_lik=TRUE , cmdstan=T, messages=FALSE)
  
  results[1,4:5]<-precis(m_continuous_continuous,depth=2)[3,3:4]
  results[2,4:5]<-precis(m_continuous_continuous,depth=2)[4,3:4]
  
  
  
  
  
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
    nestingvariable = as.integer(as.factor(dstan_continuous$nestingvariable)),
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
      logit(p) <- k[species]+a[nestingvariable] +b[nestingvariable]*continuouspredictor,
      vector[N_spp]:k ~ multi_normal( 0 , SIGMA ),
      matrix[N_spp,N_spp]: SIGMA <- cov_GPL2( Dmat, etasq, rhosq, 0.01),
      ## adaptive priors
      a[nestingvariable] ~ dnorm(0,1),
      b[nestingvariable] ~ dnorm(0,1),
      etasq~exponential(1),
      rhosq~exponential(1)
    ) , data=dat_list_continuous_phylogenetic , chains=4 , cores=4 , log_lik=TRUE , cmdstan=T, messages=FALSE)
  
  posterior<-extract.samples(m_continuous_continuous_phylogenetic)
  results[3,4:5]<-precis(posterior$b[,1])[1,3:4]
  results[4,4:5]<-precis(posterior$b[,2])[1,3:4]
  
  
  
  
  # Strict categorical classification of intersexual dominance
  
  dat_list_strict <- list(
    R = as.integer(as.factor(dstan_strict$strictfdom)),
    continuouspredictor = standardize(dstan_strict$continuouspredictor),
    nestingvariable = as.integer(as.factor(dstan_strict$nestingvariable)),
    species = as.integer(as.factor(dstan_strict$corrected_species_id))
  )
  
  m_strict_continuous <- ulam(
    alist(
      R ~ dordlogit( phi , cutpoints ),
      phi <-a[nestingvariable] + b[nestingvariable]*continuouspredictor,
      a[nestingvariable] ~ dnorm( 0 , 1 ),
      b[nestingvariable] ~ dnorm(0,1),
      cutpoints ~ dnorm( 0 , 1.5 )
    ) , data=dat_list_strict , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)
  
  posterior<-extract.samples(m_strict_continuous)
  results[5,4:5]<-precis(posterior$b[,1])[1,3:4]
  results[6,4:5]<-precis(posterior$b[,2])[1,3:4]
  
  # Strict categorical classification of intersexual dominance phylogenetic
  
  
  dat_list_strict_phylogenetic <- list(
    R = as.integer(as.factor(dstan_strict$strictfdom)),
    continuouspredictor = standardize(dstan_strict$continuouspredictor),
    nestingvariable = as.integer(as.factor(dstan_strict$nestingvariable)),
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
      phi <-k[species]+a[nestingvariable] + b[nestingvariable]*continuouspredictor,
      vector[N_spp]:k ~ multi_normal( 0 , SIGMA ),
      matrix[N_spp,N_spp]: SIGMA <- cov_GPL2( Dmat, etasq, rhosq, 0.01),
      a[nestingvariable] ~ dnorm(0,1),
      b[nestingvariable] ~ dnorm(0,1),
      cutpoints ~ dnorm( 0 , 1.5 ),
      etasq~exponential(1),
      rhosq~exponential(1)
    ) , data=dat_list_strict_phylogenetic , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)
  
  posterior<-extract.samples(m_strict_continuous_phylogenetic)
  results[7,4:5]<-precis(posterior$b[,1])[1,3:4]
  results[8,4:5]<-precis(posterior$b[,2])[1,3:4]
  
  
  
  
  # Relaxed categorical classification of intersexual dominance
  dat_list_relaxed <- list(
    R = as.integer(as.factor(dstan_relaxed$mostlyfdom))-1,
    continuouspredictor = standardize(dstan_relaxed$continuouspredictor),
    nestingvariable = as.integer(as.factor(dstan_relaxed$nestingvariable)),
    species = as.integer(as.factor(dstan_relaxed$corrected_species_id))
  )
  
  
  m_relaxed_continuouspredictor <- ulam(
    alist(
      R ~ dbinom(1,p),
      logit(p) <- a[nestingvariable] + b[nestingvariable]*continuouspredictor ,
      a[nestingvariable] ~ dnorm( 0 , 1 ),
      b[nestingvariable] ~ dnorm( 0 , 0.5 ),
      cutpoints ~ dnorm( 0 , 1.5 )
    ) , data=dat_list_relaxed , chains=4 , cores=4 , cmdstan = T, messages=FALSE )
  
  posterior<-extract.samples(m_relaxed_continuouspredictor)
  results[9,4:5]<-precis(posterior$b[,1])[1,3:4]
  results[10,4:5]<-precis(posterior$b[,2])[1,3:4]
  
  # Relaxed categorical classification of intersexual dominance phylogenetic
  
  spp_obs<-unique(dstan_relaxed$corrected_species_id)
  spp_obs<-matrix(nrow=length(spp_obs))
  rownames(spp_obs)<-unique(dstan_relaxed$corrected_species_id)
  
  # We match the species in the dataset to the species in the phylogenetic tree
  missing<-treedata(inputtree,data=spp_obs,warnings=FALSE)
  # We remove species with no data from the tree
  mtree<-missing$phy
  
  Dmat<-cophenetic(mtree)
  
  dat_list_relaxed_phylogenetic <- list(
    N_spp = nrow(dstan_relaxed),
    R = as.integer(as.factor(dstan_relaxed$mostlyfdom)),
    continuouspredictor = standardize(dstan_relaxed$continuouspredictor),
    nestingvariable = as.integer(as.factor(dstan_relaxed$nestingvariable)),
    species = as.integer(as.factor(dstan_relaxed$corrected_species_id))
  )
  
  dat_list_relaxed_phylogenetic$Dmat<-Dmat[ unique(dstan_relaxed$corrected_species_id),unique(dstan_relaxed$corrected_species_id) ]/max(Dmat)
  
  colnames(dat_list_relaxed_phylogenetic$Dmat)<-as.integer(as.factor(colnames(dat_list_relaxed_phylogenetic$Dmat)))
  rownames(dat_list_relaxed_phylogenetic$Dmat)<-as.integer(as.factor(rownames(dat_list_relaxed_phylogenetic$Dmat)))
  
  dat_list_relaxed_phylogenetic$N_spp<-length(unique(dstan_relaxed$corrected_species_id))
  
  
  
  m_relaxed_continuous_phylogenetic <- ulam(
    alist(
      R ~ dordlogit( phi , cutpoints ),
      phi <-k[species] + b[nestingvariable]*continuouspredictor,
      vector[N_spp]:k ~ multi_normal( 0 , SIGMA ),
      matrix[N_spp,N_spp]: SIGMA <- cov_GPL2( Dmat, etasq, rhosq, 0.01),
      b[nestingvariable] ~ dnorm(0,1),
      cutpoints ~ dnorm( 0 , 1.5 ),
      etasq~exponential(1),
      rhosq~exponential(1)
    ) , data=dat_list_relaxed_phylogenetic , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)
  
  posterior<-extract.samples(m_relaxed_continuous_phylogenetic)
  results[11,4:5]<-precis(posterior$b[,1])[1,3:4]
  results[12,4:5]<-precis(posterior$b[,2])[1,3:4]
  
  colnames(combined)[colnames(combined) %in% "continuouspredictor"]<-continuouspredictor
  colnames(combined)[colnames(combined) %in% "nestingvariable"]<-nestingvariable
  
  
  return(results)
  
}






run_analyses_discretepredictor_nested <- function(discretepredictor){
  
  colnames(combined)[colnames(combined) %in% discretepredictor]<-c("discretepredictor")
  colnames(combined)[colnames(combined) %in% nestingvariable]<-c("nestingvariable")
  
  dstan_continuous <- combined[ complete.cases(combined$perc_won_females,combined$discretepredictor,combined$corrected_species_id),]
  
  dstan_strict <- combined[ complete.cases(combined$strictfdom,combined$discretepredictor,combined$corrected_species_id),]
  
  dstan_relaxed <- combined[ complete.cases(combined$mostlyfdom,combined$discretepredictor,combined$corrected_species_id),]  
  
  
  totalvariants<-length(unique(dstan_continuous$discretepredictor))
  numberofcontrasts<-totalvariants*(totalvariants-1)/2
  
  results<-matrix(nrow=numberofcontrasts*3*2*2,ncol=5)    
  results<-as.data.frame(results)
  colnames(results)<-c("outcome", "nestingvariable","phylogeny", "contrast lower","contrast upper")
  results$outcome<-c(rep("perc_won", numberofcontrasts*4),rep("strict", numberofcontrasts*4),rep("relaxed",numberofcontrasts*4))
  results$nestingvariable<-rep(c(paste("no",nestingvariable,discretepredictor),paste(nestingvariable, discretepredictor)),numberofcontrasts*6)
  
  count<-1
  
  
  
  # Continuous outcome: percentage of fights won by females
  dat_list_continuous <- list(
    perc_won_females = as.integer(dstan_continuous$perc_won_females),
    discretepredictor = as.integer(as.factor(dstan_continuous$discretepredictor)),
    nestingvariable = as.integer(as.factor(dstan_continuous$nestingvariable)),
    species = as.integer(as.factor(dstan_continuous$corrected_species_id)),
    total = rep(100,nrow(dstan_continuous))
  )
  m_continuous_discrete <- ulam(
    alist(
      perc_won_females ~ dbinom(total,p),
      logit(p) <- a[nestingvariable] +b[nestingvariable]*discretepredictor,
      ## adaptive priors
      a[nestingvariable] ~ dnorm( 0 , 1 ),
      b[nestingvariable] ~ dnorm(0,1)
    ) , data=dat_list_continuous , chains=4 , cores=4 , log_lik=TRUE , cmdstan=T, messages=FALSE)
  
  posterior<-extract.samples(m_continuous_discrete)
  
  
  results[count,4]<-precis(posterior$b[,1])[1,3]
  results[count,5]<-precis(posterior$b[,1])[1,4]
  results[count,3]<-"No"
  count<-count+1
  results[count,4]<-precis(posterior$b[,2])[1,3]
  results[count,5]<-precis(posterior$b[,2])[1,4]
  results[count,3]<-"No"
  count<-count+1
  
  
  
  
  
  
  
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
    nestingvariable = as.integer(as.factor(dstan_continuous$nestingvariable)),
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
      logit(p) <- k[species]+a[nestingvariable] +b[nestingvariable]*discretepredictor,
      vector[N_spp]:k ~ multi_normal( 0 , SIGMA ),
      matrix[N_spp,N_spp]: SIGMA <- cov_GPL2( Dmat, etasq, rhosq, 0.01),
      a[nestingvariable] ~ dnorm(0,1),
      b[nestingvariable] ~ dnorm(0,1),
      etasq~exponential(1),
      rhosq~exponential(1)
    ) , data=dat_list_continuous_discrete_phylogenetic , chains=4 , cores=4 , log_lik=TRUE , cmdstan=T, messages=FALSE)
  
  posterior<-extract.samples(m_continuous_discrete_phylogenetic)
  
  
  results[count,4]<-precis(posterior$b[,1])[1,3]
  results[count,5]<-precis(posterior$b[,1])[1,4]
  results[count,3]<-"Yes"
  count<-count+1
  results[count,4]<-precis(posterior$b[,2])[1,3]
  results[count,5]<-precis(posterior$b[,2])[1,4]
  results[count,3]<-"Yes"
  count<-count+1
  
  
  # Strict categorical classification of intersexual dominance
  
  dat_list <- list(
    R = as.integer(as.factor(dstan_strict$strictfdom)),
    origin = as.integer(as.factor(dstan_strict$discretepredictor)),
    nestingvariable = as.integer(as.factor(dstan_strict$nestingvariable)),
    species = as.integer(as.factor(dstan_strict$corrected_species_id))
  )
  
  m_strict_discrete <- ulam(
    alist(
      R ~ dordlogit( phi , cutpoints ),
      phi <- a[nestingvariable]+b[nestingvariable]*origin ,
      a[nestingvariable] ~ dnorm( 0 , 0.5 ),
      b[nestingvariable] ~ dnorm( 0 , 0.5 ),
      cutpoints ~ dnorm( 0 , 1.5 )
    ) , data=dat_list , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)
  
  posterior<-extract.samples(m_strict_discrete)
  
  
  results[count,4]<-precis(posterior$b[,1])[1,3]
  results[count,5]<-precis(posterior$b[,1])[1,4]
  results[count,3]<-"No"
  count<-count+1
  results[count,4]<-precis(posterior$b[,2])[1,3]
  results[count,5]<-precis(posterior$b[,2])[1,4]
  results[count,3]<-"No"
  count<-count+1
  
  
  
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
    nestingvariable = as.integer(as.factor(dstan_strict$nestingvariable)),
    species = as.integer(as.factor(dstan_strict$corrected_species_id))
  )
  
  dat_list_strict_phylogenetic$Dmat<-Dmat[ unique(dstan_strict$corrected_species_id),unique(dstan_strict$corrected_species_id) ]/max(Dmat)
  
  colnames(dat_list_strict_phylogenetic$Dmat)<-as.integer(as.factor(colnames(dat_list_strict_phylogenetic$Dmat)))
  rownames(dat_list_strict_phylogenetic$Dmat)<-as.integer(as.factor(rownames(dat_list_strict_phylogenetic$Dmat)))
  
  dat_list_strict_phylogenetic$N_spp<-length(unique(dstan_strict$corrected_species_id))
  
  
  
  m_strict_discrete_phylogenetic <- ulam(
    alist(
      R ~ dordlogit( phi , cutpoints ),
      phi <- k[species]+a[nestingvariable]+b[nestingvariable]*origin ,
      a[nestingvariable] ~ dnorm( 0 , 0.5 ),
      b[nestingvariable] ~ dnorm( 0 , 0.5 ),
      cutpoints ~ dnorm( 0 , 1.5 ),
      vector[N_spp]:k ~ multi_normal( 0 , SIGMA ),
      matrix[N_spp,N_spp]: SIGMA <- cov_GPL2( Dmat, etasq, rhosq, 0.01),
      etasq~exponential(1),
      rhosq~exponential(1)
    ) , data=dat_list_strict_phylogenetic , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)
  
  
  posterior<-extract.samples(m_strict_discrete_phylogenetic)
  
  
  results[count,4]<-precis(posterior$b[,1])[1,3]
  results[count,5]<-precis(posterior$b[,1])[1,4]
  results[count,3]<-"Yes"
  count<-count+1
  results[count,4]<-precis(posterior$b[,2])[1,3]
  results[count,5]<-precis(posterior$b[,2])[1,4]
  results[count,3]<-"Yes"
  count<-count+1
  
  
  
  
  # Relaxed categorical classification of intersexual dominance
  dat_list_relaxed <- list(
    R = as.integer(as.factor(dstan_relaxed$mostlyfdom)),
    origin = as.integer(as.factor(dstan_relaxed$discretepredictor)),
    nestingvariable = as.integer(as.factor(dstan_relaxed$nestingvariable)),
    species = as.integer(as.factor(dstan_relaxed$corrected_species_id))
  )
  
  
  m_relaxed_discrete <- ulam(
    alist(
      R ~ dordlogit( phi , cutpoints ),
      phi <- a[nestingvariable]+b[nestingvariable]*origin ,
      a[nestingvariable] ~ dnorm( 0 , 0.5 ),    
      b[nestingvariable] ~ dnorm( 0 , 0.5 ),
      cutpoints ~ dnorm( 0 , 1.5 )
    ) , data=dat_list_relaxed , chains=4 , cores=4 , cmdstan = T, messages=FALSE )
  
  posterior<-extract.samples(m_relaxed_discrete)
  
  results[count,4]<-precis(posterior$b[,1])[1,3]
  results[count,5]<-precis(posterior$b[,1])[1,4]
  results[count,3]<-"No"
  count<-count+1
  results[count,4]<-precis(posterior$b[,2])[1,3]
  results[count,5]<-precis(posterior$b[,2])[1,4]
  results[count,3]<-"No"
  count<-count+1
  
  # Relaxed categorical classification of intersexual dominance phylogenetic
  spp_obs<-unique(dstan_relaxed$corrected_species_id)
  spp_obs<-matrix(nrow=length(spp_obs))
  rownames(spp_obs)<-unique(dstan_relaxed$corrected_species_id)
  
  # We match the species in the dataset to the species in the phylogenetic tree
  missing<-treedata(inputtree,data=spp_obs,warnings=FALSE)
  # We remove species with no data from the tree
  mtree<-missing$phy
  # We remove species not included in the tree from the dataset
  speciesnames<-mtree$tip.label
  mdata<-dstan_relaxed[dstan_relaxed$corrected_species_id %in% speciesnames,]
  
  # We calculate the pair-wise phylogenetic distances among the species in the tree; this gives a matrix where values are the total branch length needed to connect two species (diagonal is the distance of a species to itself so zero)
  Dmat<-cophenetic(mtree)
  
  dat_list_relaxed_phylogenetic <- list(
    N_spp = nrow(dstan_relaxed),
    R = as.integer(as.factor(dstan_relaxed$mostlyfdom)),
    origin = as.integer(as.factor(dstan_relaxed$discretepredictor)),
    nestingvariable = as.integer(as.factor(dstan_relaxed$nestingvariable)),
    species = as.integer(as.factor(dstan_relaxed$corrected_species_id))
  )
  
  dat_list_relaxed_phylogenetic$Dmat<-Dmat[ unique(dstan_relaxed$corrected_species_id),unique(dstan_relaxed$corrected_species_id) ]/max(Dmat)
  
  colnames(dat_list_relaxed_phylogenetic$Dmat)<-as.integer(as.factor(colnames(dat_list_relaxed_phylogenetic$Dmat)))
  rownames(dat_list_relaxed_phylogenetic$Dmat)<-as.integer(as.factor(rownames(dat_list_relaxed_phylogenetic$Dmat)))
  dat_list_relaxed_phylogenetic$N_spp<-length(unique(dstan_relaxed$corrected_species_id))
  
  
  m_relaxed_discrete_phylogenetic <- ulam(
    alist(
      R ~ dordlogit( phi , cutpoints ),
      phi <- k[species]+a[nestingvariable]+b[nestingvariable]*origin ,
      a[nestingvariable] ~ dnorm( 0 , 0.5 ),
      b[nestingvariable] ~ dnorm( 0 , 0.5 ),
      cutpoints ~ dnorm( 0 , 1.5 ),
      vector[N_spp]:k ~ multi_normal( 0 , SIGMA ),
      matrix[N_spp,N_spp]: SIGMA <- cov_GPL2( Dmat, etasq, rhosq, 0.01),
      etasq~exponential(1),
      rhosq~exponential(1)
    ) , data=dat_list_relaxed_phylogenetic , chains=4 , cores=4, cmdstan=T, messages=FALSE )
  
  posterior<-extract.samples(m_relaxed_discrete_phylogenetic)
  
  results[count,4]<-precis(posterior$b[,1])[1,3]
  results[count,5]<-precis(posterior$b[,1])[1,4]
  results[count,3]<-"Yes"
  count<-count+1
  results[count,4]<-precis(posterior$b[,2])[1,3]
  results[count,5]<-precis(posterior$b[,2])[1,4]
  results[count,3]<-"Yes"
  count<-count+1
  
  colnames(combined)[colnames(combined) %in% "discretepredictor"]<-c(discretepredictor)
  colnames(combined)[colnames(combined) %in% "nestingvariable"]<-c(nestingvariable)
  
  
  return(results)
  
}
