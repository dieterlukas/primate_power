###### Analyses of variation within species versus variation between species


###### 1. Females

# females nested in species

filtereddata<-select(combined,corrected_species_id,strictfdom,females)

filtereddata<-filtereddata[complete.cases(filtereddata),]

filtereddata<-filter(filtereddata,corrected_species_id %in% names(table(filtereddata$corrected_species_id)[table(filtereddata$corrected_species_id)>1]))


filtereddata<-filter(filtereddata,corrected_species_id %in% as.data.frame(filtereddata %>% group_by(corrected_species_id) %>% summarize(unique(strictfdom)))[as.data.frame(filtereddata %>% group_by(corrected_species_id) %>% summarize(n_distinct(strictfdom)))[,2]>1,1])




dat_list_strict <- list(
  R = as.integer(as.factor(filtereddata$strictfdom)),
  females = standardize(filtereddata$females),
  species = as.integer(as.factor(filtereddata$corrected_species_id))
)

m_withinspecies_females <- ulam(
  alist(
    R ~ dordlogit( phi , cutpoints ),
    phi <-a_species[species] + b_species[species]*females ,
    c(a_species,b_species)[species]~multi_normal(c(a,b),Rho,sigma_cafe),
    a ~ normal( 0 , 5 ),
    b ~ dnorm(0,5),
    sigma_cafe ~  dexp(1),
    Rho ~ lkj_corr(2),
    cutpoints ~ dnorm( 0 , 5 )
  ) , data=dat_list_strict , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)

precis(m_withinspecies_females)


# Females across species
speciesaverage<-as.data.frame(combined %>% group_by(corrected_species_id) %>% summarise(mean(females,na.rm=T),mean(as.numeric(strictfdom))))

colnames(speciesaverage)<-c("corrected_species_id","females","strictfdom")
speciesaverage<-speciesaverage[is.na(speciesaverage$females)==F,]
speciesaverage$strictfdom<-round(speciesaverage$strictfdom,0)

dat_list_strict <- list(
  R = as.integer(as.factor(speciesaverage$strictfdom)),
  females = standardize(speciesaverage$females),
  species = as.integer(as.factor(speciesaverage$corrected_species_id))
)

m_speciesaverage_females <- ulam(
  alist(
    R ~ dordlogit( phi , cutpoints ),
    phi <-a + b*females ,
    a ~ normal( 0 , 5 ),
    b ~ dnorm(0,5),
    cutpoints ~ dnorm( 0 , 5 )
  ) , data=dat_list_strict , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)

precis(m_speciesaverage_females)


# standardized number of females
females <- seq(from=-1.5,to=3,by=0.5)
pdat <- data.frame(females=females)

plot( NULL , type="n" , xlab="sex ratio" , ylab="probability" ,
      xlim=c(-2.5,2.5) , ylim=c(0,1) ,  yaxp=c(0,1,2) )
phi_speciesaverages <- link( m_speciesaverage_females , data=pdat )
post_speciesaverages <- extract.samples( m_speciesaverage_females )
for ( s in 1:50 ) {
  pk_averages <- pordlogit( 1:3 , phi_speciesaverages[s,] , post_speciesaverages$cutpoints[s,] )
  lines(pk_averages[,1]~females,col="grey80")
  lines(pk_averages[,2]~females,col="grey80")
}

overallphi_speciesaverages<-matrix(ncol=1, nrow=length(females))
for(i in 1:length(females)){
  overallphi_speciesaverages[i,]<-precis(m_speciesaverage_females,depth=2)[1,1]+precis(m_speciesaverage_females,depth=2)[2,1]*females[i]
}

overallprobs_speciesaverages<-pordlogit( 1:3 , overallphi_speciesaverages , precis(m_speciesaverage_females,depth=2)[3:4,1] )
lines(overallprobs_speciesaverages[,1]~females,col="black",lwd=6)
lines(overallprobs_speciesaverages[,2]~females,col="black",lwd=6)


post_withinspecies <- extract.samples( m_withinspecies_females )
phis<-matrix(nrow=50,ncol=length(females))
for(s in 1:50){
  for(i in 1:length(females)){
  phis[s,i]<-post_withinspecies$a[s]+post_withinspecies$b[s]*females[i]
  }
}

for ( s in 1:50 ) {
  pk_within <- pordlogit( 1:3 , phis[s,] , post_withinspecies$cutpoints[s,] )
  lines(pk_within[,1]~females,col="yellow")
  lines(pk_within[,2]~females,col="yellow")
}

overallphi_withinspecies<-matrix(ncol=1, nrow=length(females))
for(i in 1:length(females)){
  overallphi_withinspecies[i,]<-mean(phis[,i])
}

overallprobs_withinspecies<-pordlogit( 1:3 , overallphi_withinspecies , c(mean(post_withinspecies$cutpoints[1:50,1]),mean(post_withinspecies$cutpoints[1:50,2])) )
lines(overallprobs_withinspecies[,1]~females,col="gold",lwd=6)
lines(overallprobs_withinspecies[,2]~females,col="gold",lwd=6)


library(viridis)

plot( NULL , type="n" , xlab="sex ratio in social group" ,
      xlim=c(min(standardize(filtereddata$females),na.rm=T),max(standardize(filtereddata$females),na.rm=T)) , ylim=c(0,1) ,  yaxp=c(0,1,4) ,bty="n",xaxt="n",ylab="")
mtext("Cumulative probability for",side=2,at=0.5,line=3)
mtext("the three dominance systems",side=2,at=0.5,line=2)
axis(side=1,at=c(min(standardize(filtereddata$females),na.rm=T),-1.165669,-0.3334161,0.7069473,1.74729,max(standardize(filtereddata$females),na.rm=T)),labels=FALSE)
mtext("1 male per",side=1,at=min(standardize(filtereddata$females)),line=1)
mtext("15 females",side=1,at=min(standardize(filtereddata$females)),line=2)
mtext("1 male per",side=1,at=-1.165669,line=1)
mtext("4 females",side=1,at=-1.165669,line=2)
mtext("1 male per",side=1,at=-0.3334161,line=1)
mtext("2 females",side=1,at=-0.3334161,line=2)
mtext("1 male per",side=1,at=0.7069473,line=1)
mtext("1 female",side=1,at=0.7069473,line=2)
mtext("2 males per",side=1,at=1.74729,line=1)
mtext("1 female",side=1,at=1.74729,line=2)
mtext("3 males per",side=1,at=max(standardize(filtereddata$females),na.rm=T),line=1)
mtext("1 female",side=1,at=max(standardize(filtereddata$females),na.rm=T),line=2)

polygon(x=c(females,rev(females)),y=c(overallprobs_speciesaverages[,2],rev(overallprobs_speciesaverages[,1])),col=viridis(3)[2],border=NA)

polygon(x=c(females,rev(females)),y=c(overallprobs_speciesaverages[,1],rep(0,11)),col=viridis(3)[3],border=NA)

polygon(x=c(females,rev(females)),y=c(rep(1,11),rev(overallprobs_speciesaverages[,2])),col=viridis(3)[1],border=NA)

mtext("Male dominance",side=1,line=-5,at=-1.75)
mtext("Co dominance",side=1,line=-8,at=0.2)
mtext("Female dominance",side=1,line=-13,at=1.7)





######### 2. Males

# Males nested in species

filtereddata<-select(combined,corrected_species_id,strictfdom,males)

filtereddata<-filtereddata[complete.cases(filtereddata),]

filtereddata<-filter(filtereddata,corrected_species_id %in% names(table(filtereddata$corrected_species_id)[table(filtereddata$corrected_species_id)>1]))


filtereddata<-filter(filtereddata,corrected_species_id %in% as.data.frame(filtereddata %>% group_by(corrected_species_id) %>% summarize(unique(strictfdom)))[as.data.frame(filtereddata %>% group_by(corrected_species_id) %>% summarize(n_distinct(strictfdom)))[,2]>1,1])




dat_list_strict <- list(
  R = as.integer(as.factor(filtereddata$strictfdom)),
  males = standardize(filtereddata$males),
  species = as.integer(as.factor(filtereddata$corrected_species_id))
)

m_withinspecies_males <- ulam(
  alist(
    R ~ dordlogit( phi , cutpoints ),
    phi <-a_species[species] + b_species[species]*males ,
    c(a_species,b_species)[species]~multi_normal(c(a,b),Rho,sigma_cafe),
    a ~ normal( 0 , 5 ),
    b ~ dnorm(0,5),
    sigma_cafe ~  dexp(1),
    Rho ~ lkj_corr(2),
    cutpoints ~ dnorm( 0 , 5 )
  ) , data=dat_list_strict , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)

precis(m_withinspecies_males)


# Number of males across species
speciesaverage<-as.data.frame(combined %>% group_by(corrected_species_id) %>% summarise(mean(males,na.rm=T),mean(as.numeric(strictfdom))))

colnames(speciesaverage)<-c("corrected_species_id","males","strictfdom")
speciesaverage<-speciesaverage[is.na(speciesaverage$males)==F,]
speciesaverage$strictfdom<-round(speciesaverage$strictfdom,0)

dat_list_strict <- list(
  R = as.integer(as.factor(speciesaverage$strictfdom)),
  males = standardize(speciesaverage$males),
  species = as.integer(as.factor(speciesaverage$corrected_species_id))
)

m_speciesaverage_males <- ulam(
  alist(
    R ~ dordlogit( phi , cutpoints ),
    phi <-a + b*males ,
    a ~ normal( 0 , 5 ),
    b ~ dnorm(0,5),
    cutpoints ~ dnorm( 0 , 5 )
  ) , data=dat_list_strict , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)

precis(m_speciesaverage_males)


# standardized number of males
males <- seq(from=-1.5,to=3,by=0.5)
pdat <- data.frame(males=males)

plot( NULL , type="n" , xlab="sex ratio" , ylab="probability" ,
      xlim=c(-2.5,2.5) , ylim=c(0,1) ,  yaxp=c(0,1,2) )
phi_speciesaverages <- link( m_speciesaverage_males , data=pdat )
post_speciesaverages <- extract.samples( m_speciesaverage_males )
for ( s in 1:50 ) {
  pk_averages <- pordlogit( 1:3 , phi_speciesaverages[s,] , post_speciesaverages$cutpoints[s,] )
  lines(pk_averages[,1]~males,col="grey80")
  lines(pk_averages[,2]~males,col="grey80")
}

overallphi_speciesaverages<-matrix(ncol=1, nrow=length(males))
for(i in 1:length(males)){
  overallphi_speciesaverages[i,]<-precis(m_speciesaverage_males,depth=2)[1,1]+precis(m_speciesaverage_males,depth=2)[2,1]*males[i]
}

overallprobs_speciesaverages<-pordlogit( 1:3 , overallphi_speciesaverages , precis(m_speciesaverage_males,depth=2)[3:4,1] )
lines(overallprobs_speciesaverages[,1]~males,col="black",lwd=6)
lines(overallprobs_speciesaverages[,2]~males,col="black",lwd=6)


post_withinspecies <- extract.samples( m_withinspecies_males )
phis<-matrix(nrow=50,ncol=length(males))
for(s in 1:50){
  for(i in 1:length(males)){
    phis[s,i]<-post_withinspecies$a[s]+post_withinspecies$b[s]*males[i]
  }
}

for ( s in 1:50 ) {
  pk_within <- pordlogit( 1:3 , phis[s,] , post_withinspecies$cutpoints[s,] )
  lines(pk_within[,1]~males,col="yellow")
  lines(pk_within[,2]~males,col="yellow")
}

overallphi_withinspecies<-matrix(ncol=1, nrow=length(males))
for(i in 1:length(males)){
  overallphi_withinspecies[i,]<-mean(phis[,i])
}

overallprobs_withinspecies<-pordlogit( 1:3 , overallphi_withinspecies , c(mean(post_withinspecies$cutpoints[1:50,1]),mean(post_withinspecies$cutpoints[1:50,2])) )
lines(overallprobs_withinspecies[,1]~males,col="gold",lwd=6)
lines(overallprobs_withinspecies[,2]~males,col="gold",lwd=6)


library(viridis)

plot( NULL , type="n" , xlab="sex ratio in social group" ,
      xlim=c(min(standardize(filtereddata$males),na.rm=T),max(standardize(filtereddata$males),na.rm=T)) , ylim=c(0,1) ,  yaxp=c(0,1,4) ,bty="n",xaxt="n",ylab="")
mtext("Cumulative probability for",side=2,at=0.5,line=3)
mtext("the three dominance systems",side=2,at=0.5,line=2)
axis(side=1,at=c(min(standardize(filtereddata$males),na.rm=T),-1.165669,-0.3334161,0.7069473,1.74729,max(standardize(filtereddata$males),na.rm=T)),labels=FALSE)
mtext("1 male per",side=1,at=min(standardize(filtereddata$males)),line=1)
mtext("15 males",side=1,at=min(standardize(filtereddata$males)),line=2)
mtext("1 male per",side=1,at=-1.165669,line=1)
mtext("4 males",side=1,at=-1.165669,line=2)
mtext("1 male per",side=1,at=-0.3334161,line=1)
mtext("2 males",side=1,at=-0.3334161,line=2)
mtext("1 male per",side=1,at=0.7069473,line=1)
mtext("1 female",side=1,at=0.7069473,line=2)
mtext("2 males per",side=1,at=1.74729,line=1)
mtext("1 female",side=1,at=1.74729,line=2)
mtext("3 males per",side=1,at=max(standardize(filtereddata$males),na.rm=T),line=1)
mtext("1 female",side=1,at=max(standardize(filtereddata$males),na.rm=T),line=2)

polygon(x=c(males,rev(males)),y=c(overallprobs_speciesaverages[,2],rev(overallprobs_speciesaverages[,1])),col=viridis(3)[2],border=NA)

polygon(x=c(males,rev(males)),y=c(overallprobs_speciesaverages[,1],rep(0,11)),col=viridis(3)[3],border=NA)

polygon(x=c(males,rev(males)),y=c(rep(1,11),rev(overallprobs_speciesaverages[,2])),col=viridis(3)[1],border=NA)

mtext("Male dominance",side=1,line=-5,at=-1.75)
mtext("Co dominance",side=1,line=-8,at=0.2)
mtext("Female dominance",side=1,line=-13,at=1.7)








####### 3. Sex ratio

# Sex ratio nested in species

filtereddata<-combined %>% dplyr::select(corrected_species_id,strictfdom,sexratio)

filtereddata<-filtereddata[complete.cases(filtereddata),]

filtereddata<-filter(filtereddata,corrected_species_id %in% names(table(filtereddata$corrected_species_id)[table(filtereddata$corrected_species_id)>1]))


filtereddata<-filter(filtereddata,corrected_species_id %in% as.data.frame(filtereddata %>% group_by(corrected_species_id) %>% summarize(unique(strictfdom)))[as.data.frame(filtereddata %>% group_by(corrected_species_id) %>% summarize(n_distinct(strictfdom)))[,2]>1,1])




dat_list_strict <- list(
  R = as.integer(as.factor(filtereddata$strictfdom)),
  sexratio = standardize(filtereddata$sexratio),
  species = as.integer(as.factor(filtereddata$corrected_species_id))
)

m_withinspecies_sexratio <- ulam(
  alist(
    R ~ dordlogit( phi , cutpoints ),
    phi <-a_species[species] + b_species[species]*sexratio ,
    c(a_species,b_species)[species]~multi_normal(c(a,b),Rho,sigma_cafe),
    a ~ normal( 0 , 5 ),
    b ~ dnorm(0,5),
    sigma_cafe ~  dexp(1),
    Rho ~ lkj_corr(2),
    cutpoints ~ dnorm( 0 , 5 )
  ) , data=dat_list_strict , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)

precis(m_withinspecies_sexratio)


# Sex ratio across species
speciesaverage<-as.data.frame(combined %>% group_by(corrected_species_id) %>% summarise(mean(sexratio,na.rm=T),mean(as.numeric(strictfdom))))

colnames(speciesaverage)<-c("corrected_species_id","sexratio","strictfdom")
speciesaverage<-speciesaverage[is.na(speciesaverage$sexratio)==F,]
speciesaverage$strictfdom<-round(speciesaverage$strictfdom,0)

dat_list_strict <- list(
  R = as.integer(as.factor(speciesaverage$strictfdom)),
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
  ) , data=dat_list_strict , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)

precis(m_speciesaverage_sexratio)


# standardized sex ratio
sexratio <- seq(from=-2.5,to=2.5,by=0.5)
pdat <- data.frame(sexratio=sexratio)

plot( NULL , type="n" , xlab="sex ratio" , ylab="probability" ,
      xlim=c(-2.5,2.5) , ylim=c(0,1) ,  yaxp=c(0,1,2) )
phi_speciesaverages <- link( m_speciesaverage_sexratio , data=pdat )
post_speciesaverages <- extract.samples( m_speciesaverage_sexratio )
for ( s in 1:50 ) {
  pk_averages <- pordlogit( 1:3 , phi_speciesaverages[s,] , post_speciesaverages$cutpoints[s,] )
  lines(pk_averages[,1]~sexratio,col="grey80")
  lines(pk_averages[,2]~sexratio,col="grey80")
}

overallphi_speciesaverages<-matrix(ncol=1, nrow=11)
for(i in 1:length(sexratio)){
  overallphi_speciesaverages[i,]<-precis(m_speciesaverage_sexratio,depth=2)[1,1]+precis(m_speciesaverage_sexratio,depth=2)[2,1]*sexratio[i]
}

overallprobs_speciesaverages<-pordlogit( 1:3 , overallphi_speciesaverages , precis(m_speciesaverage_sexratio,depth=2)[3:4,1] )
lines(overallprobs_speciesaverages[,1]~sexratio,col="black",lwd=6)
lines(overallprobs_speciesaverages[,2]~sexratio,col="black",lwd=6)


post_withinspecies <- extract.samples( m_withinspecies_sexratio )
phis<-matrix(nrow=50,ncol=11)
for(s in 1:50){
  phis[s,1]<-post_withinspecies$a[s]+post_withinspecies$b[s]*sexratio[1]
  phis[s,2]<-post_withinspecies$a[s]+post_withinspecies$b[s]*sexratio[2]
  phis[s,3]<-post_withinspecies$a[s]+post_withinspecies$b[s]*sexratio[3]
  phis[s,4]<-post_withinspecies$a[s]+post_withinspecies$b[s]*sexratio[4]
  phis[s,5]<-post_withinspecies$a[s]+post_withinspecies$b[s]*sexratio[5]
  phis[s,6]<-post_withinspecies$a[s]+post_withinspecies$b[s]*sexratio[6]
  phis[s,7]<-post_withinspecies$a[s]+post_withinspecies$b[s]*sexratio[7]
  phis[s,8]<-post_withinspecies$a[s]+post_withinspecies$b[s]*sexratio[8]
  phis[s,9]<-post_withinspecies$a[s]+post_withinspecies$b[s]*sexratio[9]
  phis[s,10]<-post_withinspecies$a[s]+post_withinspecies$b[s]*sexratio[10]
  phis[s,11]<-post_withinspecies$a[s]+post_withinspecies$b[s]*sexratio[11]
}

for ( s in 1:50 ) {
  pk_within <- pordlogit( 1:3 , phis[s,] , post_withinspecies$cutpoints[s,] )
  lines(pk_within[,1]~sexratio,col="yellow")
  lines(pk_within[,2]~sexratio,col="yellow")
}

overallphi_withinspecies<-matrix(ncol=1, nrow=11)
for(i in 1:11){
  overallphi_withinspecies[i,]<-mean(phis[,i])
}

overallprobs_withinspecies<-pordlogit( 1:3 , overallphi_withinspecies , c(mean(post_withinspecies$cutpoints[1:50,1]),mean(post_withinspecies$cutpoints[1:50,2])) )
lines(overallprobs_withinspecies[,1]~sexratio,col="gold",lwd=6)
lines(overallprobs_withinspecies[,2]~sexratio,col="gold",lwd=6)


library(viridis)

plot( NULL , type="n" , xlab="sex ratio in social group" ,
      xlim=c(min(standardize(filtereddata$sexratio),na.rm=T),2.579564) , ylim=c(0,1) ,  yaxp=c(0,1,4) ,bty="n",xaxt="n",ylab="")
mtext("Cumulative probability for",side=2,at=0.5,line=3)
mtext("the three dominance systems",side=2,at=0.5,line=2)
axis(side=1,at=c(min(standardize(filtereddata$sexratio),na.rm=T),-1.165669,-0.3334161,0.7069473,1.74729,max(standardize(filtereddata$sexratio),na.rm=T)),labels=FALSE)
mtext("1 male per",side=1,at=min(standardize(filtereddata$sexratio)),line=1)
mtext("15 females",side=1,at=min(standardize(filtereddata$sexratio)),line=2)
mtext("1 male per",side=1,at=-1.165669,line=1)
mtext("4 females",side=1,at=-1.165669,line=2)
mtext("1 male per",side=1,at=-0.3334161,line=1)
mtext("2 females",side=1,at=-0.3334161,line=2)
mtext("1 male per",side=1,at=0.7069473,line=1)
mtext("1 female",side=1,at=0.7069473,line=2)
mtext("2 males per",side=1,at=1.74729,line=1)
mtext("1 female",side=1,at=1.74729,line=2)
mtext("4 males per",side=1,at=2.579564,line=1)
mtext("1 female",side=1,at=2.579564,line=2)

polygon(x=c(sexratio,rev(sexratio)),y=c(overallprobs_speciesaverages[,2],rev(overallprobs_speciesaverages[,1])),col=co_dominance_color,border=NA)

polygon(x=c(sexratio,rev(sexratio)),y=c(overallprobs_speciesaverages[,1],rep(0,11)),col=male_dominance_color,border=NA)

polygon(x=c(sexratio,rev(sexratio)),y=c(rep(1,11),rev(overallprobs_speciesaverages[,2])),col=female_dominance_color,border=NA)

mtext("Male dominance",side=1,line=-5,at=-1.75)
mtext("Co dominance",side=1,line=-8,at=0.2)
mtext("Female dominance",side=1,line=-13,at=1.7)

hist(exp(filtereddata$sexratio),breaks=40)


