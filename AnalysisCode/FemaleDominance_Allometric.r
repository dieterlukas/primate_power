# Analyses of variables that have allometric scaling

# Testis size ~ intersexual dominance + body size

filtereddata<-select(combined,corrected_species_id,perc_won_females,testes_mass,male_mass)

filtereddata<-filtereddata[complete.cases(filtereddata),]


dat_list_strict <- list(
  R = as.integer(as.factor(filtereddata$strictfdom)),
  testismass = standardize(filtereddata$testes_mass),
  bodymass = standardize(filtereddata$male_mass),
  perc_won_females =standardize(filtereddata$perc_won_females)
)

m_relativetestes <- ulam(
  alist(
    testismass ~ dnorm( mu , sigma ),
    mu <-a + b*bodymass + c*perc_won_females,
    a ~ dnorm( 0 , 1 ),
    c(b,c) ~ dnorm(0,2),
    sigma ~  dexp(1)
  ) , data=dat_list_strict , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)

precis(m_relativetestes,depth=2)
# Negative relationship: -0.32  -0.01, meaning that females win more fights when males have smaller testes.


tree_trimmed<-keep.tip(inputtree,filtereddata$corrected_species_id)
Dmat<-cophenetic(tree_trimmed)

dat_list_strict$Dmat<-Dmat[filtereddata$corrected_species_id,filtereddata$corrected_species_id]
dat_list_strict$N_spp=nrow(filtereddata)

m_relativetestes <- ulam(
  alist(
    testismass ~ multi_normal( mu , SIGMA ),
    mu <-a + b*bodymass + c*perc_won_females,
    matrix[N_spp,N_spp]:SIGMA<-cov_GPL2(Dmat,etasq,rhosq,0.01),
    a ~ dnorm( 0 , 1 ),
    c(b,c) ~ dnorm(0,2),
    etasq~half_normal(1,0.25),
    rhosq~half_normal(3,0.25)
  ) , data=dat_list_strict , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)

precis(m_relativetestes,depth=2)
# No relationship -0.08  0.06



filtereddata<-select(combined,corrected_species_id,strictfdom,testes_mass,male_mass)

filtereddata<-filtereddata[complete.cases(filtereddata),]


dat_list_strict <- list(
  R = as.integer(as.factor(filtereddata$strictfdom)),
  testismass = standardize(filtereddata$testes_mass),
  bodymass = standardize(filtereddata$male_mass),
  strictfdom = as.integer(as.factor(filtereddata$strictfdom))
)

m_relativetestes <- ulam(
  alist(
    testismass ~ dnorm( mu , sigma ),
    mu <-a + b*bodymass + c[strictfdom],
    a ~ dnorm( 0 , 1 ),
    b ~ dnorm(0,2),
    c[strictfdom] ~ dnorm(0,2),
    sigma ~  dexp(1)
  ) , data=dat_list_strict , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)

precis(m_relativetestes,depth=2)

posterior_relativetestes<-extract.samples(m_relativetestes)
fdom_vs_no<-posterior_relativetestes$c[,3]-posterior_relativetestes$c[,2]
mdom_vs_no<-posterior_relativetestes$c[,1]-posterior_relativetestes$c[,2]
fdom_vs_mdom<-posterior_relativetestes$c[,3]-posterior_relativetestes$c[,1]
# No relationship





# Female canine size ~ intersexual dominance + body size

filtereddata<-select(combined,corrected_species_id,perc_won_females,female_canine_height,female_mass)

filtereddata<-filtereddata[complete.cases(filtereddata),]


dat_list_strict <- list(
  R = as.integer(as.factor(filtereddata$strictfdom)),
  female_canine_height = standardize(log(filtereddata$female_canine_height)),
  bodymass = standardize(filtereddata$female_mass)),
  perc_won_females =standardize(filtereddata$perc_won_females)
)

m_relativetestes <- ulam(
  alist(
    testismass ~ dnorm( mu , sigma ),
    mu <-a + b*bodymass + c*perc_won_females,
    a ~ dnorm( 0 , 1 ),
    c(b,c) ~ dnorm(0,2),
    sigma ~  dexp(1)
  ) , data=dat_list_strict , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)

precis(m_relativetestes,depth=2)
# No relationship: -0.23  0.08


filtereddata<-select(combined,corrected_species_id,strictfdom,testes_mass,body_mass)

filtereddata<-filtereddata[complete.cases(filtereddata),]


dat_list_strict <- list(
  R = as.integer(as.factor(filtereddata$strictfdom)),
  testismass = standardize(filtereddata$testes_mass),
  bodymass = as.integer(as.factor(log(filtereddata$body_mass))),
  strictfdom = as.integer(as.factor(filtereddata$strictfdom))
)

m_relativetestes <- ulam(
  alist(
    testismass ~ dnorm( mu , sigma ),
    mu <-a + b*bodymass + c[strictfdom],
    a ~ dnorm( 0 , 1 ),
    b ~ dnorm(0,2),
    c[strictfdom] ~ dnorm(0,2),
    sigma ~  dexp(1)
  ) , data=dat_list_strict , chains=4 , cores=4 ,cmdstan=T, messages=FALSE)

precis(m_relativetestes,depth=2)

posterior_relativetestes<-extract.samples(m_relativetestes)
fdom_vs_no<-posterior_relativetestes$c[,3]-posterior_relativetestes$c[,2]
mdom_vs_no<-posterior_relativetestes$c[,1]-posterior_relativetestes$c[,2]
fdom_vs_mdom<-posterior_relativetestes$c[,3]-posterior_relativetestes$c[,1]
# No relationship