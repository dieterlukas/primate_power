

# Loading the data

# We load the phylogenetic tree
ifelse(fileorigin=="github",inputtree<-read.tree("https://raw.githubusercontent.com/dieterlukas/primate_power/main/data/PrimatePhylogeny_IntersexDominance.tre"),inputtree<-read.tree("data/PrimatePhylogeny_IntersexDominance.tre"))


# We load the population data with the response variables
ifelse(fileorigin=="github",inputdata<-read.table("https://raw.githubusercontent.com/dieterlukas/primate_power/main/data/PopulationData_IntersexDominance.csv",sep=",",header=T),inputdata<-read.table("data/PopulationData_IntersexDominance.csv",sep=",",header=T))


# We load the species-level predictor variables
ifelse(fileorigin=="github",specieslevelpredictors<-read.table("https://raw.githubusercontent.com/dieterlukas/primate_power/main/data/SpeciesLevelPredictors_IntersexDominance.txt"),specieslevelpredictors<-read.table("data/SpeciesLevelPredictors_IntersexDominance.txt"))

# We add the species names as a separate variable in a separate column for matching
specieslevelpredictors$corrected_species_id<-row.names(specieslevelpredictors)

# We add the species level data to the population level data
combined<-left_join(inputdata,specieslevelpredictors,by="corrected_species_id")

# We filter out those data points that are not of sufficient quality for a comparative analysis
# This combined file will now be the main dataset we refer to for the analyses
combined<-combined[is.na(combined$goodquality)==F,]
combined<-combined[combined$goodquality=="y",]

# We also exclude all those entries where we could not clearly define a general pattern of intersexual interactions
combined<-combined[is.na(combined$strictfdom)==F,]

# We sort the strict intersexual dominance relationships into the correct order (1=male dominance, 2=co-dominance, 3=female dominance)
combined[combined$strictfdom=="f",]$strictfdom<-3
combined[combined$strictfdom=="c",]$strictfdom<-2
combined[combined$strictfdom=="m",]$strictfdom<-1

# We create a new variable that separates lemur from non-lemur species
combined$lemur<-0
combined[combined$corrected_species_id %in% c("Daubentonia_madagascariensis","Eulemur_coronatus","Eulemur_flavifrons","Eulemur_fulvus","Eulemur_macaco","Eulemur_rubriventer","Eulemur_rufifrons","Hapalemur_alaotrensis","Hapalemur_griseus","Hapalemur_meridionalis","Lemur_catta","Leontocebus_tripartitus","Leontocebus_weddelli","Lepilemur_leucopus","Lepilemur_ruficaudatus","Loris_lydekkerianus","Phaner_pallescens","Varecia_rubra","Varecia_variegata","Propithecus_diadema","Propithecus_coquereli","Propithecus_edwardsi","Propithecus_coronatus","Avahi_occidentalis","Indri_indri","Propithecus_verreauxi","Microcebus_myoxinus","Microcebus_ravelobensis","Microcebus_bongolavensis","Microcebus_danfossi","Microcebus_margotmarshae","Microcebus_mamiratra","Microcebus_murinus","Microcebus_lehilahytsara","Microcebus_griseorufus","Lepilemur_edwardsi"),]$lemur<-1

# We create a new variable that splits species depending on the extent of sexual dimorphism in body size
combined$dimorphic<-0
combined[combined$SexualDimorphism_MaleWeight_over_FemaleWeight>1.1,]$dimorphic<-1


# Calculate proportion of dyads among females, among males, and between females and males
combined$ff_dyads<-combined$females*(combined$females-1)/2
combined$mm_dyads<-combined$males*(combined$males-1)/2
combined$fm_dyads<-combined$females*combined$males

combined$perc_ff_dyads<-combined$ff_dyads/(combined$ff_dyads+combined$mm_dyads+combined$fm_dyads)
combined$perc_mm_dyads<-combined$mm_dyads/(combined$ff_dyads+combined$mm_dyads+combined$fm_dyads)
combined$perc_fm_dyads<-combined$fm_dyads/(combined$ff_dyads+combined$mm_dyads+combined$fm_dyads)


# Filter male reproductive skew to only include groups with multiple males
combined[combined$males %in% c(1,1.25,1.3,1.31,1.5),]$male_skew<-NA

# Convert traits where processes can lead to potential long tails in the distribution
combined$females<-log(combined$females)
combined$males<-log(combined$males)
combined$sexualreceptivity_hours<-log(combined$sexualreceptivity_hours)
combined$SexualDimorphism_MaleWeight_over_FemaleWeight<-log(combined$SexualDimorphism_MaleWeight_over_FemaleWeight)
combined$CanineDimorphism<-log(combined$CanineDimorphism)
combined$female_mass<-log(combined$female_mass)
combined$male_mass<-log(combined$male_mass)
combined$testes_mass<-log(combined$testes_mass)
combined$rainfall_annualvariation<-log(combined$rainfall_annualvariation)
combined$relative_testes_mass<-NA
combined[is.na(combined$testes_mass)==F,]$relative_testes_mass<-lm(testes_mass~male_mass,data=combined)$residuals
combined$relative_femalecaninesize<-NA
combined[is.na(combined$female_canine_height)==F,]$relative_femalecaninesize<-lm(female_canine_height~female_mass,data=combined)$residuals



# For the variable ovulation signs, reclassify the one instance of 'slight' as 'present'
combined[combined$ovulation_signs %in% "Slight",]$ovulation_signs <- "Present"


# For the variable origin, reclassify the four instances of 'provisioned' as 'wild'
combined[combined$origin %in% "provisioned",]$origin <- "wild"


names(combined)
# Variables in the dataset

# "X"                        unique identifier for each entry                    
# "corrected_species_id"     Latin species name, matching taxonomy to the names in Burgin et al https://academic.oup.com/jmammal/article/99/1/1/4834091                    
# "entry_id"                  identifier linked to the original datasheet                   
# "origin"                    whether the studied group was in captivity, in the wild, or in the wild but provisioned; data from the original studies                   
# "females"                   the average number of adult females in the study groups; data from the original studies             
# "males"                     the average number of adult males in the study groups; data from the original studies                 
# "goodquality"               whether the observation are direct and of good quality (y) or inferred ant therefore less reliable (n); data from the original studies                  
# "permanentbisexual"         whether groups in this population usually consist of both females and males year-round; data from the original studies                   
# "mostlyfdom"                whether females can dominate males in this population (f if females win more than 50% of intersexual fights) or whether males are dominant over females (m); data from the original studies      
# "strictfdom"                whether dominants in this population are usually females (f if females win more than 90% of intersexual fights) or males (m if males win more than 90% of intersexual fights) or whether there is no clear sex-bias in aggressive interactions between the sexes (c if neither sex wins more than 90% of fights); data from the original studies              
# "fem_male_dominate"         whether both females and males can dominate over individuals of the opposite sex (yy if there is at least one individual of each sex that wins fights against at least one individual of the opposite sex), whether only females (yn) or only males (ny) can dominate individuals of the opposite sex, or whether there is no clear aggressive hierarchy between the sexes (nn); data from the original studies                   
# "perc_won_females"         percentage of observed intersexual aggressive interactions that have been won by females; data from the original studies                    
# "perc_initiated_fem"      percentage of observed intersexual dominance interactions (including those that lead to submissive journals or are undecided) that have been initiated by females; data from the original studies
# "perc_males_dominatedbyfemales"           percentage of all males in the group over which the average female in the group is dominant; data from the original studies     
# "sexratio"                  proportion of adult males of the total group size; data from the original studies                   
# "experimental"              whether the observations where taken during an experimental intervention; data from the original studies                   
# "perc_aggression_fm"   percentage of all observed aggression that occurred between females and males; data from the original studies 
# "perc_aggression_ff"   percentage of all observed aggression that occurred between females; data from the original studies 
# "perc_aggression_mm"   percentage of all observed aggression that occurred between males; data from the original studies 
# "SexualDimorphism_MaleWeight_over_FemaleWeight"   sexual dimorphism in body mass in this species, calculated using two step approach after Smith 1999: M/F if males are larger; 2-F/M if females are larger. References are:
#Jarman P. (1983). Mating system and sexual dimorphism in large terrestrial, mammalian herbivores. Biological Reviews, 58(4), 485-520.
#Loison, A., Gaillard, J. M., Pélabon, C., & Yoccoz, N. G. (1999). What factors shape sexual size dimorphism in ungulates?. Evolutionary Ecology Research, 1(5), 611-633.
#Smith, R. J., & Cheverud, J. M. (2002). Scaling of sexual dimorphism in body mass: a phylogenetic analysis of Rensch's rule in primates. International Journal of Primatology, 23(5), 1095-1135.
#Isaac, J. L. (2005). Potential causes and life‐history consequences of sexual size dimorphism in mammals. Mammal Review, 35(1), 101-115.
#Kappeler, P. M., Nunn, C. L., Vining, A. Q., & Goodman, S. M. (2019). Evolutionary dynamics of sexual size dimorphism in non-volant mammals following their independent colonization of Madagascar. Scientific reports, 9(1), 1-14.
# Heldstab, S. A., van Schaik, C. P., Müller, D. W., Rensch, E., Lackey, L. B., Zerbe, P., ... & Matsuda, I. (2021). Reproductive seasonality in primates: patterns, concepts and unsolved questions. Biological reviews, 96(1), 66-88.
# "SocOrgPMK"                  social organisation of this species (G = group living, P = pair living, S = solitary living); data are from Kappeler, P. M. & Pozzi, L. (2019)                 
# "CanineDimorphism"        sexual dimorphism in canine size calculated using two step approach after Smith 1999: M/F if males are larger; 2-F/M if females are larger. Main references are:
# Lüpold, S., Simmons, L. W., & Grueter, C. C. (2019). Sexual ornaments but not weapons trade off against testes size in primates. Proceedings of the Royal Society B, 286(1900), 20182542.
# Plavcan, J. M., van Schaik, C. P., & Kappeler, P. M. (1995). Competition, coalitions and canine size in primates. Journal of human evolution, 28(3), 245-276.
# "MatSysPMK"                   mating system of this species (POL = polygynous, PRO = promiscuous, MON = monogamous, PAN = polyandrous);  data from 
# Kappeler, P. M., Nunn, C. L., Vining, A. Q., & Goodman, S. M. (2019). Evolutionary dynamics of sexual size dimorphism in non-volant mammals following their independent colonization of Madagascar. Scientific reports, 9(1), 1-14.
# Lüpold, S., Simmons, L. W., & Grueter, C. C. (2019). Sexual ornaments but not weapons trade off against testes size in primates. Proceedings of the Royal Society B, 286(1900), 20182542.
# "female_average_relatedness"   levels of average relatedness among adult female group members (r~0 unrelated, r=0.5 mother-daughter/full-siblings) reported for this species. Data are from
# Lukas, D., & Clutton‐Brock, T. (2018). Social complexity and kinship in animal societies. Ecology letters, 21(8), 1129-1134.
# "female_dispersal"           whether in this species females leave their natal group to breed in another group (Yes)  or whether they remain in their natal group to breed (No); data are from Barsbai, T., Lukas, D., & Pondorfer, A. (2021). Local convergence of behavior across species. Science, 371(6526), 292-295.
# "male_dispersal"           whether in this species males leave their natal group to breed in another group (Yes)  or whether they remain in their natal group to breed (No); data are from Barsbai, T., Lukas, D., & Pondorfer, A. (2021). Local convergence of behavior across species. Science, 371(6526), 292-295.
# "sexbias_dispersal"           whether in this species only females leave their natal group to breed in another group (Female), only males leave their natal group to breed in another group (Male), or whether most individuals of both sexes leave to breed in another group (Both); data are from Barsbai, T., Lukas, D., & Pondorfer, A. (2021). Local convergence of behavior across species. Science, 371(6526), 292-295.
# "jointaggression_females"       whether aggressive interactions sometimes involve female coalitions. Data are from: Lukas, D., & Clutton‐Brock, T. (2018). Social complexity and kinship in animal societies. Ecology letters, 21(8), 1129-1134.
# "jointaggression_males"         whether aggressive interactions sometimes involve male coalitions. Data are from Lukas, D., & Clutton‐Brock, T. (2018). Social complexity and kinship in animal societies. Ecology letters, 21(8), 1129-1134.
# "Synchrony"                     the extent to which females synchronize their mating periods (0 = not at all, 1 = completely); data are from Burtschell et al. in press                
# "r_seasonality_value"           the extent to which the environment is seasonal (0 = similar conditions year round, 1 = clearly distinct seasons); data are from Burtschell et al. in press             
# "monopolization"              whether males limit access to females during their receptive periods
# "testes_mass"                 testes mass in gram for this species; Lukas, D., & Huchard, E. (2014). The evolution of infanticide by males in mammalian societies. Science, 346(6211), 841-844.                 
# "body_mass"                   body mass in gram for this species; Lukas, D., & Huchard, E. (2014). The evolution of infanticide by males in mammalian societies. Science, 346(6211), 841-844.
# "male_skew"                   proportion of offspring sired by dominant male in one breeding season
# Kutsukake & Nunn
# Miller, C. M., Snyder-Mackler, N., Nguyen, N., Fashing, P. J., Tung, J., Wroblewski, E. E., ... & Wilson, M. L. (2021). Extragroup paternity in gelada monkeys, Theropithecus gelada, at Guassa, Ethiopia and a comparison with other primates. Animal Behaviour, 177, 277-301.
# "cooperative_breeder"         whether species is a cooperative breeder or not; Lukas, D., & Clutton-Brock, T. (2012). Cooperative breeding and monogamy in mammalian societies. Proceedings of the Royal Society B: Biological Sciences, 279(1736), 2151-2156.
# "env_harshness"               average degree of climatic harshness across the range of the species; data are from Botero, C. A., Dor, R., McCain, C. M., & Safran, R. J. (2014). Environmental harshness is positively correlated with intraspecific divergence in mammals and birds. Molecular ecology, 23(2), 259-268.
# "rainfall_unpredictability"  average degree of unpredictability in rainfall patterns across years for that species; data are from Botero, C. A., Dor, R., McCain, C. M., & Safran, R. J. (2014). Environmental harshness is positively correlated with intraspecific divergence in mammals and birds. Molecular ecology, 23(2), 259-268.
# "female_canine_height"  height of the canine of females in cm; data are from Lüpold S, Simmons LW, Grueter CC. 2019. Sexual ornaments but not weapons trade off against testes size in primates. Proc. R. Soc. B 20182542. (doi:10.1098/rspb.2018.2542)
# "male_canine_height"  height of the canine of males in cm; data are from Lüpold S, Simmons LW, Grueter CC. 2019. Sexual ornaments but not weapons trade off against testes size in primates. Proc. R. Soc. B 20182542. (doi:10.1098/rspb.2018.2542)
# "female_evictions"  whether females forcibly evict other females from their social groups. Data are from Lukas, D., & Huchard, E. (2019). The evolution of infanticide by females in mammals. Philosophical Transactions of the Royal Society B, 374(1780), 20180075.
# "homerange_overlap" the percentage to which the home ranges of neighboring home ranges of solitary individuals/groups overlap; data are from Pearce, F., Carbone, C., Cowlishaw, G., & Isaac, N. J. (2013). Space-use scaling and home range overlap in primates. Proceedings of the Royal Society B: Biological Sciences, 280(1751), 20122122.
# "female_infanticide"  whether females kill the dependent offspring of other females. Data are from Lukas, D., & Huchard, E. (2019). The evolution of infanticide by females in mammals. Philosophical Transactions of the Royal Society B, 374(1780), 20180075.
# "receptive_synchrony": the probability that two or more females are receptive/mating on the same day.   Data are from: Carnes, L. M., Nunn, C. L., & Lewis, R. J. (2011). Effects of the distribution of female primates on the number of males. PLoS One, 6(5), e19853.; and Gogarten, J. F., & Koenig, A. (2013). Reproductive seasonality is a poor predictor of receptive synchrony and male reproductive skew among nonhuman primates. Behavioral Ecology and Sociobiology, 67(1), 123-134.
# "sexualornamentdimorphism_sum": the number of body parts in which males differ in their ornamentation from females. Data are from: 
# Lüpold, S., Simmons, L. W., & Grueter, C. C. (2019). Sexual ornaments but not weapons trade off against testes size in primates. Proceedings of the Royal Society B, 286(1900), 20182542.
# "ovulation_signs": whether females have no (0), some (0.5), slight (1), or clear (2) signs of ovulation. Data are from: 
# Rooker, K., & Gavrilets, S. (2018). On the evolution of visual female sexual signalling. Proceedings of the Royal Society B: Biological Sciences, 285(1879), 20172875.
# "sexualreceptivity_hours": the number of hours females are receptive each month. Data are from:
# Kutsukake, N., & Nunn, C. L. (2006). Comparative tests of reproductive skew in male primates: the roles of demographic factors and incomplete control. Behavioral Ecology and Sociobiology, 60(5), 695-706.
# Stockley, P. (2002). Sperm competition risk and male genital anatomy: comparative evidence for reduced duration of female sexual receptivity in primates with penile spines. Evolutionary Ecology, 16(2), 123-137.
# "Strata_Wilman": whether species are primarily found on the ground (g), below ground (s), or in trees (ar). Data are from: Wilman, H., Belmaker, J., Simpson, J., de la Rosa, C., Rivadeneira, M. M., & Jetz, W. (2014). EltonTraits 1.0: Species‐level foraging attributes of the world's birds and mammals: Ecological Archives E095‐178. Ecology, 95(7), 2027-2027.
# "StrictlyNocturnal_Wilman": whether species are strictly nocturnal (1) or also active during the day (0). Data are from: Wilman, H., Belmaker, J., Simpson, J., de la Rosa, C., Rivadeneira, M. M., & Jetz, W. (2014). EltonTraits 1.0: Species‐level foraging attributes of the world's birds and mammals: Ecological Archives E095‐178. Ecology, 95(7), 2027-2027.