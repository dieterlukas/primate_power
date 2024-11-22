###### Take the complete results and split the output by the different hypotheses, formatting them according to whether a given prediction is supported or not



combinedresults<-read.csv("../results/intersexualdominance_combinedresults.csv")

combinedresults<-combinedresults[combinedresults$phylogeny=="Yes",]

combinedresults$estimate_lower<-round(combinedresults$estimate_lower,3)
combinedresults$estimate_upper<-round(combinedresults$estimate_upper,3)


r_percwon<-combinedresults[combinedresults$outcome %in% "perc_won",]
r_strictdominance<-combinedresults[combinedresults$outcome %in% combinedresults$outcome[2],]
r_strictfemale<-combinedresults[combinedresults$outcome %in% combinedresults$outcome[3],]
r_strictmale<-combinedresults[combinedresults$outcome %in% combinedresults$outcome[4],]



femalecompetition<-c(discrete_femalecompetition,continuous_femalecompetition)

r_femalecompetition_percwon<-r_percwon[unique(grep(paste(femalecompetition,collapse="|"),r_percwon$continuouspredictor)),]
r_femalecompetition_strictdominance<-r_strictdominance[unique(grep(paste(femalecompetition,collapse="|"),r_strictdominance$continuouspredictor)),]
r_femalecompetition_strictfemale<-r_strictfemale[unique(grep(paste(femalecompetition,collapse="|"),r_strictfemale$continuouspredictor)),]
r_femalecompetition_strictmale<-r_strictmale[unique(grep(paste(femalecompetition,collapse="|"),r_strictmale$continuouspredictor)),]

write.csv(r_femalecompetition_percwon,file="../results/femalecompetition_percwon_results.csv")
write.csv(r_femalecompetition_strictdominance,file="../results/femalecompetition_strictthreeway_results.csv")
write.csv(r_femalecompetition_strictfemale,file="../results/femalecompetition_strictfemale_results.csv")
write.csv(r_femalecompetition_strictmale,file="../results/femalecompetition_strictmale_results.csv")

femalecompetition_comparison<-combinedresults[unique(grep(paste(femalecompetition,collapse="|"),combinedresults$continuouspredictor)),]
femalecompetition_comparison<-femalecompetition_comparison[,c(1,2,6)]
femalecompetition_comparison<-femalecompetition_comparison %>% pivot_wider(names_from = outcome,values_from = association_present)
femalecompetition_comparison<-as.data.frame(femalecompetition_comparison)
write.csv(femalecompetition_comparison,file="../results/femalecompetition_resultcomparison.csv")




sexualselection<-c(discrete_sexualselection,continuous_sexualselection)

r_sexualselection_percwon<-r_percwon[unique(grep(paste(sexualselection,collapse="|"),r_percwon$continuouspredictor)),]
r_sexualselection_strictdominance<-r_strictdominance[unique(grep(paste(sexualselection,collapse="|"),r_strictdominance$continuouspredictor)),]
r_sexualselection_strictfemale<-r_strictfemale[unique(grep(paste(sexualselection,collapse="|"),r_strictfemale$continuouspredictor)),]
r_sexualselection_strictmale<-r_strictmale[unique(grep(paste(sexualselection,collapse="|"),r_strictmale$continuouspredictor)),]

write.csv(r_sexualselection_percwon,file="../results/sexualselection_percwon_results.csv")
write.csv(r_sexualselection_strictdominance,file="../results/sexualselection_strictthreeway_results.csv")
write.csv(r_sexualselection_strictfemale,file="../results/sexualselection_strictfemale_results.csv")
write.csv(r_sexualselection_strictmale,file="../results/sexualselection_strictmale_results.csv")

sexualselection_comparison<-combinedresults[unique(grep(paste(sexualselection,collapse="|"),combinedresults$continuouspredictor)),]
sexualselection_comparison<-sexualselection_comparison[,c(1,2,6)]
sexualselection_comparison<-sexualselection_comparison %>% pivot_wider(names_from = outcome,values_from = association_present,values_fn = {unique})
sexualselection_comparison<-as.data.frame(sexualselection_comparison)
write.csv(sexualselection_comparison,file="../results/sexualselection_resultcomparison.csv")



femalesociality<-c(discrete_femalesociality,continuous_femalesociality)

r_femalesociality_percwon<-r_percwon[unique(grep(paste(femalesociality,collapse="|"),r_percwon$continuouspredictor)),]
r_femalesociality_strictdominance<-r_strictdominance[unique(grep(paste(femalesociality,collapse="|"),r_strictdominance$continuouspredictor)),]
r_femalesociality_strictfemale<-r_strictfemale[unique(grep(paste(femalesociality,collapse="|"),r_strictfemale$continuouspredictor)),]
r_femalesociality_strictmale<-r_strictmale[unique(grep(paste(femalesociality,collapse="|"),r_strictmale$continuouspredictor)),]

write.csv(r_femalesociality_percwon,file="../results/femalesociality_percwon_results.csv")
write.csv(r_femalesociality_strictdominance,file="../results/femalesociality_strictthreeway_results.csv")
write.csv(r_femalesociality_strictfemale,file="../results/femalesociality_strictfemale_results.csv")
write.csv(r_femalesociality_strictmale,file="../results/femalesociality_strictmale_results.csv")

femalesociality_comparison<-combinedresults[unique(grep(paste(femalesociality,collapse="|"),combinedresults$continuouspredictor)),]
femalesociality_comparison<-femalesociality_comparison[,c(1,2,6)]
femalesociality_comparison<-femalesociality_comparison %>% pivot_wider(names_from = outcome,values_from = association_present)
femalesociality_comparison<-as.data.frame(femalesociality_comparison)
write.csv(femalesociality_comparison,file="../results/femalesociality_resultcomparison.csv")



offspringloss<-c(discrete_offspringloss,continuous_offspringloss)

r_offspringloss_percwon<-r_percwon[unique(grep(paste(offspringloss,collapse="|"),r_percwon$continuouspredictor)),]
r_offspringloss_strictdominance<-r_strictdominance[unique(grep(paste(offspringloss,collapse="|"),r_strictdominance$continuouspredictor)),]
r_offspringloss_strictfemale<-r_strictfemale[unique(grep(paste(offspringloss,collapse="|"),r_strictfemale$continuouspredictor)),]
r_offspringloss_strictmale<-r_strictmale[unique(grep(paste(offspringloss,collapse="|"),r_strictmale$continuouspredictor)),]

write.csv(r_offspringloss_percwon,file="../results/offspringloss_percwon_results.csv")
write.csv(r_offspringloss_strictdominance,file="../results/offspringloss_strictthreeway_results.csv")
write.csv(r_offspringloss_strictfemale,file="../results/offspringloss_strictfemale_results.csv")
write.csv(r_offspringloss_strictmale,file="../results/offspringloss_strictmale_results.csv")

offspringloss_comparison<-combinedresults[unique(grep(paste(offspringloss,collapse="|"),combinedresults$continuouspredictor)),]
offspringloss_comparison<-offspringloss_comparison[,c(1,2,6)]
offspringloss_comparison<-offspringloss_comparison %>% pivot_wider(names_from = outcome,values_from = association_present)
offspringloss_comparison<-as.data.frame(offspringloss_comparison)
write.csv(offspringloss_comparison,file="../results/offspringloss_resultcomparison.csv")




selforganisation<-c(continuous_selforganisation)

r_selforganisation_percwon<-r_percwon[unique(grep(paste(selforganisation,collapse="|"),r_percwon$continuouspredictor)),]
r_selforganisation_strictdominance<-r_strictdominance[unique(grep(paste(selforganisation,collapse="|"),r_strictdominance$continuouspredictor)),]
r_selforganisation_strictfemale<-r_strictfemale[unique(grep(paste(selforganisation,collapse="|"),r_strictfemale$continuouspredictor)),]
r_selforganisation_strictmale<-r_strictmale[unique(grep(paste(selforganisation,collapse="|"),r_strictmale$continuouspredictor)),]

write.csv(r_selforganisation_percwon,file="../results/selforganisation_percwon_results.csv")
write.csv(r_selforganisation_strictdominance,file="../results/selforganisation_strictthreeway_results.csv")
write.csv(r_selforganisation_strictfemale,file="../results/selforganisation_strictfemale_results.csv")
write.csv(r_selforganisation_strictmale,file="../results/selforganisation_strictmale_results.csv")

selforganisation_comparison<-combinedresults[unique(grep(paste(selforganisation,collapse="|"),combinedresults$continuouspredictor)),]
selforganisation_comparison<-selforganisation_comparison[,c(1,2,6)]
selforganisation_comparison<-selforganisation_comparison %>% pivot_wider(names_from = outcome,values_from = association_present,values_fn = {unique})
selforganisation_comparison<-as.data.frame(selforganisation_comparison)
write.csv(selforganisation_comparison,file="../results/selforganisation_resultcomparison.csv")


