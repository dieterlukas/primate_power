
aggression<-combined
aggression$females<-exp(aggression$females)
aggression$males<-exp(aggression$males)
aggression$ff_dyads<-(aggression$females*(aggression$females-1))/2
aggression$mm_dyads<-(aggression$males*(aggression$males-1))/2
aggression$fm_dyads<-aggression$males*aggression$females
aggression$perc_fm_dyads<-aggression$fm_dyads*100/(aggression$fm_dyads+aggression$ff_dyads+aggression$mm_dyads)
aggression$perc_ff_dyads<-aggression$ff_dyads*100/(aggression$fm_dyads+aggression$ff_dyads+aggression$mm_dyads)
aggression$perc_mm_dyads<-aggression$mm_dyads*100/(aggression$fm_dyads+aggression$ff_dyads+aggression$mm_dyads)


plot(aggression$perc_aggression_fm~aggression$perc_fm_dyads,pch=20,col="#009E73",bty="n",cex=2,xlim=c(0,100),ylim=c(0,100),xlab="Percent of dyads in group that are between females and males",ylab="Percent of aggression in group that is between females and males")
lines(x=c(10,70),y=c(48,47),col="#009E73",lwd=4)

points(aggression$perc_aggression_ff~aggression$perc_ff_dyads,pch=20,col="#CC79A7",cex=2)
lines(x=c(10,90),y=c(30,40),col="#CC79A7",lwd=4)

points(aggression$perc_aggression_mm~aggression$perc_mm_dyads,pch=20,col="#56B4E9",cex=2)
lines(x=c(0,40),y=c(16,22),col="#56B4E9",lwd=4)
