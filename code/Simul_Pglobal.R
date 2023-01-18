


dat <- read.csv2("data/dataPop.csv", header=TRUE, stringsAsFactors = FALSE, comment.char = "#")
#attach(dat)
#pops <- c("Leff", "Trieux", "Jaudy", "Leguer", "Yar", "Douron", "Penze", "Elorn", "Aulne", "Goyen", "Odet", "Aven", "Laita", "Scorff", "Blavet")
#dat <- dat[which(dat$Population %in% pops),] # remove Couesnon
dat <- dat[-which(dat$Population == "Couesnon"),] # remove Couesnon


prop = dat$Area/sum(dat$Area)
P=array(,dim=c(length(vlocal),3))

vlocal=c(0, 0.1, 0.2, 0.3, 0.4, 0.5)
for (j in 1:length(vlocal)){

plocal=plocal1=plocal2=NULL
for (i in 1:npop){
  plocal[i] <-  vlocal[j]
  plocal1[i] <- ifelse(Type[i]=="sink", 0, vlocal[j])
  plocal2[i] <- ifelse(Type[i]=="source", 0, vlocal[j])
}

P[j,1] = sum(prop * plocal)
P[j,2] = sum(prop * plocal1)
P[j,3] = sum(prop * plocal2)

}


plot(P[,1],vlocal,type='l')
lines(P[,2],vlocal,col=2)
lines(P[,3],vlocal,col=3)




prop_rand <- sample(prop)
P=array(,dim=c(length(vlocal),3))

vlocal=c(0, 0.1, 0.2, 0.3, 0.4, 0.5)
for (j in 1:length(vlocal)){
  
  plocal=plocal1=plocal2=NULL
  for (i in 1:npop){
    plocal[i] <-  vlocal[j]
    plocal1[i] <- ifelse(Type[i]=="sink", 0, vlocal[j])
    plocal2[i] <- ifelse(Type[i]=="source", 0, vlocal[j])
  }
  
  P[j,1] = sum(prop_rand * plocal)
  P[j,2] = sum(prop_rand * plocal1)
  P[j,3] = sum(prop_rand * plocal2)
  
}

lines(P[,1],vlocal,type='l',lty=2)
lines(P[,2],vlocal,col=2,lty=2)
lines(P[,3],vlocal,col=3,lty=2)



Type_rand=NULL
for (i in 1:npop){
  Type_rand[i]<-sample(c("source","sink","neutral"),1)
}

P=array(,dim=c(length(vlocal),3))

vlocal=c(0, 0.1, 0.2, 0.3, 0.4, 0.5)
for (j in 1:length(vlocal)){
  
  plocal=plocal1=plocal2=NULL
  for (i in 1:npop){
    plocal[i] <-  vlocal[j]
    plocal1[i] <- ifelse(Type_rand[i]=="sink", 0, vlocal[j])
    plocal2[i] <- ifelse(Type_rand[i]=="source", 0, vlocal[j])
  }
  
  P[j,1] = sum(prop_rand * plocal)
  P[j,2] = sum(prop_rand * plocal1)
  P[j,3] = sum(prop_rand * plocal2)
  
}

lines(P[,1],vlocal,type='l',lty=3)
lines(P[,2],vlocal,col=2,lty=3)
lines(P[,3],vlocal,col=3,lty=3)
