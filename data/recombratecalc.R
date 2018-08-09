
pop<-initInds(100,10, 10, 10, f=0.05)# initalise a pop

x<-c()
y<-c()
cuts<-c(0:60) #cuts 0-60

compair<-function(colDam){
  o<-outer(colDam,colDam,"==") # pairwise comparison, matrix of does x[1]=x[2] etc etc 
  d<-diag(length(colDam)) #creates matrix with correct dimentions and 1s in a diaginal
  d<-ifelse(d==1, NA, o) # turns 1s to NAs for diaginal 
  d<-ifelse(d==FALSE, 1, 0) # all falses are 1s (recombination between pairs), everything else is 0
  mean(d, na.rm=T) # mean of that
}

for(j in cuts){
  for(i in 1:100){
    lLoci<- runif(30,0,1) #assign random loci locations
    nCuts<-rpois(nrow(pop), j) #randomly assign number of cuts
    lCut<-lapply(nCuts, runif, min=0, max=1) # position cuts
    lCut<-mapply(sort, lCut) # sort cuts by position
    damSet<-mapply(findInterval, lCut, MoreArgs = list(x = lLoci)) %% 2 # sort based on cut locations
    x[i]<-mean(apply(damSet, 2, compair)) # use compair function to work out recomb rate
    # 
  }  
  y<-c(y, mean(x))
}

plot(y~cuts, type="l")

approx(y, cuts, 0.25)