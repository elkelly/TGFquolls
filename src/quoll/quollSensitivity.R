source("quollTGFFunctions.R")
source("parameters.R")

# test a range of nIntr values for TGF introduction the year that toads arrive
nI<-seq(from=0,to=50,by=5) #values of nIntr
yr<-seq(21,41,1) # year of introduction
n<-100 # number of iterations (test)

##################
## h2 and K Sensisity ##
##################

## High ##
K<- 1000
h2<- 0.3
W0<- 0.03

#output
e <-vector() 
g <-vector()
g2<-vector()
popFreq<-matrix(ncol=10)
Ip<-vector()
temp<-vector()
out<-matrix(ncol=6)

for(k in yr) {
  pr<-vector()
  r.current<-vector()
  r.new<-vector()
  Ipp<-vector()
  for(j in nI) {
    for (i in 1:n) {
      test<-TGF(demPars = pars, nP = 10, nN = 10, nI = 10, W0 = W0, f0=0.05,K = K, s0=0.38, h2=h2, burnIn=10, noSel=20, selGens=50, TGFtime=k, f1=0.9, nIntr=j, outbreeding=FALSE, nl=10, H2=1, alpha=1, recombRate=50)    
      e[i]<-ifelse(tail(test$n,1)<=1, TRUE, FALSE)
      g[i]<-tail(test$genome[,1], 1) #for whole pop, pRecI = [,2]
      g2[i]<-tail(test$genome[,2], 1)
      popFreq<<-rbind(popFreq, as.numeric(tail(test$genome[,3], 1)[[1]]))
      Ip[i]<-ifelse(is.null(test$IntProp), NA, test$IntProp)
    }
    pr<- c(pr, length(e[e==TRUE])/n)
    r.current<-c(r.current, mean(as.numeric(g)))
    r.new<-c(r.new, mean(as.numeric(g2)))
    Ipp<-c(Ipp, mean(Ip, na.rm=TRUE))
  }
  temp<-cbind(rep(k, length(nI)), nI, pr, r.current, r.new, Ipp)
  out<- rbind(out, temp)
}
out<-out[2:nrow(out),]

save(out, file="../../out/quollh20.3K1000.RData")

##################

## LowK ##
K<- 125
h2<- 0.3
W0<- 0.08

#output
e <-vector() 
g <-vector()
g2<-vector()
popFreq<-matrix(ncol=10)
Ip<-vector()
temp<-vector()
out<-matrix(ncol=6)

for(k in yr) {
  pr<-vector()
  r.current<-vector()
  r.new<-vector()
  Ipp<-vector()
  for(j in nI) {
    for (i in 1:n) {
      test<-TGF(demPars = pars, nP = 10, nN = 10, nI = 10, W0 = W0, f0=0.05,K = K, s0=0.38, h2=h2, burnIn=10, noSel=20, selGens=50, TGFtime=k, f1=0.9, nIntr=j, outbreeding=FALSE, nl=10, H2=1, alpha=1, recombRate=50)    
      e[i]<-ifelse(tail(test$n,1)<=1, TRUE, FALSE)
      g[i]<-tail(test$genome[,1], 1) #for whole pop, pRecI = [,2]
      g2[i]<-tail(test$genome[,2], 1)
      popFreq<<-rbind(popFreq, as.numeric(tail(test$genome[,3], 1)[[1]]))
      Ip[i]<-ifelse(is.null(test$IntProp), NA, test$IntProp)
    }
    pr<- c(pr, length(e[e==TRUE])/n)
    r.current<-c(r.current, mean(as.numeric(g)))
    r.new<-c(r.new, mean(as.numeric(g2)))
    Ipp<-c(Ipp, mean(Ip, na.rm=TRUE))
  }
  temp<-cbind(rep(k, length(nI)), nI, pr, r.current, r.new, Ipp)
  out<- rbind(out, temp)
}
out<-out[2:nrow(out),]

save(out, file="../../out/quollh20.3K125.RData")

##################

## Lowh2 ##
K<- 1000
h2<- 0.1
W0<- 0.15

e <-vector() 
g <-vector()
g2<-vector()
popFreq<-matrix(ncol=10)
Ip<-vector()
temp<-vector()
out<-matrix(ncol=6)

for(k in yr) {
  pr<-vector()
  r.current<-vector()
  r.new<-vector()
  Ipp<-vector()
  for(j in nI) {
    for (i in 1:n) {
      test<-TGF(demPars = pars, nP = 10, nN = 10, nI = 10, W0 = W0, f0=0.05,K = K, s0=0.38, h2=h2, burnIn=10, noSel=20, selGens=50, TGFtime=k, f1=0.9, nIntr=j, outbreeding=FALSE, nl=10, H2=1, alpha=1, recombRate=50)    
      e[i]<-ifelse(tail(test$n,1)<=1, TRUE, FALSE)
      g[i]<-tail(test$genome[,1], 1) #for whole pop, pRecI = [,2]
      g2[i]<-tail(test$genome[,2], 1)
      popFreq<<-rbind(popFreq, as.numeric(tail(test$genome[,3], 1)[[1]]))
      Ip[i]<-ifelse(is.null(test$IntProp), NA, test$IntProp)
    }
    pr<- c(pr, length(e[e==TRUE])/n)
    r.current<-c(r.current, mean(as.numeric(g)))
    r.new<-c(r.new, mean(as.numeric(g2)))
    Ipp<-c(Ipp, mean(Ip, na.rm=TRUE))
  }
  temp<-cbind(rep(k, length(nI)), nI, pr, r.current, r.new, Ipp)
  out<- rbind(out, temp)
}
out<-out[2:nrow(out),]

save(out, file="../../out/quollh20.1K1000.RData")

##################

## Low ##
K<- 125
h2<- 0.1
W0<- 0.3

#output
e <-vector() 
g <-vector()
g2<-vector()
popFreq<-matrix(ncol=10)
Ip<-vector()
temp<-vector()

out<-matrix(ncol=6)

for(k in yr) {
  pr<-vector()
  r.current<-vector()
  r.new<-vector()
  Ipp<-vector()
  for(j in nI) {
    for (i in 1:n) {
      test<-TGF(demPars = pars, nP = 10, nN = 10, nI = 10, W0 = W0, f0=0.05,K = K, s0=0.38, h2=h2, burnIn=10, noSel=20, selGens=50, TGFtime=k, f1=0.9, nIntr=j, outbreeding=FALSE, nl=10, H2=1, alpha=1, recombRate=50)    
      e[i]<-ifelse(tail(test$n,1)<=1, TRUE, FALSE)
      g[i]<-tail(test$genome[,1], 1) #for whole pop, pRecI = [,2]
      g2[i]<-tail(test$genome[,2], 1)
      popFreq<<-rbind(popFreq, as.numeric(tail(test$genome[,3], 1)[[1]]))
      Ip[i]<-ifelse(is.null(test$IntProp), NA, test$IntProp)
    }
    pr<- c(pr, length(e[e==TRUE])/n)
    r.current<-c(r.current, mean(as.numeric(g)))
    r.new<-c(r.new, mean(as.numeric(g2)))
    Ipp<-c(Ipp, mean(Ip, na.rm=TRUE))
  }
  temp<-cbind(rep(k, length(nI)), nI, pr, r.current, r.new, Ipp)
  out<- rbind(out, temp)
}
out<-out[2:nrow(out),]

save(out, file="../../out/quollh20.1K125.RData")

##################

##################
## Recomb rate ##
##################

## High ##
R <- 1.329258 # recomb rate of 0.25
W0<- 0.08
n<-100

  #output
  e <-vector() 
g <-vector()
g2<-vector()
popFreq<-matrix(ncol=10)
Ip<-vector()
temp<-vector()
out<-matrix(ncol=6)

for(k in yr) {
  pr<-vector()
  r.current<-vector()
  r.new<-vector()
  Ipp<-vector()
  for(j in nI) {
    for (i in 1:n) {
      test<-TGF(demPars = pars, nP = 10, nN = 10, nI = 10, W0 = W0, f0=0.05,K = 500, s0=0.38, h2=0.2, burnIn=10, noSel=20, selGens=50, TGFtime=k, f1=0.9, nIntr=j, outbreeding=FALSE, nl=10, H2=1, alpha=1, recombRate=R)    
      e[i]<-ifelse(tail(test$n,1)<=1, TRUE, FALSE)
      g[i]<-tail(test$genome[,1], 1) #for whole pop, pRecI = [,2]
      g2[i]<-tail(test$genome[,2], 1)
      popFreq<<-rbind(popFreq, as.numeric(tail(test$genome[,3], 1)[[1]]))
      Ip[i]<-ifelse(is.null(test$IntProp), NA, test$IntProp)
    }
    pr<- c(pr, length(e[e==TRUE])/n)
    r.current<-c(r.current, mean(as.numeric(g)))
    r.new<-c(r.new, mean(as.numeric(g2)))
    Ipp<-c(Ipp, mean(Ip, na.rm=TRUE))
  }
  temp<-cbind(rep(k, length(nI)), nI, pr, r.current, r.new, Ipp)
  out<- rbind(out, temp)
}
out<-out[2:nrow(out),]

save(out, file="../../out/quollrecomb.RData")
