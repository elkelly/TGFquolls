# Run genetic population function with standard parameters and no outbreeding scenario of targeted gene flow
source("genericModelFunctions.R")

# test a range of nIntr values for TGF introduction the year that toads arrive
nI<-seq(from=0,to=50,by=5) #values of nIntr
yr<-seq(21,41,1) # year of introduction
n<-100 # number of iterations

#output
e <-vector() 
g <-vector()
g2<-vector()
popFreq<-matrix(ncol=10)
Ip<-vector()
temp<-vector()
out<-matrix(ncol=6)


w0<-0.005

for(k in yr) {
  pr<-vector()
  r.current<-vector()
  r.new<-vector()
  Ipp<-vector()
  for(j in nI) {
    for (i in 1:n) {
      print(c(j, k))
      test<-TGF(nP = 10, nN = 10, nI = 10, W0 = w0, f0 = 0.05, h2=0.2, Nstar=500, burnIn=10, noSel=20, selGens=50, TGFtime=k, f1=0.9, nIntr=j, outbreeding=FALSE, nl=10, H2=1, alpha=1, Rmax=3,recombRate=50)    
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

save(out, file="../../out/generic.Rdata")
