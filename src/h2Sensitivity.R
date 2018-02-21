source("quollTGFFunctions.R")
source("parameters.R")

# test a range of nIntr values for TGF introduction the year that toads arrive
nI<-seq(from=0,to=50,by=5) #values of nIntr
yr<-seq(21,41,1) # year of introduction
n<-100 # number of iterations (test)

#output
e <-vector() 
g <-vector()
temp<-vector()
out<-matrix(ncol=4)

for(k in yr) {
  pr<-vector()
  prop.g<-vector()
  for(j in nI) {
    for (i in 1:n) {
      test<-TGF(demPars = pars, nP = 10, nN = 10, W0 = 0.09, s0=0.38, h2=0.1, burnIn=10, noSel=20, selGens=50, TGFtime=k, W1=0.9, nIntr=j, outbreeding=FALSE, nl=10, H2=1, alpha=1)    
      e[i]<-ifelse(tail(test$n,1)<=1, TRUE, FALSE)
      g[i]<-tail(test$genome, 1)
    }
    pr<- c(pr, length(e[e==TRUE])/n)
    prop.g<-c(prop.g, mean(g))
  }
  temp<-cbind(rep(k, length(nI)), nI, pr, prop.g)
  out<- rbind(out, temp)
}
out<-out[2:232,]

save(out, file="../out/h2=0.1.RData")
