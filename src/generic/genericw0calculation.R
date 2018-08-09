#### test W0 probability of extinction with no TGF #######
source("genericModelFunctions.R")


n<-50
W0<-c(0.001, 0.003, 0.005, 0.008, 0.01, 0.05, 0.1, 0.5)
e <-vector() 
temp<-vector()
out2<- vector()

for(k in W0) {
  pr<-vector()
  print(c(k))
  for (i in 1:n) {
    test<-TGF(nP = 10, nN = 10, nI = 10, W0 = k, f0 = 0.05, h2=0.2, Nstar=500, burnIn=10, noSel=20, selGens=50, TGFtime=30, f1=0.9, nIntr=0, outbreeding=FALSE, nl=10, H2=1, alpha=1, Rmax=3,recombRate=50)    
    e[i]<-ifelse(tail(test$n,1)<=1, TRUE, FALSE)
  }
  pr<- c(pr, length(e[e==TRUE])/n)
  out2<- c(out2, pr)
}


# establish W0 for pr=0.95 for different versions of K and h2 
# course scale 
W0<-c(0.00001, 0.00005, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.03, 0.05, 0.08, 0.1, 0.5) # broad sweap of w0 values to catch drop in pe for each scenario
R<-c(1.5, 6)
NS<-c(125, 1000)
h<-c(0.1, 0.3)
n<-50

e <-vector() 
temp<-vector()
out1<- matrix()

for(k in W0) {
  pr<-vector()
  for(j in R) {
    for(m in NS) {
      for(l in h) {
        print(c(j, m, l, k))
        for (i in 1:n) {
          test<-TGF(nP = 10, nN = 10, nI = 10, W0 = k, f0=0.05, h2=l, Nstar=m, burnIn=10, noSel=20, selGens=50, TGFtime=100, f1=0.9, nIntr=0, outbreeding=FALSE, nl=10, H2=1, alpha=1, Rmax=j,recombRate=50)    
          e[i]<-ifelse(tail(test$n,1)<=1, TRUE, FALSE)
        }
        pr<- c(pr, length(e[e==TRUE])/n)
      }
    }
  }
  out1<- c(out1, pr)
}


h<-rep(h, 2)
NS<-rep(NS, each=2)
R<-rep(R, each=4)
results<-cbind(h, NS, R)
results<-results[rep(seq(nrow(results)), length(W0)),]
W02<-rep(W0, each=length(R))

out1<-out1[2:length(out1)]
results<-cbind(results, W02, out1)


out2<-out2[2:length(out1)]
NS<-rep(500, length(W0))
h<-rep(0.2, length(W0))
R<-rep(3, length(W0))
results2<-cbind(h, NS, R, W0, out2)

w095<-rbind(results, results2)

w095<-w095[!is.na(w095[,"out1"]),]

save(w095, file="../../out/genericW0.Rdata")
