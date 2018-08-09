#### test W0 probability of extinction with no TGF #######
source("quollTGFFunctions.R")
source("parameters.R")


# establish W0 for pr=0.95 for different versions of K and h2 
W0<-c(0.01, 0.03, 0.05, 0.08, 0.1, 0.2, 0.3, 0.4) # broad sweap of w0 values to catch drop in pe for each scenario
K<-c(125, 1000)
h2<-c(0.1, 0.3)
n<-100

e <-vector() 
temp<-vector()
out1<- vector()

for(j in W0) {
  pr<-vector()
    for(k in K) {
      for(l in h2) {
        print(c(j, l, k))
        for (i in 1:n) {
          test<-TGF(demPars = pars, nP = 10, nN = 10, nI = 10, W0 = j, f0=0.05, K = k, s0=0.38, h2=l, burnIn=10, noSel=20, selGens=50, TGFtime=30, f1=0.9, nIntr=0, outbreeding=FALSE, nl=10, H2=1, alpha=1, recombRate=50)    
          e[i]<-ifelse(tail(test$n,1)<=1, TRUE, FALSE)
        }
        pr<- c(pr, length(e[e==TRUE])/n)
      }
    }
  out1<- c(out1, pr)
}

## standard model 
e <-vector() 
temp<-vector()
out2<- vector()

for(j in W0) {
  pr<-vector()
      for (i in 1:n) {
        test<-TGF(demPars = pars, nP = 10, nN = 10, nI = 10, W0 = j, f0=0.05, K = 500, s0=0.38, h2=0.2, burnIn=10, noSel=20, selGens=50, TGFtime=30, f1=0.9, nIntr=0, outbreeding=FALSE, nl=10, H2=1, alpha=1, recombRate=50)    
        e[i]<-ifelse(tail(test$n,1)<=1, TRUE, FALSE)
      }
      pr<- c(pr, length(e[e==TRUE])/n)
      out2<- c(out2, pr)
}

# edited recomb rate model 
W0<-c(0.05, 0.08, 0.1) # broad sweap of w0 values to catch drop in pe for each scenario
n<-100

e <-vector() 
temp<-vector()
out3<- vector()

for(j in W0) {
  pr<-vector()
  for (i in 1:n) {
    test<-TGF(demPars = pars, nP = 10, nN = 10, nI = 10, W0 = j, f0=0.05, K = 500, s0=0.38, h2=0.2, burnIn=10, noSel=20, selGens=50, TGFtime=30, f1=0.9, nIntr=0, outbreeding=FALSE, nl=10, H2=1, alpha=1, recombRate=25)    
    e[i]<-ifelse(tail(test$n,1)<=1, TRUE, FALSE)
  }
  pr<- c(pr, length(e[e==TRUE])/n)
  out3<- c(out3, pr)
}

