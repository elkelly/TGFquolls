# targeted gene flow on generic population model with Outbreeding Depression
## Setup for outbreeding scenarios
#comuting H2 
H2<-1
h1<-H2*.5
h0<-H2*.25

#worst case ps (all hetrozygotes)
p1<-0
p2<-0
pH<-1
nl<-10

#calclate Es
es<-nl*(p1*p2*H2+(p1+p2)*pH*h1+pH*h0)
# es = 2.5

#scenarios for alpha
alpha50<- (-log(0.5, base=exp(1))/es) # 50% reduction
alpha10<- (-log(0.9, base=exp(1))/es) # 10% reduction

#run model
source("genericModelFunctions.R")

# test a range of nIntr values for TGF introduction the year that toads arrive
nI<-seq(from=0,to=50,by=5) #values of nIntr
yr<-seq(21,41,1) # year of introduction
n<-100  # number of iterations

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
        test<-TGF(nP = 10, nN = 10, nI = 10, W0 = w0, f0 = 0.05, h2=0.2, Nstar=500, burnIn=10, noSel=20, selGens=50, TGFtime=k, f1=0.9, nIntr=j, outbreeding=TRUE, nl=10, H2=1, alpha=alpha50, Rmax=3, recombRate=50)    
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

out_alpha50<-out[2:nrow(out),]
out_alpha50<-cbind(out_alpha50, "alpha"=rep("alpha50", nrow(out_alpha50)))

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
        test<-TGF(nP = 10, nN = 10, nI = 10, W0 = w0, f0 = 0.05, h2=0.2, Nstar=500, burnIn=10, noSel=20, selGens=50, TGFtime=k, f1=0.9, nIntr=j, outbreeding=TRUE, nl=10, H2=1, alpha=alpha10, Rmax=3, recombRate=50)    
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

out_alpha10<-out[2:nrow(out),]
out_alpha10<-cbind(out_alpha10, "alpha"=rep("alpha10", nrow(out_alpha10)))


OBout<-rbind(out_alpha10, out_alpha50)
save(OBout, file="../../out/genericOutB.RData")






#recomb rate
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
      test<-TGF(nP = 10, nN = 10, nI = 10, W0 = w0, f0 = 0.05, h2=0.2, Nstar=500, burnIn=10, noSel=20, selGens=50, TGFtime=k, f1=0.9, nIntr=j, outbreeding=TRUE, nl=10, H2=1, alpha=alpha50, Rmax=3, recombRate=1.329258)    
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

out_alpha50r<-out[2:nrow(out),]
out_alpha50r<-cbind(out_alpha50r, "alpha"=rep("alpha50", nrow(out_alpha50r)))

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
      test<-TGF(nP = 10, nN = 10, nI = 10, W0 = w0, f0 = 0.05, h2=0.2, Nstar=500, burnIn=10, noSel=20, selGens=50, TGFtime=k, f1=0.9, nIntr=j, outbreeding=TRUE, nl=10, H2=1, alpha=alpha10, Rmax=3, recombRate=1.329258)    
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

out_alpha10r<-out[2:nrow(out),]
out_alpha10r<-cbind(out_alpha10r, "alpha"=rep("alpha10", nrow(out_alpha10r)))


OBoutr<-rbind(out_alpha10r, out_alpha50r)
save(OBoutr, file="../../out/genericOutBrecomb.RData")
