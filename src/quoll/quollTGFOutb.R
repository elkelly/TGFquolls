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

#################

# Run scenarios with outbreeding depression
source("quollTGFFunctions.R")
source("parameters.R")

nI<-seq(from=0,to=50,by=5) #values of nIntr
yr<-seq(21,41,1) # year of introduction
n<-100 # number of iterations 


# 50% reduction
e <-vector() 
g <-vector()
g2<-vector()
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
      test<-TGF(demPars = pars, nP = 10, nN = 10, nI = 10, W0 = 0.05, f0=0.05, s0=0.38, K = 500, h2=0.2, burnIn=10, noSel=20, selGens=50, TGFtime=k, f1=0.9, nIntr=j, outbreeding=TRUE, nl=10, H2=1, alpha=alpha50, recombRate=50)    
      e[i]<-ifelse(tail(test$n,1)<=1, TRUE, FALSE)
      g[i]<-tail(test$genome[,1], 1)      
      g2[i]<-tail(test$genome[,2], 1)
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

out_alpha50<-out[2:232,]
out_alpha50<-cbind(out_alpha50, "alpha"=rep("alpha50", nrow(out_alpha50)))


# 10% reduction
e <-vector() 
g <-vector()
g2<-vector()
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
      test<-TGF(demPars = pars, nP = 10, nN = 10, nI = 10, W0 = 0.05, f0=0.05, s0=0.38, K = 500, h2=0.2, burnIn=10, noSel=20, selGens=50, TGFtime=k, f1=0.9, nIntr=j, outbreeding=TRUE, nl=10, H2=1, alpha=alpha10, recombRate=50)    
      e[i]<-ifelse(tail(test$n,1)<=1, TRUE, FALSE)
      g[i]<-tail(test$genome[,1], 1)      
      g2[i]<-tail(test$genome[,2], 1)
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
out_alpha10<-out[2:232,]

out_alpha10<-cbind(out_alpha10, "alpha"=rep("alpha10", nrow(out_alpha10)))

OBout<-rbind(out_alpha10, out_alpha50)
save(OBout, file="../../out/OutB.RData")
