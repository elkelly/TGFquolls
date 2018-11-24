# targeted gene flow on generic population model with Outbreeding Depression

source("genericModelFunctions.R")

##################################
###   Command line arguments and   ###
###   defining indices             ###
##################################

## Read in the command line arguments

command_args <- commandArgs(trailingOnly = TRUE)

recb_index <- as.numeric(command_args[1])
alpha_index <- as.numeric(command_args[2])

## Possible ID options

# calculate values of alpha
## Setup for outbreeding scenarios
#computing H2 
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

recb_options <- c(50, 1.33)
alpha_options <- c(0, alpha10, alpha50)

## Extract IDs corresponding to command line argument indices

recb <- recb_options[recb_index]
alpha <- alpha_options[alpha_index]

################
# a run through management space..
managSpaceOBRR <- function(alpha, recb, Nstar=500, h2=0.2, Rmax=3, W0=0.005, props=seq(0, 0.3, 0.02), yr=seq(from=21,to=41,by=1), n=100) {
  #input
  nI<-floor(Nstar*props) # number to introduce
  fName<-paste0("../out/genericNewOB", round(alpha*100, digits = 0), "RR", recb*100, ".RData") # filename to save to
  
  
  #output
  e <-vector() 
  g <-vector()
  g2<-vector()
  popFreq<-matrix(ncol=10)
  Ip<-vector()
  out<-matrix(ncol=7)
  
  for(k in yr) {
    pr<-vector()
    r.current<-vector()
    r.new<-vector()
    Ipp<-vector()
    for(j in nI) {
      cat("Calculating extinction probability for generic model at:",
          "\n Nstar = ", Nstar,
          "\n h2 = ", h2,
          "\n Rmax = ", Rmax,
          "\n Alpha = ", alpha,
          "\n Recombination = ", recb,
          "\n Number Introduced = ", j,
          "\n Year of introduction = ", k-31,
          "\n")
      pb <- txtProgressBar(min = 0, max = n, style = 3)
      for (i in 1:n) {
        setTxtProgressBar(pb, i)
        test<-TGF(nP = 10, nN = 10, nI = 10, W0 = W0, f0 = 0.05, h2=h2, Nstar=Nstar, burnIn=10, noSel=20, selGens=50, TGFtime=k, f1=0.9, nIntr=j, outbreeding=TRUE, nl=10, H2=1, alpha=alpha, Rmax=Rmax, recombRate=recb)    
        e[i]<-ifelse(tail(test$n,1)<=1, TRUE, FALSE)
        g[i]<-tail(test$genome[,1], 1) #for whole pop, pRecI = [,2]
        g2[i]<-tail(test$genome[,2], 1)
        popFreq<<-rbind(popFreq, as.numeric(tail(test$genome[,3], 1)[[1]]))
        Ip[i]<-ifelse(is.null(test$IntProp), NA, test$IntProp)
      }
      close(pb)
      pr<- c(pr, length(e[e==TRUE])/n)
      r.current<-c(r.current, mean(as.numeric(g)))
      r.new<-c(r.new, mean(as.numeric(g2)))
      Ipp<-c(Ipp, mean(Ip, na.rm=TRUE))
    }
    temp<-cbind(rep(k, length(nI)), nI, pr, r.current, r.new, Ipp, props)
    out<- rbind(out, temp)
  }
  out<-out[-1,]
  save(out, file=fName)
}

# Run sensitivity analysis

managSpaceOBRR(alpha, recb)










