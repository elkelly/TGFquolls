#sensitivty analysis for targeted gene flow on northern quoll population model 

# Run no outbreeding scenario of targeted gene flow
source("quollTGFFunctions.R")
source("parameters.R")

##################################
###   Command line arguments and   ###
###   defining indices             ###
##################################

## Read in the command line arguments

command_args <- commandArgs(trailingOnly = TRUE)

Nstar_index <- as.numeric(command_args[1])
h2_index <- as.numeric(command_args[2])

## Possible ID options

Nstar_options <- c(125,500,1000)
h2_options <- c(0.1,0.2, 0.3)

## Extract IDs corresponding to command line argument indices

Nstar <- Nstar_options[Nstar_index]
h2 <- h2_options[h2_index]

##################
# a run through management space..
managSpace <- function(Nstar, h2, props=seq(0, 0.3, 0.02), yr=seq(from=21,to=41,by=1), n=100) {
  #input
  nI<-floor(Nstar*props) # number to introduce
  fName<-paste0("../out/quollNewK", Nstar, "h", h2*10, ".RData") # filename to save to
  W0Tab<-read.csv("../data/Sensitivity_quoll.csv")
  W0<-W0Tab$w0[W0Tab$population.size==Nstar & 
                  W0Tab$heritability==h2 & 
                  W0Tab$outbreeding.depression==0 &
                  W0Tab$recombination.rate==0.5]
  #output
  e <-vector() 
  g <-vector()
  g2<-vector()
  popFreq<-matrix(ncol=10)
  Ip<-vector()
  temp<-vector()
  out<-matrix(ncol=7)
  
  for(k in yr) {
    pr<-vector()
    r.current<-vector()
    r.new<-vector()
    Ipp<-vector()
    for(j in nI) {
      cat("Calculating extinction probability for quoll model at:",
          "\n K = ", Nstar,
          "\n h2 = ", h2,
          "\n Number Introduced = ", j,
          "\n Year of introduction = ", k-31,
          "\n")
      pb <- txtProgressBar(min = 0, max = n, style = 3)
      for (i in 1:n) {
        setTxtProgressBar(pb, i)
        test<-TGF(demPars = pars, nP = 10, nN = 10, nI = 10, W0 = W0, f0=0.05, s0=0.38, K = Nstar, h2=h2, burnIn=10, noSel=20, selGens=50, TGFtime=k, f1=0.9, nIntr=j, outbreeding=FALSE, nl=10, H2=1, alpha=1, recombRate=50)    
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

managSpace(Nstar, h2)

