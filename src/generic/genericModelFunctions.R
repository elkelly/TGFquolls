# Modified from northern quoll popualtion model to be generic
# Mods by Ella Kelly and Ben Phillips
# For publication in Kelly & Phillips 2018. "How many, and when? Optimising targeted gene flow for a step change in the environment"

# Calculates effect size and selection threshold for turning gtypes into phenotypes
  # W = initial fitness 
  # f = initial proportion of "1" alleles
  # h2 = heritability
  # Ve = environmental variation
  # nP, number of phenotype loci (equal effects at each locus)
gSetup<-function(W, f, h2, Ve, nP){
  Vg<-h2*Ve/(1-h2)
  eSize<-sqrt(Vg/(2*nP*f*(1-f)))
  sThresh<-qnorm(p = 1-W, mean = 2*nP*f*eSize, sd = sqrt(Ve+Vg))
  list(eSize=eSize, s=sThresh, Vg=Vg)
}
# 
# # check w0 and f0 values
# # make sure e*nP*2< sThresh then no point in running
# w0<-seq(0.00001, 0.1, 0.005)
# for(i in w0){
#   set <- gSetup(i, f = 0.05, h2 = 0.2, Ve = 1, nP = 10)
#   x <- set$eSize*10*2
#   if (x<set$s) print(c(i, "fail"))
# }


# Initialises a matrix with n individuals
  # W is initial proportion of "1" alleles (population fitness)
  # nP, nI, nN are number of phnotypic, incompatibility and neutral loci, respectively
initInds<-function(Nstar, nP, nI, nN, f, Ve, h2, recipient=TRUE){
  n<- Nstar
  S<-rep(1, n) #all hermaphrodites
  A<-rep(1, n) #age 
  # allocate P alleles at random depending on W
  allelesP<-matrix(rbinom(2*n*nP, 1, f), 
                   nrow=n, dimnames=list(NULL, paste("a", 1:(2*nP), sep="")))
  # allocate neutral alleles, fixed for recipient (=0).
  if (recipient) gCode<-0 else gCode<-1
  allelesN<-matrix(rep(gCode, 2*n*nN), nrow=n,
                   dimnames=list(NULL, paste("c", 1:(2*nN), sep="")))
  # allocate incompatibility alleles, fixed for recipient (=0).
  allelesI<-matrix(rep(gCode, 2*n*nI), nrow=n,
                   dimnames=list(NULL, paste("b", 1:(2*nI), sep="")))
  cbind(S, A, allelesP, allelesN, allelesI)	
}

# calculates survival phenotype given alleles, effect size and Ve
  # effect size needs to be scaled carefully against Ve to generate initial fitness. (?????)
aPhen<-function(pop, eSize, Ve, sel){
  gtypes<-pop[,grepl("a", colnames(pop))]
  gtypes<-apply(gtypes, 1, sum)*eSize
  ptypes<-rnorm(length(gtypes), mean=gtypes, sd=sqrt(Ve)) #phenotypes
  as.numeric(ptypes>sel) #thresholded
}

# Beverton Holt population growth
  # Rmax = max number of offspring an individual can have
  # Nstar = carrying capasity given Rmax
bevHolt<-function(N, Rmax, Nstar){
  a<-(Rmax-1)/(Nstar*1) #calculates alpha for respective Rmax and Nstar values
  Rmax/(1+a*N)
}

# generates gametes for each parental individual in the population
recom<-function(pop, lLoci, gCols, recombRate){
  nCols<-sum(gCols) # number of alleles
  genos<-pop[, gCols, drop=FALSE] # genotypes
  even<- which(1:nCols %% 2 == 0) 
  odd<- which(1:nCols %% 2 != 0) 
  nLoci<-nCols/2 
  nCuts<-rpois(nrow(pop), recombRate) #randomly assign number of cuts
  
  lCut<-lapply(nCuts, runif, min=0, max=1) 
  lCut<-mapply(sort, lCut) 
  damSet<-mapply(findInterval, lCut, MoreArgs = list(x = lLoci)) %% 2 == 0 # intervals between cuts
  if (rbinom(n = 1, size = 1, prob = 0.5)==1) damSet<-!damSet #randomly select loci for selection from 
  ifelse(t(damSet), genos[,odd], genos[,even]) # select recombinated genes 
}

# Function that reproduces individuals
  # density dependence acting on number of surviving offspring according bevHolt model
repro<-function(pop, nP, nN, nI, Rmax, Nstar, lLoci, gCols, recombRate){
  if (nrow(pop)<=2) return(pop)
  nPar<-nrow(pop) #number parents/pairs (every individual breeds)
  dams<-pop[sample(nrow(pop), nPar, replace=F),] #randomly order "females"
  sires<-pop[sample(nrow(pop), nPar, replace=F),] #randomly order "males"
  pSurv<- bevHolt(nPar, Rmax, Nstar) #density dependance growth 
  nOff<- rpois(nPar, pSurv) #number of offspring per dam
  if (sum(nOff)==0) return(pop)
  dams<-dams[rep(seq_len(nrow(dams)), nOff),] #rep dams depending on number of offspring
  sires<-sires[rep(seq_len(nrow(sires)), nOff),] #rep mates depending on number of offspring  
  if (is.null(nrow(dams))) return(pop)
  offDam<-recom(dams, lLoci, gCols, recombRate) #select alleles to be passed on 
  offSire<-recom(sires, lLoci, gCols, recombRate) #select alleles to be passed on 
  if (nrow(offDam)!=nrow(offSire)) {offDam<-offDam[sample(nrow(offDam), min(nrow(offDam), nrow(offSire))), ]; offSire<-offSire[sample(nrow(offSire), min(nrow(offDam), nrow(offSire))), ]} #error catch
  orderCols<-order(c(seq(1, 2*(nP+nN+nI), 2), seq(2, 2*(nP+nN+nI), 2))) #reorganise alleles
  alleles<-cbind(offDam, offSire)[, orderCols, drop=FALSE] #bind dam and sire alleles
  nOff<-nrow(alleles) 
  S<-rep(1, nOff)
  A<-rep(0, nOff)
  rbind(pop, cbind(S, A, alleles)) # add offspring to population 
}

# Survival and age (based on age and phenotype)
surv<-function(pop, selection, eSize, Ve, s, outbreeding, nl, h2, alpha){
  if (nrow(pop)==0) return(pop) 
  if (selection==FALSE) phen<-1 else phen<-aPhen(pop, eSize, Ve, s) # threat based selection
  # sum mutually exclusive survival probabilities
  pSurv<-(pop[,"A"]==0)*phen +
    (pop[,"A"]==1)*0 # one year breeding cycle (complete mortality for age=1)
  if (outbreeding==TRUE) pSurv<-pSurv*Vmod(pop, nl, h2, alpha) else pSurv<-pSurv # outbreeding mortality
  surv<-as.logical(rbinom(nrow(pop), size = 1, prob = pSurv)) # survival based on probably (pSurv)
  pop<-pop[surv, drop=FALSE,]
  #if (is.null(nrow(pop))) pop<-rbind(pop)
  pop[,"A"]<-pop[,"A"]+1 # age everyone
  pop
}

# function to score proportion of genome that is from recipient population
  # r.current = proportion of genome from recipient across the population
  # r.new = proportion of loci with recipient alleles
neutral.genomeTracker<-function(pop){
  gCols<-grepl("c", colnames(pop)) # genotype columns
  gCols<-pop[, gCols] # genotypes
  r.current<-1-sum(gCols)/length(gCols) # proportion of genome from recipient across the population
  x<-nrow(gCols)-colSums(gCols)
  nL<-ncol(gCols)/2 #number of loci
  l<-rep(1:nL, each=2) #indexes for each locus
  xAgg<-tapply(x, INDEX = l, FUN = sum)
  r.new<-sum(xAgg>0)/nL # is there at least one individual with recipient allele at each locus?
  list(r.current=r.current, r.new=r.new, popFreq=xAgg)
}

# outputs homozygosity scores for Turelli and Orr model
genomeTracker<-function(pop){
  gCols<-grepl("b", colnames(pop)) # genotype columns
  nCols<-sum(gCols) # number of alleles
  nLoci<-nCols/2 
  nInds<-nrow(pop)
  gCols<-pop[, gCols] # genotypes
  damAllele<-seq(1, nCols, 2) # relevant columns
  sireAllele<-seq(2, nCols, 2)
  damAllele<-gCols[, damAllele] # relevant subset
  sireAllele<-gCols[, sireAllele]
  homozygous<-damAllele==sireAllele # conditions
  homP1<-homozygous & damAllele==0
  homP2<-homozygous & damAllele==1
  p1<-rowSums(homP1)/nLoci
  p2<-rowSums(homP2)/nLoci
  pH<-1-(p1+p2)
  list(p1=p1, p2=p2, pH=pH)
}

# given a breakdown score for double homozygote, h1 and h0 based on simple dosage.
hWeight<-function(H2){
  h1<-0.5*H2
  h0<-0.25*H2
  list(H2, h1, h0)
}

# Turelli and Orr model of D-M incompatibility
  # n is number of loci potentially involved in incompatibility ## is this nN???? 
  # hList is list of breakdown scores for H2, H1, H0 respectively
TOmod<-function(pop, nl, H2){
  pList<-genomeTracker(pop)[1:3]
  hList<-hWeight(H2)
  ExpS<-nl*(pList[[1]]*pList[[2]]*hList[[1]]+(pList[[1]]+pList[[2]])*pList[[3]]*hList[[1]]+pList[[3]]^2*hList[[3]])
  ExpS
}

# Use a simple negative exponential as our link between breakdown score and fitness
  # alpha is some positive constant
Vmod<-function(pop, nl, H2, alpha){
  exp(-alpha*TOmod(pop, nl, H2))
}


# Mother function
  # with Targeted gene flow; step change in environment
TGF<-function(nP, nN, nI, W0, f0, h2, Ve=1, Nstar, burnIn=10, noSel=20, selGens=50, TGFtime, outbreeding, nl, H2, alpha=1, f1, nIntr, Rmax,recombRate){
  # global parameters
  selPars<-gSetup(W0, f0, h2, Ve, nP) # selection parameters
  lLoci<-runif(nP+nN+nI,0,1) #assign random loci locations 
  
  # initialise the population
  pop<-initInds(Nstar,nP, nN, nI, f=f0)
  gCols<-grepl("b", colnames(pop)) | grepl("a", colnames(pop)) | grepl("c", colnames(pop))   # genotype columns
  # bins for collecting stuff
  n<-c()
  phen<-c()
  genome<-c()
  IntProp<-c()
  x<-c()
  
  # function for implementing TGF
  sourceInject<-function(){
    if (nIntr==0) return(NULL)
    sourceSamp<-initInds(nIntr, nP, nN, nI, f=f1, recipient = FALSE)
    pop<<-rbind(pop, sourceSamp)
  }
  
  # function for iterating over generations
  iterate<-function(start, finish, collect=TRUE, selection=FALSE, outbreeding){
    for (gg in start:finish){
      if (gg==TGFtime & nrow(pop)>1) IntProp<<- c(IntProp, nIntr/length(pop))  # record proportion of TGF inds to recipient inds
      if (gg==TGFtime) sourceInject() # add in TGF individuals 
      pop<<-repro(pop, nP, nN, nI, Rmax, Nstar, lLoci, gCols, recombRate) #reproduction
      pop<<-surv(pop, selection, selPars$eSize, Ve, selPars$s, outbreeding, nl, H2, alpha) #survival
      if (collect){
        nr<-nrow(pop)
        n<<-c(n, nr) # new population size appended to popsize vector
        if (nr<=1) {phen<-c(phen, NA); break}
        phen<<-c(phen, mean(aPhen(pop, eSize=selPars$eSize, Ve, sel=selPars$s))) 
        genome<<-rbind(genome, neutral.genomeTracker(pop)) #record genome results
      }
    }
  }
  
  iterate(1, burnIn, collect = FALSE, selection = FALSE, outbreeding) # burn in stage
  iterate(burnIn+1, burnIn+noSel, collect = TRUE, selection = FALSE, outbreeding) # no selection
  iterate(burnIn+noSel+1, burnIn+noSel+selGens, collect = TRUE, selection = TRUE, outbreeding) # selection
  
  list(n=n, phen=phen, pop=pop, genome=genome, IntProp=IntProp)
}

