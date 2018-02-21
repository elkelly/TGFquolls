# Heavily modified from IBM reported in (Kelly & Phillips 2016)
# Mods by Ben Phillips and Ella Kelly
# For publication in Kelly & Phillips 2018. "How many, and when? Optimising targeted gene flow for a step change in the environment"

# Calculates effect size and selection threshold for turning gtypes into phenotypes
  # W = initial fitness and initial proportion of "1" alleles
  # h2 = heritability
  # Ve = environmental variation
  # nP, number of phenotype loci (equal effects at each locus)
gSetup<-function(W, h2, Ve, nP){
  Vg<-h2*Ve/(1-h2)
  eSize<-sqrt(Vg/(2*nP*W*(1-W)))
  sThresh<-qnorm(p = 1-W, mean = 2*nP*W*eSize, sd = sqrt(Ve+Vg))
  list(eSize=eSize, s=sThresh, Vg=Vg)
}

# Initialises a matrix with n individuals
  # W sets initial population fitness (0-1) wrt toad-smart phenotype
  # nP, nI, nN are number of phnotypic, incompatibility and neutral loci, respectively
initInds<-function(K, nP, nI, nN, W, Ve, h2, recipient=TRUE){
  n<- K #ddep on females only
  S<-rbinom(n, 1, 0.5) #equal sex ratio: 0=Female; 1=Male
  A<-rep(1, n) #age 
  # allocate P alleles at random
  allelesP<-matrix(rbinom(2*n*nP, 1, W), 
                   nrow=n, dimnames=list(NULL, paste("a", 1:(2*nP), sep="")))
  # allocate neutral alleles, fixed for recipient (=0).
  if (recipient) gCode<-0 else gCode<-1
  allelesN<-matrix(rep(gCode, 2*n*nN), nrow=n,
                   dimnames=list(NULL, paste("c", 1:(2*nN), sep="")))
  # allocate incompatibility alleles, fixed for recipient (=0).
  if (recipient) gCode<-0 else gCode<-1
  allelesI<-matrix(rep(gCode, 2*n*nI), nrow=n,
                   dimnames=list(NULL, paste("b", 1:(2*nI), sep="")))
  cbind(S, A, allelesP, allelesN, allelesI)	
}

# calculates attack phenotype given alleles, effect size, Ve, and mu
  # effect size needs to be scaled carefully against Ve and mu to generate initial fitness.
aPhen<-function(pop, eSize, Ve, sel){
  gtypes<-pop[,grepl("a", colnames(pop))]
  nAll<-ncol(gtypes)
  gtypes<-apply(gtypes, 1, sum)*eSize
  ptypes<-rnorm(length(gtypes), mean=gtypes, sd=sqrt(Ve)) #phenotypes
  as.numeric(ptypes>sel) #thresholded
}

# Density Dependence function
  # survival of babies
  # s0 is base survival rate
  # hits p=0.3 at K
dDep<- function(nF, K, s0){
  s0*exp(nF*(log(0.3)/K)) # calculates shape of the curve 
}

# Extract mated females
matedF<- function(pop){
  nF<-sum(pop[,"S"]==0)
  nM<-sum(pop[,"S"]==1)
  if (nF==0 | nM==0) return(NULL) # one sex missing, no good
  nMates<-nM
  if (nMates>nF) {nMates<-nF}
  matedFem<-pop[sample(which(pop[,"S"]==0), nMates, replace=TRUE),, drop=FALSE]
  #if (is.null(nrow(matedFem))) matedFem<-rbind(matedFem)
  matedFem
}  

# function to extract a mate for each of our mating females
mates<-function(pop, matedFem){
  mRows<-which(pop[,"S"]==1) #males
  nF<-nrow(matedFem)
  males<-pop[sample(mRows, nF, replace=TRUE),, drop=FALSE] #mates for each female
  #if (is.null(nrow(males))) males<-rbind(males)
  males
}

# generates gametes for a population
gametes2<-function(pop){
  #browser()
  gCols<-grepl("b", colnames(pop)) | grepl("a", colnames(pop)) | grepl("c", colnames(pop))   # genotype columns
  nCols<-sum(gCols) # number of alleles
  nLoci<-nCols/2 
  nInds<-nrow(pop)
  gCols<-pop[, gCols, drop=FALSE] # genotypes
  
  damAllele<-seq(1, nCols, 2) # relevant columns
  sireAllele<-seq(2, nCols, 2)
  damAllele<-gCols[, damAllele] # relevant subset
  sireAllele<-gCols[, sireAllele]
  
  z<-matrix(rbinom(nLoci*nInds, 1, 0.5), nrow=nInds)
  ifelse(z==0, damAllele, sireAllele)
}

# Function that reproduces individuals
  # density dependence acting on number of surviving offspring according to K
repro<-function(pop, K, nP, nN, nI, s0){
  # collect breeding females and find their mates
  females<-matedF(pop)
    if (is.null(females)) return(pop)
  males<-mates(pop, females)
  nF<-nrow(females)
  # surviving offspring based on density dependence:
  pSurv<-dDep(nF, K, s0)
  nOff<- rbinom(nF, 8, pSurv) # 8 max no. of offspring
  if (sum(nOff)==0) return(pop)
  offRows<-rep(1:nF, times=nOff)
  offDam<-gametes2(females[offRows, , drop=FALSE])
  offSire<-gametes2(males[offRows, , drop=FALSE])
  orderCols<-order(c(seq(1, 2*(nP+nN+nI), 2), seq(2, 2*(nP+nN+nI), 2)))
  alleles<-cbind(offDam, offSire)[, orderCols, drop=FALSE]
  nOff<-nrow(offDam)
  S<-rbinom(nOff, 1, 0.5) #equal sex ratio: 0=Female; 1=Male
  A<-rep(0, nOff)
  rbind(pop, cbind(S, A, alleles))
}

# Survival and age (based on age and phenotype)
surv<-function(pop, fsurv1, fsurv2, msurv1, selection, eSize, Ve, s, outbreeding, nl, h2, alpha){
  if (nrow(pop)==0) return(pop)
  if (selection==FALSE) phen<-1 else phen<-aPhen(pop, eSize, Ve, s)
  # sum mutually exclusive survival probabilities
  pSurv<-(pop[,"A"]==0)*phen+
    (pop[,"S"]==0 & pop[,"A"]==1)*fsurv1+
    (pop[,"S"]==0 & pop[,"A"]==2)*fsurv2+
    (pop[,"S"]==0 & pop[,"A"]==3)*0+
    (pop[,"S"]==1 & pop[,"A"]==1)*msurv1+
    (pop[,"S"]==1 & pop[,"A"]==2)*0
  if (outbreeding==TRUE) pSurv<-pSurv*Vmod(pop, nl, h2, alpha) else pSurv<-pSurv 
  surv<-as.logical(rbinom(nrow(pop), size = 1, prob = pSurv))
  pop<-pop[surv, drop=FALSE,]
  #if (is.null(nrow(pop))) pop<-rbind(pop)
  pop[,"A"]<-pop[,"A"]+1
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
  r.new<-sum(xAgg>0)/nL
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

# Mother function I
  # dem.pars is a vector with fs1, fs2, msurv1
  # no Targeted gene flow; step change in environment
noTGF<-function(demPars, nP, nN, nI, W, selGens, h2=0.3, Ve=1, s0=1, K=1000, burnIn=20){
  fsurv1<-demPars["fsurv1"]
  fsurv2<-demPars["fsurv2"]
  msurv1<-demPars["msurv"]
  selPars<-gSetup(W, h2, Ve, nP)
  pop<-initInds(K, nP, nN, nI, W)
  n<-nrow(pop)
  phen<-mean(aPhen(pop, eSize=selPars$eSize, Ve, sel=selPars$s))
  for (bb in 1:burnIn){
    pop<-repro(pop, K, nP, nN, nI, s0) # females reproduce (density dep)
    pop<-surv(pop, fsurv1, fsurv2, msurv1, selection=FALSE, selPars$eSize, Ve, selPars$s, outbreeding=FALSE)
    nr<-nrow(pop)
    n<-c(n, nr) # new population size appended to popsize vector
    if (nr<=1) {phen<-c(phen, NA); break}
    phen<-c(phen, mean(aPhen(pop, eSize=selPars$eSize, Ve, sel=selPars$s)))
  }
  for (ss in 1:selGens){
    pop<-repro(pop, K,  nP, nN, nI, s0) # females reproduce (density dep)
    pop<-surv(pop, fsurv1, fsurv2, msurv1, selection=TRUE, selPars$eSize, Ve, selPars$s)
    nr<-nrow(pop)
    n<-c(n, nr) # new population size appended to popsize vector
    if (nr<=1) {phen<-c(phen, NA); break}
    phen<-c(phen, mean(aPhen(pop, eSize=selPars$eSize, Ve, sel=selPars$s)))
  }
  list(n=n, phen=phen, pop=pop)
}

# Mother function II
  # dem.pars is a vector with fs1, fs2, msurv1
  # with Targeted gene flow; step change in environment
TGF<-function(demPars, nP, nN, nI, W0, h2=0.3, Ve=1, s0=1, K=1000, burnIn=10, noSel=20, selGens=50, TGFtime=500,outbreeding, nl, H2, alpha=1, W1, nIntr){
  # global parameters
  fsurv1<-demPars["fsurv1"]
  fsurv2<-demPars["fsurv2"]
  msurv1<-demPars["msurv"]
  selPars<-gSetup(W0, h2, Ve, nP)
  
  # initialise the population
  pop<-initInds(K, nP, nN, nI, W0)
  
  # bins for collecting stuff
  n<-c()
  phen<-c()
  genome<-c()
  IntProp<-c()
  x<-c()
  
  # function for implementing TGF
  sourceInject<-function(){
    if (nIntr==0) return(NULL)
    sourceSamp<-initInds(nIntr, nP, nN, nI, W1, recipient = FALSE)
    pop<<-rbind(pop, sourceSamp)
  }
  
  # function for iterating over generations
  iterate<-function(start, finish, collect=TRUE, selection=FALSE, outbreeding){
    for (gg in start:finish){
      if (gg==TGFtime & nrow(pop)>1) IntProp<<- c(IntProp, nIntr/length(pop))
      if (gg==TGFtime) sourceInject()
      pop<<-repro(pop, K, nP, nN, nI, s0)
      pop<<-surv(pop, fsurv1, fsurv2, msurv1, selection, selPars$eSize, Ve, selPars$s, outbreeding, nl, H2, alpha)
      if (collect){
        nr<-nrow(pop)
        n<<-c(n, nr) # new population size appended to popsize vector
        if (nr<=1) {phen<-c(phen, NA); break}
        phen<<-c(phen, mean(aPhen(pop, eSize=selPars$eSize, Ve, sel=selPars$s)))
        genome<<-rbind(genome, neutral.genomeTracker(pop))
      }
    }
  }
  
  iterate(1, burnIn, collect = FALSE, selection = FALSE, outbreeding)
  iterate(burnIn+1, burnIn+noSel, collect = TRUE, selection = FALSE, outbreeding)
  iterate(burnIn+noSel+1, burnIn+noSel+selGens, collect = TRUE, selection = TRUE, outbreeding)
  
  list(n=n, phen=phen, pop=pop, genome=genome, IntProp=IntProp)
}

