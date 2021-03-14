source('FK_setup.R')
 
# run simulations needed to produce figures
# specify figure number: 1,2,3,s1,s2,s5,s6,s7,s8,s9
plotnum <- "1"

## default inputs for simulation. Adjust for each plot as needed.
equilFin <- 2000 #number of equilibrium females
tendin <- 1000   #number of days to run simulation
determ <- 1      # 1 to run deterministic simulation, 0 not to
stochsims <- 0   # number of stochastic simulations to run

# specify parameters as lists
rinlist <- as.list(c(1))    # release ratio
hinlist <- as.list(c(0.5))  # fitness cost degree of dominance
cAlist <- as.list(c(0.55))  # proportion of cost to allele A in 2-locus
betainlist <- as.list(c(3)) # strength of density dependence
initinlist <- as.list(c(1)) 
sHinlist <- as.list(c(0.2)) # hatching fitness cost
sMlist <- as.list(c(0.1))   # male mating competitiveness cost

# outtype: "all", "lastrow" for last time point only, "extinct" for extinct
outtype <- "all"

if (plotnum == "1"){
  #simple time series
  tendin <- 400
  betainlist <- as.list(c(2,3,4))

  allcombos <- produce_combos()

} else if(plotnum == "2a" |plotnum == "s6a"){
  # numequilib plot - specify r granularity for each mech 
  tendin <- 1000
  
  lociinlist <- as.list(c(1))
  mechlist <- as.list(c("LBS"))
  rinlist <- as.list(c(seq(0,1.4,.01),seq(2,12,1)))
  if (plotnum == "s6a"){
    betainlist <- as.list(c(2))
    rinlist <- as.list(c(seq(0,0.65,.01),seq(0.75,3,.25)))
  }
  allcombos1 <- expand.grid(r = rinlist,beta=betainlist,loci = lociinlist,initin=initinlist, hin=hinlist,
                            sMin=sMlist,sHin=sHinlist,cAin=cAlist,mech=mechlist)
  rinlist <- as.list(c(seq(0,5.5,0.05),seq(5.51,6.1,0.01),seq(7,12,1)))
  if (plotnum == "s6a"){
    rinlist <- as.list(c(seq(0,1.4,.1),seq(1.41,1.65,.01),seq(1.75,3,.25)))
  }
  mechlist <- as.list(c("EBS"))
  allcombos2 <- expand.grid(r = rinlist,beta=betainlist,loci = lociinlist,initin=initinlist,hin=hinlist,
                            sMin=sMlist,sHin=sHinlist,cAin=cAlist,mech=mechlist)
  
  lociinlist <- as.list(c(1,2))
  mechlist <- as.list(c("LFK"))
  rinlist <- as.list(c(seq(0,3.2,.02),seq(4,12,1)))
  if (plotnum == "s6a"){
    rinlist <- as.list(c(seq(0,0.4,.1),seq(0.45,0.85,.01),seq(1,3,.25)))
  }
  allcombos3 <- expand.grid(r = rinlist,beta=betainlist,loci = lociinlist,initin=initinlist,hin=hinlist,
                            sMin=sMlist,sHin=sHinlist,cAin=cAlist,mech=mechlist)
  mechlist <- as.list(c("EFK"))
  rinlist <- as.list(c(seq(0,1,0.02),seq(1.1,4,0.1),seq(4.1,6,0.01),seq(7,12,1)))
  if (plotnum == "s6a"){
    rinlist <- as.list(c(seq(0,0.6,.1),seq(.66,1,.02),seq(1.01,1.3,.01),seq(1.5,3,.25)))
  }
  allcombos4 <- expand.grid(r = rinlist,beta=betainlist,loci = lociinlist,initin=initinlist,hin=hinlist,
                            sMin=sMlist,sHin=sHinlist,cAin=cAlist,mech=mechlist)
  allcombos <- rbind(allcombos1,allcombos2,allcombos3,allcombos4)
  outtype <- "equilib"
}  else if(plotnum == "2b" | plotnum == "s6b"){
  # time to near extinction - specify r granularity for each mech 
  lociinlist <- as.list(c(1))
  mechlist <- as.list(c("LBS"))
  rinlist <- as.list(c(seq(0.5,1,.01),seq(1.1,2,0.1),seq(2.25,12,0.25)))
  if (plotnum == "s6b"){
    betainlist <- as.list(c(2))
    rinlist <- as.list(c(seq(0.4,0.7,.01),seq(0.72,1.2,.02),seq(1.25,3,.25)))
  }
  allcombos1 <- expand.grid(r = rinlist,beta=betainlist,loci = lociinlist,initin=initinlist, hin=hinlist,
                            sMin=sMlist,sHin=sHinlist,cAin=cAlist,mech=mechlist)
  rinlist <- as.list(c(seq(5.9,6.3,0.01),seq(6.4,12,.2)))
  if (plotnum == "s6b"){
    rinlist <- as.list(c(seq(1.3,1.9,.01),seq(2,3,.25)))
  }
  mechlist <- as.list(c("EBS"))
  allcombos2 <- expand.grid(r = rinlist,beta=betainlist,loci = lociinlist,initin=initinlist,hin=hinlist,
                            sMin=sMlist,sHin=sHinlist,cAin=cAlist,mech=mechlist)
  
  lociinlist <- as.list(c(1,2))
  mechlist <- as.list(c("LFK"))
  rinlist <- as.list(c(seq(2.3,3.1,.01),seq(3.2,12,0.2)))
  if (plotnum == "s6b"){
    rinlist <- as.list(c(seq(0.5,0.9,.005),seq(0.92,1.2,.02),seq(1.3,1.4,.1),seq(1.5,3,.25)))
  }
  allcombos3 <- expand.grid(r = rinlist,beta=betainlist,loci = lociinlist,initin=initinlist,hin=hinlist,
                            sMin=sMlist,sHin=sHinlist,cAin=cAlist,mech=mechlist)
  mechlist <- as.list(c("EFK"))
  rinlist <- as.list(c(seq(4.2,5.9,0.01),seq(6,12,.2)))
  if (plotnum == "s6b"){
    betainlist <- as.list(c(2))
    rinlist <- as.list(c(seq(0.75,1.4,.01),seq(1.45,2,.05),seq(2.25,3,.25)))
  }
  allcombos4 <- expand.grid(r = rinlist,beta=betainlist,loci = lociinlist,initin=initinlist,hin=hinlist,
                            sMin=sMlist,sHin=sHinlist,cAin=cAlist,mech=mechlist)
  allcombos <- rbind(allcombos1,allcombos2,allcombos3,allcombos4)
  outtype <- "extinct"
} else if(plotnum == "3"){
  #for deterministic time to near extinction comparison
  tendin <- 500
  rinlist <- as.list(c(7))
  
  hinlist <- as.list(c(0,0.5,1))
  sHinlist <- as.list(c(seq(0,1,0.01)))
  sMlist <- sHinlist
  
  allcombos <- produce_combos()
  outtype <- "extinct"

} else if(plotnum=="s1"){
  #Figure 1 with multiple beta and r, and 2 locus twice as costly 
  tendin <- 400
  rinlist <- as.list(c(0.3,1,7))
  betainlist <- as.list(c(2,3,4))
  lociinlist <- as.list(c(1))
  mechlist <- as.list(c("LBS","EBS","LFK","EFK"))
  
  allcombos1 <- expand.grid(r = rinlist,beta=betainlist,loci = lociinlist,initin=initinlist, hin=hinlist,
                            sMin=sMlist,sHin=sHinlist,cAin=cAlist,mech=mechlist)
  lociinlist <- as.list(c(2))

  mechlist <- as.list(c("LFK","EFK"))
  allcombos2 <- expand.grid(r = rinlist,beta=betainlist,loci = lociinlist,initin=initinlist,hin=hinlist,
                            sMin=sMlist,sHin=sHinlist,cAin=cAlist,mech=mechlist)
  sHinlist[[1]] <- sHinlist[[1]]*2
  sMlist[[1]] <- sMlist[[1]]*2
  allcombos3 <- expand.grid(r = rinlist,beta=betainlist,loci = lociinlist,initin=initinlist,hin=hinlist,
                            sMin=sMlist,sHin=sHinlist,cAin=cAlist,mech=mechlist)
  
  allcombos <- rbind(allcombos1,allcombos2,allcombos3)

} else if(plotnum == "s2"){
  #bistability
  tendin <- 400
  initinlist <- as.list(c(1e-3,2e-3,2.5e-3,2.8e-3,3e-3,3.5e-3,4e-3,6e-3,1e-2,5e-2,.1,1))
  mechlist <- as.list(c("LFK"))
  lociinlist <- as.list(c(1))
  rinlist <- as.list(c(2))
  allcombos <- expand.grid(r = rinlist,beta=betainlist,loci = lociinlist,initin=initinlist,hin=hinlist,
                            sMin=sMlist,sHin=sHinlist,cAin=cAlist,mech=mechlist)
} else if(plotnum == "s5"){
  # time to near extinction - specify r granularity for each mech 
  determ <- 0
  stochsims <- 300
  lociinlist <- as.list(c(1))
  mechlist <- as.list(c("LBS"))
  rinlist <- as.list(c(seq(0.5,1,.1),seq(1.1,2,0.1),seq(2.25,12,0.25)))
  allcombos1 <- expand.grid(r = rinlist,beta=betainlist,loci = lociinlist,initin=initinlist, hin=hinlist,
                            sMin=sMlist,sHin=sHinlist,cAin=cAlist,mech=mechlist)
  rinlist <- as.list(c(seq(5.9,6.3,0.1),seq(6.4,12,.2)))
  mechlist <- as.list(c("EBS"))
  allcombos2 <- expand.grid(r = rinlist,beta=betainlist,loci = lociinlist,initin=initinlist,hin=hinlist,
                            sMin=sMlist,sHin=sHinlist,cAin=cAlist,mech=mechlist)
  
  lociinlist <- as.list(c(1,2))
  mechlist <- as.list(c("LFK"))
  rinlist <- as.list(c(seq(2.3,3.1,.1),seq(3.2,12,0.2)))
  allcombos3 <- expand.grid(r = rinlist,beta=betainlist,loci = lociinlist,initin=initinlist,hin=hinlist,
                            sMin=sMlist,sHin=sHinlist,cAin=cAlist,mech=mechlist)
  mechlist <- as.list(c("EFK"))
  rinlist <- as.list(c(seq(4.2,5.9,0.1),seq(6,12,.2)))
  allcombos4 <- expand.grid(r = rinlist,beta=betainlist,loci = lociinlist,initin=initinlist,hin=hinlist,
                            sMin=sMlist,sHin=sHinlist,cAin=cAlist,mech=mechlist)
  allcombos <- rbind(allcombos1,allcombos2,allcombos3,allcombos4)
  outtype <- "extinct"

} else if(plotnum == "s7"){
  #for loss of symmetry stochastic simulations
  tendin <- 1000
  stochsims <- 10
  stochsims <- 2
  cAlist <- as.list(c(0.75))
  sHinlist <- as.list(c(0.1,0.5,.9))
  sMinlist <- as.list(c(0.25))
  rinlist <- as.list(c(2))
  mechlist <- as.list(c("LFK"))
  lociinlist <- as.list(c(1,2))
  allcombos <- expand.grid(r = rinlist,beta=betainlist,loci = lociinlist,initin=initinlist,hin=hinlist,
                          sMin=sMlist,sHin=sHinlist,cAin=cAlist,mech=mechlist)

}else if(plotnum == "s8" | plotnum == "s9"){
  #for loss of symmetry starting from perturbed population
  #need s9 to produce complete s8 figure
  #s20 and s21 are omitted from supplemental infromation:
  #s20 uses EFK, and s21 starts with the A allele in the population instead of B
  initinlist <- as.list(c(2.22))
  mechlist <- as.list(c("LFK"))
  if (plotnum == "s9"){
    initinlist <- as.list(c(3.33))
  }
  sHinlist <- as.list(seq(0,1,0.01))
  rinlist <- as.list(seq(0,5,0.1))
  sMlist <- as.list(seq(0,0.5,0.25))
  cAlist <- as.list(c(0.5,0.75,1))
  lociinlist <- as.list(c(2))
  allcombos <- expand.grid(r = rinlist,beta=betainlist,loci = lociinlist,initin=initinlist,hin=hinlist,
                            sMin=sMlist,sHin=sHinlist,cAin=cAlist,mech=mechlist)
  
  tendin <- 4000
  outtype <- "lastrow"
}

# myseed <- round(runif(1,min=1,max=100000))
# myseed
myseed <- 22
logfile <- paste0("mylog",plotnum,".txt")
writeLines(c(""),logfile)

library(foreach)
library(doParallel)
library(doRNG)
cores = detectCores() - 1
cores
# cores=3
myclust <- makeCluster(cores[1],outfile=paste0("clustout",plotnum))
registerDoParallel(myclust)
numberruns <- nrow(allcombos)
allcombos$runnum <- seq(1,numberruns)
registerDoRNG(myseed)

myout <- foreach(x=iter(allcombos,by='row'),
                 .packages=c("plyr","deSolve","magrittr","adaptivetau")) %do% {
                   
                   cat(paste("Starting iteration", x[["runnum"]][[1]],"out of", numberruns, "at", Sys.time(),"\n"),file=logfile,append=TRUE)
                   
                   out <- runsims(output=outtype,mech=x[["mech"]][[1]],loci=x[["loci"]][[1]],
                                  determ=determ,stochsims=stochsims,r=x[["r"]][[1]],
                                  equilFin=equilFin,tendin=tendin,betain=x[["beta"]][[1]],
                                  initin=x[["initin"]][[1]], sMin=x[["sMin"]][[1]],hin=x[["hin"]][[1]],
                                  sHin=x[["sHin"]][[1]],cAin=x[["cAin"]][[1]])
                   out
                 }

cat("Run completed.\n",file=logfile,append=TRUE)
stopCluster(myclust)

mydat <- rbind.fill(myout)

saveRDS(mydat,file=paste0(getwd(), "/dat_files/fig",plotnum))
