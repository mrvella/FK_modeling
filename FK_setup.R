library(magrittr)
library(deSolve)
library(doParallel)
library(adaptivetau)
library(plyr)

#########################################################

FKode <- function(t,state,params){ #right-hand side of system of ODEs
  with(as.list(c(state,params)),{
    if(loci==2){
      statevals <- c(w1_w1_w2_w2Jf, a_a_w2_w2Jf, w1_a_w2_w2Jf, w1_w1_b_bJf, 
                     a_a_b_bJf, w1_a_b_bJf, w1_w1_w2_bJf, a_a_w2_bJf, w1_a_w2_bJf,
                     w1_w1_w2_w2Jm, a_a_w2_w2Jm, w1_a_w2_w2Jm, w1_w1_b_bJm, 
                     a_a_b_bJm, w1_a_b_bJm, w1_w1_w2_bJm, a_a_w2_bJm, w1_a_w2_bJm,
                     w1_w1_w2_w2F, a_a_w2_w2F, w1_a_w2_w2F, w1_w1_b_bF, a_a_b_bF, 
                     w1_a_b_bF, w1_w1_w2_bF, a_a_w2_bF, w1_a_w2_bF, w1_w1_w2_w2M, 
                     a_a_w2_w2M, w1_a_w2_w2M, w1_w1_b_bM, a_a_b_bM, w1_a_b_bM, 
                     w1_w1_w2_bM, a_a_w2_bM, w1_a_w2_bM)
    } else{statevals <- c(w1_w1Jf,a_aJf,w1_aJf,w1_w1Jm,a_aJm,w1_aJm,w1_w1F,a_aF,w1_aF,w1_w1M,a_aM,w1_aM)}
    
    juvsum <- sum(statevals[J_index])
    #find birth rates
    sum_M_mate <- sum(statevals[male_index]*mateM)
    
    # use predefined list of probs instead of looping through pijk
    myprob <-  lapply(off_combos_list,function(y)
      sum(statevals[y[,1]]*statevals[y[,2]]*y[,5]*y[,4]/(sum_M_mate)))
    myprob <- unlist(myprob)
    
    B <- 0.5* lambda * myprob
    BF <- B * wF
    BM <- B * wM
    
    outvals <- rep(0,3*genonum)
    for (i in 1:length(genotype_names)){
      outvals[i] <- BF[i]-statevals[i]*(alpha*juvsum)^(beta-1)-muJ*statevals[i]-v*statevals[i]
      outvals[i+genonum] <- BM[i]-statevals[i+genonum]*(alpha*juvsum)^(beta-1)-muJ*statevals[i+genonum]-v*statevals[i+genonum]
      outvals[i+2*genonum] <-  v*statevals[i]*gammaF[i]-muF*statevals[i+2*genonum]
      outvals[i+3*genonum] <- v*statevals[i+genonum]*gammaM[i]-muM*statevals[i+3*genonum]+r*equilM/7*rel[i]
    }

    list(outvals)
  })
}

FKstoch2 <- function(params,transitions){ #run stochastic simulations using adaptivetau package
    simdat <- ssa.adaptivetau(init.values=round(params$states),transitions=transitions,rateFunc=ratefunc,params=params,tf=params$tend)
    simdat <- data.frame(simdat)
    
    simdat$sum_F <- rowSums(simdat[grep("F",names(simdat))])
    simdat$sum_J <- rowSums(simdat[grep("J",names(simdat))])
    sum_Jm <- rowSums(simdat[grep("Jm",names(simdat))])
    simdat$sum_M <- rowSums(simdat[grep("M",names(simdat))])
    
    simdat <- within(simdat,{
        simtype <- "stoch"
        loci <- as.factor(params$loci)
        init <- params$init
        r <- params$r
        beta <- as.factor(params$beta)
        feq <- params$equilF
        sH <- as.factor(params$sH)
        sM <- as.factor(params$sM)
        cA <- as.factor(params$cA)
        h <- as.factor(params$h)
        mech <-params$mech
        if (params$loci==1){
            a_freq <- (a_aJm + 0.5*w1_aJm)/sum_Jm
            b_freq <- 1 - a_freq
        } else {
            a_freq <- (a_a_w2_w2Jm + a_a_b_bJm + a_a_w2_bJm + 0.5* (w1_a_w2_w2Jm + w1_a_b_bJm + w1_a_w2_bJm))/sum_Jm
            b_freq <- (w1_w1_b_bJm + a_a_b_bJm + w1_a_b_bJm + 0.5* (w1_w1_w2_bJm + a_a_w2_bJm + w1_a_w2_bJm))/sum_Jm
        }
    })
    
    simdat
}

ratefunc <- function(x,params,t){ #calculation of rates for states x at time t for adaptivetau
    genonum <- params$genonum
    state_names <- names(params$states)
    sum_J <- sum(x[params$J_index]) #use names(x) instead??? also, could put grepl into params
    sum_F <- sum(x[params$fem_index]) #grepl("F",state_names)])
    sum_M <- sum(x[params$male_index])
    sum_M_mate <- sum(x[params$male_index]*params$mateM)
    death_rates <- x*c(rep(params$muJ,2*genonum),rep(params$muF,genonum),rep(params$muM,genonum))
    
    if (sum_M==0 | sum_F==0){ 
        birth_ratesF <- rep(0,genonum) #check if no adult males or females and set births to 0
        birth_ratesM <- rep(0,genonum)
    } else{ 
        total_birth <- sum_F*params$lambda

        myprob <-  lapply(params$off_combos_list,function(y)
            sum(x[y[,1]]*x[y[,2]]*y[,5]*y[,4]/(sum_F*sum_M_mate)))
        myprob <- unlist(myprob)
        
        births <- 0.5*total_birth * myprob
        birth_ratesF <- births * params$wF
        birth_ratesM <- births * params$wM
    }
    emerge_ratesM <- rep(params$v,genonum)*x[(1+genonum):(2*genonum)]#gamma is in transition
    emerge_ratesF <- rep(params$v,genonum)*x[1:genonum]
    
    dens_rates <- rep((params$alpha*sum_J)^(params$beta-1),2*genonum)*x[1:(2*genonum)]
    rel_rates <- params$r*params$equilM/7
    allrates <- c(death_rates,birth_ratesF,birth_ratesM,emerge_ratesM,emerge_ratesF,dens_rates,rel_rates)

    allrates
}

runsims <- function(output="all",mech="LFK",loci,determ,stochsims,rin,equilFin,tendin,betain,initin=1,sHin=0,sMin=0,cAin=0.5,hin=0.5){
    genonum <- 3^loci
    if(loci == 1){
      offprobs <- readRDS("offprobs_1locus")
    } else{
      offprobs <- readRDS("offprobs_2locus")
    }

    Parray <- offprobs[[1]]
    genotypes <- offprobs[[2]]
    genotype_names <- offprobs[[3]]
    
    #make off_combos_list to speed up birth rate calculations
    nonzeros <- data.frame(which(Parray!=0,arr.ind = T))
    names(nonzeros) <- c("f","m","j")
    nonzeros$prob <- apply(nonzeros,1,function(x) Parray[x[1],x[2],x[3]]) #new column with P(ijk)
    nonzeros$f <- nonzeros$f + genonum*2   #match index of states
    nonzeros$m <- nonzeros$m + genonum*3

    rel <- rep(0,genonum)   #initialize release
    gammaF <- rep(1,genonum) #initialize viability of females
    gammaM <- rep(1,genonum) #initialize viability of males
    mateM <- rep(1,genonum)
    wF <- rep(1,genonum)     #initialize fitness of females
    
    sH <- sHin
    sM <- sMin
    h <- hin                #assume additive...
    cA <- cAin #proportion of the fitness cost of A allele in 2-locus setting
    mycostpheno <- 0        #phenotype cost
    if (loci==2){
        costs <- list(a=sH*cA,b=sH*(1-cA),pheno=mycostpheno,
                      aM=sM*cA,bM=sM*(1-cA)) 
        rel[grepl("a_a_b_b",genotype_names)] <- 1   #release this genotype. (Multiplied by equilib and r in rhs)
        gammaF[grepl("a",genotype_names) & grepl("b",genotype_names)] <- 0 #unviable if have both copies
        
        wF[grepl("a_a_b_b",genotype_names)] <- (1-(costs$a+costs$b)) 
        wF[grepl("w1_a_b_b",genotype_names)] <- (1-(costs$a*h+costs$b)) 
        wF[grepl("a_a_w2_b",genotype_names)] <- (1-(costs$a+costs$b*h)) 
        wF[grepl("w1_a_w2_b",genotype_names)] <- (1-(costs$a*h+costs$b*h))
        wF[grepl("w1_w1_w2_b",genotype_names)] <- (1-(costs$b*h))
        wF[grepl("w1_w1_b_b",genotype_names)] <- (1-(costs$b))
        wF[grepl("a_a_w2_w2",genotype_names)] <- (1-(costs$a))
        wF[grepl("w1_a_w2_w2",genotype_names)] <- (1-(costs$a*h))
        
        # wF[sum(grep("a",unlist(strsplit(genotype_names),NULL)))]
        
        wF[grepl("a",genotype_names) & grepl("b",genotype_names)] <-
            wF[grepl("a",genotype_names) & grepl("b",genotype_names)] - costs$pheno 
        
        # mateM[grepl("a",genotype_names) | grepl("b",genotype_names)] <- 1-sM
        mateM[grepl("a_a_b_b",genotype_names)] <- 1-(costs$aM+costs$bM)
        mateM[grepl("w1_a_b_b",genotype_names)] <- (1-(costs$aM*h+costs$bM)) 
        mateM[grepl("a_a_w2_b",genotype_names)] <- (1-(costs$aM+costs$bM*h)) 
        mateM[grepl("w1_a_w2_b",genotype_names)] <- (1-(costs$aM*h+costs$bM*h))
        mateM[grepl("w1_w1_w2_b",genotype_names)] <- (1-(costs$bM*h))
        mateM[grepl("w1_w1_b_b",genotype_names)] <- (1-(costs$bM))
        mateM[grepl("a_a_w2_w2",genotype_names)] <- (1-(costs$aM))
        mateM[grepl("w1_a_w2_w2",genotype_names)] <- (1-(costs$aM*h))
    } else{
        costs <- list(a=sH,b=0,pheno=mycostpheno,aM=sM)
        rel[grepl("a_a",genotype_names)] <- 1 #release this genotype. (Multiplied by equilib and r in rhs)
        gammaF[grepl("a",genotype_names)] <- 0 #unviable if have both copies
        
        wF[grepl("a_a",genotype_names)] <- (1-costs$a)
        wF[grepl("w1_a",genotype_names)] <- (1-costs$a*h)
        
        wF[grepl("a",genotype_names)] <- wF[grepl("a",genotype_names)] - costs$pheno 
        
        mateM[grepl("a_a",genotype_names)] <- (1-costs$aM)
        mateM[grepl("w1_a",genotype_names)] <- (1-costs$aM*h)
    }
    wM <- wF
    if(mech=="EFK"){
        wF <- wF * gammaF
    }
    if(mech=="EBS"){
        wF <- wF * gammaF
        wM <- wF
    }
    if(mech=="LBS"){ #make juveniles also equal to zero fitness if it acts early
        gammaM <- gammaF
    }
    
    states <- rep(0,genonum*4)
    statenames <- as.vector(vapply(c("Jf","Jm","F","M"),function(y) 
        vapply(genotype_names,function(x) paste0(x,y),""),rep("",3^loci))
    )    
    names(states) <- statenames
    
    #finish making off_combos_list
    
    nonzeros$mateM <- unlist(lapply(nonzeros$m,function(x) mateM[x-genonum*3]))
    off_combos_list <- lapply(seq(1,genonum),function(x) nonzeros %>% subset(j==x)) #list with elements being combos making that offspring
    off_combos_list <- lapply(off_combos_list,as.matrix)

    params <- list(muJ=0.03,muM=0.28,muF=0.10,
        lambda=8,v=0.14,
        loci = loci,
        genonum=genonum,
        Parray = Parray,
        genotype_names=genotype_names,
        genotypes = genotypes,
        gammaF=gammaF,gammaM=gammaM,
        sH=sH,
        sM=sM,
        h=h,
        mateM=mateM,
        cA = cA,
        wF=wF,wM=wM,rel=rel,
        equilF=equilFin,
        beta=betain,
        r=rin,
        tend=tendin,
        init = initin,
        mech = mech,
        J_index = 1:(2*genonum),
        off_combos_list=off_combos_list,
        male_index=(3*genonum+1):(4*genonum),
        fem_index=(2*genonum+1):(3*genonum)
    )
    params <- within(params,{
        R0=v*lambda/(2*muF*(muJ+v)) #R0 and equilibrium wt population calculations (Roberts et al supplements)
        equilJ=1/v*equilF*(2*muF)
        alpha=1/equilJ*((muJ+v)*(R0-1))^(1/(beta-1))
        equilM=v*equilJ/(2*muM)
        # alpha = 2e-4 #use this code if specifying alpha instead of equilF
        # equilJ=1/alpha*((muJ+v)*(R0-1))^(1/(beta-1))
        # equilF=v*equilJ/(2*muF)
        # equilM=v*equilJ/(2*muM)
    })

    if (loci==2){
      states[c("w1_w1_w2_w2Jf","w1_w1_w2_w2Jm","w1_w1_w2_w2F","w1_w1_w2_w2M")] <- c(rep(0.5*params$equilJ,2),params$equilF,params$equilM)
      ##### watch out with these lines, producing different initial allele frequencies!
      if(initin == 2.22){
        states[c("w1_w1_b_bF")] <- c(params$equilF)
        states[c("w1_w1_b_bM")] <- c(params$equilM)
        initin <- 1
      }
      if(initin == 3.33){
        states[c("a_a_w2_w2F")] <- c(params$equilF)
        states[c("a_a_w2_w2M")] <- c(params$equilM)
        initin <- 1
      }
      #####
    } else{ 
      states[c("w1_w1Jf","w1_w1Jm","w1_w1F","w1_w1M")] <- c(rep(0.5*params$equilJ,2),params$equilF,params$equilM)
    }   
    
    states <- states * initin #if starting with suppressed population
    params$states <- states
    
    times <- seq(0,params$tend)
    
    ########## run deterministic simulations
    if(determ==0){datlist <- list()
    } else{
        simdat <- data.frame(ode(y=states,times=times,func <- FKode,parms=params))#,atol = 1e-14, rtol = 1e-14))
        simdat$sum_F <- rowSums(simdat[names(states[grep("F",names(states))])])
        if(output == "equilib") {
          # repeatedly run until there is little change in total number of females
          while(simdat$sum_F[nrow(simdat)-1]-simdat$sum_F[nrow(simdat)] > 1e-5){
            times <- seq(simdat[nrow(simdat),"time"],simdat[nrow(simdat),"time"]+params$tend/2)
            states <- as.double(simdat[nrow(simdat),] %>% subset(select=-c(time,sum_F)))
            names(states) <- statenames
            simdat <- data.frame(ode(y=states,times=times,func <- FKode,parms=params))
            simdat$sum_F <- rowSums(simdat[names(states[grep("F",names(states))])])
          }
          output = "lastrow"
        }

        simdat$sum_J <- rowSums(simdat[names(simdat[grep("J",names(simdat))])])
        sum_Jm <- rowSums(simdat[names(simdat[grep("Jm",names(simdat))])])
        simdat$sum_M <- rowSums(simdat[names(simdat[grep("M",names(simdat))])])

        simdat <- within(simdat,{
            time <- time
            sum_F <- sum_F
            simtype <- "determ"
            loci <- as.factor(loci)
            init <- params$init
            simnum <- 0
            r <- params$r
            beta <- as.factor(params$beta)
            feq <- params$equilF
            sH <- as.factor(params$sH)
            sM <- as.factor(params$sM)
            cA <- as.factor(params$cA)
            h <- as.factor(params$h)
            mech <- mech
            if (loci==1){ #calculate frequencies for later use
                a_freq <- (a_aJm + 0.5*w1_aJm)/sum_Jm
                b_freq <- 1 - a_freq
            } else {
                a_freq <- (a_a_w2_w2Jm + a_a_b_bJm + a_a_w2_bJm + 0.5* (w1_a_w2_w2Jm + w1_a_b_bJm + w1_a_w2_bJm))/sum_Jm
                b_freq <- (w1_w1_b_bJm + a_a_b_bJm + w1_a_b_bJm + 0.5* (w1_w1_w2_bJm + a_a_w2_bJm + w1_a_w2_bJm))/sum_Jm
            }
        })
        
        if(output == "lastrow"){ # to output last row only for numerical equilibrium
            simdat <- simdat[nrow(simdat),]
        }
        
        simdat <- simdat %>% subset(
          select=c("time","sum_F","a_freq","b_freq","mech","loci","beta","r","sH","sM","h","cA","feq","init","simtype","simnum")
        )
        datlist <- list(simdat)}
    
    ############# stochastic simulations
    totalsims <- stochsims
    simnum <- stochsims
    
    ## set up transitions for adaptive tau
    {   state_names <- names(params$states)
        genonum <- 3^params$loci
        num_states <- genonum*4
        loci <- params$loci
        minustrans <- rep(-1,num_states)
        plustrans <- minustrans * -1
        
        death_trans <- Map(setNames, as.list(minustrans),state_names)
        birth_transF <- Map(setNames, as.list(plustrans[1:genonum]),state_names[1:genonum])
        birth_transM <- Map(setNames, as.list(plustrans[(genonum+1):(2*genonum)]),state_names[(genonum+1):(2*genonum)])
        
        emerge_transF <- lapply(as.list(1:genonum),function(x) {vec <- c(-1,+1)
                        names(vec) <- c(state_names[x],state_names[x+2*genonum])
                        vec})
        emerge_transF[gammaF==0] <- lapply(emerge_transF[gammaF==0],function(x) x[1]) #if gamma is 0, only subtract from juveniles
        
        emerge_transM <- lapply(as.list((genonum+1):(2*genonum)),function(x) {vec <- c(-1,+1)
                        names(vec) <- c(state_names[x],state_names[x+2*genonum])
                        vec})
        emerge_transM[gammaM==0] <- lapply(emerge_transM[gammaM==0],function(x) x[1]) #if gamma is 0, only subtract from juveniles
        
        dens_trans <- Map(setNames, as.list(minustrans[1:(2*genonum)]),state_names[1:(2*genonum)])
        
        if(loci==1){
            rel_trans <- list(c(a_aM = +1))
        }else{ rel_trans <- list(c(a_a_b_bM = +1))
        }
        
        transitions <- unlist(list( #list of named vectors of deltas
            death_trans, #birth into every genotype (based on sumF*lambda, males, and pijks)
            birth_transF,
            birth_transM,
            emerge_transM,
            emerge_transF,
            dens_trans,
            rel_trans
        ),recursive=F)
    }
    
    while(simnum > 0){
        simdat <- FKstoch2(params,transitions)
        #keep only daily time points
        simdat$time <- floor(simdat$time)
        simdat <- simdat[!duplicated(simdat$time),]
        simdat <- simdat %>% subset(
              select=c("time","sum_F","a_freq","b_freq","mech","loci","beta","r","sH","sM","h","cA","feq","init","simtype")
        )
          
        simdat$simnum <- simnum
        
        mynum <- as.character(length(datlist))
        datlist[[mynum]] <- simdat
        simnum <- simnum - 1
    }

    dat <- rbind.fill(datlist)
    
    ####### for outputting extinct info instead of sim data
    if(output == "extinct"){
      splitdf <- lapply(levels(factor(dat$simnum)),function(z) subset(dat,simnum==z))
      lastentries <- unlist(lapply(splitdf,function(z) z[nrow(z),]["sum_F"]))
      ext_times <- unlist(lapply(splitdf,function(z) ifelse(nrow(z[z$sum_F==0,]["time"])>0,
                                                            min(z[z$sum_F==0,]["time"]),NA)))
      Jext_times <- unlist(lapply(splitdf,function(z) ifelse(nrow(z[z$sum_J==0,]["time"])>0,
                                                             min(z[z$sum_J==0,]["time"]),NA)))
      t_05 <- unlist(lapply(splitdf,function(z) {
        under <- which(z$sum_F<(z$sum_F[1]*.05))
        if(length(under)>0){out <- z$time[min(under)]
        } else{out <- NA}
        out
      }))
      t_01 <- unlist(lapply(splitdf,function(z) {
        under <- which(z$sum_F<(z$sum_F[1]*.01))
        if(length(under)>0){out <- z$time[min(under)]
        } else{out <- NA}
        out
      }))
      t_001 <- unlist(lapply(splitdf,function(z) {
        under <- which(z$sum_F<(z$sum_F[1]*.001))
        if(length(under)>0){out <- z$time[min(under)]
        } else{out <- NA}
        out
      }))
      t_0005 <- unlist(lapply(splitdf,function(z) {
        under <- which(z$sum_F<(z$sum_F[1]*.0005))
        if(length(under)>0){out <- z$time[min(under)]
        } else{out <- NA}
        out
      }
      ))
        
    dat <- data.frame(last = lastentries,t_ext = ext_times,J_ext=Jext_times,
                      t_05=t_05,t_01=t_01,t_001=t_001,t_0005=t_0005,
                      beta=params$beta,r=params$r,loci=params$loci,sH=params$sH,sM=params$sM,cA=params$cA,
                      h=params$h,mech=mech,feq=params$equilF)
    }

    dat
}

produce_combos <- function(){ # function for generating data frame of simulations for input combinations
  lociinlist <- as.list(c(1))
  mechlist <- as.list(c("LBS","EBS"))
  allcombos1 <- expand.grid(r = rinlist,beta=betainlist,loci = lociinlist,initin=initinlist,hin=hinlist,
                            sMin=sMlist,sHin=sHinlist,cAin=cAlist,mech=mechlist)
  lociinlist <- as.list(c(1,2))
  mechlist <- as.list(c("LFK","EFK"))
  allcombos2 <- expand.grid(r = rinlist,beta=betainlist,loci = lociinlist,initin=initinlist,hin=hinlist,
                            sMin=sMlist,sHin=sHinlist,cAin=cAlist,mech=mechlist)
  rbind(allcombos1,allcombos2)
}
