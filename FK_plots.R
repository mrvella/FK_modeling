library(grid)
library(ggthemes)
library(scales)
library(gridExtra)
library(dplyr)
library(reshape2)
library(ggplot2)
library(cowplot)
# setwd("~/Documents/FK_final")

dir.create("figs", showWarnings = FALSE)

my_theme <- function() {
    (theme_foundation(base_size=16)
        + theme(panel.background = element_rect(colour = NA),
                plot.background = element_rect(colour = NA),
                panel.border = element_rect(colour = NA),
                axis.title.y = element_text(angle=90,vjust =2),
                axis.title.x = element_text(vjust = -0.2),
                axis.line = element_line(colour="black"),
                axis.ticks = element_line(),
                panel.grid.major = element_line(colour="#f0f0f0"),
                panel.grid.minor = element_blank(),
                legend.key = element_rect(colour = NA),
                plot.margin=unit(c(10,5,5,5),"mm"),
                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                strip.text = element_text(face="bold")
        ))
}

plotmanip <- function(x){ # make approach match manuscript format
  x$mortality <- factor(grepl("E",x$mech))
  x$kill <- factor(grepl("F",x$mech))
  levels(x$mortality)[levels(x$mortality)==F] <- c("L-")
  levels(x$mortality)[levels(x$mortality)==T] <- c("E-")
  levels(x$kill)[levels(x$kill)==F] <- c("BK")
  levels(x$kill)[levels(x$kill)==T] <- c("FK")
  x$construct <- factor(paste0(x$kill,x$loci))
  
  levels(x$construct)[levels(x$construct)=="BK1"] <- "BK"
  x$approach <- factor(paste0(x$mortality,x$construct))
  x}

allelemaker <- function(alldat){ # set up for allele frequency plots
  holdthis1 <- alldat %>% subset(loci==2,select=-c(b_freq))
  holdthis1$allele <- "L-FK2: A"
  names(holdthis1)[names(holdthis1)=="a_freq"] <- "freq"
  holdthis2 <- alldat %>% subset(loci==2,select=-c(a_freq))
  holdthis2$allele <- "L-FK2: B"
  names(holdthis2)[names(holdthis2)=="b_freq"] <- "freq"
  holdthis3 <- alldat %>% subset(loci==1,select=-c(b_freq))
  if (nrow(holdthis3)>0){
  holdthis3$allele <- "L-FK1: K"
  names(holdthis3)[names(holdthis3)=="a_freq"] <- "freq"
  }
  newdata <- rbind(holdthis1,holdthis2,holdthis3)
  newdata
}

colors <- c('#c51214','#377eb8','#59bf56')

# make a base time series plot
plot_approaches <- function(dat,xfac = NULL,yfac=NULL){
  ggplot(dat,aes(x=time,y=(sum_F/feq[1]),group=interaction(approach,simnum),color=approach,linetype=approach)) +
        scale_linetype_manual(values = c("E-BK" = 2,"E-FK1" = 2,
                                         "E-FK2" = 2,"L-BK"= 1,
                                         "L-FK1"= 1,"L-FK2"= 1))+
        scale_color_manual(values = c("E-BK" = colors[1],"E-FK1"= colors[2],
                                      "E-FK2" = colors[3],"L-BK" = colors[1],
                                      "L-FK1" = colors[2],"L-FK2" = colors[3])) +
    facet_grid(rows=yfac,cols=xfac,labeller=label_parsed)+
    guides(linetype=guide_legend(keywidth = 4, keyheight = 1)) +
    my_theme() 
}

labelize_factor <- function(dat,myvar){ # make labels for facets using parameter values
  if (typeof(dat[,myvar]) != "factor"){
    dat[,myvar] <- factor(dat[,myvar])
  }
  myvarexpr <- myvar
  if(myvar=="pcosta" | myvar=="cA"){
    myvarexpr <- "c[A]"
  }
  if(myvar=="smate" | myvar=="sM"){
    myvarexpr <- "s[M]"
  }
  levels(dat[,myvar])<- paste(myvarexpr,"==",as.character(levels(dat[,myvar])),sep="")
  dat
}

####################### Figure 1
dat1 <- readRDS(file=paste0(getwd(),"/dat_files/fig1"))
dat1 <- plotmanip(dat1)
dat1 <- labelize_factor(dat1,"beta")
dat1 <- dat1 %>% subset(select=c("time","sum_F","beta","simtype","simnum","r","approach","feq")) %>%
  melt(measure.vars = c("sum_F"),value.name="sum_F")

figure1 <- plot_approaches(dat1,xfac=vars(beta)) +
  geom_line(size=1.5) +
  scale_x_continuous(name="time (days)",expand = c(0.05, 0.05))+
  scale_y_continuous(name="adult females (relative to baseline)",limits=c(0,1),expand = c(0.0025, 0.0025))+
  facet_grid(.~beta,labeller=label_parsed)
figure1
pdf(file = paste0(getwd(),"/figs/fig1.pdf"),width = 11.62,height=8.07)
figure1
dev.off()

####################### Figure 2 and figure S6
fignum <- "2"
fignum <- "s6"
dat2 <- readRDS(file=paste0(getwd(),"/dat_files/fig",fignum,"a"))

bifpoints <- dat2 %>%
  group_by(mech,loci) %>%
  slice(max(which(sum_F>0.001)))
bifpoints <- plotmanip(bifpoints)
bifpoints %>% subset(select=c("approach","r"))
fig1points <- dat2 %>% subset(r==1)
fig1points <- plotmanip(fig1points)
dat2 <- plotmanip(dat2)
dat2_1 <- dat2 %>% subset(sum_F>0.001) %>%
  subset(select=c("r","sum_F","loci","beta","approach","feq","simnum")) %>%
  melt(measure.vars = c("sum_F"),value.name="sum_F") 
dat2_2 <- dat2 %>% subset(sum_F<0.001) %>%
  subset(select=c("r","sum_F","loci","beta","approach","feq","simnum")) %>%
  melt(measure.vars = c("sum_F"),value.name="sum_F") 
dat2_2_2 <- dat2_2 %>% subset(approach == "L-BK")
dat2_2_1 <- dat2_2 %>% subset(approach != "L-BK")
mybreaks <- if(fignum == "2"){
  mybreaks <- seq(0,12,2)} else{
    mybreaks <- seq(0,3,0.5)
  }

labeldf <- data.frame(time=3,sum_F=1.04,letter=c("A","B"))

figure2a <- plot_approaches(dat2_1) +
  scale_x_continuous(name="weekly release ratio (r)",breaks=mybreaks,expand = c(0, 0))+
  scale_y_continuous(name="equilibrium adult females\n (relative to baseline)",expand = c(0.01, 0.0025))+
  geom_line(aes(x=r),size=1.1) +
  geom_line(data=dat2_2_2,aes(x=r),size=1.1) +
  geom_line(data=dat2_2_1,aes(x=r),size=1.1,color = "darkgrey",show.legend=FALSE) +
  guides(linetype=guide_legend(keywidth = 4, keyheight = 1))+
  geom_point(data=bifpoints,aes(x=r),size=3,shape=16,show.legend = FALSE)+
  geom_point(data=fig1points,aes(x=r),size=3,shape=8,show.legend = FALSE)+
  theme(legend.position = c(0.8,0.6)) +
  theme(axis.title.x = element_blank())

dat2b <- readRDS(file=paste0(getwd(),"/dat_files/fig",fignum,"b"))
dat2b <- plotmanip(dat2b)
dat2b$simnum <- 0

figure2b <- plot_approaches(dat2b %>% subset(!is.na(t_0005))) +
  geom_line(aes(x=r,y=t_0005),size=1.3) +
  scale_x_continuous(name=expression(paste("weekly release ratio (",italic("r"),")")), 
                                     limits = c(0,NA), breaks=mybreaks,expand = c(0, 0))+
  scale_y_continuous(name="days until >99.95% reduction",limits=c(0,1000),expand = c(0.0025, 0.0025))+
  theme(legend.position = "none")

myfig <- plot_grid(figure2a,figure2b,ncol = 1, align = "v", labels=c("A","B"),label_size = 18)
pdf(file = paste0(getwd(),"/figs/fig",fignum,".pdf"),width = 11.51,height=8.54)
myfig
dev.off()

####################### Figure 3
dat3 <- readRDS(file=paste0(getwd(),"/dat_files/fig3"))
dat3 <- plotmanip(dat3)
fig2points <- dat3 %>% subset(cA==0.55 & h == 0.5 & sH == 0.2 & sM == 0.1)
fig2points <- labelize_factor(fig2points,"h")
dat3 <- labelize_factor(dat3,"h")

figure3 <- dat3 %>% 
  ggplot(aes(x=sH,y=sM)) + 
  geom_tile(aes(fill=(t_0005))) +
  scale_x_continuous(name = expression("hatching fitness cost ("~italic(s^H)~")"),breaks = c(0,0.5,1),                   
                     labels = c("0","0.5","1")) +
  scale_y_continuous(limits = c(0,1),name = expression("male mating competitiveness fitness cost ("~italic(s^M)~")"), breaks = c(0,0.5,1),
                     labels = c("0","0.5","1")) +
  scale_fill_gradientn(name="days until\n99.95%\nreduction",na.value="white",limits=c(0,500),
                       colours = c("#361813","#61354A","#5B688C","#14A1A8","#61D18A","#EDEF5D"))+
  geom_point(data=fig2points, aes(color=approach), size=3.5,shape=15, show.legend = F) +
  geom_point(data=fig2points, color="black", size=3.5,shape=0, show.legend = F) +
  scale_color_manual(values = c(colors,colors)) +
  coord_fixed() +
  my_theme() +
  facet_grid(h~approach,labeller = label_parsed)
figure3
pdf(file = paste0(getwd(),"/figs/fig",3,".pdf"),width = 13.57,height=7.03)
figure3
dev.off()

####################### Figure S1
dats1 <- readRDS(file=paste0(getwd(),"/dat_files/figs1"))
dats1 <- plotmanip(dats1)
dats1 <- labelize_factor(dats1,"beta")
dats1 <- labelize_factor(dats1,"r")
dats1$cost <- factor(dats1$sM)
levels(dats1$cost) <- c("equal cost","double cost")
dats1 <- dats1 %>% subset(select=c("time","loci","sum_F","beta","cost","sH","simtype","simnum","r","approach","feq")) %>%
  melt(measure.vars = c("sum_F"),value.name="sum_F")
dats1_2 <- dats1 %>% subset(sH==0.4 & time%%40==0)
dats1_2$approach2 <- factor(dats1_2$approach)
dats1 <- dats1 %>% subset(sH==0.2)
dats1$approach <- factor(dats1$approach)
figureS1 <- 
  plot_approaches(dats1)+
  geom_line(size=1.1)+
  scale_x_continuous(name="time (days)",expand = c(0.05, 0.05))+
  geom_point(data=dats1_2,aes(shape=interaction(approach2,cost)),color=colors[3]) +
  scale_shape_manual(values=c(15,17),name="double costs",labels=c("EFK 2","LFK 2")) +
  scale_y_continuous(name="adult females (relative to baseline)",limits=c(0,1),expand = c(0.0025, 0.0025))+
  facet_grid(r~beta,labeller=label_parsed) +
  theme(panel.spacing = unit(1, "lines")) 
figureS1
pdf(file = paste0(getwd(),"/figs/figS1.pdf"),width = 12.92,height=9.19)
figureS1
dev.off()

####################### Figure S2
dats2 <- readRDS(file=paste0(getwd(),"/dat_files/figs2"))
my_feq <- dats2$feq[1]
dats2 <- plotmanip(dats2)
minval <- 1e-4
dats2 <- dats2 %>% subset(sum_F/dats2$feq[1]>=minval,select=c("time","sum_F","approach","init")) %>% 
  melt(measure.vars = "sum_F",value.name="sum_F")

figureS2 <- ggplot(dats2) +
  geom_line(aes(x=time,y=(sum_F/my_feq),group=init),color="black",linetype=1,size=0.5,show.legend=F) +
  scale_x_continuous(name="time (days)",expand = c(0, 0))+
  scale_y_continuous(name="adult females (relative to baseline)",trans = "log10",
                     labels = c("0.0001","0.001","0.01","0.1","1"),
                     limits=c(minval,1),expand = c(0, 0)) +
  my_theme()

figureS2
pdf(file = paste0(getwd(),"/figs/figS2.pdf"),width = 12.92,height=9.19)
figureS2
dev.off()

###################### Figure S5
dats5 <- readRDS(file=paste0(getwd(),"/dat_files/figs5"))
dats5 <- dats5 %>% plotmanip

# find mean and quantiles of times to extinction
grouped_dat <- dats5 %>% group_by(beta,approach,r,sH,sM,cA)
summaryext <- with(grouped_dat,{
  do.call(data.frame,aggregate(t_ext ~ beta + r + loci+approach,FUN=function(x)
      c(mn = mean(x),quantile(x,c(0.025,0.975)),
      len=length(x)))#percent of simulations that went extinct
      ) 
})

summaryext <- summaryext %>% subset(t_ext.len >=200)
#for plotting function only
summaryext$sum_F <- 0
summaryext$time <- 0
summaryext$feq <- 0

ribdat <- summaryext 

figureS5 <- summaryext %>%
  plot_approaches() +
  geom_line(aes(x=r,y=t_ext.mn,group=approach),size=1.5) +
  scale_x_continuous(name=expression(paste("weekly release ratio (",italic(r),")")), limits=c(0,12),breaks=seq(0,12,2),expand = c(0, 0))+
  scale_y_continuous(name="days until extinction",limits=c(0,1000),expand = c(0.0025, 0.0025))+
  geom_ribbon(data=ribdat,aes(x=r,ymin=t_ext.2.5.,ymax=t_ext.97.5.,fill=approach,group=approach),alpha=0.15,show.legend=F) +
  scale_fill_manual(values = c("E-BK" = colors[1],"E-FK1"=colors[2],
                               "E-FK2" = colors[3],"L-BK" = colors[1],
                               "L-FK1" = colors[2],"L-FK2" = colors[3])) +
  theme(legend.position = c(0.8,0.6))
figureS5
pdf(file = paste0(getwd(),"/figs/figS5.pdf"),width = 12.92,height=8.69)
figureS5
dev.off()

############### Figure S7
dats7 <- readRDS(file=paste0(getwd(),"/dat_files/figs7"))
dats7 <- plotmanip(dats7)
# dats6 <- plotmanip(mydat)
dats7 <- labelize_factor(dats7,"sH")
dats7_2 <- allelemaker(dats7)

s7_1 <- dats7 %>%
  melt(measure.vars = c("sum_F"),value.name="sum_F")  %>%
  plot_approaches(xfac = vars(sH)) +
  geom_line(aes(alpha=simtype,size=simtype)) +
  scale_alpha_manual(values=c(1,0.3))+
  scale_size_manual(values=c(1.5,0.4))+
  scale_y_continuous("adult females\n (relative to baseline)",limits=c(0,1),expand = c(0.0025, 0.0025)) +
  scale_x_continuous("time (days)",expand = c(0.05, 0.05))+
  guides(size = F,alpha=F, color=guide_legend(keywidth = 2, keyheight = 1))+
  theme(strip.background = element_blank(),
        strip.text.x = element_blank())+
  theme(legend.position = NULL)

s7_2 <- dats7_2 %>% subset(simnum < 21) %>%
  ggplot() +
  geom_line(aes(x=time,y=freq,alpha=simtype,size=simtype,
                color=allele,group=interaction(simnum,allele,simtype))) +
  scale_size_manual(name="simulation type", values=c(1.5,0.4),labels=c("deterministic","stochastic"))+
  scale_alpha_manual(name="simulation type", values=c(1,0.3),labels=c("deterministic","stochastic"))+
  scale_color_manual(values=c("L-FK2: B"="#00ff00","L-FK2: A"="#006600","L-FK1: K"=colors[2]))+
  scale_y_continuous("allele frequency",expand = c(0.0025, 0.0025)) +
  scale_x_continuous(name=NULL,
                     expand = c(0.05, 0.05))+
  facet_grid(.~sH,labeller = label_parsed) +
  my_theme() +
  guides(size = F,alpha=F,color=guide_legend(keywidth = 2, keyheight = 1))+
  theme(legend.position = NULL)

library(cowplot)
figures7 <- plot_grid(s7_2, s7_1, ncol=1, align="v",labels = c("A","B"),label_size = 18)
pdf(file = paste0(getwd(),"/figs/figS7.pdf"),width = 12.92,height=8.69)
figures7
dev.off()

############  Figure S8
alldat <- readRDS(file=paste0(getwd(),"/dat_files/figs8"))
alldat2 <- readRDS(file=paste0(getwd(),"/dat_files/figs9")) 
alldat <- plotmanip(alldat)
alldat$sH <- alldat$s
alldat$sM <- alldat$smate
alldat$cA <- alldat$pcosta
alldat$sH <- as.numeric(levels(alldat$sH))[alldat$sH]
alldat <- labelize_factor(alldat,"sM")
alldat <- labelize_factor(alldat,"cA")
alldat$b_freq[alldat$sum_F<.001] <- NA
alldat$a_freq[alldat$sum_F<.001] <- NA

# find where the A allele becomes fixed or extint when starting from different initial conditions
afixed <- alldat2
afixed[is.na(alldat$b_freq),"a_freq"] <- 0
afixed <- afixed %>% subset(a_freq>0.99)
afixed$sH <- as.numeric(levels(afixed$sH))[afixed$sH]
afixed <- labelize_factor(afixed,"sM")
afixed <- labelize_factor(afixed,"cA")

aextinct <- alldat2
aextinct[is.na(alldat$b_freq),"a_freq"] <- 0
aextinct <- aextinct %>% subset(sum_F<0.001 & alldat$sum_F >0.001)
aextinct$sH <- as.numeric(levels(aextinct$sH))[aextinct$sH]
aextinct <- labelize_factor(aextinct,"sM")
aextinct <- labelize_factor(aextinct,"cA")

alldat <- anti_join(alldat,afixed,by=c("cA","sM","r","sH"))
alldat <- anti_join(alldat,aextinct,by=c("cA","sM","r","sH"))

figs7_points <- data.frame("cA" = 0.75,"sM" = 0.25,"r"=2,"sH"=c(0.1,0.5,0.9))
figs7_points <- labelize_factor(figs7_points,"sM")
figs7_points <- labelize_factor(figs7_points,"cA")
figures8 <- alldat %>% 
  ggplot(aes(x=r,y=sH)) + 
  geom_tile(aes(fill=b_freq)) +
  scale_fill_gradientn(name="B allele\nfrequency",colours = c("white", "darkgreen","black"), values = c(0,0.99,1),limits=c(0,1))+
  geom_tile(data=afixed,fill="red4")+
  geom_tile(data=aextinct,fill="red2")+
  ylab(expression(paste("homozygote hatching fitness cost (",italic(s^H),")"))) +
  xlab(expression(paste("weekly release ratio (",italic(r),")")))+
  geom_point(data=figs7_points,color="yellow",size=2,shape=15) +
  geom_point(data=figs7_points,color="black", size=2,shape=0) +
  my_theme() +
  facet_grid(sM~cA,labeller = label_parsed)+
  theme(legend.position = NULL)
figures8
pdf(file = paste0(getwd(),"/figs/figS8.pdf"),width = 11.92,height=8.92)
figures8
dev.off()
