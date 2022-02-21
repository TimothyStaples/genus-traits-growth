# ################################################################################## ####
# Code goal: Compare tree biomass acquisition based on environment, traits and genus ####
# Author:    Tim Staples                                                             ####
# Date:      21/02/2022                                                              ####  
# ################################################################################## ####

# ####
rm(list=ls())

setwd("ENTER LOCATION OF THIS FILE HERE")

# Libraries ####
library(nlme) # mixed-effect modelling
library(vegan) # for rarefied richness etc
library(MuMIn) # model selection and averaging
library(usdm) # VIF
library(parallel) # for parallel processing
library(maptools) # for extracting env values from spatial layers
library(lme4) # for mixed-modelling
library(FD)
library(vioplot)
library(rworldmap)

# Global functions ####
sapply(list.files(path="./Functions", pattern=".R", full.names=TRUE),
       source)

# center continuous columns ("predictors") to mean=0 and sd=1 for modelling
center.data<-function(data, predictors){
  
  model.data<-data[apply(is.na(as.matrix(data[,!is.na(match(colnames(data), 
                                                            predictors))])), 
                         MARGIN=1, 
                         function(x){sum(x)==0}),]

  # re-center variables (as we are working with a data-subset here). This shouldn't affect the data
  # as a mean of 0 and sd of 1 won't change the values if variable is already centered
  
  center.data<-model.data[,colnames(model.data) %in% predictors]

  # now center data
  if(!is.null(dim(center.data))){
    center.data.output<-sapply(center.data, function(x){(x-mean(x))/sd(x)})
  } else {center.data.output<- (center.data-mean(center.data))/sd(center.data)}
  
  model.data[,colnames(model.data) %in% predictors]=center.data.output
  
  return(model.data)
}

# create all possible two-way interactions from a vector of strings (representing covariates)
two.way.interactions<-function(variables){
  interactions.temp<-combn(unique(c(variables, variables)), 2)
  return(paste(interactions.temp[1,], interactions.temp[2,], sep=":"))
}

# convert 2 long-lat pairs into a distance in km
gc.dist<-function(lat1, lon1, lat2, lon2){
  
  deg2rad <- function(deg) {(deg * pi) / (180)} # custom function to turn degrees to radians
  
  R<-6371 # Radius of the earth in km
  dLat<-deg2rad(lat2-lat1) # deg2rad below
  dLon<-deg2rad(lon2-lon1); 
  a<-sin(dLat/2) * sin(dLat/2) + cos(deg2rad(lat1)) * cos(deg2rad(lat2)) * sin(dLon/2) * sin(dLon/2)
  
  c = 2 * atan2(as.numeric(sqrt(a)), as.numeric(sqrt(1-a)))
  d = R * c # Distance in km
  return(d)
}

# see if eucalyptus and acacia CIs are non-overlapping
are.inters.nonoverlapping<-function(list.of.data.summaries){
  
  t(sapply(list.of.data.summaries, function(x){
    
    # get variable names
    x<-as.data.frame(x)
    x$var<-gsub("genusEucalyptus|:", "", rownames(x))
    
    # work out whether interaction terms are overlapping using 95% CIs
    x$lowerCI<-x[,"Estimate"]-1.96*x[,"Std. Error"]
    x$upperCI<-x[,"Estimate"]+1.96*x[,"Std. Error"]
    
    x.to.split<-x[x$var %in% names(table(x$var))[table(x$var)==2],]
    
    sapply(split(x.to.split, f=x.to.split$var), function(y){
      y[1, "lowerCI"] > y[2, "upperCI"] | 
        y[1, "upperCI"] < y[2, "lowerCI"]
    })
  }))
}

# see if eucalyptus and acacia estimates are on different sides of 0
are.marginals.overlapping<-function(list.of.data.summaries){
  
  t(sapply(list.of.data.summaries, function(x){
    
    # get variable names
    x<-as.data.frame(x)
    x$var<-gsub("genusEucalyptus|:", "", rownames(x))
    
    x.to.split<-x[x$var %in% names(table(x$var))[table(x$var)==2],]
    
    sapply(split(x.to.split, f=x.to.split$var), function(y){
      (y[1, "Estimate"] > 0 &  y[2, "Estimate"] < 0) | 
      (y[1, "Estimate"] < 0 &  y[2, "Estimate"] > 0)
      })
        
  }))
}

# what's the probability of seeing a more extreme value from the all model?
prob.of.tail<-function(estimate, dmean, dsd){
  
  temp.p<-pnorm(q=estimate,
                mean=dmean,
                sd=dsd, lower.tail=TRUE)
  
  if(temp.p>0.5){
    temp.p<-pnorm(q=estimate,
                  mean=dmean,
                  sd=dsd, lower.tail=FALSE)
  }
  
  return(temp.p)
}

# Calculate slopes at particular points on interaction plane
inter.slope.calc <- function(model, var, genus){
  
  terms<-summary(model)$coefficients[grepl(var, rownames(summary(model)$coefficients)),]
  
  if(genus == "Acacia"){
    arid.bin<-grepl("aridity", rownames(terms))
    euc.bin<-grepl("Eucalyptus", rownames(terms))
    
    slopes <- terms[!arid.bin & !euc.bin, 1] + terms[arid.bin & !euc.bin, 1] * arid.quantile
    
  }
  
  if(genus == "Eucalyptus"){
    terms<-summary(model)$coefficients[grepl(var, 
                                                      rownames(summary(model)$coefficients)),]
    arid.bin<-grepl("aridity", rownames(terms))
    euc.bin<-grepl("Acacia", rownames(terms))
    
    slopes = terms[!arid.bin & !euc.bin, 1] + terms[arid.bin & !euc.bin, 1] * arid.quantile
    
  }
  
  return(slopes)

}

# Calculate the standard error of interaction sloeps slopes
interaction.slope.se.fun<-function(model, var){
  
  model.vcov<-as.data.frame(as.matrix(vcov(model)))
  
  var.b1<-model.vcov[var,var]
  var.b3<-model.vcov[rownames(model.vcov) %in% 
                       c(paste0(var,":aridity.index"), 
                         paste0("aridity.index:", var)),
                     colnames(model.vcov) %in% 
                       c(paste0(var,":aridity.index"), 
                         paste0("aridity.index:", var))]
  cov.b1.b3<-model.vcov[var, 
                        colnames(model.vcov) %in%
                        c(paste0(var,":aridity.index"), 
                        paste0("aridity.index:", var))]
  SE<-sqrt(var.b1 + var.b3 * arid.quantile^2 + 2*arid.quantile * cov.b1.b3)
  
  # at aridity = 0, SE should be for the same as for partial regression slope of
  # var
  #test<-sqrt(var.b1 + var.b3 * 0^2 + 2*0* cov.b1.b3)
  
  return(SE)
  
}

# ####
# MODELLING ####
#         Prep data for modelling ####
mixed.plant<-read.csv("./Data/mixed plant tree.csv", stringsAsFactors = TRUE)
mixed.site<-read.csv("./Data/mixed site tree.csv", stringsAsFactors = TRUE)

mixed.plant$planting <- as.factor(mixed.plant$planting)
mixed.site$planting <- as.factor(mixed.site$planting)

# drop dead plants
mixed.plant.alive<-mixed.plant[mixed.plant$dead==0,]

# drop small plants (less than 100g)
mixed.plant.alive<-mixed.plant.alive[mixed.plant.alive$tot.ag.mass>=0.1,]

mixed.plant.alive<-merge(mixed.plant.alive,
                         mixed.site[,c("pl.p","plot.area", "age", "IBRA", "IBRAsub",
                                       "solar.annual", "aridity.index", "density")], all.x=TRUE)

colnames(mixed.plant.alive)
species<-read.csv("./Data/species.list.csv", header=TRUE)

mixed.plant.alive<-merge(mixed.plant.alive,
                         species[,c("species","genus", "fam")])

# drop large plots and small plots with few stems. Max area is 0.04 ha (20m^2)
# and min alive stems is 10
plant.count<-as.data.frame(table(mixed.plant.alive$pl.p))
head(plant.count)

mixed.site.small<-droplevels(mixed.site[mixed.site$plot.area<=0.04 &
                                          mixed.site$plant.count>=10,])

mixed.site.small<-droplevels(mixed.site[mixed.site$plot.area<=0.04 &
                                          mixed.site$pl.p %in%
                                          plant.count$Var1[plant.count$Freq>=10],])

mixed.plant.small<-droplevels(mixed.plant.alive[mixed.plant.alive$pl.p %in% 
                                                  mixed.site.small$pl.p,])

# calculate total plot biomass
plot.mass<-with(mixed.plant.small, tapply(tot.ag.mass, pl.p, sum))

# for each row get the total plot biomass minus focal biomass = nhood biomass
mixed.plant.small$nhood.mass<-plot.mass[match(mixed.plant.small$pl.p,
                                              names(plot.mass))] -
                                         mixed.plant.small$tot.ag.mass
 
# split data into list of plot dataframes for easier calculations
plant.by.pl.p<-split(mixed.plant.small, f=mixed.plant.small$pl.p)

# calculate nhood stats if we don't have a completed file
if(!"centered data for modelling.csv" %in% list.files("./Data")){
  
  plant.small.nhood<-lapply(plant.by.pl.p, function(x){
    
    # calculate stats for each tree in each plot separately
    nhood.stats<-do.call("rbind", lapply(1:dim(x)[1], function(y){
      
      # some pre work for the FD package
      x.mat<-x[-y, c("SLA","WD","MH")][!duplicated(x[-y, "species"]),]
      rownames(x.mat)=x[-y, "species"][!duplicated(x[-y, "species"])]
      x.mat<-x.mat[complete.cases(x.mat),]
      abun.mat<-t(data.matrix(table(droplevels(x[-y,"species"]))))
      abun.mat<-abun.mat[,colnames(abun.mat) %in% rownames(x.mat)]
      
      if(dim(x.mat)[1]>1){
        n.FD<-dbFD(x=x.mat, a=abun.mat, messages=FALSE)
      } else{ n.FD<-data.frame(FRic=0, FEve=NA, FDis=NA)}
      
      
      # now make single-row dataframe of values
      data.frame(
        # intraspecific density
        intra.sp.dens=sum(x[-y,"species"]==x[y,"species"], na.rm=TRUE)/x[y,"plot.area"], 
        
        # interspecific density
        inter.sp.dens=sum(x[-y,"species"]!=x[y,"species"], na.rm=TRUE)/x[y,"plot.area"],
        
        # intrageneric density
        intra.gen.dens=sum(x[-y,"genus"]==x[y,"genus"], na.rm=TRUE)/x[y,"plot.area"], 
        
        # intergeneric density
        inter.gen.dens=sum(x[-y,"genus"]!=x[y,"genus"], na.rm=TRUE)/x[y,"plot.area"],
        
        # intraspecific proportion
        intra.sp.prop = sum(x[-y,"species"]==x[y,"species"], na.rm=TRUE)/ length(x[-y, "species"]) ,
        
        inter.sp.prop = sum(x[-y,"species"]!=x[y,"species"], na.rm=TRUE)/ length(x[-y, "species"]) ,
        
        # mean nhood SLA
        nSLA.mean=mean(x[-y, "SLA"], na.rm=TRUE),
        # - SLA delta
        dnSLA.f=mean(x[y, "SLA"] - x[-y, "SLA"], na.rm=TRUE),
        # - SLA abs.delta
        abs.dnSLA.f=mean(abs(x[-y, "SLA"] - x[y, "SLA"]), na.rm=TRUE),
        
        # mean nhood WD
        nWD.mean=mean(x[-y, "WD"], na.rm=TRUE),
        # - delta
        dnWD.f=mean(x[-y, "WD"] - x[y, "WD"], na.rm=TRUE),
        # - abs.delta
        abs.dnWD.f=mean(abs(x[-y, "WD"] - x[y, "WD"]), na.rm=TRUE),
        
        # mean nhood SM
        nSM.mean=mean(x[-y, "SM"], na.rm=TRUE),
        # - delta
        dnSM.f=mean(x[-y, "SM"] - x[y, "SM"], na.rm=TRUE),
        # - abs.delta
        abs.dnSM.f=mean(abs(x[-y, "SM"] - x[y, "SM"]), na.rm=TRUE),
        
        # mean nhood MH
        nMH.mean=mean(x[-y, "MH"], na.rm=TRUE),
        # - delta
        dnMH.f=mean(x[-y, "MH"] - x[y, "MH"], na.rm=TRUE),
        # - abs.delta
        abs.dnMH.f=mean(abs(x[-y, "MH"] - x[y, "MH"]), na.rm=TRUE),
        
        # nhood richness
        n.rarerich=rarefy(table(droplevels(x[-y,"species"])), 9),
        
        # neighbourhood trait range
        n.FRr=n.FD["FRic"],
        
        # nhood trait evenness
        n.FEm=n.FD["FEve"],
        
        # nhood trait divergence
        n.FDm=n.FD["FDis"],
        
        # nhood mass compared to focal mass
        nhood.mass = mean(x[-y, "tot.ag.mass"]),
        n.mass.ratio = mean(x[y, "tot.ag.mass"] / x[-y, "tot.ag.mass"]),
        n.ln.mass.ratio = mean(log(x[y, "tot.ag.mass"]) / log(x[-y, "tot.ag.mass"]))
      
      )
    }))
    
    return(as.data.frame(cbind(x, nhood.stats)))
  })
  # bind lists back together
  plant.small.nhood<-do.call("rbind", plant.small.nhood)
  
  write.csv(plant.small.nhood, "./Data/neighbourhood stats df.csv")
} else {
  plant.small.nhood<-read.csv("./Data/neighbourhood stats df.csv")
}

plant.small.nhood.save <- plant.small.nhood

summary(plant.small.nhood$tot.ag.mass>=0.1)

  # drop plants where traits were estimated so they aren't selected as focals.
  plant.small.focals<-droplevels(plant.small.nhood[rowSums(plant.small.nhood[,c("SLA.n", "WD.n",
                                                                                "SM.n", "MH.n")])==0 &
                                                     plant.small.nhood$genus %in% c("Acacia",
                                                                                    "Eucalyptus"),])
  
  # convert biomass to biomass per year
  hist(plant.small.focals$tot.ag.mass/plant.small.focals$age)
  
  plant.small.focals$mass.per.year<-log(plant.small.focals$tot.ag.mass)/
    plant.small.focals$age
  hist(plant.small.focals$mass.per.year)
  
  plant.small.focals$raw.MH<-plant.small.focals$MH
  plant.small.focals$MH<-log(plant.small.focals$MH)
  plant.small.focals$solar.annual<-log(plant.small.focals$solar.annual)
  plant.small.focals$intra.sp.dens<-sqrt(plant.small.focals$intra.sp.dens)
  plant.small.focals$inter.sp.dens<-sqrt(plant.small.focals$inter.sp.dens) 
  plant.small.focals$intra.gen.dens<-sqrt(plant.small.focals$intra.gen.dens)
  plant.small.focals$inter.gen.dens<-sqrt(plant.small.focals$inter.gen.dens)
  plant.small.focals$density<-sqrt(plant.small.focals$density)
  plant.small.focals$log.age<-log(plant.small.focals$age)

  summary(plant.small.focals$genus)
  
  temp <- plant.small.focals$pl.p

  temp1<-plant.small.nhood.save[plant.small.nhood.save$pl.p %in% temp,]  
  
  temp2 <- mixed.plant[mixed.plant$pl.p %in% temp,]
  tapply(temp2$species, temp2$species, function(x){length(unique(x))})
  dim(plant.small.nhood.save)
  temp2$genus <- substr(temp2$species,
                        1,
                        regexpr(" ", temp2$species)-1)
  length(unique(temp2$species[temp2$unknown == 0 & temp2$dead ==0]))
  
  dim(temp2[temp2$genus %in% c("Acacia", "Eucalyptus") &
                temp2$tot.ag.mass < 0.1 &
                rowSums(temp2[,c("SLA.n","WD.n","MH.n")]==0)==3, ])
  
  head(temp2[,c("MH","MH.n","MH.dist")])

  temp3 <- temp2[!temp2$species %in% plant.small.focals$species &
                  temp2$unknown == 0 & temp2$dead ==0,]

  sum(!duplicated(temp3$species))
  summary(temp3$SLA.n[!duplicated(temp3$species)]==0)
  
  table(is.na(temp3$MH[!duplicated(temp3$species)]))

#         Center variables ####
  
  # before centering we are going to remove some plots
  plant.small.focals<-plant.small.focals[!is.na(plant.small.focals$FEve),]
  
  # centre all variables
  write.csv(plant.small.focals, "./Data/uncentered data for modelling.csv")
  
  predictors<-c("aridity.index", "solar.annual", "log.age", "density", "intra.sp.dens",
                "inter.sp.dens", "intra.gen.dens", "inter.gen.dens", "SLA",
                "WD", "MH", "n.rarerich", "FRic", "FEve", "FDis", "intra.sp.prop", "plot.area")
  
  center.data.output<-sapply(plant.small.focals[,predictors],
                             function(x){
                               (x-mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE)
                             })
  
  plant.cent<-as.data.frame(cbind(plant.small.focals[,!colnames(plant.small.focals) %in%
                                                       predictors],
                                  center.data.output))
  sapply(plant.cent[,predictors], function(x){mean(x, na.rm=TRUE)})
  sapply(plant.cent[,predictors], function(x){sd(x, na.rm=TRUE)})
  
  write.csv(plant.cent, "./Data/centered data for modelling.csv")
  
#         Base models ####
  
plant.cent<-read.csv("./Data/centered data for modelling.csv")
plant.uncent<-read.csv("./Data/uncentered data for modelling.csv")

vif(plant.cent[,c("SLA","WD","MH","aridity.index","density", "intra.sp.prop",
                  "plot.area", "solar.annual", "log.age", "FEve", "FDis", "n.rarerich")])

#                       Null Model ####

null.m <- lmer(mass.per.year ~ 1 +
                  (1|IBRAsub/planting/pl.p) + (1|species), data=plant.cent,
               control = lmerControl(optimizer="bobyqa"))

summary(null.m)
write.csv(summary(null.m)$coefficients, "./Outputs/null model coefs.csv")

#                       Genera-less Model ####

genless.m <- lmer(mass.per.year ~ aridity.index + solar.annual + 
                                  SLA + WD + MH + log.age +
                                  density + (intra.sp.prop * plot.area) + 
                                  n.rarerich + FEve + FDis +
                                  (1|IBRAsub/planting/pl.p) + (1|species), data=plant.cent,
                  control = lmerControl(optimizer="bobyqa"))

write.csv(summary(genless.m)$coefficients, "./Outputs/genus-less model coefs.csv")

#                       Genus model (ints) ####

genera.int.m <- lmer(mass.per.year ~ aridity.index + solar.annual + 
                    SLA + WD + MH + log.age +
                    density + (intra.sp.prop * plot.area) + 
                    n.rarerich + FEve + FDis + genus +
                    (1|IBRAsub/planting/pl.p) + (1|species), data=plant.cent,
                    control = lmerControl(optimizer="bobyqa"))

summary(genera.int.m)
write.csv(summary(genera.int.m)$coefficients, "./Outputs/genus (int only) model coefs.csv")

#                       Genus model (ints & slopes) ####

genera.m <- lmer(mass.per.year ~ (aridity.index + solar.annual + 
                                     SLA + WD + MH + log.age +
                                     density + (intra.sp.prop * plot.area) + 
                                     n.rarerich + FEve + FDis ) * genus +
                  (1|IBRAsub/planting/pl.p) + (1|species), data=plant.cent,
                 control = lmerControl(optimizer="bobyqa"))
r.squaredGLMM(genera.m)
summary(genera.m)

write.csv(summary(genera.m)$coefficients, "./Outputs/genus model coefs.csv")

#         Base random models ####
#                       Generate random data ####

set.seed(97531)
if(!"random data.rds" %in% list.files("./Data")){
  
  data.list<-plant.cent
  
  generate.random.data<-function(data.list, iterations, 
                                 constrain.scale = "plot"){
  
    
    lapply(1:iterations, function(n){
      
    print(n)
    
    # make sure rownames are unique and integers
    rownames(data.list) <- 1:dim(data.list)[1]
    
    # shuffle rownames within each plot
    if(constrain.scale == "plot"){
    random.assorting <- unlist(tapply(rownames(data.list), data.list$pl.p,
                               function(x){
                                 
                                 x[sample(1:length(x), length(x))]
                                 
                               }))
    random.assorting <- as.numeric(random.assorting)
    }
    
    if(constrain.scale == "unconstrain"){
      
      random.assorting = sample(1:nrow(data.list), nrow(data.list))
      
    }
    
    # re-assign growth rates based on shuffled rownames
    data.list$mass.per.year <- data.list$mass.per.year[random.assorting]
    
    return(data.list)
    
    })
      
    }
    

  random.data<-generate.random.data(plant.cent, 1000, 
                                    constrain.scale="unconstrain")
  
  saveRDS(random.data, "./Data/random data unconstrained.rds")
} else {random.data<-readRDS("./Data/random data unconstrained.rds")}
#                       Run models ####

# So returning 1000 model objects results in a large object size. Instead, I'm going to run
# a null model, genus-less model and genus model for each subset, compute anova comparisons, and
# return model coefficients and AIC etc.
random.data.sub <- lapply(random.data, function(x){
  
  x[,c("mass.per.year", "IBRAsub", "planting", "pl.p", "species", "genus",
       "aridity.index", "solar.annual", "SLA", "WD", "MH", "log.age", "density",
       "intra.sp.prop", "plot.area", "n.rarerich", "FEve", "FDis")]
  
})

if(!"random model coefs.rds" %in% list.files("./Outputs")){
library(parallel)
cl <- makeCluster(4)
setDefaultCluster(cl)
# export data and load libraries in cluster
clusterExport(varlist=c("random.data.sub"), envir=environment())
clusterCall(cl, "library", "lme4", character.only = TRUE)
clusterCall(cl, "library", "MuMIn", character.only = TRUE)
options(na.action = "na.fail")
random.model.list <- parLapply(cl, 1:length(random.data.sub), function(random.subset){
  
  print(random.subset)
  
  random.subset <- random.data.sub[[random.subset]]
  
  # run null model
  temp.null <- lmer(mass.per.year ~ 1 +
                      (1|IBRAsub/planting/pl.p) + (1|species), data=random.subset)
  
  temp.genless <- lmer(mass.per.year ~  aridity.index + solar.annual + 
                                        SLA + WD + MH + log.age +
                                        density + (intra.sp.prop * plot.area) + 
                                        n.rarerich + FEve + FDis +
                      (1|IBRAsub/planting/pl.p) + (1|species), data=random.subset)
  
  temp.gen <- lmer(mass.per.year ~ (aridity.index + solar.annual + 
                                    SLA + WD + MH + log.age +
                                    density + (intra.sp.prop * plot.area) + 
                                    n.rarerich + FEve + FDis ) * genus +
                     (1|IBRAsub/planting/pl.p) + (1|species), data=random.subset)
              
  gen.coefs <- summary(temp.gen)$coefficients[,1]
  
  comp.table <- data.frame(AIC = c(AIC(temp.null), 
                                   AIC(temp.genless), 
                                   AIC(temp.gen)),
                           loglik = c(logLik(temp.null), logLik(temp.genless), logLik(temp.gen)),
                           chi.null = c(NA, 
                                        anova(temp.null, temp.genless)$Chisq[2],
                                        anova(temp.null, temp.gen)$Chisq[2]),
                           chi.genless = c(NA, NA,
                                           anova(temp.null, temp.gen)$Chisq[2]),
                           P.chi.null = c(NA, 
                                          anova(temp.null, temp.genless)[2,8],
                                          anova(temp.null, temp.gen)[2,8]),
                           P.chi.genless = c(NA, 
                                             NA,
                                             anova(temp.genless, temp.gen)[2,8]),
                           marg.r2 = c(r.squaredGLMM(temp.null)[1],
                                       r.squaredGLMM(temp.genless)[1],
                                       r.squaredGLMM(temp.gen)[1]),
                           cond.r2 = c(r.squaredGLMM(temp.null)[2],
                                       r.squaredGLMM(temp.genless)[2],
                                       r.squaredGLMM(temp.gen)[2]))
                           
  rownames(comp.table) = c("null", "genless", "gen")
  
  return(list(gen.coefs, comp.table))
  
})
stopCluster(cl=cl)

saveRDS(random.model.list, "./Outputs/random model coefs.rds")
} else {
random.model.list <- readRDS("./Outputs/random model coefs.rds")
}

random.coefs <- do.call("rbind", lapply(random.model.list, function(x){x[[1]]}))
random.table <- array(unlist(lapply(random.model.list, function(x){x[[2]]})),
                      dim=c(dim(random.model.list[[1]][[2]])[1],
                            dim(random.model.list[[1]][[2]])[2],
                            length(random.model.list)),
                      dimnames = list(rownames(random.model.list[[1]][[2]]),
                                      colnames(random.model.list[[1]][[2]]),
                                      1:length(random.model.list)))

# For comparison, we can plot differences in AIC in random models,
# And we can rejig our Fig 3 to a simpler version that doesn't need the model selection bit.

#         Aridity models ####

# Similar to previous modelling process, let's add on moisture availability interactions and
# compare genus-less and genera models, then model randoms

#                       Genus-less aridity model ####

genless.m.arid <- lmer(mass.per.year ~ (aridity.index + solar.annual + 
                                       SLA + WD + MH + log.age +
                                       density + (intra.sp.prop * plot.area) + 
                                       n.rarerich + FEve + FDis) * aridity.index +
                    (1|IBRAsub/planting/pl.p) + (1|species), data=plant.cent,
                    control = lmerControl(optimizer="bobyqa"))

write.csv(summary(genless.m.arid)$coefficients, "./Outputs/aridity genus-less model coefs.csv")

#                       Genus aridity model (ints) ####

genera.int.arid <- lmer(mass.per.year ~ (aridity.index + solar.annual + 
                                          SLA + WD + MH + log.age +
                                          density + (intra.sp.prop * plot.area) + 
                                          n.rarerich + FEve + FDis + genus) * aridity.index +
                         (1|IBRAsub/planting/pl.p) + (1|species), data=plant.cent,
                       control = lmerControl(optimizer="bobyqa"))

write.csv(summary(genless.m.arid)$coefficients, "./Outputs/aridity genus-less model coefs.csv")

#                       Genera aridity model ####

genera.m.arid <- lmer(mass.per.year ~ ((aridity.index + solar.annual + 
                                        SLA + WD + MH + log.age +
                                        density + (intra.sp.prop * plot.area) + 
                                        n.rarerich + FEve + FDis ) * genus) * aridity.index +
                   (1|IBRAsub/planting/pl.p) + (1|species), data=plant.cent)

summary(genera.m.arid)
anova(genless.m.arid, genera.m.arid)

write.csv(summary(genera.m.arid)$coefficients, "./Outputs/aridity genus model coefs.csv")

#         Aridity random models ####

# Similar to the before process, I'm going to run the model and return the values. I'm also
# going to make predictions at particular values of moisture availability
arid.quantile<-quantile(plant.cent$aridity.index, c(0.10,0.90))

if(!"random arid models coefs.rds" %in% list.files("./Outputs")){
library(parallel)
cl <- makeCluster(4)
setDefaultCluster(cl)
# export data and load libraries in cluster
clusterExport(varlist=c("random.data.sub", "inter.slope.calc", 
                        "predictors", "arid.quantile"), envir=environment())
clusterCall(cl, "library", "lme4", character.only = TRUE)
clusterCall(cl, "library", "MuMIn", character.only = TRUE)
options(na.action = "na.fail")
arid.random.list <- parLapply(cl, 1:length(random.data.sub), function(random.subset){
  
  print(random.subset)
  
  random.subset <- random.data.sub[[random.subset]]
  
  temp.null <- lmer(mass.per.year ~ 1 +
                   (1|IBRAsub/planting/pl.p) + (1|species), data=random.subset)
  
  temp.genless <- lmer(mass.per.year ~ (aridity.index + solar.annual + 
                                          SLA + WD + MH + log.age +
                                          density + (intra.sp.prop * plot.area) + 
                                          n.rarerich + FEve + FDis) * aridity.index +
                         (1|IBRAsub/planting/pl.p) + (1|species), data=random.subset)
  
  temp.gen <- lmer(mass.per.year ~ ((aridity.index + solar.annual + 
                                            SLA + WD + MH + log.age +
                                            density + (intra.sp.prop * plot.area) + 
                                            n.rarerich + FEve + FDis ) * genus) * aridity.index +
                          (1|IBRAsub/planting/pl.p) + (1|species), data=random.subset)
  
  gen.coefs <- summary(temp.gen)$coefficients[,1]
  
  comp.table <- data.frame(AIC = c(AIC(temp.null), 
                                   AIC(temp.genless), 
                                   AIC(temp.gen)),
                           loglik = c(logLik(temp.null), logLik(temp.genless), logLik(temp.gen)),
                           chi.null = c(NA, 
                                        anova(temp.null, temp.genless)$Chisq[2],
                                        anova(temp.null, temp.gen)$Chisq[2]),
                           chi.genless = c(NA, NA,
                                           anova(temp.genless, temp.gen)$Chisq[2]),
                           P.chi.null = c(NA, 
                                          anova(temp.null, temp.genless)[2,8],
                                          anova(temp.null, temp.gen)[2,8]),
                           P.chi.genless = c(NA, 
                                             NA,
                                             anova(temp.genless, temp.gen)[2,8]),
                           marg.r2 = c(r.squaredGLMM(temp.null)[1],
                                       r.squaredGLMM(temp.genless)[1],
                                       r.squaredGLMM(temp.gen)[1]),
                           cond.r2 = c(r.squaredGLMM(temp.null)[2],
                                       r.squaredGLMM(temp.genless)[2],
                                       r.squaredGLMM(temp.gen)[2]))
  
  rownames(comp.table) = c("null", "genless", "gen")
  
  # now predict slopes at particular values of aridity, just for Acacia
  
  inter.slopes <- sapply(predictors[predictors %in% colnames(temp.gen@frame) &
                    predictors != "aridity.index"], function(x){
           
    inter.slope.calc(model = temp.gen, var=x, genus = "Acacia")
                    })
  
  return(list(gen.coefs, comp.table, inter.slopes))

})
stopCluster(cl=cl)

saveRDS(arid.random.list, "./Outputs/random arid models coefs.rds")
} else {
arid.random.list <- readRDS("./Outputs/random arid models coefs.rds")
}

arid.random.coefs <- do.call("rbind", lapply(arid.random.list, function(x){x[[1]]}))
arid.random.table <- array(unlist(lapply(arid.random.list, function(x){x[[2]]})),
                      dim=c(dim(arid.random.list[[1]][[2]])[1],
                            dim(arid.random.list[[1]][[2]])[2],
                            length(arid.random.list)),
                      dimnames = list(rownames(arid.random.list[[1]][[2]]),
                                      colnames(arid.random.list[[1]][[2]]),
                                      1:length(arid.random.list)))

arid.inter.table <-array(unlist(lapply(arid.random.list, function(x){x[[3]]})),
                         dim=c(dim(arid.random.list[[1]][[3]])[1],
                               dim(arid.random.list[[1]][[3]])[2],
                               length(arid.random.list)),
                         dimnames = list(rownames(arid.random.list[[1]][[3]]),
                                         colnames(arid.random.list[[1]][[3]]),
                                         1:length(arid.random.list)))

#                       Model comparison ####

lapply(list(null.m, genless.m, genera.int.m, genera.m, 
            genless.m.arid, genera.int.arid, genera.m.arid), AIC)

lapply(list(null.m, genless.m, genera.int.m, genera.m, 
            genless.m.arid, genera.int.arid, genera.m.arid), r.squaredGLMM)

var.table <- do.call("rbind", lapply(list(null.m, genless.m, genera.int.m, genera.m, 
                         genless.m.arid, genera.int.arid, genera.m.arid),
                    function(m){
                      
                    temp <- VarCorr(m, comp="Variance")
                    temp <- c(unlist(temp), attr(temp, "sc"))
                    names(temp)[length(temp)] <- "residual"
                    return(temp)
                    }))


comp.table <- rbind(anova(genless.m, genera.int.m)[,c("Chisq", "Pr(>Chisq)")],
                    anova(genless.m, genera.m)[,c("Chisq", "Pr(>Chisq)")])
                    
comp.table.arid <- rbind(anova(genless.m.arid, genera.int.arid)[,c("Chisq", "Pr(>Chisq)")],
                    anova(genless.m.arid, genera.m.arid)[,c("Chisq", "Pr(>Chisq)")])

comp.table.genus <- rbind(anova(genera.int.m, genera.m)[,c("Chisq", "Pr(>Chisq)")],
                          anova(genera.int.m, genera.int.arid)[,c("Chisq", "Pr(>Chisq)")],
                          anova(genera.m, genera.m.arid)[,c("Chisq", "Pr(>Chisq)")],
                          anova(genera.int.arid, genera.m.arid)[,c("Chisq", "Pr(>Chisq)")])


#         Raw partial regression plots ####

summary(genera.m.arid)

plant.cent.mm <- genera.m.arid@frame
head(plant.cent.mm)

library(merTools)

slopes.pred.list <- lapply(1:11, function(n){
  print(n)
all.predictors <- c("solar.annual", "log.age", "plot.area", 
                    "SLA", "WD", "MH", 
                    "density", "intra.sp.prop",
                    "n.rarerich", "FEve", "FDis")
target.pred <- all.predictors[n]

pred.df <- as.data.frame(matrix(0, ncol=length(all.predictors)+6, nrow=200))
colnames(pred.df) = c("aridity.index", "genus", all.predictors, 
                      "pl.p", "planting", "IBRAsub", "species")
pred.df[,(ncol(pred.df)-3):ncol(pred.df)] = "a"

# set up 4 dfs, two for each genus and two for wet and dry

arid.test <- quantile(plant.cent.mm$aridity.index, probs=c(0.1,0.9))

euc.pred <- pred.df
euc.pred$genus <- factor("Eucalyptus", levels=levels(plant.cent.mm$genus))
euc.pred[,target.pred] = seq(min(plant.cent.mm[plant.cent.mm$genus=="Eucalyptus",
                                               target.pred]),
                             max(plant.cent.mm[plant.cent.mm$genus=="Eucalyptus",
                                               target.pred]), len=200)
euc.pred.dry <- euc.pred
euc.pred.wet <- euc.pred
euc.pred.dry$aridity.index = arid.test[1]
euc.pred.wet$aridity.index = arid.test[2]

acac.pred <- pred.df
acac.pred$genus <- factor("Acacia", levels=levels(plant.cent.mm$genus))
acac.pred[,target.pred] = seq(min(plant.cent.mm[plant.cent.mm$genus=="Acacia",
                                               target.pred]),
                             max(plant.cent.mm[plant.cent.mm$genus=="Acacia",
                                               target.pred]), len=200)
acac.pred.dry <- acac.pred
acac.pred.wet <- acac.pred
acac.pred.dry$aridity.index = arid.test[1]
acac.pred.wet$aridity.index = arid.test[2]

pred.list <- list(acac.pred.dry, acac.pred.wet,
                  euc.pred.dry, euc.pred.wet)

pred.list <- lapply(pred.list, function(x){
  
  temp.pred <- cbind(x[,target.pred],
                  predictInterval(genera.m.arid,
                  newdata = x,
                  n.sims = 3999,
                  which="fixed",
                  level=0.95,
                  include.resid.var=FALSE,
                  type="linear.prediction"))
  colnames(temp.pred)[1] = target.pred

  # what were the scaling for the target pred?
  target.scale <- scale(plant.uncent[,target.pred])
  
  temp.pred$raw.pred <- (temp.pred[,target.pred] * attr(target.scale, "scaled:scale")) +
                                                   attr(target.scale, "scaled:center")
  
  return(temp.pred)
    
})

return(pred.list)
})

cols <- c(rgb(0.616,0.569,0), # dark Acacia
          rgb(0.157,0.31,0.184), # dark Euc
          rgb(0.996,0.941,0.332, 0.7), # light Acacia
          rgb(0.38,0.525,0.404, 0.4)) # light Euc

pdf("./Plots/partial regression slopes.pdf", height=14, width=6, useDingbats = FALSE)

par(mfrow=c(11,4), mar=c(3,0,0,0), oma=c(0.5,3.5,4,1), 
    las=1, ps=10, mgp=c(3,0.1,0), las=1, tcl=-0.25)

sapply(1:length(slopes.pred.list), function(n){
  
  temp.pred.list <- slopes.pred.list[[n]]

  target.pred = colnames(temp.pred.list[[1]])[1]
  
  y.axis=ifelse(n == 6, TRUE, FALSE)
  top.labels = ifelse(n==1, TRUE, FALSE)
  
sapply(1:length(temp.pred.list), function(n){

  all.predictors <- c("solar.annual", "log.age", "plot.area", 
                      "SLA", "WD", "MH", 
                      "density", "intra.sp.prop",
                      "n.rarerich", "FEve", "FDis")
  preds<-temp.pred.list[[n]]
  
  arid.quants <- quantile(plant.uncent$aridity.index, probs=c(0.5, 0.5))
  if(n %in% c(2,4)){
    target.arid = plant.uncent$aridity.index <= arid.quants[1]
  } else {
    target.arid = plant.uncent$aridity.index >= arid.quants[2]
  }
    
  target.points <- plant.uncent[(plant.uncent$genus == c("Acacia", "Eucalyptus")[(n+1) %/% 2]) &
                                target.arid,
                                c(target.pred, "mass.per.year")]
  
  if(n %in% 1:2){temp.col = cols[c(1,3)]} else {temp.col = cols[c(2,4)]}
  
  plot(x=NULL,y=NULL, xlim=range(plant.uncent[,target.pred]),
       ylim=range(plant.uncent$mass.per.year), axes=FALSE, xlab="", ylab="")
  
  if(n==1){
    axis(side=2, mgp=c(3,0.5,0), at=seq(-1,1.5,0.5))
    if(y.axis){
    mtext(side=2, las=0, line=1.75, text=expression("Annual biomass increment (kg kg"^-1*" year"^-1*")"))
    }
  } else {axis(side=2, labels=NA, at=seq(-1,1.5,0.5))}
  
  axis(side=1, at=pretty(plant.uncent[,target.pred], n=4))
  
  if(n==3){
    mtext(side=1, line=1.5, at=par("usr")[1], 
          text= c(expression("Solar radiation (ln(MJ m"^-2*" day"^-1*"))"),
                  "Planting age (ln(years))",
                  expression("Plot area (m"^2*")"),
                  expression("Specific leaf area (mm"^2* " g"^-1*")"),
                  expression("Wood density (g cm"^-3*")"),
                  "Maximum height (ln(m))",
                  expression("Neighbor density ("*sqrt("neighbors m"^-2)*")"),
                  "Proportion of intraspecific neighbors",
                  "Rarefied species richness",
                  "Functional evenness",
                  "Functional dispersion")[which(all.predictors == target.pred)],
          cex=0.9)
  }
  
  if(top.labels){
  mtext(side=3, line=0.1, 
        text=c(expression(italic("Acacia")*" (dry)"), 
               expression(italic("Acacia")*" (mesic)"),
               expression(italic("Eucalyptus")*" (dry)"),
               expression(italic("Eucalyptus")*" (mesic)"))[n],
        font=2, col=temp.col[1])
  }
  
  points(y=target.points$mass.per.year, x=target.points[,1], pch=16, col=temp.col[1])
  rect(xleft=par("usr")[1], xright=par("usr")[2], ybottom=par("usr")[3], ytop=par("usr")[4],
       border=NA, col=rgb(1,1,1,0.65))
  polygon(x=c(preds$raw.pred, rev(preds$raw.pred)),
          y=c(preds$upr, rev(preds$lwr)), col=temp.col[2],
          border=temp.col[1],
          lty=ifelse(n %in% c(1,3), "31", "solid"))
  
  lines(x=preds$raw.pred, y=preds$fit, lwd=2, col=temp.col[1])
  box()
})
})
dev.off()


  






# ####
# PLOTS ####
#         Comparison of plot slopes (Fig 3 re-done) ####
#                       Data prep ####
#                                 Marginal slopes (all data) ####

# To plot marginal slopes, we need the standard error of marginal slopes for Acacia and
# Eucalyptus. As Acacia is the reference level, we can use the slopes and SEs for the
# covariate partial regression coefficients. For Eucalyptus, we can either compute them by hand
# or run a second model with Eucalyptus as the reference level, and use those. I'm going to 
# use the second method.

acac.slopes <- as.data.frame(summary(genera.m)$coefficients)
acac.slopes <- acac.slopes[!grepl("Eucalyptus", rownames(acac.slopes)),]
acac.slopes$genus <- "Acacia"
acac.slopes$var = rownames(acac.slopes)

# re-run model for Eucalyptus
euc.cent <- plant.cent
euc.cent$genus <- relevel(as.factor(euc.cent$genus), ref="Eucalyptus")
genera.m.euc <- update(genera.m, data=euc.cent)

euc.slopes <- as.data.frame(summary(genera.m.euc)$coefficients)
euc.slopes <- euc.slopes[!grepl("Acacia", rownames(euc.slopes)),]
euc.slopes$genus <- "Eucalyptus"
euc.slopes$var = rownames(euc.slopes)

marg.slopes <- rbind(acac.slopes, euc.slopes)
marg.slopes$upperCI = marg.slopes$Estimate + 1.96 * marg.slopes$`Std. Error`
marg.slopes$lowerCI = marg.slopes$Estimate - 1.96 * marg.slopes$`Std. Error`

#                                 Marginal slopes (random data) ####

# We need to do a similar process to our random models, re-running each using Eucalyptus as the
# reference level

acac.slopes.rand <- random.coefs[, !grepl("Eucalyptus", colnames(random.coefs))]

if(!"euc slopes.rds" %in% list.files("./Outputs")){
euc.random.coefs <- do.call("rbind", lapply(1:length(random.data), function(x){
  print(x)
  x<-random.data[[x]]
  x$genus <- relevel(as.factor(x$genus), ref = "Eucalyptus")
  
  temp.gen <- lmer(mass.per.year ~ (aridity.index + solar.annual + 
                                      SLA + WD + MH + log.age +
                                      density + (intra.sp.prop * plot.area) + 
                                      n.rarerich + FEve + FDis ) * genus +
                     (1|IBRAsub/planting/pl.p) + (1|species), data=x)
  
  return(summary(temp.gen)$coefficients[,"Estimate"])
}))
saveRDS(euc.random.coefs, "./Outputs/euc slopes.rds")
} else {
euc.random.coefs <- readRDS("./Outputs/euc slopes.rds")
}

euc.slopes.rand <- euc.random.coefs[, !grepl("Acacia", colnames(euc.random.coefs))]

#                       Plot ####
pdf(paste0("./Plots/tree level coefficient plot - 2018",
           Sys.Date(), ".pdf"), height=4, width=6,
    useDingbats=FALSE)

framemat<-rbind(c(0.125,0.5,0.1,0.975), # focal effects
                c(0.625, 0.95, 0.53, 0.77), # density
                c(0.625, 0.95, 0.1, 0.5), # diversity
                c(0.625, 0.95, 0.75, 1)) # arrows/legend
                
split.screen(framemat)

#                       Screens ####
#                               Focal coefs ####

screen(1) # focal coefficients

par(mar=c(0,0,0,0), ps=8, tcl=-0.25, mgp=c(3,0.15,0))

focal.effects <- c("aridity.index", "solar.annual", "log.age", "plot.area",
                   "SLA","WD","MH")

# extract coefs to plot, and re-order
plot.frame <- marg.slopes[marg.slopes$var %in% focal.effects,]
plot.frame <- plot.frame[order(match(plot.frame$var, focal.effects),
                               plot.frame$genus),]
plot.frame$order <- (length(focal.effects)+1) - 
                    match(plot.frame$var, focal.effects) + 
                    c(0.15, -0.15)[as.factor(plot.frame$genus)]

plot(x=plot.frame$Estimate, y=plot.frame$order,
     xlim=c(-0.11, 0.03), ylim=c(0.75,7.35),
     axes=FALSE, xlab="", ylab="", type="n")

# set up frame
abline(v=0, lty="dashed", col="grey60")
axis(side=1, at=c(-0.1,-0.05,0), mgp=c(3,0,0))
axis(side=1, at=seq(-0.11,0.03, 0.01), tcl = -0.125, labels=NA)

axis(side=2, at=(length(focal.effects)):1, las=1, mgp=c(3,0.3,0),
     labels=c("Moisture\navailability"," Solar\nradiation", "Planting\nage",
              "Plot area", "Specific\nleaf area", "Wood\ndensity", "Maximum\nheight"))
mtext(side=1, text="Standardized slope estimate", line=0.75)
text(x=relative.axis.point(0.015, "x"),
     y=relative.axis.point(0.97, "y"), adj=0,
     labels="(a) Focal effects", font=2)

# plot violins
width=0.125

sapply(focal.effects, function(x){
  
  orders <- plot.frame$order[plot.frame$var == x]
  
  temp <- acac.slopes.rand[,colnames(acac.slopes.rand)==x]
  temp.s <- qnorm(c(0.025, 0.5, 0.975), mean=mean(temp), sd=sd(temp))
  
  rect(xleft=temp.s[1], xright=temp.s[3],
       ybottom=orders[1]-width, ytop=orders[1]+width,
       col=rgb(0.996,0.941,0.332, 0.7))
  # segments(x0=temp.s[1], x1=temp.s[1],
  #          y0=orders[1]-width +0.0125, 
  #          y1=orders[1]+width -0.0125,
  #          lend=2, lwd=2)
 
  temp <- euc.slopes.rand[,colnames(euc.slopes.rand)==x]
  temp.s <- qnorm(c(0.025, 0.5, 0.975), mean=mean(temp), sd=sd(temp))
  
  rect(xleft=temp.s[1], xright=temp.s[3],
       ybottom=orders[2]-width, ytop=orders[2]+width,
       col=rgb(0.38,0.525,0.404, 0.4))
  # segments(x0=temp.s[1], x1=temp.s[1],
  #          y0=orders[2]-width +0.0125, 
  #          y1=orders[2]+width -0.0125,
  #          lend=2, lwd=2)
  
  })      

# plot points
points(x=plot.frame$Estimate, 
       y=plot.frame$order + c(0, -0.025)[as.factor(plot.frame$genus)],
       col=c(rgb(0.616,0.569,0), rgb(0.157,0.31,0.184))[as.factor(plot.frame$genus)], 
       pch=c(16,17)[as.factor(plot.frame$genus)])

box()

close.screen(1)

#                               Density coefs ####

screen(2) # density coefficients 
par(mar=c(0,0,0,0), ps=8, tcl=-0.25, mgp=c(3,0.15,0))

focal.effects <- c("density", "intra.sp.prop")

# extract coefs to plot, and re-order
plot.frame <- marg.slopes[marg.slopes$var %in% focal.effects,]
plot.frame <- plot.frame[order(match(plot.frame$var, focal.effects),
                               plot.frame$genus),]
plot.frame$order <- (length(focal.effects)+1) - 
  match(plot.frame$var, focal.effects) + 
  c(0.15, -0.15)[as.factor(plot.frame$genus)]

plot(x=plot.frame$Estimate, y=plot.frame$order,
     xlim=c(-0.06, 0.02), ylim=c(0.6,2.75),
     axes=FALSE, xlab="", ylab="", type="n")

# set up frame
abline(v=0, lty="dashed", col="grey60")
axis(side=1, mgp=c(3,0,0), labels=NA)
axis(side=1, at=seq(-0.06,0.02,0.01), tcl = -0.125, labels=NA)

axis(side=2, at=(length(focal.effects)):1, las=1, mgp=c(3,0.3,0),
     labels=c("Neighbor\ndensity","Proportion\nintraspecific\nneighbours"))

rect(xleft=par("usr")[1],
     xright=par("usr")[2],
     ybottom=relative.axis.point(0.85, "y"),
     ytop=par("usr")[4],
     col="white", border=NA)

text(x=relative.axis.point(0.01, "x"),
     y=relative.axis.point(0.9, "y"), adj=0,
     labels="(b) Neighbor density", font=2)

# plot violins
width=0.125
sapply(focal.effects, function(x){
  
  orders <- plot.frame$order[plot.frame$var == x]
  
  temp <- acac.slopes.rand[,colnames(acac.slopes.rand)==x]
  temp.s <- qnorm(c(0.025, 0.5, 0.975), mean=mean(temp), sd=sd(temp))
  
  rect(xleft=temp.s[1], xright=temp.s[3],
       ybottom=orders[1]-width, ytop=orders[1]+width,
       col=rgb(0.996,0.941,0.332, 0.7))
  # segments(x0=temp.s[1], x1=temp.s[1],
  #          y0=orders[1]-width +0.0125, 
  #          y1=orders[1]+width -0.0125,
  #          lend=2, lwd=2)
  
  temp <- euc.slopes.rand[,colnames(euc.slopes.rand)==x]
  temp.s <- qnorm(c(0.025, 0.5, 0.975), mean=mean(temp), sd=sd(temp))
  
  rect(xleft=temp.s[1], xright=temp.s[3],
       ybottom=orders[2]-width, ytop=orders[2]+width,
       col=rgb(0.38,0.525,0.404, 0.4))
  # segments(x0=temp.s[1], x1=temp.s[1],
  #          y0=orders[2]-width +0.0125, 
  #          y1=orders[2]+width -0.0125,
  #          lend=2, lwd=2)
  
})    
# plot points
points(x=plot.frame$Estimate, 
       y=plot.frame$order + c(0, -0.025)[as.factor(plot.frame$genus)],
       col=c(rgb(0.616,0.569,0), rgb(0.157,0.31,0.184))[as.factor(plot.frame$genus)], 
       pch=c(16,17)[as.factor(plot.frame$genus)])

box()
close.screen(2)

#                               Diversity coefs ####

screen(3) # Diversity coefficients 
par(mar=c(0,0,0,0), ps=8, tcl=-0.25, mgp=c(3,0.15,0))

focal.effects <- c("n.rarerich", "FEve", "FDis")

# extract coefs to plot, and re-order
plot.frame <- marg.slopes[marg.slopes$var %in% focal.effects,]
plot.frame <- plot.frame[order(match(plot.frame$var, focal.effects),
                               plot.frame$genus),]
plot.frame$order <- (length(focal.effects)+1) - 
  match(plot.frame$var, focal.effects) + 
  c(0.15, -0.15)[as.factor(plot.frame$genus)]

plot(x=plot.frame$Estimate, y=plot.frame$order,
     xlim=c(-0.06, 0.02), ylim=c(0.75, 3.5),
     axes=FALSE, xlab="", ylab="", type="n")

# set up frame
abline(v=0, lty="dashed", col="grey60")
axis(side=1, mgp=c(3,0,0))
axis(side=1, at=seq(-0.06,0.02,0.01), tcl = -0.125, labels=NA)

axis(side=2, at=(length(focal.effects)):1, las=1, mgp=c(3,0.3,0),
     labels=c("Species\nrichness","Functional\nevenness",
              "Functional\ndispersion"))
mtext(side=1, text="Standardized slope estimate", line=0.75)

rect(xleft=par("usr")[1],
     xright=par("usr")[2],
     ybottom=relative.axis.point(0.9, "y"),
     ytop=par("usr")[4],
     col="white", border=NA)

text(x=relative.axis.point(0.01, "x"),
     y=relative.axis.point(0.94, "y"), adj=0,
     labels="(c) Neighborhood diversity", font=2)

# plot violins
width=0.12
sapply(focal.effects, function(x){
  
  orders <- plot.frame$order[plot.frame$var == x]
  
  temp <- acac.slopes.rand[,colnames(acac.slopes.rand)==x]
  temp.s <- qnorm(c(0.025, 0.5, 0.975), mean=mean(temp), sd=sd(temp))
  
  rect(xleft=temp.s[1], xright=temp.s[3],
       ybottom=orders[1]-width, ytop=orders[1]+width,
       col=rgb(0.996,0.941,0.332, 0.7))
  # segments(x0=temp.s[1], x1=temp.s[1],
  #          y0=orders[1]-width +0.0125, 
  #          y1=orders[1]+width -0.0125,
  #          lend=2, lwd=2)
  
  temp <- euc.slopes.rand[,colnames(euc.slopes.rand)==x]
  temp.s <- qnorm(c(0.025, 0.5, 0.975), mean=mean(temp), sd=sd(temp))
  
  rect(xleft=temp.s[1], xright=temp.s[3],
       ybottom=orders[2]-width, ytop=orders[2]+width,
       col=rgb(0.38,0.525,0.404, 0.4))
  # segments(x0=temp.s[1], x1=temp.s[1],
  #          y0=orders[2]-width +0.0125, 
  #          y1=orders[2]+width -0.0125,
  #          lend=2, lwd=2)
  
})       

# plot points
points(x=plot.frame$Estimate, 
       y=plot.frame$order + c(0, -0.025)[as.factor(plot.frame$genus)],
       col=c(rgb(0.616,0.569,0), rgb(0.157,0.31,0.184))[as.factor(plot.frame$genus)], 
       pch=c(16,17)[as.factor(plot.frame$genus)])

box()
close.screen(3)

#                               Legend ####

screen(4)
par(mar=c(0,0,0,0), ps=8, tck=-0.02, mgp=c(3,0.15,0))
plot.new()

acac.y<- 0.525
euc.y<-  0.325
all.x<-  0.485
random.x<- 0.785

label.y<- 0.7
label.x<- 0.355


rect(ybottom=acac.y - 0.06,
     ytop = acac.y + 0.06,
     xleft=random.x - 0.085,
     xright = random.x + 0.085,
     col=rgb(0.996,0.941,0.332, 0.7))

rect(ybottom=euc.y - 0.06,
     ytop = euc.y + 0.06,
     xleft=random.x - 0.085,
     xright = random.x + 0.085,
     col=rgb(0.38,0.525,0.404, 0.4))

points(y=c(acac.y, euc.y), x=rep(all.x, 2),
       pch=c(16,17),
       col=c(rgb(0.616,0.569,0), rgb(0.157,0.31,0.184)))

text(y=c(acac.y, euc.y), x=rep(label.x, 2),
     labels=c("Acacia", "Eucalyptus"), font=3, adj=1)

text(x=c(all.x, random.x), y= label.y - c(0,0.01),
     labels=c("Observed", "Null"))

rect(xleft=0.5 - 0.465, 
     xright=0.5 + 0.465, 
     ybottom= 0.5 - 0.275, 
     ytop= 0.5 + 0.275,
     lwd=0.5)

close.screen(4)

dev.off()
#         Aridity interactions plot ####
#                             Data prep - All data model ####

# We need to extract a value for each of our interesting interacting variables, at particular
# levels of moisture availability -perhaps the 10th and 90th quantile.

main.inter.slopes <- data.frame(t(sapply(predictors[predictors %in% colnames(genera.m.arid@frame) &
                                        predictors != "aridity.index"],
                           function(x){
                             slopes <- inter.slope.calc(model = genera.m.arid,
                                     var = x,
                                     genus = "Acacia")
                             
                             SEs <- interaction.slope.se.fun(model = genera.m.arid,
                                                      var = x)
                             
                             return(cbind(slopes, SEs))
                           })))

colnames(main.inter.slopes) <- c("mean10", "mean90", "se10", "se90")
main.inter.slopes$genus = "Acacia"
main.inter.slopes$var = rownames(main.inter.slopes)

# Now do the same for Eucalyptus, based on a re-leveled model
genera.m.arid.euc <- update(genera.m.arid, data=euc.cent)

euc.inter.slopes <- data.frame(t(sapply(predictors[predictors %in% colnames(genera.m.arid@frame) &
                                                     predictors != "aridity.index"],
                                        function(x){
                                          slopes <- inter.slope.calc(model = genera.m.arid.euc,
                                                                     var = x,
                                                                     genus = "Eucalyptus")
                                          
                                          SEs <- interaction.slope.se.fun(model = genera.m.arid.euc,
                                                                          var = x)
                                          
                                          if(x %in% c("intra.sp.prop", "plot.area")){
                                            slopes <- slopes[c(1,3)]
                                          }
                                          
                                          return(cbind(slopes, SEs))
                                        })))

colnames(euc.inter.slopes) <- c("mean10", "mean90", "se10", "se90")
euc.inter.slopes$genus = "Eucalyptus"
euc.inter.slopes$var = rownames(euc.inter.slopes)

arid.marg.slopes <- rbind(main.inter.slopes, euc.inter.slopes)
arid.marg.slopes$upperCI10 <- arid.marg.slopes$mean10 + 1.96 * arid.marg.slopes$se10
arid.marg.slopes$lowerCI10 <- arid.marg.slopes$mean10 - 1.96 * arid.marg.slopes$se10
arid.marg.slopes$upperCI90 <- arid.marg.slopes$mean90 + 1.96 * arid.marg.slopes$se90
arid.marg.slopes$lowerCI90 <- arid.marg.slopes$mean90 - 1.96 * arid.marg.slopes$se90

#                             Data prep - Random subset models ####

# we only need slopes for our random models, not standard errors, which makes the process much
# easier. We already computed the 10% and 90% slopes for Acacia, we just need to re-run the random
# models and extract the Eucalyptus slopes
random.data.sub <- lapply(random.data, function(x){
  
  x[,c("mass.per.year", "IBRAsub", "planting", "pl.p", "species", "genus",
       "aridity.index", "solar.annual", "SLA", "WD", "MH", "log.age", "density",
       "intra.sp.prop", "plot.area", "n.rarerich", "FEve", "FDis")]
  
})

if(!"euc.arid.slopes.rds" %in% list.files("./Outputs")){
# temporary coefficient generator
library(parallel)
cl <- makeCluster(4)
setDefaultCluster(cl)
# export data and load libraries in cluster
clusterExport(varlist=c("random.data.sub", "inter.slope.calc", 
                        "predictors", "arid.quantile"), envir=environment())
clusterCall(cl, "library", "lme4", character.only = TRUE)
clusterCall(cl, "library", "MuMIn", character.only = TRUE)
options(na.action = "na.fail")
random.coefs <- parLapply(cl, 1:length(random.data.sub), function(x){
  
  print(x)
  x <- random.data.sub[[x]]
  
  acac.gen <- lmer(mass.per.year ~ ((aridity.index + solar.annual + 
                                       SLA + WD + MH + log.age +
                                       density + (intra.sp.prop * plot.area) + 
                                       n.rarerich + FEve + FDis ) * genus) * aridity.index +
                     (1|IBRAsub/planting/pl.p) + (1|species), data=x)
  
  
  x$genus <- relevel(as.factor(x$genus), ref = "Eucalyptus")
  
  temp.gen <- lmer(mass.per.year ~ ((aridity.index + solar.annual + 
                                       SLA + WD + MH + log.age +
                                       density + (intra.sp.prop * plot.area) + 
                                       n.rarerich + FEve + FDis ) * genus) * aridity.index +
                     (1|IBRAsub/planting/pl.p) + (1|species), data=x)
  
  acac.inter.slopes <- sapply(predictors[predictors %in% colnames(acac.gen@frame) &
                                      predictors != "aridity.index"], function(x){
                                        
                                        inter.slope.calc(model = acac.gen, var=x, genus = "Acacia")
                                        
                                      })
  
  
  euc.inter.slopes <- sapply(predictors[predictors %in% colnames(temp.gen@frame) &
                                      predictors != "aridity.index"], function(x){
                                        
                                        inter.slope.calc(model = temp.gen, var=x, genus = "Eucalyptus")
                                        
                                      })
  
  return(list(acac.inter.slopes, euc.inter.slopes))
})
stopCluster(cl)

saveRDS(random.coefs, "./Outputs/combined.random.coefs.rds")

acac.random.coefs <- lapply(random.coefs, function(x){x[[1]]})
euc.random.coefs <- lapply(random.coefs, function(x){x[[2]]})

euc.random.coefs <- lapply(1:length(random.data.sub), function(x){
  
  print(x)
  x <- random.data.sub[[x]]
  x$genus <- relevel(x$genus, ref = "Eucalyptus")
  
  temp.gen <- lmer(mass.per.year ~ ((aridity.index + solar.annual + 
                                       SLA + WD + MH + log.age +
                                       density + (intra.sp.prop * plot.area) + 
                                       n.rarerich + FEve + FDis ) * genus) * aridity.index +
                     (1|IBRAsub/planting/pl.p) + (1|species), data=x)
  
  inter.slopes <- sapply(predictors[predictors %in% colnames(temp.gen@frame) &
                                      predictors != "aridity.index"], function(x){
                                        
                                        inter.slope.calc(model = temp.gen, var=x, genus = "Eucalyptus")
                                        
                                      })
  
  return(inter.slopes)
})

saveRDS(euc.random.coefs, "./Outputs/euc.arid.slopes.rds")
} else {
acac.random.coefs <- readRDS("./Outputs/combined.random.coefs.rds")
acac.random.coefs <- lapply(acac.random.coefs, function(x){x[[1]]})
euc.random.coefs <- readRDS("./Outputs/euc.arid.slopes.rds")
}

arid.inter.table <-array(unlist(acac.random.coefs),
                             dim=c(dim(acac.random.coefs[[1]])[1],
                                   dim(acac.random.coefs[[1]])[2],
                                   length(acac.random.coefs)),
                             dimnames = list(rownames(acac.random.coefs[[1]]),
                                             colnames(acac.random.coefs[[1]]),
                                             1:length(acac.random.coefs)))

arid.inter.table.euc <-array(unlist(euc.random.coefs),
                         dim=c(dim(euc.random.coefs[[1]])[1],
                               dim(euc.random.coefs[[1]])[2],
                               length(euc.random.coefs)),
                         dimnames = list(rownames(euc.random.coefs[[1]]),
                                         colnames(euc.random.coefs[[1]]),
                                         1:length(euc.random.coefs)))


#                       PLOT ####

pdf(paste0("./Plots/tree level aridity interaction - 2018",
           Sys.Date(), ".pdf"), height=6, width=6,
    useDingbats=FALSE)

# framemat<-rbind(c(0.09,0.34,0.1,0.975), # focal effects
#                 c(0.44, 0.7, 0.53, 0.77), # density
#                 c(0.44, 0.7, 0.1, 0.5), # diversity
#                 c(0.44, 0.7, 0.75, 1), # legend
#                 c(0.77,0.99,0.6,0.975), # Acacia callout
#                 c(0.77,0.99,0.1, 0.475)) # Eucalyptus callout

framemat<-rbind(c(0.1,0.48,0.1,0.975), # focal effects
                c(0.60, 0.98, 0.565, 0.85), # density
                c(0.60, 0.98, 0.1, 0.55), # diversity
                c(0.60, 0.98, 0.85, 1)) # legend

split.screen(framemat)

#                               Focal coefs ####

screen(1) # focal coefficients

par(mar=c(0,0,0,0), ps=8, tcl=-0.25, mgp=c(3,0.55,0))

focal.effects <- c("solar.annual", "log.age", "plot.area",
                   "SLA","WD","MH")

# extract coefs to plot, and re-order
plot.frame <- arid.marg.slopes[arid.marg.slopes$var %in% focal.effects,]
plot.frame <- plot.frame[order(match(plot.frame$var, focal.effects),
                               plot.frame$genus),]

# offset genera points
plot.frame$order <- (length(focal.effects)+1) - 
  match(plot.frame$var, focal.effects) + 
  c(0.175, -0.175)[as.factor(plot.frame$genus)]

plot(x=plot.frame$mean10, y=plot.frame$order,
     xlim=c(-0.16, 0.05), ylim=c(0.75,6.35),
     axes=FALSE, xlab="", ylab="", type="n")

# set up frame
abline(v=0, lty="dashed", col="grey60")
axis(side=1, mgp=c(3,0,0))
axis(side=1, at=seq(-0.125,0.025,0.025), tcl= -0.125, labels=NA)
#axis(side=1, mgp=c(3,0,0), at=c(-0.1,0,0.1))

axis(side=2, at=(length(focal.effects)):1, las=1,
     labels=c("Solar\nradiation", "Planting\nage",
              "Plot area", "Specific\nleaf area", "Wood\ndensity", "Maximum\nheight"))
mtext(side=1, text="Standardized slope estimate", line=0.75)
text(x=relative.axis.point(0.015, "x"),
     y=relative.axis.point(0.98, "y"), adj=0,
     labels="(a) Focal effects", font=2)

# plot violins
wetdryoffset <- 0.09
width<-0.075
sapply(focal.effects, function(x){

  orders <- plot.frame$order[plot.frame$var == x]

  temp <- arid.inter.table["90%",x,]
  temp.s <- qnorm(c(0.025, 0.5, 0.975), mean=mean(temp), sd=sd(temp))
  
  # Acacia violin plot - wet
  rect(xleft=temp.s[1], 
       xright=temp.s[3],
       ybottom=orders[1]-width - wetdryoffset , 
       ytop=orders[1]+width  - wetdryoffset,
       col=rgb(0.996,0.941,0.332, 0.7))

  temp <- arid.inter.table["10%",x,]
  temp.s <- qnorm(c(0.025, 0.5, 0.975), mean=mean(temp), sd=sd(temp))
  
  # Acacia violin plot - dry
  rect(xleft=temp.s[1], xright=temp.s[3],
       ybottom=orders[1]-width + wetdryoffset, 
       ytop=orders[1]+width + wetdryoffset,
       col=rgb(0.996,0.941,0.332, 0.4), lty="21")

  # Eucalyptus violin plot - wet
  temp <- arid.inter.table.euc["90%",x,]  
  temp.s <- qnorm(c(0.025, 0.5, 0.975), mean=mean(temp), sd=sd(temp))
  rect(xleft=temp.s[1], 
       xright=temp.s[3],
       ybottom=orders[2]-width - wetdryoffset, ytop=orders[2]+width  - wetdryoffset,
       col=rgb(0.38,0.525,0.404, 0.4))

  # Eucalyptus violin plot - dry
  temp <- arid.inter.table.euc["10%",x,]
  temp.s <- qnorm(c(0.025, 0.5, 0.975), mean=mean(temp), sd=sd(temp))
  
  # Acacia violin plot - dry
  rect(xleft=temp.s[1], xright=temp.s[3],
       ybottom=orders[2]-width + wetdryoffset, 
       ytop=orders[2]+width + wetdryoffset,
       col=rgb(0.38,0.525,0.404, 0.2), lty="21")

})

# # callout boxes for (e) and (f)
# x.limits <- unlist(plot.frame[plot.frame$var=="SLA" & plot.frame$genus == "Acacia", 1:2])
# y.centre <- plot.frame[plot.frame$var=="SLA" & plot.frame$genus == "Acacia", "order"]
# 
# rect(ybottom=y.centre + 0.135, 
#      ytop=y.centre - 0.135, 
#      xleft= min(x.limits) - 0.015, 
#      xright= max(x.limits) + 0.016, 
#      lwd=0.5, col="white")
# text(x=min(x.limits) - 0.015, y=y.centre, pos=2, labels="(e)", offset=0.1)
# 
# # MH Euc
# x.limits <- unlist(plot.frame[plot.frame$var=="MH" & plot.frame$genus == "Eucalyptus", 1:2])
# y.centre <- plot.frame[plot.frame$var=="MH" & plot.frame$genus == "Eucalyptus", "order"]
# 
# rect(ybottom=y.centre + 0.175, 
#      ytop=y.centre - 0.175, 
#      xleft= min(x.limits) - 0.015, 
#      xright= max(x.limits) + 0.015, 
#      lwd=0.5, col="white")
# text(x=max(x.limits) + 0.015, y=y.centre, pos=4, labels="(f)", offset=0.1)

# plot points - wet
points(x=plot.frame$mean90, 
       y=plot.frame$order - wetdryoffset - c(0, 0.025)[as.factor(plot.frame$genus)],
       col=c(rgb(0.616,0.569,0), rgb(0.157,0.31,0.184))[as.factor(plot.frame$genus)], 
       pch=c(16,17)[as.factor(plot.frame$genus)])

# plot points - dry
points(x=plot.frame$mean10, 
       y=plot.frame$order + wetdryoffset - c(0, 0.025)[as.factor(plot.frame$genus)],
       col=c(rgb(0.616,0.569,0), rgb(0.157,0.31,0.184))[as.factor(plot.frame$genus)], 
       bg=rgb(1,1,1,1),
       pch=c(21,24)[as.factor(plot.frame$genus)])

box()

close.screen(1)

#                               Density coefs ####

screen(2) # density coefficients
par(mar=c(0,0,0,0), ps=8, tcl=-0.25, mgp=c(3,0.15,0))

focal.effects <- c("density", "intra.sp.prop")

# extract coefs to plot, and re-order
plot.frame <- arid.marg.slopes[arid.marg.slopes$var %in% focal.effects,]
plot.frame <- plot.frame[order(match(plot.frame$var, focal.effects),
                               plot.frame$genus),]
plot.frame$order <- (length(focal.effects)+1) -
  match(plot.frame$var, focal.effects) +
  c(0.18, -0.18)[as.factor(plot.frame$genus)]

plot(x=plot.frame$mean10, y=plot.frame$order,
     xlim=c(-0.09, 0.03), ylim=c(0.6,2.75),
     axes=FALSE, xlab="", ylab="", type="n")

# set up frame
abline(v=0, lty="dashed", col="grey60")
axis(side=1, labels=NA)

axis(side=2, at=(length(focal.effects)):1, las=1, mgp=c(3,0.3,0),
     labels=c("Neighbor\ndensity","Proportion\nintraspecific\nneighbours"))

rect(xleft=par("usr")[1],
     xright=par("usr")[2],
     ybottom=relative.axis.point(0.85, "y"),
     ytop=par("usr")[4],
     col="white", border=NA)

text(x=relative.axis.point(0.01, "x"),
     y=relative.axis.point(0.95, "y"), adj=0,
     labels="(b) Neighbor density", font=2)

# plot violins
wetdryoffset <- 0.09
width<-0.08
sapply(focal.effects, function(x){
  
  orders <- plot.frame$order[plot.frame$var == x]
  
  
  temp <- arid.inter.table["90%",x,]
  temp.s <- qnorm(c(0.025, 0.5, 0.975), mean=mean(temp), sd=sd(temp))
  
  # Acacia violin plot - wet
  rect(xleft=temp.s[1], 
       xright=temp.s[3],
       ybottom=orders[1]-width - wetdryoffset , 
       ytop=orders[1]+width  - wetdryoffset,
       col=rgb(0.996,0.941,0.332, 0.7))
  
  temp <- arid.inter.table["10%",x,]
  temp.s <- qnorm(c(0.025, 0.5, 0.975), mean=mean(temp), sd=sd(temp))
  
  # Acacia violin plot - dry
  rect(xleft=temp.s[1], xright=temp.s[3],
       ybottom=orders[1]-width + wetdryoffset, 
       ytop=orders[1]+width + wetdryoffset,
       col=rgb(0.996,0.941,0.332, 0.4), lty="21")
  
  # Eucalyptus violin plot - wet
  temp <- arid.inter.table.euc["90%",x,]  
  temp.s <- qnorm(c(0.025, 0.5, 0.975), mean=mean(temp), sd=sd(temp))
  rect(xleft=temp.s[1], 
       xright=temp.s[3],
       ybottom=orders[2]-width - wetdryoffset, ytop=orders[2]+width  - wetdryoffset,
       col=rgb(0.38,0.525,0.404, 0.4))
  
  # Eucalyptus violin plot - dry
  temp <- arid.inter.table.euc["10%",x,]
  temp.s <- qnorm(c(0.025, 0.5, 0.975), mean=mean(temp), sd=sd(temp))
  
  # Acacia violin plot - dry
  rect(xleft=temp.s[1], xright=temp.s[3],
       ybottom=orders[2]-width + wetdryoffset, 
       ytop=orders[2]+width + wetdryoffset,
       col=rgb(0.38,0.525,0.404, 0.2), lty="21")
  
})

# plot points - wet
points(x=plot.frame$mean90, 
       y=plot.frame$order - wetdryoffset - c(0, 0.025)[as.factor(plot.frame$genus)],
       col=c(rgb(0.616,0.569,0), rgb(0.157,0.31,0.184))[as.factor(plot.frame$genus)], 
       pch=c(16,17)[as.factor(plot.frame$genus)])

# plot points - dry
points(x=plot.frame$mean10, 
       y=plot.frame$order + wetdryoffset - c(0, 0.025)[as.factor(plot.frame$genus)],
       col=c(rgb(0.616,0.569,0), rgb(0.157,0.31,0.184))[as.factor(plot.frame$genus)], 
       bg=rgb(1,1,1,1),
       pch=c(21,24)[as.factor(plot.frame$genus)])

box()
close.screen(2)
# 
#                               Diversity coefs ####

screen(3) # Diversity coefficients
par(mar=c(0,0,0,0), ps=8, tcl = -0.25, mgp=c(3,0.15,0))

focal.effects <- c("n.rarerich", "FEve", "FDis")

# extract coefs to plot, and re-order
plot.frame <- arid.marg.slopes[arid.marg.slopes$var %in% focal.effects,]
plot.frame <- plot.frame[order(match(plot.frame$var, focal.effects),
                               plot.frame$genus),]
plot.frame$order <- (length(focal.effects)+1) -
  match(plot.frame$var, focal.effects) +
  c(0.175, -0.175)[as.factor(plot.frame$genus)]

plot(x=plot.frame$mean10, y=plot.frame$order,
     xlim=c(-0.08, 0.03), ylim=c(0.75, 3.5),
     axes=FALSE, xlab="", ylab="", type="n")

# set up frame
abline(v=0, lty="dashed", col="grey60")
axis(side=1, mgp=c(3,0,0))
#axis(side=1, at=c(-0.15, -0.1,-0.05,0,0.05,0.1), labels=NA)
#axis(side=1, mgp=c(3,0,0), at=c(-0.1,0,0.1))

axis(side=2, at=(length(focal.effects)):1, las=1, mgp=c(3,0.3,0),
     labels=c("Species\nrichness","Functional\nevenness",
              "Functional\ndispersion"))
mtext(side=1, text="Standardized slope estimate", line=0.75)

rect(xleft=par("usr")[1],
     xright=par("usr")[2],
     ybottom=relative.axis.point(0.9, "y"),
     ytop=par("usr")[4],
     col="white", border=NA)

text(x=relative.axis.point(0.01, "x"),
     y=relative.axis.point(0.96, "y"), adj=0,
     labels="(c) Neighborhood diversity", font=2)

# plot violins
wetdryoffset <- 0.09
width<-0.065
sapply(focal.effects, function(x){
  
  orders <- plot.frame$order[plot.frame$var == x]
  
  
  temp <- arid.inter.table["90%",x,]
  temp.s <- qnorm(c(0.025, 0.5, 0.975), mean=mean(temp), sd=sd(temp))
  
  # Acacia violin plot - wet
  rect(xleft=temp.s[1], 
       xright=temp.s[3],
       ybottom=orders[1]-width - wetdryoffset , 
       ytop=orders[1]+width  - wetdryoffset,
       col=rgb(0.996,0.941,0.332, 0.7))
  
  temp <- arid.inter.table["10%",x,]
  temp.s <- qnorm(c(0.025, 0.5, 0.975), mean=mean(temp), sd=sd(temp))
  
  # Acacia violin plot - dry
  rect(xleft=temp.s[1], xright=temp.s[3],
       ybottom=orders[1]-width + wetdryoffset, 
       ytop=orders[1]+width + wetdryoffset,
       col=rgb(0.996,0.941,0.332, 0.7), lty="21")
  
  # Eucalyptus violin plot - wet
  temp <- arid.inter.table.euc["90%",x,]  
  temp.s <- qnorm(c(0.025, 0.5, 0.975), mean=mean(temp), sd=sd(temp))
  rect(xleft=temp.s[1], 
       xright=temp.s[3],
       ybottom=orders[2]-width - wetdryoffset, ytop=orders[2]+width  - wetdryoffset,
       col=rgb(0.38,0.525,0.404, 0.4))
  
  # Eucalyptus violin plot - dry
  temp <- arid.inter.table.euc["10%",x,]
  temp.s <- qnorm(c(0.025, 0.5, 0.975), mean=mean(temp), sd=sd(temp))
  
  # Acacia violin plot - dry
  rect(xleft=temp.s[1], xright=temp.s[3],
       ybottom=orders[2]-width + wetdryoffset, 
       ytop=orders[2]+width + wetdryoffset,
       col=rgb(0.38,0.525,0.404, 0.4), lty="21")
  
})

# plot points - wet
points(x=plot.frame$mean90, 
       y=plot.frame$order - wetdryoffset - c(0, 0.02)[as.factor(plot.frame$genus)],
       col=c(rgb(0.616,0.569,0), rgb(0.157,0.31,0.184))[as.factor(plot.frame$genus)], 
       pch=c(16,17)[as.factor(plot.frame$genus)])

# plot points - dry
points(x=plot.frame$mean10, 
       y=plot.frame$order + wetdryoffset - c(0, 0.02)[as.factor(plot.frame$genus)],
       col=c(rgb(0.616,0.569,0), rgb(0.157,0.31,0.184))[as.factor(plot.frame$genus)], 
       bg=rgb(1,1,1,1),
       pch=c(21,24)[as.factor(plot.frame$genus)])

box()
close.screen(3)

#                               Legend ####

screen(4)
par(mar=c(0,0,0,0), ps=8, tck=-0.02, mgp=c(3,0.15,0))
plot.new()

acac.y<- 0.4
euc.y <- 0.2

moist.center <- 0.4
moist.exp.x <-  moist.center + 0.1
moist.obs.x <-  moist.center - 0.05

dry.center <- 0.75
dry.exp.x <- dry.center + 0.1
dry.obs.x <- dry.center - 0.05

label.y<- 0.6
label.x<- 0.28

box.width = 0.08
box.height = 0.065

points(y=c(acac.y, euc.y), x=rep(moist.obs.x, 2),
       pch=c(16,17),
       col=c(rgb(0.616,0.569,0), rgb(0.157,0.31,0.184)))

points(y=c(acac.y, euc.y), x=rep(dry.obs.x, 2),
       pch=c(21,24),
       bg="white",
       col=c(rgb(0.616,0.569,0), rgb(0.157,0.31,0.184)))

rect(xleft=moist.exp.x - box.width,
     xright = moist.exp.x + box.width,
     ytop = acac.y + box.height,
     ybottom = acac.y - box.height,
     col=rgb(0.996,0.941,0.332, 0.7))

rect(xleft=moist.exp.x - box.width,
     xright = moist.exp.x + box.width,
     ytop = euc.y + box.height,
     ybottom = euc.y - box.height,
     col=rgb(0.38,0.525,0.404, 0.4))

rect(xleft=dry.exp.x - box.width,
     xright = dry.exp.x + box.width,
     ytop = acac.y + box.height,
     ybottom = acac.y - box.height,
     col=rgb(0.996,0.941,0.332, 0.4), lty="21")

rect(xleft=dry.exp.x - box.width,
     xright = dry.exp.x + box.width,
     ytop = euc.y + box.height,
     ybottom = euc.y - box.height,
     col=rgb(0.38,0.525,0.404, 0.2), lty="21")

text(y=c(acac.y, euc.y), x=rep(label.x, 2),
     labels=c("Acacia", "Eucalyptus"), font=3, adj=1)

text(x=c(moist.exp.x, 
         moist.obs.x, 
         dry.exp.x, 
         dry.obs.x), 
     y=label.y,
     labels=c("Null", "Obs", "Null", "Obs"), adj=0.5)

text(x=c(moist.center + 0.5*box.width, 
         dry.center + 0.5*box.width), y=label.y+0.15,
     labels=c("Mesic", "Dry"), font=2)

rect(xleft=0.5 - 0.495, 
     xright=0.5 + 0.495, 
     ybottom= 0.5 - 0.45, 
     ytop= 0.5 + 0.375,
     lwd=0.5)

close.screen(4)

dev.off()
#         Biomass diameter correlation plot ####

cor(log(mixed.plant.small$dbh), log(mixed.plant.small$tot.ag.mass), use="complete.obs")

summary(plant.uncent$tot.ag.mass)

dbh.cor <- cor(log(plant.uncent$dbh),
               log(plant.uncent$tot.ag.mass), use="complete.obs")

d10.cor <- cor(log(plant.uncent$d10.d50),
               log(plant.uncent$tot.ag.mass), use="complete.obs")

pdf("./Plots/diameter biomass comparison.pdf", height=3, width=6, useDingbats=FALSE)

par(mfrow=c(1,2), mar=c(0,0,0,0), oma=c(2.5,3,1,1), mgp=c(3,0.5,0), ps=8, tck= -0.025, las=1)

plot(y=log(plant.uncent$tot.ag.mass),
     x=(log(plant.uncent$dbh)),
     col=rgb(0.2,0.2,0.2,0.025), pch=16, cex=0.8,
     axes=FALSE, xlab="", ylab="")
text(x = relative.axis.point(0.01, "x"),
     y = relative.axis.point(0.95, "y"),
     adj = 0, labels = paste0("(a) DBH: r = ", round(dbh.cor, 3)),
     font = 2)
axis(side=2, at=log(c(0.1, 1,10,100,1000,10000)), labels=c(0.1, 1,10,100,1000,10000))
axis(side=2, at=log(c(seq(0.1,1,0.1), seq(1,10,1), seq(10,100,10),seq(100,1000,100),
                      seq(1000,10000,1000))), labels=NA, tck= -0.01)
axis(side=1, at=log(c(0.1, 1,10,50,100,1000,10000)), 
     labels=c(0.1, 1,10,50,100,1000,10000), mgp=c(3,0,0))
axis(side=1, at=log(c(seq(0.1,1,0.1), seq(1,10,1), seq(10,100,10),seq(100,1000,100),
                      seq(1000,10000,1000))), labels=NA, tck= -0.01)
mtext(side=1, line=1, text="Diameter at breast height (cm)")
mtext(side=2, line=2, text="Above-ground biomass (kg)", las=0)
box()


plot(y=log(plant.uncent$tot.ag.mass),
     x=(log(plant.uncent$d10.d50)),
     col=rgb(0.2,0.2,0.2,0.025), pch=16, cex=0.8,
     axes=FALSE, xlab="", ylab="")

text(x = relative.axis.point(0.01, "x"),
     y = relative.axis.point(0.95, "y"),
     adj = 0, labels = paste0("(b) D10/D50: r = ", round(d10.cor, 3)),
     font = 2)
box()

axis(side=2, at=log(c(0.1, 1,10,100,1000,10000)), labels=NA)
axis(side=2, at=log(c(seq(0.1,1,0.1), seq(1,10,1), seq(10,100,10),seq(100,1000,100),
                      seq(1000,10000,1000))), labels=NA, tck= -0.01)
axis(side=1, at=log(c(0.1, 1,10,50,100,1000,10000)), 
     labels=c(0.1, 1,10,50,100,1000,10000), mgp=c(3,0,0))
axis(side=1, at=log(c(seq(0.1,1,0.1), seq(1,10,1), seq(10,100,10),seq(100,1000,100),
                      seq(1000,10000,1000))), labels=NA, tck= -0.01)

mtext(side=1, line=1, text="Diameter at 10cm or 50cm aboveground (cm)")
dev.off()

#         Site map and trait plots ####

#                         Data prep ####

plant.uncent <- read.csv("./Data/uncentered data for modelling.csv")
ausmap<-readShapeSpatial("/home/timothy/University files - offline/PhD - offline/Shape files/Australia outline/AUS_adm/AUS_adm1.shp")

data(countriesLow)

ausmap <- countriesLow[countriesLow$SOVEREIGNT == "Australia",]

plot(ausmap)
aridity <- raster("/home/timothy/University files - offline/PhD - offline/Shape files/aridity.extent.tif")
proj4string(ausmap) <- proj4string(aridity)
log.aridity <- log(aridity)

mixed.planting<-mixed.site[!duplicated(mixed.site$planting) &
                             mixed.site$pl.p %in% plant.cent$pl.p,]
mixed.planting<-merge(mixed.planting, mean.planting.dens, all.x=TRUE)

#                         Acacia occurrence raster ####

# species occurrence records
acac.occur <- read.csv('./Data/Acacia_occurrence_records.csv')
acac.occur <- acac.occur[,c("Matched.Scientific.Name",
                            "Latitude...processed",
                            "Longitude...processed")]
colnames(acac.occur) <- c("species", "lat", "long")
acac.occur <- acac.occur[complete.cases(acac.occur),]
acac.occur <- droplevels(acac.occur[acac.occur$species %in% plant.uncent$species,])
acac.occur <- acac.occur[acac.occur$lat > -45 & acac.occur$lat < -9 &
                           acac.occur$long > 110 & acac.occur$long < 155, ]

# remnove points that aren't inside Australia polygon
acac.coord <- acac.occur[,c("long", "lat")]
coordinates(acac.coord) <- c("long", "lat")
proj4string(acac.coord) <- proj4string(ausmap)
acac.occur <- acac.occur[!is.na(over(acac.coord, ausmap)$ScaleRank),]

# turn into raster of species richness at particular coordinates
raster.res <- c(1,1)
acac.occur$lat.bin <- cut(acac.occur$lat, breaks = seq(min(acac.occur$lat),
                                                       max(acac.occur$lat),
                                                       raster.res[1]))
acac.occur$long.bin <- cut(acac.occur$long, breaks = seq(min(acac.occur$long),
                                                         max(acac.occur$long),
                                                         raster.res[2]))


acac.occur.table <- as.array(table(acac.occur$lat.bin, acac.occur$long.bin, acac.occur$species))
acac.occur.table <- ifelse(acac.occur.table > 0, 1 , 0)
acac.occur.table <- apply(acac.occur.table, c(1,2), sum)
acac.occur.table <- acac.occur.table[dim(acac.occur.table)[1]:1,]

acac.occur.rast <- raster(acac.occur.table,
                          ymn = as.numeric(substr(levels(acac.occur$lat.bin)[1],2,
                                                  regexpr(",", levels(acac.occur$lat.bin)[1])-1)),
                          ymx = as.numeric(substr(rev(levels(acac.occur$lat.bin))[1],
                                                  regexpr(",", rev(levels(acac.occur$lat.bin))[1])+1,
                                                  nchar(rev(levels(acac.occur$lat.bin))[1])-1)),
                          xmn = as.numeric(substr(levels(acac.occur$long.bin)[1],2,
                                                  regexpr(",", levels(acac.occur$long.bin)[1])-1)),
                          xmx = as.numeric(substr(rev(levels(acac.occur$long.bin))[1],
                                                  regexpr(",", rev(levels(acac.occur$long.bin))[1])+1,
                                                  nchar(rev(levels(acac.occur$long.bin))[1])-1)),
                          crs = proj4string(aridity))

acac.occur.rast[acac.occur.rast == 0] = NA
acac.occur.rast[acac.occur.rast > 1 & acac.occur.rast <=5] = 2
acac.occur.rast[acac.occur.rast > 5 & acac.occur.rast <=10] = 3
acac.occur.rast[acac.occur.rast > 10 & acac.occur.rast <=20] = 4
acac.occur.rast[acac.occur.rast > 20] = 5

maxValue(acac.occur.rast)

dim(acac.ramp)
acac.ramp <- data.frame(vals=c(0,1,5,10,20,50),
                        ramp=colorRampPalette(c("white", rgb(0.616,0.569,0)))(6))

#                         Eucalyptus occurrence raster ####

# species occurrence records
euc.occur <- read.csv('./Data/Eucalyptus_occurrence_records.csv')
euc.occur <- euc.occur[,c("Matched.Scientific.Name",
                          "Latitude...processed",
                          "Longitude...processed")]

colnames(euc.occur) <- c("species", "lat", "long")
euc.occur <- euc.occur[complete.cases(euc.occur),]
euc.occur <- droplevels(euc.occur[euc.occur$species %in% plant.uncent$species,])

euc.coord <- euc.occur[,c("long", "lat")]
coordinates(euc.coord) <- c("long", "lat")
proj4string(euc.coord) <- proj4string(ausmap)
euc.occur <- euc.occur[!is.na(over(euc.coord, ausmap)$ScaleRank),]

# turn into raster of species richness at particular coordinates

euc.occur$lat.bin <- cut(euc.occur$lat, breaks = seq(min(euc.occur$lat),
                                                     max(euc.occur$lat),
                                                     raster.res[1]))
euc.occur$long.bin <- cut(euc.occur$long, breaks = seq(min(euc.occur$long),
                                                       max(euc.occur$long),
                                                       raster.res[2]))

euc.occur.table <- as.array(table(euc.occur$lat.bin, euc.occur$long.bin, euc.occur$species))
euc.occur.table <- ifelse(euc.occur.table > 0, 1 , 0)
euc.occur.table <- apply(euc.occur.table, c(1,2), sum)
euc.occur.table <- euc.occur.table[dim(euc.occur.table)[1]:1,]

euc.occur.rast <- raster(euc.occur.table,
                         ymn = as.numeric(substr(levels(euc.occur$lat.bin)[1],2,
                                                 regexpr(",", levels(euc.occur$lat.bin)[1])-1)),
                         ymx = as.numeric(substr(rev(levels(euc.occur$lat.bin))[1],
                                                 regexpr(",", rev(levels(euc.occur$lat.bin))[1])+1,
                                                 nchar(rev(levels(euc.occur$lat.bin))[1])-1)),
                         xmn = as.numeric(substr(levels(euc.occur$long.bin)[1],2,
                                                 regexpr(",", levels(euc.occur$long.bin)[1])-1)),
                         xmx = as.numeric(substr(rev(levels(euc.occur$long.bin))[1],
                                                 regexpr(",", rev(levels(euc.occur$long.bin))[1])+1,
                                                 nchar(rev(levels(euc.occur$long.bin))[1])-1)),
                         crs = proj4string(aridity))

euc.occur.rast[euc.occur.rast == 0] = NA
euc.occur.rast[euc.occur.rast > 1 & euc.occur.rast <=5] = 2
euc.occur.rast[euc.occur.rast > 5 & euc.occur.rast <=10] = 3
euc.occur.rast[euc.occur.rast > 10 & euc.occur.rast <=20] = 4
euc.occur.rast[euc.occur.rast > 20] = 5

euc.ramp <- data.frame(vals=c(0,1,5,10,20,50),
                       ramp=colorRampPalette(c("white", rgb(0.157,0.31,0.184)))(6))

#                           Plot setup ####

pdf(paste0("./Plots/combined map ", Sys.Date(), ".pdf"), height=6, width=5)


# split.screen(rbind(c(0.11, 0.99, 0.52, 0.96), # raster
#                    c(0.11, 0.99, 0.52, 0.96), # map outline
#                    
#                    c(0.25,0.55, 0.57,0.62), # aridity density plot
#                    
#                    c(0.11, 0.36, 0.19, 0.31), # Max height
#                    c(0.425, 0.675, 0.19, 0.31), # Wood dens
#                    c(0.74, 0.99, 0.19, 0.31), # SLA
#                    
#                    c(0.11, 0.99, 0.52, 0.96), # map box
#                    
#                    c(0.075, 0.475, 0.31, 0.51), # Acac occur map
#                    c(0.075, 0.475, 0.31, 0.51),
#                    
#                    c(0.61, 1, 0.31, 0.51), # Euc occur map
#                    c(0.61, 1, 0.31, 0.51),
#                    
#                    c(0.475, 0.625, 0.36, 0.46), # Occur legend
#                    c(0.11, 0.99, 0.31, 0.52)) # occur box


split.screen(rbind(c(0.11, 0.99, 0.51, 0.96), # raster
                   c(0.11, 0.99, 0.51, 0.96), # map outline
                   
                   c(0.25,0.55, 0.57,0.62), # aridity density plot
                   
                   c(0.11, 0.36, 0.06, 0.23), # Max height
                   c(0.425, 0.675, 0.06, 0.23), # Wood dens
                   c(0.74, 0.99, 0.06, 0.23), # SLA
                   
                   c(0.11, 0.99, 0.51, 0.96), # map box
                   
                   c(0.075, 0.475, 0.23, 0.51), # Acac occur map
                   c(0.075, 0.475, 0.23, 0.51),
                   
                   c(0.61, 1, 0.23, 0.51), # Euc occur map
                   c(0.61, 1, 0.23, 0.51),
                   
                   c(0.475, 0.625, 0.28, 0.42), # Occur legend
                   c(0.11, 0.99, 0.23, 0.51)) # occur box
             
) 

#                                     Map ####

screen(1)
par(mar=c(0,0,0,0), ps=8, tcl= -0.25, mgp=c(3,0.5,0), las=1)

arid.ramp <- data.frame(log.aridity = seq(minValue(log.aridity),
                                          maxValue(log.aridity), length.out=200),
                        ramp = colorRampPalette(c("white",rgb(0.9,0.9,1,1), rgb(0.15,0.15,1,1)))(200))

plot(x=NULL,y=NULL, xlim=c(118,148), ylim=c(-45,-9), axes=FALSE,
     xlab = "", ylab = "", asp=1)
plot(log.aridity, col=as.character(arid.ramp$ramp), 
     box=FALSE, add=TRUE, legend=FALSE, maxpixels=1000000, interpolate=TRUE)
mappar<-par("usr")
close.screen(1)

screen(2)
par(mar=c(0,0,0,0), ps=8, tcl= -0.25, mgp=c(3,0.5,0), las=1)
plot(x=NULL,y=NULL, xlim=c(118,148), ylim=c(-45,-9), axes=FALSE,
     ylab = "", xlab = "", asp=1)
plot(ausmap, lwd=1, axes=FALSE, border="grey40", add=TRUE)

text(x=relative.axis.point(0.01, "x"),
     y=relative.axis.point(0.965, "y"),
     labels="(a)", font=2, cex=1.25, adj=0)
axis(side=3, at=seq(110,150,10), labels=parse(text=paste(seq(110,150,10), "*degree~E", sep="")),
     mgp=c(3,0.2,0))
axis(side=2, at=seq(-45,-10,5), labels=parse(text=paste(seq(-45,-10,5), "*degree~S", sep="")),
     mgp=c(3,0.5,0))

points(mixed.planting$lat ~ mixed.planting$long, 
       cex=0.8, pch=21, bg="white")

close.screen(2)

#                                     Aridity density plot ####

screen(3)
par(mar=c(0,0,0,0), ps=8, tcl=-0.25, mgp=c(3,0.15,0), las=1)

arid<-plant.uncent$aridity.index[!duplicated(plant.uncent$pl.p)]

arid.dens <- density(arid, adjust=0.75)

plot(x=NULL, y=NULL, yaxs="i", xaxs="i",
     xlim=c(summary(arid.dens$x)[c(1)], 1.4),
     ylim=c(-0.05,max(arid.dens$y)+0.1),
     axes=FALSE, main="", xlab= "", ylab = "")

rect.dist <- arid.ramp[2,1] - arid.ramp[1,1]
for(i in 1:dim(arid.ramp)[1]){
  
  rect(xleft = exp(arid.ramp$log.aridity[i] - rect.dist),
       xright = exp(arid.ramp$log.aridity[i] + rect.dist),
       ybottom= par("usr")[3],
       ytop= par("usr")[4],
       border=NA,
       col=as.character(arid.ramp$ramp[i]))
  
}

polygon(x=c(arid.dens$x, 
            arid.dens$x[length(arid.dens$x)],
            par("usr")[2],
            par("usr")[2],
            arid.dens$x[1]),
        y=c(arid.dens$y,
            0,
            0,
            par("usr")[4],
            par("usr")[4]), 
        col="white",
        border="black")

box(col="white", lwd=2)
axis(side=1, mgp=c(3,0,0), line= -0.05)
axis(side=1, at=seq(200,1400,100), labels=NA, tcl=-0.25, cex=0.8, line = -0.0)
segments(x0=par("usr")[1], x1=par("usr")[2], y0=0, y1=0,
         lwd=1)
axis(side=2, at=c(0,1,2,3), mgp=c(3,0.5,0))
axis(side=2, at=c(0.5,1.5,2.5), tcl = -0.125, labels=NA)
mtext(side=2, line=1, text="Density", las=0)

mtext(side=1, line=0.65, text="Moisture availability (unitless)")
close.screen(3)

#                                     Trait density plots ####

euc.sp<-plant.uncent[plant.uncent$genus=="Eucalyptus" &
                       !duplicated(plant.uncent$species), ]
acac.sp<-plant.uncent[plant.uncent$genus=="Acacia" &
                        !duplicated(plant.uncent$species), ]

labels<-list(expression("Specific leaf area (mm"^2*" g"^-1*")"),
             expression("Wood density (g cm"^-3*")"),
             "Seed mass (mg)",
             "Maximum height (m)")

text.pos<-list(c(log(11), 0.18, log(4), 0.8),
               c(0.7, 2.5, 0.875, 0.4),
               c(log(5.5), 0.4, log(20), 0.5))
names(text.pos)<-c("SLA","WD","MH")

lma.lims <- c(14,1400)
sla.lims <- lma.lims*1000 # g/m2 -> mg/m2
sla.lims <- sla.lims/1e6 # mg/m2 -> mg/mm2
sla.lims <- 1/sla.lims # mg/mm2 -> mm2/mg

xlimits <- rbind(SLA = log(c(rev(sla.lims))),
                 WD = c(0.1,1.5),
                 MH = log(c(0.5, 110)))

sapply(3:1,
       function(n){
         
         screen(n+3)
         
         x <- c("SLA","WD","MH")[n]
         par(mar=c(0,0,0,0), ps=8, tcl=-0.25, mgp=c(3,0.5,0), las=1)
         
         if(x == "SLA"){
           
           euc.dens<-density(log(euc.sp[,x]), 
                             from=log(min(euc.sp[,x])), 
                             to=log(max(euc.sp[,x])))
           
           acac.dens<-density(log(acac.sp[,x]), 
                              from=log(min(acac.sp[,x])), 
                              to=log(max(acac.sp[,x])))
           
         } else {
           
           euc.dens<-density(euc.sp[,x], from=min(euc.sp[,x]), to=max(euc.sp[,x]))
           acac.dens<-density(acac.sp[,x], from=min(acac.sp[,x]), to=max(acac.sp[,x]))
           
         }
         
         y.max<-max(c(euc.dens$y, acac.dens$y))
         
         plot(x=NULL,y=NULL, 
              xlim=c(xlimits[n,1],
                     xlimits[n,2]), 
              ylim=c(0, y.max + 0.05*y.max), 
              yaxs="i", axes=FALSE, ylab="", xlab="")
         
         axis(side=2, at = pretty(par("usr")[3:4], 3))
         box()
         
         polygon(x=c(acac.dens$x, 
                     rev(acac.dens$x)[1],
                     acac.dens$x[1]),
                 y=c(acac.dens$y, 
                     par("usr")[3],
                     par("usr")[3]), 
                 col=rgb(0.996,0.941,0.332, 0.7),
                 border=rgb(0.616,0.569,0))
         
         polygon(x=c(euc.dens$x, 
                     rev(euc.dens$x)[1],
                     euc.dens$x[1]),
                 y=c(euc.dens$y, 
                     par("usr")[3],
                     par("usr")[3]), 
                 col=rgb(0.38,0.525,0.404, 0.4),
                 border=rgb(0.157,0.31,0.184))
         
         # Ac and Euc labels
         temp.text<-text.pos[[which(names(text.pos)==x)]]
         text(x=temp.text[c(1,3)],
              y=temp.text[c(2,4)],
              labels=c("A", "E"), font=2,
              col=c(rgb(0.616,0.569,0), rgb(0.157,0.31,0.184)))
         
         # SLA axes
         if(x == "SLA"){
           axis(side=1, at=log(c(0.1,1,10,50,100)), labels=c(0.1,1,10,50,100),
                tcl=-0.25, mgp=c(3,0,0))
           axis(side=1, at=log(c(seq(0.2,1,0.1),
                                 seq(2,10,1),
                                 seq(10,100,10))),
                labels=NA,
                tcl=-0.125)
           mtext(side=2, text="Density", line=1.5, las=0)
           
           text(x=relative.axis.point(0.03, "x"),
                y=relative.axis.point(0.9, "y"),
                labels="(c)", font=2, cex=1.25, adj=0)
         }
         
         # max height axes
         if(x == "MH"){
           axis(side=1, at=log(c(1,10,100)), labels=c(1,10,100),
                tcl=-0.25, mgp=c(3,0,0))
           axis(side=1, at=log(c(seq(2,10,1),
                                 seq(10,100,10))),
                labels=NA,
                tcl=-0.125)
           text(x=relative.axis.point(0.03, "x"),
                y=relative.axis.point(0.9, "y"),
                labels="(e)", font=2, cex=1.25, adj=0)
         }
         
         # other axes
         if(x %in% c("WD")){
           axis(side=1, tck=-0.04, mgp=c(3,0,0))
           text(x=relative.axis.point(0.03, "x"),
                y=relative.axis.point(0.9, "y"),
                labels="(d)", font=2, cex=1.25, adj=0)
         }
         
         mtext(side=1, line=0.9, text=labels[[which(c("SLA","WD","SM","MH") == x)]])
         
         close.screen(n+3)
         
       })

screen(7)
par(mar=c(0,0,0,0), ps=8, tcl=-0.25, mgp=c(3,0.15,0), las=1)
plot.new()
box()
close.screen(7)

#                                     Species occurrence rasters ####
screen(8)
par(mar=c(0,0,0,0), ps=8, tcl= -0.25, mgp=c(3,0.5,0), las=1)
plot(x=NULL,y=NULL, xlim=c(110,155), ylim=c(-45,-9), axes=FALSE,
     xlab = "", ylab = "", asp=1)
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
     border=NA, col=rgb(1,1,1,0.9))
plot(acac.occur.rast, axes=FALSE, add=TRUE, legend=FALSE,
     col=as.character(acac.ramp$ramp[2:6]))
close.screen(8)

screen(9)
par(mar=c(0,0,0,0), ps=8, tcl= -0.25, mgp=c(3,0.5,0), las=1)
plot(x=NULL,y=NULL, xlim=c(110,155), ylim=c(-45,-9), axes=FALSE,
     xlab = "", ylab = "", asp=1)
plot(ausmap, axes=FALSE, add=TRUE, lwd=0.5)

par(xpd=NA)
text(x=relative.axis.point(0.1, "x"),
     y=relative.axis.point(0.95, "y"),
     labels="(b)", font=2, cex=1.25, adj=0)
par(xpd=FALSE)

close.screen(9)

screen(10)
par(mar=c(0,0,0,0), ps=8, tcl= -0.25, mgp=c(3,0.5,0), las=1)
plot(x=NULL,y=NULL, xlim=c(110,155), ylim=c(-45,-9), axes=FALSE,
     xlab = "", ylab = "", asp=1)
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
     border=NA, col=rgb(1,1,1,0.9))
plot(euc.occur.rast, axes=FALSE, add=TRUE, legend=FALSE,
     col=as.character(euc.ramp$ramp[2:6]))
close.screen(10)

screen(11)
par(mar=c(0,0,0,0), ps=8, tcl= -0.25, mgp=c(3,0.5,0), las=1)
plot(x=NULL,y=NULL, xlim=c(110,155), ylim=c(-45,-9), axes=FALSE,
     xlab = "", ylab = "", asp=1, lwd=0.5)
plot(ausmap, axes=FALSE, add=TRUE)
close.screen(11)

screen(12)
par(mar=c(0,0,0,0), ps=8, tcl= -0.25, mgp=c(3,0.5,0), las=1)
y.sections <-seq(par("usr")[3], par("usr")[4], length.out=dim(acac.ramp)[1]+1)
y.dist <- 0.5* (y.sections[2] - y.sections[1])

plot.new()
rect(xleft=rep(0, dim(acac.ramp)[1]), 
     xright=rep(0.15, dim(acac.ramp)[1]),
     ybottom=rev(rev(y.sections)[-1]),
     ytop = y.sections[-1],
     col=as.character(acac.ramp$ramp))

rect(xleft=rep(0.85, dim(acac.ramp)[1]), 
     xright=rep(1, dim(acac.ramp)[1]),
     ybottom=rev(rev(y.sections)[-1]),
     ytop = y.sections[-1],
     col=as.character(euc.ramp$ramp))

segments(x0=0.15, x1=0.2,
         y0=y.sections[-1] - y.dist,
         y1 =y.sections[-1] - y.dist)

segments(x0=0.85, x1=0.8,
         y0=y.sections[-1] - y.dist,
         y1 =y.sections[-1] - y.dist)

text(x=0.5, y = y.sections[-1] - y.dist,
     labels = c("0", "1", "2 - 5", "6 - 10", "11 - 20", "> 20"))

mtext(side=3, at=0.075, text="A", font=2)
mtext(side=3, at=0.925, text="E", font=2)
mtext(side=3, at=0.5, text="# sp.", font=2)

close.screen(12)     

screen(13)
par(mar=c(0,0,0,0), ps=8, tcl= -0.25, mgp=c(3,0.5,0), las=1)
plot.new()
box()
close.screen(13)

dev.off()


#         Random models Likelihood ratio plot ####

pdf("./Plots/random likelihood comparison.pdf", height=2.4, width=5, useDingbats=FALSE)

par(mfrow=c(1,2), mar=c(0,0,0,0), oma=c(2.5,3,1,1), ps=8, 
    mgp=c(3,0.5,0), las=1, tck= -0.025)

xaxis <- log(c(1e-15, 1e-10, 1e-5, 1e-2, 0.1, 1))

temp <- density(log(random.table[3,6,]), to=0)
plot(x=NULL, y=NULL, xlim=c(log(min(arid.random.table[3,6,])), 2), ylim=c(0,0.1), xaxt="n", yaxs="i",
     xlab="", ylab="")
polygon(y=c(temp$y, 0),
        x=c(temp$x, rev(temp$x)[1]), col="grey70")

temp.non <- data.frame(x=temp$x[temp$x >= log(0.05)],
                       y=temp$y[temp$x >= log(0.05)])

polygon(y=c(0, temp.non$y, 0),
        x=c(log(0.05), temp.non$x, rev(temp.non$x)[1]), col="red")

axis(side=1, at = xaxis,
     labels= NA)
mtext(side=2, line=2, text="Density", las=0)
axis(side=1, at = log(as.numeric(paste0("1e-", 15:1))), labels=NA, tck = -0.01)
text(x=xaxis, y=relative.axis.point(-0.05, "y"), srt=55, adj=1,
     labels = c("1e-15", "1e-10", "1e-5", "0.01", "0.1", 1))
abline(v=log(0.05))
text(x=relative.axis.point(0.01, "x"),
     y=relative.axis.point(0.95, "y"), labels="(a)", font=2, adj=0)

sig.sum<-sum(random.table[3,6,]<=0.05) / 1000
text(x=log(0.05), y=relative.axis.point(0.95, "y"),
     labels=paste0(sig.sum*100, "%"), pos=2, offset=0.1)

text(x=log(0.05), y=relative.axis.point(0.95, "y"),
     labels=paste0((1-sig.sum)*100, "%"), pos=4, offset=0.1,
     col="red")

par(xpd=NA)
text(x=xaxis[c(1:3,6)], y=relative.axis.point(-0.06, "y"), 
     srt=0,
     adj=0.5,
     labels = c("1e-15", "1e-10", "1e-5", "1"))
text(x=xaxis[4:5], y=relative.axis.point(-0.04, "y"), 
     srt=60,
     adj=1,
     labels = c("0.01", "0.1"))
mtext(side=1, at=par("usr")[2], line=1.25,
      text=expression(italic(P)*"-value from likelihood ratio test"))

par(xpd=FALSE)

temp <- density(log(arid.random.table[3,6,]), to=0, from = log(min(arid.random.table[3,6,])))
plot(x=NULL, y=NULL, xlim=c(log(min(arid.random.table[3,6,])), 2), 
     ylim=c(0,0.1), axes=FALSE, yaxs="i",
     xlab="", ylab="")
polygon(y=c(temp$y, 0),
        x=c(temp$x, rev(temp$x)[1]), col="grey70")

temp.non <- data.frame(x=temp$x[temp$x >= log(0.05)],
                       y=temp$y[temp$x >= log(0.05)])

polygon(y=c(0, temp.non$y, 0),
        x=c(log(0.05), temp.non$x, rev(temp.non$x)[1]), col="red")

axis(side=1, at = xaxis, labels= NA)
axis(side=1, at = log(as.numeric(paste0("1e-", 15:1))), labels=NA, tck = -0.01)

abline(v=log(0.05))
box()
text(x=relative.axis.point(0.01, "x"),
     y=relative.axis.point(0.95, "y"), labels="(b)", font=2, adj=0)

sig.sum<-sum(arid.random.table[3,6,]<=0.05) / 1000
text(x=log(0.05), y=relative.axis.point(0.95, "y"),
     labels=paste0(sig.sum*100, "%"), pos=2, offset=0.1)

text(x=log(0.05), y=relative.axis.point(0.95, "y"),
     labels=paste0((1-sig.sum)*100, "%"), pos=4, offset=0.1,
     col="red")


par(xpd=NA)
text(x=xaxis[c(1:3,6)], y=relative.axis.point(-0.06, "y"), 
     srt=0,
     adj=0.5,
     labels = c("1e-15", "1e-10", "1e-5", "1"))
text(x=xaxis[4:5], y=relative.axis.point(-0.04, "y"), 
     srt=60,
     adj=1,
     labels = c("0.01", "0.1"))
par(xpd=FALSE)

dev.off()

# Obs vs Exp table ####

acac.comp.table <- do.call("rbind", lapply(1:dim(acac.slopes)[1], function(n){
  
  temp<-c(mean(acac.slopes.rand[,n]),
          sd(acac.slopes.rand[,n]))
  
  lower.tail <- ifelse(acac.slopes[n,1] > temp[1], FALSE, TRUE)
  
  data.frame(var = rownames(acac.slopes)[n],
             genus = "Acacia",
             obs = acac.slopes[n,1],
             exp.lower = qnorm(0.025, mean=temp[1], sd=temp[2]),
             exp.upper = qnorm(0.975, mean=temp[1], sd=temp[2]),
             p.value = pnorm(acac.slopes[n,1], mean=temp[1],  sd=temp[2], lower.tail=lower.tail))
}))

euc.comp.table <- do.call("rbind", lapply(1:dim(euc.slopes)[1], function(n){
  
  temp<-c(mean(euc.slopes.rand[,n]),
          sd(euc.slopes.rand[,n]))
  
  lower.tail <- ifelse(euc.slopes[n,1] > temp[1], FALSE, TRUE)
  
  data.frame(var = rownames(euc.slopes)[n],
             genus = "Eucalyptus",
             obs = euc.slopes[n,1],
             exp.lower = qnorm(0.025, mean=temp[1], sd=temp[2]),
             exp.upper = qnorm(0.975, mean=temp[1], sd=temp[2]),
             p.value = pnorm(euc.slopes[n,1], mean=temp[1],  sd=temp[2], lower.tail=lower.tail))
}))

main.table <- rbind(acac.comp.table, euc.comp.table)

acac.arid.slopes <- arid.marg.slopes[arid.marg.slopes$genus == "Acacia", ]
euc.arid.slopes <- arid.marg.slopes[arid.marg.slopes$genus == "Eucalyptus", ]
rownames(euc.arid.slopes) <- gsub("1", "", rownames(euc.arid.slopes))

acac.rand.table <-do.call("rbind", lapply(rownames(acac.arid.slopes), function(var){
  
  dry <- arid.inter.table["10%",var,]
  mesic <- arid.inter.table["90%",var,]
  
  obs <- acac.arid.slopes[var,]
  
  #Dry test
  dry.lower.tail <- ifelse(obs$mean10 > mean(dry), FALSE, TRUE)
  mesic.lower.tail <- ifelse(obs$mean90 > mean(mesic), FALSE, TRUE)
  
  # mesic test
  data.frame(var = var,
             genus = "Acacia",
             conditions = c("Dry", "Mesic"),
             obs = obs[,c("mean10","mean90")],
             exp.lower = c(qnorm(0.025, mean=mean(dry), sd=sd(dry)),
                           qnorm(0.025, mean=mean(mesic), sd=sd(mesic))),
             exp.upper = c(qnorm(0.975, mean=mean(dry), sd=sd(dry)),
                           qnorm(0.975, mean=mean(mesic), sd=sd(mesic))),
             p.value = c(pnorm(obs$mean10, mean=mean(dry),  sd=sd(dry), lower.tail=dry.lower.tail),
                         pnorm(obs$mean90, mean=mean(mesic),  sd=sd(mesic), lower.tail=mesic.lower.tail)))
}))
euc.rand.table <-do.call("rbind", lapply(rownames(euc.arid.slopes), function(var){
  
  print(var)
  dry <- arid.inter.table.euc["10%",var,]
  mesic <- arid.inter.table.euc["90%",var,]
  
  obs <- euc.arid.slopes[var,]
  
  #Dry test
  dry.lower.tail <- ifelse(obs$mean10 > mean(dry), FALSE, TRUE)
  mesic.lower.tail <- ifelse(obs$mean90 > mean(mesic), FALSE, TRUE)
  
  # mesic test
  data.frame(var = var,
             genus = "Eucalyptus",
             conditions = c("Dry", "Mesic"),
             obs = obs[,c("mean10","mean90")],
             exp.lower = c(qnorm(0.025, mean=mean(dry), sd=sd(dry)),
                           qnorm(0.025, mean=mean(mesic), sd=sd(mesic))),
             exp.upper = c(qnorm(0.975, mean=mean(dry), sd=sd(dry)),
                           qnorm(0.975, mean=mean(mesic), sd=sd(mesic))),
             p.value = c(pnorm(obs$mean10, mean=mean(dry),  sd=sd(dry), lower.tail=dry.lower.tail),
                         pnorm(obs$mean90, mean=mean(mesic),  sd=sd(mesic), lower.tail=mesic.lower.tail)))
}))

arid.table <- rbind(acac.rand.table,
                    euc.rand.table)

write.csv(main.table, "./Outputs/genus model table.csv")
write.csv(arid.table, "./Outputs/genus moisture model table.csv")

# Summary stats for paper ####

mean(with(plant.cent[!duplicated(plant.cent$pl.p),],
        tapply(pl.p, planting, length)))


sd(with(plant.cent[!duplicated(plant.cent$pl.p),],
     tapply(pl.p, planting, length)))


length(unique(plant.cent$pl.p))

summary(is.na(plant.cent$dbh))



cor(log(plant.cent$dbh) / plant.cent$age, 
    log(plant.cent$tot.ag.mass) /plant.cent$age, use="complete.obs")



# temp reviewer stats ####

envMat <- sapply(mixed.site[,c("tmax.annual", "tmin.annual", 
                               "solar.annual", "aridity.index")], scale)
envMat <- envMat[complete.cases(envMat),]
tmaxPrin <- princomp(envMat, cor=TRUE)
tmaxPrin$loadings
tmaxPrin$sdev / sum(tmaxPrin$sdev)

library(usdm)
cor(mixed.site[,c("tmax.annual", "solar.annual", "aridity.index")], use="complete.obs")
cor(mixed.site[,c("tmax.annual", "tmin.annual", 
                  "solar.annual", "aridity.index")], use="complete.obs")


vif(mixed.site[,c("tmax.annual", "solar.annual", "aridity.index")])

# range of species by climate
solarRange <- do.call("rbind", tapply(plant.cent$solar.annual, plant.cent$species, range))
solarRange <- solarRange[,2] - solarRange[,1]
hist(solarRange)
summary(solarRange==0)
summary(solarRange<1)

summary(plant.uncent$solar.annual)

aridityRange <- do.call("rbind", tapply(plant.cent$aridity.index, plant.cent$species, range))
aridityRange <- aridityRange[,2] - aridityRange[,1]
hist(aridityRange)
summary(aridityRange<1)
summary(plant.uncent$aridity.index)

summary(aridityRange==0)

# add in random slopes

plantBiggies <- droplevels(plant.cent[plant.cent$species %in% names(aridityRange)[aridityRange > 1],])

genera.mRS <- lmer(mass.per.year ~ (aridity.index + solar.annual + 
                                    SLA + WD + MH + log.age +
                                    density + (intra.sp.prop * plot.area) + 
                                    n.rarerich + FEve + FDis ) * genus +
                   (1|IBRAsub/planting/pl.p) + (aridity.index|species), data=plantBiggies,
                 control = lmerControl(optimizer="bobyqa"))

genera.mRI <- lmer(mass.per.year ~ (aridity.index + solar.annual + 
                                      SLA + WD + MH + log.age +
                                      density + (intra.sp.prop * plot.area) + 
                                      n.rarerich + FEve + FDis ) * genus +
                     (1|IBRAsub/planting/pl.p) + (1|species), data=plantBiggies,
                   control = lmerControl(optimizer="bobyqa"))

anova(genera.mRI, genera.mRS)

RSran <- ranef(genera.mRS)

aridityMean <- tapply(plant.uncent$aridity.index, plant.uncent$species, mean)

RSran <- cbind(RSran$species, aridityRange =aridityMean[aridityRange > 1])

pdf("./Plots/random slopes.pdf", height=5, width=5.5, useDingbats=FALSE)
par(mar=c(3,3,1,1), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)
plot(RSran$`(Intercept)` ~ RSran$aridityRange, 
     col=c(rgb(0.616,0.569,0), rgb(0.157,0.31,0.184))[as.factor(grepl("Euc", rownames(RSran)))],
pch=16, xlim=c(0.25,1.32))
text(y=RSran$`(Intercept)`, x=RSran$aridityRange,
     col=c(rgb(0.616,0.569,0), rgb(0.157,0.31,0.184))[as.factor(grepl("Euc", rownames(RSran)))],
     labels=rownames(RSran), font=3, cex=0.75, pos=4, offset=0.25)
mtext(side=1, line=1.5, text="Species mean moisture availability")
mtext(side=2, line=2, text="Random slope", las=0)
dev.off()