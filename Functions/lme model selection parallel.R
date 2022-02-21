parallel.lme.model.selection<-function(response,
                                       main.fixed,
                                       covariate.fixed=NULL,
                                       interactions=NULL,
                                       random.effects,
                                       dataframe,
                                       AIC.threshold=2,
                                       cores=detectCores() - 1){
  
  ## AUTHOR: Timothy Staples
  
  ## DATE: 01/07/2016
  
  ## FUNCTION PURPOSE: Run two-step lme model simplification using AIC. 
                       # Features: # Fixed terms that are retained regardless of
                       #             importance to model, for hypothesis testing
  #           # Customise AIC threshold to retain variables
  #           # Standardise numeric variables to allow
  #             comparison of partial regression slopes.
  
  ## ARGUMENTS: response: String representing data frame column heading for 
  ##                      response variable of models
  ##          main.fixed: String vector of data frame column headings to act
  ##                      as fixed terms in model. These terms are kept in
  ##                      all models even if deemed unimportant in explaining
  ##                      response variation
  ##     covariate.fixed: String vector of data frame column headings to act
  ##                      as fixed terms in model. These terms can be removed in
  ##                      second step of model simplification.
  ##        interactions: String vector of interaction terms to fit to model. 
  ##                      Terms must be listed in correct format ("var1:var2"),
  ##                      variables must match data frame column headings and be
  ##                      present as main terms in the model.
  ##      random.effects: String of random effects to feed into lme. Must be
  ##                      specified in correct format to be read by as.formula().
  ##                      See ?lme for more details on format.
  ##           dataframe: Data frame containing terms listed in response,
  ##                      main.fixed and covariate.fixed. Other columns are
  ##                      allowed, as they are subsetted in the model.
  ##       AIC.threshold: Custom AIC threshold to select variables. Compared to
  ##                      the increase in AIC when variable is removed from model.
  
  
  # FUNCTION CHECKS ####
  
  # make sure packages are loaded
  package.temp<-c(library(nlme, logical.return=TRUE),
                  library(MuMIn, logical.return=TRUE))
  if(sum(!package.temp)>0){
    stop(paste("\"", c("nlme","MuMIn")[!package.temp], "\"",
               "is required for this function and is not installed."))
  }
  
  # make sure all interactions have main terms
  if(!is.null(interactions)){
    main.terms<-c(main.fixed, covariate.fixed)
    interaction.terms<-unique(c(substr(interactions,
                                       1,
                                       regexpr(':',interactions)-1),
                                substr(interactions,
                                       regexpr(':',interactions)+1,
                                       100)))
    missing.main<-interaction.terms[!interaction.terms %in% main.terms]
    
    if(length(missing.main)>0){
      stop(paste("\"", missing.main, "\"", 
                 " is listed as part of an interaction, but is not present as a main effect.",
                 sep=""))
    }
  } # end interaction check
  
  
# DATA PREP ####
 
  # subset rows with NAs
  model.data<-dataframe[apply(is.na(dataframe), 1, function(x){sum(x)==0}),]
  
  if(dim(model.data)[1] != dim(dataframe)[1]){
    warning(paste(dim(dataframe)[1] - dim(model.data)[1],
                  " rows were removed due to NAs."))
  }
  
  # Center predictors to mean=0 and sd=1
  # remove non-numeric variables
  center.vect<-c(main.fixed, covariate.fixed)
  center.class<-sapply(model.data[colnames(model.data) %in% center.vect], class)
  center.numeric<-names(center.class)[center.class=="numeric"]
  center.data<-model.data[,colnames(model.data) %in% center.numeric]
  
  # Now center data if there's at least 1 numeric variable
  if(sum(center.class=="numeric")>1){
    center.data.output<-sapply(center.data, function(x){(x-mean(x))/sd(x)})
    model.data[,colnames(model.data) %in% center.numeric]=center.data.output
  }
  if(sum(center.class=="numeric")==1){
    center.data.output=(center.data-mean(center.data))/sd(center.data)
    model.data[,colnames(model.data) %in% center.numeric]=center.data.output
  }
  
  # Set up model formula
  f.fixed<- paste(response, "~", paste(c(main.fixed, 
                                         covariate.fixed, 
                                         interactions), collapse=" + "))
  
# INITIAL MODEL ####
  
      # First check for error
  error.trap<-tryCatch(
    do.call("lme", 
            list(as.formula(f.fixed), 
                 data=as.name("model.data"),
                 random=as.formula(random.effects))),
    error=function(e) e)
  
      # Run model if no error, else return error message
  if(!inherits(error.trap, "error")){
    m1<-do.call("lme", list(as.formula(f.fixed), 
                            data=as.name("model.data"),
                            random=as.formula(random.effects)))
  } else {
    stop(paste("Full model has failed. Error message from lme:", error.trap))
  }
  
  
  
# SUBSET AND TEST INTERACTIONS ####
  
  if(length(interactions)>0){
        
      # create list of models removing one interaction at a time
    no_cores <- cores
    cl <- makeCluster(no_cores)
    setDefaultCluster(cl)
    # export data and load libraries in cluster
    clusterExport(varlist=c("model.data", "m1"),
                  envir=environment())
    clusterCall(cl, "library", "nlme", "MuMIn", character.only = TRUE)
    options(na.action = "na.fail")
    int.list<-parLapply(cl, interactions, function(x){update(m1, as.formula(paste0(".~. -",x)))})
    stopCluster(cl=cl)
  
      # compare models using AIC
  m.sel.temp<-model.sel(int.list, rank="AIC")[order(rownames(model.sel(int.list)))]
  
      # record which interactions to retain
  interactions.to.retain<-interactions[as.numeric(rownames(m.sel.temp)[m.sel.temp$AIC-AIC(m1)>AIC.threshold])]
  
  } else {interactions.to.retain=NULL} # End interaction IF
  
# REDUCED MODEL ####
  
      # update fixed effects formula
  f.fixed <- paste(response, " ~ ", paste(c(main.fixed, covariate.fixed,
                                          interactions.to.retain), collapse=" + "))
  
  # First check for error
  error.trap<-tryCatch(
    do.call("lme", 
            list(as.formula(f.fixed), 
                 data=as.name("model.data"),
                 random=as.formula(random.effects))),
    error=function(e) e)
  
  # Run model if no error, else return error message
  if(!inherits(error.trap, "error")){
    m2<-do.call("lme", list(as.formula(f.fixed), 
                            data=as.name("model.data"),
                            random=as.formula(random.effects)))
  } else {
    stop(paste("Interaction subset model has failed. Error message from lme:", 
               error.trap))
  }
  

# SUBSET AND TEST NON-INTERACTING MAIN EFFECTS ####
  
      # Subset out any remaining interactions
  int.temp<-rownames(coefTable(m2))[grepl(":", rownames(coefTable(m2)))]
  
      # get main terms in these interactions (subsetting string before and after
      # ":"). If there were no important interactions, contingency creates an empty vector
  if(length(int.temp)>0){
    interacting.main.effects<-unique(c(substr(int.temp,
                                              1,
                                              regexpr(pattern =':',int.temp)-1),
                                       substr(int.temp,
                                              regexpr(pattern =':',int.temp)+1,
                                              100)))
  } else {interacting.main.effects<-vector()}
  
      # get non-interacting main terms
  main.effects<-rownames(coefTable(m2))[
    !grepl(paste(c(interacting.main.effects, 
                   main.fixed, 
                   "Intercept"), collapse="|"), 
           rownames(coefTable(m2)))]
  
      # As long as there's at least one non-interacting main term,
      # run models subsetting one main effect at a time
  if(length(main.effects)!=0){
    
      # Get list of models
    no_cores <- cores
    cl <- makeCluster(no_cores)
    setDefaultCluster(cl)
    # export data and load libraries in cluster
    clusterExport(varlist=c("model.data", "m2"),
                  envir=environment())
    clusterCall(cl, "library", "nlme", "MuMIn", character.only = TRUE)
    options(na.action = "na.fail")
    main.list<-parLapply(cl, main.effects, function(x){update(m2, as.formula(paste0(".~. -",x)))})
    stopCluster(cl=cl)
    
      # Compare with AIC
    m.sel.temp<-model.sel(main.list)[order(rownames(model.sel(main.list)))]
    
      # Record main effects to retain
    non.interacting.main.effects<-main.effects[as.numeric(rownames(m.sel.temp)[m.sel.temp$AICc-AICc(m2)>AIC.threshold])]
    
    # where there are no non-interacting main terms, 
    # keep them all out as they'll be put in in as "interacting.main.effects"
    } else { non.interacting.main.effects=NULL}
# FINAL MODEL ####
  
  formula.terms<-c(main.fixed, 
                   interacting.main.effects, 
                   non.interacting.main.effects, 
                   interactions.to.retain)
  
  f.fixed <- paste(response, " ~ ", 
                   paste(formula.terms, collapse=" + "))
  
  # First check for error
  error.trap<-tryCatch(
    do.call("lme", 
            list(as.formula(f.fixed), 
                 data=as.name("model.data"),
                 random=as.formula(random.effects))),
    error=function(e) e)
  
  # Run model if no error, else return error message
  if(!inherits(error.trap, "error")){
    m3<-do.call("lme", list(as.formula(f.fixed), 
                            data=as.name("model.data"),
                            random=as.formula(random.effects)))
  } else {
    stop(paste("Final model has failed. Error message from lme:", 
               error.trap))
  }

  return(m3)
  
  } # end Function
