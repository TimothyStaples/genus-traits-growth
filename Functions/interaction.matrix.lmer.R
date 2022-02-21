interaction.matrix.fun.lmer<-function(model, 
                                 x.var, 
                                 y.var, 
                                 grid.size, 
                                 quantile){
  ## AUTHOR: Timothy Staples
  
  ## DATE: 01/07/2016
  
  ## FUNCTION PURPOSE: Create surface plot matrix for two-way interaction
  ##                   between two continuous variables.
  
  ## ARGUMENTS: model: model of a type with a valid 'predict' function.
  ##             x.var: String, matches model term. Interacting variable to
  ##                    be plotted along x-axis of surface.
  ##             y.var: String, matches model term. Interacting variable to
  ##                    be plotted along x-axis of surface.
  ##         grid.size: number of grid cells along each axis.
  ##          quantile: Quantile to use as axis limits 
  ##                    e.g., 0 = max and min value,
  ##                       0.01 = 1% and 99% quantile).
  
  ## OUTPUT: A list containing 4 elements:
  ##                x: Vector of grid cells along the x axis.
  ##                y: Vector of grid cells along the y axis.
  ##                z: Matrix of grid cell value from model prediction.
  ##               se: Matrix of standard error grid cell value from model
  ##                   prediction.
  
# FUNCTION CHECKS ####

  # make sure interaction term exists
  # get mean values for covariates
  
# DATA PREP ####
  
  # extract model data
  model.data<-model@frame
  
  # extract x and y data columns
  x.var.data<-model.data[,colnames(model.data) %in% x.var]
  y.var.data<-model.data[,colnames(model.data) %in% y.var]
  
  # generate x and y values along each axis
  x.points<-seq(quantile(x.var.data, probs=quantile, na.rm=TRUE),
                quantile(x.var.data, probs=1-quantile, na.rm=TRUE), 
                length.out=grid.size)
  
  y.points<-seq(quantile(y.var.data, probs=quantile, na.rm=TRUE),
               quantile(y.var.data, probs=1-quantile, na.rm=TRUE), 
               length.out=grid.size)
  
  # extract out main terms from model to feed into predict function
  model.main.effects<-rownames(coefTable(model))[-1][!grepl(":", rownames(coefTable(model))[-1])]
  
  # create pair-wise combination of x and y points, one for each grid cell
  interaction.grid<-as.matrix(expand.grid(x.points, y.points))
  
  interaction.mat<-matrix(0, nrow=length(interaction.grid[,1]),
                          ncol=length(model.main.effects), 
                          dimnames=list(NULL, model.main.effects))
  interaction.mat[,colnames(interaction.mat)==x.var]=interaction.grid[,1]
  interaction.mat[,colnames(interaction.mat)==y.var]=interaction.grid[,2]
  
  predictions<-predict(object=model, newdata=as.data.frame(interaction.mat), re.form=~0, se.fit=TRUE)
  inter.list<-list(x=x.points, y=y.points,
                   z=matrix(predictions$fit, 
                            nrow=length(x.points), ncol=length(y.points)),
                   se=matrix(predictions$se.fit, 
                             nrow=length(x.points), ncol=length(y.points)))
  return(inter.list)
  
}
