

# ================== Script containing setup priors and functions for INLA spatiotemporal models ====================



# ----------------- 1. Project hyperpriors -------------------

# iid model 
hyper.iid = list(theta = list(prior="pc.prec", param=c(0.5, 0.01)))

# ar1 model
hyper.ar1 = list(theta1 = list(prior='pc.prec', param=c(0.5, 0.01)),
                  rho = list(prior='pc.cor0', param = c(0.5, 0.75)))

# bym model
hyper.bym = list(theta1 = list(prior="pc.prec", param=c(1, 0.01)),
                 theta2 = list(prior="pc.prec", param=c(1, 0.01)))

# bym2 model
hyper.bym2 = list(theta1 = list(prior="pc.prec", param=c(1, 0.01)),
                  theta2 = list(prior="pc", param=c(0.5, 0.5)))

# hyperpriors for model grouping (iid / ar1) if used
# group.control.iid = list(model='iid', hyper = list(prec = list(prior='pc.prec',param=c(1, 0.01))))
# group.control.ar1 = list(model='ar1', hyper = list(theta1 = list(prior='pc.prec', param=c(1, 0.01)), rho = list(prior='pc.cor0', param = c(0.5, 0.75))))

# rw1/rw2 model: three levels of constraint on precision parameter 
# (puts more or less prior probability density on more or less wiggly)
hyper1.rw = list(prec = list(prior='pc.prec', param=c(0.1, 0.01))) # strictest smoothing; sd constrained to be low
hyper2.rw = list(prec = list(prior='pc.prec', param=c(0.3, 0.01))) # medium
hyper3.rw = list(prec = list(prior='pc.prec', param=c(1, 0.01))) # weaker (suggested INLA default) 
hyper4.rw = list(prec = list(prior='pc.prec', param=c(2, 0.01))) # weakest; sd can be quite wide 




# ------------------ inla functions --------------------

### fitINLAModel: fit and return INLA model with specified formula and family

#' @param formx inla formula object; i.e. created using formula(y ~ x + f())
#' @param family likelihood
#' @param stack inla stack object for model
#' @param spde spde object for fitting
#' @param config boolean (default FALSE): set config in compute to TRUE for inla.posterior.sample()
#' @param verbose verbose reporting on or off? default FALSE
#' @param return.marginals specify whether model should save and return marginals

fitINLAModel = function(formx, family, stack, spde, verbose=FALSE, config=FALSE, return.marginals=FALSE, inla.mode="experimental"){
  return(
    inla(formx,
         verbose = verbose,
         data = inla.stack.data(stack, spde=spde),
         family=family,
         control.fixed = control.fixed1, 
         control.predictor=list(A=inla.stack.A(stack), 
                                compute=TRUE, 
                                link=1),
         control.compute=list(cpo=TRUE, waic=TRUE, dic=TRUE, 
                              config=config, 
                              return.marginals=return.marginals),
         control.inla = list(strategy='adaptive', # adaptive gaussian
                             cmin=0), # fixing Q factorisation issue https://groups.google.com/g/r-inla-discussion-group/c/hDboQsJ1Mls)
         inla.mode = inla.mode)
  )
}







# ------------------ 3. Wrapper functions to extract useful stuff from fitted INLA models -----------------


### fitMetricsINLA: extract and report on metrics of fit

#' @param mod fitted INLA model
#' @param data data used to fit the model; n.b. assumes that response variable column is called "y"
#' @param modname name to give model in resulting dataframe

fitMetricsINLA = function(mod, data, modname="mod"){
  
  dx = data[ , c("y"), drop=FALSE] %>%
    dplyr::mutate(fitted = mod$summary.fitted.values$mean,
                  abs_err = abs(y - fitted),
                  cpo = mod$cpo$cpo,
                  cpo_fail = mod$cpo$failure)
  fit = data.frame(modname = modname,
                   dic = mod$dic$dic, 
                   waic = mod$waic$waic,
                   waic_neffp = mod$waic$p.eff,
                   mae = mean(dx$abs_err, na.rm=TRUE),
                   logscore = -mean(log(dx$cpo), na.rm=TRUE),
                   cpo_fail = sum(dx$cpo_fail == 1 & !is.na(dx$cpo_fail)))
  return(fit)
}


### extractRandomINLA: extract random effect and rename columns

# if effect is grouped/replicated by a factor, automatically assign each subgroup to its grouping factor (labelled 1:n) 
# if BYM model, further partition into u and v components

#' @param summary_random points to model$summary.random$effect_of_interest
#' @param effect_name name to assign to fitted effect in dataframe (can be anything)
#' @param model_is_bym boolean; to specify if model is joint Besag-York-Mollie
#' @param transform specify whether to exponentiate coefficients (i.e. back transform to relative risk)
extractRandomINLA = function(summary_random, effect_name, model_is_bym=FALSE, transform=FALSE){
  
  # extract model effect
  rf = summary_random %>%
    dplyr::rename("value"=1, "lower"=4, "median"=5, "upper"=6)

  # label by grouping factor (if not replicated, group is 1 for all observations)
  rf$group = rep(1:as.vector(table(rf$value)[1]), each=n_distinct(rf$value))
  
  # partition BYM into u and v components
  if(model_is_bym){
    rf$component = rep(c("uv_joint", "u_besag"), each=n_distinct(rf$value)/2)
    rf$value = rep(1:(n_distinct(rf$value)/2), n_distinct(rf$group)*2)
  }
  
  # back transform if specified
  if(transform == TRUE){
    rf[ , 2:7 ] = exp(rf[ , 2:7])
  }
  
  # name and return
  rf$effect = effect_name
  return(rf)
}


## extractFixedINLA: extracts fixed effects from specified model

#' @param model fitted INLA object
#' @param model_name modelname to provide in dataframe
#' @param transform exponentiate coefficients to RR/OR scale YN?

extractFixedINLA = function(model, model_name="mod", transform=FALSE){
  ff = model$summary.fixed
  ff$param = row.names(ff)
  ff$param[ ff$param == "(Intercept)" ] = "Intercept"
  names(ff)[3:5] = c("lower", "median", "upper")
  if(transform == TRUE){
    ff[ 1:5 ] = exp(ff[ 1:5 ])
  }
  ff$model = model_name
  ff
}


### rasteriseSPDE: generates a raster of the fitted spatial field masked by actual study area

#' @param fitted_mesh_x = model$summary_random$spde$mean
#' @param mesh_x = mesh object used for model
#' @param shp_x = shapefile of study boundaries used for defining mesh (to mask to correct extent)

rasteriseSPDE = function(fitted_mesh_x, mesh_x, shp_x){
  
  proj = inla.mesh.projector(mesh_x, dims=c(500, 500))
  full.proj = expand.grid(x = proj$x, y=proj$y)
  full.proj$field = c(inla.mesh.project(proj, fitted_mesh_x))
  pras = rasterFromXYZ(full.proj)
  mask = fasterize::fasterize(shp_x, pras)
  pras = mask(pras, mask)
  pras = crop(pras, shp_x)
  return(pras)
}

# extract info criteria from model
getIC = function(m, mod_name){
  data.frame(
    model=mod_name,
    waic = m$waic$waic,
    dic = m$dic$dic
  )
}
