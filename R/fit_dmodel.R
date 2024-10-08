#####################################################################################
### This function fit a dissagregation model Based on the ELGM approach of        ###
### Stringer and others 2021 (DOI: https://doi.org/10.1080/10618600.2022.2099403) ###
### This code uses source code of the disaggregation package:                     ###
### https://cran.r-project.org/web/packages/disaggregation/index.html             ###
#####################################################################################
#' Fit a dissagregation model based on the Extended Latent Gaussian Process
#' approach of Stringer and others 2021
#'
#' This function fit a dissagregation model
#' @param mortality_sf a ` "sf"  "data.frame"` object with mortality counts.
#' @param offset_sr a `SpatRaster` with the population offset
#' @param id_var a `character` with the name of polygon variable
#' @param response_var a `character` with name of response variable
#' @param priors a `list` with priors
#' @param covariate_raster a `SpatRaster` with covariates
#' @param k `integer`, the number of quadrature points to use.  Default k=5.
#' k= 1 corresponds to a Laplace approximation.
#' @param startingvalue Value to start the optimization
#' @references Stringer, A., Brown, P., & Stafford, J. (2022). Fast, Scalable
#' Approximations to Posterior Distributions in Extended Latent Gaussian Models.
#' Journal of Computational and Graphical Statistics, 32(1), 84–98.
#' https://doi.org/10.1080/10618600.2022.2099403
#' @examples
#'result = fit_dmodel(mortality_data, population_offset, id_var='code_muni',
#'                    response_var='nCounts', priors,
#'                    covariate_rasters = covariate_raster,
#'                    startingvalue=c(0,0,20),
#'                    k=7)
#'summary(result$aghqmodel)
#' @export
fit_dmodel = function(mortality_sf,offset_sr,id_var,response_var,priors,
                    covariate_rasters = NULL, k=5,startingvalue=rep(0,3)){

  # Fit model with no covariates
  if(is.null(covariate_rasters)){
    data = disaggregation::prepare_data(
    polygon_shapefile = mortality_sf,
    covariate_rasters = offset_sr,
    aggregation_raster = offset_sr,
    id_var = id_var,
    response_var = response_var,
    na.action = TRUE)


    spde = (INLA::inla.spde2.matern(data$mesh, alpha = 1 + 1)$param.inla)[c("M0", "M1", "M2")]
    Apix = fmesher::fm_evaluator(data$mesh, loc = data$coordsForFit)$proj$A
    n_s = nrow(spde$M0)

    field <- iid <- TRUE

    cov_matrix = as.matrix(data$covariate_data[, -c(1:2)])
    if (ncol(cov_matrix) == 1) {
      cov_matrix = as.matrix(base::apply(cov_matrix,1, as.numeric))
    } else {
      cov_matrix = t(base::apply(cov_matrix, 1, as.numeric))
      }

    parameters = list(intercept = -10,
                       slope =  rep(0, ncol(cov_matrix)),
                       log_tau_gaussian = 8,
                       iideffect = rep(0, nrow(data$polygon_data)),
                       iideffect_log_tau = 1,
                       log_sigma = 0,
                       log_rho = 4,
                       nodemean = rep(0, n_s))

    input_data = list(x = cov_matrix,
                       aggregation_values = data$aggregation_pixels,
                       Apixel = Apix,
                       spde = spde,
                       startendindex = data$startendindex,
                       polygon_response_data = data$polygon_data$response,
                       response_sample_size = data$polygon_data$N,
                       family = 2, # Poisson distribution
                       link = 1, # Log link
                       nu = 1,
                       field = as.integer(field),
                       iid = as.integer(iid))

    input_data = c(input_data, priors)
    tmb_map = list(log_tau_gaussian = as.factor(NA), slope=as.factor(NA))
    ## END: taken directly from disaggregate::make_model_object()

    random_effects = c('nodemean','iideffect','intercept')

    obj = TMB::MakeADFun(
      data = input_data,
      parameters = parameters,
      map = tmb_map,
      random = random_effects,
      hessian = TRUE,
      silent = TRUE, # Always do silent = TRUE
      DLL = "disaggregation")

    ## Fit the model using AGQH ----
    # theta = log(tau), log(sigma), log(rho)
    cat("Fitting using AGHQ.\n")
    aghqmodel = aghq::marginal_laplace_tmb(obj,k=k,startingvalue = startingvalue,
                                           control = aghq::default_control_tmb())

    data[[11]] = 0
    names(data)[11]='Ncov'

  } else {  # Fit model with covariates
    data = disaggregation::prepare_data(
    polygon_shapefile = mortality_sf,
    covariate_rasters = covariate_rasters,
    aggregation_raster = offset_sr,
    id_var = id_var,
    response_var = response_var,
    na.action = TRUE)

    # This was added to avoid NA in covariates. However, when 2 or more covariates
    # are included this should not be necessary.
    # Is this a bug in disaggregation::prepare_data ?

    n_covs = dim(data$covariate_data)[2]-2
    for(i in 1:n_covs){
      which_nones = which(is.na(data$covariate_data[,i]))
      if (sum(length(which_nones))>0){
        data$covariate_data[which_nones,i] =
          median(data$covariate_data[,i],na.rm = TRUE)
      }
    }

    # End of the added

    spde = (INLA::inla.spde2.matern(data$mesh, alpha = 1 + 1)$param.inla)[c("M0", "M1", "M2")]
    Apix = fmesher::fm_evaluator(data$mesh, loc = data$coordsForFit)$proj$A
    n_s = nrow(spde$M0)

    field <- iid <- TRUE

    # Set covariates. This depend on the number of covarites
    ndim = dim(data$covariate_data)[2]
    cov_matrix = as.matrix(data$covariate_data[, -c((ndim-1):ndim)])

    if (ncol(cov_matrix) == 1) {cov_matrix <- as.matrix(base::apply(cov_matrix, 1, as.numeric))
    }else {cov_matrix = t(base::apply(cov_matrix, 1, as.numeric))}

    parameters = list(intercept = -5,
                       slope =  rep(0, ncol(cov_matrix)),
                       log_tau_gaussian = 8,
                       iideffect = rep(0, nrow(data$polygon_data)),
                       iideffect_log_tau = 1,
                       log_sigma = 0,
                       log_rho = 4,
                       nodemean = rep(0, n_s))


    input_data = list(x = cov_matrix,
                       aggregation_values = data$aggregation_pixels,
                       Apixel = Apix,
                       spde = spde,
                       startendindex = data$startendindex,
                       polygon_response_data = data$polygon_data$response,
                       response_sample_size = data$polygon_data$N,
                       family = 2, #Poisson Model
                       link = 1, # log link
                       nu = 1,
                       field = as.integer(field),
                       iid = as.integer(iid))
    input_data = c(input_data, priors)
    tmb_map = list(log_tau_gaussian = as.factor(NA))

    random_effects = c('nodemean','iideffect','intercept','slope')

    obj = TMB::MakeADFun(
      data = input_data,
      parameters = parameters,
      map = tmb_map,
      random = random_effects,
      hessian = TRUE,
      silent = TRUE, # Always do silent = TRUE
      DLL = "disaggregation")

    cat("Fitting using AGHQ.\n")
    aghqmodel = aghq::marginal_laplace_tmb(obj,k=k,startingvalue = startingvalue, control = aghq::default_control_tmb())

    data[[11]] = ndim-2
    names(data)[11] = 'Ncov'
    }

    return(list(aghqmodel=aghqmodel,pre_data=data))
}










