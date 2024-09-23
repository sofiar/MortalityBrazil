################################### Example ####################################

# Analysis of intentional aggression (X85-Y09) and from aggression of undermined
# intent (Y10-Y34) of women  between 10 and 49 years old
# Period: 2007-2010

devtools::load_all("./")
devtools::document()

############################## 1. Load data ####################################

# Brazil limits
data(brazil_limits)

# Spatial Covariates
name_covs =  c('p_he','p_urban','p_wome_resp','p_ind', 'p_black_brown','p_ev')
covariate_raster = load_covariate_raster(name_covs)

# Population offset
population_offset = load_population_offset()

# Mortality data
data(mortality_data)

############################## 2. Plot data ####################################

# Population offset
plot_raster (population_offset,mortality_data,brazil_limits)
# Spatial covariates
plot_raster(covariate_raster,mortality_data)
# Mortality counts
ggplot(mortality_data) + geom_sf(aes(fill = nCounts))+ theme_minimal()

####################### 3. Fit model and prediction ############################

#set priors
priors = list(priormean_intercept = 0,
              priorsd_intercept = 10,
              priormean_slope = 0.0,
              priorsd_slope = 10,
              prior_rho_min = 1e5,
              prior_rho_prob =.01,
              prior_sigma_max = 5,
              prior_sigma_prob =.01,
              prior_iideffect_sd_max=.1,
              prior_iideffect_sd_prob =.01)

# Fit model
result = fit_dmodel(mortality_data, population_offset, id_var='code_muni',
                    response_var='nCounts', priors,
                    covariate_rasters = covariate_raster,
                    startingvalue=c(0,0,20),
                    k=7)

summary(result$aghqmodel)

# Make predictions
model_pred = predict_dmodel(result,brazil_limits)
plot_prediction(model_pred, result, limits = brazil_limits,
                breaks = seq(0,3,by=0.35))



