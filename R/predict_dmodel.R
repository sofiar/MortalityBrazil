make_eta0 = function(W,dis_data,ndmean_index) {

  coords = dis_data$coordsForPrediction
  Amatrix = fmesher::fm_evaluator(dis_data$mesh, loc = as.matrix(coords))$proj$A
  eta = 0

  # Covariate effect
  eta = eta + W[1]

  # spatial field
  nodemean = W[ndmean_index]
  field = Amatrix %*% nodemean
  field_ras = raster::rasterFromXYZ(cbind(coords,field))
  eta = eta + field_ras
  eta
}


make_eta = function(W,dis_data,ndmean_index) {
  coords = dis_data$coordsForPrediction
  Amatrix = fmesher::fm_evaluator(dis_data$mesh, loc = as.matrix(coords))$proj$A
  eta = 0
  cov_to_predict = stack(dis_data$covariate_rasters)
  ncov = dim(cov_to_predict)[3]

  # Covariate effect
  eta = eta + W[1]
  for(i in 1:ncov){
    eta = eta+cov_to_predict[[i]] * W[i+1]
  }

  # spatial field
  nodemean = W[ndmean_index]
  field = Amatrix %*% nodemean
  field_ras = raster::rasterFromXYZ(cbind(coords,field))
  eta = eta + field_ras
  eta
  }

predict_dmodel = function(out_fit,limits,all_samples = FALSE){

aghqmodel = out_fit$aghqmodel
data = out_fit$pre_data
ncovs = as.numeric(data[[11]])

# Regression coefficients
samps = aghq::sample_marginal(aghqmodel,1e03)
beta_samps =  samps$samps[1:(ncovs+1), ]
ndmean_index = which(names(samps$samps[,1])=='nodemean')

cat("Making additive predictor rasters, AGHQ.\n")
if (ncovs==0){
  ee = base::apply(samps$samps,2,make_eta0,dis_data=data,ndmean_index=ndmean_index)
}else{ee = base::apply(samps$samps,2,make_eta,dis_data=data,ndmean_index=ndmean_index)}

#make_eta(samps$samps[1,],data,ndmean_index)
eb = brick(ee)
predmean = mean(exp(eb))
#predmean = calc(exp(eb), fun = function(x) {mean(x ,na.rm=TRUE)} )
q10 = calc(exp(eb), fun = function(x) {quantile(x,probs = c(.1,.5,.9),na.rm=TRUE)} )

coords = data$coordsForPrediction
Amatrix = fmesher::fm_evaluator(data$mesh, loc = as.matrix(coords))$proj$A
uu = samps$samps[ndmean_index,]

uest = base::apply(uu,1,mean)
field = Amatrix%*%uest
uest = raster::rasterFromXYZ(cbind(coords,field))
covcontrib = mean(eb) - uest

# Get U fild and mean incidence rate: E(lamda|Y)
cat("Getting raster.\n")

MGBorder = spTransform(limits, projection(data$polygon_shapefile))
U.rast = raster::mask(uest,MGBorder)
EL.rast = raster::mask(predmean,MGBorder)
EL1.rast = raster::mask(q10$layer.1,MGBorder)
#EL.rast = raster::mask(q10$layer.2,MGBorder)
EL2.rast = raster::mask(q10$layer.3,MGBorder)


if(all_samples){

    return(list(SpatialField=U.rast,Predmean=EL.rast,PredQ1=EL1.rast,PredQ3=EL2.rast,
              posterior_beta=beta_samps, Covcontrib = covcontrib,
              all_samples = eb))
} else {

  return(list(SpatialField=U.rast,Predmean=EL.rast,PredQ1=EL1.rast,
              PredQ3=EL2.rast, posterior_beta=beta_samps,
              Covcontrib = covcontrib))

}


}
