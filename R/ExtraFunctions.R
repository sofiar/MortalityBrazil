### Extra Functions ###
#' Plot raster object
#'
#' This function plot a raster object
#' @param rater a `SpatRaster` object to be plotted.
#' @param mortality_map a ` sf "data.frame` with the mortality data structure
#' @param limits a `SpatialPolygonsDataFrame` with boundary limits to add to the
#' final plot. If it is not given, it is taken as NULL.
#' @return plot
#' @examples
#' population_offset = load_population_offset()
#' data(femicides_data)
#' data(brazil_limits)
#' plot_raster (population_offset,femicides_data,brazil_limits)
#' @export
plot_raster = function(raster,mortality_map,limits=NULL) {

  # limits:  "SpatialPolygonsDataFrame"
  plotraster = terra::mask(raster,mortality_map)
  breaks = quantile(na.omit(values(plotraster)),
                    probs=c(0.05,0.15,0.25,0.35,0.5,0.75,0.85,0.95,1))
  terra::plot(plotraster, breaks =breaks,
              col = (RColorBrewer::brewer.pal(length(breaks),"Blues")))
  if(!is.null(limits)) sp::plot(limits,add=TRUE)

}
################################################################################


################################################################################

plot_prediction = function(out_list,fit_list,limits=NULL,breaks=NULL){

  # out_list = output object from predict_dmodel
  # fit_list = output object from fit_dmodel
  # limits = optional, boundaries to add
  # breacks = optional, breacks to plot raster objects

  mortality_sf = fit_list$pre_data$polygon_shapefile
  if(is.null(breaks)){

    max_val = range(values(out_list$PredQ3),na.rm=TRUE)[2]
    min_val = range(values(out_list$PredQ1),na.rm=TRUE)[1]
    theCol = mapmisc::colourScale(breaks = seq(round(min_val,digits = 1),
                                               round(max_val,digits = 1),
                                               length.out=9), style='fixed',
                                  col='YlOrRd', rev=FALSE)
  } else {
    theCol = mapmisc::colourScale(breaks = breaks,
             style='fixed', col='YlOrRd', rev=FALSE)}

  betas = as.data.frame(out_list$posterior_beta)

  if (dim(betas)[2]!=1) {
    betas = t(betas)
    nbetas = dim(betas)[2]-1
  } else {nbetas = 0}

  colnames(betas) = paste('beta_',seq(0,nbetas),sep='')
  betas = as.data.frame(betas)
  betas = betas %>% tidyr::pivot_longer(cols=colnames(betas),
                               names_to='beta',
                               values_to='sample')

  beta_plot = ggplot(betas, aes(x=sample))+
  geom_histogram(color="darkblue", fill="lightblue") + facet_wrap(~beta) +
  geom_vline(xintercept=0, linetype="dashed",color = "darkgrey", size=0.5) +
    theme_bw()


  plots = list(beta_plot)
  print(cowplot::plot_grid(plotlist = plots))

  readline("Press Enter to next plot")
  mapmisc::map.new(mortality_sf,mar = c(2,2,2,2))
  raster::plot(out_list$PredQ1, add=TRUE, col=theCol$col,
               breaks=theCol$breaks, legend=FALSE, maxcells=1e12)

  if(!is.null(limits)){sp::plot(limits,add=TRUE)}

  graphics::title(main = "Lower CI")
  mapmisc::legendBreaks("bottomright", theCol, bty='n', inset=0, cex=1)
  plots[[2]] = recordPlot()

  readline("Press Enter to next plot")
  mapmisc::map.new(mortality_sf,mar = c(2,2,2,2))
  raster::plot(out_list$Predmean, add=TRUE, col=theCol$col,
               breaks=theCol$breaks, legend=FALSE, maxcells=1e12)

  if(!is.null(limits)){sp::plot(limits,add=TRUE,asp=1)}

  graphics::title(main = "Mean")
  mapmisc::legendBreaks("bottomright", theCol, bty='n', inset=0, cex=1)
  plots[[3]] = recordPlot()

  readline("Press Enter to next plot")
  mapmisc::map.new(mortality_sf,mar = c(2,2,2,2))
  raster::plot(out_list$PredQ3, add=TRUE, col=theCol$col, breaks=theCol$breaks,
               legend=FALSE, maxcells=1e12)

  if(!is.null(limits))sp::plot(limits,add=TRUE)

  graphics::title(main = "Upper CI")
  mapmisc::legendBreaks("bottomright", theCol, bty='n', inset=0, cex=1)
  plots[[4]] = recordPlot()

  names(plots) = c('posterior_beta','Lower_CI','Mean_Pred','Upper_CI')
  return(invisible(plots))


}

################################################################################
#' Load Covariate Raster
#'
#' This function loads a raster file containing spatial covariates and
#'  subsets it based on the specified covariates.The covariates are based
#'  on census 2010.
#' @param names A character vector indicating the covariates to include.
#'  If `NULL`, all covariates are included.
#' @return A `SpatRaster` object containing spatial covariates
#'  \describe{
#'   \item{p_he}{Percentage of population with High School completed.}
#'   \item{p_ca}{Percentage of Catholics in the total population}
#'   \item{p_ev}{Percentage of Evangelicals in the total population}
#'   \item{p_oth}{Percentage of other religions in the total population.
#'   This includes #'   all religions except Catholicism,, Evangelical
#'   Christianity and Non Christianity.}
#'  \item{p_nr}{Percentage of individuals with no religious affiliation in the
#'  total population}
#'  \item{p_le}{Percentage of population less that High School completed.}
#' \item{p_rural}{Percentage of rurality.}
#' \item{p_lwp}{Percentage of women over 18 living with their partners.}
#' \item{p_nelec}{Percentage of population with no electricity service.}
#' \item{p_whe}{Percentage of female population with Hish School completed.}
#' \item{p_urban}{Percentage urbanity.}
#' \item{p_nonchirstian}{Percentage of the population identifying as
#' non-Christian.}
#' \item{p_ind}{Percentage of indigenous individuals in the total population.}
#' \item{p_women_resp}{Percentage of women who are responsible for the
#' household}
#' \item{p_work}{Percentage of the total population that is working}
#' \item{m_area}{Population density}
#' \item{p_black}{Percentage of the population identifying as black.}
#' \item{p_black_brown}{Percentage of the population identifying as black and
#' brown.}
#' }
#' @examples
#' covariate_raster <- load_covariate_raster('c(p_he','p_rural'))
#' covariate_raster
#' @export
load_covariate_raster = function(names = NULL) {
  output  = terra::rast(system.file("extdata", "covariates_raster.tif",
                          package = "MortalityBrazil"))
  if(!is.null(names)){
    which_covs = match(names, names(output))
    output = subset(output, which_covs)
    }

  return(output)
}



################################################################################
#' Load Population offset
#'
#' This function loads the population offset of femicides for period 2007-2010.
#' The raster is based on 2010 census data.
#' @return A `SpatRaster` object containing population offset
#' @examples
#' population_offset <- load_population_offset()
#' population_offset
#' @export
load_population_offset = function() {
  output  = terra::rast(system.file("extdata", "population_offset.tif",
                                    package = "MortalityBrazil"))
  return(output)
}
