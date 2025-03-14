
# if CO2, selection of start stop
# if CH4, selection of start stop + selection of diffusive pattern


clickflux <- function(dataframe, myauxfile,
                      shoulder = 0,
                      gastype, plot.lim,
                      fluxSeparation,
                      displayPlots){

  # ----------------------- Check arguments -------------------------###

  ## Check dataframe ####
  if(missing(dataframe)) stop("'dataframe' is required")
  if(!is.null(dataframe) & !is.data.frame(dataframe)){
    stop("'dataframe' must be of class data.frame")}

  ## Check myauxfile ####
  if(missing(myauxfile)) stop("'myauxfile' is required")
  if(!is.null(myauxfile) & !is.data.frame(myauxfile)){
    stop("'myauxfile' must be of class data.frame")}

  ### Is there any match between dataframe and myauxfile?
  if(!any(unique(dataframe$UniqueID) %in% unique(myauxfile$UniqueID))){
    stop("'UniqueID' in 'myauxfile' has no match for 'UniqueID' in 'dataframe'")}

  ## Check shoulder ####
  if(is.null(shoulder)) stop("'shoulder' is required") else{
    if(!is.numeric(shoulder)) stop("'shoulder' must be of class numeric") else{
      if(shoulder < 0) stop("'shoulder' cannot be a negative value")}}

  ### gastype and match in dataframe ####
  if(missing(gastype)) stop("'gastype' is required")
  if(!is.null(gastype) & !is.character(gastype)) stop("'gastype' must be a character string")
  if(!any(grepl(paste("\\<", gastype, "\\>", sep = ""),
                c("CO2dry_ppm", "COdry_ppb", "CH4dry_ppb", "N2Odry_ppb", "NH3dry_ppb", "H2O_ppm")))){
    stop("'gastype' must be of class character and one of the following: 'CO2dry_ppm', 'COdry_ppm', 'CH4dry_ppb', 'N2Odry_ppb', 'NH3dry_ppb' or 'H2O_ppm'")}
  if(!any(grepl(paste("\\<", gastype, "\\>", sep = ""), names(dataframe)))){
    stop("'dataframe' must contain a column that matches 'gastype'")}
  if(any(grepl(paste("\\<", gastype, "\\>", sep = ""), names(dataframe))) &
     !is.numeric(dataframe[,gastype][[1]])){
    stop("The column that matches 'gastype' in 'dataframe' must be of class character")}

  if(!is.numeric(plot.lim) | length(plot.lim) != 2){
    stop("'plot.lim' must be numeric and of length 2")}
  if(!is.numeric(warn.length)) {stop("'warn.length' must be of class numeric")
  } else {if(warn.length <= 0) stop("'warn.length' must be greater than 0")}
  if(!is.null(save.plots)){
    if(!is.character(save.plots)) stop("'save.plot' must be a character string")}

  ## Check fluxSeparation ####
  if(missing(fluxSeparation)) {
    fluxSeparation = FALSE
    warning("in automaticflux(), 'fluxSeparation' was not provided and automatically set to FALSE")
  }
  if(fluxSeparation != TRUE & fluxSeparation != FALSE){
    stop("'fluxSeparation' must be TRUE or FALSE")}

  ## Check displayPlots ####
  if(missing(displayPlots)) {
    displayPlots = FALSE
    warning("in automaticflux(), 'displayPlots' was not provided and automatically set to FALSE")
  }
  if(displayPlots != TRUE & displayPlots != FALSE){
    stop("'displayPlots' must be TRUE or FALSE")}



  # ----------------------- Function starts here -------------------------###

  n_incubations <- length(unique(dataframe$UniqueID))

  mydata_ow <- obs.win(inputfile = dataframe, auxfile = myauxfile, shoulder = shoulder)

  message("... Please click on start and end points of VALID OBSERVATIONS")
  mydata_manID <- click.peak2(mydata_ow, seq = seq(1,n_incubations),
                              gastype = gastype,
                              plot.lim = plot.lim)

  # Calculate fluxes
  flux_manID <- goFlux(mydata_manID, gastype)

  # Use best.flux to select the best flux estimates (LM or HM)
  # based on a list of criteria
  criteria <- c("g.factor", "kappa", "MDF", "R2", "SE.rel")

  best.flux_manID <- best.flux(flux_manID, criteria)

  if(!fluxSeparation){

    if(displayPlots){
      p <- flux.plot(
        flux.results = best.flux_manID, dataframe = mydata_manID,
        gastype = gastype, quality.check = TRUE,
        plot.legend = c("MAE", "AICc", "k.ratio", "g.factor"),
        plot.display = c("Ci", "C0", "MDF", "prec", "nb.obs", "flux.term"))
      print(p)
    }

    return(best.flux_manID)
  } else { # in  that case we proceed with the separation between diffusion and ebullition
    message("... Please now select DIFFUSION patterns only")

    # for each incubation, extract data selected at previous step
    myauxfile_corr <- NULL
    for(id in unique(myauxfile$UniqueID)){
      myauxfile_corr.tmp <- myauxfile[myauxfile$UniqueID==id,]
      myauxfile_corr.tmp$start.time <- first(mydata_manID$start.time_corr[mydata_manID$UniqueID==id])
      myauxfile_corr.tmp$obs.length <- first(mydata_manID$obs.length_corr[mydata_manID$UniqueID==id])
      myauxfile_corr.tmp$end.time <- first(mydata_manID$end.time_corr[mydata_manID$UniqueID==id])
      myauxfile_corr <- rbind(myauxfile_corr, myauxfile_corr.tmp)
    }
    mydata_ow_corr <- obs.win(inputfile = dataframe, auxfile = myauxfile_corr, shoulder = 0)

    mydata_diffusionID <- click.peak2(mydata_ow_corr, seq = seq(1,n_incubations),
                                      gastype = gastype,
                                      plot.lim = plot.lim)

    # for each incubation, extract data selected at previous step
    linear_chunk <- NULL
    for(id in unique(mydata_diffusionID$UniqueID)){
      linear_chunk.tmp <- data.frame(UniqueID = id,
                                 start.time = unique(mydata_diffusionID$start.time[mydata_diffusionID$UniqueID==id]),
                                 end = unique(mydata_diffusionID$end.time_corr[mydata_diffusionID$UniqueID==id]))
      linear_chunk.tmp$obs.length <- as.numeric(linear_chunk.tmp$end) - as.numeric(linear_chunk.tmp$start)
      linear_chunk <- rbind(linear_chunk, linear_chunk.tmp)
    }

    best.flux_manID <- lapply(seq_along(mydata_ow_corr), flux.separator.loop,
                             list_of_dataframes = mydata_ow_corr, gastype = gastype, auxfile = myauxfile, linear_chunk = linear_chunk, criteria)%>%
      map_df(., ~as.data.frame(.x))

    if(displayPlots){
      p <- separated.flux.plot(
        flux.results = best.flux_manID, dataframe = mydata_manID,
        gastype = gastype, quality.check = TRUE,
        plot.legend = c("SD"),
        plot.display = c("nb.obs"))
      print(p)
    }


    # # Calculate fluxes
    # flux_diffusion <- goFlux(mydata_diffusionID, gastype)
    #
    # best.flux_diffusion <- best.flux(flux_diffusion, criteria)
    #
    # best.flux_manID <- best.flux_diffusion
    # names(best.flux_manID)[-1] <- paste0(names(best.flux_manID)[-1],"_diffusion")
    #
    # # Estimating ebullition component
    # for (id in unique(best.flux_manID$UniqueID)){
    #   i <- which(best.flux_manID$UniqueID == id)
    #
    #   CH4_initial <-  flux_manID$C0[i]
    #   CH4_final <- flux_manID$Ct[i]
    #   incubation_time <- myauxfile_corr$obs.length[which(myauxfile_corr$UniqueID == id)]
    #   best.flux_manID$total.flux[i] <- (CH4_final-CH4_initial)/incubation_time*flux_manID$flux.term[i] # nmol/m2/s
    #   best.flux_manID$ebullition.flux[i] <- best.flux_manID$total.flux[i]  - best.flux_manID$LM.flux_diffusion[i] # total flux - diffusive term
    #   best.flux_manID$diffusion.flux[i] <- best.flux_manID$LM.flux_diffusion[i]
    # }
    #
    # # Error propagation (expressed as SD)
    # best.flux_manID$total.flux.SD <- best.flux_manID$diffusion.flux.SD <- best.flux_manID$ebullition.flux.SD <- NULL
    #
    # SD_C0 <- SD_Cf <- best.flux_manID$prec
    # deltaconcs = flux_manID$Ct-flux_manID$C0
    # SD_deltaconcs <- sqrt(SD_C0^2+SD_Cf^2)
    # SD_total.flux <- abs(best.flux_manID$total.flux) * SD_deltaconcs/deltaconcs
    # if(best.flux_manID$model == "LM"){SE_diffusion.flux = best.flux_manID$LM.SE} else {SE_diffusion.flux = best.flux_manID$HM.SE}
    # SD_diffusion.flux <- best.flux_manID$LM.SE*sqrt(best.flux_manID$nb.obs)
    # SD_ebullition.flux <- sqrt(SD_diffusion.flux^2+SD_total.flux^2)
    #
    # best.flux_manID$total.flux.SD <- SD_total.flux
    # best.flux_manID$diffusion.flux.SD <- SD_diffusion.flux
    # best.flux_manID$ebullition.flux.SD <- SD_ebullition.flux
    #
    # # adding information of data used for diffusion model
    # best.flux_manID$obs.length_diffusion <- best.flux_diffusion$obs.length[match(best.flux_diffusion$UniqueID, best.flux_manID$UniqueID)]
    #
    # if(displayPlots){
    #   p <- separated.flux.plot(
    #     flux.results = best.flux_manID, dataframe = mydata_manID,
    #     gastype = gastype, quality.check = TRUE,
    #     plot.legend = c("SD"),
    #     plot.display = c("nb.obs"))
    #   print(p)
    # }

    return(best.flux_manID)
  }
}





