# Copyright (c) 2020, ETH Zurich

#' Internal function that saves the current simulation state
#'
#' @param val the internal state
#' @param save_val a tri-state flag: NA to not save, "all" to save the current state,
#' "last" to remove any previous save state after saving the current state
#'
#' @noRd
save_val <- function(val, save_val = NA){
  ### save current simulation state
  #save the !!global!! state of the rng
  if(is.na(save_val)){
    #save_val not set: don't save anything
    return()
  } else if(save_val == "all"){
    val$config$seed <- .GlobalEnv$.Random.seed
    val$data$distance_matrix <- NULL
    saveRDS(object = val, file = file.path(val$config$directories$output_val, paste0("val_t_", val$vars$ti, ".rds")))
  } else if(save_val == "last"){
    files <- list.files(val$config$directories$output_val, full.names = TRUE)
    val$config$seed <- .GlobalEnv$.Random.seed
    val$data$distance_matrix <- NULL
    saveRDS(object = val, file = file.path(val$config$directories$output_val, paste0("val_t_", val$vars$ti, ".rds")))
    file.remove(files)
  }
}


#' Saves the current species and landscape objects
#'
#' @param val the current simulation state
#' @noRd
save_ecogengeo <- function(val){
  # save eco at time ti
  if(is.null(val$vars$ti)){
    ti <- val$config$gen3sis$general$start_time
  } else {
    ti <- val$vars$ti
  }

  species <- val$data$all_species
  saveRDS(object = species,
          file = file.path(val$config$directories$output_species, paste0("species_t_", val$vars$ti, ".rds")))

  landscape <- val$data$landscape
  saveRDS(object = landscape,
          file = file.path(val$config$directories$output_landscapes, paste0("landscape_t_", val$vars$ti, ".rds")))

}


#' Restores the simulation form a previously saved state
#'
#' @param val a semi valid simulation state
#' @param restart_timestep the time-step to restart from
#' @noRd
restore_state <- function(val, restart_timestep){
  ### val contains a populated config, required for directory information
  ### restore previous simulation state
  if(restart_timestep == "ti"){
    #look for the most recent completed time-step
    regex <- "\\d+"
    files <- list.files(val$config$directories$output_val)
    if(length(files) == 0){
      print("no val found, starting from the initial time-step")
      return(val)
    }
    numbers <- as.integer(regmatches(files, regexpr(regex, files)))
    timestep <- min(numbers)

  } else {
    timestep <- as.integer(restart_timestep)
  }

  val <- readRDS(file.path(val$config$directories$output_val, paste0("val_t_", timestep, ".rds")))
  .GlobalEnv$.Random.seed <- val$config$seed

  if(timestep > 0){
    val$vars$save_steps <- (timestep-1):(val$config$gen3sis$general$end_time)
    val$vars$steps <- (timestep-1):(val$config$gen3sis$general$end_time)
    print(paste("restarting at time-step:", timestep))
  } else {
    val$vars$save_steps <- NULL
    val$vars$steps <- NULL
    print("simulation already completed")
  }

  return(val)
}
