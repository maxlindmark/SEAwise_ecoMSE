temporal.SOM = function (x, time, plot = TRUE, seed = NULL, parallel = list(parallel = FALSE, 
                                                             cores = NULL), reps = 100, return.all = FALSE, ...) 
{
  if (is.null(seed)) {
    warning("No internal seed set! Solution might be different for each run...")
  }
  else {
    if (exists(".Random.seed")) {
      old.seed = .Random.seed
    }
    else {
      set.seed(NULL)
      old.seed = .Random.seed
    }
    on.exit(assign(".Random.seed", value = old.seed, envir = .GlobalEnv))
    set.seed(seed)
    sim.seeds = stats::runif(reps, min = 1, max = 1e+09)
  }
  if (class(x) == "list") {
    if (all(lapply(x, class) %in% c("raster", "RasterBrick", 
                                    "RasterStack"))) {
      geo.pos = vector("list", length(x))
      for (i in 1:length(geo.pos)) {
        geo.pos[[i]] = raster::xyFromCell(x[[i]], 1:raster::ncell(x[[i]]))
      }
      if (do.call(all.equal.numeric, geo.pos)) {
        geo.pos = geo.pos[[1]]
      }
      else {
        stop("Rasters in list have not the same spatial resolution!")
      }
      print("1. convert raster to matrix...")
      tmp.mat = vector("list", length(x))
      for (i in 1:length(x)) {
        tmp.mat[[i]] = raster::as.matrix(x[[i]])
      }
      mat = do.call(cbind, tmp.mat)
      rm(tmp.mat)
    }
    else {
      stop("Not all elements of list(x) are raster-files!")
    }
  }
  else {
    print("1. convert raster to matrix...")
    mat = raster::as.matrix(x)
    geo.pos = raster::xyFromCell(x, 1:raster::ncell(x))
  }
  NA.rows = which(apply(mat, 1, function(x) all(is.na(x))))
  if (length(NA.rows > 0)) {
    mat = mat[-NA.rows,]
    geo.pos = geo.pos[-NA.rows, ]
  }
  print("2. perform SOM...")
  if (parallel[[1]] == TRUE) {
    library(parallel)
    library(doParallel)
    library(foreach)
    if (is.null(parallel[[2]])) {
      n.cores = parallel::detectCores() - 1
    }
    else {
      n.cores = parallel[[2]]
    }
    cl = parallel::makeCluster(n.cores)
    doParallel::registerDoParallel(cl)
    cat("Running in parallel mode...", "\n")
    tmp.store = foreach::foreach(i = 1:reps, .packages = "kohonen") %dopar% 
    {
      if (!is.null(seed)) {
        set.seed(sim.seeds[i])
      }
      # temporal mode SOM
      SOM = kohonen::som(X = mat, ...)
    }
    parallel::stopCluster(cl)
  }
  else {
    cat("Running in single core mode...", "\n")
    tmp.store = rep(list(NA), reps)
    for (i in 1:reps) {
      if (!is.null(seed)) {
        set.seed(sim.seeds[i])
      }
      tmp.store[[i]] = kohonen::som(X = mat, ...)
    }
  }
  if (reps == 1) {
    best.SOM = tmp.store[[1]]
    SOM.quality = NULL
    warning("Robustness of SOM was not evaluated!")
  }
  else {
    cat("Checking quality of the solutions...", "\n")
    QE = lapply(tmp.store, marmalaid::SOM.quant.error)
    TE = lapply(tmp.store, marmalaid::SOM.topo.error)
    all.TEs = sapply(TE, function(x) x$TE)
    all.QEs = do.call(c, QE)
    indx.TE = which(all.TEs == min(all.TEs))
    indx.QE = which(all.QEs == min(all.QEs[indx.TE]))
    best.SOM = tmp.store[[indx.QE]]
    SOM.quality = data.frame(QE = all.QEs, TE = all.TEs)
  }
  print("3. convert back to raster-format...")
  #SOM.space = best.SOM$unit.classif
  SOM.space.df = data.frame(geo.pos, best.SOM$unit.classif)
  SOM.raster = raster::rasterFromXYZ(SOM.space.df)
  # if (class(x) == "list") {
  #   split.df = split(SOM.space.df, f = rep(1:length(x), 
  #                                          each = nrow(geo.pos)))
  #   SOM.raster = vector("list", length = length(x))
  #   for (i in 1:length(x)) {
  #     SOM.raster[[i]] = raster::rasterFromXYZ(split.df[[i]])
  #   }
  # }
  # else {
  #   SOM.raster = raster::rasterFromXYZ(SOM.space.df)
  # }
  
  if (class(x) == "list") {
    split.df = split(data.frame(t(best.SOM$codes[[1]])), f = rep(1:length(x),
                                           each = length(time)))
    SOM.ts = vector("list", length = length(x))
    for (i in 1:length(x)) {
      SOM.ts[[i]] = data.frame(time = time, ts = split.df[[i]])
    }
    names(SOM.ts) = names(x)
  }
  else {
    SOM.ts = data.frame(time = time, ts = t(best.SOM$codes[[1]]))
  }
  
  
  if (plot == TRUE) {
    opar = graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(opar))
    graphics::par(mfrow = c(2, 2),mar = c(2,2,3,2))
    graphics::plot(best.SOM, type = "changes")
    graphics::plot(best.SOM, type = "codes")
    graphics::plot(best.SOM, type = "mapping")
    image(SOM.raster,main = "map")
    maps::map(add = T)
  }
  if (return.all == FALSE) {
    return.all.SOMs = NULL
  }
  else {
    return.all.SOMs = tmp.store
  }
  output = list(SOM.out = best.SOM, SOM.quality = SOM.quality, 
                SOM.raster = SOM.raster, 
                Clustered.ts = SOM.ts, 
                return.all.SOMs = return.all.SOMs)
  return(output)
}
