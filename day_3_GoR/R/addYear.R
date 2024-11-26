addYear <- function(dat, nyears, method = 'last', avg_years = 3){

  if(is.data.frame(dat)){

    if(method == 'last'){
      tmp <- dat[nyears, ]
      tmp$years <- tmp$years+1
      dat.out <- rbind(dat, tmp)
    }

    if(method == 'average'){

      tmp <- (colMeans(dat[(nyears-avg_years):nyears,]))
      tmp <- as.data.frame(t(tmp))
      tmp$years <- max(dat$years)+1
      dat.out <- rbind(dat, tmp)

    }

  }

  if(is.array(dat)){

    if(method == 'last'){
      nage <- dim(dat)[1]

      tmp <- array(dat[,nyears,1], dim = c(nage, 1, 1))
      dat.out <- abind::abind(dat, tmp, along = 2)

    }

    if(method == 'average'){

      if(dim(dat)[3] > 1){
        stop('only 1 season supported currently')
      }

      nage <- dim(dat)[1]
      tmp <- array(((rowMeans(dat[,(nyears-avg_years):nyears,1]))), dim = c(nage, 1, 1))
      dat.out <- abind::abind(dat, tmp, along = 2)

    }



  }



  return(dat.out)
}
