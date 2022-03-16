require(Distance)
require(MASS)
require(tidyverse)

integrate2 <- function(f,lower,upper,...,subs=100){
  h <- (upper-lower)/subs
  sq <- lower+(1:(subs-1))*h
  (upper-lower) * (f(lower,...)+2*sum(f(sq,...))+f(upper,...)) / (2*subs)
}

dfpdf <- function(x, prm, ispoint, key, series, order, exp, w){
  res <- switch(key,
                "unif" = rep(1/w, length(x)),
                "hn" = mrds:::keyfct.hn(x, prm[1]),
                "hr" = mrds:::keyfct.hz(x, prm[2], prm[1])
                )
  aprm <- switch(key,
                 "unif" = prm,
                 "hn" = prm[-1],
                 "hr" = prm[-(1:2)]
                 )
  if(length(aprm)>0){
    adj <- switch(series,
                  "herm" = mrds:::adjfct.herm(x, scaling=w, order, aprm, exp),
                  "cos" = mrds:::adjfct.cos(x, scaling=w, order, aprm, exp),
                  "poly" = mrds:::adjfct.poly(x, scaling=w, order, aprm, exp)
                  )
    res <- res*(1+adj)
  }
  if(ispoint) res <- x*res
  res
} 

esw <- function(prm, ispoint, key, series, order, exp, w){
  res <- integrate2(dfpdf, 0, w, prm, ispoint, key, series, order, exp, w)
  if(ispoint) res <- sqrt(2*res)
  res
}

fitdf <- function(formula, data, reps=999, ...){
  args <- list(...)
  depname <- all.vars(formula)[1]
  classes <- dplyr::summarise_all(data, class)
  if(classes[depname]=="numeric") 
    data <- dplyr::rename(data, distance=all_of(depname)) else{
      cats <- strsplit(as.character(dplyr::pull(data, depname)), "-")
      data$distbegin <- unlist(lapply(cats, function(x) as.numeric(x[1])))
      data$distend <- unlist(lapply(cats, function(x) as.numeric(x[2])))
      data$distance <- (data$distbegin + data$distend) / 2
    }
  data <- as.data.frame(data)
  if("quiet" %in% names(args))
    args <- c(data=list(data), formula=formula[-2], ...) else
      args <- c(data=list(data), formula=formula[-2], quiet=TRUE, ...)
  mod <- do.call(ds, args)$ddf
  prdn <- predict_esw(mod, reps=reps)
  list(ddf=mod, edd=prdn)
}

predict_esw <- function(mod, newdata=NULL, reps=999){
  data <- as.data.frame(unclass(mod$data), stringsAsFactors=T)
  key <- mod$ds$aux$ddfobj$type
  series <- mod$ds$aux$ddfobj$adjustment$series
  order <- mod$ds$aux$ddfobj$adjustment$order
  exp <- mod$ds$aux$ddfobj$adjustment$exp
  w <- mod$meta.data$width
  ispoint <- mod$ds$aux$point
  cfs <- mod$par
  names(cfs)[grep("Intercept", names(cfs))] <- "(Intercept)"
  names(cfs) <- gsub("\\.", ":", names(cfs))
  formula <- formula(mod$ds$aux$ddfobj$scale$formula)
  covnames <- all.vars(formula)
  vcov <- abs(solve(-mod$hessian))
  if(is.na(vcov[1])){
    message("Model failed to converge, no edd calculated")
    return(list(ddf=mod, edd=NULL))
  }
  
  scfs <- mvrnorm(reps, cfs, vcov)
  if(length(covnames)==0){
    ESW <- esw(exp(cfs), ispoint, key, series, order, exp, w)
    SE <- sd(apply(exp(scfs), 1, esw, ispoint, key, series, order, exp, w))
    return(data.frame(estimate=ESW, se=SE))
  } else{
    if(is.null(newdata)){
      newdata <- data %>% dplyr::select(all_of(covnames)) %>% 
        lapply(function(x) if(is.numeric(x)) mean(x, na.rm=T) else sort(unique(x)))  %>% 
        expand.grid()
    }
    for(nm in names(newdata)) 
      if(!is.numeric(newdata[[nm]]))
        newdata[[nm]] <- factor(newdata[[nm]], levels=levels(data[[nm]]))
    mat <- model.matrix(formula, newdata)
    prmat <- mat %*% cfs[colnames(mat)]
    if(key=="hr") prmat <- cbind(cfs[1], prmat)
    ESW <- apply(exp(prmat), 1, esw, ispoint, key, series, order, exp, w)
    
    prmat <- matrix(scfs[, colnames(mat)] %*% t(mat), ncol=1)
    if(key=="hr") prmat <- cbind(scfs[, 1], prmat)
    esws <- matrix(apply(prmat, 1, esw, ispoint, key, series, order, exp, w), ncol=nrow(newdata))
    SE <- apply(esws, 2, sd, na.rm=TRUE)
    data.frame(newdata, estimate=ESW, se=SE)
  }
}

AICdf <- function(mods){
  getf <- function(m){
    ff <- strsplit(m$ddf$dsmodel[2], "formula = ")[[1]][2]
    substr(ff, 1, nchar(ff)-1)
  }
  ff <- unlist(lapply(mods, getf))
  AIC <- unlist(lapply(mods, function(m) m$ddf$criterion))
  dAIC <- round(AIC-min(AIC), 2)
  AICw <- exp(-0.5*dAIC)
  AICw <- round(AICw/sum(AICw), 3)
  data.frame(ff,AIC, dAIC, AICw)[order(AIC), ]
}
