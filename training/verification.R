scatter <- function(fit, x) {
   fcst = pG(0*x$OBS+0.5, fit, x)
   plot(x$P0, fcst)
}

pithist <- function(fit, x, x0=NULL, x1=NULL, threshold=NULL, ...) {
   hist(pit(fit, x, x0=x0, x1=x1, threshold=threshold), ...)
}

# above: Use probabilities above threshold, instead of below
reliability <- function(fit, x, threshold, above=FALSE) {
   p = pG(0*x$OBS+threshold, fit, x)
   breaks = seq(0,1,0.1)
   I = cut(p, breaks, labels=FALSE)
   X = matrix(0,nrow=length(breaks)-1,ncol=1)
   Y = matrix(0,nrow=length(breaks)-1,ncol=1)
   for(i in 1:(length(breaks)-1)) {
      II = which(I == i)
      X[i] = mean(p[II])
      Y[i] = mean(x$OBS[II] <= threshold)
   }
   clim = mean(x$OBS <= threshold)
   if(above) {
      X = 1-X
      Y = 1-Y
      clim = 1 - clim
   }
   plot(X, Y, 'b', pch=19, xlim=c(0,1), ylim=c(0,1))
   lines(c(0,1),c(0,1))
   lines(c(0,1), c(clim,clim))
   lines(c(clim,clim), c(0,1))
   lines(c(0,1), c(clim/2,1-(1-clim)/2))

   # Draw confidence
   z = 1.96
   n = length(p)
   lower = Y - z*sqrt(Y*(1-Y)/n)
   upper = Y + z*sqrt(Y*(1-Y)/n)
   lines(X, lower, lty=2)
   lines(X, upper, lty=2)
}

ign <- function(fit, x, threshold=NULL) {
   if(!is.null(threshold)) {
      pdf = pG(0*x$OBS+threshold, fit, x)
      I = which(x$OBS > threshold)
      pdf[I] = 1-pdf[I]
   }
   else {
      pdf = dG(x$OBS, fit, x)
   }
   value = -log2(pdf)
   return(value)
}

brier <- function(fit, x, threshold) {
   p = pG(0*x$OBS+threshold, fit, x)
   obs = x$OBS <= threshold
   value = (obs-p)**2
   return(value)
}

pit <- function(fit, x, x0=NULL, x1=NULL, threshold=NULL) {
   if(is.null(threshold)) {
      cdf = pG(x$OBS, fit, x)
      if(!is.null(x0)) {
         I = which(x$OBS <= x0)
         cdf = runif(length(I), 0, 1)*cdf[I]
      }
      if(!is.null(x1)) {
         I = which(x$OBS >= x1)
         cdf = 1 - runif(length(I), 0, 1)*(1 - cdf[I])
      }

   }
   else {
      cdf = pG(x$OBS*0+threshold, fit, x)
      I = which(x$OBS < threshold)
      cdf[I]  = runif(length(I), 0, 1)*cdf[I]
      cdf[!I] = 1 - runif(length(!I), 0, 1)*(1-cdf[!I])
   }
   return(cdf)
}

pG <- function(p, fit, x) {
   if(length(p) == 1)
      p = 0*x$OBS + p
   return(getValues(p, fit, x, "p"))
}
qG <- function(q, fit, x) {
   if(length(q) == 1)
      q = 0*x$OBS + q
   values = getValues(q, fit, x, "q")
   return(values)
}
dG <- function(d, fit, x) {
   if(length(d) == 1)
      d = 0*x$OBS + d
   return(getValues(d, fit, x, "d"))
}
mG <- function(fit, x) {
   return(getValues(0, fit, x, "m"))
}
getValues <- function(q, fit, x, type="q") {
   if(length(q) != dim(x)[1]) {
      #stop("getValues: q must be the same size as x")
   }
   if(is.character(fit)) {
      if(fit == "bpe") {
         values = matrix(NA, nrow=dim(x)[1], ncol=1)
         kens = grep("ENS", names(x))
         if(type == "p") {
            values = apply(x[,kens] <= q, 1, mean)
         }
         else if(type == "q") {
            temp  = t(apply(x[,kens], 1, sort))
            Nrows = dim(temp)[1]
            Nens  = dim(temp)[2]
            k     = floor(Nens*q)
            ind   = (k-1)*Nrows + (1:Nrows)
            values[,1] = temp[ind]
         }
         else if(type == "d") {
            warning("PDF not implemented for bpe")
            values[,1] = -999
         }
         else if(type == "m") {
            values[,1] = apply(x[,kens], 1, mean)
         }
      }
   }
   else {
      family = fit$family[1]
      par = predictAll(fit, x)
      mu = par$mu
      sigma = par$sigma
      nu = par$nu
      tau = par$tau
      if(type == "m") {
         values = matrix(NA, nrow=dim(x)[1], ncol=1)
         values[,1] = mu*(1-nu)
         return(values)
      }
      # print(paste(mu,sigma,nu,tau))
      string = paste(type, family, "(q, mu=mu, sigma=sigma", sep="")
      if(!is.null(nu)) {
         string = paste(string, ",nu=nu", sep="")
      }
      if(!is.null(tau))
         string = paste(string, ",tau=tau", sep="")
      string = paste(string, ")", sep="")
      values = eval(parse(text=string))
   }
   return(values)
}

mapplot <- function(x, y, z, ...) {
   loc = unique(data.frame(x=x,y=y))
   X = loc$x
   Y = loc$y
   Z = matrix(0, nrow=dim(loc)[1], ncol=1)
   for(i in 1:dim(loc)[1]) {
      I = which((x == X[i]) & (y == Y[i]))
      Z[i] = mean(z[I])
   }
   breaks = quantile(Z, seq(0,1,0.1))
   cols = rainbow(length(breaks))
   I = cut(Z, breaks)

   plot(X, Y, pch=19, "p", col=cols[I], ..., xlab="", ylab="")
   legend("topleft", legend=breaks, col=cols, pch=19)
}

cutplot <- function(x, y, ...) {
   edges = unique(quantile(x, seq(0,1,0.01), labels=FALSE))
   X = (edges[2:length(edges)]+edges[1:(length(edges)-1)])/2
   #X = seq(0,100,5)
   #X = unique(x)
   I = cut(x, edges, labels=FALSE)
   Y = matrix(0, nrow=length(X), ncol=1)
   for(i in 1:length(X)) {
      II = which(I == i)
      print(length(II))
      Y[i] = mean(y[II])
   }

   plot(x,y,'p', pch=19)
   lines(X, Y, 'b', ...)
}
