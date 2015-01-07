source("rutil/dates.R", chdir=TRUE)
source("verification.R")

getSites <- function(type="default")  {
   x = loadData(type)
   sites = unique(x[,c("SITE", "LAT", "LON", "ELEV")])
   row.names(sites) <- NULL

   return(sites)
}


loadData <- function(type="default") {
   file = getDataFile(type)
   x = read.table(file, header=TRUE)
   kens = grep("ENS", names(x))
   kens0 = kens[c(6,12,18,24,30,36,42)]
   p.thr  = 0.5

   # Add extra fields
   x$MEAN   = apply(x[,kens], 1, mean)
   x$MEDIAN = apply(x[,kens], 1, median)
   x$SD     = apply(x[,kens], 1, sd)
     N0     = apply(x[,kens] <= 0, 1, sum)
     N05    = apply(x[,kens] <= 0.5, 1, sum)
     N5     = apply(x[,kens] <= 5, 1, sum)
     NV     = apply(!is.na(x[,kens]), 1, sum)
   x$P0     = N0  / NV
   x$P05    = N05 / NV
   x$P5     = N5  / NV
   # Subensemble
   x$MEAN0  = apply(x[,kens0], 1, mean)
     N050   = apply(x[,kens0] <= 0.5, 1, sum)
     N50    = apply(x[,kens0] <= 5, 1, sum)
     NV0    = apply(!is.na(x[,kens0]), 1, sum)
   x$P050   = N050 / NV0
   # Climate fields
   sites    = unique(x$SITE)
   x$CMEAN  = mean(x$MEAN)
   x$CP05   = mean(x$P05)
   for(site in sites) {
      I = which(x$SITE == site)
      x$CMEAN[I] = mean(x$MEAN[I])
      x$CP05[I]  = mean(x$CP05[I])
   }

   return(x)
}

getDataFile <- function(type="default") {
   if(type == "default") {
      file = "data/compsData.txt"
   }
   else if(type == "arome") {
      file = "data/compsDataArome.txt"
   }
   else if(type == "manyoffsets") {
      file = "data/compsDataManyOffsets.txt"
   }
   return(file)
}

# Writes to verification file
toComps <- function(fit, x, filename) {
   library("ncdf")
   MV  = -1e30
   MVL = -999
   cdfx  <- c(0,0.2,0.5,1,5,10,20)
   cdfp  <- c(0.1,0.5,0.9)

   # Dimensions
   date <- sort(unique(x$DATE))
   offset <- sort(unique(x$OFFSET))
   location <- sort(unique(x$SITE))
   dDate <- dim.def.ncdf("Date", "", date)
   dOffset <- dim.def.ncdf("Offset", "", offset)
   dLocation <- dim.def.ncdf("Location", "", location)

   # Variables
   Iloc <- match(location, x$SITE)
   vLat <- var.def.ncdf("Lat", "degrees", dLocation, MV)#x$LAT[Iloc])
   vLon <- var.def.ncdf("Lon", "degrees", dLocation, MV)#x$LON[Iloc])
   vElev <- var.def.ncdf("Elev", "m", dLocation, MV)#x$Elev[Iloc])

   vObs  <- var.def.ncdf("obs", "", list(dLocation, dOffset, dDate), MV)
   vFcst <- var.def.ncdf("fcst", "", list(dLocation, dOffset, dDate), MV)
   vPit  <- var.def.ncdf("pit", "", list(dLocation, dOffset, dDate), MV)
   vIgn  <- var.def.ncdf("ign", "", list(dLocation, dOffset, dDate), MV)

   varList <- list(vLat, vLon, vElev, vObs, vFcst, vPit, vIgn)
   vcdf    <- NULL
   for(c in 1:length(cdfx)) {
      varname = getPVarName(cdfx[c])
      v <- var.def.ncdf(varname, "", list(dLocation, dOffset, dDate), MV)
      varList[length(varList)+1] <- list(v)
   }
   for(c in 1:length(cdfp)) {
      varname = getQVarName(cdfp[c])
      v <- var.def.ncdf(varname, "", list(dLocation, dOffset, dDate), MV)
      varList[length(varList)+1] <- list(v)
   }

   fid <- create.ncdf(filename, varList)
   close.ncdf(fid)

   # Set up data

   fid <- open.ncdf(filename, write=TRUE)
   put.var.ncdf(fid, vLat, x$LAT[Iloc])
   put.var.ncdf(fid, vLon, x$LON[Iloc])
   put.var.ncdf(fid, vElev, x$ELEV[Iloc])

   xfcst <- mG(fit, x)
   xpit  <- pG(x$OBS, fit, x)
   xign  <- -log2(dG(x$OBS, fit, x))
   xp    <- array(0, c(length(xfcst), length(cdfx)))
   xq    <- array(0, c(length(xfcst), length(cdfp)))
   for(c in 1:length(cdfx)) {
      xp[,c] <- pG(cdfx[c], fit, x)
   }
   for(c in 1:length(cdfp)) {
      xq[,c] <- qG(cdfp[c], fit, x)
   }
   #stopifnot(length(xfcst) == length(x$OBS))
   obs   <- array(MVL, c(length(location), length(offset), length(date)))
   fcst  <- array(MVL, c(length(location), length(offset), length(date)))
   pit   <- array(MVL, c(length(location), length(offset), length(date)))
   ign   <- array(MVL, c(length(location), length(offset), length(date)))
   p     <- array(MVL, c(length(location), length(offset), length(date), length(cdfx)))
   q     <- array(MVL, c(length(location), length(offset), length(date), length(cdfp)))
   for(d in 1:length(date)) {
      print(paste(d, "/", length(date), sep=""))
      Id = which(x$DATE == date[d])
      for(o in 1:length(offset)) {
         Io = which(x$OFFSET == offset[o])
         I0 = intersect(Id, Io)
         I  = match(x$SITE[I0], location)
         obs[I,o,d]  = x$OBS[I0]
         fcst[I,o,d] = xfcst[I0]
         pit[I,o,d]  = xpit[I0]
         ign[I,o,d]  = xign[I0]
         for(c in 1:length(cdfx)) {
            p[I,o,d,c] = xp[I0,c]
         }
         for(c in 1:length(cdfp)) {
            q[I,o,d,c] = xq[I0,c]
         }
      }
   }
   put.var.ncdf(fid, vObs, obs)
   put.var.ncdf(fid, vFcst, fcst)
   put.var.ncdf(fid, vPit, pit)
   put.var.ncdf(fid, vIgn, ign)

   att.put.ncdf(fid, 0, "Variable", "Precip")
   att.put.ncdf(fid, 0, "Units", "mm")

   for(c in 1:length(cdfx)) {
      varname = getPVarName(cdfx[c])
      put.var.ncdf(fid, varname, p[,,,c])
   }
   for(c in 1:length(cdfp)) {
      varname = getQVarName(cdfp[c])
      put.var.ncdf(fid, varname, q[,,,c])
   }
   close.ncdf(fid)
}
getPVarName <- function(cdfx) {
   if(cdfx < 1 & cdfx > 0) {
      varname = paste("p0", cdfx*10, sep="")
   }
   else {
      varname = paste("p", cdfx, sep="")
   }
   return(varname)
}

getQVarName <- function(cdfp) {
   varname = paste("q", floor(cdfp*100), sep="")
   return(varname)
}
writeTrainingFile <- function(x, filename) {
   offsets = seq(18,234,24)
   offsets = sort(unique(x$OFFSET))
   npar = 9
   pars    = array(0, dim=c(length(offsets), npar))
   for(i in 1:length(offsets)) {
      offset = offsets[i]
      I = which(x$OFFSET == offset)
      fit = gamlss(OBS~curt(MEAN), sigma.formula=~MEAN, nu.formula=~MEAN+P05+curt(MEAN),
                   family=ZAGA(), data=x[I,])
      pars[i, 1]   = i-1
      pars[i, 2:3] = fit$mu.coefficients
      pars[i, 4:5] = fit$sigma.coefficients
      pars[i, 6:9] = fit$nu.coefficients
   }
   # Change to missing value flag
   pars[which(is.na(pars))] = -999

   write(t(pars), filename, ncolumns=npar)
}
