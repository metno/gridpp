source("rutil/dates.R", chdir=TRUE)
# From sites
# LAT_DEC, LON_DEC, 
# From model
# VEG_FRAC, LAI, ORO_ROUGH, SST_BIN, Zs, LAND_FRAC

getRegionId <- function(regionName) {
   if(regionName == "va_ytrekyst")
      return(4)
   if(regionName == "va_kyst")
      return(5)
   if(regionName == "vv_ytrekyst")
      return(6)
   if(regionName == "vv_kyst")
      return(7)
   if(regionName == "vnn_kyst")
      return(9)
   if(regionName == "fjell")
      return(10)
   if(regionName == "kystfjell")
      return(11)
   if(regionName == "utsatte_fjellomrader")
      return(12)
   if(regionName == "andre")
      return(0)
   else
      warning(paste("Region name", regionName, "is invalid"))
   return("")

}

getSites <- function()  {

   # Cached sites
   if(exists("AROME_SITES"))
      return(AROME_SITES)

   options(stringsAsFactors=FALSE)

   ## Forecast fields

   ##  observation type (FF, FX_1, ...)
   obs.type   <- "FX_1"

   ##  list of stations
   #file.sites <- list("stations_both.txt", "stations_SST.txt", "stations_TG.txt")
   #file.names <- list("Both", "Coast", "Terrain")

   folder <- "stations"
   file.sites  <- list.files(path=folder, full.names=1, pattern=".txt$")
   file.names  <- list.files(path=folder, pattern=".txt$")
   file.names  <- gsub(".txt", "", file.names)

   ##  lead times
   prg        <- seq(0, 42, by=6)

   ##  get sites
   sites           <- lapply(file.sites, read.table, header=TRUE, sep=";", strip.white=TRUE)
   print(length(sites))
   for(i in 1:length(file.names)) {
      print(file.names[[i]])
      tmp <- file.names[[i]]
      sites[[i]]$TYPE <- tmp
      sites[[i]]$TYPEID <- getRegionId(tmp)
   }
   sites           <- do.call("rbind", sites)
   sites           <- sites[order(sites$STNR),]
   rownames(sites) <- sites$STNR

   ##  get surface information
   if(1) {
      require(miIO)
      sfc      <- list()
      sfc[[1]] <- miReadFelt("/opdata/arome_norway25/AROME_Norway_files_moved_here/preprod/AROME_Norway25_00_sfx.dat",
                             sites=sites[,c("STNR","LAT_DEC","LON_DEC")], prg=cbind(1,0),
                             prm=data.frame(VC=2,PRM=c(13196,13198,10083),LEV1=1000,LEV2=0),
                             interpolation="nearest", df=TRUE, use.names=FALSE)
      names(sfc[[1]]) <- c("SITE", "TIME", "PRG", "VEG_FRAC", "LAI", "ORO_ROUGH")
      sfc[[2]] <- miReadFelt("/opdata/arome_norway25/AROME_Norway_files_moved_here/preprod/AROME_Norway25_00_sfx.dat",
                             sites=sites[,c("STNR","LAT_DEC","LON_DEC")], prg=cbind(4,0),
                             prm=data.frame(VC=2,PRM=103,LEV1=1000,LEV2=0),
                             interpolation="nearest", df=TRUE, use.names=FALSE)
      sfc[[2]]$SST_BIN <- !is.na(sfc[[2]]$SST)
      sfc[[3]] <- miReadFelt("/opdata/arome_norway25/AROME_Norway_files_moved_here/preprod/AROME_Norway25_00.dat",
                             sites=sites[,c("STNR","LAT_DEC","LON_DEC")], prg=cbind(4,0),
                             prm=data.frame(VC=2,PRM=c(101,181),LEV1=1000,LEV2=0),
                             interpolation="nearest", df=TRUE, use.names=TRUE)
      names(sfc[[3]]) <- c("SITE", "TIME", "PRG", "Zs", "LAND_FRAC")
      sfc[[4]] <- miReadFelt("/opdata/arome_norway25/AROME_Norway_files_moved_here/preprod/AROME_Norway25_00.dat",
                             sites=sites[,c("STNR","LAT_DEC","LON_DEC")], prg=cbind(3,0),
                             prm=data.frame(VC=2,PRM=30,LEV1=1000,LEV2=0),
                             interpolation="nearest", df=TRUE, use.names=TRUE)
      tmp <- merge(sfc[[1]], sfc[[2]])
      for (i in 3:length(sfc)) 
         tmp <- merge(tmp, sfc[[i]])
      sfc <- tmp
      sfc <- sfc[,-(which(names(sfc) %in% c("TIME", "PRG")))]
      rm(tmp)

      ##  add position
      sfc <- cbind(sfc[,"SITE",drop=FALSE], sites[match(sfc$SITE, sites$STNR),], sfc[,-1])
      sfc <- sfc[order(as.numeric(sfc$SITE)),]
   }
   else {
      sfc <- sites
   }
   sfc$SITE = sfc$STNR

   # Add representativeness field
   sfc$REP <- 1
   I <- which(sfc$TYPE=="andre")
   sfc$REP[I] <- 0

   # Cache sites
   AROME_SITES <<- sfc

   return(sfc)
}

addInitTime <- function(data) {
   times <- data$TIME
   progs <- data$PRG

   dates <- as.Date(as.character(floor(data$TIME / 100)), "%Y%m%d") - progs/24 + (data$TIME %% 100)/24
   init <- as.integer(format(dates, "%Y%m%d"))
   data$INIT <- init

   return(data)
}


loadData <- function(type="default", pp=FALSE) {
   ############
   # Settings #
   ############
   options(stringsAsFactors=FALSE)
   obs.type   <- "RR_1"            # Field used as observation (FF, FX_1, ...)
   prm        <- list(RR=list(prm="precipitation_amount_acc"))
   prg        <- seq(0, 240, by=6)
   obsHours   <- 0:23
   sites      <- getSites()
   #sites = sites[51,]

   ##################
   # surface fields #
   ##################
   require(miIO)
   if(type == "default") {
      # forecastFolder <- "../../aromeWind/src/metcoopRotated/" 
      # files.sfx   <- list.files(path=forecastFolder, pattern="*2012.*.nc",recursive=TRUE,
      #                           full.names=TRUE) 
      forecastFolder <- "/vol/fou/nwparc2/eps25/"
      files <- list.files(path=forecastFolder, pattern="eps25_201401*00Z.nc",
                                recursive=TRUE, full.names=TRUE)
      forecastFolder <- "/disk1/comps/eps/"
      files <- list.files(path=forecastFolder, pattern="20140101.nc",
                                recursive=TRUE, full.names=TRUE)
   }
   else {
      error("Incorrect type")
   }
   for(file in files) {
      fc  <- miReadNetCDF(files=files, lt.format="hours", sites=sites[,c("STNR","LAT_DEC","LON_DEC")],
                               prm=prm, df=TRUE, time.var="valid", time.format="%Y%m%d%H")
      fc$TIME = as.integer(fc$TIME)
   }

   # Deaccumulate

   ################
   # observations #
   ################
   print(range(fc$TIME)%/%100)
   r = range(fc$TIME)%/%100
   obs      <- miDVH(stnr=sites$STNR, period=r, prm=obs.type, UTC=TRUE,
                     hour=seq(0,23), collapse.time=TRUE)
   obs$SITE <- obs$STNR
   obs      <- obs[,-match("STNR",names(obs))]

   # Accumulate precipitation

   # merge forecasts and observations
   x        <- merge(fc, obs, by=c("SITE","TIME"))

   #  remove missing data
   #x        <- x[rowSums(is.na(x))==0,]
   #x        <- merge(x, sites, by="SITE")#[match(x$SITE, sites$SITE),])

   #x <- addInitTime(x)

   return(x)
}

getObs <- function(dates, site, prgs) {
   require(miIO)
   source("rutil/dates.R", chdir=TRUE)
   stopifnot(max(prgs) < 24)
   minDate <- min(dates)
   maxDate <- getNextDate(max(dates))
   obs <- miDVH(stnr=site$SITE, period=c(minDate, maxDate), prm="FX_1", UTC=TRUE,
                hour=prgs, collapse.time=TRUE)
   obs$SITE <- obs$STNR
   obs$STNR <- NULL
   obs$PRG  <- obs$TIME %% 100
   obs$DATE <- floor(obs$TIME/100)
   obs$TIME <- NULL
   return(obs)
}

getElevs500 <- function() {
   lats  <- as.numeric(read.table("lats500.dat")$V1)
   lons  <- as.numeric(read.table("lons500.dat")$V1)
   elevs <- as.numeric(read.table("elevs500.dat")$V1)
   return(data.frame(LAT=lats,LON=lons,Z=elevs))
}

getAromeLatLon <- function() {
   # Get lat/lon
   lats <- as.numeric(read.table("lats.dat")$V1)
   lons <- as.numeric(read.table("lons.dat")$V1)
   return(data.frame(LAT=lats,LON=lons))
}

loadField <- function(date, prgs=seq(0,42,6)) {
   require(miIO)
   source("miReadFelt.R")
   year  <- floor(date / 10000)
   month <- floor(date / 100) %% 100
   day   <- date %% 100
   if(date >= 20131101) {
      forecastFolder <- "/vol/fou/atmos4/arome_norway25/"
      fileprefix <- "AROME_Norway25_00.dat_"
      file <- sprintf("%s%04d/%02d/%02d/%s%04d%02d%02d", forecastFolder, year, month, day, fileprefix, year, month, day, sep="")
   }
   else {
      forecastFolder <- "/vol/fou/atmos2/jakobks/HM25_spring/"
      fileprefix <- "HM25_spring_"
      file <- sprintf("%s%s%04d%02d%02d%s", forecastFolder,
                      fileprefix, year, month, day, "00.dat", sep="")
   }

   winds   <- miReadFelt(files=file, prm=cbind(2,c(33,34),1000,0), prg=prgs,
                         uv2ffdd=FALSE, df=TRUE,
                         collapse.time=TRUE)$data
   # Use patched version of miReadFelt to get elevations
   z       <- miReadFeltTN(files=file, prm=data.frame(2,c(101),1000,0), prg=data.frame(4,c(0)),
                         df=TRUE)$data
   X <- nrow(winds[,,1])   # 739
   Y <- ncol(winds[,,1])   # 949
   T <- length(prgs)

   # Get lat/lon
   lats <- as.numeric(read.table("lats.dat")$V1)
   lons <- as.numeric(read.table("lons.dat")$V1)

   data <- array(0, c(length(lats)*T, 7))
   stopifnot(length(lats) == X*Y)
   I <- 1:Y
   # Create dataframe out of matrices
   for(t in 1:T) {
      prg <- prgs[t]
      for(i in 1:X) {
         LL <- seq(i,X*Y,X)
         data[I,1] <- date
         data[I,2] <- prg
         data[I,3] <- lats[LL]
         data[I,4] <- lons[LL]
         data[I,5] <- winds[i,,2*t-1]
         data[I,6] <- winds[i,,2*t]
         data[I,7] <- z[i,,1]
         I <- I + Y
      }
   }
   data <- as.data.frame(data)
   names(data) <- c("TIME", "PRG", "LAT", "LON", "U", "V", "Z")

   # Only use an excerpt of the domain
   #I    <- which(data$LAT < 64 & data$LAT > 56.5 & data$LON >  4 & data$LON < 12)
   #data <- data[I,]
   gc()
   #save(data, file=paste("griddedAromeR/", date, ".RData", sep=""))
   return(data)
}

loadInterpolationData2 <- function(sites) {
   source("interpolation.R")
   startDate <- 20130301
   endDate   <- 20130413
   date      <- startDate
   methods   <- c("elev", "nearest")
   prgs      <- 0:23

   allData   <- NULL
   index     <- array(NA, c(nrow(sites), length(methods)))
   while(date <= endDate) {
      print(date)
      # Load data
      data <- loadField(date, prgs)

      # Interpolate
      for(s in 1:nrow(sites)) {
         site <- sites[s,]
         temp <- as.data.frame(array(NA, c(length(prgs), length(methods))))
         obs  <- array(NA,length(prgs))
         obs0 <- getObs(date, site, prgs)
         obs[match(obs0$PRG, prgs)] <- obs0$FX_1
         for(p in 1:length(prgs)) {
            prg <- prgs[p]
            I <- which(data$PRG == prg)
            currData <- data[I,]
            for(i in 1:length(methods)) {
               method <- methods[i]
               if(is.na(index[s,i])) {
                  index[s,i] <- interpolateIndex(currData, site, method)
               }
               stopifnot(!is.na(index[s,i]))
               temp[p,i] <- ws(data[I[index[s,i]],])
            }
         }
         names(temp) <- methods
         temp$DATE <- date
         temp$PRG  <- prgs
         temp$SITE <- site$SITE
         temp$OBS  <- obs
         allData <- rbind(allData, temp)
      }
      date <- getNextDate(date)
   }
   allData <- unique(allData)
   return(allData)

}

# Retrive gridded forecasts of gridpoints near sites
loadNearbyForecasts <- function(sites) {
   source("interpolation.R")
   startDate <- 20130301
   endDate   <- 20130331
   date      <- startDate
   prgs      <- 0:23

   allData   <- NULL
   S <- nrow(sites)

   I <- NULL
   while(date <= endDate) {
      print(date)
      # Load data
      data <- loadField(date, prgs)
      print("Done loading")
      if(is.null(I)) {
         print("Finding indices")
         for(s in 1:S) {
            I <- c(I,which(abs(data$LAT - sites$LAT_DEC[s]) < 0.15 & abs(data$LON - sites$LON_DEC[s]) < 0.3))
         }
      }
      allData <- rbind(allData, data[I,])
      date    <- getNextDate(date)
   }
   return(allData)
}

plotCorrelationData <- function(data, elevs500=getElevs500()) {
   require(akima)
   source("rutil/metrics.R")
   latlons <- unique(data[,c("LAT","LON", "Z")])
   corrs <- array(0, c(nrow(latlons),1))
   for(s in 1:nrow(latlons)) {
      I <- which(data$LAT == latlons$LAT[s] & data$LON == latlons$LON[s])
      #corrs[s] <- as.numeric(etsMean(data[I,"FCST"], data[I, "OBS"]))
      corrs[s] <- cor(data[I,"FCST"], data[I, "OBS"], use="complete.obs")
   }
   breaks <- seq(0.4,1,0.05)
   N <- length(breaks)
   colours <- rainbow(N)
   colI  <- cut(corrs, breaks=breaks)
   cols <- colours[colI]
   plot(latlons$LON, latlons$LAT, col=cols, pch=16)
   I <- which(elevs500$LON > min(latlons$LON) &
              elevs500$LON < max(latlons$LON) &
              elevs500$LAT > min(latlons$LAT) &
              elevs500$LAT < max(latlons$LAT))
   elevInterp <- interp(elevs500$LON[I], elevs500$LAT[I], elevs500$Z[I])
   contour(elevInterp$x, elevInterp$y, elevInterp$z, add=TRUE, col="yellow")
   elevInterp <- interp(latlons$LON, latlons$LAT, latlons$Z)
   contour(elevInterp$x, elevInterp$y, elevInterp$z, add=TRUE, col="blue", levels=seq(0,2000,100))
   legend("topleft", legend=breaks, pch=16, col=colours, bg="white")

   return(corrs)
}


plotInterpolationData <- function(data) {
   savePlot = FALSE
   sites <- unique(data$SITE)
   cex = 1
   if(savePlot)
      postscript(paste("../img/", data$SITE[1], ".ps", sep=""), height=6, width=20, pointsize=12, paper="special")

   par(mfrow=c(2,3))
   for(site in sites[1:5]) {
      I <- which(data$SITE == site)
      ylim <- c(0, 5*ceiling(max(data$elev[I])/5))
      x <- as.Date(as.character(data$DATE[I]), "%Y%m%d") + data$PRG[I]/24
      plot(x,data$elev[I], type="o", pch=16, col="green", ylim=ylim, ylab="Wind speed (m/s)", xlab="",
           xaxs="i",yaxs="i", cex=cex, cex.lab=cex,cex.axis=cex, xaxt="n")
      points(x,data$nearest[I], type="o", pch=16, col="red", cex=cex)
      points(x,data$OBS[I], type="o", pch=16, col="black", cex=cex)
      dates <- unique(as.Date(as.character(data$DATE[I]), "%Y%m%d"))
      abline(v=dates, lty=3)
      axis.Date(1, Day, at=seq(min(dates), max(dates), by="days"), cex=cex)
      legend("top", legend=c("Elevation based", "Nearest neighbour", "Observations"), lty=1, pch=16,
             col=c("green", "red", "black"), cex=1.5, bg="white")
      box()
   }
   if(savePlot)
      dev.off()
}
evalInterpolationData <- function(data) {
   savePlot = FALSE
   siteIds <- unique(data$SITE)
   sites < getSites()
   if(savePlot)
      postscript(paste("../img/", data$SITE[1], ".ps", sep=""), height=6, width=20, pointsize=12, paper="special")

   names <- array("", length(siteIds))

   par(mfrow=c(2,ceiling(length(siteIds)/2)))
   par(mar=c(5,5,5,5))
   corrs <- as.matrix(array(0,c(length(siteIds), 2)))
   for(s in 1:length(siteIds)) {
      site <- siteIds[s]
      ylim <- c(0,1)
      I <- which(data$SITE == site)
      plot(ets(data$OBS[I], data$elev[I]), type="o", pch=16, col="green", ylim=ylim, ylab="ETS",
               xlab="Wind speed (m/s)",
           xaxs="i",yaxs="i", cex=1.5, cex.lab=1.5,cex.axis=1.5, xaxt="n")
      points(ets(data$OBS[I], data$nearest[I]), type="o", pch=16, col="red", cex=1.5)
      legend("top", legend=c("Elevation based", "Nearest neighbour"), lty=1, pch=16,
             col=c("green", "red"), cex=1.5, bg="white")
      names[s] <- sites$ST_NAME[which(sites$SITE == site)]
      title(paste("Station:", names[s]))
      corrs[s,1] <- as.numeric(correlation(data$OBS[I], data$elev[I]))
      corrs[s,2] <- as.numeric(correlation(data$OBS[I], data$nearest[I]))
      box()
   }
   par(mar=c(5,11,0,1))
   barplot(t(corrs), names.arg=names, beside=TRUE, xlab="Correlation", xlim=c(0,1),
           col=c("green", "red"), horiz=TRUE, las=1)
   if(savePlot)
      dev.off()
}


loadInterpolationDataOld <- function() {
   #startDate <- 20131101
   #endDate   <- 20131113
   startDate <- 20130301
   endDate   <- 20130413
   date      <- startDate

   allData <- NULL
   while(date <= endDate) {
      print(date)
      filename <- paste("griddedAromeR/", date, ".RData", sep="")
      if(!file.exists(filename))
         loadField(date, prgs=16)#0,23,1))
      gc()
      load(filename)
      allData <- rbind(allData, data)
      date <- getNextDate(date)
   }
   return(allData)
}

getNearbySites <- function(sites) {
   lats <- sites$LAT
   lons <- sites$LON
   startDate <- 20130301
   endDate   <- 20130413
   date      <- startDate

   # Get grid points near sites
   latlon <- getAromeLatLon()
   allI <- NULL
   for(i in 1:nrow(sites)) {
      I <- which(abs(latlon$LAT - sites$LAT[i]) < 0.1 & abs(latlon$LON - sites$LON[i])< 0.1)
      allI <- rbind(allI, I)
   }
   gridSites <- data.frame(STNR=seq(1,length(allI)), LAT_DEC=latlon$LAT[allI], LON_DEC=latlon$LON[allI])
   return(gridSites)
}

loadInterpolationData <- function(sites, date) {
   require(miIO)
   #gridSites <- getNearbySites(sites)
   gridSites <- sites
   prgs = seq(0,23)

   year  <- floor(date / 10000)
   month <- floor(date / 100) %% 100
   day   <- date %% 100
   if(date >= 20131101) {
      forecastFolder <- "/vol/fou/atmos4/arome_norway25/"
      fileprefix <- "AROME_Norway25_00.dat_"
      file <- sprintf("%s%04d/%02d/%02d/%s%04d%02d%02d", forecastFolder, year, month, day, fileprefix, year, month, day, sep="")
   }
   else {
      forecastFolder <- "/vol/fou/atmos2/jakobks/HM25_spring/"
      fileprefix <- "HM25_spring_"
      file <- sprintf("%s%s%04d%02d%02d%s", forecastFolder,
                      fileprefix, year, month, day, "00.dat", sep="")
   }

   winds   <- miReadFelt(files=file, sites=gridSites, prm=cbind(2,c(33,34),1000,0), prg=prgs,
                         uv2ffdd=TRUE, df=TRUE,
                         collapse.time=TRUE)
   elev   <- miReadFelt(files=file, sites=gridSites, prm=cbind(2,c(101),1000,0), prg=data.frame(4,c(0)),
                         df=TRUE,
                         collapse.time=TRUE)
   print(elev)
   winds$LAT = gridSites$LAT_DEC[as.numeric(winds$SITE)]
   winds$LON = gridSites$LON_DEC[as.numeric(winds$SITE)]
   return(winds)
}

getCoastLine <- function() {
   sites <- getSites()
   #I         <- c(194,193,192,191,189,187,186,183,182,175,174,173,172,171,170,168,165,164,162,160,161,147,146,136,135,134,133,131,132,127,126,114,112,111,108,102,99,97,94,84,79,75,74,73,71,72,67,63,60,58,53,48,47,46,37,32,14,12)
   coastIds <- c(98790, 98580, 98550, 98400, 98090,
            96400, 96310, 94680, 94500, 92750,
            91740, 91380, 90800, 90490, 90450,
            88690, 87110, 86740, 86070, 85840,
            85890, 80610, 80102, 76530, 76450,
            76330, 75550, 75220, 75410, 71990,
            71850, 65940, 65310, 64330, 62480,
            60990, 59800, 59110, 57770, 52535,
            50500, 48330, 48120, 47350, 47260,
            47300, 44610, 43350, 42160, 41770,
            39100, 36200, 35860, 34130, 29950,
            27500, 17280, 17000)
   #I <- which(sites$SITE %in% ids)
   lats  <- NULL
   lons  <- NULL
   ids   <- NULL
   types <- NULL
   for(i in 1:length(coastIds)) {
      I <- which(sites$SITE == coastIds[i])
      if(length(I) == 0)
         print(coastIds[i])
      if(sites$REP[I] == 1) {
         lats <- rbind(lats, sites$LAT_DEC[I])
         lons <- rbind(lons, sites$LON_DEC[I])
         ids  <- rbind(ids, sites$SITE[I])
         types<- rbind(types, sites$TYPE[I])
      }
   }

   return(data.frame(SITE=ids,LAT_DEC=lats,LON_DEC=lons,TYPE=types))
}

getBroadType <- function(type) {
   if(grepl("_", type) == TRUE) {
      type  <- gsub("^.*_", "", type)
   }
   return(type)
}

