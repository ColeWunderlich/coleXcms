#had to modify MPI
xcmsSet2 <-
function(files = NULL, snames = NULL, sclass = NULL, phenoData = NULL,
                    profmethod = "bin", profparam = list(),
                    polarity = NULL, lockMassFreq=FALSE,
                    mslevel=NULL, nSlaves=0, progressCallback=NULL,
                    scanrange=NULL, ...) {

    object <- new("xcmsSet2")

    ## initialise progress information
    xcms.options <- getOption("BioC")$xcms
    xcms.methods <- c(paste("group", xcms.options$group.methods,sep="."), paste("findPeaks", xcms.options$findPeaks.methods,sep="."),
                      paste("retcor", xcms.options$retcor.methods,sep="."), paste("fillPeaks", xcms.options$fillPeaks.methods,sep="."))
    eval(parse(text=paste("object@progressInfo <- list(",paste(xcms.methods,"=0",sep="",collapse=","),")") ))

    if (is.function(progressCallback))
        object@progressCallback <- progressCallback

    filepattern <- c("[Cc][Dd][Ff]", "[Nn][Cc]", "([Mm][Zz])?[Xx][Mm][Ll]",
                     "[Mm][Zz][Dd][Aa][Tt][Aa]", "[Mm][Zz][Mm][Ll]")
    filepattern <- paste(paste("\\.", filepattern, "$", sep = ""), collapse = "|")

    if (is.null(files))
        files <- getwd()
    info <- file.info(files)
    listed <- list.files(files[info$isdir], pattern = filepattern,
                         recursive = TRUE, full.names = TRUE)
    files <- c(files[!info$isdir], listed)

    ## try making paths absolute
    files_abs <- file.path(getwd(), files)
    exists <- file.exists(files_abs)
    files[exists] <- files_abs[exists]

    if(lockMassFreq==TRUE){
        ## remove the 02 files if there here
        lockMass.files<-grep("02.CDF", files)
        if(length(lockMass.files) > 0){
            files<-files[-lockMass.files]
        }
    }

    filepaths(object) <- files

    if (length(files) == 0)
        stop("No NetCDF/mzXML/mzData/mzML files were found.\n")

    ## determine experimental design
    fromPaths <- phenoDataFromPaths(files)
    if (is.null(snames)) {
        snames <- rownames(fromPaths)
    } else {
        rownames(fromPaths) <- snames
    }

    pdata <- phenoData
    if (is.null(pdata)) {
        pdata <- sclass
        if (is.null(pdata))
            pdata <- fromPaths
    }
    phenoData(object) <- pdata
    if (is.null(phenoData))
        rownames(phenoData(object)) <- snames

    rtlist <- list(raw = vector("list", length(snames)),
                   corrected = vector("list", length(snames)))

    if ("step" %in% names(profparam)) {
        if ("step" %in% names(list(...)) && profparam$step != list(...)$step) {
            stop("different step values defined in profparam and step arguments")
        }
        profstep <- profparam$step
        profparam <- profparam[names(profparam) != "step"]
    } else if ("step" %in% names(list(...))) {
        profstep <- list(...)$step
    } else {
        profstep <- 0.1
    }

    if ("method" %in% names(profparam)) {
        if (profparam$method != profmethod) {
            stop("different method values defined in profparam and profmethod arguments")
        }
        profmethod <- profparam$method
        profparam <- profparam[names(profparam) != "method"]
    }

    profinfo(object) <- c(list(method = profmethod, step = profstep), profparam)

    object@polarity <- as.character(polarity)
    includeMSn=FALSE

    ## implicitely TRUE if selecting MSn
    includeMSn <- !is.null(mslevel) &&  mslevel>1

    ## implicitely TRUE if MS1 parent peak picking
    xcmsSetArgs <- as.list(match.call())
    if (!is.null(xcmsSetArgs$method)) {
        if (xcmsSetArgs$method=="MS1") {
            includeMSn=TRUE
        }
    }

    parmode <- xcmsParallelSetup(nSlaves=nSlaves)
    runParallel <- parmode$runParallel
    parMode <- parmode$parMode    
    snowclust <- parmode$snowclust
    
        params <- list(...);
        params$profmethod <- profmethod;
        params$profparam <- profparam;
        params$includeMSn <- includeMSn;
        params$scanrange <- scanrange;

        params$mslevel <- mslevel; ## Actually, this is 
        params$lockMassFreq <- lockMassFreq;

        ft <- cbind(file=files,id=1:length(files))
        argList <- apply(ft,1,function(x) list(file=x["file"],id=as.numeric(x["id"]),params=params))

        if (parMode == "MPI") {
            res <- xcmsPapply(argList, findPeaksPar)
            mpi.close.Rslaves()
        } else if (parMode == "SOCK") {
                res <- xcmsClusterApply(cl=snowclust, x=argList, fun=findPeaksPar, msgfun=msgfun.featureDetection)
                stopCluster(snowclust)
            } else {
              ## serial mode
              res <- lapply(argList, findPeaksPar)
            }

        peaklist <- lapply(res, function(x) x$peaks)
        rtlist$raw <-  rtlist$corrected <-  lapply(res, function(x) x$scantime)
        if(lockMassFreq){
            object@dataCorrection[1:length(files)]<-1
        }

    lapply(1:length(peaklist), function(i) {
        if (is.null(peaklist[[i]]))
            warning("No peaks found in sample ", snames[i], call. = FALSE)
        else  if (nrow(peaklist[[i]]) == 0)
            warning("No peaks found in sample ", snames[i], call. = FALSE)
        else if (nrow(peaklist[[i]]) == 1)
            warning("Only 1 peak found in sample ", snames[i], call. = FALSE)
        else if (nrow(peaklist[[i]]) < 10)
            warning("Only ", nrow(peaklist[[i]]), " peaks found in sample",
                    snames[i], call. = FALSE)
    })

    peaks(object) <- do.call(rbind, peaklist)
    object@rt <- rtlist
	object@xcmsRaw <- lapply(res, function(x) x$xRaw) #store xcmsRaw objects

    object
}


setMethod("group.density", "xcmsSet2", function(object, bw = 30, minfrac = 0.5, minsamp = 1,
                                               mzwid = 0.25, max = 50, sleep = 0) {

    samples <- sampnames(object)
    classlabel <- sampclass(object)
    classnames <- as.character(unique(sampclass(object)))
    classlabel <- as.vector(unclass(classlabel))
    classnum <- table(classlabel)

    peakmat <- peaks(object)
    porder <- order(peakmat[,"mz"])
    peakmat <- peakmat[porder,, drop=FALSE]
    rownames(peakmat) <- NULL
    retrange <- range(peakmat[,"rt"])

    mass <- seq(peakmat[1,"mz"], peakmat[nrow(peakmat),"mz"] + mzwid, by = mzwid/2)
    masspos <- findEqualGreaterM(peakmat[,"mz"], mass)

    groupmat <- matrix(nrow = 512, ncol = 7 + length(classnum))
    groupindex <- vector("list", 512)

    endidx <- 0
    num <- 0
    gcount <- integer(length(classnum))
    for (i in seq(length = length(mass)-2)) {
        if (i %% 500 == 0) {
            cat(round(mass[i]), "")
            flush.console()
        }
        startidx <- masspos[i]
        endidx <- masspos[i+2]-1
        if (endidx - startidx < 0)
            next
        speakmat <- peakmat[startidx:endidx,,drop=FALSE]
        den <- density(speakmat[,"rt"], bw, from = retrange[1]-3*bw, to = retrange[2]+3*bw)
        maxden <- max(den$y)
        deny <- den$y
        gmat <- matrix(nrow = 5, ncol = 2+length(classnum))
        snum <- 0
        while (deny[maxy <- which.max(deny)] > maxden/20 && snum < max) {
            grange <- descendMin(deny, maxy)
            deny[grange[1]:grange[2]] <- 0
            gidx <- which(speakmat[,"rt"] >= den$x[grange[1]] & speakmat[,"rt"] <= den$x[grange[2]])
            gnum <- classlabel[unique(speakmat[gidx,"sample"])]
            for (j in seq(along = gcount))
                gcount[j] <- sum(gnum == j)
            if (! any(gcount >= classnum*minfrac & gcount >= minsamp))
                next
            snum <- snum + 1
            num <- num + 1
### Double the size of the output containers if they're full
            if (num > nrow(groupmat)) {
                groupmat <- rbind(groupmat, matrix(nrow = nrow(groupmat), ncol = ncol(groupmat)))
                groupindex <- c(groupindex, vector("list", length(groupindex)))
            }
            groupmat[num, 1] <- median(speakmat[gidx, "mz"])
            groupmat[num, 2:3] <- range(speakmat[gidx, "mz"])
            groupmat[num, 4] <- median(speakmat[gidx, "rt"])
            groupmat[num, 5:6] <- range(speakmat[gidx, "rt"])
            groupmat[num, 7] <- length(gidx)
            groupmat[num, 7+seq(along = gcount)] <- gcount
            groupindex[[num]] <- sort(porder[(startidx:endidx)[gidx]])
        }
        if (sleep > 0) {
            plot(den, main = paste(round(min(speakmat[,"mz"]), 2), "-", round(max(speakmat[,"mz"]), 2)))
            for (i in seq(along = classnum)) {
                idx <- classlabel[speakmat[,"sample"]] == i
                points(speakmat[idx,"rt"], speakmat[idx,"into"]/max(speakmat[,"into"])*maxden, col = i, pch=20)
            }
            for (i in seq(length = snum))
                abline(v = groupmat[num-snum+i, 5:6], lty = "dashed", col = i)
            Sys.sleep(sleep)
        }
    }
    cat("\n")

    colnames(groupmat) <- c("mzmed", "mzmin", "mzmax", "rtmed", "rtmin", "rtmax",
                            "npeaks", classnames)

    groupmat <- groupmat[seq(length = num),,drop=FALSE]
    groupindex <- groupindex[seq(length = num)]

    ## Remove groups that overlap with more "well-behaved" groups
    numsamp <- rowSums(groupmat[,(match("npeaks", colnames(groupmat))+1):ncol(groupmat),drop=FALSE])
    uorder <- order(-numsamp, groupmat[,"npeaks"])

    uindex <- rectUnique(groupmat[,c("mzmin","mzmax","rtmin","rtmax"),drop=FALSE],
                         uorder)

    groups(object) <- groupmat[uindex,,drop=FALSE]
    groupidx(object) <- groupindex[uindex]

    object
})


setMethod("group", "xcmsSet2", function(object, method=getOption("BioC")$xcms$group.method,
                                       ...) {

    method <- match.arg(method, getOption("BioC")$xcms$group.methods)
    if (is.na(method))
        stop("unknown method : ", method)
    method <- paste("group", method, sep=".")
    invisible(do.call(method, alist(object, ...)))
})

#had to modify MPI
setMethod("fillPeaks.chrom", "xcmsSet2", function(object, nSlaves=NULL) {
  ## development mockup:
  if (FALSE) {
    library(xcms)
    library(faahKO)
    object <- group(faahko)
    gf <- fillPeaks(object)
    pkgEnv = getNamespace("xcms")
    attach(pkgEnv)
  }
  
    peakmat <- peaks(object)
    groupmat <- groups(object)
    if (length(groupmat) == 0)
        stop("No group information found")
    files <- filepaths(object)
    samp <- sampnames(object)
    classlabel <- as.vector(unclass(sampclass(object)))
    prof <- profinfo(object)
    rtcor <- object@rt$corrected

    ## Remove groups that overlap with more "well-behaved" groups
    numsamp <- rowSums(groupmat[,(match("npeaks", colnames(groupmat))+1):ncol(groupmat),drop=FALSE])
    uorder <- order(-numsamp, groupmat[,"npeaks"])
    uindex <- rectUnique(groupmat[,c("mzmin","mzmax","rtmin","rtmax"),drop=FALSE],
                         uorder)
    groupmat <- groupmat[uindex,]
    groupindex <- groupidx(object)[uindex]
    gvals <- groupval(object)[uindex,]

    peakrange <- matrix(nrow = nrow(gvals), ncol = 4)
    colnames(peakrange) <- c("mzmin","mzmax","rtmin","rtmax")

    mzmin <- peakmat[gvals,"mzmin"]
    dim(mzmin) <- c(nrow(gvals), ncol(gvals))
    peakrange[,"mzmin"] <- apply(mzmin, 1, median, na.rm = TRUE)
    mzmax <- peakmat[gvals,"mzmax"]
    dim(mzmax) <- c(nrow(gvals), ncol(gvals))
    peakrange[,"mzmax"] <- apply(mzmax, 1, median, na.rm = TRUE)
    retmin <- peakmat[gvals,"rtmin"]
    dim(retmin) <- c(nrow(gvals), ncol(gvals))
    peakrange[,"rtmin"] <- apply(retmin, 1, median, na.rm = TRUE)
    retmax <- peakmat[gvals,"rtmax"]
    dim(retmax) <- c(nrow(gvals), ncol(gvals))
    peakrange[,"rtmax"] <- apply(retmax, 1, median, na.rm = TRUE)

    lastpeak <- nrow(peakmat)
    lastpeakOrig <- lastpeak

##    peakmat <- rbind(peakmat, matrix(nrow = sum(is.na(gvals)), ncol = ncol(peakmat)))

    cnames <- colnames(object@peaks)

# Making gvals environment so that when it is repeated for each file it only uses the memory one time
gvals_env <- new.env(parent=baseenv())
assign("gvals", gvals, envir = gvals_env)


    ft <- cbind(file=files,id=1:length(files))
    argList <- apply(ft,1,function(x) {
      ## Add only those samples which actually have NA in them
      if (!any(is.na(gvals[,as.numeric(x["id"])]))) {
        ## nothing to do.
        list()
      } else {
        list(file=x["file"],id=as.numeric(x["id"]), object = object,
             params=list(method="chrom",
               gvals=gvals_env, 
               prof=prof,
               dataCorrection=object@dataCorrection,
               polarity=object@polarity,
               rtcor=object@rt$corrected[[as.numeric(x["id"])]],
               peakrange=peakrange))
      }
    })

  nonemptyIdx <- (sapply(argList, length) > 0)

  if (!any(nonemptyIdx)) {
    ## Nothing to do
    return(invisible(object))
  }
    
  argList <- argList[nonemptyIdx]
    
  parmode <- xcmsParallelSetup(nSlaves=nSlaves)
  runParallel <- parmode$runParallel
  parMode <- parmode$parMode    
  snowclust <- parmode$snowclust

  if (parMode == "MPI") {
    newpeakslist <- xcmsPapply(argList, fillPeaksChromPar)
    mpi.close.Rslaves()
  } else if (parMode == "SOCK") {
    newpeakslist <- xcmsClusterApply(cl=snowclust, x=argList,
                                     fun=fillPeaksChromPar,
                                     msgfun=msgfunGeneric)    
    stopCluster(snowclust)
  } else {
    ## serial mode
    newpeakslist <- lapply(argList, fillPeaksChromPar)
  }



  o <- order(sapply(newpeakslist, function(x) x$myID))
  newpeaks <- do.call(rbind, lapply(newpeakslist[o], function(x) x$newpeaks))

  ## Make sure colnames are compatible 
  newpeaks <- newpeaks[, match(cnames, colnames(newpeaks)), drop = FALSE]
  colnames(newpeaks) <- cnames

  peakmat <- rbind(peakmat, newpeaks)

  for (i in seq(along = files)) {
    naidx <- which(is.na(gvals[,i]))
    
    for (j in seq(along = naidx))
      groupindex[[naidx[j]]] <- c(groupindex[[naidx[j]]], lastpeak+j)

    lastpeak <- lastpeak + length(naidx)
  }

    peaks(object) <- peakmat
    object@filled <- seq((lastpeakOrig+1),nrow(peakmat))
    groups(object) <- groupmat
    groupidx(object) <- groupindex

    invisible(object)
})


setMethod("fillPeaks", "xcmsSet2", function(object, method=getOption("BioC")$xcms$fillPeaks.method,...) {
    method <- match.arg(method, getOption("BioC")$xcms$fillPeaks.methods)
    if (is.na(method))
        stop("unknown method : ", method)
    method <- paste("fillPeaks", method, sep=".")
    invisible(do.call(method, alist(object, ...)))
})


setMethod("getEIC", "xcmsSet2", function(object, mzrange, rtrange = 200,
                                        groupidx, sampleidx = sampnames(object),
                                        rt = c("corrected", "raw")) {

    files <- filepaths(object)
    grp <- groups(object)
    samp <- sampnames(object)
    prof <- profinfo(object)

    rt <- match.arg(rt)

    if (is.numeric(sampleidx))
        sampleidx <- sampnames(object)[sampleidx]
    sampidx <- match(sampleidx, sampnames(object))

    if (!missing(groupidx)) {
        if (is.numeric(groupidx))
            groupidx <- groupnames(object)[unique(as.integer(groupidx))]
        grpidx <- match(groupidx, groupnames(object, template = groupidx))
    }

    if (missing(mzrange)) {
        if (missing(groupidx))
            stop("No m/z range or groups specified")
        if (any(is.na(groupval(object, value = "mz"))))
            stop('Please use fillPeaks() to fill up NA values !')
        mzmin <- -rowMax(-groupval(object, value = "mzmin"))
        mzmax <- rowMax(groupval(object, value = "mzmax"))
        mzrange <- matrix(c(mzmin[grpidx], mzmax[grpidx]), ncol = 2)
    } else if (all(c("mzmin","mzmax") %in% colnames(mzrange)))
        mzrange <- mzrange[,c("mzmin", "mzmax"),drop=FALSE]
    else if (is.null(dim(mzrange)))
        stop("mzrange must be a matrix")
    colnames(mzrange) <- c("mzmin", "mzmax")

    if (length(rtrange) == 1) {
        if (missing(groupidx))
            rtrange <- matrix(rep(range(object@rt[[rt]][sampidx]), nrow(mzrange)),
                              ncol = 2, byrow = TRUE)
        else {
            rtrange <- retexp(grp[grpidx,c("rtmin","rtmax"),drop=FALSE], rtrange)
        }
    } else if (is.null(dim(rtrange)))
        stop("rtrange must be a matrix or single number")
    colnames(rtrange) <- c("rtmin", "rtmax")

    if (missing(groupidx))
        gnames <- character(0)
    else
        gnames <- groupidx

    eic <- vector("list", length(sampleidx))
    names(eic) <- sampleidx

    for (i in seq(along = sampidx)) {

        cat(sampleidx[i], "")
        flush.console()
        lcraw <- object@xcmsRaw[[i]]  #Need to make sure order is preserved!!! (ie. this is referencing the sample it should)
        if(length(object@dataCorrection) > 1){
            if(object@dataCorrection[i] == 1)
                lcraw<-stitch(lcraw, AutoLockMass(lcraw))
        }
        if (rt == "corrected")
            lcraw@scantime <- object@rt$corrected[[sampidx[i]]]
        if (length(prof) > 2)
            lcraw@profparam <- prof[seq(3, length(prof))]
        currenteic <- getEIC(lcraw, mzrange, rtrange, step = prof$step)
        eic[[i]] <- currenteic@eic[[1]]
        rm(lcraw)
        gc()
    }
    cat("\n")

    invisible(new("xcmsEIC", eic = eic, mzrange = mzrange, rtrange = rtrange,
                  rt = rt, groupnames = gnames))
})




setMethod("retcor.obiwarp", "xcmsSet2", function(object, plottype = c("none", "deviation"),
                                                profStep=1, center=NULL,
                                                col = NULL, ty = NULL,
                                                response=1, distFunc="cor_opt",
                                                gapInit=NULL, gapExtend=NULL,
                                                factorDiag=2, factorGap=1,
                                                localAlignment=0, initPenalty=0) {

    if (is.null(gapInit)) {
        if (distFunc=="cor") {gapInit=0.3}
        if (distFunc=="cor_opt") {gapInit=0.3}
        if (distFunc=="cov") {gapInit=0.0}
        if (distFunc=="euc") {gapInit=0.9}
        if (distFunc=="prd") {gapInit=0.0}
    }

    if (is.null(gapExtend)) {
        if (distFunc=="cor") {gapExtend=2.4}
        if (distFunc=="cor_opt") {gapExtend=2.4}
        if (distFunc=="cov") {gapExtend= 11.7}
        if (distFunc=="euc") {gapExtend= 1.8}
        if (distFunc=="prd") {gapExtend= 7.8}
    }

    peakmat <- peaks(object)
    samples <- sampnames(object)
    classlabel <- as.vector(unclass(sampclass(object)))
    N <- length(samples)
    corpeaks <- peakmat
    plottype <- match.arg(plottype)

    if (length(object@rt) == 2) {
        rtcor <- object@rt$corrected
    } else {
        fnames <- filepaths(object)
        rtcor <- vector("list", length(fnames))
        for (i in seq(along = fnames)) {
            xraw <- object@xcmsRaw[[i]]
            rtcor[[i]] <- xraw@scantime
        }
        object@rt <- list(raw = rtcor, corrected = rtcor)
    }

    rtimecor <- vector("list", N)
    rtdevsmo <- vector("list", N)
    plength <- rep(0, N)

    if (missing(center)) {
        for(i in 1:N){
            plength[i] <- length(which(peakmat[,"sample"]==i))
        }
        center <- which.max(plength)
    }

    cat("center sample: ", samples[center], "\nProcessing: ")
    idx <- which(seq(1,N) != center)
    obj1 <- object@xcmsRaw[[center]] #xcmsRaw(object@filepaths[center], profmethod="bin", profstep=0)
	
	## added t automatically find the correct scan range from the xcmsSet object
	if(length(obj1@scantime) != length(object@rt$raw[[center]])){
		##figure out the scan time range
		scantime.start	<-object@rt$raw[[center]][1]
		scantime.end	<-object@rt$raw[[center]][length(object@rt$raw[[center]])]
		
		scanrange.start	<-which.min(abs(obj1@scantime - scantime.start)) 
		scanrange.end	<-which.min(abs(obj1@scantime - scantime.end))
		scanrange<-c(scanrange.start, scanrange.end)
		obj1 <- xcmsRaw(object@filepaths[center], profmethod="bin", profstep=0, scanrange=scanrange) ###May need to change this 
	} else{
		scanrange<-NULL	
	}

    for (si in 1:length(idx)) {
        s <- idx[si]
        cat(samples[s], " ")

        profStepPad(obj1) <- profStep ## (re-)generate profile matrix, since it might have been modified during previous iteration
		if(is.null(scanrange)){
			obj2 <- object@xcmsRaw[[s]] #xcmsRaw(object@filepaths[s], profmethod="bin", profstep=0)	
		} else{
			obj2 <- xcmsRaw(object@filepaths[s], profmethod="bin", profstep=0, scanrange=scanrange)
		}
        profStepPad(obj2) <- profStep ## generate profile matrix

        mzmin <-  min(obj1@mzrange[1], obj2@mzrange[1])
        mzmax <-  max(obj1@mzrange[2], obj2@mzrange[2])

        mz <- seq(mzmin,mzmax, by=profStep)
        mz <- as.double(mz)
        mzval <- length(mz)

        scantime1 <- obj1@scantime
        scantime2 <- obj2@scantime

        mstdiff <- median(c(diff(scantime1), diff(scantime2)))

        rtup1 <- c(1:length(scantime1))
        rtup2 <- c(1:length(scantime2))

        mst1 <- which(diff(scantime1)>5*mstdiff)[1]
        if(!is.na(mst1)) {
            rtup1 <- which(rtup1<=mst1)
            cat("Found gaps: cut scantime-vector at ", scantime1[mst1],"seconds", "\n")
        }

        mst2 <- which(diff(scantime2)>5*mstdiff)[1]
        if(!is.na(mst2)) {
            rtup2 <- which(rtup2<=mst2)
            cat("Found gaps: cut scantime-vector at ", scantime2[mst2],"seconds", "\n")
        }

        scantime1 <- scantime1[rtup1]
        scantime2 <- scantime2[rtup2]

        rtmaxdiff <- abs(diff(c(scantime1[length(scantime1)],
                                scantime2[length(scantime2)])))
        if(rtmaxdiff>(5*mstdiff)){
            rtmax <- min(scantime1[length(scantime1)],
                         scantime2[length(scantime2)])
            rtup1 <- which(scantime1<=rtmax)
            rtup2 <- which(scantime2<=rtmax)
        }

        scantime1 <- scantime1[rtup1]
        scantime2 <- scantime2[rtup2]
        valscantime1 <- length(scantime1)
        valscantime2 <- length(scantime2)

        if(length(obj1@scantime)>valscantime1) {
            obj1@env$profile <- obj1@env$profile[,-c((valscantime1+1):length(obj1@scantime))]
        }
        if(length(obj2@scantime)>valscantime2) {
            obj2@env$profile <- obj2@env$profile[,-c((valscantime2+1):length(obj2@scantime))]
        }

        if(mzmin < obj1@mzrange[1]) {
            seqlen <- length(seq(mzmin, obj1@mzrange[1], profStep))-1
            x <- matrix(0, seqlen,dim(obj1@env$profile)[2])
            obj1@env$profile <- rbind(x, obj1@env$profile)
        }
        if(mzmax > obj1@mzrange[2]){
            seqlen <- length(seq(obj1@mzrange[2], mzmax, profStep))-1
            x <- matrix(0, seqlen, dim(obj1@env$profile)[2])
            obj1@env$profile <- rbind(obj1@env$profile, x)
        }
        if(mzmin < obj2@mzrange[1]){
            seqlen <- length(seq(mzmin, obj2@mzrange[1], profStep))-1
            x <- matrix(0, seqlen, dim(obj2@env$profile)[2])
            obj2@env$profile <- rbind(x, obj2@env$profile)
        }
        if(mzmax > obj2@mzrange[2]){
            seqlen <- length(seq(obj2@mzrange[2], mzmax, profStep))-1
            x <- matrix(0, seqlen, dim(obj2@env$profile)[2])
            obj2@env$profile <- rbind(obj2@env$profile, x)
        }

        intensity1 <- obj1@env$profile
        intensity2 <- obj2@env$profile

        if ((mzval * valscantime1 != length(intensity1)) ||  (mzval * valscantime2 != length(intensity2)))
            stop("Dimensions of profile matrices do not match !\n")

        rtimecor[[s]] <-.Call("R_set_from_xcms",
                              valscantime1,scantime1,mzval,mz,intensity1,
                              valscantime2,scantime2,mzval,mz,intensity2,
                              response, distFunc,
                              gapInit, gapExtend,
                              factorDiag, factorGap,
                              localAlignment, initPenalty)

        if(length(obj2@scantime) > valscantime2) {
            object@rt$corrected[[s]] <- c(rtimecor[[s]],
                                          obj2@scantime[(max(rtup2)+1):length(obj2@scantime)])
        } else {
            object@rt$corrected[[s]] <- rtimecor[[s]]
        }

        rtdevsmo[[s]] <- round(rtcor[[s]]-object@rt$corrected[[s]],2)

        rm(obj2)
        gc()

        ## updateProgressInfo
        object@progressInfo$retcor.obiwarp <-  si / length(idx)
        xcms:::progressInfoUpdate(object)

    }

    cat("\n")
    rtdevsmo[[center]] <- round(rtcor[[center]] - object@rt$corrected[[center]], 2)

    if (plottype == "deviation") {

        ## Set up the colors and line type
        if (missing(col)) {
            col <- integer(N)
            for (i in 1:max(classlabel))
                col[classlabel == i] <- 1:sum(classlabel == i)
        }
        if (missing(ty)) {
            ty <- integer(N)
            for (i in 1:max(col))
                ty[col == i] <- 1:sum(col == i)
        }
        if (length(palette()) < max(col))
            mypal <- rainbow(max(col), end = 0.85)
        else
            mypal <- palette()[1:max(col)]

        rtrange <- range(do.call("c", rtcor))
        devrange <- range(do.call("c", rtdevsmo))

        layout(matrix(c(1, 2),ncol=2,  byrow=F),widths=c(1,0.3))
        par(mar=c(4,4,2,0))

        plot(0, 0, type="n", xlim = rtrange, ylim = devrange,
             main = "Retention Time Deviation vs. Retention Time",
             xlab = "Retention Time", ylab = "Retention Time Deviation")

        for (i in 1:N) {
            points(rtcor[[i]], rtdevsmo[[i]], type="l", col = mypal[col[i]], lty = ty[i])
        }

        plot.new() ;  par(mar= c(2, 0, 2, 0))
        plot.window(c(0,1), c(0,1))
        legend(0,1.04, basename(samples), col = mypal[col], lty = ty)
    }

    for (i in 1:N) {
        cfun <- stepfun(rtcor[[i]][-1] - diff(rtcor[[i]])/2, rtcor[[i]] - rtdevsmo[[i]])
        rtcor[[i]] <- rtcor[[i]] - rtdevsmo[[i]]

        sidx <- which(corpeaks[,"sample"] == i)
        corpeaks[sidx, c("rt", "rtmin", "rtmax")] <- cfun(corpeaks[sidx, c("rt", "rtmin", "rtmax")])
    }

    peaks(object) <- corpeaks
    groups(object) <- matrix(nrow = 0, ncol = 0)
    groupidx(object) <- list()
    invisible(object)

})


setMethod("retcor", "xcmsSet2", function(object, method=getOption("BioC")$xcms$retcor.method,
                                        ...) {

    ## Backward compatibility for old "methods"
    if (method == "linear" || method == "loess") {
        return(invisible(do.call(retcor.peakgroups, alist(object, smooth=method, ...))))
    }

    method <- match.arg(method, getOption("BioC")$xcms$retcor.methods)
    if (is.na(method))
        stop("unknown method : ", method)
    method <- paste("retcor", method, sep=".")

    invisible(do.call(method, alist(object, ...)))
})


setMethod("peakTable", "xcmsSet2", function(object, filebase = character(), ...) {

    if (length(sampnames(object)) == 1) {
        return(object@peaks)
    }

    if (nrow(object@groups) < 1) {
        stop ('First argument must be an xcmsSet with group information or contain only one sample.')
    }

    groupmat <- groups(object)


    if (! "value" %in% names(list(...))) {
        ts <- data.frame(cbind(groupmat,groupval(object, value="into",  ...)), row.names = NULL)
    } else {
        ts <- data.frame(cbind(groupmat,groupval(object, ...)), row.names = NULL)
    }

    cnames <- colnames(ts)

    if (cnames[1] == 'mzmed') {
        cnames[1] <- 'mz'
    } else {
        stop ('mzmed column missing')
    }
    if (cnames[4] == 'rtmed') {
        cnames[4] <- 'rt'
    } else {
        stop ('mzmed column missing')
    }

    colnames(ts) <- cnames

    if (length(filebase))
        write.table(ts, paste(filebase, ".tsv", sep = ""), quote = FALSE, sep = "\t", col.names = NA)

    ts
})


##Dependencies for obiwarp
#setGeneric("profStepPad<-", function(object, value) standardGeneric("profStepPad<-"))

#setReplaceMethod("profStepPad", "xcmsRaw", function(object, value) {

   # if ("profile" %in% ls(object@env))
   #     rm("profile", envir = object@env)
  #  if (!value)
   #     return(object)

  #  if (length(object@env$mz)==0) {
  #      warning("MS1 scans empty. Skipping profile matrix calculation.")
   #     return(object)
  #  }

  #  mzrange <- range(object@env$mz)

    ## calculate the profile matrix with whole-number limits
   # minmass <- floor(mzrange[1])
   # maxmass <- ceiling(mzrange[2])

   # num <- (maxmass - minmass)/value + 1
   # profFun <- match.profFun(object)
   # object@env$profile <- profFun(object@env$mz, object@env$intensity,
   #                               object@scanindex, num, minmass, maxmass,
   #                               FALSE, object@profparam)

   # object@mzrange <- c(minmass, maxmass)
   # return(object)
#})


