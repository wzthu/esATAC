################################################################################
##                                                                            ##
##                                                                            ##
##            Functions related with qQCReport are from QuasR.                ##
##     We integrate these function in order to support the pipeline           ##
##     in our package. We recommend using qQCReport in QuauR if you           ##
##          want do quality control for your sequencing solely.               ##
##                                                                            ##
##                                                                            ##
################################################################################
#' @importFrom ShortRead FastqSampler
#' @importFrom ShortRead yield
#' @importFrom ShortRead readFasta
#' @importFrom ShortRead sread
#' @importFrom ShortRead alphabetByCycle
#' @importFrom tools file_ext
#' @importFrom IRanges width
#' @importFrom graphics strwidth
#' @importFrom BiocGenerics which
#' @importFrom S4Vectors runValue
#' @importFrom S4Vectors runLength
#' @importFrom Biostrings alphabetFrequency
#' @importFrom graphics layout
#' @importFrom ShortRead tables
#' @importFrom graphics par
#' @importFrom graphics box
#' @importFrom graphics strheight
#' @importFrom graphics text
#' @importFrom graphics matplot
#' @importFrom graphics rect
#' @importFrom grDevices dev.new
qQCReport <- function(input, pdfFilename=NULL, chunkSize=1e6L, useSampleNames=FALSE, ...){
    if(is.character(input)){
        filetype <- unique(consolidateFileExtensions(input, compressed=TRUE))
        if(length(filetype) > 1L)
            stop("parameter 'input' must consist of unique filetype")
        if(useSampleNames && !is.null(names(input))) {
            label <- sprintf("%i. %s", 1:length(input), names(input))
        } else {
            label <- sprintf("%i. %s", 1:length(input), basename(input))
        }
        mapLabel <- label
        if(all(file.exists(input))){
            if(filetype != "bam"){
                readFilename <- input
                alnFilename <- NULL
                genome <- NULL
            } else {
                readFilename <- input
                alnFilename <- input
                genome <- NULL
            }
        } else {
            stop("could not find the files '", paste(input[!file.exists(input)], collapse=" "), "'")
        }
    } else {
        stop("'input' must be filenames")
    }

    qc1L <- mapply(FUN = calcQaInformation, filename = as.list(readFilename),
                   label = as.list(label), MoreArgs = list(filetype=filetype, chunkSize=chunkSize))

    qa <- do.call(rbind, qc1L)

    # open/close pdf stream
    if(!is.null(pdfFilename)){
        pdf(pdfFilename, paper="default", onefile=TRUE, width=0, height=0)
        on.exit(dev.off())
    }

    unique <- NULL
    mm <- NULL
    frag <- NULL

    message("creating QC plots")
    plotdata <- list(raw=list(qa=qa, mm=mm, frag=frag, unique=unique, mapdata=NULL))
    if(!is.null(qa)){
        if(filetype == "fastq" || filetype == "bam"){
            if(is.null(pdfFilename))
                dev.new()
            plotdata[['qualByCycle']] <- plotQualByCycle(qa)
        }

        if(is.null(pdfFilename))
            dev.new()
        plotdata[['nuclByCycle']] <- plotNuclByCycle(qa)

        if(is.null(pdfFilename))
            dev.new()
        plotdata[['duplicated']] <- plotDuplicated(qa)

    } else if(!is.null(mm)){
        if(is.null(pdfFilename))
            dev.new()
        plotdata[['nuclByCycle']] <- plotNuclByCycle(mm)
    }

    invisible(plotdata)
}

consolidateFileExtensions <- function(filename, compressed=FALSE) {
    # convert various sequence file extension to one representative
    if(compressed)
        fileExtension <- tolower(tools::file_ext(sub("[.](gz|bz2|xz)$","",filename)))
    else
        fileExtension <- tolower(tools::file_ext(filename))
    fileExtension[ fileExtension %in% c("fa", "fna", "fasta") ] <- "fasta"
    fileExtension[ fileExtension %in% c("fq", "fastq") ] <- "fastq"
    return(fileExtension)
}

compressedFileFormat <- function(filename) {
    ifelse(grepl("[.](gz|bz2|xz)$",filename),
           c("gz"="gzip", "bz2"="bzip2", "xz"="xz")[tools::file_ext(filename)],
           "none")
}

calcQaInformation <- function(filename, label, filetype, chunkSize){
    if(any(filetype == "fasta") && any(compressedFileFormat(filename) != "none")){
        qa <- NULL
        warning("compressed 'fasta' input is not yet supported")
    }else{
        reads <- switch(as.character(filetype),
                        fastq = {
                            f <- ShortRead::FastqSampler(filename, n=chunkSize)
                            reads <- ShortRead::yield(f)
                            close(f)
                            reads[IRanges::width(reads)>0]
                        },
                        fasta = {
                            reads <- ShortRead::readFasta(as.character(filename), nrec=chunkSize)
                            reads[IRanges::width(reads)>0]
                        }
                        # bam are not in consideration
                        # bam = {
                        #     bf <- Rsamtools::BamFile(filename, yieldSize=1e6)
                        #     myyield <- function(x) {
                        #         tmp <- Rsamtools::scanBam(x, param=ScanBamParam(what=c("seq", "qual", "strand")))[[1]]
                        #         minusStrand <- !is.na(tmp$strand) & tmp$strand == "-"
                        #         ShortRead::ShortReadQ(sread=c(tmp$seq[!minusStrand],Biostrings::reverseComplement(tmp$seq[minusStrand])),
                        #                               quality=c(tmp$qual[!minusStrand],Biostrings::reverse(tmp$qual[minusStrand])))
                        #     }
                        #     reads <- reduceByYield(X=bf, YIELD=myyield, MAP=identity,
                        #                            REDUCE=REDUCEsampler(sampleSize=chunkSize, verbose=FALSE),
                        #                            parallel=FALSE)
                        #     reads[width(reads)>0]
                        # }
        )
    }
    return(qa(reads, label))
}

truncStringToPlotWidth <- function(s, plotwidth) {
    sw <- graphics::strwidth(s)
    if(any(sw > plotwidth)) {
        w <- 10 # number of character to replace with ".."
        l <- nchar(s)
        news <- s
        i <- sw > plotwidth
        while(w < l && any(i)) {
            news <- ifelse(i, paste(substr(s, 1, ceiling((l-w)/2)+5), substr(s, floor((l+w)/2)+5, l), sep="..."), news)
            sw <- graphics::strwidth(news)
            i <- sw > plotwidth
            w <- w + 2
        }
        return(news)
    } else {
        return(s)
    }
}

plotQualByCycle <- function(qcdata, lmat=matrix(1:12, nrow=6, byrow=TRUE)) {
    data <- qcdata[['perCycle']][['quality']]
    qtiles <- by(list(data$Score, data$Count), list(data$lane, data$Cycle), function(x) {
        coef <- 1.5
        scoreRle <- Rle(x[[1]], x[[2]])
        n <- length(scoreRle)
        nna <- !is.na(scoreRle)
        stats <- c(min(scoreRle), quantile(scoreRle, c(0.25, 0.5, 0.75)), max(scoreRle))
        iqr <- diff(stats[c(2,4)])
        out <- if (!is.na(iqr)) {
            scoreRle < (stats[2L] - coef * iqr) | scoreRle > (stats[4L] + coef * iqr)
        } else !is.finite(scoreRle)
        if (any(out[nna], na.rm = TRUE))
            stats[c(1, 5)] <- range(scoreRle[!out], na.rm = TRUE)
        conf <- stats[3L] + c(-1.58, 1.58) * iqr/sqrt(n)
        list(stats=stats, n=n, conf=conf, out=BiocGenerics::which(out))
    }, simplify=FALSE)
    ns <- nrow(qtiles)

    qtilesL <- lapply(1:ns, function(i) {
        tmpconf <- do.call(cbind,lapply(qtiles[i,], "[[", 'conf'))
        tmpxn <- ncol(tmpconf)
        list(stats=do.call(cbind,lapply(qtiles[i,][1:tmpxn], "[[", 'stats')),
             n=sapply(qtiles[i,][1:tmpxn], "[[", 'n'),
             conf=tmpconf,
             out=numeric(0), #sapply(qtiles[i,][1:tmpxn], "[[", 'out'),
             group=numeric(0), #rep(1:ncol(qtiles), sapply(qtiles[i,], function(x) length(x$out))),
             names=colnames(tmpconf))
    })
    names(qtilesL) <- rownames(qtiles)

    graphics::layout(lmat)
    for(i in 1:ns) {
        xn <- length(qtilesL[[i]]$names)
        ym <- max(35,max(qtilesL[[i]]$stats))
        par(mar=c(5-1,4-1,4-4,2-1)+.1, mgp=c(3-1,1-0.25,0))
        plot(0:1,0:1,type="n", xlab="Position in read (bp)", ylab="Quality score", xlim=c(0,xn)+0.5, xaxs="i", ylim=c(0,ym))
        rect(xleft=seq.int(xn)-0.5, ybottom=-10, xright=seq.int(xn)+0.5, ytop=20,    col=c("#e6afaf","#e6c3c3"), border=NA)
        rect(xleft=seq.int(xn)-0.5, ybottom=20,  xright=seq.int(xn)+0.5, ytop=28,    col=c("#e6d7af","#e6dcc3"), border=NA)
        rect(xleft=seq.int(xn)-0.5, ybottom=28,  xright=seq.int(xn)+0.5, ytop=ym+10, col=c("#afe6af","#c3e6c3"), border=NA)
        do.call("bxp", c(list(qtilesL[[i]], notch=FALSE, width=NULL, varwidth=FALSE, log="", border=par('fg'),
                              pars=list(boxwex=0.8, staplewex=0.5,  outwex=0.5, boxfill="#99999944"),
                              outline=FALSE, horizontal=FALSE, add=TRUE, at=1:xn, axes=FALSE)))
        cxy <- par('cxy')
        text(x=par('usr')[1]+cxy[1]/4, y=par('usr')[3]+cxy[2]/4, adj=c(0,0),
             labels=truncStringToPlotWidth(rownames(qtiles)[i], diff(par("usr")[1:2]) - cxy[1]/2))
        box()
    }

    invisible(qtilesL)
}

plotNuclByCycle <- function(qcdata, lmat=matrix(1:12, nrow=6, byrow=TRUE)) {
    if(!is.null(qcdata[['perCycle']])){
        data <- qcdata[['perCycle']][['baseCall']]
        nfreq <- by(list(data$Base, data$Count), list(data$lane, data$Cycle), function(x) {
            y <- x[[2]] /sum(x[[2]]) *100
            names(y) <- x[[1]]
            y
        }, simplify=TRUE)
        ns <- nrow(nfreq)

        nfreqL <- lapply(1:ns, function(i) {
            tmp <- do.call(cbind, lapply(nfreq[i,], function(x) x[c('A','C','G','T','N')]))
            tmp[is.na(tmp)] <- 0
            rownames(tmp) <- c('A','C','G','T','N')
            tmp
        })
        names(nfreqL) <- rownames(nfreq)
    } else {
        nfreqL <- lapply(qcdata, function(x){
            nfreq <- do.call(rbind, lapply(1:dim(x)[3], function(j){
                c(A=sum(x[,"A",j]),
                  C=sum(x[,"C",j]),
                  G=sum(x[,"G",j]),
                  T=sum(x[,"T",j]),
                  N=sum(x[,"N",j]))
            }))
            rownames(nfreq) <- 1:dim(x)[3]
            t(nfreq / rowSums(nfreq) * 100)
        })
        ns <- length(nfreqL)
        #         names(nfreqL) <- names(qcdata)
    }
    nfreqL[is.na.data.frame(nfreqL)] <- 0

    mycols <- c("#5050ff","#e00000","#00c000","#e6e600","darkgray")
    graphics::layout(lmat)
    for(i in 1:ns) {
        xn <- ncol(nfreqL[[i]])
        ym <- max(50,ceiling(max(nfreqL[[i]], na.rm = TRUE) /5) *5 +5)
        par(mar=c(5-1,4-1,4-4,2-1)+.1, mgp=c(3-1,1-0.25,0))
        matplot(1:xn, t(nfreqL[[i]]), type="o", xlab="Position in read (bp)", ylab="Nucleotide frequency (%)",
                xlim=c(0,xn)+0.5, xaxs="i", ylim=c(0,ym), lwd=2, lty=1, pch=20, cex=0.6, col=mycols)
        abline(h=0,lty=2,col="gray")
        cxy <- par('cxy')
        text(x=par('usr')[1]+cxy[1]/4, y=par('usr')[4]-cxy[2]/4, adj=c(0,1),
             labels=truncStringToPlotWidth(names(nfreqL)[i], diff(par('usr')[1:2])-1.5*cxy[1]-strwidth(paste(rownames(nfreqL[[i]]), collapse=" "))))
        text(x=par('usr')[2]-cxy[1]/4-(4:0)*cxy[1]*0.8, y=par('usr')[4]-cxy[2]/4, adj=c(1,1), col=mycols, labels=rownames(nfreqL[[i]]))
    }

    invisible(nfreqL)
}

plotDuplicated <- function(qcdata, breaks=c(1:10), lmat=matrix(1:6, nrow=3, byrow=TRUE)) {
    if(breaks[length(breaks)]<Inf)
        breaks <- c(breaks,breaks[length(breaks)]+1,Inf)
    breakNames <- c(as.character(breaks[1:(length(breaks)-2)]),paste(">",breaks[length(breaks)-2],sep=""))
    data <- qcdata[['sequenceDistribution']]
    nocc <- by(list(data$nOccurrences, data$nReads), list(data$lane),
               function(x) Rle(values=x[[1]], lengths=x[[2]]), simplify=FALSE)
    ns <- nrow(nocc)

    nbin <- length(breaks)-1
    bocc <- do.call(rbind, lapply(1:ns, function(i) {
        bin <- findInterval(S4Vectors::runValue(nocc[[i]]), breaks)
        tmp <- tapply(S4Vectors::runLength(nocc[[i]]),bin,sum) /length(nocc[[i]]) *100
        tmp2 <- numeric(nbin)
        tmp2[ as.integer(names(tmp)) ] <- tmp
        tmp2
    }))
    dimnames(bocc) <- list(sampleName=rownames(nocc), duplicationLevel=breakNames)

    graphics::layout(lmat)
    for(i in 1:ns) {
        nm <- rownames(nocc)[i]
        fs <- qcdata[['frequentSequences']]
        fs <- fs[ fs$lane==nm, ]
        xn <- length(breaks)-1
        ym <- max(50,ceiling(max(bocc[i,]) /5) *5)
        graphics::par(mar=c(5-1,4-1,4-4,2-1)+.1, mgp=c(3-1,1-0.25,0))
        plot(1:xn, bocc[i,], type="o", xlab="Sequence duplication level", ylab="Percent of unique sequences",
             xlim=c(0,xn)+0.5, xaxs="i", ylim=c(0,ym), lwd=2, lty=1, pch=20, cex=0.6, axes=FALSE,
             panel.first=abline(h=0, lty=2, col="gray"))
        axis(1, at=1:xn, labels=breakNames)
        axis(2)
        box()
        frqseqcex <- 0.8
        frqseqS <- sprintf("%-*s",max(nchar(as.character(fs[,'sequence']))),fs[,'sequence'])
        frqseqF <- sprintf("(%6.0f)",fs[,'count']/sum(runValue(nocc[[i]])*as.numeric(runLength(nocc[[i]])))*1e6)
        frqseqJ <- " "
        frqseqW <- max(nchar(as.character(fs[,'sequence'])))
        xleft <- par('usr')[2] - max(strwidth(paste(frqseqS," ",frqseqF,sep=""), cex=frqseqcex, family="mono"))
        while(xleft < 1 && frqseqW > 7) {
            frqseqJ <- ".. "
            frqseqW <- frqseqW - 2
            frqseqS <- strtrim(frqseqS,frqseqW)
            xleft <- par('usr')[2] - max(strwidth(paste(frqseqS,frqseqJ,frqseqF,sep=""), cex=frqseqcex, family="mono"))
        }
        if(xleft >= 1 && frqseqW > 5) {
            cxy <- par('cxy')
            ytop <- par('usr')[4] - 2.0*cxy[2]
            yoff <- ytop - 1.8*cumsum(strheight(frqseqS, cex=frqseqcex, family="mono"))
            ii <- yoff+diff(yoff[1:2]) > max(bocc[i, 1:xn > floor(xleft)])
            if(any(ii)) {
                text(x=xleft, y=ytop,     adj=c(0,0),
                     labels=paste(truncStringToPlotWidth(nm,par("usr")[2]-xleft-cxy[1]/2),"frequent sequences (per Mio.):",sep="\n"))
                text(x=xleft, y=yoff[ii], adj=c(0,1),
                     labels=paste(frqseqS,frqseqJ,frqseqF,sep="")[ii], family="mono", cex=frqseqcex)
            } else {
                text(x=xleft, y=ytop,     adj=c(0,0),
                     labels=truncStringToPlotWidth(nm,par("usr")[2]-xleft-cxy[1]/2))
            }
        }
    }

    invisible(bocc)
}

.qa_ShortRead <- function(dirPath, lane, ..., verbose=FALSE) {
    if (missing(lane))
        stop("Paramteter 'lane' is missing.")
    obj <- dirPath
    alf <- Biostrings::alphabetFrequency(ShortRead::sread(obj), baseOnly=TRUE, collapse=TRUE)
    #     bqtbl <- alphabetFrequency(quality(obj), collapse=TRUE)
    #     rqs <- .qa_qdensity(quality(obj))
    freqtbl <- ShortRead::tables(ShortRead::sread(obj))
    abc <- ShortRead::alphabetByCycle(obj)
    names(dimnames(abc)) <- c("base", "cycle")
    dimnames(abc)$cycle <- as.character(1:dim(abc)[2])
    ac <- ShortRead:::.qa_adapterContamination(obj, lane, ...)
    perCycleBaseCall <- data.frame(Cycle = as.integer(colnames(abc)[col(abc)]),
                                   Base = factor(rownames(abc)[row(abc)]), Count = as.vector(abc),
                                   lane = lane, row.names = NULL)
    perCycleBaseCall <- perCycleBaseCall[perCycleBaseCall$Count != 0, ]
    #     perCycleBaseCall <- ShortRead:::.qa_perCycleBaseCall(abc, lane)
    #     perCycleQuality <- .qa_perCycleQuality(abc, quality(obj), lane)
    lst <- list(readCounts=data.frame(
        read=length(obj), filter=NA, aligned=NA,
        row.names=lane),
        baseCalls=data.frame(
            A=alf[["A"]], C=alf[["C"]], G=alf[["G"]], T=alf[["T"]],
            N=alf[["other"]], row.names=lane),
        readQualityScore=data.frame(
            quality=NULL, #rqs$x,
            density=NULL, #rqs$y,
            lane=NULL, #lane,
            type=NULL #"read"
        ),
        baseQuality=data.frame(
            score=NULL, #names(bqtbl),
            count=NULL, #as.vector(bqtbl),
            lane=NULL #lane
        ),
        alignQuality=data.frame(
            score=as.numeric(NA),
            count=as.numeric(NA),
            lane=lane, row.names=NULL),
        frequentSequences=data.frame(
            sequence=names(freqtbl$top),
            count=as.integer(freqtbl$top),
            type="read",
            lane=lane),
        sequenceDistribution=cbind(
            freqtbl$distribution,
            type="read",
            lane=lane),
        perCycle=list(
            baseCall=perCycleBaseCall,
            quality=NULL #perCycleQuality
        ),
        perTile=list(
            readCounts=data.frame(
                count=integer(0), type=character(0),
                tile=integer(0), lane=character(0)),
            medianReadQualityScore=data.frame(
                score=integer(), type=character(), tile=integer(),
                lane=integer(), row.names=NULL)),
        adapterContamination=ac

    )

    ShortRead:::.ShortReadQQA(lst)
}

setMethod(qa, "ShortRead", .qa_ShortRead)



############################get PWM matrix from user's file#####################
#' @name getMotifInfo
#' @title Generate PFMatrix or PFMatrixList from file.
#' @description
#' atacMotifScan and atacMotifScanPair accept PFM in a \code{list}, this
#' function convert JASPAR PFM file to \code{\link{PFMatrix}} or \code{\link{PFMatrixList}}.
#' @param motif.file Motif PFM file downloaded from JASPAR.
#' @details Generate \code{\link{PFMatrix}} or \code{\link{PFMatrixList}}.
#' @return \code{\link{PFMatrix}} or \code{\link{PFMatrixList}}.
#' @author Wei Zhang
#' @importFrom TFBSTools PFMatrix
#' @importFrom TFBSTools as.matrix
#' @importFrom TFBSTools PFMatrixList
#' @examples
#'
#' motif_file <- system.file("extdata", "CustomizedMotif.txt", package="esATAC")
#' pfm <- getMotifInfo(motif.file = motif_file)
#'
#' @export
getMotifInfo <- function(motif.file = NULL){
    PWMList <- list()
    name.flag <- FALSE
    A.flag <- FALSE
    C.flag <- FALSE
    G.flag <- FALSE
    T.flag <- FALSE

    con <- file(motif.file, "r")
    line <- readLines(con, n = 1)
    while( length(line) != 0 ){

        if(substring(line, 1, 1) == ">"){
            name_str <- unlist(strsplit(x = sub(pattern = ">", replacement = "", x = line), split = "\\s+"))
            motif_name <- tail(name_str, n = 1)
            name.flag <- TRUE
        }else if(substring(line, 1, 1) == "A"){
            A_str_num <- regmatches(line, gregexpr("[[:digit:]]+", line))
            A_num <- as.numeric(unlist(A_str_num))
            A.flag <- TRUE
        }else if(substring(line, 1, 1) == "C"){
            C_str_num <- regmatches(line, gregexpr("[[:digit:]]+", line))
            C_num <- as.numeric(unlist(C_str_num))
            C.flag <- TRUE
        }else if(substring(line, 1, 1) == "G"){
            G_str_num <- regmatches(line, gregexpr("[[:digit:]]+", line))
            G_num <- as.numeric(unlist(G_str_num))
            G.flag <- TRUE
        }else if(substring(line, 1, 1) == "T"){
            T_str_num <- regmatches(line, gregexpr("[[:digit:]]+", line))
            T_num <- as.numeric(unlist(T_str_num))
            T.flag <- TRUE
        }

        if(name.flag & A.flag & C.flag & G.flag & T.flag){
            p_matrix <- matrix(data = c(A_num, C_num, G_num, T_num), nrow = 4,
                               byrow = TRUE,  dimnames=list(c("A", "C", "G", "T")))

            p_matrix <- TFBSTools::PFMatrix(profileMatrix = p_matrix)

            PWMList[[motif_name]] <- p_matrix

            name.flag <- FALSE
            A.flag <- FALSE
            C.flag <- FALSE
            G.flag <- FALSE
            T.flag <- FALSE
        }

        line <- readLines(con, n = 1)
    }
    close(con)

    PWMList <- do.call(PFMatrixList, PWMList)

    return(PWMList)
}

checkParam <- function(paramlist,paramPattern,...){
    rs<-grepl(paramPattern, paramlist)
    if(sum(rs)>0){
        banp=paste(paramlist[rs], collapse = "'/'")
        stop(sprintf("Parameter(s) '%s' are not acceptable in paramList. it should be set as fix parameter.",banp))
    }
}

