CutSiteCountR <- R6::R6Class(
  classname = "CutSiteCountR",
  inherit = BaseProc,

  public = list(
    initialize = function(atacProc, CutSiteFile = NULL, MotifFile = NULL, MatrixPath = NULL,
                          motif_length = NULL, strand_length = NULL, FootPrint = FALSE, editable = FALSE){
      super$initialize("CutSiteCountR",editable,list(arg1=atacProc))
      print("CutSiteCountRInitCall")

      # necessary and unchanged parameters, this should be tested
      if(!is.null(atacProc)){
        print("Parameter atacProc is not using now! We will add more functions in the future!")
      }
      private$paramlist[["CutSiteFile"]] <- CutSiteFile
      private$paramlist[["MotifFile"]] <- MotifFile
      private$paramlist[["MatrixPath"]] <- MatrixPath
      private$paramlist[["motif_length"]] <- motif_length
      private$paramlist[["strand_length"]] <- strand_length
      private$paramlist[["FootPrint"]] <- FootPrint
      # check parameter
      private$checkRequireParam()
      private$checkPathExist(private$paramlist[["CutSiteFile"]])
      private$checkPathExist(private$paramlist[["MotifFile"]])
      private$checkPathExist(private$paramlist[["MatrixPath"]])
      print("finishCutSiteCountRInitCall")
    }, # initialization end

    processing = function(){
      if(!super$processing()){
        return()
      }
      chr <- as.list(1:22)
      chr[[23]] <- "X"
      chr[[24]] <- "Y"
      for(i in seq(1:24)){
        CutSiteInput <- paste0(private$paramlist[["CutSiteFile"]], "_chr", chr[[i]], ".bed", collapse = "")
        MotifInput <- paste0(private$paramlist[["MotifFile"]], "_chr", chr[[i]], ".bed", collapse = "")
        MatrixOutput <- paste0(private$paramlist[["MatrixPath"]], "_chr", chr[[i]], ".matrix", collapse = "")
        .CutSiteCount(readsfile = CutSiteInput, motiffile = MotifInput, matrixfile = MatrixOutput,
                      motif_len = private$paramlist[["motif_length"]], strand_len = private$paramlist[["strand_length"]])
        if(i == 1){
          data <- read.table(MatrixOutput)
        }else{
          temp <- read.table(MatrixOutput)
          data <- rbind(data, temp)
        }
      }
      OverallMatrix <- paste0(private$paramlist[["MatrixPath"]], "data", ".matrix", collapse = "")
      write.table(data, file = OverallMatrix, row.names = FALSE, col.names = FALSE, quote = FALSE)
      if(private$paramlist[["FootPrint"]]){
        fp <- apply(data, 2, sum)
        plot(fp, type = "l", col = "blue", lwd = 2, xlab = "Relative Distance From Motif (bp)", ylab = "Cut Site Count", xaxt = "n", yaxt = "n")
        axis(1, at = seq(1, private$paramlist[["strand_length"]], len = 3),
              labels = -(private$paramlist[["strand_length"]] + 1 - seq(1, private$paramlist[["strand_length"]] + 1, len = 3)),
              padj = -1.0, tck = -0.01)
        axis(1, at = private$paramlist[["strand_length"]] + private$paramlist[["motif_length"]] + seq(1, private$paramlist[["strand_length"]], len = 3),
             labels = seq(0, private$paramlist[["strand_length"]], len = 3),
             padj = -1.0, tck = -0.01)
        axis(2, padj = 1.0,tck = -0.02)
        abline(v = c(private$paramlist[["strand_length"]], private$paramlist[["strand_length"]] + private$paramlist[["motif_length"]] + 1),
               lty = 2)
      }
      private$setFinish()
    }, # processing end

    setResultParam = function(){

    }  # setResultParam end

  ), # public end


  private = list(
    checkRequireParam = function(){
      if(private$editable){
        return();
      }
      if(is.null(private$paramlist[["CutSiteFile"]])){
        stop("Parameter CutSiteFile is required!")
      }
      if(is.null(private$paramlist[["MotifFile"]])){
        stop("Parameter MotifFile is required!")
      }
      if(is.null(private$paramlist[["MatrixPath"]])){
        stop("Parameter MatrixPath is required!")
      }
      if(is.null(private$paramlist[["motif_length"]])){
        stop("Parameter motif_length is required!")
      }
      if(is.null(private$paramlist[["strand_length"]])){
        stop("Parameter strand_length is required!")
      }
    }
  ) # private end

) # class end
