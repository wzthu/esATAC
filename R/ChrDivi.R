ChrDivi <- R6::R6Class(
  classname = "ChrDivi",
  inherit = BaseProc,
  public = list(
    initialize = function(atacProc, ReadsIfile = NULL, ReadsOpath = NULL, name = NULL,
                          sort = NULL, sort_col = NULL, editable = FALSE){
      super$initialize("ChrDivi",editable,list(arg1=atacProc))
      print("ChrDiviInitCall")
      # add "/" in the output path
      if(substr(ReadsOpath, nchar(ReadsOpath), nchar(ReadsOpath)) != "/"){
        ReadsOpath <- paste0(ReadsOpath, "/", collapse = "")
      }
      # necessary and unchanged parameters, this should be tested
      private$paramlist[["sort"]] <- sort
      private$paramlist[["sort_col"]] <- sort_col
      private$paramlist[["ReadsIfile"]] <- ReadsIfile
      private$paramlist[["ReadsOpath"]] <- ReadsOpath
      if(is.null(name)){
        private$paramlist[["name"]] <- "Output"
      }else{
        private$paramlist[["name"]] <- name
      }
      # parameter check
      private$checkRequireParam()
      private$checkFileExist(private$paramlist[["ReadsIfile"]])
      if(!dir.exists(private$paramlist[["ReadsOpath"]])){
        stop(paste("error, path does not exist:",private$paramlist[["ReadsOpath"]]))
      }
      print("finishChrDiviInitCall")
    },

    processing = function(){
      if(!super$processing()){
        return()
      }
      .chr_separate_call(ReadsIfile = private$paramlist[["ReadsIfile"]], ReadsOpath = private$paramlist[["ReadsOpath"]],
                         Name = private$paramlist[["name"]])
      if(!is.null(private$paramlist[["sort"]]) && (private$paramlist[["sort"]])){
        chr <- as.list(1:22)
        chr[[23]] <- "X"
        chr[[24]] <- "Y"
        for(i in seq(1:24)){
          sfile <- paste(private$paramlist[["ReadsOpath"]], private$paramlist[["name"]],
                        "_chr", chr[[i]], ".bed", sep = "", collapse = "")
          data <- read.table(sfile)
          data <- data[order(data[, private$paramlist[["sort_col"]]]), ]
          write.table(data, file = sfile, row.names = FALSE, col.names = FALSE, quote = FALSE)
        }
      }
      private$setFinish()
    },

    setResultParam = function(){
      print("This function is not using!")
    }

  ), #public end

  private = list(
    checkRequireParam = function(){
      if(private$editable){
        return();
      }
      if(is.null(private$paramlist[["ReadsIfile"]])){
        stop("Parameter ReadsIfile is required!")
      }
      if(is.null(private$paramlist[["ReadsOpath"]])){
        stop("Parameter ReadsOpath is required!")
      }
      if(!is.null(private$paramlist[["sort"]]) && (private$paramlist[["sort"]])){
        if(is.null(private$paramlist[["sort_col"]])){
          stop("Please specify which column you want sort!")
        }
      }
    }
  ) # private end

)
