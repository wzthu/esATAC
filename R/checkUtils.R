checkFileExist <- function(filePath,argname){
    if(!is.null(filePath)){
        for(i in 1:length(filePath)){
            if(!file.exists(filePath[i])){
                stop(sprintf("For argument `%s`, file does not exist: `%s`",
                             argname,filePath[i]))
            }
        }
    }
}
checkPathExist <- function(filePath,argname){
    if(!is.null(filePath)){
        if(!dir.exists(dirname(filePath))){
            stop(sprintf("For argument `%s`,path does not exist: `%s`",
                         argname,filePath))
        }
    }
}
checkFileCreatable <- function(filePath,argname,overwrite){
    if(!is.null(filePath)){
        if(file.exists(filePath)){
            if(overwrite){
                warning(sprintf(paste0("For argument `%s`, file exist:%s. ",
                                       "It will be overwrited"),
                                       argname,filePath));
            }else{
                stop(sprintf(paste0("For argument `%s`,file exist: %s. ",
                                    "Use 'overwrite=TRUE' to overwrite"),
                             argname,filePath));
            }
        }else if(!file.create(filePath)){
            stop(sprintf(paste0("For argument `%s`, cannot create file `%s`.",
                                "\nNo such directory or permission denied."),
                         argname,filePath));
            stop("")
        }else{
            unlink(filePath)
        }
    }
}

checkAddArgus<- function(pattern,...){
    paramlist<-trimws(as.character(list(...)))
    paramArray<-c()
    if(length(paramlist)>0){
        for(i in 1:length(paramlist)){
            paramArray<-c(paramArray,strsplit(paramlist[i],"\\s+")[[1]])
        }
    }
    fixed<-grepl(pattern,paramArray)
    if(sum(fixed)>0){
        invalidp<-paste0(paramArray[fixed],collapse = " ")
        stop(sprintf(paste0("Argument(s) `%s` are invalid for additional ",
                            "argument. Please set them as fixed arguments."),
                     invalidp))
    }
    invisible(paramArray)
}


