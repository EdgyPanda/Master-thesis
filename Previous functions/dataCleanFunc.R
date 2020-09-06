#
# 
#-------------------------------DATA CLEANING------------------------------------------
# DESCRIPTION: data.csv contains Nan-values in each columns that needs to be
# removed.  This implies columns of different row-lengths and had to be saved in
# a list.  We need to construct a timestamp of Date-Time object, such that we can
# work with xts matrices in R. This implies better functionality since R
# recognizes xts files as multivariate time-series objects.



#' nanRemover
#'
#' Takes data and removes the  missing data (Nan-values) in each column in each list component in data.
#' @param data of type matrix.
#' @return data of type list with removed Nan-values for each day.
#' @import xts stats
#' @export
nanRemover <- function(data) {

    datatest <- list()

    counter <- seq(1, length(data), 0.5)

    for (i in seq(1, length(data), 2)) {

        datatest[[counter[i]]] <- data[!is.na(data[seq(2, length(data), 2)][counter[i]]),
            i:(i + 1)]

    }

    return(datatest)

}


#' getDates
#'
#' Isolates dates from the header, when the first header-value is a date.
#' @param data of type matrix.
#' @return list with all of the dates in data.
#' @export
getDates <- function(data) {

    dates <- colnames(data)[seq(1, length(data), 2)]
    dates <- gsub(".*_", "", dates)

    return(dates)
}


#' datesToPOSIXct
#' 
#' Turns daily dates into POSIXct dates used in xts objects.
#' @param dates of type character.
#' @return dates as POSIXct elements.
#' @export 
datesToPOSIXct <- function(dates) {

    for (i in 1:length(dates)){

        dates[i] <- gsub('^([0-9]{4})([0-9]+)$', '\\1-\\2', dates[i])
        dates[i] <- gsub('([[:digit:]]{2,2})$', '-\\1\\2', dates[i])

    }

    dates <- as.POSIXct((dates))

    return(dates)
}


#' timeStamp
#'
#' Create usable timestamps for xts object.
#' @param cleanedData A List object with removed NAN-values from the nanRemover-function
#' @param Dates A List object with the scraped dates following from the getDates function.
#' @return list with all timestamps. List components are days and list-elements are intraday timestamps.
#' @export
timestamp <- function(cleanedData, Dates) {

    timestamps <- list()

    for (i in 1:length(cleanedData)) {

        timestamps[[i]] <- as.character(unlist(cleanedData[[i]][1:dim(cleanedData[[i]])[1],
            1]))
        timestamps[[i]] <- as.POSIXct(paste(Dates[i], timestamps[[i]]), format = "%Y%m%d %H:%M:%S",
            tz = "")
        print(sprintf("current iteration: %d out of: %d", i, length(cleanedData)))
    }

    return(timestamps)
}


#' timeseriesList
#'
#' Create list of xts objects for all days in data.
#' @param data The original data used to get the columnnames of the timeseriesList.
#' @param cleanedData A List object with removed NAN-values from the nanRemover-function.
#' @param timestamps A List object with the timestamps from the timestamp function.
#' @return list with all intraday data using the timestamps created from the timestamp function as index.
#' @export
timeseriesList <- function(data, cleanedData, timestamps) {

    timeseriesList <- list()

    for (i in 1:length(cleanedData)) {

        timeseriesList[[i]] <- xts(cleanedData[[i]][1:dim(cleanedData[[i]])[1], 2],
            order.by = timestamps[[i]])
        colnames(timeseriesList[[i]]) <- gsub("\\_.*", "", colnames(data)[1])
    }

    return(timeseriesList)

}


#' mergedTimeseries
#'
#' Merges two timeseries objects from the function timeseriesList, and removes the removes the rows where each timeseries contains missing data.
#' @param timeseries1 A List object containing XTS objects as list elements.
#' @param timeseries2 A List object containing XTS objects as list elements..
#' @return list object containing the merged timeseries objects.
#' @export
mergedTimeseries <- function(timeseries1, timeseries2) {

    mergedList <- list()

    if ((length(timeseries1) != length(timeseries2))) {

        stop("List ojects not the same length")

    }

    for (i in 1:length(timeseries1)) {

        mergedList[[i]] <- na.omit(merge.xts(timeseries1[[i]], timeseries2[[i]]))
    }

    return(mergedList)

}


