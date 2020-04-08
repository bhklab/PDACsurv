#' Define an S4 Generic for the methods::as function
#'
#' This will allow creation of new definitons for object conversions
#'
#' @export
setGeneric('as', function(object, value) methods::as(object, value, ...))

#' Convert factor vector to numeric vector
#' 
#' @param f A \code{factor} vector with numeric levels
#' 
#' @return A \code{numeric} vector of the factor levels in the same order
#' 
#' @export
setMethod('as',
          signature('factor'),
          function(object, value) {
            switch(value,
                   'numeric'={as.numeric(levels(object))[object]},
                   'character'={as.character(object)},
                   'list'={as.list(object) },
                   stop(paste0('Cannot coerce to class ', value))
            )
          })