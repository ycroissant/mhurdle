#' @importFrom generics glance
#' @export
generics::glance

#' @importFrom generics tidy
#' @export
generics::tidy


#' @importFrom prediction prediction
#' @export
prediction::prediction


#' @importFrom margins margins
#' @export
margins::margins


#' broom's methods
#'
#' Methods to compute extract in a tidy way the elements of a fitted
#' model
#'
#' @name broom
#' @param x a model fitted with mhurdle
#' @param conf.int,conf.level current see `generics::tidy` (currently
#'     unused)
#' @param ... further arguments, currently unused
#' @details `mhurdle` exports the `generics::tidy` and
#'     `generics::glance` functions. The specific method provided for
#'     `mhurdle` objects enables the use of some package that relies
#'     on these functions (`modelsummary` for example)
NULL

#' @rdname broom
#' @method tidy mhurdle
#' @export
tidy.mhurdle <- function(x, conf.int = FALSE, conf.level = 0.95, ...){
    result <- summary(x)$coefficients
    nms <- rownames(result)
    rownames(result) <- NULL
    result <- data.frame(term = nms, result)
    names(result) <- c("term", "estimate", "std.error", "statistic", "p.value")
    result
}

#' @rdname broom
#' @method glance mhurdle
#' @export
glance.mhurdle <- function(x, ...){
    N <- nobs(x)
    N.zero <- nobs(x, "null")
    N.pos <- nobs(x, "pos")
    logLik = logLik(x)
    R2.zero = x$R2[1]
    R2.pos = x$R2[2]
    dpar = x$dpar
    data.frame(nobs = N, dpar = dpar, nobs.zero = N.zero, nobs.pos = N.pos, logLik = logLik, R2.zero = R2.zero, R2.pos = R2.pos)
}

# copied from the prediction package
# internal function that overrides the defaults of `data.frame()`
make_data_frame <- function(...) {
    data.frame(..., check.rows = FALSE, check.names = FALSE, fix.empty.names = FALSE, stringsAsFactors = FALSE)
}


#' prediction methods
#'
#' Methods to compute the predictions and the marginal effects for
#' tobit1 objects
#'
#' @name prediction_margins
#' @param model a model fitted using `mhurdle`
#' @param data,at,vcov,calculate_se see `prediction::prediction`
#' @param what see `mhurdle:::predict.mhurdle`
#' @param ... further arguments, especially, a `what` argument can be
#'     provided and will be passed to `predict`
#' @details `tobit1` exports the `prediction::prediction` and
#'     `margins::margins` functions. `prediction` use the `predict`
#'     method to compute the predictions in a "tidy way", it returns
#'     the data frame provided for the predictions augmented by the
#'     predictions. `margins` compute the average marginal effect of
#'     every covariate. It uses the numerical derivatives of the
#'     predictions using the `prediction` function.
#' @importFrom prediction prediction find_data build_datalist
#' @importFrom stats na.pass delete.response
#' @method prediction mhurdle
#' @export
prediction.mhurdle <- function (model, data = find_data(model, parent.frame()), at = NULL, 
                               what = c("E", "Ep", "p"), vcov = stats::vcov(model), calculate_se = FALSE, 
                               ...){
    what <- match.arg(what)
    datalist <- build_datalist(data, at = at, as.data.frame = TRUE)
    at_specification <- attr(datalist, "at_specification")
    pred <- predict(model, newdata = datalist, what = what, se.fit = FALSE, ...)
    pred <- make_data_frame(datalist, fitted = pred, se.fitted = rep(NA_real_, nrow(datalist)))
    if (length(at)) vc <- rep(NA_real_, nrow(at_specification)) else vc <- NA_real_
    structure(pred, class = c("prediction", "data.frame"), at = if (is.null(at)) 
                                                                    at
    else at_specification, what = what, call = if ("call" %in% 
        names(model)) 
        model[["call"]]
    else NULL, model_class = class(model), row.names = seq_len(nrow(pred)), 
        vcov = vc, jacobian = NULL, weighted = FALSE)
}


## old.prediction.mhurdle <- function (model, data = find_data(model, parent.frame()), at = NULL, 
##                                type = "response", vcov = stats::vcov(model), calculate_se = FALSE, 
##                                ...){
##     type <- match.arg(type)
##     data <- data
##     if (missing(data) || is.null(data)) {
##         if (isTRUE(calculate_se)) {
##             pred <- predict(model, type = type, se.fit = TRUE, 
##                 ...)
##             pred <- make_data_frame(fitted = pred[["fit"]], se.fitted = pred[["se.fit"]])
##         }
##         else {
##             pred <- predict(model, type = type, se.fit = FALSE, 
##                 ...)
##             pred <- make_data_frame(fitted = pred, se.fitted = rep(NA_real_, 
##                                                                    length(pred)))
##         }
##     }
##     else {
##         cat("ICI\n")
##         model[["model"]] <- NULL
##         datalist <- build_datalist(data, at = at, as.data.frame = TRUE)
##         at_specification <- attr(datalist, "at_specification")
##         if (isTRUE(calculate_se)) {
##             tmp <- predict(model, newdata = datalist, type = type, 
##                 se.fit = TRUE, ...)
##             pred <- make_data_frame(datalist, fitted = tmp[["fit"]], 
##                 se.fitted = tmp[["se.fit"]])
##         }
##         else {
##             tmp <- predict(model, newdata = datalist, type = type, 
##                            se.fit = FALSE, ...)
##             print(head(tmp))
##             pred <- make_data_frame(datalist, fitted = tmp, se.fitted = rep(NA_real_, 
##                                                                             nrow(datalist)))
##             print(head(pred));stop()
##         }
##     }
##     if (isTRUE(calculate_se)) {
##         J <- NULL
##         model_terms <- delete.response(terms(model))
##         if (is.null(at)) {
##             model_frame <- model.frame(model_terms, data, na.action = na.pass, 
##                 xlev = model$xlevels)
##             model_mat <- model.matrix(model_terms, model_frame, 
##                 contrasts.arg = model$contrasts)
##             means_for_prediction <- colMeans(model_mat)
##             vc <- (means_for_prediction %*% vcov %*% means_for_prediction)[1L, 
##                 1L, drop = TRUE]
##         }
##         else {
##             datalist <- build_datalist(data, at = at, as.data.frame = FALSE)
##             vc <- unlist(lapply(datalist, function(one) {
##                 model_frame <- model.frame(model_terms, one, 
##                   na.action = na.pass, xlev = model$xlevels)
##                 model_mat <- model.matrix(model_terms, model_frame, 
##                   contrasts.arg = model$contrasts)
##                 means_for_prediction <- colMeans(model_mat)
##                 means_for_prediction %*% vcov %*% means_for_prediction
##             }))
##         }
##     }
##     else {
##         J <- NULL
##         if (length(at)) {
##             vc <- rep(NA_real_, nrow(at_specification))
##         }
##         else {
##             vc <- NA_real_
##         }
##     }
##     structure(pred, class = c("prediction", "data.frame"), at = if (is.null(at)) 
##         at
##     else at_specification, type = type, call = if ("call" %in% 
##         names(model)) 
##         model[["call"]]
##     else NULL, model_class = class(model), row.names = seq_len(nrow(pred)), 
##         vcov = vc, jacobian = J, weighted = FALSE)
## }
