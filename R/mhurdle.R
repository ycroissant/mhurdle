
#' Estimation of limited dependent variable models
#' 
#' mhurdle fits a large set of models relevant when the dependent variable is 0
#' for a part of the sample.
#' 
#' 
#' `mhurdle` fits models for which the dependent variable is zero for
#' a part of the sample. Null values of the dependent variable may
#' occurs because of one or several mechanisms : good rejection, lack
#' of ressources and purchase infrequency. The model is described
#' using a three-parts formula : the first part describes the
#' selection process if any, the second part the regression equation
#' and the third part the purchase infrequency process.  `y ~ 0 | x1 +
#' x2 | z1 + z2` means that there is no selection process.  `y ~ w1 +
#' w2 | x1 + x2 | 0` and `y ~ w1 + w2 | x1 + x2` describe the same
#' model with no purchase infrequency process. The second part is
#' mandatory, it explains the positive values of the dependant
#' variable. The `dist` argument indicates the distribution of the
#' error term. If `dist = "n"`, the error term is normal and (at least
#' part of) the zero observations are also explained by the second
#' part as the result of a corner solution. Several models described
#' in the litterature are obtained as special cases :
#' 
#' A model with a formula like `y~0|x1+x2` and `dist="n"` is the Tobit
#' model proposed by \insertCite{TOBIN/58}{mhurdle}.
#' 
#' `y~w1+w2|x1+x2` and `dist="l"` or `dist="t"` is the single hurdle
#' model proposed by \insertCite{CRAGG/71}{mhurdle}. With `dist="n"`,
#' the double hurdle model also proposed by
#' \insertCite{CRAGG/71}{mhurdle} is obtained. With `corr="h1"` we get
#' the correlated version of this model described by
#' \insertCite{BLUNDELL/87}{mhurdle}.
#' 
#' `y~0|x1+x2|z1+z2` is the P-Tobit model of
#' \insertCite{DEATO/IRISH/84}{mhurdle}, which can be a single hurdle
#' model if `dist="t"` or `dist="l"` or a double hurdle model if
#' `dist="n"`.
#' 
#' @name mhurdle
#' @aliases mhurdle
#' @param formula a symbolic description of the model to be fitted,
#' @param data a `data.frame`,
#' @param subset see [stats::lm()],
#' @param weights see [stats::lm()],
#' @param na.action see [stats::lm()],
#' @param start starting values,
#' @param dist the distribution of the error of the consumption
#'     equation: one of `"n"` (normal), `"ln"` (log-normal) `"bc"`
#'     (box-cox normal) and `"ihs"` (inverse hyperbolic sinus
#'     transformation),
#' @param h2 if `TRUE` the second hurdle is effective, it is not
#'     otherwise,
#' @param scaled if `TRUE`, the dependent variable is divided by its
#'     geometric mean,
#' @param corr a boolean indicating whether the errors of the
#'     different equations are correlated or not,
#' @param robust transformation of the structural parameters in order
#'     to avoid numerical problems,
#' @param check_gradient if `TRUE`, a matrix containing the analytical and
#'     the numerical gradient for the starting values are returned,
#' @param \dots further arguments.
#' @return
#' #' an object of class `c("mhurdle", "maxLik")`.
#' 
#' A `mhurdle` object has the following elements :
#' 
#' - coefficients: the vector of coefficients,
#' - vcov: the covariance matrix of the coefficients,
#' - fitted.values: a matrix of fitted.values, the first column being
#' the probability of 0 and the second one the mean values for the
#' positive observations,
#' - logLik: the log-likelihood,
#' - gradient: the gradient at convergence,
#' - model: a data.frame containing the variables used for the estimation,
#' - coef.names: a list containing the names of the coefficients in the
#' selection equation, the regression equation, the infrequency of purchase
#' equation and the other coefficients (the standard deviation of the error
#' term and the coefficient of correlation if `corr = TRUE`,
#' - formula: the model formula, an object of class `Formula`
#' - call: the call,
#' - rho: the lagrange multiplier test of no correlation.
#' 
#' @references
#'
#' \insertRef{BLUNDELL/87}{mhurdle}
#'
#' \insertRef{CRAGG/71}{mhurdle}
#'
#' \insertRef{DEATO/IRISH/84}{mhurdle}
#'
#' \insertRef{TOBIN/58}{mhurdle}
#'
#' @keywords regression
#' @examples
#' 
#' data("Interview", package = "mhurdle")
#' 
#' # independent double hurdle model
#' idhm <- mhurdle(vacations ~ car + size | linc + linc2 | 0, Interview,
#'               dist = "ln", h2 = TRUE, method = "bfgs")
#' 
#' # dependent double hurdle model
#' ddhm <- mhurdle(vacations ~ car + size | linc + linc2  | 0, Interview,
#'               dist = "ln", h2 = TRUE, method = "bfgs", corr = TRUE)
#' 
#' # a double hurdle p-tobit model
#' ptm <- mhurdle(vacations ~ 0 | linc + linc2 | car + size, Interview,
#'               dist = "ln", h2 = TRUE, method = "bfgs", corr = TRUE)
#' @importFrom Formula Formula
#' @importFrom survival survreg Surv
#' @importFrom truncreg truncreg
#' @importFrom stats binomial cor df.residual dnorm ecdf formula glm
#'     integrate lm lm.fit model.frame model.matrix model.response
#'     optimize pchisq pnorm printCoefmat qnorm quantile rnorm terms
#'     uniroot var .getXlevels deviance
#' @importFrom maxLik maxLik activePar
#' @export
mhurdle <- function(formula, data, subset, weights, na.action,
                    start = NULL, dist = c("ln", "n", "bc", "ihs"), h2 = FALSE,
                    scaled = TRUE, corr = FALSE, robust = TRUE,
                    check_gradient = FALSE, ...){
    fitted = TRUE
    dots <- list(...)
    oldoptions <- options(warn = -1)
    on.exit(options(oldoptions))
    cl <- match.call()
    posT <- as.list(cl) == "T" ; posF <- as.list(cl) == "F"
    cl[posT] <- TRUE           ; cl[posF] <- FALSE
    dist <- match.arg(dist)
    if (dist == "ln" & h2) dist <- "ln2"
    if (dist == "n" & ! h2) dist <- "tn"
    if (dist == "bc" & h2) dist <- "bc2"
    # 1. Compute the model.frame and the model.matrix

    if (! inherits(formula, "Formula")) formula <- Formula(formula)
    if (length(formula)[2] > 4) stop("at most 4 rhs should be provided in the formula")
    mf <- match.call(expand.dots = FALSE)
    mc <- match.call(expand.dots = TRUE)
    
    m <- match(c("formula", "data", "subset", "na.action", "weights"),
               names(mc), 0L)
    mf <- mc[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf$formula <- formula
    mf <- eval(mf, parent.frame())
    .terms <- terms(mf)
    X1 <- model.matrix(formula, data = mf, rhs = 1)
    X2 <- model.matrix(formula, data = mf, rhs = 2)
    X3 <- model.matrix(formula, data = mf, rhs = 3)
    X4 <- model.matrix(formula, data = mf, rhs = 4)
    # Remove the intercept if any for the heteroscedastic model
    if (ncol(X4)){
        int4 <- attr(terms(formula, data = mf, rhs = 4), "intercept")
        X4 <- X4[, - 1, drop = FALSE]
    }
    y <- model.response(mf)
    if (scaled){
        geomean <- exp(mean(log(y[y > 0])))
        y <- y / geomean
        attr(y, "geomean") <- attr(mf, "geomean") <- geomean
    }
    N <- length(y)
    if (length(X1) == 0) X1 <- NULL
    if (length(X2) == 0) stop("the second hurdle (consumption equation) is mandatory")
    if (length(X3) == 0) X3 <- NULL
    if (length(X4) == 0) X4 <- NULL
    h1 <- ! is.null(X1)
    h3 <- ! is.null(X3)
    
    #  2. One equation models
    if (! h1 && ! h3 && is.null(X4)){
        result <- one_equation_model(X2, y, dist, check_gradient = check_gradient, start = start)
        if (check_gradient){
            result$call <- cl
            return(result)
        }
        result$naive <- NULL#naive
        result$call <- cl
        result$model <- mf
        result$formula <- formula
        result$terms <- .terms
        result$R2 <- c(null = NA, positive = NA)
        return(result)
    }
    
    # 3. Selection single hurdle models without correlation that can
    # be estimated simply in two parts using seperate.mhurdle()
    if (h1 && !h3  && ! corr && !h2){
        result <- two_parts_model(X1, X2, y, dist = dist, check_gradient = check_gradient, start = start)
        if (check_gradient){
            result$call <- cl
            return(result)
        }

        P0 <- mean(y == 0)
        m <- mean(log(y[y > 0]))
        s2 <- mean( log(y[y > 0]) ^ 2) - m ^ 2
        parts <- c(log(P0), log(1 - P0), log(1 - P0) - m - 1 / 2 * log(s2) - 1 / 2 * log(2 * pi) - 1 /2)
        parts <- parts * c(P0, 1 - P0, 1 - P0) * N
        logLik.naive <- structure((parts[1] + parts[3]), nobs = length(y),
                                  parts = parts, 
                                  df = 3, class = "logLik")
        naive <- list(logLik = logLik.naive)

        lnL1c <- attr( naive$logLik, "parts")[1] + attr( naive$logLik, "parts")[2]
        lnL1u <- attr(result$logLik, "parts")[1] + attr(result$logLik, "parts")[2]
        lnL2c <- attr( naive$logLik, "parts")[3] - attr( naive$logLik, "parts")[2]
        lnL2u <- attr(result$logLik, "parts")[3] - attr(result$logLik, "parts")[2]

        R1 <- 1 - (lnL1c / lnL1u) ^ (2 / N * lnL1c)
        R2old <- 1 - (exp(lnL2c) / exp(lnL2u)) ^ (2 / (N * (1 - P0)))
        R2 <- 1 - exp(- 2 * (lnL2u - lnL2c) / (N * (1 - P0)))
        result$R2 <- c(null = unname(R1), positive = unname(R2), old = unname(R2old))
        result$naive <- naive
        result$call <- cl
        result$model <- mf
        result$formula <- formula
        result$terms <- .terms
        return(result)
    }
    # 5. Compute the starting values if not provided (use the linear
    # specification as the starting value for ihs and the log-linear
    # specification for Box-Cox)
    if (is.null(start)){
        dist.start <- dist
        if (dist %in% c("bc", "bc2", "ln2")) dist.start <- "ln"
        if (dist == "ihs") dist.start <- "n"
        start <- start.mhurdle(X1, X2, X3, y, dist.start)
        # in case of heteroscedasctic model, add K4 zeros to the start
        # vector and the intercept should be ln(sigma_o) (not sigma_o)
        # because of the exp form
        sd.pos <- getindex(X1, X2, X3, X4, corr, dist, which = "sd")
        if (robust) start[sd.pos] <- log(start[sd.pos])
        if (! is.null(X4)) start <- c(start[1:sd.pos], rep(0, ncol(X4)))
        
        # add shape and/or scale parameters
        if (corr){
            if (robust) rhoinit <- tan(0.1 * pi / 2) else rhoinit <- 0.1
            if (h1 + h3 == 2) start <- c(start, rho12 = rhoinit, rho13 = rhoinit, rho23 = rhoinit)
            else start <- c(start, rho = rhoinit)
        }
        if (dist %in% c("bc", "bc2", "ihs")) start <- c(start, tr =  0.1)
        if (dist %in% c("bc2", "ln2")){
            pos = 1
            if (robust) pos <- max(log(pos), -100)
            start <- c(start, pos = pos)
        }
    }
    else{
        if (robust){
            sd.pos <- getindex(X1, X2, X3, X4, corr, dist, which = "sd")
            if (dist %in% c("ln2", "bc2")) mu.pos <- getindex(X1, X2, X3, X4, corr, dist, which = "pos")
            rho.pos <- getindex(X1, X2, X3, X4, corr, dist, which = "corr")
            start[sd.pos] <- log(start[sd.pos])
            if (dist %in% c("ln2", "bc2")) start[mu.pos] <- max(log(start[mu.pos]), - 50)
            if (corr) start[rho.pos] <- tan(start[rho.pos] * pi / 2)
        }
    }

    result <- mhurdle.fit(start, X1, X2, X3, X4, y, gradient = TRUE,
                          dist = dist, corr = corr, robust = robust,
                          fitted = fitted,
                          check_gradient = check_gradient,
                          ...)

    if (check_gradient){
        result$call <- cl
        return(result)
    }
    
    # 3. Compute the naive model

    P0 <- mean(y == 0)
    dist.naive <- dist
    if (dist %in% c("ihs")) dist.naive <- "n"
    if (dist %in% c("bc", "bc2", "ln2")) dist.naive <- "ln"
    if (dist.naive != "ln"){
        Ec <- mean(y[y > 0])
        Vc <- var(y[y > 0])}
    else{
        Ec <- mean(log(y[y > 0]))
        Vc <- var(log(y[y > 0]))
    }

    start.naive <- c(rep(0.1, 1 + h1 + h3), 1)
    moments <- c(P0, Ec, Vc)
    naive <- maxLik::maxLik(lnl.naive, start = start.naive,
                    dist = dist.naive, moments = moments,
                    h1 = h1, h3 = h3)
    coef.naive <- naive$est
    parts <- attr(naive$max, "parts") * c(P0, 1 - P0, 1 - P0) * N
    logLik.naive <- structure(as.numeric(naive$max) * N, nobs = length(y),
                              parts = parts, 
                              df = length(coef.naive), class = "logLik")
    naive <- list(logLik = logLik.naive)
    # QDF for log-normal models

#    cat("sans h2\n")
    if (dist == "ln"){
        P0 <- mean(y == 0)
        m <- mean(log(y[y > 0]))
        s2 <- mean( log(y[y > 0]) ^ 2) - m ^ 2
        parts <- c(log(P0), log(1 - P0), log(1 - P0) - m - 1 / 2 * log(s2) - 1 / 2 * log(2 * pi) - 1 /2)
        parts <- parts * c(P0, 1 - P0, 1 - P0) * N
        logLik.naive <- structure((parts[1] + parts[3]), nobs = length(y),
                                  parts = parts, 
                                  df = 3, class = "logLik")
        naive <- list(logLik = logLik.naive)
    }

    ## cat("avec h2\n")
    ## param <- c(Np / N, mean(log(y[y > 0])), sd(log(y[y > 0])), 0.5)
    ## No <- sum(y == 0)
    ## Np <- sum(y > 0)
    ## f <- function(param){
    ##     Phio <- param[1]
    ##     m <- param[2]
    ##     sig <- param[3]
    ##     a <- param[4]
    ##     v <- No * log(1 - Phio * pnorm((m - log(a)) / sig)) + Np * (log(Phio) - log(sig) - 1 / 2  * log(2 * pi)) -
    ##         sum( log(y[y > 0] + a) + (log(y[y > 0] + a) - m) ^ 2 / (2 * sig ^ 2))
    ##     v
    ## }
    
    ## print(maxLik::maxLik(f, start = param, method = "bfgs"))
    
    ## g <- function(param){
    ##     m <- param[1]
    ##     sig <- param[2]
    ##     a <- param[3]
    ##     Pp <- Np / N
    ##         v <- No * log(1 - Pp) + Np * log(Pp) +
    ##             Np * (- pnorm((m - log(a)) / sig, log.p = TRUE) - log(sig) - 1 / 2  * log(2 * pi)) -
    ##             sum( log(y[y > 0] + a) + (log(y[y > 0] + a) - m) ^ 2 / (2 * sig ^ 2))
    ##     v
    ## }

    ## print(maxLik::maxLik(g, start = param[-1], method = "bfgs"))

    ## a <- 0.00000001
    ## h <- function(param){
    ##     m <- param[1]
    ##     sig <- param[2]
    ##     Pp <- Np / N
    ##         v <- No * log(1 - Pp) + Np * log(Pp) +
    ##             Np * (- pnorm((m - log(a)) / sig, log.p = TRUE) - log(sig) - 1 / 2  * log(2 * pi)) -
    ##             sum( log(y[y > 0] + a) + (log(y[y > 0] + a) - m) ^ 2 / (2 * sig ^ 2))
    ##     v
    ## }
    ## print(maxLik::maxLik(h, start = param[2:3], method = "bfgs"))

    ## stop()



    result$naive <- naive
    result$call <- cl
    result$formula <- formula
    result$model <- mf
    result$terms <- .terms
    result$xlevels = .getXlevels(.terms, mf)
    
    lnL1c <- attr( naive$logLik, "parts")[1] + attr( naive$logLik, "parts")[2]
    lnL1u <- attr(result$logLik, "parts")[1] + attr(result$logLik, "parts")[2]
    lnL2c <- attr( naive$logLik, "parts")[3] - attr( naive$logLik, "parts")[2]
    lnL2u <- attr(result$logLik, "parts")[3] - attr(result$logLik, "parts")[2]

    R1 <- 1 - (lnL1c / lnL1u) ^ (2 / N * lnL1c)
    R2old <- 1 - (exp(lnL2c) / exp(lnL2u)) ^ (2 / (N * (1 - P0)))
    R2 <- 1 - exp(- 2 * (lnL2u - lnL2c) / (N * (1 - P0)))
    result$R2 <- c(null = unname(R1), positive = unname(R2), old = unname(R2old))

    if (result$type_hessian == "neg_definite")
        result$dpar <- - sum(diag(solve(result$hessian, crossprod(result$gradient))))
    if (robust){
        result <- update(result, start = coef(result), iterlim = 0, robust = FALSE)
        result$call$robust <- TRUE
        result$call$iterlim <- NULL
        result$call$start <- NULL
    }
    result
}

mhurdle.fit <- function(start, X1, X2, X3, X4, y, gradient = FALSE, fit = FALSE,
                        dist = c("ln", "n", "tn", "bc", "ihs", "bc2", "ln2"),
                        corr = FALSE, robust = TRUE,  fitted = FALSE,
                        check_gradient = FALSE, ...){
    start.time <- proc.time()
    h1 <- ! is.null(X1)
    h3 <- ! is.null(X3)
    h4 <- ! is.null(X4)
    
    # fancy coefficients names
    sd.names <- "sd"
    
    if (corr){
        if (h1 & h3) rho.names <- c("corr12", "corr13", "corr23")
        else rho.names <- ifelse(h1, "corr12", "corr23")
    }
    else rho.names <- NULL
    if (dist %in% c("bc", "bc2", "ihs")) tr.names <- "tr" else tr.names <- NULL
    if (dist %in% c("ln2", "bc2")) mu.names <- "pos" else mu.names <- NULL

    coef.names <- list(h1   = colnames(X1),
                       h2   = colnames(X2),
                       h3   = colnames(X3),
                       sd   = sd.names,
                       h4   = colnames(X4),
                       corr = rho.names,
                       tr   = tr.names,
                       pos   = mu.names)

    start.names <- coef.names
    if (h1) start.names$h1 <- paste("h1", start.names$h1, sep = ".")
    start.names$h2 <- paste("h2", start.names$h2, sep = ".")
    if (h3) start.names$h3 <- paste("h3", start.names$h3, sep = ".")
    if (h4) start.names$h4 <- paste("h4", start.names$h4, sep = ".")
    names(start) <- Reduce("c", start.names)
    f <- function(param) mhurdle.lnl(param, X1 = X1, X2 = X2, X3 = X3, X4 = X4, y = y,
                                     gradient = TRUE, fitted = FALSE,
                                     dist = dist, corr = corr,
                                     robust = robust)

    if (check_gradient) return(compare_gradient(f, start))
    maxl <- maxLik::maxLik(f, start = start, control = list(lambdatol = 1E-20),  ...)
    nb.iter <- maxl$iterations
    convergence.OK <- maxl$code <= 2
    coefficients <- maxl$estimate

    h <- function(param) sum(mhurdle.lnl(param, X1 = X1, X2 = X2, X3 = X3, X4 = X4, y = y,
                                     gradient = FALSE, fitted = FALSE,
                                     dist = dist, corr = corr,
                                     robust = FALSE))
#print(robust)
#    H <- numDeriv::hessian(h, coefficients)
    if (fitted) fitted.values <- attr(mhurdle.lnl(coefficients, X1 = X1, X2 = X2, X3 = X3, X4 = X4, y = y,
                                                  gradient = FALSE, fitted = TRUE, robust = robust,
                                                  dist = dist, corr = corr), "fitted")
    else fitted.values <- NULL

    # contribution of every single observation to the likelihood and
    # its gradient (as an attribute)
    logLik <- f(coefficients)
    gradi <- attr(logLik, "gradi")
    logLik <- structure(as.numeric(logLik), df = length(coefficients),
                        parts = attr(logLik, "parts"),
                        nobs = length(y), class = "logLik")
    hessian <- maxl$hessian

    ev <- eigen(- hessian)$values

    type_hessian <- "neg_definite"
    if (any(is.complex(ev))){
        type_hessian <- "complex"
    } else {
        if (any(ev < 0)){
            type_hessian <- "not_neg_definite"
        } else {
            if (any(ev < 1E-07)){
                type_hessian <- "singular"
            }
        }
    }
    opg <- crossprod(gradi)
    ev_opg <- eigen(opg)$values
    
    elaps.time <- proc.time() - start.time

    if (type_hessian == "neg_definite") information <- - hessian else information <- opg

#    eps <- drop(crossprod(maxl$gradient, solve(information, maxl$gradient)))
    eps <- NA
    est.stat <- list(elaps.time = elaps.time,
                     nb.iter = nb.iter,
                     eps = eps,
                     method = maxl$type,
                     message = maxl$message
                     )
    class(est.stat) <- "est.stat"
    gtheta <- rep(1, length(coefficients))
    ## if (robust){
    ##     if (corr){
    ##         poscor <- sub.mhurdle(coef.names, "corr")
    ##         gtheta[poscor] <- 2 / pi / (1 + coefficients[poscor] ^ 2)
    ##         coefficients[poscor] <- atan(coefficients[poscor]) * 2 / pi
    ##     }
    ##     if (dist %in% c("bc2", "ln2")){
    ##         posmu <- sub.mhurdle(coef.names, "pos")
    ##         gtheta[posmu] <- exp(coefficients[posmu])
    ##         coefficients[posmu] <- exp(coefficients[posmu])
    ##     }
    ##     possd <- sub.mhurdle(coef.names, "sd")
    ##     gtheta[possd] <-  exp(coefficients[possd])
    ##     coefficients[possd] <- exp(coefficients[possd])
    ##     print(gtheta)
    ##     g1 <- apply(gradi, 2, sum)
    ##     gradi <- t(t(gradi) / gtheta)
    ##     g2 <- apply(gradi, 2, sum)
    ##     print(cbind(g1, g2))
    ## }
    if (robust){
        if (corr){
            poscor <- sub.mhurdle(coef.names, "corr")
            coefficients[poscor] <- atan(coefficients[poscor]) * 2 / pi
        }
        if (dist %in% c("bc2", "ln2")){
            posmu <- sub.mhurdle(coef.names, "pos")
            coefficients[posmu] <- exp(coefficients[posmu])
        }
        possd <- sub.mhurdle(coef.names, "sd")
        coefficients[possd] <- exp(coefficients[possd])
    }

    result <- list(coefficients  = coefficients,
                   fitted.values = fitted.values,
                   logLik        = logLik,
                   gradient      = gradi,
                   hessian       = hessian,
                   formula       = NULL,
                   terms         = NULL,
                   xlevels       = NULL,
                   model         = NULL,
                   coef.names    = coef.names,
                   call          = NULL,
                   est.stat      = est.stat,
                   naive         = NULL,
                   type_hessian  = type_hessian
                   )
    class(result) <- c("mhurdle", class(result))
    result
}


sanitize <- function(x, myeps = 1E-07, mymax = 1E02, string = c("", ""), replace = TRUE, verbal = TRUE){
    string <- paste("of", string[1], "in", string[2])
    if (replace){
        if (any(is.na(x))){
            if (verbal) cat(paste(sum(is.na(x)), "NA values", string, "replaced by 0\n"))
            x[is.na(x)] <- 0
        }
        if ( any(x > 0 & x < myeps)){
            if (verbal) cat(paste(sum(x > 0 & x < myeps), "values", string, "lower than",  myeps,"replaced by", myeps, "\n"))
            x[x > 0 & x < myeps] <- myeps
        }
        if (any(x < - mymax)){
            if (verbal) cat(paste(sum(x < - mymax), "values", string, "lower than", - mymax, "replaced by", - mymax, "\n"))
            x[x < - mymax] <- - mymax
        }
        if (any(x > mymax)){
            if (verbal) cat(paste(sum(x >  mymax), "values", string, "greater than",  mymax, "replaced by", mymax, "\n"))
            x[x > mymax] <- mymax
        }
    }
    x
}
    
