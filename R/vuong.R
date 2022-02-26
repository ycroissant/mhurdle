#' Vuoung test for non-nested models
#' 
#' The Vuong test is suitable to discriminate between two non-nested models.
#' 
#' 
#' @aliases vuongtest
#' @param x a first fitted model of class \code{"mhurdle"},
#' @param y a second fitted model of class \code{"mhurdle"},
#' @param type the kind of test to be computed,
#' @param true_model a boolean, \code{TRUE} if one of the models is asumed to be the
#' true model,
#' @param variance the variance is estimated using the \code{centered} or
#' \code{uncentered} expression,
#' @param matrix the W matrix can be computed using the general expression
#' \code{large} or the reduced matrix \code{reduced} (only relevant for the
#' nested case),
#' @return an object of class \code{"htest"}
#' @seealso \code{vuong} in package \code{pscl}.
#' @references Vuong Q.H. (1989) Likelihood ratio tests for model selection and
#' non-nested hypothesis, Econometrica, vol.57(2), pp.307-33.
#' @keywords htest
#' @importFrom CompQuadForm davies
#' @examples
#' 
#' data("Interview", package = "mhurdle")
#' # dependent double hurdle model
#' dhm <- mhurdle(vacations ~ car + size | linc + linc2 | 0, Interview,
#'               dist = "ln", h2 = TRUE, method = "bhhh", corr = TRUE)
#' 
#' # a double hurdle p-tobit model
#' ptm <- mhurdle(vacations ~ 0 | linc + linc2 | car + size, Interview,
#'               dist = "ln", h2 = TRUE, method = "bhhh", corr = TRUE)
#' vuongtest(dhm, ptm)
#' @export
vuongtest <- function(x, y,
                      type = c("non-nested", "nested", "overlapping"),
                      true_model = FALSE,
                      variance = c("centered", "uncentered"),
                      matrix = c("large", "reduced")){# match the arguments
    type <- match.arg(type)
    variance <- match.arg(variance)
    matrix <- match.arg(matrix)
    pds <- function(x) paste(deparse(substitute(x)))
    data.name <- paste(paste(deparse(substitute(x))),
                       paste(deparse(substitute(x))),
                       sep = " - ")
    
    # for convenience, call f the larger model and g the other one if
    # the models are nested, else call the first model f and g for
    # consistency with Vuong paper. Check also that the models are
    # really nested, otherwise stop
    f <- x
    ;   g <- y
    if (type == "nested" && (length(coef(x)) < length(coef(y)))){
        f <- y ;  g <- x
    }
    
    Kf <- length(coef(f))      ;  Kg <- length(coef(g))
    N <- nobs(f)

    # irrelevant or not implemented tests
    if (nobs(g) != nobs(f)) stop("the number of observations of the two models differ")
    if (matrix == "reduced" && type != "nested")
        stop("the reduced matrix is only relevant for nested models")
    if (type == "nested")
        if (! all(names(coef(g)) %in% names(coef(f))) || (Kf == Kg))
            stop("the two models are not nested")
    if (type == "non-nested" & true_model)
        stop("non-nested test with the hypothese of a true model is not implemented")
    
    # extract the characteristics of the models, using the meat,
    # bread, estfun (sandwich) and llcont (nonnest2)
    lf <- llcont(f)    ; lg <- llcont(g)
    Kf <- length(lf)   ; Kg <- length(lg)
    Bf <- meat(f)      ; Bg <- meat(g)
    Af1 <- - bread(f)  ; Ag1 <- - bread(g)
    Bfg <- crossprod(estfun(f), estfun(g)) / N
    
    # compute the likelihood ratio
    LR <- as.numeric(logLik(f) - logLik(g))
    # and its variance (if centered, eq 4.2 else eq 4.3 in Vuong 1989,
    # p. 314)
    w2 <- ifelse(variance == "centered",
                 1 / N * sum( (lf - lg) ^ 2) - (1 / N * LR) ^ 2,
                 1 / N * sum( (lf - lg) ^ 2))

    # construct the large or reduced matrix and its eigen values (see
    # Vuong (1989) eq 3.6 (for the large) and eq 7.4 for the reduced
    # matrix
    if (matrix == "large"){
        W <- rbind(cbind( -  Bf   %*% Af1, - Bfg %*% Ag1),
                   cbind(  t(Bfg) %*% Af1,   Bg  %*% Ag1))
#        same as:        
#        B <- cbind(estfun(f), - estfun(g))
#        B <- crossprod(B) / nrow(B)
#        Ainv <- bdiag(- bread(f), bread(g))
#        W <- - crossprod(B, Ainv)
        print(W[1:3, 1:10])
    }
    else{
        common.coef <- names(coef(f)) %in% names(coef(g))
        D <- t(diag(1, Kf)[common.coef, ])
        W <- Bf %*% (D %*% Ag1 %*% t(D) - Af1)
    }
    Z <- eigen(W)$values
    print(sum(Z))
    
    # non nested test
    if (type == "non-nested"){
        method <- "Vuong Test (non-nested)"
        statistic <- c(z = LR / sqrt(N * w2))
        # if the stat is negative, the p-value is the lower tail,
        # otherwise it is the upper tail
        pval <- pnorm(statistic, lower.tail = statistic < 0)
        parameter <- NULL
    }

    # nested test
    if (type == "nested"){
        method <- "Vuong Test (nested)"
        statistic <- 2 * LR
        if (true_model){
            # the classical likelihood ratio test
            parameter <- c(df = Kf - Kg)
            names(statistic) <- "chisq"
            pval <- pchisq(statistic, parameter, lower.tail = FALSE)
        }
        else{
            parameter <- c(sev = sum(Z))
            names(statistic) <- "wchisq"
            pval <- CompQuadForm::davies(statistic, Z)$Qq
        }
    }
    
    # overlapping test
    if (type == "overlapping"){
        method <- "Vuong Test (overlapping)"
        if (true_model){
            statistic <- c(wchisq = 2 * LR)
            parameter <- c(sev = sum(Z))
            pval <- CompQuadForm::davies(statistic, Z)$Qq
            pval <- 2 * ifelse(statistic < 0, 1 - pval, pval)
        }
        else{
            # first compute the variance statistic (eq. 4.6 p. 315 of
            # Vuong, 1989)"
            statistic <- c(wchisq = N * w2)
            pval <- CompQuadForm::davies(statistic, Z ^ 2)$Qq
            method <- paste(method, ": first step (variance test), p-value = ", round(pval, 3), "\n")
            parameter <- c(sev = sum(Z ^ 2))
            variance_test <- structure(list(statistic = statistic,
                                            method = "variance test",
                                            p.value = pval,
                                            data.name = data.name,
                                            parameter = parameter))
            # if H_o is rejected, then compute the likelihood-ratio
            # like test. The LR statistic can be either positive or
            # negative, depending on the order of the models. The test
            # is two-sided so that the p-value is twice the lower tail
            # if the statistic is negative and twice the upper tail
            # otherwise
            statistic <- c(z = LR / sqrt(N * w2))
            parameter <- NULL
            pval <- 2 * ifelse(statistic < 0, pnorm(statistic), pnorm(statistic, lower.tail = FALSE))
        }
    }

    result <- list(statistic = statistic,
                   method = method,
                   p.value = pval,
                   data.name = data.name,
                   parameter = parameter)
    class(result) <- "htest"
    result$Z <- Z

    result
}

## #' Shi test for non-nested models
## #' 
## #' The Shi test correct the bias of the Vuong test
## #' 
## #' 
## #' @aliases shitest
## #' @param x a first fitted model of class \code{"mhurdle"},
## #' @param y a second fitted model of class \code{"mhurdle"},
## #' @param size the size of the test,
## #' @param pval should the p-value be computed ?
## #' @param type the kind of test to be computed,
## #' @param ndraws the number of draws for the simulations,
## #' @param diffnorm a creuser,
## #' @param seed the seed,
## #' @param print.level the level of details to be printed.
## #' @return an object of class \code{"htest"}
## #' @seealso \code{vuong} in package \code{pscl}.
## #' @references Vuong Q.H. (1989) Likelihood ratio tests for model selection and
## #' non-nested hypothesis, Econometrica, vol.57(2), pp.307-33.
## #' @keywords htest
## #' @examples
## #' 
## #' data("Interview", package = "mhurdle")
## #' # dependent double hurdle model
## #' dhm <- mhurdle(vacations ~ car + size | linc + linc2 | 0, Interview,
## #'               dist = "ln", h2 = TRUE, method = "bhhh", corr = TRUE)
## #' 
## #' # a double hurdle p-tobit model
## #' ptm <- mhurdle(vacations ~ 0 | linc + linc2 | car + size, Interview,
## #'               dist = "ln", h2 = TRUE, method = "bhhh", corr = TRUE)
## #' shitest(dhm, ptm)
## #' @export
## shitest <- function(x, y, size = 0.05, pval = TRUE,
##                     type = c("non-nested", "nested", "overlapping"),
##                     ndraws = 1E04, diffnorm = 0.1, seed = 1,
##                     numbers = NULL,
##                     print.level = 0){
##     # if pval is TRUE, size is adjusted, otherwise it is fixed

##     type <- match.arg(type)
##     eigentol <- 1E-12
##     data.name <- c(
##         paste(deparse(substitute(x))),
##         paste(deparse(substitute(y)))
##     )
##     data.name <- paste(data.name, collapse = "-")
##     set.seed(seed)

##     # the two models are renamed f and g like in Vuong paper, the
##     # generics llcont, estfun and bread are used to extract the
##     # relevant features of the log-likelihood
##     f <- x                    ;  g  <- y    
##     N <- length(llcont(f))
##     K.f <- ncol(estfun(f))    ;  K.g <- ncol(estfun(g))
##     K <- K.f + K.g

##     # Compute the LR statistic as an average and its variance
##     LR <- as.numeric(logLik(f) - logLik(g)) / N
##     w2 <- sum((llcont(f)- llcont(g)) ^ 2) / N - LR ^ 2#(LR / N) ^ 2
    
##     solveA <- bdiag(- bread(f), bread(g))
##     grad.tot <- cbind(estfun(f), - estfun(g))
##     # substract the mean of the gradient (usefull if the gradient is
##     # not null at the optimum)
##     grad.tot <- t(t(grad.tot) - apply(grad.tot, 2, mean))
    
##     B <- crossprod(grad.tot) / N# + diag(rep(eigentol, K))
##     sqrtB <- sqrtm(B)
##     V <- sqrtB %*% solveA %*% sqrtB
##     V <- eigen(V, symmetric = TRUE)$values

##     # Compute the 4 matrices that are used to compute empirically the
##     # distribution of the Vuong statistic
##     if (type != "nested"){
##         rho <- as.numeric((abs(V) - max(abs(V))) == 0)
##         rho <- rho / sqrt(sum(rho))
##     }
##     else rho <- rep(0, length(V))
##     SIGMA <- rbind(c(1 , rho),
##                    cbind(rho, diag(rep(1, K))))

##     SIGMA <- rbind(c(1, rho),
##                    cbind(rho, diag(K)))
##     SIGMA[c(FALSE, as.logical(rho)), ] <- 0

##     if (is.null(numbers)) Z <- matrix(rnorm(ndraws * (K + 1)), ndraws)
##     else Z <- numbers

    
##     Z <- Z %*% SIGMA
##     ZL <- Z[,   1]
##     ZP <- Z[, - 1]
##     trVsq <- sum(V ^ 2)
##     Vnmlzd <- V / sqrt(trVsq)
##     A1 <- ZL
##     A2 <- as.numeric(apply(ZP, 1, function(x) sum(Vnmlzd * x ^ 2) / 2) -
##                      sum(Vnmlzd) / 2)
##     A3 <- as.numeric(apply(ZP, 1, function(x) sum(Vnmlzd * rho * x)))
##     A4 <- as.numeric(apply(ZP, 1, function(x) sum(x ^ 2 * Vnmlzd ^ 2)))

##     Tmod <- function(sigma, cst){
##         # R draws in the distribution of the Vuong statistic
##         num <- sigma * A1 - A2
##         denom <- sigma ^ 2 - 2 * sigma * A3 + A4 + cst
##         num / sqrt(denom)
##     }

##     quant <- function(sigma, cst, size)
##         # compute the empirical quantile of level 1 - size in the
##         # distribution of the Vuong statistic
##         as.numeric(quantile(abs(Tmod(sigma, cst)), 1 - size))
        
##     sigstar <- function(cst, size)
##         # compute the value of sigma which maximize the empirical
##         # quantile of level 1 - size of the distribution of the Vuong
##         # statistic
##         optimize(function(x) quant(x, cst, size), c(0, 5), maximum = TRUE)$maximum
    
##     seekpval <- function(size){
##         # for a given size, compute the constant so that the critical
##         # value equals the target
##         c.value <- 0
##         cv.value <- quant(sigstar(0, size), 0, size)
##         # what is the interest of that equation
##         cv.normal <- max(qnorm(1 - size / 2), quantile(abs(ZL), 1 - size))
##         cv.target <- cv.normal + diffnorm
##         if (cv.value < cv.target){
##             cv.value <- max(cv.value, cv.normal)
##         }
##         else {
## #            asig <- sigstar(2, 0.05)
## #            afr <- quant(asig, 2, 0.05) - cv.target
## #            print(asig)
## #            print(afr)
##             froot <- function(cst) quant(sigstar(cst, size), cst, size) - cv.target
##             zo <- uniroot(froot, c(0, 10))
##             c.value <- zo$root
##             cv.value <- quant(sigstar(c.value, size), c.value, size)
##         }

##         LRmod <- LR + sum(V) / (2 * N)
##         w2mod <- w2 + c.value * sum(V ^ 2) / N
##         Tnd <-  sqrt(N) * LRmod / sqrt(w2mod)
##         pvalue <- 1 - ecdf(abs(Tmod(sigstar(c.value, size), c.value)))(abs(Tnd))
##         list(pvalue = pvalue, stat = Tnd,
##              constant = c.value, cv = cv.value)
##     }
##     if (type == "nested"){
##         LRmod <- LR + sum(V) / (2 * N)
##         Tnd <- sqrt(N) * (LRmod) / sqrt(w2)
##         if (Tnd < 0) pvalue <- ecdf(Tmod(0, 0))(Tnd)
##         else pvalue <- 1 - ecdf(Tmod(0, 0))(Tnd)
##         TndVuong <- list(statistic = c(z = Tnd),
##                          method = "Non-degenerate Vuong test",
##                          p.value = pvalue,
##                          data.name = data.name,
##                          alternative = "different models")
##     }
##     else{
##         if (! pval) results <- seekpval(size)
##         else{
##             froot <- function(alpha) seekpval(alpha)$pvalue - alpha
##             if (froot(1E-100) < 0) results <- list(pvalue = 0, stat = Inf, constant = 0, cv = 0)
##             else{
##                 pvalue <- uniroot(froot, c(1E-100, 1-1E-100))
##                 results <- seekpval(pvalue$root)
##             }
##         }
##         TndVuong <- list(statistic = c(z = as.numeric(results$stat)),
##                          method = "Non-degenerate Vuong test",
##                          p.value = results$pvalue,
##                          data.name = data.name,
##                          alternative = "different models",
##                          parameters = c(constant = results$constant,
##                                         'crit-value' = results$cv,
##                                         'sum e.v.' = sum(V)))
##         if (! pval) TndVuong$size <- size
##     }
##     structure(TndVuong, class = "htest")
## }

## sqrtm <- function(x){
##     eigenx <- eigen(x)
## #    if (any(eigenx$values < 0)) stop("negative eigen values")
##     eigenx$vectors %*% diag(sqrt(abs(eigenx$values))) %*% t(eigenx$vectors)
## #    eigenx$vectors %*% diag(sqrt(eigenx$values)) %*% t(eigenx$vectors)
## }


## ## # simulated probabilities and quantiles for the weighted chi-squares
## ## # distribution
## ## pwchisq <- function(q, weights, lower.tail = TRUE, R = 1E03, seed = 1){
## ##     set.seed(seed)
## ##     K <- length(weights)
## ##     e <- matrix(rnorm(R * K) ^ 2, R, K)
## ##     wcs <- apply(e, 1, function(x) sum(x * weights))
## ##     F <- ecdf(wcs)
## ##     ifelse(lower.tail, F(q), 1 - F(q))
## ## }
## ## qwchisq <- function(p, weights, lower.tail = TRUE, R = 1000, seed = 1){
## ##     set.seed(1)
## ##     K <- length(weights)
## ##     e <- matrix(rnorm(R * K) ^ 2, R, K)
## ##     wcs <- apply(e, 1, function(x) sum(x * weights))
## ##     ifelse(lower.tail, quantile(wcs, p), quantile(wcs, 1 - p))
## ## }  
