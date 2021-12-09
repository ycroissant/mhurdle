one_equation_model <- function(X2, y, dist = NULL, start = NULL, check_gradient = FALSE, robust = TRUE){
    if (dist %in% c("ln", "tn", "bc")) stop("the specified model doesn't allow zero observations")        
    if (dist == "n"){
        start_length <- ncol(X2) + 1
        if (is.null(start)) start <- c(rep(0.1, ncol(X2)), 1)
        names(start) <- c(colnames(X2), "sigma")

        # the tobit model, efficiently estimated using Survreg from
        # the survival package
        tobit.lnl <- function(param, fitted = FALSE){
            beta <- param[1 : (length(param) - 1)]
            sigma <- param[length(param)]
            bX <- as.numeric(crossprod(t(X2), beta))
            lnL <- (y == 0) * pnorm(- bX / sigma, log.p = TRUE) +
                (y > 0) * (- log(sigma) +  dnorm(   (y - bX) / sigma, log   = TRUE))
            gbeta <- - (y == 0) * mills(- bX / sigma) / sigma + (y > 0) * (y - bX) / sigma ^ 2
            gsigma <- (y == 0) * mills(- bX / sigma) * bX / sigma ^ 2 +
                (y > 0) * (- 1 / sigma + (y - bX) ^ 2 / sigma ^ 3)
            gradi <- cbind(gbeta * X2, sigma = gsigma)
            attr(lnL, "gradient") <- gradi
            lnL
        }
        if (check_gradient) return(compare_gradient(tobit.lnl, start))

        tobit <- survival::survreg(survival::Surv(y, y > 0, type = "left") ~ X2 - 1, dist = "gaussian")
        beta <- coef(tobit)
        sigma <- tobit$scale
        bX <- as.numeric(crossprod(t(X2), beta))
        resid <- y - bX
        z <- tobit.lnl(c(beta, sigma))
        gradi <- attr(z, "gradient")
        result <- list(coefficients = c(beta, sigma = sigma),
                       vcov = solve(crossprod(gradi)),
                       logLik = logLik(tobit),
                       gradient = gradi,
                       fitted = cbind(1 - pnorm(bX / sigma), bX + sigma * mills(bX / sigma)))
    }
    if (dist == "ln2"){
        start_length <- ncol(X2) + 2
        if (is.null(start)){
            thelm <- lm(log(y + 1) ~ X2 - 1)
            sigma <- summary(thelm)$sigma
            beta <- coef(thelm)
            start <- c(beta, sigma = sigma, alpha = 1)
        }
        else if (length(start) != start_length) stop("the starting values vector has a wrong length")
        ltobit.lnl <- function(param, fitted = FALSE){
            beta <- param[1 : (length(param) - 2)]
            sigma <- param[length(param) - 1]
            alpha <- param[length(param)]
            bX <- as.numeric(crossprod(t(X2), beta))
            lnL <- (y == 0) * pnorm( (log(    alpha) - bX) / sigma, log.p = TRUE) +
                (y > 0) * ( dnorm(   (log(y + alpha) - bX) / sigma, log   = TRUE) -
                               log(y + alpha) - log(sigma)
                           )
            mu <- mills( (log(alpha) - bX) / sigma)
            gbeta <- - (y == 0) * mu / sigma +
                (y > 0) * (log(y + alpha) - bX) / sigma ^ 2
            gsigma <- - (y == 0) * mu * (log(alpha) - bX) / sigma ^ 2 +
                (y > 0) * ( - 1 / sigma + (log(alpha + y) - bX) ^ 2 / sigma ^ 3)
            galpha <- (y == 0) * mu / (alpha * sigma) -
                (y > 0) * ( 1 / (y + alpha) + ( log(y + alpha) - bX) / (sigma ^ 2 * (y + alpha)))
            attr(lnL, "gradient") <- cbind(gbeta * X2, gsigma, galpha)
            if (fitted){
                P0 <- pnorm( (log(alpha) - bX) / sigma)
                Ep <- exp(bX + 0.5 * sigma ^ 2) * pnorm( (bX - log(alpha)) / sigma + sigma) / (1 - P0)
                attr(lnL, "fitted") <- cbind(P0, Ep)
            }
            lnL
        }
        if (check_gradient) return(compare_gradient(ltobit.lnl, start))
        result <- maxLik::maxLik(ltobit.lnl, start = start, method = "bhhh", print.level = 0)
        fitted <- attr(ltobit.lnl(coef(result), fitted = TRUE), "fitted")
        result <- list(coefficients = coef(result),
                       vcov = vcov(result),
                       logLik = logLik(result),
                       gradient = result$gradientObs,
                       fitted = fitted)

    }
    if (dist == "bc2"){
        if (is.null(start)) start <- c(rep(0.1, ncol(X2)), 0.1, 1)
        names(start) <- c(colnames(X2), "alpha", "sigma")
        f2 <- function(x) boxcox.lnl(x, X2, y, gradient = TRUE, truncated = FALSE, robust = FALSE)
        if (check_gradient) return(compare_gradient(f2, start))
        result <- boxcox.fit(X2, y, robust = robust, start = start,
                             check_gradient = check_gradient, truncated = FALSE)
    }

    coef.names <- list(h1    = NULL,
                       h2    = colnames(X2),
                       h3    = NULL,
                       sd    = "sd",
                       h4    = NULL,
                       corr  = NULL,
                       tr    = NULL,
                       pos   = NULL)
    if (dist == "bc2") coef.names$tr <- "tr"
    if (dist %in% c("bc2", "ln2")) coef.names$pos <- "pos"

    
    result <- list(coefficients = result$coefficients, 
                   vcov = result$vcov,
                   fitted.values = result$fitted,
                   logLik = result$logLik,
                   gradient = result$gradient,
                   model = NULL,
                   formula = NULL,
                   coef.names = coef.names,
                   call = NULL
                   )
    class(result) <- c("mhurdle","maxLik")
    result

    
    return(result)
}
        

# Compute the estimation of hurdle models in the cases where it can be
# done using two independent estimations (a binomial logit model and a
# normal/log-normal/truncated linear model). This is relevant for
# uncorrelated models with selection

two_parts_model <- function(X1, X2, y, dist = NULL, start = NULL, check_gradient = FALSE){
    if (is.null(start)) start <- c(rep(0.1, ncol(X1) + ncol(X2)), 1)
    names(start) <- c(paste("h1", colnames(X1), sep = "."),
                      paste("h2", colnames(X2), sep = "."),
                      "sd")
    
    yp <- y[y > 0]
    q <- 2 * (as.numeric(y > 0)) - 1
    
    lnl.probit <- function(x){
        bX1 <- as.numeric(crossprod(x, t(X1)))
        lnL <- pnorm(q * bX1, log.p = TRUE)
        attr(lnL, "gradient") <- q * mills(q * bX1) * X1
        attr(lnL, "hessian") <- crossprod(q ^ 2 * dmills(q * bX1) * X1, X1)
        attr(lnL, "parts") <- c(lnLNull = sum(pnorm(- bX1[y == 0], log.p = TRUE)),
                                lnLOne  = sum(pnorm(  bX1[y > 0], log.p = TRUE)))
        lnL
    }

    lnl.ln <- function(x){
        sigma <- x[length(x)]
        beta <- x[- length(x)]
        bX2 <- as.numeric(crossprod(beta, t(X2)))
        resid <- log2(y) - bX2
        lnL <- (- log(sigma) - log2(y) + dnorm(resid / sigma, log = TRUE)) * (y > 0)
        lnL_beta <- resid / sigma ^ 2
        lnL_sigma <- resid ^ 2 / sigma ^ 3 - 1 / sigma
        attr(lnL, "gradient") <- cbind(lnL_beta * X2, sigma = lnL_sigma) * (y > 0)
        attr(lnL, "hessian") <- bdiag(- 1 / sigma ^ 2 * crossprod(X2),
                                      - 3 * sum( (y - bX2) ^ 2) / sigma ^ 4 + 1 / sigma ^ 2)
        lnL
    }

    lnl.boxcox <- function(x){
        boxcox.lnl(x, X2, yp, truncated = TRUE, robust = FALSE, gradient = TRUE)
    }
    
    lnl.tn <- function(x){
        sigma <- x[length(x)]
        beta <- x[- length(x)]
        bX2 <- as.numeric(crossprod(beta, t(X2)))
        resid <- y - bX2
        lnL <- - log(sigma) + dnorm(resid / sigma, log = TRUE) - pnorm(bX2 / sigma, log.p = TRUE)
        lnL <- lnL * (y > 0)
        mu <- mills( bX2 / sigma)
        dmu <- dmills( bX2 / sigma)
        lnL_beta <- resid / sigma ^ 2 - mu / sigma
        lnL_sigma <- - 1 / sigma + resid ^ 2 / sigma ^ 3 + mu * bX2 / sigma ^2
        lnL_beta_beta <- crossprod(- (bX2 + dmu) / sigma ^ 2 * X2, X2)
        lnL_beta_sigma <- - 2 * apply(resid / sigma + mu + dmu * bX2 / sigma * X2, 2, sum) / sigma ^ 2
        lnL_sigma_sigma <- sum(1 - 3 * resid ^ 2 / sigma ^ 2 - dmu * bX2 / sigma ^ 2 - 2 * mu * bX2 / sigma) / sigma ^ 2
        attr(lnL, "gradient") <- cbind(lnL_beta * X2, sigma = lnL_sigma) * (y > 0)
        attr(lnL, "hessian") <- rbind(cbind(lnL_beta_beta * crossprod(X2), lnL_beta_sigma),
                                      lnL_beta_sigma, lnL_sigma_sigma)
        lnL
    }

    lnl.tn.olsen <- function(x, gradient = TRUE, hessian = TRUE){
        theta <- x[  length(x)]
        gamma <- x[- length(x)]
        gX2 <- as.numeric(crossprod(gamma, t(X2)))
        lnL <-  dnorm(theta * y - gX2, log = TRUE) + log(theta) - pnorm(gX2, log.p = TRUE)
        lnL <- - sum(lnL * (y > 0))
        lnL_gamma <- - mills(gX2) + (theta * y - gX2)
        lnL_theta <- 1 / theta - y * (theta * y - gX2)
        lnL_gamma_gamma <- - crossprod( (1 + dmills(gX2)) * X2 * (y > 0), X2 * (y > 0))
        lnL_gamma_theta <- apply(y * X2 * (y > 0), 2, sum)
        lnL_theta_theta <- - sum( (1  / theta ^ 2 + y ^ 2) * (y > 0))
        attr(lnL, "gradient") <- - apply(cbind(lnL_gamma * X2, theta = lnL_theta) * (y > 0), 2, sum)
        attr(lnL, "hessian") <- - rbind(cbind(lnL_gamma_gamma, lnL_gamma_theta),
                                        c(lnL_gamma_theta, lnL_theta_theta))
        lnL
    }
        
    
    if (dist == "tn") f_second <- function(x) lnl.tn(x)
    if (dist == "ln") f_second <- function(x) lnl.ln(x)
    if (dist == "bc") f_second <- function(x) lnl.boxcox(x)
    fs_second <- function(x) sum(as.numeric(f_second(x)))
    fs_first <- function(x) sum(as.numeric(lnl.probit(x)))


    lnl_two_parts <- function(x){
        coefs_one <- x[1:ncol(X1)]
        coefs_two <- x[(ncol(X1) + 1):(ncol(X1) + ncol(X2) + 1)]
        one <- lnl.probit(coefs_one)
        two <- f_second(coefs_two)
        lnL <- as.numeric(one) + as.numeric(two)
        grad_one <- attr(one, "gradient")
        grad_two <- attr(two, "gradient")
        grad <- cbind(grad_one, grad_two)
        attr(lnL, "gradient") <- grad
        attr(lnL, "parts") <- c(lnLNull = unname(attr(one, "parts")["lnLNull"]),
                                lnLOne  = unname(attr(one, "parts")["lnLOne"]),
                                lnLPos  = sum(as.numeric(two)) + unname(attr(one, "parts")["lnLOne"]))
        lnL
    }
    if (check_gradient) return(compare_gradient(lnl_two_parts, start))

    first <- glm(y != 0 ~ X1 - 1, family = binomial(link = "probit"))
    coefs_first <- coef(first)
    first <- lnl.probit(coefs_first)
    L.null <- as.numeric(first)
    lp <- as.numeric(X1 %*% coefs_first)
    gradient_first <- attr(first, "gradient")
    hessian_first <- attr(first, "hessian")
    
    # Computation of the likelihood for positive observations for log-normal distribution

    if (dist == "ln"){
        second <- lm(log(yp) ~ X2[y > 0, ] - 1)
        beta2 <- coef(second)
        sigma <- sqrt(deviance(second) / nobs(second))
        coefs_second <- c(beta2, sigma = sigma)
    }
    else{
        if (dist == "tn"){
            start <- rep(0.1, ncol(X2) + 1)
            # la version olsen ne marche pas pour une raison inconnue
#            za <- nlm(lnl.tn.olsen, start, print.level = 2L)
            second <- truncreg::truncreg(yp ~ X2[y > 0, ] - 1)
        }
        
        if (dist == "bc") second <- boxcoxreg(yp ~ X2[y > 0, ] - 1, method = "bhhh", print.level = 0)
        coefs_second <- coef(second)
    }
    gradient_second <- attr(f_second(coefs_second), "gradient")
    hessian_second <- attr(f_second(coefs_second), "hessian")
    gradient <- cbind(gradient_first, gradient_second)
    
    L.pos <- logLik(second)
    

    coef.names <- list(h1    = colnames(X1),
                       h2    = colnames(X2),
                       h3    = NULL,
                       sd    = "sd",
                       h4    = NULL,
                       corr  = NULL,
                       tr    = NULL,
                       pos   = NULL)
    if (dist == "bc") coef.names$tr <- "tr"
    logLik <- structure(sum(L.null) + sum(L.pos), df = length(coef), nobs = length(y), class = "logLik")
    result <- list(coefficients = c(coefs_first, coefs_second), 
                   fitted.values = NULL,
                   logLik = logLik,
                   gradient = gradient,
                   hessian = bdiag(hessian_first, hessian_second),
                   model = NULL,
                   formula = NULL,
                   coef.names = coef.names,
                   call = NULL
                   )
    class(result) <- c("mhurdle","maxLik")
    result
}


newton <- function(fun, coefs, trace, direction = c("min", "max"), ...){
    if (trace){
        cat("Initial values of the coefficients:\n")
    }
    direction <- match.arg(direction)
    i <- 1
    eps <- 10
    while (abs(eps) > 1E-07){
        f <- fun(coefs, gradient = TRUE, hessian = TRUE, ...)
        g <- attr(f, "gradient")
        if (is.matrix(g)) g <- apply(g, 2, sum)
        h <- attr(f, "hessian")
        if (direction == "max"){
            f <- - f
            g <- - g
            h <- - h
        }
        lambda <- 1
        newcoefs <- coefs - as.numeric(solve(h, g))
        as_scalar <- function(x) sum(as.numeric(x))
        while (as_scalar(- fun(newcoefs, ...)) > as_scalar(f)){
            lambda <- lambda / 2
            if(trace) cat(paste("function is increasing, lambda set to:", lambda, "\n"))
            newcoefs <- coefs - lambda * as.numeric(solve(h, g))
        }
        eps <- as.numeric(crossprod(solve(h, g), g))
                if (trace) cat(paste("iteration:", i, "criteria:", round(eps, 5), "\n"))
        i <- i + 1
        if (i > 500) stop("max iter reached")
        coefs <- newcoefs
    }
    coefs
}

