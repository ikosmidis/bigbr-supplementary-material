bigglm <- function(formula, data, family = gaussian(), ...) {
    UseMethod("bigglm", data)
}
setGeneric("bigglm", signature = c("formula", "data"))


bigglm.data.frame <- function(formula, data, ..., chunksize = 5000) {
    n <- nrow(data)
    cursor <- 0
    datafun <- function(reset = FALSE) {
        if (reset) {
            cursor <<- 0
            return(NULL)
        }
        if (cursor >= n) {
            return(NULL)
        }
        start <- cursor + 1
        cursor <<- cursor + min(chunksize, n - cursor)
        data[start:cursor, ]
    }
    rval <- bigglm(formula = formula, data = datafun, ...)
    rval$call <- sys.call()
    rval$call[[1]] <- as.name(.Generic)
    rval
}

## bigglm.function with a bias reduction option (only for family =
## binomial() and family = poisson()). If br = TRUE then two passes
## are required through the database per IWLS iteration: one to get
## the new ML iterate from the previous one, and R from the QR
## decomposition of W^{1/2} X (and relevant quantities, like deviance,
## etc), and one to get the first-order bias at the rpevious ML
## iterate using the current value of R. This implementats the
## algorithms in Kosmidis et al (2021, StCo)
bigglm.function <- function(formula, data, family = gaussian(),
                            weights = NULL, sandwich = FALSE,
                            maxit = 8, tolerance = 1e-7, start = NULL,
                            quiet = FALSE,
                            type = c("ML", "BRASE", "MJPL"),
                            implementation = c("1pass", "2pass"),
                            verbose = FALSE, ...) {

    ## From glm to allow for the various ways family may be specified
    if (is.character(family))
        family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family))
        family <- family()
    if (is.null(family$family)) {
        print(family)
        stop("'family' not recognized")
    }

    type <- match.arg(type)
    br <- isTRUE(type == "BRASE") ## bias-reducing adjusted score equations
    mjpl <- isTRUE(type == "MJPL") ## maximum Jefrreys-penalized likelihood
    adj <- br | mjpl
    implementation <- match.arg(implementation)
    one_pass <- implementation == "1pass" & adj
    two_pass <- implementation == "2pass" & adj

    ## Enrich family
    if (adj) {
        if (br & !(family$family %in% c("binomial", "poisson"))) {
            ## Because we need
            stop("Reduced-bias estimation is currently implemented only for binomial and poisson responses. If your data is of moderate size and in-memory consider using the `brglm_fit` method from the brglm2 R package.")
        }
        linkglm <- make.link(family$link)
        linkglm <- enrichwith::enrich(linkglm, with = c("d2mu.deta"))
        ## Put everything into the family object
        family[names(linkglm)] <- linkglm
        family <- enrichwith::enrich(family, with = c("d1variance"))
    }

    tt <- terms(formula)
    beta <- betaprev <- start

    ## Starting values are set to rep(0, nrow(x)) here if start = NULL
    etafun <- function(x) if (is.null(beta)) rep(0, nrow(x)) else drop(x %*% beta)
    etafun_prev <- function(x) if (is.null(betaprev)) rep(0, nrow(x)) else drop(x %*% betaprev)

    get_hats <- function(R1inv, wow) {
        xRi <- mm %*% R1inv
        rowSums(xRi * xRi * wow)
    }

    converged <- FALSE
    for (i in 1:maxit) {

        ## BEGIN: IWLS for ML iterate
        firstchunk <- TRUE
        deviance <- 0
        rss <- 0
        data(reset = TRUE)
        n <- 0
        sh <- 0
        while (!is.null(chunk <- data(reset = FALSE))) {
            n <- n + nrow(chunk)
            mf <- model.frame(tt, chunk)
            mm <- model.matrix(tt, mf)
            p <- NCOL(mm)
            if (!is.null(weights)) {
                if (!inherits(weights, "formula")) {
                    stop("`weights' must be a formula")
                }
                w <- model.frame(weights, chunk)[[1]]
            } else {
                w <- rep(1, nrow(mm))
            }
            if (firstchunk) {
                qr <- bigqr.init(p)
                assn <- attr(mm, "assign")
                if (sandwich) {
                    xyqr <- bigqr.init(p * (p + 1))
                }
            }
            if (!identical(assn, attr(mm, "assign"))) {
                stop("model matrices incompatible")
            }
            y <- model.response(mf)
            if (is.null(off <- model.offset(mf))) off <- 0
            eta <- etafun(mm) + off
            mu <- family$linkinv(eta)
            dmu <- family$mu.eta(eta)
            ww <- w * dmu * dmu / (family$variance(mu))
            z <- eta + (y - mu) / dmu
            if (one_pass & i > 1) {
                eta_prev <- etafun_prev(mm)
                mu_prev <- family$linkinv(eta_prev)
                dmu_prev <- family$mu.eta(eta_prev)
                ww_prev <- w * dmu_prev * dmu_prev / (family$variance(mu_prev))
                h <- get_hats(R1inv, ww_prev)
                xi <- h * family$d2mu.deta(eta) / (2 * dmu * ww)
                z <- z + xi
                if (mjpl) {
                    vd <- family$d1variance(mu_prev)
                    z <- z + xi - h * vd / (2 * w * dmu)
                }
                sh <- sh + sum(h)
            }
            qr <- update(qr, mm, z - off, ww)
            if (!is.null(beta)) {
                deviance <- deviance + sum(family$dev.resids(y, mu, w))
                rss <- rss + sum((y - mu)^2 / (w * family$variance(mu))) * (sum(w) / length(w))
                if (sandwich) {
                    xx <- matrix(nrow = nrow(mm), ncol = p * (p + 1))
                    xx[, 1:p] <- mm * (drop(z) - off)
                    for (j in 1:p) {
                        xx[, p * j + (1:p)] <- mm * mm[, j]
                    }
                    if (verbose > 0) cat("Updating QR for sandwich, n = ", n, "...")
                    xyqr <- update(xyqr, xx, rep(0, nrow(mm)), ww * ww)
                    if (verbose > 0) cat("Done!\n")
                }
            }
            firstchunk <- FALSE
        }
        iwlm <- list(
            call = sys.call(-1), qr = qr, iterations = i,
            assign = attr(mm, "assign"), terms = tt, converged = FALSE,
            n = n, names = colnames(mm), weights = weights, rss = rss
        )
        if (sandwich) {
            iwlm$sandwich <- list(xy = xyqr)
        }
        class(iwlm) <- "biglm"
        ## END: IWLS for ML iterate

        ## BEGIN: IWLS for first-order bias iterate
        ##
        ## The checks for early termination in the ML loop are not
        ## necessary for the bias loop because bigglm.function would
        ## have already failed at this stage
        ## bias <- 0
        if (adj) {
            R1inv <- backsolve(getR(qr), diag(p))
        }
        if (two_pass) {
            firstchunk <- TRUE
            data(reset = TRUE)
            while (!is.null(chunk <- data(reset = FALSE))) {
                mf <- model.frame(tt, chunk)
                mm <- model.matrix(tt, mf)
                p <- NCOL(mm)
                if (!is.null(weights)) {
                    w <- model.frame(weights, chunk)[[1]]
                } else {
                    w <- rep(1, nrow(mm))
                }
                if (firstchunk) {
                    qr <- bigqr.init(p)
                }
                if (is.null(off <- model.offset(mf))) off <- 0
                eta <- etafun(mm) + off
                mu <- family$linkinv(eta)
                dmu <- family$mu.eta(eta)
                ww <- w * dmu * dmu / (family$variance(mu))
                h <- get_hats(R1inv, ww)
                xi <- h * family$d2mu.deta(eta) / (2 * dmu * ww)
                if (mjpl) {
                    vd <- family$d1variance(mu)
                    xi <- 2 * xi - h * vd / (2 * w * dmu)
                }
                qr <- update(qr, mm, xi, ww)
                sh <- sh + sum(h)
                firstchunk <- FALSE
            }
            ## bias <- -coef(qr)
            ## Adjust thetab in iwlm by thetab in iwlm
            iwlm$qr$thetab <- iwlm$qr$thetab + qr$thetab
        }
        ## END: IWLS for first-order bias iterate

        betaprev <- beta
        beta <- coef(iwlm)

        if (isTRUE(i %% verbose == 0)) {
            delta <- betaprev - beta
            cat("Type:", type, "\n")
            if (adj) {
                cat("Implementation:", if (one_pass) "1pass\n" else "2pass\n")
            }
            cat("Iteration:", paste0(sprintf(paste0("%0", nchar(maxit), "d"),  i)), "/", maxit, ":\n")
            cat("||beta||_infinity =", max(abs(beta)), "\n")
            if (length(delta)) {
                linf <- max(abs(delta))
                l1 <- sum(abs(delta))
                l2 <- sqrt(sum(delta^2))
            } else {
                linf <- l1 <- l2 <- NA
            }
            cat("||beta - beta_prev||_infinity =", linf, "\n")
            cat("||beta - beta_prev||_1 =", l1, "\n")
            cat("||beta - beta_prev||_2 =", l2, "\n")
            ## cat("sum(hats) =", if (adj) sh else p, "\n")
        }

        if (i >= maxit) {
            if (!quiet) warning("ran out of iterations and failed to converge")
            break
        }

        if (is.null(betaprev)) {
            delta <- NA
        } else {
            ## delta <- (betaprev - beta) / sqrt(diag(vcov(iwlm)))
            delta <- betaprev - beta
            ## Use Linf norm
            if (max(abs(delta)) < tolerance) {
                iwlm$converged <- TRUE
                break
            }
        }
    }

    ## iwlm is a good object to use even if adj = TRUE since
    ## quantities are evaluated at the beta iterate just before the
    ## the beta value that signals convergence
    rval <- iwlm
    rval$type <- type
    rval$delta <- delta
    rval$impementation <- if (adj) implementation else "1pass"
    rval$family <- family
    rval$deviance <- deviance
    rval$df.resid <- rval$n - length(rval$qr$D)
    class(rval) <- c("bigglm", "biglm")
    rval
}



print.bigglm <- function(x, ...) {
    cat("Large data regression model: ")
    print(x$call)
    cat("Sample size = ", x$n, "\n")
    if (is.null(x$converged) || !x$converged) {
        cat("failed to converge after", x$iterations, "iterations\n")
    }
    invisible(x)
}
