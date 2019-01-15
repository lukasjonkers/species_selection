# modified from rioja:::randomPTF
# performs MAT instead of WA
# no normalistaion of relative abundances during bootstrapping

randomMAT <- function (spec, env, ncol = 1, nVar, nTF, verbose = TRUE, do.parallel = FALSE, ...) 
{
	    if (do.parallel) {
        if (!requireNamespace("foreach", quietly = TRUE)) 
            stop("Package foreach is required for parallel operation. Please install and try again.")
    }
	
    do.TF <- function(y, x, nsp, nsam, fun, ncol, nVar, ...) {
        res2 <- vector("numeric", length = nsp) # create vector for result
        sel.sp <- sample.int(nsp, nVar)  # randomly select 1/3 of species
        boot <- sample.int(nsam, replace = TRUE) # select samples for boot
        test <- setdiff(1:nsam, boot) # which sample does not occur in boot
        x.b <- x[boot] # select env[boot]
        y.b <- y[boot, sel.sp] # select species in boot and sel.sp
        x.t <- x[test] # select env[!boot]: out of bag
        y.t <- y[test, sel.sp] # select species[!boot, sel.sp]: out of bag
        sel <- apply(y.b, 2, sum) > 0.001 # test if summed species occur in subsample
        y.b <- y.b[, sel]
        y.t <- y.t[, sel]
        sel.sp <- sel.sp[sel]
        sel <- apply(y.b, 1, sum) > 0.01
        y.b <- y.b[sel, ]
        x.b <- x.b[sel]
       	mod <- MAT(y.b, x.b) # build TF model
        pred.oob <- predict(mod, y.t, verbose = FALSE)$fit[,1] # get mean prediction of oob
        mseOOB <- mean((pred.oob - x.t)^2) # get error of oob prediction
        nv <- ncol(y.t)
        nsam.t <- nrow(y.t)
        for (j in 1:nv) {
            y.t2 <- y.t
            y.t2[, j] <- y.t[sample.int(nsam.t), j]
            pred.oob_h <- predict(mod, y.t2, verbose = FALSE)$fit[,1]
            mseOOB_H <- mean((pred.oob_h - x.t)^2)
            res2[sel.sp[j]] <- mseOOB_H - mseOOB
        }
        res2
    }
    y <- as.matrix(spec)
    x <- as.matrix(env)
    nsam <- nrow(y)
    nsp <- ncol(spec)
    if (missing(nVar)) 
        nVar <- max(as.integer(nsp/3), 1)
	if (do.parallel) {
        .paropts <- NULL
        .paropts$.combine <- function(...) NULL
        i <- seq_len(nTF)
        fe_call <- as.call(c(list(quote(foreach::foreach), i = i), 
            .paropts))
        fe <- eval(fe_call)
        `%dopar%` <- foreach::`%dopar%`
        res <- foreach::foreach(1:nTF, .packages = c("rioja")) %dopar% 
            {
                do.TF(y, x, nsp, nsam, fun, ncol, nVar, ...)
            }
        res <- t(sapply(res, "["))
    }
    else {
        if (verbose) {
            writeLines("Running:")
            pb <- txtProgressBar(min = 0, max = 1, style = 3)
            on.exit(close(pb))
        }
        res <- matrix(ncol = nsp, nrow = nTF)
        for (i in 1:nTF) {
            if (verbose) {
                setTxtProgressBar(pb, i/nTF)
            }
            res[i, ] <- do.TF(y, x, nsp, nsam, fun, ncol, nVar, 
                ...)
        }
    }
    #print(class(res))
    colnames(res) <- colnames(y)
    SumErr <- colSums(res, na.rm = TRUE)
    nTree <- colSums(!is.na(res)) 
    VI <- SumErr/nTree # mean difference
    names(VI) <- colnames(y)
    VI <- data.frame(VI = sort(VI, decreasing = TRUE))
    res <- list(VI = VI, spec = y, env = x)
    class(res) <- "randomPTF"
    res
}
