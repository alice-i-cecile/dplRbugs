# Libraries ####
library(dplR)
library(ggplot2)
library(reshape2)

# Reproducting the bug ####

# Load bugged series
ngxp_bugged <- read.rwl("ngxp_bugged.rwl")

# Attempt to detrend, showing plots
std_series <- detrend.series(ngxp_bugged[,1], make.plot=TRUE)

# Reconstruct values of the detrending curve
# Raw data / detrending curve = standardized
# Raw data /  standardized = detrending curve
detrending_curve <- ngxp_bugged / std_series$ModNegExp

# Show raw data and detrending curve in a clearer graph
plot_df <- melt(data.frame(Raw=ngxp_bugged[[1]], Detrending=detrending_curve[[1]], Year=as.numeric(rownames(ngxp_bugged))), id.vars="Year")
ggplot(data=plot_df,aes(x=Year, y=value, colour=variable)) + geom_line(size=2) + theme_bw() + geom_hline(x=0) + ylab("Ring Width (mm)") + scale_colour_discrete(name="Curve")


# Fixing detrend.series ####

# Opening up the relevant sections of detrend.series
# Relevant line of code for model-fitting
# nec <- nls(formula = Y ~ a * exp(b * seq_along(Y)) + k, start = list(a = a, b = b, k = k))

# k and a must be greater than 0 
# to ensure that the predicted values used in the detrending curve stay positive as well
# But no constraint is given! Let's add some

# Specified lower constraints
# Also set b < 0 to ensure correct shape (manually checked for just below in original code)
# Had to switch algorithm to "port" to use constraints
# lower=c(a=0, b=-Inf, k=0), upper=c(a=Inf, b=0, k=Inf), algorithm="port"

# Redefining detrend.series to limit allowable range of values
detrend.series2 <- function (y, y.name = "", make.plot = TRUE, method = c("Spline",                                                        "ModNegExp", "Mean"), nyrs = NULL, f = 0.5, pos.slope = FALSE) 
{
    known.methods <- c("Spline", "ModNegExp", "Mean")
    method2 <- match.arg(arg = method, choices = known.methods, 
                         several.ok = TRUE)
    good.y <- which(!is.na(y))
    if (length(good.y) == 0) {
        stop("all values are 'NA'")
    }
    else if (any(diff(good.y) != 1)) {
        stop("'NA's are not allowed in the middle of the series")
    }
    y2 <- y[good.y]
    y2[y2 == 0] <- 0.001
    resids <- list()
    if ("ModNegExp" %in% method2) {
        nec.func <- function(Y) {
            a <- mean(Y[seq_len(floor(length(Y) * 0.1))])
            b <- -0.01
            k <- mean(Y[floor(length(Y) * 0.9):length(Y)])
            nec <- nls(formula = Y ~ a * exp(b * seq_along(Y)) + 
                       k, start = list(a = a, b = b, k = k),
                       lower=c(a=0, b=-Inf, k=0), upper=c(a=Inf, b=0, k=Inf),
                       algorithm="port")
            if (coef(nec)[2] >= 0) 
                stop()
            fits <- predict(nec)
            if (fits[1] < fits[length(fits)]) 
                stop()
            if (fits[length(fits)] < 0) 
                stop()
            fits
        }
        ModNegExp <- try(nec.func(y2), silent = TRUE)
        if (class(ModNegExp) == "try-error") {
            tm <- seq_along(y2)
            lm1 <- lm(y2 ~ tm)
            ModNegExp <- predict(lm1)
            if (coef(lm1)[2] > 0 && !pos.slope) 
                ModNegExp <- rep(mean(y2), length(y2))
        }
        resids$ModNegExp <- y2/ModNegExp
        do.mne <- TRUE
    }
    else {
        do.mne <- FALSE
    }
    if ("Spline" %in% method2) {
        if (is.null(nyrs)) 
            nyrs2 <- floor(length(y2) * 0.67)
        else nyrs2 <- nyrs
        Spline <- ffcsaps(y = y2, x = seq_along(y2), nyrs = nyrs2, 
                          f = f)
        resids$Spline <- y2/Spline
        do.spline <- TRUE
    }
    else {
        do.spline <- FALSE
    }
    if ("Mean" %in% method2) {
        Mean <- rep(mean(y2), length(y2))
        resids$Mean <- y2/Mean
        do.mean <- TRUE
    }
    else {
        do.mean <- FALSE
    }
    resids <- data.frame(resids)
    if (make.plot) {
        op <- par(no.readonly = TRUE)
        on.exit(par(op))
        par(mar = c(2.5, 2.5, 2.5, 0.5) + 0.1, mgp = c(1.5, 0.5, 
                                                       0))
        n.rows <- 1 + ncol(resids)
        mat <- matrix(seq_len(n.rows), n.rows, 1)
        layout(mat, widths = rep(0.5, ncol(mat)), heights = rep(1, 
                                                                nrow(mat)))
        plot(y2, type = "l", ylab = "mm", xlab = gettext("Age (Yrs)", 
                                                         domain = "R-dplR"), main = gettextf("Raw Series %s", 
                                                                                             y.name, domain = "R-dplR"))
        if (do.spline) 
            lines(Spline, col = "green", lwd = 2)
        if (do.mne) 
            lines(ModNegExp, col = "red", lwd = 2)
        if (do.mean) 
            lines(Mean, col = "blue", lwd = 2)
        if (do.spline) {
            plot(resids$Spline, type = "l", col = "green", main = gettext("Spline", 
                                                                          domain = "R-dplR"), xlab = gettext("Age (Yrs)", 
                                                                                                             domain = "R-dplR"), ylab = gettext("RWI", domain = "R-dplR"))
            abline(h = 1)
        }
        if (do.mne) {
            plot(resids$ModNegExp, type = "l", col = "red", main = gettext("Neg. Exp. Curve or Straight Line", 
                                                                           domain = "R-dplR"), xlab = gettext("Age (Yrs)", 
                                                                                                              domain = "R-dplR"), ylab = gettext("RWI", domain = "R-dplR"))
            abline(h = 1)
        }
        if (do.mean) {
            plot(resids$Mean, type = "l", col = "blue", main = gettext("Horizontal Line (Mean)", 
                                                                       domain = "R-dplR"), xlab = gettext("Age (Yrs)", 
                                                                                                          domain = "R-dplR"), ylab = gettext("RWI", domain = "R-dplR"))
            abline(h = 1)
        }
    }
    resids2 <- matrix(NA, ncol = ncol(resids), nrow = length(y))
    resids2 <- data.frame(resids2)
    names(resids2) <- names(resids)
    if (!is.null(names(y))) 
        row.names(resids2) <- names(y)
    resids2[good.y, ] <- resids
    resids2 <- resids2[, method2]
    if (!is.data.frame(resids2)) 
        names(resids2) <- names(y)
    resids2
}

# Verifying bug fix ####
# Attempt to detrend, showing plots
std_series2 <- detrend.series2(ngxp_bugged[,1], make.plot=TRUE)

# Reconstruct values of the detrending curve
# Raw data / detrending curve = standardized
# Raw data /  standardized = detrending curve
detrending_curve2 <- ngxp_bugged / std_series2$ModNegExp

# Show raw data and detrending curve in a clearer graph
plot_df2 <- melt(data.frame(Raw=ngxp_bugged[[1]], Detrending=detrending_curve2[[1]], Year=as.numeric(rownames(ngxp_bugged))), id.vars="Year")
ggplot(data=plot_df2,aes(x=Year, y=value, colour=variable)) + geom_line(size=2) + theme_bw() + geom_hline(x=0) + ylab("Ring Width (mm)") + scale_colour_discrete(name="Curve")

# Success!