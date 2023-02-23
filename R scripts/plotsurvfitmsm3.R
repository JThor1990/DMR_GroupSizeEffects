# I have added a very minor tweak to the function within msm in order to plot a shaded polygon
# and to invert the plot to plot the probability of disappearance

plot.survfit.msm3 <- function (x, from = 1, to = NULL, range = NULL, covariates = "mean", 
                               interp = c("start", "midpoint"), 
                               ci = c("none", "normal", "bootstrap"), B = 100, 
                               legend.pos = NULL, 
                               xlab = "Time", ylab = "Probability of disappearance", 
                               lty = 1, lwd = 1, col = "red", lty.ci = 2, lwd.ci = 1, 
                               col.ci = "red", mark.time = TRUE, col.surv = "blue", 
                               lty.surv = 2, lwd.surv = 1, survdata = FALSE, 
                               col.poly = "red", alpha.poly = 0.2, ...) # I have added these two arguments (see line 52)
{
  require(survival)
  
  if (!inherits(x, "msm")) 
    stop("expected \"x\" to be a msm model")
  if (is.null(to)) 
    to <- max(absorbing.msm(x))
  else {
    if (!is.numeric(to)) 
      stop("\"to\" must be numeric")
    if (!(to %in% absorbing.msm(x))) 
      stop("\"to\" must be an absorbing state")
  }
  if (!(from %in% transient.msm(x))) 
    stop("\"from\" must be a non-absorbing state")
  if (is.null(range)) 
    rg <- range(model.extract(x$data$mf, "time"))
  else {
    if (!is.numeric(range) || length(range) != 2) 
      stop("\"range\" must be a numeric vector of two elements")
    rg <- range
  }
  interp <- match.arg(interp)
  ci <- match.arg(ci)
  timediff <- (rg[2] - rg[1])/50
  times <- seq(rg[1], rg[2], timediff)
  pr <- lower <- upper <- numeric()
  for (t in times) {
    P <- pmatrix.msm(x, t, t1 = times[1], covariates = covariates, 
                     ci = ci, B = B)
    if (ci != "none") {
      pr <- c(pr, P$estimates[from, to])
      lower <- c(lower, P$L[from, to])
      upper <- c(upper, P$U[from, to])
    }
    else pr <- c(pr, P[from, to])
  }
  plot(times, pr, type = "l", xlab = xlab, ylab = ylab, 
       lwd = lwd, ylim = c(0, 1), lty = lty, col = col, ...)
  if (ci != "none") {
    lines(times, lower, lty = lty.ci, col = col.ci, lwd = lwd.ci)
    lines(times, upper, lty = lty.ci, col = col.ci, lwd = lwd.ci)
    polygon(c(times, rev(times)), c(lower, rev(upper)), col = adjustcolor(col.poly, alpha.f = alpha.poly), border=NA) # added
  }
  dat <- x$data$mf[, c("(subject)", "(time)", "(state)")]
  dat$"(subject)" <- match(dat$"(subject)", unique(dat$"(subject)"))
  dat$subjstate <- paste(dat$"(subject)", dat$"(state)")
  anyfrom <- tapply(dat$"(state)", dat$"(subject)", 
                    function(x) any(x == from))[as.character(dat$"(subject)")]
  dat <- dat[anyfrom, ]
  obspt <- sequence(table(dat$"(subject)"))
  minfrom <- rep(which(dat$"(state)" == from & !duplicated(dat$subjstate)) - 
                   which(!duplicated(dat$"(subject)")), table(dat$"(subject)")) +   1
  dat <- dat[obspt >= minfrom, ]
  first <- !duplicated(dat$"(subject)")
  last <- !duplicated(dat$"(subject)", fromLast = TRUE)
  anyabs <- tapply(dat$"(state)", dat$"(subject)", 
                   function(x) any(x == to))[as.character(dat$"(subject)")]
  subjstate <- paste(dat$"(subject)", dat$"(state)")
  minabs <- dat$"(state)" == to & !duplicated(subjstate)
  dtime <- dat$"(time)" - tapply(dat$"(time)", 
                                 dat$"(subject)", min)[as.character(dat$"(subject)")]
  if (interp == "midpoint") {
    prevminabs <- c(minabs[-1], FALSE)
    dtime[minabs] <- 0.5 * (dtime[minabs] + dtime[prevminabs])
  }
  minabs[!anyabs] <- last[!anyabs]
  survdat <- data.frame(survtime = dtime[minabs], died = as.numeric(anyabs[minabs]))
  lines(survival::survfit(Surv(survdat$survtime, survdat$died) ~ 1), 
        mark.time = mark.time, col = col.surv, lty = lty.surv, fun = "event",
        lwd = lwd.surv)
  timediff <- (rg[2] - rg[1])/50
  if (!is.numeric(legend.pos) || length(legend.pos) != 2) 
    legend.pos <- c(max(x$data$mf$"(time)") - 25 * 
                      timediff, 1)
  if (ci == "none") 
    legend(legend.pos[1], legend.pos[2], lty = c(lty, lty.surv), 
           lwd = c(lwd, lwd.surv), col = c(col, col.surv), legend = c("Fitted", "Empirical"))
  else legend(legend.pos[1], legend.pos[2], lty = c(lty, lty.ci, lty.surv), 
              lwd = c(lwd, lwd.ci, lwd.surv), 
              col = c(col, col.ci, col.surv), 
              legend = c("Fitted", "Fitted (confidence interval)",  "Empirical"))
  if (survdata) 
    survdat
  else invisible()
}
