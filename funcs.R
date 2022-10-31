#############################
## Useful custom functions ##
#############################

#Shade function - from rethinking package
shade <- function (object, lim, label = NULL, 
          col = alpha("black", 0.15), 
          border = NA, ...) 
{
  if (missing(lim)) 
    stop("Interval limits missing.")
  if (missing(object)) 
    stop("No density or formula object.")
  from <- lim[1]
  to <- lim[2]
  if (class(object)[1] == "formula") {
    x1 <- eval(object[[3]])
    y1 <- eval(object[[2]])
    x <- x1[x1 >= from & x1 <= to]
    y <- y1[x1 >= from & x1 <= to]
  }
  if (class(object)[1] == "density") {
    x <- object$x[object$x >= from & object$x <= to]
    y <- object$y[object$x >= from & object$x <= to]
  }
  if (class(object)[1] == "matrix" & length(dim(object)) == 
      2) {
    y <- c(object[1, ], object[2, ][ncol(object):1])
    x <- c(lim, lim[length(lim):1])
  }
  if (class(object)[1] == "matrix") {
    polygon(x, y, col = col, border = border, ...)
  }
  else {
    polygon(c(x, to, from), c(y, 0, 0), col = col, border = border, 
            ...)
  }
  if (!is.null(label)) {
    lx <- mean(x)
    ly <- max(y)/2
    text(lx, ly, label)
  }
}

#Round function
round2 = function(x, n=0) {
  posneg = sign(x)
  z = abs(x)*10^n
  z = z + 0.5 + sqrt(.Machine$double.eps)
  z = trunc(z)
  z = z/10^n
  z*posneg
}

#Estimate k from prevalence and mean worm burden
k_estim <- function(p,m){
  numer = -m*log(1-p)
  denom1 = log(1-p)
  denom2 = (1-p)^(1/m)*log(1-p)/m
  lambert = lambertWm1(denom2)
  k = numer/(denom1-m*lambert)
  return(k)
}

#examine model output
getDist <- function(x){
  m = mean(x)
  int95 = quantile(x, probs=c(0.025, 0.975))
  return(c(m, int95))
}

#Power Law function
PL_func <- function(x, y1=53, gamma=0.96) y1*x^gamma

#Zero inflated negative binomial sensitivity
nbh_sens<- function(x, b) x/(x+b)

#observed prevalecne given true prevalence, sensivity and specificity
obs_prev_func <- function(prev, se, sp=1) prev*se + (1-prev)*(1-sp)

#True prevalence based on M and k from negative binomial
prev <- function(x, k) 1-pnbinom(q=0, mu=x, size=k)

#Population level sensitivity
pop_sensitivity <- function(x, k, b){
  max_worms = 1:35000
  pop_prob_dist = dnbinom(max_worms, mu=x, size=k) #Matrix with worm prob dist given NB population params
  sens_x = nbh_sens(max_worms, b=b) #Individual level sensitivity given worm burden x
  pop_sens = sum(pop_prob_dist*sens_x)/(1-dnbinom(0, mu=x, size=k))
  return(pop_sens)
}