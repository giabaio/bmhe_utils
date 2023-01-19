#' Computes and prints summary statistics for a vector or matrix of
#' simulated values
#'
#' @param x A vector or a matrix containing simulations from, eg, BUGS
#' @return tab A table with some specific summary statistics
#' @examples
#' x=rnorm(1000)
#' stats(x)
#'
stats <- function(x,dim=2){
  # This is a utility function that is used by 'stats' to produce the summaries
  bugs.stats <- function(x) {
    ## Used by the function stats
    c(mean(x),sd(x),quantile(x,.025),median(x),quantile(x,.975))
  }

  if(is.null(dim(x))==TRUE) {
	  tab <- bugs.stats(x)
	  names(tab) <- c("mean","sd","2.5%","median","97.5%")
  }
  ## for a matrix, it's a bit more complex --- need to specify the dimension along which want to take stats
  if(is.null(dim(x))==FALSE) {
  	tab <- t(apply(x,dim,function(x) bugs.stats(x)))
  	colnames(tab) <- c("mean","sd","2.5%","median","97.5%")
  }
  list(tab=tab)
}


#' Computes and prints summary statistics for a vector or matrix of
#' simulated values - tidyverse style
#'
#' @param x A vector or a matrix containing simulations from, eg, BUGS
#' @param digits The number of significant digits shown (default = 3)
#' @param na.rm A logical value (default TRUE) to indicate whether NA should be
#' removed
#' @examples
#' x=rnorm(1000)
#' stats2(x)
#'
# Tidyverse version of 'stats'
stats2 <- function(x,digits=3,na.rm=TRUE) {
 # Makes sure tidyverse is installed
 required_packages=c("tidyverse")
 for (pkg in required_packages) {
   if (!requireNamespace(pkg, quietly = TRUE)) {
     stop("`", pkg, "` is required: install.packages('", pkg, "')")
   }
   if (requireNamespace(pkg, quietly = TRUE)) {
     if (!is.element(pkg, (.packages()))) {
       suppressMessages(suppressWarnings(attachNamespace(pkg)))
     }
   }
 }
 # Specifies the number of significant digits
 options(pillar.sigfig=digits)
 # If a single vector needs to fiddle to keep the name
 if(is.null(dim(x))==TRUE) {
   nm=deparse(substitute(x))
   x %>% as_tibble() %>% mutate(variable=nm) %>%
     gather("variable", "value") %>% group_by(variable) %>%
     rename("Parameter"="variable") %>% summarise(
       mean = mean(value,na.rm=na.rm),
       sd = sd(value,na.rm=na.rm),
       `2.5%` = quantile(value, probs = .025,na.rm=na.rm),
       median = quantile(value, probs = .5,na.rm=na.rm),
       `97.5%` = quantile(value, probs = .975,na.rm=na.rm)
     )
 } else {
   # If a named matrix then can do this
   x %>% as_tibble() %>% gather("variable", "value") %>% group_by(variable) %>%
     rename("Parameter"="variable") %>% summarise(
       mean = mean(value,na.rm=na.rm),
       sd = sd(value,na.rm=na.rm),
       `2.5%` = quantile(value, probs = .025,na.rm=na.rm),
       median = quantile(value, probs = .5,na.rm=na.rm),
       `97.5%` = quantile(value, probs = .975,na.rm=na.rm)
     )
 }
 # As the result is a tibble, it can be further post-processed, eg to fix the number
 # of total digits (not the significant digits!) to, say, 5, can do
 # stats2(object) %>% mutate(across(where(is.numeric), num, digits=5))
 # Or to print the table, can use 'kable'/'kableExtra', eg
 # stats2(object) %>% knitr::kable(digits=4) %>% kableExtra::kable_styling()
}


#' Compute the parameters of a Beta distribution, given a prior guess for key
#' parameters. Based on "Bayesian ideas and data analysis", page 100.
#' Optimisation method to identify the values of a,b that give required
#' conditions on the Beta distribution
#'
#' @param mode The implied mode of the distribution
#' @param upp An upper bound value for the distribution
#' @param prob The estimated probability that theta <= upp
#' @return The list of relevant output including the values for the
#' parameters of the Beta distribution and some underlying summary statistics
#' of the resulting variable
#' @examples
#' res=betaPar2(.6,.7,.9)
#'
betaPar2 <- function(mode,upp,prob){

N <- 10000
b <- 1:N
a <- (1+mode*(b-2))/(1-mode)
sim <- qbeta(prob,a,b)
m <- ifelse(prob>=.5,max(which(sim>=upp)),min(which(sim>=upp)))
M <- ifelse(prob>=.5,min(which(sim<=upp)),max(which(sim<=upp)))

b <- min(m,M)+(b/N)
a <- (1+mode*(b-2))/(1-mode)
sim <- qbeta(prob,a,b)
m <- ifelse(prob>=.5,max(which(sim>=upp)),min(which(sim>=upp)))
M <- ifelse(prob>=.5,min(which(sim<=upp)),max(which(sim<=upp)))
a <- ifelse(m==M,a[m],mean(a[m],a[M]))
b <- ifelse(m==M,b[m],mean(b[m],b[M]))

step <- 0.001
theta <- seq(0,1,step)
density <- dbeta(theta,a,b)

norm.dens <- density/sum(density)
cdf <- cumsum(norm.dens)
M <- min(which(cdf>=.5))
m <- max(which(cdf<=.5))

theta.mode <- theta[which(density==max(density))]
theta.mean <- a/(a+b)
theta.median <- mean(theta[m],theta[M])
theta.sd <- sqrt((a*b)/(((a+b)^2)*(a+b+1)))
beta.params <- c(a,b,theta.mode,theta.mean,theta.median,theta.sd)
res1 <- beta.params[1]
res2 <- beta.params[2]
theta.mode <- beta.params[3]
theta.mean <- beta.params[4]
theta.median <- beta.params[5]
theta.sd <- beta.params[6]
list(
alpha=res1,beta=res2,theta.mode=theta.mode,theta.mean=theta.mean,theta.median=theta.median,theta.sd=theta.sd)
}


#' Computes the parameters of a Beta distribution so that the mean and
#' standard dev are the input (m,s)
#'
#' @param m The implied mean for the underlying Beta distribution
#' @param s The implied standard deviation for the underlying Beta distribution
#' @return The list of relevant output including the values for the
#' parameters of the Beta distribution (alpha and beta)
#' @examples
#' betaPar(.5,.15)
#'
betaPar <- function(m,s){
	a <- m*( (m*(1-m)/s^2) -1 )
	b <- (1-m)*( (m*(1-m)/s^2) -1 )
	list(alpha=a,beta=b)
}


#' Computes mean and variance of a logNormal distribution
#' so that the parameters on the natural scale are mu and sigma
#'
#' @param m The implied mean for the underlying Beta distribution
#' @param s The implied standard deviation for the underlying Beta distribution
#' @return The list of relevant output including the values for the
#' parameters of the logNormal distribution in terms of the mean on the log
#' scale (mulog) and the sd on the log scale (sigmalog)
#' @examples
#' lognPar(3,.15)
lognPar <- function(m,s) {
	s2 <- s^2
	mulog <- log(m) - .5*log(1+s2/m^2)
	s2log <- log(1+(s2/m^2))
	sigmalog <- sqrt(s2log)
	list(mulog=mulog,sigmalog=sigmalog)
}

#' Computes the parameters of a Gamma distribution so that the mean and
#' standard dev are the input (m,s)
#'
#' @param m The implied mean for the underlying Beta distribution
#' @param s The implied standard deviation for the underlying Beta distribution
#' @return The list of relevant output including the values for the
#' parameters of the Beta distribution (alpha and beta)
#' @examples
#' gammaPar(12,3)
#'
gammaPar <- function(m,s){
  b <- m/s^2
  a <- m*b
  list(shape=a,rate=b)
}


#' Makes a traceplot (eg to visualise MCMC simulations from multiple chains)
#'
#' @param node a *string* with the name of the node to be plotted, eg "theta"
#' (in quotes)
#' @param model the name of the object containing the MCMC simulations
#' @param title the title of the graph (defaults to nothing)
#' @param lab the label to write on the y-axis (defaults to nothing)
#' @return the graph with the traceplot
#' @examples
#' \dontrun{
#' }
#'
mytraceplot <- function(node,model=m,title="",lab=""){
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("`", pkg, "` is required: install.packages('", pkg, "')")
    }
    if (requireNamespace(pkg, quietly = TRUE)) {
      if (!is.element(pkg, (.packages()))) {
        suppressMessages(suppressWarnings(attachNamespace(pkg)))
      }
    }
  }
	xlab <- "Iteration"
## this way works with R2jags as well as with R2WinBUGS
	cmd <- ifelse(class(model)=="rjags",mdl <- model$BUGSoutput,mdl <- model)
	if (mdl$n.chains==1) {
		plot(mdl$sims.array[,1,node],t="l",col="blue",xlab=xlab,
			ylab=lab,main=title)
	}
	if (mdl$n.chains==2) {
		plot(mdl$sims.array[,1,node],t="l",col="blue",xlab=xlab,
			ylab=lab,main=title,ylim=range(mdl$sims.array[,1:2,node]))
		points(mdl$sims.array[,2,node],t="l",col="red")
	}
	if (mdl$n.chains>2) {
		plot(mdl$sims.array[,1,node],t="l",col="blue",xlab=xlab,
			ylab=lab,main=title,
			ylim=range(mdl$sims.array[,1:mdl$n.chains,node]))
		col <- c("red","green","magenta","yellow")
		for (i in 2:mdl$n.chains) {
			points(mdl$sims.array[,i,node],t="l",col=col[i])
		}
	}
}

#' Produces a plot of the values of the Gelman Rubin stats to determine
#' visually convergence (and see clearly which node has reached it)
#'
#' @param m is an object in the class jags or bugs (the output of the MCMC run)
#' @return the graph with the Gelmn Rubin statistics plot
#' @examples
#' \dontrun{
#' }
#'
plotGR <- function(m) {
cmd <- ifelse(class(m)=="rjags",mdl <- m$BUGSoutput,mdl <- m)
plot(mdl$summary[,"Rhat"],xlab="Saved parameters",ylab="Gelman-Rubin diagnostic",
	axes=F,col="white")
thresh <- 1.1
abline(h=thresh,col="dark grey",lwd=2)
text(1:dim(mdl$summary)[1],mdl$summary[,"Rhat"],
	row.names(mdl$summary),cex=.8)
axis(2)
box()
}


#' Computes the logit of a number
#'
#' @param x a number between 0 and 1
#' @return logit(x)=log(x/(1-x))
#' @examples
#' logit(.2)
#'
logit <- function(x){log(x/(1-x))}


#' Computes the inverse logit of a number between -infinity and +infinity
#'
#' @param x a real number
#' @return inverse-logit(x) = exp(x)/(1+exp(x))
#' @examples
#' ilogit(2)
#'
ilogit <- function(x){exp(x)/(1+exp(x))}


#' Maps from odds to probabilities
#'
#' @param odds the odds ratio *against* p: OR=(1-p)/p
#' @return the value of the underlying probability, p
#' @examples
#' odds2probs(4)
#'
odds2probs <- function(odds) {
  p <- 1/(odds+1)
  return(p)
}


#' Computes the odds ratio between two probabilities
#'
#' @param p1 a probability
#' @param p2 another probability
#' @return OR=(p1/(1-p1))/(p2/(1-p2))
#' @examples
#' OR(.5,.2)
#'
OR <- function(p1,p2) {
  OR <- (p1/(1-p1))/(p2/(1-p2))
  return(OR)
}
