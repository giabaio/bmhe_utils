#' Tidyverse based function to do traceplots
#'
#' Traceplot for a \code{bugs} or \code{jags} object
#'
#'
#' @param x an object of class `bugs', see \code{\link{bugs}}, or of class
#' 'jags', see \code{\link{jags}} for details
#' @param parameter a string with the name of the parameter for which to show
#' the traceplot
#' @param ... further arguments to \code{\link{traceplot}}
#' @author Gianluca Baio
#' @seealso \code{\link{bugs}}, \code{\link{jags}}
#' @export traceplot
traceplot=function(x,parameter=NULL,...) {
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
  # If the object is in the class JAGS,then selects the relevant list
  if(class(x)=="JAGS") {x=x$BUGSoutput}
  #
  if(is.null(parameter)) {
    x$sims.array %>%
      as_tibble(.name_repair = ~paste0("Chain",x$sims.array %>% as_tibble() %>% colnames())) %>%
      mutate(iteration=row_number()) %>%
      gather(variable,value,c(-iteration)) %>%
      separate(variable,c("chain","parameter"),extra = "merge") %>%
      ggplot(aes(x=iteration,y=value,color=chain))+
      geom_line()+facet_wrap(~parameter,scales="free")+
      labs(title="Traceplot for all model parameters")+
      theme_bw()
  } else {
    x$sims.array[,,parameter] %>%
      as_tibble(.name_repair = ~paste("Chain",1:x$n.chains)) %>% mutate(iteration=row_number()) %>%
      gather(chain,value,contains("Chain"),-iteration) %>% ggplot(aes(x=iteration,y=value,color=chain))+geom_line() +
      theme_bw() + labs(title=paste("Traceplot for",parameter))
  }
}


#' Various plots for the posteriors in a \code{bugs} or \code{jags} object
#'
#'
#' @param x an object of class `bugs', see \code{\link{bugs}}, or of class
#' 'jags', see \code{\link{jags}} for details
#' @param parameter a string with the name of the parameter for which to show
#' the density plot
#' @param plot the type of plot (options are 'density' (default) or
#' 'hist' for a binned barplot of the posterior)
#' @param ... further arguments to \code{\link{densityplot}}
#' @author Gianluca Baio
#' @seealso \code{\link{bugs}}, \code{\link{jags}}
#' @export posteriorplot
posteriorplot=function(x,parameter=NULL,plot="density",...) {
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
  # If the object is in the class JAGS,then selects the relevant list
  if(class(x)=="JAGS") {x=x$BUGSoutput}
  #

  if(is.null(parameter)) {
    if(plot=="density") {
      p=x$sims.matrix %>% as_tibble() %>% gather(variable,value,1:ncol(.)) %>%
        ggplot(aes(value))+geom_density()+facet_wrap(~variable,scales="free") +
        theme_bw()
    }
    if(plot=="hist") {
      p=x$sims.matrix %>% as_tibble() %>% gather(variable,value,1:ncol(.)) %>%
        ggplot(aes(value))+geom_bar()+scale_x_binned() + facet_wrap(~variable,scales="free") +
        theme_bw()
    }
    # This would do a barplot of all the variables
    #
  } else {
    if(plot=="density") {
      p=x$sims.matrix[,parameter] %>% as_tibble() %>% ggplot(aes(value))+geom_density() +
        labs(x=parameter,title=paste("Density plot for",parameter)) + theme_bw()
    }
    if (plot=="hist") {
      p=x$sims.matrix[,parameter] %>% as_tibble() %>% ggplot(aes(value))+geom_bar() +
        scale_x_binned()+labs(x=parameter,title=paste("Density plot for",parameter)) +
        theme_bw()
    }
  }
  p
}


#' Specialised diagnostic plots
#'
#' Creates a plot showing the output of convergence indicators, such as
#' the Potential Scale Reduction and the effective sample size
#'
#' @param x an object of class `bugs', see \code{\link{bugs}}, or of class
#' 'jags', see \code{\link{jags}} for details
#' @param what A string indicating what diagnostic measure should be plotted.
#' Options are 'Rhat' (default), indicating the PSR statistic, or 'n.eff',
#' indicating the effective sample size
#' @param ...  Additional options
#' @author Gianluca Baio
#' @seealso \code{\link{bugs}}, \code{\link{jags}}
#' @keywords Diagnostic plots
#' @examples
#' \dontrun{
#' }
#' @export diagplot
#'
diagplot=function(x,what="Rhat",...) {

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
  # If the object is in the class JAGS,then selects the relevant list
  if(class(x)=="JAGS") {x=x$BUGSoutput}
  #

  x$summary %>% as_tibble() %>% ggplot(aes(1:nrow(.),!!sym(what))) +
    geom_point(color="red",size=2) + geom_hline(yintercept=ifelse(what=="Rhat",1.1,x$n.sims),linetype="dashed",size=.5) +
    theme_bw() + labs(x="Parameters",title=ifelse(what=="Rhat","Potential scale reduction","Effective sample size"))
}

#' Coefplot for the parameters in the model
#'
#' Creates a plot showing the mean and an interval estimate for the posterior
#' distributions in a given model.
#'
#' @param x an object of class `bugs', see \code{\link{bugs}}, or of class
#' 'jags', see \code{\link{jags}} for details
#' @param low the lower quantile to consider (default 2.5 percentile)
#' @param upp the upper quantile to consider (default 97.5 percentile)
#' @param params a vector of strings with the names of the parameters to be
#' included. Defaults to all those in the original model
#' @param deviance a logical value (defaults to FALSE) to indicate whether
#' the model deviance should be considered in the plot.
#' @param ...  Additional options
#' @author Gianluca Baio
#' @seealso \code{\link{bugs}}, \code{\link{jags}}
#' @keywords Diagnostic plots
#' @examples
#' \dontrun{
#' }
#' @export coefplot
#'
coefplot=function(x,low=.025,upp=.975,params=NULL,deviance=FALSE,...) {

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
  # If the object is in the class JAGS,then selects the relevant list
  if(class(x)=="JAGS") {x=x$BUGSoutput}
  #

  if(is.null(params)) {
    params=x$sims.matrix %>% colnames()
  }
  if((any(grepl("deviance",x$sims.matrix %>% colnames()))) & (deviance==FALSE)) {
    params=params[-grep("deviance",x$sims.matrix %>% colnames())]
  }
  x$sims.matrix %>% apply(2,function(x) c(mean(x,na.rm=T),sd(x,na.rm=T),quantile(x,low,na.rm=T),quantile(x,upp,na.rm=T))) %>%
    t() %>% as_tibble(.name_repair=~c("mean","sd",paste0(low*100,"%"),paste0(upp*100,"%"))) %>%
    mutate(Parameter=x$sims.matrix %>% colnames()) %>% select(Parameter,everything()) %>%
    filter(Parameter%in%params) %>%
    ggplot(aes(mean,Parameter))+
    geom_linerange(aes(xmin=`2.5%`,xmax=`97.5%`),position=position_dodge(.3)) +
    geom_point(position = position_dodge(0.3)) + theme_bw() + geom_vline(xintercept=0,linetype="dashed") +
    labs(x="Interval estimate",title="Coefplot")
}
