#' Tidyverse based function to do traceplots
#'
#' Traceplot for a \code{bugs} or \code{jags} object
#'
#'
#' @param x an object of class `bugs', see \code{\link{bugs}}, or of class
#' 'jags', see \code{\link{jags}} for details
#' @param parameter a string with the name of the parameter for which to show
#' the traceplot. Can be a vector, eg \code{c("par1","par2")}
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
  if(any(grepl("rjags",class(x)))) {x=x$BUGSoutput}
  #
  x$sims.array %>%
    as_tibble(.name_repair = ~paste0("Chain",x$sims.array %>% as_tibble() %>% colnames())) %>%
    # This is only run if 'parameter' is not null (so some parameters are selected by the user)
    { if(!is.null(parameter)) select(., contains(parameter)) else . } %>%
    mutate(iteration=row_number()) %>%
    gather(variable,value,c(-iteration)) %>%
    separate(variable,c("chain","parameter"),extra = "merge") %>%
    ggplot(aes(x=iteration,y=value,color=chain))+
    geom_line()+facet_wrap(~parameter,scales="free")+
    labs(title="Traceplot for all model parameters")+
    theme_bw()
}


#' Various plots for the posteriors in a \code{bugs} or \code{jags} object
#'
#'
#' @param x an object of class `bugs', see \code{\link{bugs}}, or of class
#' 'jags', see \code{\link{jags}} for details
#' @param parameter a string with the name of the parameter for which to show
#' the density plot. Can be a vector, eg \code{c("par1","par2")}
#' @param plot the type of plot (options are 'density' (default) or
#' 'bar' for a binned barplot of the posterior) or
#' 'hist' for a histogram
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
  if(any(grepl("rjags",class(x)))) {x=x$BUGSoutput}
  #

  if(plot=="density") {
    p=x$sims.matrix %>% as_tibble() %>%
      { if(!is.null(parameter)) select(.,contains(parameter)) else . } %>%
      gather(variable,value,1:ncol(.)) %>%
      ggplot(aes(value))+geom_density()+facet_wrap(~variable,scales="free") +
      theme_bw()
  }
  if(plot=="bar") {
    p=x$sims.matrix %>% as_tibble() %>%
      { if(!is.null(parameter)) select(.,contains(parameter)) else . } %>%
      gather(variable,value,1:ncol(.)) %>%
      ggplot(aes(value))+geom_bar()+scale_x_binned() + facet_wrap(~variable,scales="free") +
      theme_bw()
  }
  if(plot=="hist") {
    p=x$sims.matrix %>% as_tibble() %>%
      { if(!is.null(parameter)) select(.,contains(parameter)) else . } %>%
      gather(variable,value,1:ncol(.)) %>%
      ggplot(aes(value))+geom_histogram()+ facet_wrap(~variable,scales="free") +theme_bw()
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
  if(any(grepl("rjags",class(x)))) {x=x$BUGSoutput}
  #
  if(x$n.chains==1) {
    cat("You need to run 2 or more parallel chains to be able to make this plot")
  } else {
    x$summary %>% as_tibble() %>% ggplot(aes(1:nrow(.),!!sym(what))) +
      geom_point(color="red",size=2) + geom_hline(yintercept=ifelse(what=="Rhat",1.1,x$n.sims),linetype="dashed",size=.5) +
      theme_bw() + labs(x="Parameters",title=ifelse(what=="Rhat","Potential scale reduction","Effective sample size"))
  }
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
#' @param parameter a vector of strings with the names of the parameters to be
#' included. Defaults to all those in the original model, but can be a
#' vector eg \code{c("par1","par2")}
#' @param ...  Additional options
#' @author Gianluca Baio
#' @seealso \code{\link{bugs}}, \code{\link{jags}}
#' @keywords Coefficient plots
#' @examples
#' \dontrun{
#' }
#' @export coefplot
#'
coefplot=function(x,low=.025,upp=.975,params=NULL,...) {

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
  if(any(grepl("rjags",class(x)))) {x=x$BUGSoutput}
  #

  x$sims.matrix %>% as_tibble() %>%
    { if(grepl("deviance",x$sims.list %>% names()) %>% any()) select(.,-deviance) else . } %>%
    { if(!is.null(parameter)) select(.,contains(parameter)) else . } %>%
    apply(2,function(x) c(mean(x,na.rm=T),sd(x,na.rm=T),quantile(x,low,na.rm=T),quantile(x,upp,na.rm=T))) %>%
    t() %>% as_tibble(.name_repair=~c("mean","sd",paste0(low*100,"%"),paste0(upp*100,"%"))) %>%
    mutate(
      Parameter=x$sims.matrix %>% as_tibble() %>%
        { if(grepl("deviance",x$sims.list %>% names()) %>% any()) select(.,-deviance) else . } %>%
        { if(!is.null(parameter)) select(.,contains(parameter)) else . } %>%
        colnames()
    ) %>%
    select(Parameter,everything()) %>%
    ggplot(aes(mean,Parameter))+
    geom_linerange(aes(xmin=`2.5%`,xmax=`97.5%`),position=position_dodge(.3)) +
    geom_point(position = position_dodge(0.3)) + theme_bw() + geom_vline(xintercept=0,linetype="dashed") +
    labs(x="Interval estimate",title="Coefplot")
}
