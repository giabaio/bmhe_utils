#' Tidyverse based function to do traceplots
#'
#' Traceplot for a \code{bugs} or \code{jags} object
#'
#'
#' @param x an object of class `bugs', see \code{BUGS}, or of class
#' 'jags', see \code{JAGS} for details. It can also be a vector with
#' simulations for a single variable
#' @param parameter a string with the name of the parameter for which to show
#' the traceplot. Can be a vector, eg \code{c("par1","par2")}
#' @param ... further arguments to \code{\link{traceplot}}
#' @author Gianluca Baio
#' @seealso \code{BUGS}, \code{JAGS}
#' @export traceplot
traceplot=function(x,parameter=NULL,...) {
  # Makes sure tidyverse is installed
  required_packages=c("dplyr","ggplot2")
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

  # If the input is a `BUGS`/`JAGS` object than use the relevant simulations
  if(inherits(x,what=c("rjags","bugs"))) {
    # If the object is in the class JAGS,then selects the relevant list
    if(any(grepl("rjags",class(x)))) {x=x$BUGSoutput}
    #
    if(is.null(parameter)) {
      title="Traceplot for all model parameters"
    } else {
      title=paste0("Traceplot for selected model parameters")
    }
    p=x$sims.array %>%
      as_tibble(.name_repair = ~paste0("Chain",x$sims.array %>% as_tibble() %>% colnames())) %>%
      # This is only run if 'parameter' is not null (so some parameters are selected by the user)
      { if(!is.null(parameter)) select(., contains(parameter)) else . } %>%
      mutate(iteration=row_number()) %>%
      tidyr::gather(variable,value,c(-iteration)) %>%
      tidyr::separate(variable,c("chain","parameter"),extra = "merge") %>%
      ggplot(aes(x=iteration,y=value,color=chain))+
      geom_line()+facet_wrap(~parameter,scales="free")+
      labs(title=title)+ theme_bw()
  }

  # If the input is a vector it will still work
  if(is.logical(x)) {x=as.numeric(x)}
  if(is.numeric(x) || is.integer(x)) {
    p=x |> as_tibble() |> ggplot(aes(1:length(x),x)) +
      geom_line() + theme_bw() + xlab("Iterations") + ylab("")
  }
  p
}


#' Various plots for the posteriors in a \code{bugs} or \code{jags} object
#'
#'
#' @param x an object of class `bugs', see \code{BUGS}, or of class
#' 'jags', see \code{JAGS} for details
#' @param parameter a string with the name of the parameter for which to show
#' the density plot. Can be a vector, eg \code{c("par1","par2")}
#' @param plot the type of plot (options are 'density' (default) or
#' 'bar' for a binned barplot of the posterior) or
#' 'hist' for a histogram
#' @param add_deviance a logical argument to determine whether the `deviance`
#' should be added to the plot (in case it is monitored). Defaults to `FALSE`
#' @param ... further arguments
#' @author Gianluca Baio
#' @seealso \code{BUGS}, \code{JAGS}
#' @export posteriorplot
posteriorplot=function(x,parameter=NULL,plot="density",add_deviance=FALSE,...) {
  # Makes sure tidyverse is installed
  required_packages=c("dplyr","ggplot2")
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
      { if(grepl("deviance",x$sims.list %>% names()) %>% any() && !add_deviance) select(.,-deviance) else . } %>%
      { if(!is.null(parameter)) select(.,contains(parameter)) else . } %>%
      tidyr::gather(variable,value,1:ncol(.)) %>%
      ggplot(aes(value))+geom_density()+facet_wrap(~variable,scales="free") +
      theme_bw()
  }
  if(plot=="bar") {
    p=x$sims.matrix %>% as_tibble() %>%
      { if(grepl("deviance",x$sims.list %>% names()) %>% any() && !add_deviance) select(.,-deviance) else . } %>%
      { if(!is.null(parameter)) select(.,contains(parameter)) else . } %>%
      tidyr::gather(variable,value,1:ncol(.)) %>%
      ggplot(aes(value))+geom_bar()+scale_x_binned() + facet_wrap(~variable,scales="free") +
      theme_bw()
  }
  if(plot=="hist") {
    p=x$sims.matrix %>% as_tibble() %>%
      { if(grepl("deviance",x$sims.list %>% names()) %>% any() && !add_deviance) select(.,-deviance) else . } %>%
      { if(!is.null(parameter)) select(.,contains(parameter)) else . } %>%
      tidyr::gather(variable,value,1:ncol(.)) %>%
      ggplot(aes(value))+geom_histogram()+ facet_wrap(~variable,scales="free") +theme_bw()
  }
  p
}


#' Specialised diagnostic plots
#'
#' Creates a plot showing the output of convergence indicators, such as
#' the Potential Scale Reduction and the effective sample size
#'
#' @param x an object of class `bugs', see \code{BUGS}, or of class
#' 'jags', see \code{JAGS} for details
#' @param what A string indicating what diagnostic measure should be plotted.
#' Options are 'Rhat' (default), indicating the PSR statistic, or 'n.eff',
#' indicating the effective sample size
#' @param label A logical input. If set to 'FALSE' (default), then does not
#' include text labels next to each node
#' @param ...  Additional options
#' @author Gianluca Baio
#' @seealso \code{BUGS}, \code{JAGS}
#' @keywords Diagnostic plots
#' @examples
#' \dontrun{
#' diagplot(m)
#' }
#' @export diagplot
#'
diagplot=function(x,what="Rhat",label=FALSE,...) {

  required_packages=c("dplyr","ggplot2")
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
    pl=x$summary %>% as_tibble() %>% ggplot(aes(1:nrow(.),!!sym(what))) +
      geom_point(color="red",size=2) +
      geom_hline(
        yintercept=ifelse(what=="Rhat",1.1,x$n.sims),linetype="dashed",size=.5
      ) +
      theme_bw() +
      labs(
        x="Parameters",
        title=ifelse(
          what=="Rhat","Potential scale reduction","Effective sample size"
        )
      )
  }
  # May add labels to the plot
  if (label) {
    pl=pl + geom_text(
      aes(
        label=
          as.character(x$summary |> rownames())
    ),vjust=-1)
  }
  # Now produces the plot
  pl
}

#' Coefplot for the parameters in the model
#'
#' Creates a plot showing the mean and an interval estimate for the posterior
#' distributions in a given model.
#'
#' @param x an object of class `bugs', see \code{BUGS}, or of class
#' 'jags', see \code{JAGS} for details
#' @param low the lower quantile to consider (default 2.5 percentile)
#' @param upp the upper quantile to consider (default 97.5 percentile)
#' @param parameter a vector of strings with the names of the parameters to be
#' included. Defaults to all those in the original model, but can be a
#' vector eg \code{c("par1","par2")}
#' @param ...  Additional options
#' @author Gianluca Baio
#' @seealso \code{BUGS}, \code{JAGS}
#' @keywords Coefficient plots
#' @examples
#' \dontrun{
#' coefplot(m)
#' }
#' @export coefplot
#'
coefplot=function(x,low=.025,upp=.975,parameter=NULL,...) {

  required_packages=c("dplyr","ggplot2")
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
  xmin=sym(paste0(low*100,"%"))
  xmax=sym(paste0(upp*100,"%"))
  toplot=x$sims.matrix %>% as_tibble() %>%
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
    select(Parameter,everything()) |> mutate(num=row_number())
  # This ensures the correct ordering of variables
  toplot |> ggplot(aes(mean,num)) +
    geom_linerange(aes(xmin=!!xmin,xmax=!!xmax),position=position_dodge(.3)) +
    geom_point(position = position_dodge(0.3)) + theme_bw() + geom_vline(xintercept=0,linetype="dashed") +
    labs(x="Interval estimate",title="Coefplot",y="Parameter") +
    scale_y_continuous(breaks=1:nrow(toplot),labels=toplot$Parameter)
}

#' Trial-and-error Beta plot
#'
#' Provides a quick and dirty, trial-and-error tool to identify suitable values for the the
#' parameters of a Beta distribution to match set properties (eg mean, sd, 95% interval)
#'
#' @param a_max The maximum value for the parameter \code{a} of the Beta distribution
#' @param b_max The maximum value for the parameter \code{b} of the Beta distribution
#' @param step The increment in the grid of values for \code{a} and \code{b}
#' @author Gianluca Baio
#' @keywords Beta distribution
#' @examples
#' \dontrun{
#' betaplot()
#' }
#' @export betaplot
#'
betaplot=function(a_max=30,b_max=30,step=.01) {
  # Checks and loads the necessary packages
  required_packages=c("manipulate")
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
  # Needs to create a "fake" plot (using base-R) so that 'manipulate' works OK with 'ggplot'
  a=b=1 # Initialise a,b to avoid issue with no-visible binding
  manipulate::manipulate(plot(a,b,col="white",axes=F,xlab="",ylab=""),a=manipulate::picker(1),b=manipulate::picker(1))

  # Utility functions to get stats for the resulting distribution
  mbeta=function(a,b){a/(a+b)}
  sdbeta=function(a,b){sqrt((a*b)/((a+b)^2*(a+b+1)))}
  intbeta=function(a,b){c(stats::qbeta(.025,a,b),stats::qbeta(.975,a,b))}

  # Create the actual plot using 'ggplot'
  manipulate::manipulate(
    # Creates the plot
    ggplot()+stat_function(fun=stats::dbeta,args=list(shape1=a,shape2=b))+
      xlim(0,1)+scale_y_continuous(expand=expansion(mult=c(0,.11)))+theme_bw()+
      annotate(
        "text",x=0,y=Inf,
        label=paste0(
          "a=",a,"\nb=",b,"\nmean=",format(mbeta(a,b),digits=3,nsmall=3),
          "\nsd=",format(sdbeta(a,b),digits=3,nsmall=3),
          "\n95% interval=[",format(intbeta(a,b)[1],digits=3,nsmall=3)," - ",format(intbeta(a,b)[2],digits=3,nsmall=3),"]"
        ),
        vjust=1.5,hjust=0,size=3
      )+geom_segment(aes(x=intbeta(a,b)[1],xend=intbeta(a,b)[2],y=0,yend=0),size=1.5),
    # Customise the values for the parameters
    a=manipulate::slider(0,a_max,step=step,label="Value for a",initial=0.01),
    b=manipulate::slider(0,b_max,step=step,label="Value for b",initial=0.01)
  )
}


#' Trial-and-error Gamma plot
#'
#' Provides a quick and dirty, trial-and-error tool to identify suitable values for the the
#' parameters of a Gamma distribution to match set properties (eg mean, sd, 95% interval)
#'
#' @param shape_max The maximum value for the parameter \code{shape} of the Gamma distribution
#' @param rate_max The maximum value for the parameter \code{rate} of the Gamma distribution
#' @param step The increment in the grid of values for \code{shape} and \code{rate}
#' @author Gianluca Baio
#' @keywords Gamma distribution
#' @examples
#' \dontrun{
#' gammaplot()
#' }
#' @export gammaplot
#'
gammaplot=function(shape_max=30,rate_max=30,step=.01) {
  # Checks and loads the necessary packages
  required_packages=c("manipulate")
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
  # Needs to create a "fake" plot (using base-R) so that 'manipulate' works OK with 'ggplot'
  manipulate::manipulate(plot(a,b,col="white",axes=F,xlab="",ylab=""),a=picker(1),b=picker(1))

  # Utility functions to get stats for the resulting distribution
  mgamma=function(shape,rate){shape/rate}
  sdgamma=function(shape,rate){sqrt(shape/rate^2)}
  intgamma=function(shape,rate){c(qgamma(.025,shape=shape,rate=rate),qgamma(.975,shape=shape,rate=rate))}

  # Create the actual plot using 'ggplot'
  manipulate::manipulate(
    # Creates the plot
    ggplot()+stat_function(fun=dgamma,args=list(shape=shape,rate=rate))+
      xlim(0,30)+theme_bw()+scale_y_continuous(expand=expansion(mult=c(0,.11)))+ #c(0, 0), limits = c(0, NA))+
      annotate(
        "text",x=0,y=Inf,
        label=paste0(
          "shape=",shape,"\nrate=",rate,"\nmean=",format(mgamma(shape,rate),digits=3,nsmall=3),
          "\nsd=",format(sdgamma(shape,rate),digits=3,nsmall=3),
          "\n95% interval=[",format(intgamma(shape,rate)[1],digits=3,nsmall=3)," - ",format(intgamma(shape,rate)[2],digits=3,nsmall=3),"]"
        ),
        vjust=1.5,hjust=0,size=3
      )+geom_segment(aes(x=intgamma(shape,rate)[1],xend=intgamma(shape,rate)[2],y=0,yend=0),size=1.5),
    # Customise the values for the parameters
    shape=manipulate::slider(0,shape_max,step=step,label="Value for a",initial=0.01),
    rate=manipulate::slider(0,rate_max,step=step,label="Value for b",initial=0.01)
  )
}


#' Autocorrelation plot
#'
#' Plots the ACF function
#'
#' @param x A vector with simulations from a MCMC process (eg from a \code{BUGS}
#' or \code{JAGS} run)
#' @param col The color with which to plot the ACF (default to \code{"black"})
#' @param parameter A text string to select a named parameter (eg if using
#' a \code{BUGS} or \code{JAGS} object, that would be one of the monitored
#' parameters)
#' @param add_deviance a logical argument to determine whether the `deviance`
#' should be added to the plot (in case it is monitored). Defaults to `TRUE`,but
#' is only relevant if the input object `x` is a \code{BUGS} or \code{JAGS}
#' object
#' @param ... Extra arguments
#' @author Gianluca Baio
#' @keywords Autocorrelation function
#' @examples
#' \dontrun{
#' acfplot(m)
#' }
#' @export acfplot
#'
acfplot=function(x,col="black",parameter=NULL,add_deviance=TRUE,...) {
  # Needs to add options to customise
  # a. calling the `stats::acf` function
  # b. the `ggplot` graph
  #
  required_packages=c("dplyr","ggplot2")
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
  if(class(x) %in% c("rjags","bugs")) {
    if(any(grepl("rjags",class(x)))) {
      x=x$BUGSoutput
    }
    ac=x$sims.matrix %>% as_tibble() %>%
      { if(grepl("deviance",x$sims.list %>% names()) %>% any() && !add_deviance) select(.,-deviance) else . } %>%
      { if(!is.null(parameter)) select(.,contains(parameter)) else . } %>%
      lapply(.,acf,plot=FALSE)

    tab=ac[[1]] |> purrr::map_df(~.) |> mutate(series=names(ac)[1])
    if(length(ac)>1) {
      for (i in 2:length(ac)) {
        tab=tab |> bind_rows(ac[[i]] |> purrr::map_df(~.) |> mutate(series=names(ac)[i]))
      }
    }
    p=tab |>
      ggplot(aes(lag,acf))+geom_hline(aes(yintercept = 0)) +
      geom_segment(mapping = aes(xend = lag, yend = 0),linewidth=1,col=col) +
      geom_hline(aes(
        yintercept=qnorm((1 + (1 - 0.05))/2)/sqrt(n.used)
      ), linetype = 2, color = 'blue') +
      geom_hline(aes(
        yintercept=-qnorm((1 + (1 - 0.05))/2)/sqrt(n.used)
      ), linetype = 2, color = 'blue') + xlab("Lag") + ylab("ACF")+
      facet_wrap(~series,scales="free") +
      theme_bw()
  }

  # If the object is a vector, then plot the ACF for that single variable
  if(is.numeric(x)) {
    ac=acf(x,plot=F)
    p=tibble(x=ac$lag,y=ac$acf) |>
      ggplot(aes(x,y))+geom_hline(aes(yintercept = 0)) +
      geom_segment(mapping = aes(xend = x, yend = 0),linewidth=1,col=col) + theme_bw() +
      geom_hline(aes(
        yintercept=qnorm((1 + (1 - 0.05))/2)/sqrt(ac$n.used)
      ), linetype = 2, color = 'blue') +
      geom_hline(aes(
        yintercept=-qnorm((1 + (1 - 0.05))/2)/sqrt(ac$n.used)
      ), linetype = 2, color = 'blue') + xlab("Lag") + ylab("ACF")
  }
  p
}
