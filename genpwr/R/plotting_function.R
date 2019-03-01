#########################################################################################################
#Plotting functions to summarize
#########################################################################################################

#########################################################################################################
#Plots for Sample Size
#########################################################################################################

#' Function to convert to numeric with scientific notation contaning the "." character
#'
#' convert to numeric with scientific notation contaning the "." character
#'
#'
#' @param char string to be converted to numeric
#'
#' @return a number
#'
#' @export
#'
as.numeric2<-function(char){
  if(grepl("e.", char)) char <- gsub("e.", "e-", char)
  return(as.numeric(char))
}


#' Function to Plot Sample Size Results
#'
#' Plot the sample size results by MAF, OR, Alpha or Power
#'
#' @import ggplot2
#'
#' @param data The data frame result from \code{\link{ss.calc}}
#' @param x The desired variable on the y axis: "MAF", "OR", "ES","Alpha", or "Power"
#' @param panel.by A grouping variable to panel the graphs by: "True.Model", "MAF", "OR", "Alpha", or "Power"
#' @param y_limit An object specifying the minimum and maximum of the y-axis (eg c(0,4)) default is NULL, which allows the limits to be picked automatically
#' @param y_log  Logical, specifying whether the y axis should be logarithmic. Default is F
#' @param return_gg Logical, specifying whether to return the ggplot object instead of printing out the plot
#' @param linear.effect.measure Should the graphs indicate ES values, or R2 values? (default ES)
#' @param select.Alpha Only produce graphs for the specified Alpha level(s).
#' @param select.OR Only produce graphs for the specified odds ratio(s).
#' @param select.ES Only produce graphs for the specified effect size(s).
#' @param select.Power Only produce graphs for the specified power(s).
#' @param select.MAF Only produce graphs for the specified minor allele frequency(ies).
#' @param select.Case.Rate Only produce graphs for the specified case rate(s).
#' @param select.SD Only produce graphs for the specified standard deviation(s).
#' @param select.True.Model Only produce graphs for the specified true genetic model(s): "Additive", "Dominant", "Recessive".
#' @param select.Test.Model Only produce graphs for the specified testing model(s): "Additive", "Dominant", "Recessive", "2df".
#'
#' @return A series of plots with sample size on the Y axis.
#'
#' @examples
#' ss <- ss.calc(power=0.8, Case.Rate=c(0.5), k=NULL,
#'     MAF=seq(0.01, 0.05, 0.01), OR=c(3,4),Alpha=c(0.05),
#'     True.Model='All', Test.Model='All')
#' ss.plot(data=ss, x='MAF',panel.by='OR')
#'
#' @export
#'
ss.plot<-function(data=NULL,x='MAF', panel.by='True.Model', y_limit = NULL, y_log = F, return_gg = F,
                  linear.effect.measure = 'ES',
                  select.Alpha = NULL,
                  select.OR = NULL,
                  select.ES = NULL,
                  select.Power = NULL,
                  select.MAF = NULL,
                  select.Case.Rate = NULL,
                  select.SD = NULL,
                  select.True.Model = NULL,
                  select.Test.Model = NULL)
{

  #library(ggplot2)
  "%ni%" <- Negate("%in%")
  
  nstr <- "N_total_at_Alpha_"
  if(any(grepl("SampleSize_GE_at_Alpha_", colnames(data)))) nstr <- "SampleSize_GE_at_Alpha_"
 
  #Create Dataset with separate rows for each Alpha level  
  indices<-grep(nstr,colnames(data))
  Alphas<-as.numeric2(unlist( strsplit(colnames(data)[indices],nstr))[2])
  
  ss.new<-NULL
  for (i in 1:length(Alphas)){ss.new<-rbind(ss.new, data.frame(data[-indices], N_total=data[,indices[i]], Alpha=Alphas[i]))}

  #Subset data to only include those in the select. statements
  if(is.null(select.Alpha)==F){ss.new<-ss.new[ss.new$Alpha %in% select.Alpha,]}
  if(is.null(select.OR)==F){ss.new<-ss.new[ss.new$OR %in% select.OR,]}
  if(is.null(select.ES)==F){ss.new<-ss.new[ss.new$ES %in% select.ES,]}
  if(is.null(select.Power)==F){ss.new<-ss.new[ss.new$Power %in% select.Power,]}
  if(is.null(select.MAF)==F){ss.new<-ss.new[ss.new$MAF %in% select.MAF,]}
  if(is.null(select.Case.Rate)==F){ss.new<-ss.new[ss.new$Case.Rate %in% select.Case.Rate,]}
  if(is.null(select.SD)==F){ss.new<-ss.new[ss.new$SD %in% select.SD,]}
  if(is.null(select.True.Model)==F){ss.new<-ss.new[ss.new$True.Model %in% select.True.Model,]}
  if(is.null(select.Test.Model)==F){ss.new<-ss.new[ss.new$Test.Model %in% select.Test.Model,]}


  # var<-c("MAF", 'OR','ES','Power', 'Case.Rate', 'SD','Alpha', 'True.Model')
  var<-c("MAF", "OR", "OR_G", "OR_E", "OR_GE", "P_e", "sd_e", "ES", "ES_G", "ES_E", "ES_GE", "R2", "R2_E", "R2_G", "R2_GE",
      "Power", "Case.Rate", "SD_Y","Alpha", "True.Model")
  var<-var[!(var %in% c(x, panel.by))]
  var<-var[var %in% colnames(ss.new)]
  
  if(linear.effect.measure %ni% c("ES", "R2")) stop("'linear.effect.measure should be either 'ES' or 'R2'")
  if(linear.effect.measure=="ES"){var <- var[var %ni% c("R2", "R2_G", "R2_E", "R2_GE")]}
  if(linear.effect.measure=="R2"){var <- var[var %ni% c("ES", "ES_G", "ES_E", "ES_GE")]}

  if(return_gg) resplots <- list()

  graphs<-unique.data.frame(ss.new[,var])
  # graphs$Test.Model<-graphs$
  for(j in 1:nrow(graphs))
  {
    subtitle<-paste0("N by ", x, ": ",
                  paste0(sapply(1:ncol(graphs), function(ii) 
                    paste0(var[ii], "=", graphs[j,ii])), collapse =  ", "))
    temp1 <- ss.new[apply(sapply(1:ncol(graphs), function(ii) ss.new[,var[ii]]==graphs[j,ii]), 1, all),]

    temp2<-temp1[order(temp1$True.Model),]
    temp2[,panel.by]<-paste(panel.by,"=", temp2[,panel.by])

    plot_obj <- ggplot(data=temp2, aes(x=temp2[,x], y=N_total, group = Test.Model, colour = Test.Model)) +
            geom_line() +
            geom_point( size=4, shape=21, fill="white") + xlab(x)+ggtitle(subtitle)+
            facet_wrap( ~temp2[,panel.by])
    if(!is.null(y_limit)){
      plot_obj <- plot_obj + ylim(y_limit)
    }
    if(y_log){
      plot_obj <- plot_obj + scale_y_log10()
    }
    if(return_gg){ 
      resplots <- c(resplots, list(plot_obj))
      return(resplots)
    }else print(plot_obj)
  }

}


#########################################################################################################
#Plots for Power
#########################################################################################################
#' Function to Plot Power Results
#'
#' Plot the power results by MAF, OR, Alpha or N
#'
#' @import ggplot2
#'
#' @param data The data frame result from \code{\link{power.calc}}
#' @param x The desired variable on the y axis: "MAF", "OR", "Alpha", or "N_total"
#' @param panel.by A grouping variable to panel the graphs by: "True.Model", "MAF", "OR", "Alpha", or "N_total"
#' @param y_limit An object specifying the minimum and maximum of the y-axis (eg c(0,4)) default is NULL, which allows the limits to be picked automatically
#' @param y_log  Logical, specifying whether the y axis should be logarithmic. Default is F
#' @param return_gg Logical, specifying whether to return the ggplot object instead of printing out the plot
#' @param linear.effect.measure Should the graphs indicate ES values, or R2 values? (default ES)
#' @param select.Alpha Only produce graphs for the specified Alpha level(s).
#' @param select.OR Only produce graphs for the specified odds ratio(s).
#' @param select.ES Only produce graphs for the specified effect sizes(s).
#' @param select.N Only produce graphs for the specified sample size(s).
#' @param select.MAF Only produce graphs for the specified minor allele frequency(ies).
#' @param select.Case.Rate Only produce graphs for the specified case rate(s).
#' @param select.SD Only produce graphs for the specified standard deviation(s).
#' @param select.True.Model Only produce graphs for the specified true genetic model(s): "Additive", "Dominant", "Recessive".
#' @param select.Test.Model Only produce graphs for the specified testing model(s): "Additive", "Dominant", "Recessive", "2df".
#'
#' @return A series of plots with power on the Y axis.
#'
#' @examples
#' pw <- power.calc(N=c(400, 600), Case.Rate=c(0.4), k=NULL,
#'    MAF=seq(0.05, 0.1, 0.01), OR=c(3),Alpha=c(0.05),
#'    True.Model=c("Dominant", "Recessive", "Additive"), 
#'    Test.Model=c("Dominant", "Recessive", "Additive", "2df"))
#'
#' @export
#'
power.plot<-function(data=NULL,x='MAF', panel.by='True.Model', y_limit = NULL, y_log = F, return_gg = F,
                     linear.effect.measure = 'ES',
                     select.Alpha = NULL,
                     select.OR = NULL,
                     select.ES = NULL,
                     select.N = NULL,
                     select.MAF = NULL,
                     select.Case.Rate = NULL,
                     select.SD = NULL,
                     select.True.Model = NULL,
                     select.Test.Model = NULL)
{

  "%ni%" <- Negate("%in%")

  #Create Dataset with separate rows for each Alpha level
  powerstr <- "Power_at_Alpha_"
  if(any(grepl("Power_GE_at_Alpha_", colnames(data)))) powerstr <- "Power_GE_at_Alpha_"
  indices<-grep(powerstr,colnames(data))
  Alphas<-as.numeric2(unlist( strsplit(colnames(data)[indices],powerstr))[2])
    
 
  ss.new<-NULL
  for (i in 1:length(Alphas)){ss.new<-rbind(ss.new, data.frame(data[-indices], Power=data[,indices[i]], Alpha=Alphas[i]))}

  ss.new$Test.Model <- as.character(ss.new$Test.Model)
  #Subset data to only include those in the select. statements
  if(!is.null(select.Alpha)){ss.new<-ss.new[ss.new$Alpha %in% select.Alpha,]}
  if(!is.null(select.OR)){ss.new<-ss.new[ss.new$OR %in% select.OR,]}
  if(!is.null(select.ES)){ss.new<-ss.new[ss.new$ES %in% select.ES,]}
  if(!is.null(select.N)){ss.new<-ss.new[ss.new$N_total %in% select.N,]}
  if(!is.null(select.MAF)){ss.new<-ss.new[ss.new$MAF %in% select.MAF,]}
  if(!is.null(select.Case.Rate)){ss.new<-ss.new[ss.new$Case.Rate %in% select.Case.Rate,]}
  if(!is.null(select.SD)){ss.new<-ss.new[ss.new$SD_Y %in% select.SD,]}
  if(!is.null(select.True.Model)){ss.new<-ss.new[ss.new$True.Model %in% select.True.Model,]}
  if(!is.null(select.Test.Model)){ss.new<-ss.new[ss.new$Test.Model %in% select.Test.Model,]}


  var<-c("MAF", "OR", "OR_G", "OR_E", "OR_GE", "P_e", "sd_e", "ES", "ES_G", "ES_E", "ES_GE", "R2", "R2_E", "R2_G", "R2_GE",
      "N_total", "Case.Rate", "SD_Y","Alpha", "True.Model")
  if(linear.effect.measure %ni% c("ES", "R2")) stop("'linear.effect.measure should be either 'ES' or 'R2'")
  if(linear.effect.measure=="ES"){var <- var[var %ni% c("R2", "R2_G", "R2_E", "R2_GE")]}
  if(linear.effect.measure=="R2"){var <- var[var %ni% c("ES", "ES_G", "ES_E", "ES_GE")]}
  var<-var[!(var %in% c(x, panel.by))]
  var<-var[var %in% colnames(ss.new)]
  graphs<-unique.data.frame(ss.new[,var])

  if(return_gg) resplots <- list()

  for(j in 1:nrow(graphs)){
    subtitle<-paste0("Power by ", x, ": ",
                  paste0(sapply(1:ncol(graphs), function(ii) 
                    paste0(var[ii], "=", graphs[j,ii])), collapse =  ", "))
    if(length(unlist(strsplit(subtitle, split = ", "))) > 4){
      # split up the title if it's too long for one line
      st_vec0 <- unlist(strsplit(subtitle, split = ", "))
      num_items <- 5
      n_locs <- seq(from = num_items, to = length(st_vec0) + ceiling(length(st_vec0)/(num_items - 1)) - 1, by = num_items)
      st_vec <- character(length(st_vec0) + length(n_locs))
      qq = 1
      for(ii in 1:length(st_vec)){
        if(ii %in% n_locs){st_vec[ii] <- "\n"
        }else{st_vec[ii] <- st_vec0[qq]; qq = qq + 1}
      }
      subtitle <- paste0(st_vec, collapse = ", ")
      subtitle <- gsub("\n, ", "\n", subtitle)
    }

    # temp1<-ss.new[ss.new[,var[1]]==graphs[j,1] & ss.new[,var[2]]==graphs[j,2]
    #               &ss.new[,var[3]]==graphs[j,3] & ss.new[,var[4]]==graphs[j,4],]
    temp1 <- ss.new[apply(sapply(1:ncol(graphs), function(ii) ss.new[,var[ii]]==graphs[j,ii]), 1, all),]

    temp2<-temp1[order(temp1$True.Model),]
    temp2[,panel.by]<-paste(panel.by,"=", temp2[,panel.by])

    plot_obj <- ggplot(data=temp2, aes(x=temp2[,x], y=Power, group = Test.Model, colour = Test.Model)) +
      geom_line() +
      geom_point( size=4, shape=21, fill="white") + xlab(x)+ggtitle(subtitle)+
      facet_wrap( ~temp2[,panel.by])
    if(!is.null(y_limit)){
      plot_obj <- plot_obj + ylim(y_limit)
    }
    if(y_log){
      plot_obj <- plot_obj + scale_y_log10()
    }
    if(return_gg){ 
      resplots <- c(resplots, list(plot_obj))
      return(resplots)
    }else print(plot_obj)
  }
}


#########################################################################################################
#Plots for Odds Ratio
#########################################################################################################
#' Function to Plot Odds Ratio Results
#'
#' Plot the power results by MAF, Power, Alpha or N
#'
#' @import ggplot2
#'
#' @param data The data frame result from \code{\link{power.calc}}
#' @param x The desired variable on the y axis: "MAF", "OR", "Alpha", or "N_total"
#' @param panel.by A grouping variable to panel the graphs by: "True.Model", "MAF", "Power", "Alpha", or "N_total"
#' @param y_limit An object specifying the minimum and maximum of the y-axis (eg c(0,4)) default is NULL, which allows the limits to be picked automatically
#' @param y_log  Logical, specifying whether the y axis should be logarithmic. Default is F
#' @param return_gg Logical, specifying whether to return the ggplot object instead of printing out the plot
#' @param linear.effect.measure Should the graphs indicate ES values, or R2 values? (default ES)
#' @param select.Alpha Only produce graphs for the specified Alpha level(s).
#' @param select.power Only produce graphs for the specified Power(s).
#' @param select.ES Only produce graphs for the specified effect sizes(s).
#' @param select.N Only produce graphs for the specified sample size(s).
#' @param select.MAF Only produce graphs for the specified minor allele frequency(ies).
#' @param select.Case.Rate Only produce graphs for the specified case rate(s).
#' @param select.SD Only produce graphs for the specified standard deviation(s).
#' @param select.True.Model Only produce graphs for the specified true genetic model(s): "Additive", "Dominant", "Recessive".
#' @param select.Test.Model Only produce graphs for the specified testing model(s): "Additive", "Dominant", "Recessive", "2df".
#'
#' @return A series of plots with power on the Y axis.
#'
#' @examples
#' or <- odds_ratio_function(N=c(1000,2000), Case.Rate=c(0.5), k=NULL,
#'     MAF=seq(0.05, 0.1, 0.01), power=c(0.5,0.8),Alpha=c(0.05),
#'     True.Model='All', Test.Model='All')
#' or.plot(data=or, x='MAF')
#'
#' @export
#'
or.plot<-function(data=NULL,x='MAF', panel.by='True.Model', y_limit = NULL, y_log = F, return_gg = F,
                    linear.effect.measure = 'ES',
                    select.Alpha = NULL,
                    select.power = NULL,
                    select.ES = NULL,
                    select.N = NULL,
                    select.MAF = NULL,
                    select.Case.Rate = NULL,
                    select.SD = NULL,
                    select.True.Model = NULL,
                    select.Test.Model = NULL)
{

  #library(ggplot2)
  "%ni%" <- Negate("%in%")

  #Create Dataset with separate rows for each Alpha level
  indices<-grep("OR_at_Alpha_",colnames(data))
  Alphas<-as.numeric2(unlist( strsplit(colnames(data)[indices],"OR_at_Alpha_"))[2])
    

  ss.new<-NULL
  for (i in 1:length(Alphas)){ss.new<-rbind(ss.new, data.frame(data[-indices], OR=data[,indices[i]], Alpha=Alphas[i]))}

  ss.new$Test.Model <- as.character(ss.new$Test.Model)
  #Subset data to only include those in the select. statements
  if(is.null(select.Alpha)==F){ss.new<-ss.new[ss.new$Alpha %in% select.Alpha,]}
  if(is.null(select.power)==F){ss.new<-ss.new[ss.new$power %in% select.power,]}
  if(is.null(select.ES)==F){ss.new<-ss.new[ss.new$ES %in% select.ES,]}
  if(is.null(select.N)==F){ss.new<-ss.new[ss.new$N_total %in% select.N,]}
  if(is.null(select.MAF)==F){ss.new<-ss.new[ss.new$MAF %in% select.MAF,]}
  if(is.null(select.Case.Rate)==F){ss.new<-ss.new[ss.new$Case.Rate %in% select.Case.Rate,]}
  if(is.null(select.SD)==F){ss.new<-ss.new[ss.new$SD_Y %in% select.SD,]}
  if(is.null(select.True.Model)==F){ss.new<-ss.new[ss.new$True.Model %in% select.True.Model,]}
  if(is.null(select.Test.Model)==F){ss.new<-ss.new[ss.new$Test.Model %in% select.Test.Model,]}


  var<-c("MAF", 'power','ES','R2','N_total', 'Case.Rate', 'SD_Y','Alpha', 'True.Model')
  if(linear.effect.measure %ni% c("ES", "R2")) stop("'linear.effect.measure should be either 'ES' or 'R2'")
  if(linear.effect.measure=="ES"){var <- var[!(var=='R2')]}
  if(linear.effect.measure=="R2"){var <- var[!(var=='ES')]}
  var<-var[!(var %in% c(x, panel.by))]
  var<-var[var %in% colnames(ss.new)]
  graphs<-unique.data.frame(ss.new[,var])

  if(return_gg) resplots <- list()

  for(j in 1:nrow(graphs)){
    subtitle<-paste("OR by ", x, ": ",
       var[1], '=', graphs[j,1], ', ',
       var[2], '=', graphs[j,2],', ',
       var[3], '=', graphs[j,3],', ',
       var[4], '=', graphs[j,4], sep='')
    temp1<-ss.new[ss.new[,var[1]]==graphs[j,1] & ss.new[,var[2]]==graphs[j,2]
                  &ss.new[,var[3]]==graphs[j,3] & ss.new[,var[4]]==graphs[j,4],]

    temp2<-temp1[order(temp1$True.Model),]
    temp2[,panel.by]<-paste(panel.by,"=", temp2[,panel.by])

    plot_obj <- ggplot(data=temp2, aes(x=temp2[,x], y=OR, group = Test.Model, colour = Test.Model)) +
            geom_line() +
            geom_point( size=4, shape=21, fill="white") + xlab(x)+ggtitle(subtitle)+
            facet_wrap( ~temp2[,panel.by])
    if(!is.null(y_limit)){
      plot_obj <- plot_obj + ylim(y_limit)
    }
    if(y_log){
      plot_obj <- plot_obj + scale_y_log10()
    }
    if(return_gg){ 
      resplots <- c(resplots, list(plot_obj))
      return(resplots)
    }else print(plot_obj)


  }
}
