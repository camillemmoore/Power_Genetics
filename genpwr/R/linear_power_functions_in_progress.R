#' Function to Calculate Power for Linear Models
#'
#' Calculates the power to detect an difference in means/effect size/regression coefficient, at a given sample size, N, with type 1 error rate, Alpha
#'
#' @param N Vector of the desired sample size(s)
#' @param Alpha the desired type 1 error rate(s)
#' @param MAF Vector of minor allele frequencies
#' @param sd_y Standard deviation of the outcome in the population (ignoring genotype). Either sd_y_x or sd_y must be specified.
#' @param ES Vector of effect sizes (difference in means) to detect. Either ES or R2 must be specified.
#' @param R2 Vector of R-squared values to detect. Either ES or R2 must be specified.
#' @param True.Model A vector specifying the true underlying genetic model(s): 'Dominant', 'Additive1', 'Additive2', 'Recessive' or 'All'
#' @param Test.Model A vector specifying the assumed genetic model(s) used in testing: 'Dominant', 'Additive', 'Recessive' or 'All'
#'
#' @return A data frame including the power for all combinations of the specified parameters (Case.Rate, ES, Power, etc)
#'
#' @examples
#' pw <- power.calc.linear(N=c(1000,2000),
#'     MAF=seq(0.05, 0.1, 0.01), ES=c(3,4),sd_y = c(1,2),Alpha=c(0.05),
#'     True.Model='All', Test.Model='All')
#'
#' @export
#'
power.calc.linear<-function(N=NULL, MAF=NULL, ES=NULL,R2=NULL, sd_y=NULL,
           Alpha=0.05, True.Model='All', Test.Model='All'){

    ############################################################################################################
    #Error Messages for insufficient sample size information, MAF, and case vs. control ratio
    ############################################################################################################
    if(is.null(N)==T ) {
      stop("N, the total sample size, must be specified.")}

    if(is.null(MAF)==T){
      stop("MAF (minor allele frequency) must be specified.")
    }

    if(is.null(ES)==T & is.null(R2)==T){
      stop("Either ES (detectable effect size) or R2 (detectable R-squared) must be specified.")
    }

    if(is.null(ES)==F & is.null(R2)==F){
      stop("Specify either ES (detectable effect size) or R2 (detectable R-squared), not both.")
    }

    if(is.null(sd_y)==T){
      stop("sd_y, the standard deviation of the outcome in the overall population, must be specified.")
    }

    ############################################################################################################
    #Error Messages for out of range values
    ############################################################################################################
    if(sum(sd_y<=0)>0){
      stop("sd_y must be greater than 0.")
    }

    if(sum(R2>=1)>0 | sum(R2<=0)>0){
      stop("R2 must be greater than 0 and less than 1.")
    }

    if(sum(MAF>=1)>0 | sum(MAF<=0)>0){
      stop("MAF must be greater than 0 and less than 1.")
    }

    if(sum(N<=0)>0){
      stop("N must be greater than 0.")
    }

    if(sum(Alpha>=1)>0 | sum(Alpha<=0)>0){
      stop("Alpha must be greater than 0 and less than 1.")
    }

    if(sum(!(Test.Model %in% c("Dominant", "Recessive", "Additive", "2df", "All")))>0){
      stop(paste("Invalid Test.Model:",
               paste(Test.Model[!(Test.Model %in% c("Dominant", "Recessive", "Additive", "2df", "All"))], collapse=', ')))
    }

    if(sum(!(True.Model %in% c("Dominant", "Recessive", "Additive1", "Additive2", "All")))>0){
      stop(paste("Invalid True.Model:",
               paste(True.Model[!(True.Model %in% c("Dominant", "Recessive", "Additive1", "Additive2", "All"))], collapse=', ')))
    }
    ############################################################################################################
    #Create model vectors if model = 'All'
    ############################################################################################################
    #Test model vector
    if('All' %in% Test.Model){Test.Model<-c("Dominant", "Recessive", "Additive", "2df")}

    #True model vector
    if('All' %in% True.Model){True.Model<-c("Dominant", "Recessive", "Additive1", "Additive2")}


    ############################################################################################################
    # Calculate var_x (genotype) to do effect size calculations
    ############################################################################################################
      var_x_dom = (1^2)*(1-(1-MAF)^2)-(1*(1-(1-MAF)^2))^2
      var_x_add = (1^2)*(2*MAF*(1-MAF))+(2^2)*(MAF^2)-(1*(2*MAF*(1-MAF))+2*(MAF^2))^2
      var_x_rec = (1^2)*(MAF^2)-(1*(MAF^2))^2

      var_x <- data.frame(True.Model=c(rep('Dominant', length(var_x_dom)),
                                       rep('Additive1', length(var_x_add)),
                                       rep('Additive2', length(var_x_add)),
                                       rep('Recessive', length(var_x_rec))),
                          MAF = rep(MAF, 4),
                          var_x = c(var_x_dom, var_x_add, var_x_add, var_x_rec)
                          )

  ############################################################################################################
  # Create a data.frame with all possible combinations of MAF, SD and effect size
  # Will calculate power for each of these scenarios
  ############################################################################################################
      # Depending on if ES or R2 is calculated, calcuate the other effect size measure
      if (is.null(ES)==T){e.save.tab = expand.grid(True.Model, MAF, sd_y, R2)
        colnames(e.save.tab) <- c('True.Model','MAF', 'sd_y', 'R2')
        e.save.tab <- merge(e.save.tab, var_x)
        e.save.tab$ES <- sqrt(e.save.tab$R2)*e.save.tab$sd_y/sqrt(e.save.tab$var_x)
        e.save.tab <- e.save.tab[!(e.save.tab$True.Model=='Additive1'),]
      }

      if (is.null(R2)==T){e.save.tab = expand.grid(True.Model, MAF, sd_y, ES)
        colnames(e.save.tab) <- c("True.Model", 'MAF', 'sd_y', 'ES')
        e.save.tab <- merge(e.save.tab, var_x)
        e.save.tab$R2 <- (e.save.tab$ES^2)*e.save.tab$var_x/(e.save.tab$sd_y^2)
        if(sum(e.save.tab$R2>=1)>0){
          excluded <- e.save.tab[e.save.tab$R2>=1,c("True.Model", "MAF", "sd_y", "ES", "R2")]
          e.save.tab <- e.save.tab[e.save.tab$R2<1,]

          message("Some combinations of the specified ES and sd_y imply R2>1 and 0 or negative variance of y within a genotype.\n
                  Power was not calculated for the following combinations: \n", paste(capture.output(print(excluded)), collapse = "\n"))
        }
        if(nrow(e.save.tab)==0){
          stop("All combinations of the specified ES and sd_y imply R2>1 and 0 or negative variance of y within a genotype.\n
                  Power could not be calculated. Try using smaller ES and/or larger sd_y.")
        }
      }

      # For each scenario calculate the true differences in means AB-AA and BB-AA
      e.save.tab$es_ab = ifelse(e.save.tab$True.Model=='Dominant', e.save.tab$ES,
                              ifelse(e.save.tab$True.Model=='Recessive', 0,
                                     ifelse(e.save.tab$True.Model=='Additive1', 0.5*e.save.tab$ES,
                                            e.save.tab$ES)))

      e.save.tab$es_bb = ifelse(e.save.tab$True.Model=='Dominant', e.save.tab$ES,
                              ifelse(e.save.tab$True.Model=='Recessive', e.save.tab$ES,
                                     ifelse(e.save.tab$True.Model=='Additive1', e.save.tab$ES,
                                            2*e.save.tab$ES)))

      # For each scenario calculate the Likelihoof and SD of Y given X for the Null/Intercept only model
      # sd_y_x_null <- mapply(function(x){linear.sds(e.save.tab[x,'MAF'], e.save.tab[x,'es_ab'], e.save.tab[x,'es_bb'], e.save.tab[x,'sd_y'],model = "Null")},
      #                      seq(1:nrow(e.save.tab)))
      #NOTE -THE SD_Y_X UNDER THE NULL IS JUST SD_Y...WILL NO LONGER HOLD IF CONTROLLING FOR A CONFOUNDER, FOR EX.
      ll.null = mapply(function(x){calc.like.linear(linear.mles(e.save.tab[x,'MAF'], e.save.tab[x,'es_ab'], e.save.tab[x,'es_bb'], model = 'null'),
                                                    e.save.tab[x,'MAF'], e.save.tab[x,'es_ab'], e.save.tab[x,'es_bb'], e.save.tab[x, 'sd_y'], model='null')}, seq(1:nrow(e.save.tab)))


    ############################################################################################################
    # Calculate Power for each scenario in e.save.tab under the specified testing model
    ############################################################################################################
      power.tab <- NULL

    ############################################################################################################
    #Loop over sample size
    ############################################################################################################
    for (n in N){
      ################################################################################################
      #Loop over all of the testing models and calculate power for each ES, SD, and MAF scenario
      ################################################################################################
      for (mod in Test.Model){
        # Calculate SD of Y given X for each scenario, given the test model
          sd_y_x <- mapply(function(x){linear.sds(e.save.tab[x,'MAF'], e.save.tab[x,'es_ab'], e.save.tab[x,'es_bb'], e.save.tab[x,'sd_y'],model = mod)},
                         seq(1:nrow(e.save.tab)))

          ll.alt = mapply(function(x){calc.like.linear(linear.mles(e.save.tab[x,'MAF'], e.save.tab[x,'es_ab'], e.save.tab[x,'es_bb'], model = mod),
                                    e.save.tab[x,'MAF'], e.save.tab[x,'es_ab'], e.save.tab[x,'es_bb'], sd_y_x[x], model=mod)}, seq(1:nrow(e.save.tab)))

          ll.stat = 2*(ll.alt-ll.null)

        #Calculate the power for the given sample size for a range of Alpha levels
        if(mod=='2df'){pow = mapply(function(stat) 1-pchisq(qchisq(1-Alpha, df=2, ncp=0), df=2, ncp = n*stat), ll.stat)
        }else{pow = mapply(function(stat) pnorm(sqrt(n*stat) - qnorm(1-Alpha/2))+pnorm(-sqrt(n*stat) - qnorm(1-Alpha/2)), ll.stat)
        }
        if(length(Alpha)>1){pow <- t(pow)
        rownames(pow) <- seq(1:nrow(pow))}

        #Save the power calculations for each testing model in a final table for the sample size and case rate
        power.tab<-rbind(power.tab,data.frame(Test.Model=mod, True.Model = as.character(e.save.tab[, 'True.Model']),
                 MAF = e.save.tab[, 'MAF'], N=n, ES=e.save.tab[, "ES"] , R2=e.save.tab[, "R2"], SD=e.save.tab[, "sd_y"], ES_AB = e.save.tab[, "es_ab"], ES_BB = e.save.tab[, "es_bb"], pow),row.names = NULL)

      }}
      colnames(power.tab)<-c('Test.Model', 'True.Model', 'MAF', 'N_total','ES', 'R2','SD_Y','ES_AB', 'ES_BB',
                             paste("Power_at_Alpha_", Alpha, sep=''))

    return(power.tab)
  }


