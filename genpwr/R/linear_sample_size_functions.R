#' Function to Calculate Sample Size in Linear Models
#'
#' Calculates the necessary sample size to acheive the specified level of power to detect an effect size, ES or R2 value, with type 1 error rate, Alpha
#'
#' @param power Vector of the desired power(s)
#' @param Alpha the desired type 1 error rate(s)
#' @param MAF Vector of minor allele frequencies
#' @param sd_y Standard deviation of the outcome in the population (ignoring genotype). Either sd_y_x or sd_y must be specified.
#' @param ES Vector of effect sizes (difference in means) to detect. Either ES or R2 must be specified.
#' @param R2 Vector of R-squared values to detect. Either ES or R2 must be specified.
#' @param True.Model A vector specifying the true underlying genetic model(s): 'Dominant', 'Additive', 'Recessive' or 'All'
#' @param Test.Model A vector specifying the assumed genetic model(s) used in testing: 'Dominant', 'Additive', 'Recessive' or 'All'
#'
#' @return A data frame including the total number of subjects required for all combinations of the specified parameters
#'
#' @examples
#' ss <- ss.calc.linear(power=0.8,MAF=0.1,
#'     ES=3, R2=NULL, sd_y = 1,Alpha=0.05,
#'     True.Model='All', Test.Model='All')
#'
#' @export
#'
ss.calc.linear<-
  function(power=0.8, MAF=NULL, ES=NULL, R2=NULL, sd_y=NULL,
           Alpha=0.05, True.Model='All', Test.Model='All')
{

    ############################################################################################################
    #Error Messages for insufficient sample size information, MAF, and case vs. control ratio
    ############################################################################################################
    if(is.null(power)==T ) {
      stop("Power must be specified.")}

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

    if(sum(Alpha>=1)>0 | sum(Alpha<=0)>0){
      stop("Alpha must be greater than 0 and less than 1.")
    }

    if(sum(power>=1)>0 | sum(power<=0)>0){
      stop("MAF must be greater than 0 and less than 1.")
    }

    if(sum(!(Test.Model %in% c("Dominant", "Recessive", "Additive", "2df", "All")))>0){
      stop(paste("Invalid Test.Model:",
                 paste(Test.Model[!(Test.Model %in% c("Dominant", "Recessive", "Additive", "2df", "All"))], collapse=', ')))
    }

    if(sum(!(True.Model %in% c("Dominant", "Recessive", "Additive", "All")))>0){
      stop(paste("Invalid True.Model:",
                 paste(True.Model[!(True.Model %in% c("Dominant", "Recessive", "Additive", "All"))], collapse=', ')))
    }
    ############################################################################################################
    #Create model vectors if model = 'All'
    ############################################################################################################
    #Test model vector
    if('All' %in% Test.Model){Test.Model<-c("Dominant", "Recessive", "Additive", "2df")}

    #True model vector
    if('All' %in% True.Model){True.Model<-c("Dominant", "Recessive", "Additive")}

    ############################################################################################################
    # Calculate var_x (genotype) to do effect size calculations
    ############################################################################################################
    var_x_dom = (1^2)*(1-(1-MAF)^2)-(1*(1-(1-MAF)^2))^2
    var_x_add = (1^2)*(2*MAF*(1-MAF))+(2^2)*(MAF^2)-(1*(2*MAF*(1-MAF))+2*(MAF^2))^2
    var_x_rec = (1^2)*(MAF^2)-(1*(MAF^2))^2

    var_x <- data.frame(True.Model=c(rep('Dominant', length(var_x_dom)),
                                     rep('Additive', length(var_x_add)),
                                     rep('Recessive', length(var_x_rec))),
                        MAF = rep(MAF, 3),
                        var_x = c(var_x_dom, var_x_add, var_x_rec)
    )

    ############################################################################################################
    # Create a data.frame with all possible combinations of MAF, SD and effect size
    # Will calculate power for each of these scenarios
    ############################################################################################################
    # Depending on if ES or R2 is calculated, calcuate the other effect size measure
    if (is.null(ES)){
      e.save.tab = expand.grid(True.Model, MAF, sd_y, R2)
      colnames(e.save.tab) <- c('True.Model','MAF', 'sd_y', 'R2')
      e.save.tab <- merge(e.save.tab, var_x)
      e.save.tab$ES <- sqrt(e.save.tab$R2)*e.save.tab$sd_y/sqrt(e.save.tab$var_x)
      # e.save.tab <- e.save.tab[!(e.save.tab$True.Model=='Additive1'),]
    }

    if (is.null(R2)){
      e.save.tab = expand.grid(True.Model, MAF, sd_y, ES)
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
    e.save.tab$es_ab = ifelse(e.save.tab$True.Model=='Recessive', 0, e.save.tab$ES)

    e.save.tab$es_bb = ifelse(e.save.tab$True.Model=='Additive', 2*e.save.tab$ES,e.save.tab$ES)

    e.save.tab$True.Model <- as.character(e.save.tab$True.Model)

    # For each scenario calculate the SD of Y give X for the true model
    e.save.tab$sd_y_x_true = mapply(function(x){linear.sds(e.save.tab[x,'MAF'], e.save.tab[x,'es_ab'], e.save.tab[x,'es_bb'], e.save.tab[x,'sd_y'],model = e.save.tab[x,'True.Model'])},
                                    seq(1:nrow(e.save.tab)))


    # For each scenario calculate the Likelihood and SD of Y given X for the Null/Intercept only model
    # sd_y_x_null <- mapply(function(x){linear.sds(e.save.tab[x,'MAF'], e.save.tab[x,'es_ab'], e.save.tab[x,'es_bb'], e.save.tab[x,'sd_y'],model = "Null")},
    #                      seq(1:nrow(e.save.tab)))
    #NOTE -THE SD_Y_X UNDER THE NULL IS JUST SD_Y...WILL NO LONGER HOLD IF CONTROLLING FOR A CONFOUNDER, FOR EX.
    ll.null = mapply(function(x){calc.like.linear(linear.mles(e.save.tab[x,'MAF'], e.save.tab[x,'es_ab'], e.save.tab[x,'es_bb'], model = 'null'),
                                                  e.save.tab[x,'MAF'],
                                                  e.save.tab[x,'es_ab'],
                                                  e.save.tab[x,'es_bb'],
                                                  e.save.tab[x, 'sd_y'],
                                                  e.save.tab[x, 'sd_y_x_true'],
                                                  model='null')}, seq(1:nrow(e.save.tab)))

    ############################################################################################################
    # Calculate Sample Size for each scenario in e.save.tab under the specified testing model
    ############################################################################################################
    ss.tab <- NULL

    ############################################################################################################
    #Loop over power
    ############################################################################################################
    for (p in power){
      ################################################################################################
      #Loop over all of the testing models and calculate power for each ES, SD, and MAF scenario
      ################################################################################################
      for (mod in Test.Model){
        # Calculate SD of Y given X for each scenario, given the test model
        sd_y_x <- mapply(function(x){linear.sds(e.save.tab[x,'MAF'], e.save.tab[x,'es_ab'], e.save.tab[x,'es_bb'], e.save.tab[x,'sd_y'],model = mod)},
                         seq(1:nrow(e.save.tab)))

        ll.alt = mapply(function(x){calc.like.linear(linear.mles(e.save.tab[x,'MAF'], e.save.tab[x,'es_ab'], e.save.tab[x,'es_bb'], model = mod),
                                                     e.save.tab[x,'MAF'],
                                                     e.save.tab[x,'es_ab'],
                                                     e.save.tab[x,'es_bb'],
                                                     sd_y_x[x],
                                                     sd_y_x_truth = e.save.tab[x, "sd_y_x_true"],
                                                     model=mod)}, seq(1:nrow(e.save.tab)))

        ll.stat = 2*(ll.alt-ll.null)
        #Calculate the SS for the given power for a range of Alpha levels
        if(mod=='2df'){ss<-NULL
        for (q in 1:length(Alpha)){
          ss = rbind(ss, mapply(function(stat){uniroot(function(x) ncp.search(x, p, Alpha[q], df=2),
                             lower=0, upper=1000, extendInt = 'upX', tol=0.00001)$root/stat}, ll.stat))
          }
          ss<-t(ss)
        }else{ss<-NULL
        for (q in 1:length(Alpha)){
          ss = rbind(ss, mapply(function(stat){uniroot(function(x) ncp.search(x, p, Alpha[q], df=1),
                                                           lower=0, upper=1000, extendInt = 'upX', tol=0.00001)$root/stat}, ll.stat))
          #ss = mapply(function(stat)(qnorm(1-Alpha/2)+qnorm(p))^2/stat, ll.stat)
        }
          ss<-t(ss)  
        #if(length(Alpha)>1){ss<-t(ss)}
        }



        #Save the power calculations for each testing model in a final table for the sample size and case rate
        ss.tab<-rbind(ss.tab,data.frame(Test.Model=mod, True.Model = as.character(e.save.tab[, 'True.Model']),
                                              MAF = e.save.tab[, 'MAF'], power=p, ES=e.save.tab[, "ES"] , R2=e.save.tab[, "R2"], SD=e.save.tab[, "sd_y"],
                                       ES_AB = e.save.tab[, "es_ab"], ES_BB = e.save.tab[, "es_bb"], ss),row.names = NULL)

      }}
    colnames(ss.tab)<-c('Test.Model', 'True.Model', 'MAF', 'Power','ES', 'R2','SD_Y','ES_AB', 'ES_BB',
                           paste("N_total_at_Alpha_", Alpha, sep=''))

    return(ss.tab)
}

