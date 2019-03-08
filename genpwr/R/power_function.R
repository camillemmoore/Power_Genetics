#' Function to Calculate Power
#'
#' Calculates the power to detect an odds ratio, OR, at a given sample size, N, with type 1 error rate, Alpha
#'
#' @param N Vector of the desired sample size(s)
#' @param Case.Rate Vector of the proportion(s) of cases in the sample (cases/(cases + controls)). Either k or Case.Rate must be specified.
#' @param k Vector of the number of controls per case. Either k or Case.Rate must be specified.
#' @param Alpha the desired type 1 error rate(s)
#' @param MAF Vector of minor allele frequencies
#' @param OR Vector of odds ratios to detect
#' @param True.Model A vector specifying the true underlying genetic model(s): 'Dominant', 'Additive', 'Recessive' or 'All'
#' @param Test.Model A vector specifying the assumed genetic model(s) used in testing: 'Dominant', 'Additive', 'Recessive' or 'All'
#'
#' @return A data frame including the power for all combinations of the specified parameters (Case.Rate, OR, Power, etc)
#'
#' @examples
#' pw <- power.calc(N=2000, Case.Rate=0.5, k=NULL,
#'     MAF=0.2, OR=1.5,Alpha=0.05,
#'     True.Model='All', Test.Model='All')
#'
#' @export
#'
power.calc<-
  function(N=NULL, Case.Rate=NULL, k=NULL, MAF=NULL, OR=NULL,
                     Alpha=0.05, True.Model='All', Test.Model='All')
{

  ############################################################################################################
  #Error Messages for insufficient sample size information, MAF, and case vs. control ratio
  ############################################################################################################
  if(is.null(N)==T ) {
    stop("N, the total sample size, must be specified.")}

  if(is.null(k)==T & is.null(Case.Rate)==T){
    stop("k, the number of controls per case, or Case.Rate, the proportion of cases in the study sample, must be specified.")
  }

  if(is.null(k)==F & is.null(Case.Rate)==F){
    stop("Specify one of k, the number of controls per case, or Case.Rate, the proportion of cases in the study sample, not both.")
  }

  if(is.null(MAF)==T){
    stop("MAF (minor allele frequency) must be specified.")
  }

  if(is.null(OR)==T){
    stop("OR (detectable odds ratio) must be specified.")
  }

  ############################################################################################################
  #Error Messages for out of range values
  ############################################################################################################

  if(sum(Case.Rate>=1)>0 | sum(Case.Rate<=0)>0){
    stop("R2 must be greater than 0 and less than 1.")
  }

  if(sum(MAF>=1)>0 | sum(MAF<=0)>0){
    stop("MAF must be greater than 0 and less than 1.")
  }

  if(sum(N<=0)>0){
    stop("N must be greater than 0.")
  }

  if(sum(k<=0)>0){
    stop("k must be greater than 0.")
  }

  if(sum(OR<=0)>0){
    stop("OR must be greater than 0.")
  }

  if(sum(Alpha>=1)>0 | sum(Alpha<=0)>0){
    stop("Alpha must be greater than 0 and less than 1.")
  }

  if(sum(!(Test.Model %in% c("Dominant", "Recessive", "Additive", "2df", "All")))>0){
    stop(paste("Invalid Test.Model:",
               paste(Test.Model[!(Test.Model %in% c("Dominant", "Recessive", "Additive", "2df", "All"))], collapse=', ')))
  }

  if(sum(!(True.Model %in% c("Dominant", "Recessive", "Additive", "Additive", "All")))>0){
    stop(paste("Invalid True.Model:",
               paste(True.Model[!(True.Model %in% c("Dominant", "Recessive", "Additive", "All"))], collapse=', ')))
  }
  ############################################################################################################
  #Calculate needed sample size information from provided inputs
  ############################################################################################################
  #Test model vector
  if('All' %in% Test.Model){Test.Model<-c("Dominant", "Recessive", "Additive", "2df")}

  #True model vector
  # True.Model[True.Model == "Additive"] <- "Additive"
  if('All' %in% True.Model){True.Model<-c("Dominant", "Recessive", "Additive")}

  #If k is provided calculate the Case.Rate
  if(is.null(Case.Rate)==T){Case.Rate = 1/(1+k)}

  #Find all possible combinations of N and case_rate
  sample.size.tab <- expand.grid(N, Case.Rate)
  colnames(sample.size.tab) <- c('N', 'Case.Rate')
  sample.size.tab$N_cases <- sample.size.tab$N*sample.size.tab$Case.Rate
  sample.size.tab$N_controls <- sample.size.tab$N - sample.size.tab$N_cases

  iter <- nrow(sample.size.tab)

  final.pow.tab <- NULL

  ############################################################################################################
  #Loop over sample size and Case.Rate combiantions, calculate power for all OR's and MAF's
  ############################################################################################################

  for (zz in 1:iter){
    N <- sample.size.tab[zz,'N']
    Case.Rate <- sample.size.tab[zz,'Case.Rate']
    N_cases <- sample.size.tab[zz,'N_cases']
    N_controls <- sample.size.tab[zz,'N_controls']

    ############################################################################################################
    #Use OR's and MAF's to calculate true distriubiton of genotypes and disease
    ############################################################################################################
    ##############################################################################
    #For each OR
    ##############################################################################
    o.save.tab <-NULL

    for (o in OR){

      m.save.tab<-NULL

      ##########################################################################
      #For each MAF
      ##########################################################################
      for (m in MAF){
        #Temporary place to save 2x3 tables for each true model with MAF=m and OR=o
        save.tab <- NULL

        #Proportion with each genotype. This is the same for all true.models.
        P_AA <- (1-m)^2
        P_AB <- 2*m*(1-m)
        P_BB <- m^2

        ############################################################################
        #Create 2x3 tables of joint probabilities for each true model of interest
        ############################################################################

        if('Dominant' %in% True.Model){
          #Dominant Model
          #Solve Quadratic Equation to get 2x3 tables corresponding to the OR
          a <- (1-o)
          b <- (P_AA-Case.Rate+o*(Case.Rate+P_AB+P_BB))
          c <- -o*(P_AB+P_BB)*Case.Rate

          soln <- quad_roots(a,b,c)[2]

          #Under a dominant model
          #Proabilities of disease conditional on genotype
          P_AA_case_d <- (Case.Rate-soln)/P_AA
          P_AB_case_d <- P_BB_case_d <- soln/(P_AB+P_BB)

          #Joint probabilities of disease and genotype
          prob_AA_case_d <- P_AA_case_d*P_AA
          prob_AB_case_d <- P_AB_case_d*P_AB
          prob_BB_case_d <- P_BB_case_d*P_BB
          prob_AA_control_d <- (1-P_AA_case_d)*P_AA
          prob_AB_control_d <- (1-P_AB_case_d)*P_AB
          prob_BB_control_d <- (1-P_BB_case_d)*P_BB

          dom.tab <- data.frame(model=rep('Dominant',2),table=rbind(c(prob_AA_case_d, prob_AB_case_d, prob_BB_case_d),
                                                                    c(prob_AA_control_d, prob_AB_control_d,prob_BB_control_d)))


          save.tab<-rbind(save.tab, dom.tab)
        }

        # if('Additive1' %in% True.Model){
        #   a <- (sqrt(o)-1)
        #   b <- (P_AB+sqrt(o)*P_BB+Case.Rate-Case.Rate*sqrt(o))
        #   c <- -P_AB*Case.Rate
        #   soln <- quad_roots(a,b,c)[2]
        #   upper.lim<-min(soln, P_AB)

        #   fa.1<-function(x){sqrt(o)-x*(P_AA-Case.Rate+x+((sqrt(o)*x*P_BB)/(P_AB-x+sqrt(o)*x)))/((Case.Rate-x-((sqrt(o)*x*P_BB)/(P_AB-x+sqrt(o)*x)))*(P_AB-x))}

        #   add1.root<-zero_finder_nleqslv(fa.1,veclength = 1, x.start.vals = runif(1)*soln)


        #   #Proabilities of disease conditional on genotype
        #   P_AB_case_a1 <- add1.root/P_AB

        #   #Joint probabilities of disease and genotype
        #   prob_AB_case_a1 <- P_AB_case_a1*P_AB
        #   prob_AB_control_a1 <- (1-P_AB_case_a1)*P_AB
        #   prob_AA_case_a1 <-P_AA*prob_AB_case_a1/(sqrt(o)*prob_AB_control_a1+prob_AB_case_a1)
        #   prob_AA_control_a1 <- P_AA-prob_AA_case_a1
        #   prob_BB_case_a1 <- (prob_AA_case_a1*P_BB*o)/(prob_AA_case_a1*o + prob_AA_control_a1)
        #   prob_BB_control_a1 <- P_BB-prob_BB_case_a1
        #   if(length(add1.root) > 1){
        #     qqcr <- numeric(0)
        #     for(aqq in 1:length(add1.root)){
        #       qqt <- rbind(c(prob_AA_case_a1[aqq], prob_AB_case_a1[aqq], prob_BB_case_a1[aqq]),
        #         c(prob_AA_control_a1[aqq], P_AB-prob_AB_case_a1[aqq],prob_BB_control_a1[aqq]))
        #       qqcr <- c(qqcr, apply(qqt, 1, sum)[1])
        #     }
        #     qqz <- which(qqcr - Case.Rate < 1e-3)
        #     if(length(qqz) > 1) stop("trouble selecting correct zero in additive function")
        #     prob_AA_case_a2 <- prob_AA_case_a2[qqz]
        #     prob_AB_case_a2 <- prob_AB_case_a2[qqz]
        #     prob_BB_case_a2 <- prob_BB_case_a2[qqz]
        #     prob_AA_control_a2 <- prob_AA_control_a2[qqz]
        #     prob_AB_case_a2 <- prob_AB_case_a2[qqz]
        #     prob_BB_control_a2 <- prob_BB_control_a2[qqz]
        #   }

        #   add.tab1<-data.frame(model=rep('Additive1',2),table=rbind(c(prob_AA_case_a1, prob_AB_case_a1, prob_BB_case_a1),
        #                                                             c(prob_AA_control_a1, P_AB-prob_AB_case_a1,prob_BB_control_a1)))

        #   save.tab<-rbind(save.tab, add.tab1)
        # }

        if('Additive' %in% True.Model){
          a <- (o-1)
          b <- (P_AB+o*P_BB+Case.Rate-Case.Rate*o)
          c <- -P_AB*Case.Rate
          soln <- quad_roots(a,b,c)[2]
          upper.lim<-min(soln, P_AB)

          fa.2<-function(x){o-x*(P_AA-Case.Rate+x+((o*x*P_BB)/(P_AB-x+o*x)))/((Case.Rate-x-((o*x*P_BB)/(P_AB-x+o*x)))*(P_AB-x))}

          add2.root<-zero_finder_nleqslv(fa.2,veclength = 1, x.start.vals = runif(1)*soln)

          #Proabilities of disease conditional on genotype
          P_AB_case_a2 <- add2.root/P_AB

          #Joint probabilities of disease and genotype
          prob_AB_case_a2 <- P_AB_case_a2*P_AB
          prob_AB_control_a2 <- (1-P_AB_case_a2)*P_AB
          prob_AA_case_a2 <-P_AA*prob_AB_case_a2/(o*prob_AB_control_a2+prob_AB_case_a2)
          prob_AA_control_a2 <- P_AA-prob_AA_case_a2
          prob_BB_case_a2 <- (prob_AA_case_a2*P_BB*o^2)/(prob_AA_case_a2*o^2 + prob_AA_control_a2)
          prob_BB_control_a2 <- P_BB-prob_BB_case_a2
          if(length(add2.root) > 1){
            qqcr <- numeric(0)
            for(aqq in 1:length(add2.root)){
              qqt <- rbind(c(prob_AA_case_a2[aqq], prob_AB_case_a2[aqq], prob_BB_case_a2[aqq]),
                c(prob_AA_control_a2[aqq], P_AB-prob_AB_case_a2[aqq],prob_BB_control_a2[aqq]))
              qqcr <- c(qqcr, apply(qqt, 1, sum)[1])
            }
            qqz <- which(qqcr - Case.Rate < 1e-3)
            if(length(qqz) > 1) qqz <- which.min(qqcr - Case.Rate < 1e-3)
            if(length(qqz) > 1) stop("trouble selecting correct zero in additive function")
            prob_AA_case_a2 <- prob_AA_case_a2[qqz]
            prob_AB_case_a2 <- prob_AB_case_a2[qqz]
            prob_BB_case_a2 <- prob_BB_case_a2[qqz]
            prob_AA_control_a2 <- prob_AA_control_a2[qqz]
            prob_AB_case_a2 <- prob_AB_case_a2[qqz]
            prob_BB_control_a2 <- prob_BB_control_a2[qqz]
          }

          add.tab2<-data.frame(model=rep('Additive',2),table=rbind(c(prob_AA_case_a2, prob_AB_case_a2, prob_BB_case_a2),
                                                                    c(prob_AA_control_a2, P_AB-prob_AB_case_a2,prob_BB_control_a2)))

          save.tab<-rbind(save.tab, add.tab2)
        }

        if('Recessive' %in% True.Model){
          a <- (1-o)
          b <- o*P_BB+(o-1)*Case.Rate+P_AA + P_AB
          c <- -o*(P_BB)*Case.Rate

          soln <- quad_roots(a,b,c)[2]

          #Proabilities of disease conditional on genotype
          P_AB_case_r <- P_AA_case_r <- (Case.Rate-soln)/(P_AB+P_AA)
          P_BB_case_r <- soln/P_BB

          #Joint probabilities of disease and genotype
          prob_AA_case_r <- P_AA_case_r*P_AA
          prob_AB_case_r <- P_AB_case_r*P_AB
          prob_BB_case_r <- P_BB_case_r*P_BB
          prob_AA_control_r <- (1-P_AA_case_r)*P_AA
          prob_AB_control_r <- (1-P_AB_case_r)*P_AB
          prob_BB_control_r <- (1-P_BB_case_r)*P_BB

          rec.tab<-data.frame(model=rep('Recessive',2),table=rbind(c(prob_AA_case_r, prob_AB_case_r, prob_BB_case_r),
                                                                   c(prob_AA_control_r, P_AB-prob_AB_case_r,prob_BB_control_r)))

          save.tab<-rbind(save.tab, rec.tab)
        }

        m.save.tab<-rbind(m.save.tab,
                          data.frame(True.Model = save.tab[,1], MAF=m, OR = o,
                                     Disease.Status = rep(c('case', "control"),nrow(save.tab)/2),
                                     Geno.AA = save.tab[,2],Geno.AB = save.tab[,3], Geno.BB = save.tab[,4]))
      }
      o.save.tab<-rbind(o.save.tab, m.save.tab)
    }

    ############################################################################################################
    #Calculate Power for each scenario under the specified testing model
    ############################################################################################################

    power.tab<-NULL

    ################################################################################################
    #Loop over all of the testing models and calculate power for each OR and MAF scenario
    ################################################################################################
    for (mod in Test.Model){
      temp<-NULL

      #Repeat calculation for each OR/MAF combination
      for (j in seq(1, nrow(o.save.tab),2)){
        #Grab the correct 2x3 table of probabilities
        t<-o.save.tab[j:(j+1),c("Geno.AA", "Geno.AB", "Geno.BB")]

        #Calculate the null and alternative likelihoods
        ll.alt<-calc.like(logistic.mles(t, model = mod), t, model=mod)
        ll.null<-null.ll(t)

        #Calculate the LRT statistic
        stat<-2*(as.numeric(ll.alt-ll.null))

        #Calculate the power for the given sample size for a range of Alpha levels
        if(mod=='2df'){pow = 1-pchisq(qchisq(1-Alpha, df=2, ncp=0), df=2, ncp = N*stat)
        }else{pow = pnorm(sqrt(N*stat) - qnorm(1-Alpha/2))+pnorm(-sqrt(N*stat) - qnorm(1-Alpha/2))*1}#

        temp<-rbind(temp, pow)
      }

      #Save the power calculations for each testing model in a final table for the sample size and case rate
      power.tab<-rbind(power.tab,data.frame(Test.Model=mod, o.save.tab[seq(1, nrow(o.save.tab),2),1:3],
                                            N, N_cases, N_controls, Case.Rate,temp))
    }
    colnames(power.tab)<-c('Test.Model', 'True.Model', 'MAF', 'OR', 'N_total', 'N_cases', 'N_controls','Case.Rate',
                           paste("Power_at_Alpha_", Alpha, sep=''))



    final.pow.tab<-rbind(final.pow.tab, power.tab)
  }
  return(final.pow.tab)
}

