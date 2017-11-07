###############################################################################################################################
###############################################################################################################################
#Function to calculate power for a given sample size
###############################################################################################################################
###############################################################################################################################

#Inputs:
#First must specify sample size and proportion of cases
#N = Sample size (total), vector
#and
#Case.Rate = overall proportion of sample with disease; in a case control study, N cases/(N cases + N controls), vector
#or
#k = controls per case, vector

#Minor Allele Freqency
#MAF, Vector

#Alpha level
#Alpha, a vector of significance levels

#Effect Size
#OR, a vector of ORs to detect

#True Model
#True.Model, Vector of any or all of the following:
#'Dominant', 'Additive1', 'Additive2', 'Recessive' or 'All'

#Testing Model
#Test.Model, Vector of any or all of the following
#'Dominant', 'Additive', '2df', 'Recessive' or 'All'

###########################################################################################################
#Power calculation function
###########################################################################################################
power.calc<-function(N=NULL, Case.Rate=NULL, k=NULL, MAF=NULL, OR=NULL,
                     Alpha=0.05, True.Model='All', Test.Model='All'){
  
  ############################################################################################################  
  #Error Messages for insufficient sample size information, MAF, and
  ############################################################################################################
  if(is.null(N)==T ) {
    stop("N, the total sample size, must be specified.")}
  
  if(is.null(k)==T & is.null(Case.Rate)==T){
    stop("k, the number of controls per case, or Case.Rate, the proportion of cases in the study sample, must be specified.")
  }
  
  if(is.null(MAF)==T){
    stop("MAF (minor allele frequency) must be specified.")
  }
  
  if(is.null(OR)==T){
    stop("OR (detectable odds ratio) must be specified.")
  }
  
  ############################################################################################################  
  #Calculate needed sample size information from provided inputs
  ############################################################################################################  
  #Test model vector
  if(Test.Model=='All'){Test.Model<-c("Dominant", "Recessive", "Additive", "2df")}
  
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
        
        #Proportion with each genotype
        P_AA <- (1-m)^2
        P_AB <- 2*m*(1-m)
        P_BB <- m^2
        
        #Solve Quadratic Equation to get 2x3 tables corresponding to the OR
        a <- (1-o)
        b <- (P_AA-Case.Rate+o*(Case.Rate+P_AB+P_BB))
        c <- -o*(P_AB+P_BB)*Case.Rate
        
        soln <- quad_roots(a,b,c)[2]
        
        #Under a dominant model
        #Proabilities of disease conditional on genotype
        P_AA_case <- (Case.Rate-soln)/P_AA
        P_AB_case <- P_BB_case <- soln/(P_AB+P_BB)
        
        #Joint probabilities of disease and genotype
        prob_AA_case <- P_AA_case*P_AA
        prob_AB_case <- P_AB_case*P_AB
        prob_BB_case <- P_BB_case*P_BB
        prob_AA_control <- (1-P_AA_case)*P_AA 
        prob_AB_control <- (1-P_AB_case)*P_AB 
        prob_BB_control <- (1-P_BB_case)*P_BB 
        
        ############################################################################
        #Create 2x3 tables of joint probabilities for each true model of interest
        ############################################################################
        
        save.tab<-NULL
        
        if('Dominant' %in% True.Model | True.Model=='All'){
          dom.tab <- data.frame(model=rep('Dominant',2),table=rbind(c(prob_AA_case, prob_AB_case, prob_BB_case),
                                                                    c(prob_AA_control, prob_AB_control,prob_BB_control)))
          
          save.tab<-rbind(save.tab, dom.tab)
        }
        
        if('Additive1' %in% True.Model | True.Model=='All'){
          OR1=exp(0.5*log(o))
          
          prob_AB_case_a1<-(OR1*P_AA_case*P_AB/(OR1*P_AA_case+P_AA-P_AA_case))
          prob_AB_control_a1<-P_AB-prob_AB_case_a1
          
          add.tab1<-data.frame(model=rep('Additive1',2),table=rbind(c(prob_AA_case, prob_AB_case_a1, prob_BB_case),
                                                                    c(prob_AA_control, prob_AB_control_a1,prob_BB_control)))
          
          save.tab<-rbind(save.tab, add.tab1)   
        }
        
        if('Additive2' %in% True.Model | True.Model=='All'){
          OR2=exp(2*log(o))
          
          prob_BB_case_a2<-(OR2*P_AA_case*P_BB/(OR2*P_AA_case+P_AA-P_AA_case))
          prob_BB_control_a2<-P_BB - prob_BB_case_a2
          
          add.tab2<-data.frame(model=rep('Additive2',2),table=rbind(c(prob_AA_case, prob_AB_case, prob_BB_case_a2),
                                                                    c(prob_AA_control, prob_AB_control,prob_BB_control_a2)))
          save.tab<-rbind(save.tab, add.tab2) 
        }
        
        if('Recessive' %in% True.Model | True.Model=='All'){
          P_AB_case_r<-P_AA_case
          prob_AB_case_r<-(P_AB_case_r*P_AB)
          
          rec.tab<-data.frame(model=rep('Recessive',2),table=rbind(c(prob_AA_case, prob_AB_case_r, prob_BB_case),
                                                                   c(prob_AA_control, P_AB-prob_AB_case_r,prob_BB_control)))
          
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
    for (mod in Test.Model){temp<-NULL
    
    #Repeat calcualtion for each OR/MAF combination
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
      }else{pow = pnorm(sqrt(N*stat) - qnorm(1-Alpha/2))+pnorm(-sqrt(N*stat) - qnorm(1-Alpha/2))}
      
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


