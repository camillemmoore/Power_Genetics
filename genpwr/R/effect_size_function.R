#Should make this into a function called effect_size_function

#GOAL: taking input of vectors of power, alpha, sample size, MAF and case rate, true model and test model,
# return a table of detectable odds ratios
# user should specify if m is a risk allele or protective allele (OR's >1 or <1)

####################################################################################
# First identify the detectable LRT test statistic
#####################################################################################
#For 1DF Test.Models the detectable LRT test statistic is:

####################################################################################
# odds ratios for different models
####################################################################################

or_calc <- function(a,b,c,d,e,f, mod, risk_allele)
{
  if(mod == "Recessive"){
    if(risk_allele) return(c(c*e/(b*f))[c(c*e/(b*f)) > 1])
    if(!risk_allele) return(c(c*e/(b*f))[c(c*e/(b*f)) < 1])
  }
  if(mod == "Dominant"){
    if(risk_allele) return(c(b*d/(a*e))[c(b*d/(a*e)) > 1])
    if(!risk_allele) return(c(b*d/(a*e))[c(b*d/(a*e)) < 1])
  }
  if(mod == "Additive"){
    if(risk_allele) return(c(c(b*d/(a*e))[b*d/(a*e) > 1]))
    if(!risk_allele) return(c(c(b*d/(a*e))[b*d/(a*e) < 1]))
  }
}



####################################################################################
# For a dominant model
####################################################################################
dom.or.function <-function(like, Case.Rate, P_AA, P_AB, P_BB, True.Model, risk_allele)
{
  if(F){ #to delete
    ############################################################################################################
    #log likelihood stats
    ############################################################################################################

    #For 1DF Test.Models the detectable LRT test statistic is:
    stat = ((qnorm(1-alpha/2)+qnorm(power))^2)/N


    ####################################################################################
    # Calculate the null and alternative likelihoods
    ####################################################################################

    # Can calculate the null log likelihood from the case.rate alone (intercept only model)
    beta0 <- logit(Case.Rate)
    ll.null <- Case.Rate*log(exp(beta0)/(1+exp(beta0))) + (1-Case.Rate)*log(1/(1+exp(beta0)))

    # Note that stat = 2*(as.numeric(ll.alt-ll.null))
    # Then the alt.likelihood to detect is
    ll.alt <- ll.null+(stat/2)
    like <- ll.alt

    ####################################################################################
    # Calculate the probabilities of each genotype
    ####################################################################################
    P_AA <- (1-m)^2
    P_AB <- 2*(1-m)*m
    P_BB <- m^2

    ####################################################################################
    # Assign alternate variable names
    ####################################################################################

    cr <- Case.Rate
    m <- MAF
  }
  cr <- Case.Rate
  # m <- MAF
  # Solve for a, the joint prob of geno AA and disease;
  # a and d will be the same for all true genetic models
  # warning = error = 0 #For error handling

  # tryCatch(
    a <- ll_zero_finder2(function(x) {x*log(x/P_AA)+
        (cr-x)*log((cr-x)/(P_AB+P_BB))+
        (P_AA-x)*log((P_AA-x)/P_AA)+
        (1-cr-P_AA+x)*log((1-cr-P_AA+x)/(P_AB+P_BB)) - like})#,
      # lower=ifelse(risk_allele==T, max(cr+P_AA-1+delta, delta), P_AA*cr),
      # upper=ifelse(risk_allele==T, cr*P_AA-delta, min(P_AA, cr)-delta))$root, #this assumes that m is a risk allele/that OR's will be > 0; Consider adding in option for protective alleles
    # warning=function(warn){warning<<-1},
    # error=function(err){error<<-1})


  # if(warning+error==0){
    d <- P_AA - a

    if (True.Model == "Dominant"){
      y <- (cr-a)/(P_AB+P_BB)
      b <- y*P_AB
      c <- y*P_BB
      e <- P_AB - b
      f <- P_BB - c
    }

    if (True.Model == "Recessive"){
      # The conditional probability of disease in AA and AB should be the same
      P_AB_case = a/P_AA

      # Now solve for the elements of the recessive table
      b = P_AB_case*P_AB
      c = cr - a - b
      e = P_AB - b
      f = P_BB - c
    }

    if (True.Model == "Additive"){
      # warning<-error<-0

      # Calculate the OR under the dominant model
      # o = (cr-a)*d/(a*(1-cr-d))

      #Here to solve for b, solve for c in terms of b in two ways, first using the dominant OR relationship and second setting OR(AB vs. AA) equal to OR(BB vs AB)
      #subtract the two version of c from one another and find the root of the equation
      b <- sapply(a, function(aa){
        # tryCatch(
          bb <- ll_zero_finder2(function(x){
                        dd <- P_AA - aa
                        o = (cr-aa)*dd/(aa*(1-cr-dd))
                        (aa*(x-x*o+o*(P_BB+ P_AB)) - x*P_AA)/(aa*(-1 + o) + P_AA) - (P_BB*(P_AA-aa)*x^2)/(aa*P_AB*(P_AB-2*x)+P_AA*x^2)
                      })#,
                            # lower = ifelse(o>1, aa*P_AB/P_AA, 0), upper = ifelse(o<1, aa*P_AB/P_AA, P_AB))$root,
               # warning=function(warn){warning<<-1},
               # error=function(err){error<<-1})
        return(bb)
      })

      # if(warning+error==0){
        c = cr - a - b
        e = P_AB - b
        f = P_BB - c
      # }else{a=b=c=d=e=f=NA}
    }

    #For all models, calc the OR's
    #tab <- rbind(c(a,b,c), c(d,e,f)) #table
    if(all(!is.na(a) & a>0 & b>0 & c>0 & d>0 & e>0 & f>0)){
      return(or_calc(a,b,c,d,e,f, mod = True.Model, risk_allele))
    }else{stop("Could not calculate odds ratio for the given parameters")} #before returning anything, check that no negative values in table
  # }else{stop("Could not calculate odds ratio for the given parameters")}
}

####################################################################################
# For a recessive model
####################################################################################
rec.or.function <-function(like, Case.Rate, P_AA, P_AB, P_BB, True.Model, risk_allele)
{
  if(F){#to delete
    ############################################################################################################
    #log likelihood stats
    ############################################################################################################

     #For 1DF Test.Models the detectable LRT test statistic is:
    stat = ((qnorm(1-alpha/2)+qnorm(power))^2)/N

    #For the 2DF Test the detectable LRT test statistic is:
    stat_2df = uniroot(function(x) ncp.search(x, power, stat, alpha, df=2),
                   lower=0, upper=1000, extendInt = 'upX', tol=0.00001)$root/N


    ####################################################################################
    # Calculate the null and alternative likelihoods
    ####################################################################################

    # Can calculate the null log likelihood from the case.rate alone (intercept only model)
    beta0 <- logit(Case.Rate)
    ll.null <- Case.Rate*log(exp(beta0)/(1+exp(beta0))) + (1-Case.Rate)*log(1/(1+exp(beta0)))

    # Note that stat = 2*(as.numeric(ll.alt-ll.null))
    # Then the alt.likelihood to detect is
    ll.alt <- ll.null+(stat/2)
    like <- ll.alt

    ####################################################################################
    # Calculate the probabilities of each genotype
    ####################################################################################
    P_AA <- (1-m)^2
    P_AB <- 2*(1-m)*m
    P_BB <- m^2

    ####################################################################################
    # Assign alternate variable names
    ####################################################################################

    cr <- Case.Rate
    m <- MAF
  }
  cr <- Case.Rate
  # Solve for c, the prob of BB and disease
  # warning<-error<-0
  # tryCatch(
    c <- ll_zero_finder2(function(x) {
          x*log(x/P_BB)+
          (cr-x)*log((cr-x)/(P_AB+P_AA))+
          (P_BB-x)*log((P_BB-x)/P_BB)+
          (1-cr-P_BB+x)*log((1-cr-P_BB+x)/(P_AB+P_AA)) - like})#,
        # lower=ifelse(risk_allele==T, cr*P_BB+delta, delta),
        # upper=ifelse(risk_allele==T, P_BB-delta, cr*P_BB - delta)),
    # warning=function(warn){warning<<-1},
    # error=function(err){error<<-1}
    # )

  # if(warning+error==0){
    f <- P_BB - c

    if(True.Model=='Recessive'){
      y = (cr-c)/(P_AB+P_AA)
      a = y*P_AA
      d = P_AA - a
      b = y*P_AB
      e = P_AB - b
    }

    if(True.Model=='Dominant'){
      # the conditional prob of disease for genotypes AB and BB are equal
      P_AB_case = c/P_BB
      P_AB_control = f/P_BB #

      # Solve for the joint probabilities
      b = P_AB_case*P_AB
      d = P_AB - b
      a = cr - b - c
      d = P_AA - a
      e = P_AB_control*P_AB
    }

    if(True.Model=='Additive'){
      # Calculate the OR for the recessive model
      # o = c*(1-cr-f)/((cr-c)*f)

      # To find b, solve for a in two ways (see dominant.or.function)
      # warning<-error<-0
        # tryCatch(
          # b <- ll_zero_finder2(function(x) ((-(x^2)*c*P_AA + (x^2)*P_AA*P_BB)/(-2*x*c*P_AB + c*(P_AB^2) + (x^2)*P_BB)) -
          #                  ((-c*(P_AA + P_AB) + x*(c - c*o + o*P_BB))/(c*(-1 + o) - o*P_BB)))
                         # ,lower = ifelse(o>1,0, c*P_AB/P_BB), upper = ifelse(o<1, P_AB,c*P_AB/P_BB))#$root,
             # warning=function(warn){warning<<-1},
             # error=function(err){error<<-1})
      b <- sapply(c, function(cc){
        # tryCatch(
          bb <- ll_zero_finder2(function(x){
                        ff <- P_BB - cc
                        o = cc*(1-cr-ff)/((cr-cc)*ff)
                        (-(x^2)*cc*P_AA + (x^2)*P_AA*P_BB)/(-2*x*cc*P_AB + cc*(P_AB^2) + (x^2)*P_BB) -
                           (-cc*(P_AA + P_AB) + x*(cc - cc*o + o*P_BB))/(cc*(-1 + o) - o*P_BB)
                      })#,
                            # lower = ifelse(o>1, aa*P_AB/P_AA, 0), upper = ifelse(o<1, aa*P_AB/P_AA, P_AB))$root,
               # warning=function(warn){warning<<-1},
               # error=function(err){error<<-1})
        return(bb)
      })

      # if(warning+error==0){
        e = P_AB - b
        a = cr - b - c
        d = P_AA - a
      # }else{a=b=c=d=e=f=NA}
    }

    # if(is.na(a)==F & a>0 & b>0 & c>0 & d>0 & e>0 & f>0){
    if(all(!is.na(a) & a>0 & b>0 & c>0 & d>0 & e>0 & f>0)){
      return(or_calc(a,b,c,d,e,f, mod = True.Model, risk_allele))
        # if(risk_allele) return(OR = c(c*e/(b*f))[c(c*e/(b*f)) > 1])
        # if(!risk_allele) return(OR = c(c*e/(b*f))[c(c*e/(b*f)) < 1])
    }else{stop("Could not calculate odds ratio for the given parameters")} #before returning anything, check that no negative values in table
  # }else{stop("Could not calculate odds ratio for the given parameters")}
}

####################################################################################
# For an additive model
####################################################################################

#First need a helper functions
# sometimes R gives very small numbers instead of 0
# replace_0 <- function(anum){
#   ifelse(abs(anum) < 1e-4, return(0), return(anum))
# }
#Function to solve for a when b is known in an additive model; this should be moved outside the add.or function
solve_a<-function(b, cr, P_AA, P_AB, P_BB)
{
  a1<- P_AB*(2*b-P_AB)
  b1<- cr*(P_AB^2 - 2*P_AB*b)-b*P_AB^2 + (2*P_AB-P_AA+P_BB)*b^2
  c1<- -(b^2)*(P_BB-cr+b)*P_AA

  a.test<-quad_roots(a1,b1,c1)[1]
  a<-Re(polyroot(c(c1,b1,a1)))
  a <- a[a>0]
  a <- a[a<1]
  # if(length(a) == 2) a <- a[which.min(abs(a - a.test))]
  # a <- a[a>1e-5]
  # want to know if the two functions are giving different answers, and we know that a small a1 
  #  causes problems, so we want to know if something else becomes an issue
  # if(abs(a - a.test) > 1e-6 & a1 > 1e-4) warning("quad_roots and polyroot giving different answers")
  if(length(a) == 0) a <- NA
  return(a)
}

# xs<-seq(0,1,0.01)
# as<-sapply(xs, function(x){solve_a(x, cr, P_AA, P_AB, P_BB)})

#Find probability of disease in a dominant truth given AB or BB, additive test model
find.prob.dom<-function(x,m,cr, like)
{
  b<-2*m*(1-m)*x
  c<-m*m*x
  e<-2*m*(1-m)-b
  f<-m*m-c
  a<-cr-b-c
  d<-(1-m)*(1-m)-a

  tab <- rbind(c(a,b,c), c(d,e,f))
  mles<-optim(c(0,0),function(x) -additive.ll(x,tab), control=c(abstol = 0.00001))$par

  additive.ll(mles, tab)-like
}

#Find the probability of disease in a recessive truth given BB, additive test model
find.prob.rec<-function(x,m,cr, like)
{
  c<-m*m*x
  f<-m*m-c
  y<-(cr-c)/(2*m*(1-m)+(1-m)*(1-m))
  b<-2*m*(1-m)*y
  e<-2*m*(1-m)-b
  a<-(1-m)*(1-m)*y
  d<-(1-m)*(1-m)-a

  tab <- rbind(c(a,b,c), c(d,e,f))
  mles<-optim(c(0,0),function(x) -additive.ll(x,tab), control=c(abstol = 0.00001))$par

  additive.ll(mles, tab)-like
}

add.or.function <-function(like, Case.Rate, P_AA, P_AB, P_BB, True.Model, risk_allele)
{
  if(F){#to delete
    ############################################################################################################
    #log likelihood stats
    ############################################################################################################

    #For 1DF Test.Models the detectable LRT test statistic is:
    stat = ((qnorm(1-alpha/2)+qnorm(power))^2)/N

    #For the 2DF Test the detectable LRT test statistic is:
    stat_2df = uniroot(function(x) ncp.search(x, power, stat, alpha, df=2),
                   lower=0, upper=1000, extendInt = 'upX', tol=0.00001)$root/N


    ####################################################################################
    # Calculate the null and alternative likelihoods
    ####################################################################################

    # Can calculate the null log likelihood from the case.rate alone (intercept only model)
    beta0 <- logit(Case.Rate)
    ll.null <- Case.Rate*log(exp(beta0)/(1+exp(beta0))) + (1-Case.Rate)*log(1/(1+exp(beta0)))

    # Note that stat = 2*(as.numeric(ll.alt-ll.null))
    # Then the alt.likelihood to detect is
    ll.alt <- ll.null+(stat/2)
    like <- ll.alt

    ####################################################################################
    # Calculate the probabilities of each genotype
    ####################################################################################
    P_AA <- (1-m)^2
    P_AB <- 2*(1-m)*m
    P_BB <- m^2

    ####################################################################################
    # Assign alternate variable names
    ####################################################################################

    cr <- Case.Rate
    m <- MAF
  }
  cr <- Case.Rate
  # For a true additive model
  if (True.Model == "Additive"){
    # warning = error = 0 #For error handling

    # tryCatch(
      b <- ll_zero_finder2(function(x){
        a <- solve_a(x, cr, P_AA, P_AB, P_BB)
        # a <- a[length(a)]
        if(length(a) > 1){
          for(aa in a){
            d<-P_AA-aa; e<-P_AB-x; c<-cr-aa-x; f<-P_BB-c
            if(any( c( c( c(c,d,e,f) > 1), c( c(c,d,e,f) < 0) ) )) a <- a[a != aa]
          }
        }
        d<-P_AA-a; e<-P_AB-x; c<-cr-a-x; f<-P_BB-c

        ll <- a*log(a/P_AA)+
          x*log(x/P_AB)+
          c*log(c/P_BB)+
          d*log(d/P_AA)+
          e*log(e/P_AB)+
          f*log(f/P_BB)
        if(length(ll) == 0) ll <- NA
        return(like - ll)
      })#, 
      # lower=ifelse(risk_allele==T, max(c(0,cr-P_AA-P_BB, (cr-P_BB)*P_AB/(P_AA+P_AB))), max(c(0,cr-P_AA-P_BB,P_AB*(cr-P_AA)/(P_AB+P_BB)))),
      # upper=ifelse(risk_allele==T, min(P_AB,P_AB*cr/(P_AB+P_BB)),min(P_AB,P_AB*cr/(P_AB+P_AA)))
      # ), 
      # warning=function(warn){warning<<-1},
      # error=function(err){error<<-1}
      # )
    if(F){
      # tryCatch(b<-uniroot(function(x){
      #     a<-replace_0(solve_a(x, cr, P_AA, P_AB, P_BB))
      #     d<-replace_0(P_AA-a)
      #     e<-replace_0(P_AB-x)
      #     c<-replace_0(cr-a-x)
      #     f<-replace_0(P_BB-c)

      #     ll <- a*log(a/P_AA)+
      #       x*log(x/P_AB)+
      #       c*log(c/P_BB)+
      #       d*log(d/P_AA)+
      #       e*log(e/P_AB)+
      #       f*log(f/P_BB)
      #   return(like - ll)
      #   }, 
      #   lower=ifelse(risk_allele==T, max(c(0,cr-P_AA-P_BB, (cr-P_BB)*P_AB/(P_AA+P_AB)))+delta, max(c(0,cr-P_AA-P_BB,P_AB*(cr-P_AA)/(P_AB+P_BB)))+delta),
      #   upper=ifelse(risk_allele==T, min(P_AB,P_AB*cr/(P_AB+P_BB))-delta,min(P_AB,P_AB*cr/(P_AB+P_AA))-delta))$root,
      #   warning=function(warn){warning<<-1},
      #   error=function(err){error<<-1}
      #   )
    }

    # if(warning+error==0){
      # a <- solve_a(b, cr, P_AA, P_AB, P_BB)
      a <- lapply(b, solve_a, cr, P_AA, P_AB, P_BB)
      if(any(sapply(a, length) > 1)){
        # what to do if there are multiple zeros in solve_a:
        for(ix in which(sapply(a, length) > 1)){
          d<-P_AA-a[[ix]]
          e<-P_AB-b[ix]
          c<-cr-a[[ix]]-b[ix]
          f<-P_BB-c
          a[[ix]] <- a[[ix]][which(a[[ix]] > 0 & c > 0 & d > 0 & e > 0 & f > 0 & a[[ix]] < 1 & c < 1 & d < 1 & e < 1 & f < 1)]
          d<-P_AA-a[[ix]]
          e<-P_AB-b[ix]
          c<-cr-a[[ix]]-b[ix]
          f<-P_BB-c
          a[[ix]] <- a[[ix]][which.min(
            a[[ix]]*log(a[[ix]]/P_AA)+b[ix]*log(b[ix]/P_AB)+c*log(c/P_BB)+d*log(d/P_AA)+e*log(e/P_AB)+f*log(f/P_BB) - like)]
        }
        a <- unlist(a)
        if(length(a) != length(b)) stop("tried to fix multiple zeros in solve_a, but there is still a problem")
      }else a <- unlist(a)
      d<-P_AA-a
      e<-P_AB-b
      c<-cr-a-b
      f<-P_BB-c
    # }else{a<-b<-c<-d<-e<-f<-NA}
  }

  #Dominant
  if (True.Model == "Dominant"){
    # warning = error = 0 #For error handling

    #x = prob of disease given AB or BB
    # tryCatch(
    x <- ll_zero_finder2(function(x) find.prob.dom(x, m, cr, like=like))#,lower = ifelse(risk_allele==T, cr,0)+delta, upper = ifelse(risk_allele==T, 1,cr)-delta)$root,
             # warning=function(warn){warning<<-1},
             # error=function(err){error<<-1})

    # if(warning+error==0){
        a = cr-x*P_AB-x*P_BB
        b = x*P_AB
        c = x*P_BB
        d = P_AA - a
        e = P_AB - b
        f = P_BB - c
    # }else{a<-b<-c<-d<-e<-f<-NA}
  }
  #Recessive
  if(True.Model=='Recessive'){
    #Let x = prob of disease given BB
    # warning = error = 0 #For error handling

    # tryCatch(
    x <- ll_zero_finder2(function(x) find.prob.rec(x, m, cr, like=like))#,lower = cr, upper = 1)$root,
             # warning=function(warn){warning<<-1},
             # error=function(err){error<<-1})

    # if(warning+error==0){
      y<-(cr-x*P_BB)/(P_AA+P_AB)
      a = P_AA*y
      b = P_AB*y
      c = x*P_BB
      d = P_AA - a
      e = P_AB - b
      f = P_BB - c
    # }else{a<-b<-c<-d<-e<-f<-NA}
  }
  if(all(!is.na(a) & a>0 & b>0 & c>0 & d>0 & e>0 & f>0)){
    # if(risk_allele) res <- c(OR1 = c(b*d/(a*e))[b*d/(a*e) > 1], OR2 = c(c*d/(a*f))[c*d/(a*f) > 1])
    # if(!risk_allele) res <- c(OR1 = c(b*d/(a*e))[b*d/(a*e) < 1], OR2 = c(c*d/(a*f))[c*d/(a*f) < 1])
    return(or_calc(a,b,c,d,e,f, mod = True.Model, risk_allele))
  }else{stop("Could not calculate odds ratio for the given parameters")}
}

####################################################################################
# For a 2DF model
####################################################################################

or.function.2df <-function(power, N, alpha = 0.05, Case.Rate, MAF, True.Model, risk_allele)
{
  #For 1DF Test.Models the detectable LRT test statistic is:
  stat = ((qnorm(1-alpha/2)+qnorm(power))^2)/N

  #For the 2DF Test the detectable LRT test statistic is:
  stat_2df = uniroot(function(x) ncp.search(x, power, stat, alpha, df=2),
                 lower=0, upper=1000, extendInt = 'upX', tol=0.00001)$root/N

  # Can calculate the null log likelihood from the case.rate alone (intercept only model)
  cr = Case.Rate
  m = MAF
  beta0 <- logit(Case.Rate)
  ll.null <- Case.Rate*log(exp(beta0)/(1+exp(beta0))) + (1-Case.Rate)*log(1/(1+exp(beta0)))

  # Note that stat = 2*(as.numeric(ll.alt-ll.null))
  # Then the alt.likelihood to detect is
  ll.alt <- ll.null+(stat_2df/2)
  like <- ll.alt
  P_AA <- (1-m)^2
  P_AB <- 2*(1-m)*m
  P_BB <- m^2


  # If model is dominant, will estimate the same coef for AB and BB
  # Solve for a, the joint prob of geno AA and disease
  if (True.Model == "Dominant"){
    warning<-error<-0
    # tryCatch(
      a <- ll_zero_finder2(function(x) {x*log(x/P_AA)+
          (cr-x)*log((cr-x)/(P_AB+P_BB))+
          (P_AA-x)*log((P_AA-x)/P_AA)+
          (1-cr-P_AA+x)*log((1-cr-P_AA+x)/(P_AB+P_BB)) - like})#,
          # lower=max(cr+P_AA-1+0.0000000001, 0.0000000001),upper=cr*P_AA),
          #this assumes that m is a risk allele/that OR's will be > 0; Consider adding in option for protective alleles
        # warning=function(warn){warning<<-1},
        # error=function(err){error<<-1})


    if(warning+error==0){
      d <- P_AA - a
      y <- (cr-a)/(P_AB+P_BB)
      b <- y*P_AB
      c <- y*P_BB
      e <- P_AB - b
      f <- P_BB - c
    }
  }

  # If model is recessive
  if (True.Model == "Recessive"){
    # warning<-error<-0
      # tryCatch(
        c <- ll_zero_finder2(function(x){
            x*log(x/P_BB)+
            (cr-x)*log((cr-x)/(P_AB+P_AA))+
            (P_BB-x)*log((P_BB-x)/P_BB)+
            (1-cr-P_BB+x)*log((1-cr-P_BB+x)/(P_AB+P_AA)) - like})#,
            # lower=cr*P_BB, upper=P_BB-0.0000000001),#this assumes that m is a risk allele
        # warning=function(warn){warning<<-1},
        # error=function(err){error<<-1})

    # if(warning+error==0){
      f <- P_BB - c
      y = (cr-c)/(P_AB+P_AA)
      a = y*P_AA
      d = P_AA - a
      b = y*P_AB
      e = P_AB - b
    # }
  }

  # If model is additive
  if(True.Model=='Additive'){
    solve_a<-function(b){
      a1<- P_AB*(2*b-P_AB)
      b1<- cr*(P_AB^2 - 2*P_AB*b)-b*P_AB^2 + (2*P_AB-P_AA+P_BB)*b^2
      c1<- -(b^2)*(P_BB-cr+b)*P_AA

      a<-quad_roots(a1,b1,c1)[1]

      d<-P_AA-a
      e<-P_AB-b
      c<-cr-a-b
      f<-P_BB-c

      ll <- a*log(a/P_AA)+
        b*log(b/P_AB)+
        c*log(c/P_BB)+
        d*log(d/P_AA)+
        e*log(e/P_AB)+
        f*log(f/P_BB)
        return(like - ll)
    }
    #b<-uniroot(solve_a, lower=max(0,cr-P_AA)+0.00001, upper=min(P_AB,cr-P_BB)-0.000001 ) #2 rootsin this window

    b<-uniroot(solve_a, lower=P_AB*cr+0.0000000001, upper=min(P_AB,cr-P_BB)-0.000001 )$root #this gives the root for an OR>1

    a1<- P_AB*(2*b-P_AB)
    b1<- cr*(P_AB^2 - 2*P_AB*b)-b*P_AB^2 + (2*P_AB-P_AA+P_BB)*b^2
    c1<- -(b^2)*(P_BB-cr+b)*P_AA
    a<-quad_roots(a1,b1,c1)[1]
    d<-P_AA-a
    e<-P_AB-b
    c<-cr-a-b
    f<-P_BB-c
  }
  if(all(!is.na(a) & a>0 & b>0 & c>0 & d>0 & e>0 & f>0)){
    # if(risk_allele) res <- c(OR1 = c(b*d/(a*e))[b*d/(a*e) > 1], OR2 = c(c*d/(a*f))[c*d/(a*f) > 1])
    # if(!risk_allele) res <- c(OR1 = c(b*d/(a*e))[b*d/(a*e) < 1], OR2 = c(c*d/(a*f))[c*d/(a*f) < 1])
    return(or_calc(a,b,c,d,e,f, mod = True.Model, risk_allele))
  }else{stop("Could not calculate odds ratio for the given parameters")}
}

# example
# effect_size_function(N=N, Case.Rate = Case.Rate, MAF = MAF, power = power, risk_allele = risk_allele, Alpha = 0.05, True.Model = True.Model, Test.Model = Test.Model)

effect_size_function <-
  function(N=NULL, Case.Rate=NULL, k=NULL, MAF=NULL, power=NULL, risk_allele = T,
                     Alpha=0.05, True.Model='All', Test.Model='All')
{

  ############################################################################################################
  #Error Messages for insufficient sample size information, MAF, and case vs. control ratio
  ############################################################################################################
  if(is.null(N)) {
    stop("N, the total sample size, must be specified.")}

  if(is.null(k)& is.null(Case.Rate)){
    stop("k, the number of controls per case, or Case.Rate, the proportion of cases in the study sample, must be specified.")
  }

  if(is.null(MAF)){
    stop("MAF (minor allele frequency) must be specified.")
  }

  if(is.null(power)==T){
    stop("power (statistical power) must be specified.")
  }

  if(sum(Case.Rate>=1)>0 | sum(Case.Rate<=0)>0){
    stop("R2 must be greater than 0 and less than 1.")
  }
  if(sum(power>=1)>0 | sum(power<=0)>0){
    stop("power must be greater than 0 and less than 1.")
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

  if(sum(Alpha>=1)>0 | sum(Alpha<=0)>0){
    stop("Alpha must be greater than 0 and less than 1.")
  }

  if(sum(!(Test.Model %in% c("Dominant", "Recessive", "Additive", "2df", "All")))>0){
    stop(paste("Invalid Test.Model:",
               paste(Test.Model[!(Test.Model %in% c("Dominant", "Recessive", "Additive", "2df", "All"))], collapse=', ')))
  }

  if(sum(!(True.Model %in% c("Dominant", "Recessive", "Additive", "All")))>0){
    stop(paste("Invalid True.Model:",
               paste(True.Model[!(True.Model %in% c("Dominant", "Recessive", "Additive", "All"))], collapse=', ')))
  }




  ####################################################################################
  # Next, solve for the OR or table that would correspond to this alternative likelihood
  # May need to use an optimizer or root finding function, depending on model
  ####################################################################################
  # Consider a 2x3 table
  #         AA    AB    BB    TOTAL
  # Case    a     b     c     cr
  # Control d     e     f     1-cr
  # TOTAL   P_AA  P_AB  P_BB  1
  # Need to solve for joint probabilities a-f and calculate the detectable OR
  ####################################################################################
  # delta = 0.0000000001 #I'm not sure how to set this.  It is to nudge the limits of the uniroot function so we don't get NaN's

  ############################################################################################################
  #Calculate needed sample size information from provided inputs
  ############################################################################################################
  #Test model vector
  if('All' %in% Test.Model){Test.Model<-c("Dominant", "Recessive", "Additive", "2df")}

  #True model vector
  if('All' %in% True.Model){True.Model<-c("Dominant", "Recessive", "Additive1", "Additive2")}

  #If k is provided calculate the Case.Rate
  if(is.null(Case.Rate)==T){Case.Rate = 1/(1+k)}

  #Find all possible combinations of N and case_rate
  sample.size.tab <- expand.grid(N, Case.Rate)
  colnames(sample.size.tab) <- c('N', 'Case.Rate')
  sample.size.tab$N_cases <- sample.size.tab$N*sample.size.tab$Case.Rate
  sample.size.tab$N_controls <- sample.size.tab$N - sample.size.tab$N_cases

  iter <- nrow(sample.size.tab)

  final.or.tab <- NULL

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
    pow.save.tab <-NULL


    for (pow in power){
      for(alpha0 in Alpha){

        ############################################################################################################
        #log likelihood stats
        ############################################################################################################

        stat = ((qnorm(1-alpha0/2)+qnorm(pow))^2)/N

        #For the 2DF Test the detectable LRT test statistic is:
        stat_2df = uniroot(function(x) ncp.search(x, pow, stat, alpha0, df=2),
                       lower=0, upper=1000, extendInt = 'upX', tol=0.00001)$root/N

        beta0 <- logit(Case.Rate)
        ll.null <- Case.Rate*log(exp(beta0)/(1+exp(beta0))) + (1-Case.Rate)*log(1/(1+exp(beta0)))

        # Note that stat = 2*(as.numeric(ll.alt-ll.null))
        # Then the alt.likelihood to detect is
        ll.alt <- ll.null+(stat/2)
        like <- ll.alt

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

          if('Dominant' %in% Test.Model){
            o <- dom.or.function(like=like, Case.Rate=Case.Rate, P_AA=P_AA, P_AB=P_AB, P_BB=P_BB, True.Model=True.Model, risk_allele=risk_allele)
          }

          # if('Additive1' %in% Test.Model){
          # }

          if('Additive' %in% Test.Model){
            o <- add.or.function(plike=like, Case.Rate=Case.Rate, P_AA=P_AA, P_AB=P_AB, P_BB=P_BB, True.Model=True.Model, risk_allele=risk_allele)
          }

          if('Recessive' %in% Test.Model){
            o <- rec.or.function(like=like, Case.Rate=Case.Rate, P_AA=P_AA, P_AB=P_AB, P_BB=P_BB, True.Model=True.Model, risk_allele=risk_allele)
          }

          m.save.tab<-rbind(m.save.tab,
                            data.frame(Test.Model = Test.Model, True.Model = True.Model, MAF=m, power = pow, N_total = N, N_cases = N_cases, N_controls = N_controls, 
                                       Case.Rate = Case.Rate, OR = o, Alpha = alpha0))
            # Test.Model True.Model MAF OR N_total N_cases N_controls Case.Rate Power_at_Alpha_0.05
        }
        pow.save.tab<-rbind(pow.save.tab, m.save.tab)
      }
    }

    ############################################################################################################
    #Calculate Odds Ratio for each scenario under the specified testing model
    ############################################################################################################

    final.or.tab<-rbind(final.or.tab, pow.save.tab)

    # c("Test.Model", "True.Model", "MAF", "OR", "N_total", "N_cases", "N_controls", "Case.Rate", "Power_at_Alpha_0.05")
    # colnames(final.or.tab)<-c('Test.Model', 'True.Model', 'MAF', 'OR', 'N_total', 'N_cases', 'N_controls','Case.Rate',
    #                        paste("Power_at_Alpha_", Alpha, sep=''))
    # c("Test.Model", "True.Model", "MAF", "OR", "N_total", "N_cases", "N_controls", "Case.Rate", "Alpha")
  }
  # fix table so that it's in the same format as Camille's
  final.or.tab2 <- unique(final.or.tab[,1:8])
  for(analpha in unique(final.or.tab$Alpha)){
    final.or.tab2 <- cbind(final.or.tab2, NA)
    names(final.or.tab2)[ncol(final.or.tab2)] <- paste0("OR_at_Alpha_", analpha)
    final.or.tab2[,ncol(final.or.tab2)] <- final.or.tab[final.or.tab$Alpha == analpha, "OR"]
  }
  return(final.or.tab2)
}