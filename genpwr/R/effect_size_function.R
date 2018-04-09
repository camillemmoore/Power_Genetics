#Should make this into a function called effect_size_function

#GOAL: taking input of vectors of power, alpha, sample size, MAF and case rate, true model and test model,
# return a table of detectable odds ratios
# user should specify if m is a risk allele or protective allele (OR's >1 or <1)

####################################################################################
# First identify the detectable LRT test statistic
#####################################################################################
#For 1DF Test.Models the detectable LRT test statistic is:

stat = ((qnorm(1-alpha/2)+qnorm(power))^2)/N

#For the 2DF Test the detectable LRT test statistic is:
stat = uniroot(function(x) ncp.search(x, power, stat, alpha, df=2),
               lower=0, upper=1000, extendInt = 'upX', tol=0.00001)$root/N


####################################################################################
# Calculate the null and alternative likelihoods
####################################################################################

# Can calculate the null log likelihood from the case.rate alone (intercept only model)
  beta0 = logit(Case.Rate)
  ll.null = Case.Rate*log(exp(beta0)/(1+exp(beta0))) + (1-Case.Rate)*log(1/(1+exp(beta0)))

# Note that stat = 2*(as.numeric(ll.alt-ll.null))
# Then the alt.likelihood to detect is
  ll.alt = ll.null+(stat/2)

####################################################################################
# Calculate the probabilities of each genotype
####################################################################################
  P_AA = (1-m)^2
  P_AB = 2*(1-m)*m
  P_BB = m^2

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
delta = 0.0000000001 #I'm not sure how to set this.  It is to nudge the limits of the uniroot function so we don't get NaN's
####################################################################################
# For a dominant model
####################################################################################
dom.or.function <-function(like, cr, P_AA, P_AB, P_BB, True.Model, risk_allele, delta){

  # Solve for a, the joint prob of geno AA and disease;
  # a and d will be the same for all true genetic models
    warning = error = 0 #For error handling

    tryCatch(a <- uniroot(function(x) {x*log(x/P_AA)+
          (cr-x)*log((cr-x)/(P_AB+P_BB))+
          (P_AA-x)*log((P_AA-x)/P_AA)+
          (1-cr-P_AA+x)*log((1-cr-P_AA+x)/(P_AB+P_BB)) - like},
        lower=ifelse(risk_allele==T, max(cr+P_AA-1+delta, delta), P_AA*cr),
        upper=ifelse(risk_allele==T, cr*P_AA-delta, min(P_AA, cr)-delta))$root, #this assumes that m is a risk allele/that OR's will be > 0; Consider adding in option for protective alleles
      warning=function(warn){warning<<-1},
      error=function(err){error<<-1})


  if(warning+error==0){
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
      warning<-error<-0

      # Calculate the OR under the dominant model
      o = (cr-a)*d/(a*(1-cr-d))

    #Here to solve for b, solve for c in terms of b in two ways, first using the dominant OR relationship and second setting OR(AB vs. AA) equal to OR(BB vs AB)
    #subtract the two version of c from one another and find the root of the equation
    tryCatch(b <- uniroot(function(x) (a*(x-x*o+o*(P_BB+ P_AB)) - x*P_AA)/(a*(-1 + o) + P_AA) - (P_BB*(P_AA-a)*x^2)/(a*P_AB*(P_AB-2*x)+P_AA*x^2)
                        ,lower = ifelse(o>1, a*P_AB/P_AA, 0), upper = ifelse(o<1, a*P_AB/P_AA, P_AB))$root,
             warning=function(warn){warning<<-1},
             error=function(err){error<<-1})

      if(warning+error==0){
      c = cr - a - b
      e = P_AB - b
      f = P_BB - c
      }else{a=b=c=d=e=f=NA}
}

    #For all models, calc the OR's
    #tab <- rbind(c(a,b,c), c(d,e,f)) #table
    if(is.na(a)==F & a>0 & b>0 & c>0 & d>0 & e>0 & f>0){
    return(c(OR1 = b*d/(a*e), OR2 = c*d/(a*f)))
      }else{return(c(OR1 = NaN, OR2 = NaN))} #before returning anything, check that no negative values in table
  }else{return(c(OR1 = NaN, OR2 = NaN))}
}

####################################################################################
# For a recessive model
####################################################################################
rec.or.function <-function(like, cr, P_AA, P_AB, P_BB, True.Model, risk_allele, delta){

  # Solve for c, the prob of BB and disease
  warning<-error<-0
  tryCatch(
    c <- uniroot(function(x) {x*log(x/P_BB)+
          (cr-x)*log((cr-x)/(P_AB+P_AA))+
          (P_BB-x)*log((P_BB-x)/P_BB)+
          (1-cr-P_BB+x)*log((1-cr-P_BB+x)/(P_AB+P_AA)) - like},
        lower=ifelse(risk_allele==T, cr*P_BB+delta, delta),
        upper=ifelse(risk_allele==T, P_BB-delta, cr*P_BB - delta))$root,
    warning=function(warn){warning<<-1},
    error=function(err){error<<-1})

  if(warning+error==0){
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

    # Solve for the joint probabilities
      b = P_AB_case*P_AB
      d = P_AB - b
      a = cr - b - c
      d = P_AA - a
    }

    if(True.Model=='Additive'){
    # Calculate the OR for the recessive model
      o = c*(1-cr-f)/((cr-c)*f)

    # To find b, solve for a in two ways (see dominant.or.function)
      warning<-error<-0
        tryCatch(b <- uniroot(function(x) ((-(x^2)*c*P_AA + (x^2)*P_AA*P_BB)/(-2*x*c*P_AB + c*(P_AB^2) + (x^2)*P_BB)) -
                           ((-c*(P_AA + P_AB) + x*(c - c*o + o*P_BB))/(c*(-1 + o) - o*P_BB))
                         ,lower = ifelse(o>1,0, c*P_AB/P_BB), upper = ifelse(o<1, P_AB,c*P_AB/P_BB))$root,
             warning=function(warn){warning<<-1},
             error=function(err){error<<-1})


    if(warning+error==0){
      e = P_AB - b
      a = cr - b - c
      d = P_AA - a
    }else{a=b=c=d=e=f=NA}
    }

    if(is.na(a)==F & a>0 & b>0 & c>0 & d>0 & e>0 & f>0){
        return(c(OR1 = b*d/(a*e), OR2 = c*d/(a*f)))
      }else{return(c(OR1 = NaN, OR2 = NaN))} #before returning anything, check that no negative values in table
    }else{return(c(OR1 = NaN, OR2 = NaN))}
  }

####################################################################################
# For an additive model
####################################################################################

#First need a helper functions
  #Function to solve for a when b is known in an additive model; this should be moved outside the add.or function
  solve_a<-function(b, cr, P_AA, P_AB, P_BB){
    a1<- P_AB*(2*b-P_AB)
    b1<- cr*(P_AB^2 - 2*P_AB*b)-b*P_AB^2 + (2*P_AB-P_AA+P_BB)*b^2
    c1<- -(b^2)*(P_BB-cr+b)*P_AA

    a<-quad_roots(a1,b1,c1)[1]
    return(a)}

  xs<-seq(0,1,0.01)
  as<-mapply(function(x){solve_a(x, cr, P_AA, P_AB, P_BB)}, xs)

  #Find probability of disease in a dominant truth given AB or BB, additive test model
  find.prob.dom<-function(x,m,cr, like){
    b<-2*m*(1-m)*x
    c<-m*m*x
    e<-2*m*(1-m)-b
    f<-m*m-c
    a<-cr-b-c
    d<-(1-m)*(1-m)-a

    tab <- rbind(c(a,b,c), c(d,e,f))
    mles<-optim(c(0,0),function(x) -additive.ll(x,tab), control=c(abstol = 0.00001))$par

    additive.ll(mles, tab)-like}

  #Find the probability of disease in a recessive truth given BB, additive test model
  find.prob.rec<-function(x,m,cr, like){
    c<-m*m*x
    f<-m*m-c
    y<-(cr-c)/(2*m*(1-m)+(1-m)*(1-m))
    b<-2*m*(1-m)*y
    e<-2*m*(1-m)-b
    a<-(1-m)*(1-m)*y
    d<-(1-m)*(1-m)-a

    tab <- rbind(c(a,b,c), c(d,e,f))
    mles<-optim(c(0,0),function(x) -additive.ll(x,tab), control=c(abstol = 0.00001))$par

    additive.ll(mles, tab)-like}

add.or.function <-function(like, cr, P_AA, P_AB, P_BB, True.Model, risk_allele, delta){
  # For a true additive model
  if (True.Model == "Additive"){
  warning = error = 0 #For error handling

  tryCatch(b<-uniroot(function(x){
      a <- solve_a(x, cr, P_AA, P_AB, P_BB)
      d<-P_AA-a
      e<-P_AB-x
      c<-cr-a-x
      f<-P_BB-c

      ll <- a*log(a/P_AA)+
        x*log(x/P_AB)+
        c*log(c/P_BB)+
        d*log(d/P_AA)+
        e*log(e/P_AB)+
        f*log(f/P_BB)
      return(like - ll)
      }, lower=ifelse(risk_allele==T, max(c(0,cr-P_AA-P_BB, (cr-P_BB)*P_AB/(P_AA+P_AB)))+delta, max(c(0,cr-P_AA-P_BB,P_AB*(cr-P_AA)/(P_AB+P_BB)))+delta),
      upper=ifelse(risk_allele==T, min(P_AB,P_AB*cr/(P_AB+P_BB))-delta,min(P_AB,P_AB*cr/(P_AB+P_AA))-delta))$root,
      warning=function(warn){warning<<-1},
      error=function(err){error<<-1})

  if(warning+error==0){
  a <- solve_a(b, cr, P_AA, P_AB, P_BB)
  d<-P_AA-a
  e<-P_AB-b
  c<-cr-a-b
  f<-P_BB-c
  }else{a<-b<-c<-d<-e<-f<-NA}}

  #Dominant
  if (True.Model == "Dominant"){
    warning = error = 0 #For error handling

    #x = prob of disease given AB or BB
    tryCatch(x <- uniroot(function(x) find.prob.dom(x, m, cr, like=like),lower = ifelse(risk_allele==T, cr,0)+delta, upper = ifelse(risk_allele==T, 1,cr)-delta)$root,
             warning=function(warn){warning<<-1},
             error=function(err){error<<-1})

    if(warning+error==0){
        a = cr-x*P_AB-x*P_BB
        b = x*P_AB
        c = x*P_BB
        d = P_AA - a
        e = P_AB - b
        f = P_BB - c
    }else{a<-b<-c<-d<-e<-f<-NA}}

#Recessive
  if(True.Model=='Recessive'){
    #Let x = prob of disease given BB
    warning = error = 0 #For error handling

    tryCatch(x <- uniroot(function(x) find.prob.rec(x, m, cr, like=like),lower = cr, upper = 1)$root,
             warning=function(warn){warning<<-1},
             error=function(err){error<<-1})

    if(warning+error==0){
      y<-(cr-x*P_BB)/(P_AA+P_AB)
      a = P_AA*y
      b = P_AB*y
      c = x*P_BB
      d = P_AA - a
      e = P_AB - b
      f = P_BB - c
    }else{a<-b<-c<-d<-e<-f<-NA}}

  if(is.na(a)==F & a>0 & b>0 & c>0 & d>0 & e>0 & f>0){
    return(c(OR1 = b*d/(a*e), OR2 = c*d/(a*f)))
  }else{return(c(OR1 = NaN, OR2 = NaN))}
}


####################################################################################
# For a 2DF model
####################################################################################

P_AA = (1-m)^2
P_AB = 2*(1-m)*m
P_BB = m^2

# If model is dominant, will estimate the same coef for AB and BB
# Solve for a, the joint prob of geno AA and disease
tryCatch(
  a = uniroot(function(x) {x*log(x/P_AA)+
      (cr-x)*log((cr-x)/(P_AB+P_BB))+
      (P_AA-x)*log((P_AA-x)/P_AA)+
      (1-cr-P_AA+x)*log((1-cr-P_AA+x)/(P_AB+P_BB)) - like},
      lower=max(cr+P_AA-1+0.0000000001, 0.0000000001),
      upper=cr*P_AA)$root, #this assumes that m is a risk allele/that OR's will be > 0; Consider adding in option for protective alleles
  warning=function(warn){warning<<-1},
  error=function(err){error<<-1})


if(warning+error==0){
  d <- P_AA - a
  y <- (cr-a)/(P_AB+P_BB)
  b <- y*P_AB
  c <- y*P_BB
  e <- P_AB - b
  f <- P_BB - c}

# If model is recessive
warning<-error<-0
tryCatch(
  c <- uniroot(function(x) {x*log(x/P_BB)+
      (cr-x)*log((cr-x)/(P_AB+P_AA))+
      (P_BB-x)*log((P_BB-x)/P_BB)+
      (1-cr-P_BB+x)*log((1-cr-P_BB+x)/(P_AB+P_AA)) - like},
      lower=cr*P_BB, upper=P_BB-0.0000000001)$root,#this assumes that m is a risk allele
  warning=function(warn){warning<<-1},
  error=function(err){error<<-1})

if(warning+error==0){
  f <- P_BB - c
  y = (cr-c)/(P_AB+P_AA)
  a = y*P_AA
  d = P_AA - a
  b = y*P_AB
  e = P_AB - b}

# If model is additive
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



Pnorm(sqrt(N * (p1 – p2)^2 * (1 – B)) – Za * sqrt(p*(1-p)/B)/ sqrt(p1 * (1 - p1) + p2 * (1 - p2) * (1 - B)/B))

