#Should make this into a function called effect_size_function

#GOAL: taking input of vectors of power, sample size, MAF and case rate, true model and test model,
# return a table of detectable odds ratios

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
# Next, solve for the OR or table that would correspond to this alternative likelihood
# May need to use an optimizer or root finding function, depending on model
####################################################################################

# First code for when true model and test model match:
  # For a dominant model + dominant truth
    dom.or.function <-function(like, m, cr){
      # Calculate the probabilities of each genotype
        P_AA = (1-m)^2
        P_AB = 2*(1-m)*m
        P_BB = m^2

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

  #if true model == dominant
      y <- (cr-a)/(P_AB+P_BB)
      b <- y*P_AB
      c <- y*P_BB
      e <- P_AB - b
      f <- P_BB - c

  #if true model == recessive
    # The conditional probability of disease in AA and AB should be the same
      P_AB_case = P_AA_case = a/P_AA

    # Now solve for the elements of the recessive table
      b = P_AB_case*P_AB
      c = cr - a - b
      e = P_AB - b
      f = P_BB - c

  #if(true.model=='Additive1'|true.model=='Additive2')
      warning<-error<-0

      # Calculate the OR under the dominant model
      o = (cr-a)*d/(a*(1-cr-d))

      #Here to solve for b, solve for c in terms of b in two ways, first using the dominant OR relationship and second setting OR(AB vs. AA) equal to OR(BB vs AB)
      #subtract the two version of c from one another and find the root of the equation
      tryCatch(b<-uniroot(function(x) (a*(x-x*o+o*(P_BB+ P_AB)) - x*P_AA)/(a*(-1 + o) + P_AA) - (P_BB*(P_AA-a)*x^2)/(a*P_AB*(P_AB-2*x)+P_AA*x^2)
                             ,lower = ifelse(o>1, a*P_AB/P_AA, 0), upper = ifelse(o<1, a*P_AB/P_AA, P_AB))$root,
                 warning=function(warn){warning<<-1},
                 error=function(err){error<<-1})


  #For all models, save the table and calc the OR's
  tab <- rbind(c(a,b,c), c(d,e,f)) #table
  return(list(tab = tab, OR1 = b*d/(a*e), OR2 = c*d/(a*f))) #before returning anything, check that no negative values in table
  }else{return(list(tab = NaN, OR1 = NaN, OR2 = NaN))}
    }

  # For a recessive model + recessive truth
    rec.or.function <-function(like, m, cr){
      # Calculate Genotype Probabilities
      P_AA = (1-m)^2
      P_AB = 2*(1-m)*m
      P_BB = m^2

      # Solve for c, the prob of BB and disease
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

  # If the true model is recessive:
    y = (cr-c)/(P_AB+P_AA)
    a = y*P_AA
    d = P_AA - a
    b = y*P_AB
    e = P_AB - b

  # If the true model is dominant:
    # the conditional prob of disease for genotypes AB and BB are equal
      P_AB_case = P_BB_case= c/P_BB

    # Solve for the joint probabilities
      b = P_AB_case*P_AB
      d = P_AB - b
      a = cr - b - c
      d = P_AA - a

  # If the true model is additive:
    # Calculate the OR for the recessive model
      o = c*(1-cr-f)/((cr-c)*f)

      # To find b, solve for a in two ways (see dominant.or.function)
      warning<-error<-0
      tryCatch(b = uniroot(function(x) ((-(x^2)*c*P_AA + (x^2)*P_AA*P_BB)/(-2*x*c*P_AB + c*(P_AB^2) + (x^2)*P_BB)) -
                             ((-c*(P_AA + P_AB) + x*(c - c*o + o*P_BB))/(c*(-1 + o) - o*P_BB))
                           ,lower = ifelse(o>1,0, c*P_AB/P_BB), upper = ifelse(o<1, P_AB,c*P_AB/P_BB))$root,
               warning=function(warn){warning<<-1},
               error=function(err){error<<-1})


      if(warning+error==0){
        e = P_AB - b
        a = cr - b - c
        d = P_AA - a}

  tab <- rbind(c(a,b,c), c(d,e,f)) #table
  return(list(tab = tab, OR = c*d/(a*f)))
  }else{return(list(tab = NaN, OR = NaN))}
}

  # For an additive model + additive truth
    add.or.function <-function(like, m, cr){
      # Calculate genotype probabilities
        P_AA = (1-m)^2
        P_AB = 2*(1-m)*m
        P_BB = m^2

      # For a true additive model

      #Function to solve for a when b is known in an additive model; this should be moved outside the add.or function
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

      #Dominant

      #Let x = prob of disease given AB or BB
      #Then
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

      root<-uniroot(function(x) find.prob.dom(x, m, cr, like=like),lower = cr, upper = 1)$root

      dom.tab_ma2<-data.frame(model=rep('Dominant',2),table=rbind(c(Case.Rate-root*P_AB-root*P_BB, root*P_AB, root*P_BB),
                                                                  c(P_AA-(Case.Rate-root*P_AB-root*P_BB), P_AB-root*P_AB,P_BB-root*P_BB)))


      dom_OR_ma2 <- dom.tab_ma2[2,2]*(dom.tab_ma2[1,3]+dom.tab_ma2[1,4])/(dom.tab_ma2[2,3]+dom.tab_ma2[2,4])/dom.tab_ma2[1,2]
    }
    #Recessive if(true.model=='Recessive'){
      #Let x = prob of disease given BB
      #Then
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

      root<-uniroot(function(x) find.prob.rec(x, m, cr, like=like),lower = cr, upper = 1)$root

      y<-(Case.Rate-root*P_BB)/(P_AA+P_AB)

      rec.tab_ma2<-data.frame(model=rep('Recessive',2),table=rbind(c(P_AA*y, P_AB*y, root*P_BB),
                                                                   c((1-P_AA)*y, (1-P_AB)*y,P_BB-root*P_BB)))



      rec_OR_ma2 <- rec.tab_ma2[2,2]*(rec.tab_ma2[1,3]+rec.tab_ma2[1,4])/(rec.tab_ma2[2,3]+rec.tab_ma2[2,4])/rec.tab_ma2[1,2]



    tab <- rbind(c(a,b,c), c(d,e,f)) #table
    return(list(tab = tab, OR1 = b*d/(a*e), OR2 = c*d/(a*f)))
  }


# For 2 DF Test
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



# Now consider model mis-specification




    # Additive Test, Dominant Truth
    # Additive Test, Recessive Truth
    # Recessive Test, Additive Truth

    # 2df Test, Dominant Truth
    # 2df Test, Additive Truth
    # 2df Test, Recessive Truth





########################################################################################
# Suppose Using Dominant Model.  What OR would be needed for a significant result
# under other models?
########################################################################################
P_AA = (1-m)^2
P_AB = 2*(1-m)*m
P_BB = m^2

prob_AA_case = tab.assuming.test.true[1,1]
prob_AB_case = tab.assuming.test.true[1,2]
prob_BB_case = tab.assuming.test.true[1,3]
prob_AA_control = tab.assuming.test.true[2,1]
prob_AB_control = tab.assuming.test.true[2,2]
prob_BB_control = tab.assuming.test.true[2,3]

if(test.model=='Dominant'){
  P_AA_case = prob_AA_case/P_AA

  if(true.model=='Recessive'){
    #Recessive Truth
      P_AB_case_r_md<-P_AA_case
      prob_AB_case_r_md<-(P_AB_case_r_md*P_AB)
      prob_AB_control_r_md<-(1-P_AB_case_r_md)*P_AB
      prob_BB_case_r_md <- (prob_AB_case+prob_BB_case)-prob_AB_case_r_md
      prob_BB_control_r_md <- P_BB-prob_BB_case_r_md

      rec.tab_md<-data.frame(model=rep('Recessive',2),table=rbind(c(prob_AA_case, prob_AB_case_r_md, prob_BB_case_r_md),
                                                            c(prob_AA_control, P_AB-prob_AB_case_r_md,P_BB-prob_BB_case_r_md)))

if(prob_BB_case_r_md<0 | prob_BB_control_r_md<0) {rec_OR_md <- NaN
}else {rec_OR_md <- rec.tab_md[1,4]*(rec.tab_md[2,3]+rec.tab_md[2,2])/(rec.tab_md[1,3]+rec.tab_md[1,2])/rec.tab_md[2,4]
}
}
#Additive
if(true.model=='Additive1'|true.model=='Additive2'){
warning<-error<-0
tryCatch(zz<-uniroot(function(x) (prob_AA_case*(x-x*o+o*(P_BB+ P_AB)) - x*P_AA)/(prob_AA_case*(-1 + o) + P_AA) - (P_BB*(P_AA-prob_AA_case)*x^2)/(prob_AA_case*P_AB*(P_AB-2*x)+P_AA*x^2)
                     ,lower = 0, upper = 1)$root,
         warning=function(warn){warning<<-1},
         error=function(err){error<<-1})


if(warning+error==0){
  prob_AB_case_a_md<-zz
  P_AB_case_a_md<-prob_AB_case_a_md/P_AB
  prob_AB_control_a_md<-(1-P_AB_case_a_md)*P_AB
  prob_BB_case_a_md <- (prob_AB_case+prob_BB_case)-prob_AB_case_a_md
  prob_BB_control_a_md <- P_BB-prob_BB_case_a_md

  add.tab_md<-data.frame(model=rep('Additive',2),table=rbind(c(prob_AA_case, prob_AB_case_a_md, prob_BB_case_a_md),
                                                             c(prob_AA_control, P_AB-prob_AB_case_a_md,P_BB-prob_BB_case_a_md)))

  if(prob_BB_case_a_md<0 | prob_BB_control_a_md<0) {add_OR1_md<-add_OR2_md <- NaN
  }else {add_OR1_md <- add.tab_md[1,3]*(add.tab_md[2,2])/(add.tab_md[1,2])/add.tab_md[2,3]
  add_OR2_md <- add.tab_md[1,4]*(add.tab_md[2,2])/(add.tab_md[1,2])/add.tab_md[2,4]
  }
}else{add_OR1_md<-add_OR2_md<-NaN}
}}

##################################################################################
# Suppose Using Recessive Model.  What OR would be needed for a significant result
# under other models?
##################################################################################
if(test.model=='Recessive'){
  if(true.model=='Additive1'|true.model=='Additive2'){
#Additive
warning<-error<-0
tryCatch(zz<-uniroot(function(x) ((-(x^2)*prob_BB_case*P_AA + (x^2)*P_AA*P_BB)/(-2*x*prob_BB_case*P_AB + prob_BB_case*(P_AB^2) + (x^2)*P_BB)) -
                       ((-prob_BB_case*(P_AA + P_AB) + x*(prob_BB_case - prob_BB_case*o + o*P_BB))/(prob_BB_case*(-1 + o) - o*P_BB))
                     ,lower = 0, upper = 1)$root,
         warning=function(warn){warning<<-1},
         error=function(err){error<<-1})

if(warning+error==0){
  prob_AB_case_a_mr<-zz
  P_AB_case_a_mr<-prob_AB_case_a_mr/P_AB
  prob_AB_control_a_mr<-(1-P_AB_case_a_mr)*P_AB
  prob_AA_case_a_mr <- (prob_AB_case+prob_AA_case)-prob_AB_case_a_mr
  prob_AA_control_a_mr <- P_AA-prob_AA_case_a_mr

  add.tab_mr<-data.frame(model=rep('Additive',2),table=rbind(c(prob_AA_case_a_mr, prob_AB_case_a_mr, prob_BB_case),
                                                             c(prob_AA_control_a_mr, P_AB-prob_AB_case_a_mr,P_BB-prob_BB_case)))

  if(prob_AB_case_a_mr<0 | prob_AB_control_a_mr<0) {add_OR1_mr<-add_OR2_mr <- NaN
  }else {add_OR1_mr <- add.tab_mr[1,3]*(add.tab_mr[2,2])/(add.tab_mr[1,2])/add.tab_mr[2,3]
  add_OR2_mr <- add.tab_mr[1,4]*(add.tab_mr[2,2])/(add.tab_mr[1,2])/add.tab_mr[2,4]
  }
}else{add_OR1_mr<-add_OR2_mr<-NaN}
}

  if(true.model=='Dominant'){
#Dominant

P_AB_case_d_mr<-P_BB_case
prob_AB_case_d_mr<-(P_AB_case_d_mr*P_AB)
prob_AB_control_d_mr<-(1-P_AB_case_d_mr)*P_AB
prob_AA_case_d_mr <- (prob_AB_case+prob_AA_case)-prob_AB_case_d_mr
prob_AA_control_d_mr <- P_AA-prob_AA_case_d_mr

dom.tab_mr<-data.frame(model=rep('Dominant',2),table=rbind(c(prob_AA_case_d_mr, prob_AB_case_d_mr, prob_BB_case),
                                                           c(prob_AA_control_d_mr, P_AB-prob_AB_case_d_mr,P_BB-prob_BB_case)))

if(prob_AA_case_d_mr<0 | prob_AA_control_d_mr<0) {dom_OR_mr <- NaN
}else {dom_OR_mr <- dom.tab_mr[2,2]*(dom.tab_mr[1,3]+dom.tab_mr[1,4])/(dom.tab_mr[2,3]+dom.tab_mr[2,4])/dom.tab_mr[1,2]
}

}
##################################################################################
# Suppose Using Additive Model.  What OR would be needed for a significant result
# under other models?
##################################################################################
if(test.model=='Additive'){
  #First get betas and likelihood for additive2 model
  beta<-optim(c(0,0),function(x) -additive.ll(x,tab.assuming.test.true), control=c(abstol = 0.00001))$par
  like<-additive.ll(beta, tab.assuming.test.true)

    if(true.model=='Dominant'){
#Dominant

#Let x = prob of disease given AB or BB
#Then
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

root<-uniroot(function(x) find.prob.dom(x, m, cr, like=like),lower = cr, upper = 1)$root

dom.tab_ma2<-data.frame(model=rep('Dominant',2),table=rbind(c(Case.Rate-root*P_AB-root*P_BB, root*P_AB, root*P_BB),
                                                            c(P_AA-(Case.Rate-root*P_AB-root*P_BB), P_AB-root*P_AB,P_BB-root*P_BB)))


dom_OR_ma2 <- dom.tab_ma2[2,2]*(dom.tab_ma2[1,3]+dom.tab_ma2[1,4])/(dom.tab_ma2[2,3]+dom.tab_ma2[2,4])/dom.tab_ma2[1,2]
}
#Recessive
if(true.model=='Recessive'){
#Let x = prob of disease given BB
#Then
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

root<-uniroot(function(x) find.prob.rec(x, m, cr, like=like),lower = cr, upper = 1)$root

y<-(Case.Rate-root*P_BB)/(P_AA+P_AB)

rec.tab_ma2<-data.frame(model=rep('Recessive',2),table=rbind(c(P_AA*y, P_AB*y, root*P_BB),
                                                             c((1-P_AA)*y, (1-P_AB)*y,P_BB-root*P_BB)))



rec_OR_ma2 <- rec.tab_ma2[2,2]*(rec.tab_ma2[1,3]+rec.tab_ma2[1,4])/(rec.tab_ma2[2,3]+rec.tab_ma2[2,4])/rec.tab_ma2[1,2]
}}







  ########odds and ends
  #####################################################################
  # Function to create a true recessive model
  # x = prob case given BB
  # m = MAF
  # cr = case rate
  calc.rec.tab<-function(x, m, cr){
    c<-m*m*x #prob BB and case
    f<-m*m-c #prob BB and control
    y<-(cr-c)/(2*m*(1-m)+(1-m)*(1-m)) #conditional prob case in AA and AB
    b<-2*m*(1-m)*y #prob AB and case
    e<-2*m*(1-m)-b #prob AB and control
    a<-(1-m)*(1-m)*y #prob AA and case
    d<-(1-m)*(1-m)-a #prob AA and control
    tab <- rbind(c(a,b,c), c(d,e,f)) #table
    return(tab)
  }

  calc.dom.tab<-function(x, m, cr){
    a<-(1-m)*(1-m)*x #prob AA and case
    d<-(1-m)*(1-m)-a #prob AA and control
    y<-(cr-a)/(2*m*(1-m)+m*m) #conditional prob: probabilty of case give BB or AB
    c<-m*m*y #prob BB and case
    f<-m*m-c #prob BB and control
    b<-2*m*(1-m)*y #prob AB and case
    e<-2*m*(1-m)-b #prob AB and control
    tab <- rbind(c(a,b,c), c(d,e,f)) #table
    return(tab)
  }


  # Function to create a true additive model
  # x = prob case given AA
  # m = MAF
  # cr = case rate
  calc.add.tab<-function(x, m, cr){
    A<-(1-m)^2
    B<-2*m*(1-m)
    C<-m^2
    a<-A*x
    d<-A*(1-x)



    b<-uniroot(function(x) {((cr-a-x)*d /(a*(C-(cr-a-x)))) - (x*d/(a*(B-x)))^2}
               ,lower=(-(C-cr+a)), upper=B)$root
    e<-B-b
    c<-cr-b-a
    f<-C-c

    tab <- rbind(c(a,b,c), c(d,e,f)) #table
    return(tab)
  }
