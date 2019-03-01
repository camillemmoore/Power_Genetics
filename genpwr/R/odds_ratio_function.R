
####################################################################################
# For a dominant model
####################################################################################
#' Dominant Model Function
#'
#' Operates within odds_ratio_function to calculate odds ratios for a Test.Model of "Dominant"
#'
#' @param like Expected log likelihood
#' @param Case.Rate proportion of cases in the sample (cases/(cases + controls)). 
#' @param P_AA Probability the allele is homozygous for the major allele
#' @param P_AB Probability the allele is heteroygous
#' @param P_BB Probability the allele is homozygous for the minor allele
#' @param risk_allele Logical: If OR > 1, the allele is classified as a "risk allele"
#' @param True.Model A vector object specifying the true underlying genetic model(s): 'Dominant', 'Additive', or 'Recessive'
#'
#' @return: The odds ratios and their corresponding genetic model(s)
#'
#' @export
#'
dom.or.function <-function(like, Case.Rate, P_AA, P_AB, P_BB, True.Model, risk_allele)
{
  or.tab <- NULL
  cr <- Case.Rate

  a <- ll_zero_finder2(function(x) {x*log(x/P_AA)+
      (cr-x)*log((cr-x)/(P_AB+P_BB))+
      (P_AA-x)*log((P_AA-x)/P_AA)+
      (1-cr-P_AA+x)*log((1-cr-P_AA+x)/(P_AB+P_BB)) - like})#,

  d <- P_AA - a

  if ("Dominant" %in% True.Model){
    y <- (cr-a)/(P_AB+P_BB)
    b <- y*P_AB
    c <- y*P_BB
    e <- P_AB - b
    f <- P_BB - c
    mm00 <- !apply(rbind(a,b,c,d,e,f), 2, function(x) any(is.na(x)) | any(x > 1) | any(x < 0))
    a00 <- a[mm00]; b00 <- b[mm00]; c00 <- c[mm00]; d00 <- d[mm00]; e00 <- e[mm00]; f00 <- f[mm00];
    if(any(dim(rbind(a00,b00,c00,d00,e00,f00))==0)){
      or.tab <- rbind(or.tab, data.frame("Dominant", NA, fix.empty.names = F))
    }else if(all(!is.na(rbind(a00,b00,c00,d00,e00,f00)) & rbind(a00,b00,c00,d00,e00,f00) > 0 & rbind(a00,b00,c00,d00,e00,f00) < 1 & length(rbind(a00,b00,c00,d00,e00,f00)) > 0)){
      mmor <- or_calc(a00,b00,c00,d00,e00,f00, mod = "Dominant", risk_allele)
      if(length(mmor) == 0) mmor <- NA
      or.tab <- rbind(or.tab, data.frame("Dominant", mmor, fix.empty.names = F))
    }else{or.tab <- rbind(or.tab, data.frame("Dominant", NA, fix.empty.names = F))}
  }
  
  if ("Recessive" %in% True.Model){
    # The conditional probability of disease in AA and AB should be the same
    P_AB_case = a/P_AA

    # Now solve for the elements of the recessive table
    b = P_AB_case*P_AB
    c = cr - a - b
    e = P_AB - b
    f = P_BB - c
    mm00 <- !apply(rbind(a,b,c,d,e,f), 2, function(x) any(is.na(x)) | any(x > 1) | any(x < 0))
    a00 <- a[mm00]; b00 <- b[mm00]; c00 <- c[mm00]; d00 <- d[mm00]; e00 <- e[mm00]; f00 <- f[mm00];
    if(any(dim(rbind(a00,b00,c00,d00,e00,f00))==0)){
      or.tab <- rbind(or.tab, data.frame("Recessive", NA, fix.empty.names = F))
    }else if(all(!is.na(rbind(a00,b00,c00,d00,e00,f00)) & rbind(a00,b00,c00,d00,e00,f00) > 0 & rbind(a00,b00,c00,d00,e00,f00) < 1 & length(rbind(a00,b00,c00,d00,e00,f00)) > 0)){
      mmor <- or_calc(a00,b00,c00,d00,e00,f00, mod = "Recessive", risk_allele)
      if(length(mmor) == 0) mmor <- NA
      or.tab <- rbind(or.tab, data.frame("Recessive", mmor, fix.empty.names = F))
    }else{or.tab <- rbind(or.tab, data.frame("Recessive", NA, fix.empty.names = F))}
  }
  
  if ("Additive" %in% True.Model){
    #Here to solve for b, solve for c in terms of b in two ways, first using the dominant OR relationship and second setting OR(AB vs. AA) equal to OR(BB vs AB)
    #subtract the two version of c from one another and find the root of the equation
    b <- sapply(a, function(aa){
      bb <- ll_zero_finder2(function(x){
                    dd <- P_AA - aa
                    o = (cr-aa)*dd/(aa*(1-cr-dd))
                    (aa*(x-x*o+o*(P_BB+ P_AB)) - x*P_AA)/(aa*(-1 + o) + P_AA) - (P_BB*(P_AA-aa)*x^2)/(aa*P_AB*(P_AB-2*x)+P_AA*x^2)
                  })
      return(bb)
    })

    c = cr - a - b
    e = P_AB - b
    f = P_BB - c
    mm00 <- !apply(rbind(a,b,c,d,e,f), 2, function(x) any(is.na(x)) | any(x > 1) | any(x < 0))
    a00 <- a[mm00]; b00 <- b[mm00]; c00 <- c[mm00]; d00 <- d[mm00]; e00 <- e[mm00]; f00 <- f[mm00];
    if(any(dim(rbind(a00,b00,c00,d00,e00,f00))==0)){
      or.tab <- rbind(or.tab, data.frame("Additive", NA, fix.empty.names = F))
    }else if(all(!is.na(rbind(a00,b00,c00,d00,e00,f00)) & rbind(a00,b00,c00,d00,e00,f00) > 0 & rbind(a00,b00,c00,d00,e00,f00) < 1 & length(rbind(a00,b00,c00,d00,e00,f00)) > 0)){
      mmor <- or_calc(a00,b00,c00,d00,e00,f00, mod = "Additive", risk_allele)
      if(length(mmor) == 0) mmor <- NA
      or.tab <- rbind(or.tab, data.frame("Additive", mmor, fix.empty.names = F))
    }else{or.tab <- rbind(or.tab, data.frame("Additive", NA, fix.empty.names = F))}
  }
  
  return(or.tab)
}

####################################################################################
# For a recessive model
####################################################################################
#' Recessive Model Function
#'
#' Operates within odds_ratio_function to calculate odds ratios for a Test.Model of "Recessive"
#'
#' @param like Expected log likelihood
#' @param Case.Rate proportion of cases in the sample (cases/(cases + controls)). 
#' @param P_AA Probability the allele is homozygous for the major allele
#' @param P_AB Probability the allele is heteroygous
#' @param P_BB Probability the allele is homozygous for the minor allele
#' @param risk_allele Logical: If OR > 1, the allele is classified as a "risk allele"
#' @param True.Model A vector object specifying the true underlying genetic model(s): 'Dominant', 'Additive', or 'Recessive'
#'
#' @return: The odds ratios and their corresponding genetic model(s)
#'
#' @export
#'
rec.or.function <-function(like, Case.Rate, P_AA, P_AB, P_BB, True.Model, risk_allele)
{
  cr <- Case.Rate
  # Solve for c, the prob of BB and disease
  c <- ll_zero_finder2(function(x) {
        x*log(x/P_BB)+
        (cr-x)*log((cr-x)/(P_AB+P_AA))+
        (P_BB-x)*log((P_BB-x)/P_BB)+
        (1-cr-P_BB+x)*log((1-cr-P_BB+x)/(P_AB+P_AA)) - like})#,
  f <- P_BB - c

  or.tab <- NULL
  
  if('Recessive' %in% True.Model){
    y = (cr-c)/(P_AB+P_AA)
    a = y*P_AA
    d = P_AA - a
    b = y*P_AB
    e = P_AB - b
    mm00 <- !apply(rbind(a,b,c,d,e,f), 2, function(x) any(is.na(x)) | any(x > 1) | any(x < 0))
    a00 <- a[mm00]; b00 <- b[mm00]; c00 <- c[mm00]; d00 <- d[mm00]; e00 <- e[mm00]; f00 <- f[mm00];
    if(any(dim(rbind(a00,b00,c00,d00,e00,f00))==0)){
      or.tab <- rbind(or.tab, data.frame("Recessive", NA, fix.empty.names = F))
    }else if(all(!is.na(rbind(a00,b00,c00,d00,e00,f00)) & rbind(a00,b00,c00,d00,e00,f00) > 0 & rbind(a00,b00,c00,d00,e00,f00) < 1 & length(rbind(a00,b00,c00,d00,e00,f00)) > 0)){
      mmor <- or_calc(a00,b00,c00,d00,e00,f00, mod = "Recessive", risk_allele)
      if(length(mmor) == 0) mmor <- NA
      or.tab <- rbind(or.tab, data.frame("Recessive", mmor, fix.empty.names = F))
    }else{or.tab <- rbind(or.tab, data.frame("Recessive", NA, fix.empty.names = F))}
  }

  if('Dominant' %in% True.Model){
    # the conditional prob of disease for genotypes AB and BB are equal
    P_AB_case = c/P_BB
    P_AB_control = f/P_BB #

    # Solve for the joint probabilities
    b = P_AB_case*P_AB
    d = P_AB - b
    a = cr - b - c
    d = P_AA - a
    e = P_AB_control*P_AB
    mm00 <- !apply(rbind(a,b,c,d,e,f), 2, function(x) any(is.na(x)) | any(x > 1) | any(x < 0))
    a00 <- a[mm00]; b00 <- b[mm00]; c00 <- c[mm00]; d00 <- d[mm00]; e00 <- e[mm00]; f00 <- f[mm00];
    if(any(dim(rbind(a00,b00,c00,d00,e00,f00))==0)){
      or.tab <- rbind(or.tab, data.frame("Dominant", NA, fix.empty.names = F))
    }else if(all(!is.na(rbind(a00,b00,c00,d00,e00,f00)) & rbind(a00,b00,c00,d00,e00,f00) > 0 & rbind(a00,b00,c00,d00,e00,f00) < 1 & length(rbind(a00,b00,c00,d00,e00,f00)) > 0)){
      mmor <- or_calc(a00,b00,c00,d00,e00,f00, mod = "Dominant", risk_allele)
      if(length(mmor) == 0) mmor <- NA
      or.tab <- rbind(or.tab, data.frame("Dominant", mmor, fix.empty.names = F))
    }else{or.tab <- rbind(or.tab, data.frame("Dominant", NA, fix.empty.names = F))}
  }

  if('Additive' %in% True.Model){
  
    # To find b, solve for a in two ways (see dominant.or.function)
    b <- sapply(c, function(cc){
      if(is.na(cc)){
        return(NA)
      }else{
        bb <- ll_zero_finder2(function(x){
                      ff <- P_BB - cc
                      o <- cc*(1-cr-ff)/((cr-cc)*ff)
                      (-(x^2)*cc*P_AA + (x^2)*P_AA*P_BB)/(-2*x*cc*P_AB + cc*(P_AB^2) + (x^2)*P_BB) -
                         (-cc*(P_AA + P_AB) + x*(cc - cc*o + o*P_BB))/(cc*(-1 + o) - o*P_BB)
                    })
        return(bb)
      }
    })

    e = P_AB - b
    a = cr - b - c
    d = P_AA - a
    mm00 <- !apply(rbind(a,b,c,d,e,f), 2, function(x) any(is.na(x)) | any(x > 1) | any(x < 0))
    a00 <- a[mm00]; b00 <- b[mm00]; c00 <- c[mm00]; d00 <- d[mm00]; e00 <- e[mm00]; f00 <- f[mm00];
    if(any(dim(rbind(a00,b00,c00,d00,e00,f00))==0)){
      or.tab <- rbind(or.tab, data.frame("Additive", NA, fix.empty.names = F))
    }else if(all(!is.na(rbind(a00,b00,c00,d00,e00,f00)) & rbind(a00,b00,c00,d00,e00,f00) > 0 & rbind(a00,b00,c00,d00,e00,f00) < 1 & length(rbind(a00,b00,c00,d00,e00,f00)) > 0)){
      mmor <- or_calc(a00,b00,c00,d00,e00,f00, mod = "Additive", risk_allele)
      if(length(mmor) == 0) mmor <- NA
      or.tab <- rbind(or.tab, data.frame("Additive", mmor, fix.empty.names = F))
    }else{or.tab <- rbind(or.tab, data.frame("Additive", NA, fix.empty.names = F))}
  }
  return(or.tab)
}

####################################################################################
# For an additive model
####################################################################################

#' Binomial coefficient calculation
#'
#' Operates within add.or.function to solve for 'a' when 'b' is known in an additive model
#'
#' @param b The "b" in the binomial function ax^2 + bx + c that arises in solution for the additive OR functions
#' @param cr proportion of cases in the sample (cases/(cases + controls)). 
#' @param P_AA Probability the allele is homozygous for the major allele
#' @param P_AB Probability the allele is heteroygous
#' @param P_BB Probability the allele is homozygous for the minor allele
#'
#' @return: The "a" in the binomial function ax^2 + bx + c that arises in solution for the additive OR functions
#'
#' @export
#'
solve_a<-function(b, cr, P_AA, P_AB, P_BB)
{
  a1<- P_AB*(2*b-P_AB)
  b1<- cr*(P_AB^2 - 2*P_AB*b)-b*P_AB^2 + (2*P_AB-P_AA+P_BB)*b^2
  c1<- -(b^2)*(P_BB-cr+b)*P_AA

  a.test<-quad_roots(a1,b1,c1)[1]
  a<-Re(polyroot(c(c1,b1,a1)))
  a <- a[a>0]
  a <- a[a<1]
  if(length(a) == 0) a <- NA
  return(a)
}

#' Dominant probability finding function
#'
#' Operates within add.or.function to find probability of disease in a dominant truth given AB or BB, additive test model
#'
#' @param x Probability of disease given AB or BB
#' @param P_AA Probability the allele is homozygous for the major allele
#' @param P_AB Probability the allele is heteroygous
#' @param P_BB Probability the allele is homozygous for the minor allele
#' @param cr proportion of cases in the sample (cases/(cases + controls)). 
#' @param like Expected log likelihood
#'
#' @return: The "a" in the binomial function ax^2 + bx + c that arises in solution for the additive OR functions
#'
#' @export
#'
find.prob.dom<-function(x,P_AA,P_AB,P_BB,cr, like)
{
  b<-P_AB*x
  c<-P_BB*x
  e<-P_AB-b
  f<-P_BB-c
  a<-cr-b-c
  d<-P_AA-a

  tab <- rbind(c(a,b,c), c(d,e,f))
  mles<-optim(c(-1,0),function(x) -additive.ll(x,tab), control=c(abstol = 0.00001))$par
  # ll_calc <- sum(tab[1,]*log(tab[1,]/c(P_AA,P_AB,P_BB))) + sum(tab[2,]*log(tab[2,]/c(P_AA,P_AB,P_BB)))
  additive.ll(mles, tab)-like
  # ll_calc - like
}

#' Recessive probability finding function
#'
#' Operates within add.or.function to find probability of disease in a recessive truth given AB or BB, additive test model
#'
#' @param x Probability of disease given AB or BB
#' @param P_AA Probability the allele is homozygous for the major allele
#' @param P_AB Probability the allele is heteroygous
#' @param P_BB Probability the allele is homozygous for the minor allele
#' @param cr proportion of cases in the sample (cases/(cases + controls)). 
#' @param like Expected log likelihood
#'
#' @return: The "a" in the binomial function ax^2 + bx + c that arises in solution for the additive OR functions
#'
#' @export
#'
find.prob.rec<-function(x,P_AA,P_AB,P_BB,cr, like)
{
  c<-P_BB*x
  f<-P_BB-c
  y<-(cr-c)/(P_AB+P_AA)
  b<-P_AB*y
  e<-P_AB-b
  a<-P_AA*y
  d<-P_AA-a

  tab <- rbind(c(a,b,c), c(d,e,f))
  mles<-optim(c(-1,0),function(x) -additive.ll(x,tab), control=c(abstol = 0.00001))$par
  # ll_calc <- sum(tab[1,]*log(tab[1,]/c(P_AA,P_AB,P_BB))) + sum(tab[2,]*log(tab[2,]/c(P_AA,P_AB,P_BB)))
  additive.ll(mles, tab)-like
  # ll_calc - like
}


#' Additive Model Function
#'
#' Operates within odds_ratio_function to calculate odds ratios for a Test.Model of "Additive"
#'
#' @param like Expected log likelihood
#' @param Case.Rate proportion of cases in the sample (cases/(cases + controls)). 
#' @param P_AA Probability the allele is homozygous for the major allele
#' @param P_AB Probability the allele is heteroygous
#' @param P_BB Probability the allele is homozygous for the minor allele
#' @param risk_allele Logical: If OR > 1, the allele is classified as a "risk allele"
#' @param True.Model A vector object specifying the true underlying genetic model(s): 'Dominant', 'Additive', or 'Recessive'
#'
#' @return: The odds ratios and their corresponding genetic model(s)
#'
#' @export
#'
add.or.function <-function(like, Case.Rate, P_AA, P_AB, P_BB, True.Model, risk_allele)
{

  cr <- Case.Rate
  or.tab <- NULL
  # For a true additive model
  if ("Additive" %in% True.Model){
    b <- ll_zero_finder2(function(x){
      a <- solve_a(x, cr, P_AA, P_AB, P_BB)#[1]
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
    })
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
        a[[ix]] <- a[[ix]][which.min(abs(
          a[[ix]]*log(a[[ix]]/P_AA)+b[ix]*log(b[ix]/P_AB)+c*log(c/P_BB)+d*log(d/P_AA)+e*log(e/P_AB)+f*log(f/P_BB) - like))]
      }
      a <- unlist(a)
      if(length(a) != length(b)) stop("tried to fix multiple zeros in solve_a, but there is still a problem")
    }else a <- unlist(a)
    d<-P_AA-a
    e<-P_AB-b
    c<-cr-a-b
    f<-P_BB-c
    mm00 <- !apply(rbind(a,b,c,d,e,f), 2, function(x) any(is.na(x)) | any(x > 1) | any(x < 0))
    a00 <- a[mm00]; b00 <- b[mm00]; c00 <- c[mm00]; d00 <- d[mm00]; e00 <- e[mm00]; f00 <- f[mm00];
    if(any(dim(rbind(a00,b00,c00,d00,e00,f00))==0)){
      or.tab <- rbind(or.tab, data.frame("Additive", NA, fix.empty.names = F))
    }else if(all(!is.na(rbind(a00,b00,c00,d00,e00,f00)) & rbind(a00,b00,c00,d00,e00,f00) > 0 & rbind(a00,b00,c00,d00,e00,f00) < 1 & length(rbind(a00,b00,c00,d00,e00,f00)) > 0)){
      mmor <- or_calc(a00,b00,c00,d00,e00,f00, mod = "Additive", risk_allele)
      if(length(mmor) == 0) mmor <- NA
      or.tab <- rbind(or.tab, data.frame("Additive", mmor, fix.empty.names = F))
    }else{or.tab <- rbind(or.tab, data.frame("Additive", NA, fix.empty.names = F))}
  }

  #Dominant
  if ("Dominant" %in% True.Model){
    #x = prob of disease given AB or BB
    x <- ll_zero_finder2(function(x) find.prob.dom(x, P_AA,P_AB,P_BB, cr, like=like))

    a = cr-x*P_AB-x*P_BB
    b = x*P_AB
    c = x*P_BB
    d = P_AA - a
    e = P_AB - b
    f = P_BB - c
    mm00 <- !apply(rbind(a,b,c,d,e,f), 2, function(x) any(is.na(x)) | any(x > 1) | any(x < 0))
    a00 <- a[mm00]; b00 <- b[mm00]; c00 <- c[mm00]; d00 <- d[mm00]; e00 <- e[mm00]; f00 <- f[mm00];
    if(any(dim(rbind(a00,b00,c00,d00,e00,f00))==0)){
      or.tab <- rbind(or.tab, data.frame("Dominant", NA, fix.empty.names = F))
    }else if(all(!is.na(rbind(a00,b00,c00,d00,e00,f00)) & rbind(a00,b00,c00,d00,e00,f00) > 0 & rbind(a00,b00,c00,d00,e00,f00) < 1 & length(rbind(a00,b00,c00,d00,e00,f00)) > 0)){
      mmor <- or_calc(a00,b00,c00,d00,e00,f00, mod = "Dominant", risk_allele)
      if(length(mmor) == 0) mmor <- NA
      or.tab <- rbind(or.tab, data.frame("Dominant", mmor, fix.empty.names = F))
    }else{or.tab <- rbind(or.tab, data.frame("Dominant", NA, fix.empty.names = F))}
  }
  #Recessive
  if('Recessive' %in% True.Model){
    #Let x = prob of disease given BB
    x <- ll_zero_finder2(function(x) find.prob.rec(x, P_AA,P_AB,P_BB, cr, like=like))
    
    y<-(cr-x*P_BB)/(P_AA+P_AB)
    a = P_AA*y
    b = P_AB*y
    c = x*P_BB
    d = P_AA - a
    e = P_AB - b
    f = P_BB - c
    mm00 <- !apply(rbind(a,b,c,d,e,f), 2, function(x) any(is.na(x)) | any(x > 1) | any(x < 0))
    a00 <- a[mm00]; b00 <- b[mm00]; c00 <- c[mm00]; d00 <- d[mm00]; e00 <- e[mm00]; f00 <- f[mm00];
    if(any(dim(rbind(a00,b00,c00,d00,e00,f00))==0)){
      or.tab <- rbind(or.tab, data.frame("Recessive", NA, fix.empty.names = F))
    }else if(all(!is.na(rbind(a00,b00,c00,d00,e00,f00)) & rbind(a00,b00,c00,d00,e00,f00) > 0 & rbind(a00,b00,c00,d00,e00,f00) < 1 & length(rbind(a00,b00,c00,d00,e00,f00)) > 0)){
      mmor <- or_calc(a00,b00,c00,d00,e00,f00, mod = "Recessive", risk_allele)
      if(length(mmor) == 0) mmor <- NA
      or.tab <- rbind(or.tab, data.frame("Recessive", mmor, fix.empty.names = F))
    }else{or.tab <- rbind(or.tab, data.frame("Recessive", NA, fix.empty.names = F))}
  }
  return(or.tab)
}

####################################################################################
# For a 2DF model
####################################################################################
#' 2df Model Function
#'
#' Operates within odds_ratio_function to calculate odds ratios for a Test.Model of "2df"
#'
#' @param like Expected log likelihood
#' @param Case.Rate proportion of cases in the sample (cases/(cases + controls)). 
#' @param P_AA Probability the allele is homozygous for the major allele
#' @param P_AB Probability the allele is heteroygous
#' @param P_BB Probability the allele is homozygous for the minor allele
#' @param risk_allele Logical: If OR > 1, the allele is classified as a "risk allele"
#' @param True.Model A vector object specifying the true underlying genetic model(s): 'Dominant', 'Additive', or 'Recessive'
#'
#' @return: The odds ratios and their corresponding genetic model(s)
#'
#' @export
#'
or.function.2df <-function(like, Case.Rate, P_AA, P_AB, P_BB, True.Model, risk_allele)
{
  cr <- Case.Rate
  or.tab <- NULL
  # ll.alt <- ll.null+(stat_2df/2)
  # like <- ll.alt


  # If model is dominant, will estimate the same coef for AB and BB
  # Solve for a, the joint prob of geno AA and disease
  if ("Dominant" %in% True.Model){
    a <- ll_zero_finder2(function(x) {x*log(x/P_AA)+
        (cr-x)*log((cr-x)/(P_AB+P_BB))+
        (P_AA-x)*log((P_AA-x)/P_AA)+
        (1-cr-P_AA+x)*log((1-cr-P_AA+x)/(P_AB+P_BB)) - like})

    d <- P_AA - a
    y <- (cr-a)/(P_AB+P_BB)
    b <- y*P_AB
    c <- y*P_BB
    e <- P_AB - b
    f <- P_BB - c
    mm00 <- !apply(rbind(a,b,c,d,e,f), 2, function(x) any(is.na(x)) | any(x > 1) | any(x < 0))
    a00 <- a[mm00]; b00 <- b[mm00]; c00 <- c[mm00]; d00 <- d[mm00]; e00 <- e[mm00]; f00 <- f[mm00];
    if(any(dim(rbind(a00,b00,c00,d00,e00,f00))==0)){
      or.tab <- rbind(or.tab, data.frame("Dominant", NA, fix.empty.names = F))
    }else if(all(!is.na(rbind(a00,b00,c00,d00,e00,f00)) & rbind(a00,b00,c00,d00,e00,f00) > 0 & rbind(a00,b00,c00,d00,e00,f00) < 1 & length(rbind(a00,b00,c00,d00,e00,f00)) > 0)){
      mmor <- or_calc(a00,b00,c00,d00,e00,f00, mod = "Dominant", risk_allele)
      if(length(mmor) == 0) mmor <- NA
      or.tab <- rbind(or.tab, data.frame("Dominant", mmor, fix.empty.names = F))
    }else{or.tab <- rbind(or.tab, data.frame("Dominant", NA, fix.empty.names = F))}
  }

  # If model is recessive
  if ("Recessive" %in% True.Model){
    c <- ll_zero_finder2(function(x){
        x*log(x/P_BB)+
        (cr-x)*log((cr-x)/(P_AB+P_AA))+
        (P_BB-x)*log((P_BB-x)/P_BB)+
        (1-cr-P_BB+x)*log((1-cr-P_BB+x)/(P_AB+P_AA)) - like})#,

    f <- P_BB - c
    y = (cr-c)/(P_AB+P_AA)
    a = y*P_AA
    d = P_AA - a
    b = y*P_AB
    e = P_AB - b
    mm00 <- !apply(rbind(a,b,c,d,e,f), 2, function(x) any(is.na(x)) | any(x > 1) | any(x < 0))
    a00 <- a[mm00]; b00 <- b[mm00]; c00 <- c[mm00]; d00 <- d[mm00]; e00 <- e[mm00]; f00 <- f[mm00];
    if(any(dim(rbind(a00,b00,c00,d00,e00,f00))==0)){
      or.tab <- rbind(or.tab, data.frame("Recessive", NA, fix.empty.names = F))
    }else if(all(!is.na(rbind(a00,b00,c00,d00,e00,f00)) & rbind(a00,b00,c00,d00,e00,f00) > 0 & rbind(a00,b00,c00,d00,e00,f00) < 1 & length(rbind(a00,b00,c00,d00,e00,f00)) > 0)){
      mmor <- or_calc(a00,b00,c00,d00,e00,f00, mod = "Recessive", risk_allele)
      if(length(mmor) == 0) mmor <- NA
      or.tab <- rbind(or.tab, data.frame("Recessive", mmor, fix.empty.names = F))
    }else{or.tab <- rbind(or.tab, data.frame("Recessive", NA, fix.empty.names = F))}
  }

  # If model is additive
  if("Additive" %in% True.Model){
    solve_a_2df_add<-function(b){
      a1<- P_AB*(2*b-P_AB)
      b1<- cr*(P_AB^2 - 2*P_AB*b)-b*P_AB^2 + (2*P_AB-P_AA+P_BB)*b^2
      c1<- -(b^2)*(P_BB-cr+b)*P_AA

      a<-Re(polyroot(c(c1,b1,a1)))#[1]

      d<-P_AA-a
      e<-P_AB-b
      c<-cr-a-b
      f<-P_BB-c

      # getting rid of negative values
      # zznum <- unique(unlist(sapply(c("a", "b", "c", "d", "e", "f"), function(aletter) eval(parse(text=sprintf("which(%s < 0)", aletter))))))
      # if(zznum >= 2) return(NA)
      # if(zznum == 1) for(aletter in c("a", "b", "c", "d", "e", "f")) eval(parse(text=sprintf("%s <- %s[-zznum]", aletter, aletter)))
      logna <- function(x) sapply(x, function(xx) ifelse(xx < 0, return(NA), return(log(xx))))

      ll <- a*logna(a/P_AA)+
        b*logna(b/P_AB)+
        c*logna(c/P_BB)+
        d*logna(d/P_AA)+
        e*logna(e/P_AB)+
        f*logna(f/P_BB)
      ll <- ll[!is.na(ll)]
        return(like - ll)
    }
    # b<-uniroot(solve_a, lower=P_AB*cr+0.0000000001, upper=min(P_AB,cr-P_BB)-0.000001 )$root #this gives the root for an OR>1
    b<-ll_zero_finder2(solve_a_2df_add) #this gives the root for an OR>1

    a1<- P_AB*(2*b-P_AB)
    b1<- cr*(P_AB^2 - 2*P_AB*b)-b*P_AB^2 + (2*P_AB-P_AA+P_BB)*b^2
    c1<- -(b^2)*(P_BB-cr+b)*P_AA
    a <- list()
    for(ii1 in 1:length(a1)) {
      if(any(sapply(c(c1[ii1], b1[ii1], a1[ii1]), is.na))){
        a <- c(a, list(NA))
      }else{
        a <- c(a, list(Re(polyroot(c(c1[ii1], b1[ii1], a1[ii1])))))
      }
    }
    if(any(sapply(a, length) > 1)){
      # what to do if there are multiple zeros in polyroot:
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
        a[[ix]] <- a[[ix]][which.min(abs(
          a[[ix]]*log(a[[ix]]/P_AA)+b[ix]*log(b[ix]/P_AB)+c*log(c/P_BB)+d*log(d/P_AA)+e*log(e/P_AB)+f*log(f/P_BB) - like))]
      }
      a <- unlist(a)
      if(length(a) != length(b)) stop("tried to fix multiple zeros in solve_a, but there is still a problem")
    }else a <- unlist(a)
    d<-P_AA-a
    e<-P_AB-b
    c<-cr-a-b
    f<-P_BB-c
    mm00 <- !apply(rbind(a,b,c,d,e,f), 2, function(x) any(is.na(x)) | any(x > 1) | any(x < 0))
    a00 <- a[mm00]; b00 <- b[mm00]; c00 <- c[mm00]; d00 <- d[mm00]; e00 <- e[mm00]; f00 <- f[mm00];
    if(any(dim(rbind(a00,b00,c00,d00,e00,f00))==0)){
      or.tab <- rbind(or.tab, data.frame("Additive", NA, fix.empty.names = F))
    }else if(all(!is.na(rbind(a00,b00,c00,d00,e00,f00)) & rbind(a00,b00,c00,d00,e00,f00) > 0 & rbind(a00,b00,c00,d00,e00,f00) < 1 & length(rbind(a00,b00,c00,d00,e00,f00)) > 0)){
      mmor <- or_calc(a00,b00,c00,d00,e00,f00, mod = "Additive", risk_allele)
      if(length(mmor) == 0) mmor <- NA
      or.tab <- rbind(or.tab, data.frame("Additive", mmor, fix.empty.names = F))
    }else{or.tab <- rbind(or.tab, data.frame("Additive", NA, fix.empty.names = F))}
  }
  return(or.tab)
}

#' Odds Ratio Function
#'
#' Calculates the odds ratio for a given power, at a given sample size, N, with type 1 error rate, Alpha
#'
#' @param N Vector of the desired sample size(s)
#' @param Case.Rate Vector of the proportion(s) of cases in the sample (cases/(cases + controls)). Either k or Case.Rate must be specified.
#' @param k Vector of the number of controls per case. Either k or Case.Rate must be specified.
#' @param Alpha the desired type 1 error rate(s)
#' @param MAF Vector of minor allele frequencies
#' @param power Vector of powers to detect
#' @param risk_allele Logical: If OR > 1, the allele is classified as a "risk allele"
#' @param True.Model A vector vector the true underlying genetic model(s): 'Dominant', 'Additive', 'Recessive' or 'All'
#' @param Test.Model A vector specifying the assumed genetic model(s) used in testing: 'Dominant', 'Additive', 'Recessive' or 'All'
#'
#' @examples
#' or <- odds_ratio_function(N=c(100,1000), Case.Rate=seq(from=0.2,to=0.3,by=0.1), k=NULL, MAF= seq(from=0.2, to = 0.3, by = 0.025), 
#'   power=c(0.5, 0.8), Alpha = c(0.05), risk_allele = T, True.Model = 'All', Test.Model = 'All')
#'
#' @export
#'
odds_ratio_function <-
  function(N=NULL, Case.Rate=NULL, k=NULL, MAF=NULL, power=NULL, risk_allele = T,
                     Alpha=0.05, True.Model='All', Test.Model='All')
{

  ############################################################################################################
  #Error Messages for insufficient sample size information, MAF, and case vs. control ratio
  ############################################################################################################
  if(is.null(N)) {
    stop("N, the total sample size, must be specified.")
  }

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

  if(any(Test.Model == "Additive2")){
    Test.Model[Test.Model == "Additive2"] <- "Additive"
  }

  if(any(Test.Model == "Additive1")){
    Test.Model[Test.Model == "Additive1"] <- "Additive"
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

  ############################################################################################################
  #Calculate needed sample size information from provided inputs
  ############################################################################################################
  #Test model vector
  if('All' %in% Test.Model){Test.Model<-c("Dominant", "Recessive", "Additive", "2df")}

  #True model vector
  if('All' %in% True.Model){True.Model<-c("Dominant", "Recessive", "Additive")}

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
    #Use power's and MAF's to calculate true distribution of genotypes and disease
    ############################################################################################################
    ##############################################################################
    #For each power
    ##############################################################################
    pow.save.tab <-NULL

    for (pow in power){

      alph.save.tab <- NULL
      
      for(alpha0 in Alpha){

        ############################################################################################################
        #log likelihood stats
        ############################################################################################################

        #stat = ((qnorm(1-alpha0/2)+qnorm(pow))^2)/N
        stat = uniroot(function(x) ncp.search(x, pow, stat, alpha0, df=1),
                           lower=0, upper=1000, extendInt = 'upX', tol=0.00001)$root/N
        
        #For the 2DF Test the detectable LRT test statistic is:
        stat_2df = uniroot(function(x) ncp.search(x, pow, stat, alpha0, df=2),
                       lower=0, upper=1000, extendInt = 'upX', tol=0.00001)$root/N

        beta0 <- logit(Case.Rate)
        ll.null <- Case.Rate*log(exp(beta0)/(1+exp(beta0))) + (1-Case.Rate)*log(1/(1+exp(beta0)))

        # Note that stat = 2*(as.numeric(ll.alt-ll.null))
        # Then the alt.likelihood to detect is
        ll.alt <- ll.null+(stat/2)
        like <- ll.alt

        ll.alt_2df <- ll.null+(stat_2df/2)
        like_2df <- ll.alt_2df

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

          ############################################################################################################
          #Calculate Odds Ratio for each scenario under the specified testing model
          ############################################################################################################

          if('Dominant' %in% Test.Model){
            s.tab <- dom.or.function(like=like, Case.Rate=Case.Rate, P_AA=P_AA, P_AB=P_AB, P_BB=P_BB, True.Model=True.Model, risk_allele=risk_allele)
            s.tab <- data.frame(Test.Model = rep('Dominant', nrow(s.tab)), True.Model = s.tab[,1], MAF=rep(m, nrow(s.tab)), power = rep(pow, nrow(s.tab)),
              N_total = rep(N, nrow(s.tab)), N_cases = rep(N_cases, nrow(s.tab)), N_controls = rep(N_controls, nrow(s.tab)), 
                                       Case.Rate = rep(Case.Rate, nrow(s.tab)), OR = s.tab[,2], Alpha = rep(alpha0, nrow(s.tab)))
            save.tab <- rbind(save.tab, s.tab)
          }

          if('Additive' %in% Test.Model){
            s.tab <- add.or.function(like=like, Case.Rate=Case.Rate, P_AA=P_AA, P_AB=P_AB, P_BB=P_BB, True.Model=True.Model, risk_allele=risk_allele)
            s.tab <- data.frame(Test.Model = rep('Additive', nrow(s.tab)), True.Model = s.tab[,1], MAF=rep(m, nrow(s.tab)), power = rep(pow, nrow(s.tab)),
              N_total = rep(N, nrow(s.tab)), N_cases = rep(N_cases, nrow(s.tab)), N_controls = rep(N_controls, nrow(s.tab)), 
                                       Case.Rate = rep(Case.Rate, nrow(s.tab)), OR = s.tab[,2], Alpha = rep(alpha0, nrow(s.tab)))
            save.tab <- rbind(save.tab, s.tab)
          }

          if('Recessive' %in% Test.Model){
            s.tab <- rec.or.function(like=like, Case.Rate=Case.Rate, P_AA=P_AA, P_AB=P_AB, P_BB=P_BB, True.Model=True.Model, risk_allele=risk_allele)
            s.tab <- data.frame(Test.Model = rep('Recessive', nrow(s.tab)), True.Model = s.tab[,1], MAF=rep(m, nrow(s.tab)), power = rep(pow, nrow(s.tab)),
              N_total = rep(N, nrow(s.tab)), N_cases = rep(N_cases, nrow(s.tab)), N_controls = rep(N_controls, nrow(s.tab)), 
                                       Case.Rate = rep(Case.Rate, nrow(s.tab)), OR = s.tab[,2], Alpha = rep(alpha0, nrow(s.tab)))
            save.tab <- rbind(save.tab, s.tab)
          }

          if('2df' %in% Test.Model){
            s.tab <- or.function.2df(like=like_2df, Case.Rate=Case.Rate, P_AA=P_AA, P_AB=P_AB, P_BB=P_BB, True.Model=True.Model, risk_allele=risk_allele)
            s.tab <- data.frame(Test.Model = rep('2df', nrow(s.tab)), True.Model = s.tab[,1], MAF=rep(m, nrow(s.tab)), power = rep(pow, nrow(s.tab)),
              N_total = rep(N, nrow(s.tab)), N_cases = rep(N_cases, nrow(s.tab)), N_controls = rep(N_controls, nrow(s.tab)), 
                                       Case.Rate = rep(Case.Rate, nrow(s.tab)), OR = s.tab[,2], Alpha = rep(alpha0, nrow(s.tab)))
            save.tab <- rbind(save.tab, s.tab)
          }

          m.save.tab<-rbind(m.save.tab, save.tab)
            # Test.Model True.Model MAF OR N_total N_cases N_controls Case.Rate Power_at_Alpha_0.05
        }
        alph.save.tab <- rbind(alph.save.tab, m.save.tab)
      }
      pow.save.tab<-rbind(pow.save.tab, alph.save.tab)
    }

    final.or.tab<-rbind(final.or.tab, pow.save.tab)
  }

  # fix table so that it's in the same format as Camille's
  final.or.tab2 <- final.or.tab[final.or.tab$Alpha == unique(final.or.tab$Alpha)[1],][,1:8]
  for(analpha in unique(final.or.tab$Alpha)){
    final.or.tab2 <- cbind(final.or.tab2, NA)
    names(final.or.tab2)[ncol(final.or.tab2)] <- paste0("OR_at_Alpha_", analpha)
    final.or.tab2[,paste0("OR_at_Alpha_", analpha)] <- final.or.tab[final.or.tab$Alpha == analpha, "OR"]
  }
  return(final.or.tab2)
}
