####################################################################
#Search for Non-Centrality Parameter for SS calcualtion for 2Df Test
####################################################################

ncp.search<-function(x, power, stat, Alpha, df){
  (1-power)-pchisq(qchisq(1-Alpha, df, ncp=0), df=1, ncp = x) }
