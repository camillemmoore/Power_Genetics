#########################################################################################################
#Plotting functions to summarize
#########################################################################################################

#########################################################################################################
#Plots for Sample Size
#########################################################################################################


ss.plot<-function(data=NULL,x='MAF', panel.by='True.Model',
                  select.Alpha = NULL, 
                  select.OR = NULL,
                  select.Power = NULL, 
                  select.MAF = NULL, 
                  select.Case.Rate = NULL,
                  select.True.Model = NULL,
                  select.Test.Model = NULL){
  
  #Create Dataset with separate rows for each Alpha level
  indices<-grep("N_total_at_Alpha_",colnames(data))
  Alphas<-as.numeric(matrix(unlist( strsplit(colnames(data)[indices],"N_total_at_Alpha_")), ncol=2, byrow=TRUE)[,2])
  
  ss.new<-NULL
  for (i in 1:length(Alphas)){ss.new<-rbind(ss.new, data.frame(data[-indices], N_total=data[,indices[i]], Alpha=Alphas[i]))}
  
  #Subset data to only include those in the select. statements
  if(is.null(select.Alpha)==F){ss.new<-ss.new[ss.new$Alpha %in% select.Alpha,]}
  if(is.null(select.OR)==F){ss.new<-ss.new[ss.new$OR %in% select.OR,]}
  if(is.null(select.Power)==F){ss.new<-ss.new[ss.new$Power %in% select.Power,]}
  if(is.null(select.MAF)==F){ss.new<-ss.new[ss.new$MAF %in% select.MAF,]}
  if(is.null(select.Case.Rate)==F){ss.new<-ss.new[ss.new$Case.Rate %in% select.Case.Rate,]}
  if(is.null(select.True.Model)==F){ss.new<-ss.new[ss.new$True.Model %in% select.True.Model,]}
  if(is.null(select.Test.Model)==F){ss.new<-ss.new[ss.new$Test.Model %in% select.Test.Model,]}
  
  
  var<-c("MAF", 'OR','Power', 'Case.Rate', 'Alpha', 'True.Model')
  var<-var[!(var %in% c(x, panel.by))]
  
  graphs<-unique.data.frame(ss.new[,var])
  
  for(j in 1:nrow(graphs)){subtitle<-paste("N by ", x, ": ",
                                           var[1], '=', graphs[j,1], ', ',
                                           var[2], '=', graphs[j,2],', ',
                                           var[3], '=', graphs[j,3],', ',
                                           var[4], '=', graphs[j,4], sep='')
  temp1<-ss.new[ss.new[,var[1]]==graphs[j,1] & ss.new[,var[2]]==graphs[j,2] 
                &ss.new[,var[3]]==graphs[j,3] & ss.new[,var[4]]==graphs[j,4],]
  
  temp2<-temp1[order(temp1$True.Model),]
  temp2[,panel.by]<-paste(panel.by,"=", temp2[,panel.by])
  
  print(ggplot(data=temp2, aes(x=temp2[,x], y=N_total, group = Test.Model, colour = Test.Model)) +
          geom_line() +
          geom_point( size=4, shape=21, fill="white") + xlab(x)+ggtitle(subtitle)+ 
          facet_wrap( ~temp2[,panel.by])
  )
  }}


#########################################################################################################
#Plots for Power
#########################################################################################################

power.plot<-function(data=NULL,x='MAF', panel.by='True.Model',
                     select.Alpha = NULL, 
                     select.OR = NULL,
                     select.N = NULL, 
                     select.MAF = NULL, 
                     select.Case.Rate = NULL,
                     select.True.Model = NULL,
                     select.Test.Model = NULL){
  
  #Create Dataset with separate rows for each Alpha level
  indices<-grep("Power_at_Alpha_",colnames(data))
  Alphas<-as.numeric(matrix(unlist( strsplit(colnames(data)[indices],"Power_at_Alpha_")), ncol=2, byrow=TRUE)[,2])
  
  ss.new<-NULL
  for (i in 1:length(Alphas)){ss.new<-rbind(ss.new, data.frame(data[-indices], Power=data[,indices[i]], Alpha=Alphas[i]))}
  
  #Subset data to only include those in the select. statements
  if(is.null(select.Alpha)==F){ss.new<-ss.new[ss.new$Alpha %in% select.Alpha,]}
  if(is.null(select.OR)==F){ss.new<-ss.new[ss.new$OR %in% select.OR,]}
  if(is.null(select.N)==F){ss.new<-ss.new[ss.new$N_total %in% select.N,]}
  if(is.null(select.MAF)==F){ss.new<-ss.new[ss.new$MAF %in% select.MAF,]}
  if(is.null(select.Case.Rate)==F){ss.new<-ss.new[ss.new$Case.Rate %in% select.Case.Rate,]}
  if(is.null(select.True.Model)==F){ss.new<-ss.new[ss.new$True.Model %in% select.True.Model,]}
  if(is.null(select.Test.Model)==F){ss.new<-ss.new[ss.new$Test.Model %in% select.Test.Model,]}
  
  
  var<-c("MAF", 'OR','N_total', 'Case.Rate', 'Alpha', 'True.Model')
  var<-var[!(var %in% c(x, panel.by))]
  
  graphs<-unique.data.frame(ss.new[,var])
  
  for(j in 1:nrow(graphs)){subtitle<-paste("Power by ", x, ": ",
                                           var[1], '=', graphs[j,1], ', ',
                                           var[2], '=', graphs[j,2],', ',
                                           var[3], '=', graphs[j,3],', ',
                                           var[4], '=', graphs[j,4], sep='')
  temp1<-ss.new[ss.new[,var[1]]==graphs[j,1] & ss.new[,var[2]]==graphs[j,2] 
                &ss.new[,var[3]]==graphs[j,3] & ss.new[,var[4]]==graphs[j,4],]
  
  temp2<-temp1[order(temp1$True.Model),]
  temp2[,panel.by]<-paste(panel.by,"=", temp2[,panel.by])
  
  print(ggplot(data=temp2, aes(x=temp2[,x], y=Power, group = Test.Model, colour = Test.Model)) +
          geom_line() +
          geom_point( size=4, shape=21, fill="white") + xlab(x)+ggtitle(subtitle)+ 
          facet_wrap( ~temp2[,panel.by])
  )
  }}
