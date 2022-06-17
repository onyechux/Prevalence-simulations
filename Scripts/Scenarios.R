

# Reading and tidy the data
FG = read_excel(here("Data","PrevComp.xlsx"),sheet = "Sheet1",na =".")

FG$Farm = as.factor(FG$Farm)
FG$ROOM = as.factor(FG$ROOM)
FG$FOFC=as.factor(FG$FOFC)
FG$ROOM=as.factor(FG$ROOM)
#Categorical variable of FOF positive (1) or negative (0)

# Logistic regression

MF2 = glmer(FOFC~Prop_pos + (1|Farm), family = binomial(link = "logit"), data = FG, na.action = na.omit)


#################################
# Stochastic monte carlo model  #
#################################

#Read source

source(here("Data", "data.R"))


scenarios <- function(nrpl,nbase,sce_n,cbase,sce_c){ #Model as function
  
  sce_n<-c(sce_n,rep(nbase,length(sce_n)))
  sce_c<-c(rep(cbase,length(sce_c)),sce_c)
  
  
  
  ALP_me<-vector(mode = "list", length = length(sce_n))
  TLP_me<-vector(mode = "list", length = length(sce_c))
  
  
  
for (s in 1:length(sce_n)){ #Scenarios
  cat("\n\n Scenario: ", s,"n=",sce_n[s]," c=",sce_c[s],"\n")
  tempo_inicial = proc.time()
  
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = length(prevalence), # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar
  
  for (w in 1:length(prevalence)){ ##Iterations for each prevalence
    
    
    
    
    
    for(j in 1:nrpl){ #Each iteration is a different room by chance
      
      
      # Creating objects for i (litter)  
      prob<-NULL
      N<-NULL
      
      
      # Room characteristics
      
      tot<-sample(sort(unique(FG$`TOTAL PIGS`)),sce_n[s],replace=T,prob=prop.table(table(FG$`TOTAL PIGS`))) # Litters' size (stochastic in j) based on Almeida's data
      
      tot_t<-sum(tot) # Number of piglets in the room
      
      pos<-round(prevalence[w]*tot_t,0) #Expected number of positive piglets in the room
      
      
      #Initial values for the binomial recursive model
      
      prob[1]<-(tot[1]/tot_t)+(1-(tot[1]/tot_t))*sce_c[s] #Initial condition for the within litter prevalence
      
      N[1]<-min(rbinom(1,pos,prob[1]),tot[1]) #Initial condition for the number of positive piglets in the litter
      
      
      
      
      for(i in 2:sce_n[s]){ # Start the binomial recursive model
        
        prob[i]<-tot[i]/(tot_t-sum(tot[1:i-1]))+(1-(tot[i]/(tot_t-sum(tot[1:i-1]))))*sce_c[s] # Vector of within litter prevalence
        
        N[i]<-min(rbinom(1,pos-sum(N[1:i-1]),prob[i]),tot[i]) #Vector for the within litter number of positive piglets
        ###
        
      } #End of the recursive Binomial model
      
      
      prev<-data.frame(Prop_pos=N/tot)
      ff = predict(MF2, prev, type = "response", re.form = NA)
      
      pos_prop<-rbinom(length(ff),1,ff) # status of each crate in the iteration j
      
      ALP[j]<-mean(pos_prop) # Proportion of observed positive crates in the iteration j
      
      TLP[j]<-mean(N>0) #Proportion of positive litters in the iteration j
      
      
  
      
    } # End of the stochastic simulation node j
    
    ALP_me[[s]][w]<-mean(ALP)
    TLP_me[[s]][w]<-mean(TLP)
    
    
    
    
    setTxtProgressBar(pb, w)  
  
  } #End of w piglets prevalence
  close(pb)
  
    

  } # End of S Scenarios
  dir.create(here("Output","Scenarios"),showWarnings = F) #Create the directory for the output
  dir.create(here("Figures","Scenarios"),showWarnings = F) #Create the directory for the figure
  
  sce_list<-list(ALP=ALP_me,TLP=TLP_me) #Saving output
  
  
  capture.output(sce_list, file = here("Output","Scenarios","Scenarios.txt"), append=FALSE) #Saving crude outcome
  
  all_data<-cbind.data.frame(cbind("10"=ALP_me[[1]],"30"=ALP_me[[2]],"50"=ALP_me[[3]],"90"=ALP_me[[4]], prevalence))#ALP vs Crates
  all_data2<-cbind.data.frame(cbind("0.01"=ALP_me[[5]],"0.3"=ALP_me[[6]],"0.6"=ALP_me[[7]], "0.9"=ALP_me[[8]],prevalence)) #ALP vs Clustering
  all_data3<-cbind.data.frame(cbind("10"=TLP_me[[1]],"30"=TLP_me[[2]],"50"=TLP_me[[3]],"90"=TLP_me[[4]], prevalence)) #TLP vs crates
  all_data4<-cbind.data.frame(cbind("0.01"=TLP_me[[5]],"0.3"=TLP_me[[6]],"0.6"=TLP_me[[7]], "0.9"=TLP_me[[8]],prevalence)) #TLP vs Clustering
  
  ALP_crate<-all_data%>%gather(key="Number_of_crates",value="ALP",-prevalence)#vs crates
  TLP_crate<-all_data3%>%gather(key="Number_of_crates",value="TLP",-prevalence)#vs crates
  ALP_clust<-all_data2%>%gather(key="Clustering_level",value="ALP",-prevalence)#vs Clustering
  TLP_clust<-all_data4%>%gather(key="Clustering_level",value="TLP",-prevalence)#vs Clustering
  
  G1<-ggplot(ALP_crate,aes(x=prevalence,y=ALP,color=Number_of_crates))+
    geom_line()
  G2<-ggplot(TLP_crate,aes(x=prevalence,y=TLP,color=Number_of_crates))+
    geom_line()
  G3<-ggplot(ALP_clust,aes(x=prevalence,y=ALP,color=Clustering_level))+
    geom_line()
  G4<-ggplot(TLP_clust,aes(x=prevalence,y=TLP,color=Clustering_level))+
    geom_line()
    
  ggarrange(G1, G2, G3, G4, widths = c(5.5,5,5))
  ggsave(here("Figures","Scenarios","Figure2.png"),width = 7,height = 5,device = "png",dpi=300)  #save the plot
  
} # End of the function
