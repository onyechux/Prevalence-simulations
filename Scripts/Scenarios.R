

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

#Model as function
scenarios <- function(nrpl){
  
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
  
  sce_list<-list(ALP=ALP_me,TLP=TLP_me) #Saving output
  
  capture.output(sce_list, file = here("Output","Scenarios","Scenarios.txt"), append=FALSE) #Saving crude outcome
  
  
} # End of the function
