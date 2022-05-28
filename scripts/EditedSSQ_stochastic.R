#####Learning clustering####

#Packages to be used
packages<-c("here","tidyverse","ggplot2","gridExtra","lme4","lmtest","readxl")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))  


##Calling each room from the data data##
FG = read_excel(here("Data","PrevComp.xlsx"),sheet = "Sheet1",na =".")

Rm1 = subset(FG, FG$ROOM==9)  
Rm2 = subset(FG, FG$ROOM==10)
Rm3 = subset(FG, FG$ROOM==11)
Rm4 = subset(FG, FG$Farm=='EAST' & FG$ROOM==12)
Rm5 = subset(FG, FG$Farm=='EGP' & FG$ROOM==12) #All negatives
Rm6 = subset(FG, FG$ROOM==3)
Rm7 = subset(FG, FG$ROOM==4) #All negatives
Rm8 = subset(FG, FG$Farm=='G&W LEIBOL' & FG$ROOM==6)
Rm9 = subset(FG, FG$Farm=='LONGRUN' & FG$ROOM==6) #All negatives
Rm10 = subset(FG, FG$Farm=='PETE SCHNEIDER' & FG$ROOM==6)
Rm11 = subset(FG, FG$ROOM==1)#All negatives



##Monte carlo simulation to investigate clustering level##

arg_min=NULL

for (k in 1:5000){
  cat("\n\n node: ", k, "\n")
  tempo_inicial = proc.time()
  
  # Creating objects
  prob<-NULL
  N<-NULL
  pos_prop<-NULL
  prev<-NULL
  Results = NULL
  mssq = NULL
  cck = NULL
  res_fit=NULL
  Rm<-NULL
  nc=NULL
  
#Sampling the room
    
  
    Rm<-get(paste("Rm",sample(c(1:4,6,8,10),1),sep=""))
    Rm <- Rm[order(-Rm$Prop_pos),] #Rearranging rooms to have the litters with highest
    
    #Specifying prevalence and number 
    pc = sum(Rm$Positive)/sum(Rm$`TOTAL PIGS`) #Prevalence of PRRSV at the piglet level #You can repeat by changing the numbers
    nc = length(Rm$`TOTAL PIGS`) #number of crates sampled
    sortrm = Rm$Prop_pos 
    
    
    # Room characteristics
    totc<- Rm$`TOTAL PIGS` # Litters size (stochastic in j)
    tot_tc<-sum(totc) # Number of piglets in the room
    posc<-round(pc*tot_tc,0) #Expected number of positive piglets in the room
    
    
#number of positive pigs in front


  for(j in 1:300){ #Each iteration is a different room by chance
    
    
    
  #Initial values for the binomial recursive model
  cck[j]<-runif(1,0,1)
  prob[1]<-(totc[1]/tot_tc)+(1-(totc[1]/tot_tc))*cck[j] #Initial condition for the within litter prevalence 
  
  N[1]<-min(rbinom(1,posc,prob[1]),totc[1]) #Initial condition for the number of positive piglets in the litter
  
  #lit_prev[1] = N[1]/tot[1]
  
    for(i in 2:nc){ # Start the binomial recursive model
    
    prob[i]<-totc[i]/(tot_tc-sum(totc[1:i-1]))+(1-(totc[i]/(tot_tc-sum(totc[1:i-1]))))*cck[j] # Vector of within litter prevalence
    
    N[i]<-min(rbinom(1,posc-sum(N[1:i-1]),prob[i]),totc[i]) #Vector for the within litter number of positive piglets
    ###
    
  } #End of the recursive Binomial model
  prevc<-data.frame(yy=N/totc) #the proportion of positive pigs per litter in each room
  
  kc = sum((sortrm - prevc)^2)/nc  #Sum of squares difference for matched litters per room
  
  Results[j] = kc 
  
  
} # End of the stochastic simulation node j


#Visualizing results
#plot(Results~cck) #the clustering with the lowest mssq is chosen

final_data<-cbind.data.frame(Results,cck)

arg_min[k]<-final_data[which.min(final_data$Results),][2]

}

#Visualize and getting estimates
df=unlist(arg_min)
hist(df)
summary(df)
