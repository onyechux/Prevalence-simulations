##Sensitivity for room size##
# For looping

# Creating objects
prob<-NULL
N<-NULL
TLP<-NULL
teste<-NULL
ff = NULL
prev<-NULL
ALP<-NULL
ALPq =NULL
TLPq = NULL
ALPq2 = list()
TLPq2 = list()

prevalence = c(0.01,seq(0.1,0.5,0.1)) #List of possible prevalences

crates = c(10,20,28,40,56,102) #Specifying crate sizes

c = 0.63 #Fixed clusteringg

#Truncated Poisson sampling function
poisrange <- function(n, lambda, smallest, highest){
  sample(x=smallest:highest, size=n,
         prob=dpois(smallest:highest, lambda), replace=TRUE)
}

sensitivity <- function(nrp){
  
for (tt in 1:length(crates)){
  nn = crates[tt]
  
  for (pp in 1:length(prevp)){ ##Iterations for each clustering
    p = prevp[pp]
    
    for(j in 1:nrp){ #Each iteration is a different room by chance
      cat("\n\n Room: ", j, "\n")
      tempo_inicial = proc.time()
      
      # Room characteristics
      tot<-poisrange(nn,11,2,19) # Litters size (stochastic in j)
      
      tot_t<-sum(tot) # Number of piglets in the room
      
      pos<-round(p*tot_t,0) #Expected number of positive piglets in the room
      
      
      #Initial values for the binomial recursive model
      
      prob[1]<-(tot[1]/tot_t)+(1-(tot[1]/tot_t))*c #Initial condition for the within litter prevalence 
      
      N[1]<-min(rbinom(1,pos,prob[1]),tot[1]) #Initial condition for the number of positive piglets in the litter
      
      #lit_prev[1] = N[1]/tot[1]
      
      for(i in 2:nn){ # Start the binomial recursive model
        
        prob[i]<-tot[i]/(tot_t-sum(tot[1:i-1]))+(1-(tot[i]/(tot_t-sum(tot[1:i-1]))))*c # Vector of within litter prevalence
        
        N[i]<-min(rbinom(1,pos-sum(N[1:i-1]),prob[i]),tot[i]) #Vector for the within litter number of positive piglets
        ###
        
      } #End of the recursive Binomial model
      
      
      prev<-data.frame(Prop_pos=N/tot)
      ff = predict(MF2, prev, type = "response", re.form = NA)
      
      
      #Convolution 
      e = 1:nn
      g <- length(ff) + 1
      t <- c(1, rep(0, g-1))
      for (i in ff) t <- (1-i)*t + i*c(0, t[-g])
      
      
      #Expectation calculation (Expected number of FOF per room)
      t<-as.matrix(t,nrow = length(t))
      e<-as.matrix(e,ncol = nn)
      a =  t[-1,]%*%e #Expectation (Mean) of number of positive litters
      ALP[j] = a/nn #Apparent litter level prevalence in the iteration j
      
      
      TLP[j]<-mean(N>0) #Proportion of positive litters in the iteration j
    }
    
    ALPq[pp] = median(ALP)
    TLPq[pp] = median(TLP)
  } # End of the stochastic simulation node j
  ALPq2[tt] = data.frame(ALPq)
  TLPq2[tt] = data.frame(TLPq)
}

  Output = data.frame (TLPq2, ALPq2)
  
  # Plotting
  
  #For ALP##
  all_data<-cbind.data.frame(cbind("10"=ALPq2[[1]],"20"=ALPq2[[2]],"28"=ALPq2[[3]],"40"=ALPq2[[4]],
                                   "56"=ALPq2[[5]],"102"=ALPq2[[6]]),prevp)
  
  all_data<-all_data%>%gather(key="parameter","ALPq2",-prevp)
  
  ggplot(all_data,aes(x=prevp,y=ALPq2,color=parameter))+
    geom_line()
  
  ##For TLP##
  all_data<-cbind.data.frame(cbind("10"=TLPq2[[1]],"20"=TLPq2[[2]],"28"=TLPq2[[3]],"40"=TLPq2[[4]],
                                   "56"=TLPq2[[5]],"102"=TLPq2[[6]]),prevp)
  
  all_data<-all_data%>%gather(key="parameter","TLPq2",-prevp)
  
  ggplot(all_data,aes(x=prevp,y=TLPq2,color=parameter))+
    geom_line()
  
}

sensitivity(nrp = 1000)
 