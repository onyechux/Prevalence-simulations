
FG = read_excel(here("Data","PrevComp.xlsx"),sheet = "Sheet1",na =".")



FG$Farm = as.factor(FG$Farm)
FG$ROOM = as.factor(FG$ROOM)
FG$FOFC=as.factor(FG$FOFC)
FG$ROOM=as.factor(FG$ROOM)
#Categorical variable of FOF positive (1) or negative (0)


MF2 = glmer(FOFC~Prop_pos + (1|Farm), family = binomial(link = "logit"), data = FG, na.action = na.omit)



#############################
# Recursive binomial model  #
#############################

#Set parameters
prevalence = c(0.01,seq(0.05,1,0.05)) #List of possible prevalences

n = 56
c = 0.63

# Creating objects
prob<-NULL
N<-NULL
TLP<-NULL
teste<-NULL
ff = NULL
prev<-NULL
ALP<-NULL
ALP2= NULL
ALPmedian =NULL
TLPmedian = NULL
ALPlowerquantile = NULL
TLPlowerquantile = NULL
ALPupperquantile = NULL
TLPupperquantile = NULL

#Truncated Poisson sampling function
poisrange <- function(n, lambda, smallest, highest){
  sample(x=smallest:highest, size=n, 
         prob=dpois(smallest:highest, lambda), replace=TRUE)
}

#Model proper
test <- function(nrpl){
  
for (pp in 1:length(prevalence)){ ##Iterations for each prevalence
  p = prevalence[pp]
  
  for(j in 1:nrpl){ #Each iteration is a different room by chance
    cat("\n\n Room: ", j, "\n")
    tempo_inicial = proc.time()
    
    # Room characteristics
    
    tot<-poisrange(n,11,2,17) # Litters size (stochastic in j)
    
    tot_t<-sum(tot) # Number of piglets in the room
    
    pos<-round(p*tot_t,0) #Expected number of positive piglets in the room
    
    
    #Initial values for the binomial recursive model
    
    prob[1]<-(tot[1]/tot_t)+(1-(tot[1]/tot_t))*c #Initial condition for the within litter prevalence 
    
    N[1]<-min(rbinom(1,pos,prob[1]),tot[1]) #Initial condition for the number of positive piglets in the litter
    
    #lit_prev[1] = N[1]/tot[1]
    
    for(i in 2:n){ # Start the binomial recursive model
      
      prob[i]<-tot[i]/(tot_t-sum(tot[1:i-1]))+(1-(tot[i]/(tot_t-sum(tot[1:i-1]))))*c # Vector of within litter prevalence
      
      N[i]<-min(rbinom(1,pos-sum(N[1:i-1]),prob[i]),tot[i]) #Vector for the within litter number of positive piglets
      ###
      
    } #End of the recursive Binomial model
    
    
    prev<-data.frame(Prop_pos=N/tot)
    ff = predict(MF2, prev, type = "response", re.form = NA)
    
    #Convolution 
    e = 1:n
    g <- length(ff) + 1
    t <- c(1, rep(0, g-1))
    for (i in ff) t <- (1-i)*t + i*c(0, t[-g])
    
    #Expectation calculation (Expected number of FOF per room)
    t<-as.matrix(t,nrow = length(t))
    e<-as.matrix(e,ncol = n)
    a =  t[-1,]%*%e #Expectation (Mean) of number of positive litters
    ALP[j] = a/n #Apparent litter level prevalence in the iteration j
    
    
    TLP[j]<-mean(N>0) #Proportion of positive litters in the iteration j
  }
  
  ALPmedian[pp] = round(median(ALP),3)
  TLPmedian[pp] = round(median(TLP),3)
  ALPlowerquantile[pp] = round(quantile(ALP, probs = 0.025,names = FALSE),3)
  TLPlowerquantile[pp] = round(quantile(TLP, probs = 0.025,names = FALSE), 3)
  ALPupperquantile[pp] = round(quantile(ALP, probs = 0.975,names = FALSE),3)
  TLPupperquantile[pp] = round(quantile(TLP, probs = 0.975,names = FALSE),3)
  #TLPq[pp] = data.frame(quantile(TLP, probs = c(0.5,0.025,0.975),names = FALSE))
  #TLPq[pp] = data.frame(quantile(TLP, probs = c(0.5,0.025,0.975),names = FALSE))
  
} # End of the stochastic simulation node j

Output = data.frame(cbind(prevalence,TLPmedian, ALPmedian, TLPlowerquantile, ALPlowerquantile, 
                          TLPupperquantile, ALPupperquantile))
Tab = datatable(Output)

capture.output(Tab)
}

test(nrpl=5000)
