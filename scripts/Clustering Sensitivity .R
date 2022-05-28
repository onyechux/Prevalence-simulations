###Sensitivity analysis###


#Loading data
FG = read_excel("PrevComp.xlsx", sheet = "Sheet1", na = ".")

FG$Farm = as.factor(FG$Farm)
FG$ROOM = as.factor(FG$ROOM)
FG$FOFC=as.factor(FG$FOFC)
FG$ROOM=as.factor(FG$ROOM)
#Categorical variable of FOF positive (1) or negative (0)


MF2 = glmer(FOFC~Prop_pos + (1|Farm), family = binomial(link = "logit"), data = FG, na.action = na.omit)

#Set parameters
prevp = c(0.01,seq(0.05,1,0.05)) #List of possible prevalences
n = 56

#Truncated Poisson sampling function
poisrange <- function(n, lambda, smallest, highest){
  sample(x=smallest:highest, size=n,
         prob=dpois(smallest:highest, lambda), replace=TRUE)
}

##Sensitivity for clustering##

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

#Specifying clustering levels
#

clusx<-seq(-1,1,0.2) #Support vector for percentage change around the baseline

base <- 0.63
base_min<-0 #min for clustering
base_max<-1 #max for clustering

clustern<-ifelse(clusx<0,base+(base-base_min)*clusx,base-(base-base_max)*clusx) #Calculation for the uncertainty analysis

#Alternative
#clustern = c(0.2,0.3,0.4,0.55,0.8,0.9)  #0.55 is my baseline
n = 56

for (tc in 1:length(clustern)){
  c = clustern[tc]

  for (pp in 1:length(prevp)){ ##Iterations for each clustering
    p = prevp[pp]

    for(j in 1:5000){ #Each iteration is a different room by chance
      cat("\n\n Room: ", j, "\n")
      tempo_inicial = proc.time()

      # Room characteristics
      tot<-poisrange(n,11,2,19) # Litters size (stochastic in j)

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

    ALPq[pp] = median(ALP)
    TLPq[pp] = median(TLP)
  } # End of the stochastic simulation node j
  ALPq2[tc] = data.frame(ALPq)
  TLPq2[tc] = data.frame(TLPq)
}

####End of model runs#####


###Plots###

#For ALP##
all_data<-cbind.data.frame(cbind("0.00"=ALPq2[[1]],"0.126"=ALPq2[[2]],"0.252"=ALPq2[[3]],"0.378"=ALPq2[[4]],
                                 "0.504"=ALPq2[[5]],"0.630"=ALPq2[[6]],"0.704"=ALPq2[[7]],
                                 "0.778"=ALPq2[[8]],"0.852"=ALPq2[[9]],"0.926"=ALPq2[[10]], "1.000"=ALPq2[[11]]),prevp)

all_data<-all_data%>%gather(key="parameter","ALPq2",-prevp)

ggplot(all_data,aes(x=prevp,y=ALPq2,color=parameter))+
  geom_line()

##For TLP##
all_data<-cbind.data.frame(cbind("0.00"=TLPq2[[1]],"0.126"=TLPq2[[2]],"0.252"=TLPq2[[3]],"0.378"=TLPq2[[4]],
                                 "0.504"=TLPq2[[5]],"0.630"=TLPq2[[6]],"0.704"=TLPq2[[7]],
                                 "0.778"=TLPq2[[8]],"0.852"=TLPq2[[9]],"0.926"=TLPq2[[10]], "1.000"=TLPq2[[11]]),prevp)

all_data<-all_data%>%gather(key="parameter","TLPq2",-prevp)

ggplot(all_data,aes(x=prevp,y=TLPq2,color=parameter))+
  geom_line()
