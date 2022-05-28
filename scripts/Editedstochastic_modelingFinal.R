


#Packages to be used
packages<-c("here","tidyverse","ggplot2","gridExtra","lme4","lmtest","readxl")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))  



FG = read_excel("PrevComp.xlsx", sheet = "Sheet1", na = ".") 

View(FG)

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
prevp = c(0.01,seq(0.05,1,0.05)) #List of possible prevalences

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
ALPq =NULL
TLPq = NULL
ALPq2 = list()
TLPq2 = list()

#Poisson sampling function
poisrange <- function(n, lambda, smallest, highest){
  sample(x=smallest:highest, size=n, 
         prob=dpois(smallest:highest, lambda), replace=TRUE)
}

#Model proper
for (pp in 1:length(prevp)){ ##Iterations for each prevalence
  p = prevp[pp]
  
for(j in 1:5000){ #Each iteration is a different room by chance
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
  
    ##Random subsamplingsampling
    #for (ss in 1:1000){
      #hh = sample(ff, 15, replace = FALSE)
      
      #e2 = 1:length(hh)
      #g2 <- length(hh) + 1
      #t2 <- c(1, rep(0, g2-1))
      #for (i in hh) t2 <- (1-i)*t2 + i*c(0, t2[-g2])
      #Expectation calculation (Expected number of FOF per room)
      #t2<-as.matrix(t2,nrow = length(t2))
      #e2<-as.matrix(e2,ncol = length(hh))
      #a2 =  t2[-1,]%*%e2 #Expectation (Mean) of number of positive litters
      #ALP2[ss] = a2/length(hh) #Appare
    #}
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
  
  ALPq[pp] = data.frame(quantile(ALP, probs = c(0.5,0.025,0.975),names = FALSE))
  TLPq[pp] = data.frame(quantile(TLP, probs = c(0.5,0.025,0.975),names = FALSE))
} # End of the stochastic simulation node j


Output = data.frame(cbind(prevp,TLPq, ALPq))

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
  
# Plotting sensitivity

#plot(rt,TLPq2[[6]],type="l",lwd=3, col="red", ylim=c(0,max(unlist(TLPq2))), main="Sensitivity analysis of clustering", xlab="Proportion of PRRSV-positive pigs", ylab="True proportion of positive litters"
#)
#for   (tc in 1:length(clustern)){
# lines(rt,TLPq2[[tc]],type="l", col="grey")
#}
#rt = seq(0,0.5,0.05)

#2nd plot
all_data<-cbind.data.frame(cbind("0.00"=ALPq2[[1]],"0.126"=ALPq2[[2]],"0.252"=ALPq2[[3]],"0.378"=ALPq2[[4]],
                                 "0.504"=ALPq2[[5]],"0.630"=ALPq2[[6]],"0.704"=ALPq2[[7]],
                                 "0.778"=ALPq2[[8]],"0.852"=ALPq2[[9]],"0.926"=ALPq2[[10]], "1.000"=ALPq2[[11]]),prevp)

all_data<-all_data%>%gather(key="parameter","ALPq2",-prevp)

ggplot(all_data,aes(x=prevp,y=ALPq2,color=parameter))+
  geom_line()

#For TLP
all_data<-cbind.data.frame(cbind("0.00"=TLPq2[[1]],"0.126"=TLPq2[[2]],"0.252"=TLPq2[[3]],"0.378"=TLPq2[[4]],
                                 "0.504"=TLPq2[[5]],"0.630"=TLPq2[[6]],"0.704"=TLPq2[[7]],
                                 "0.778"=TLPq2[[8]],"0.852"=TLPq2[[9]],"0.926"=TLPq2[[10]], "1.000"=TLPq2[[11]]),prevp)

all_data<-all_data%>%gather(key="parameter","TLPq2",-prevp)

ggplot(all_data,aes(x=prevp,y=TLPq2,color=parameter))+
  geom_line()
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

#Specifying crate sizes
crates = c(10,20,28,40,56,102)
c = 0.63

for (tt in 1:length(crates)){
  nn = crates[tt]
  
  for (pp in 1:length(prevp)){ ##Iterations for each clustering
    p = prevp[pp]
    
    for(j in 1:1000){ #Each iteration is a different room by chance
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

# Plotting

plot(1:length(prevp),TLPq2[[3]],type="l",lwd=2,ylim=c(0,max(unlist(TLPq2))))
for   (tt in 1:length(crates)){
  lines(1:length(prevp),TLPq2[[tt]],type="l")
}
#or
plot(rt,ALPq2[[5]],type="l",lwd=3, col="red", ylim=c(0,0.7), main="Sensitivity analysis of number of litters within room", xlab="Proportion of PRRSV-positive pigs", ylab="Apparent proportion of positive litters"
)
for   (tt in 1:length(crates)){
  lines(rt,ALPq2[[tt]],type="l", col="grey")
}
###The values for predicted FOF status by WLP
ggh = data.frame(Prop_pos = seq(0,1,0.01))
fv = cbind(predict(MF2,  ggh, type = "response", re.form = NA))
uu = cbind(ggh,fv)

#Plot of the predictive model
library(esquisse)
esquisse::esquisser()
#OR
library(ggplot2)
ggplot(uu) +
 aes(x = Prop_pos, y = fv) +
 geom_line(size = 0.5, colour = "#145695") +
 labs(x = "Within-litter prevalence", 
 y = "Probability of PRRSV detection in FOF", title = "Predicted probability of PRRSV detection in FOF by Within-litter prevalence (WLP)") +
 theme_linedraw() +
 theme(plot.title = element_text(size = 20L, face = "bold", hjust = 0.5), axis.title.y = element_text(size = 18L, 
 face = "bold"), axis.title.x = element_text(size = 18L, face = "bold"))

###Cross validating the quality of the fourier transform###
rrt = data.frame(Prop_pos= c(0,0,0,0,0,0.556))
rrt2 = predict(MF2, rrt, type = "response", re.form = NA)
yu = NULL
for (rt in 1:6000){
  
  x = rbinom(6,1,rrt2)
  
  yu[rt] = sum(x)
}
summary(yu)
