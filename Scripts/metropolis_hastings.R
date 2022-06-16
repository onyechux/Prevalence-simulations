


#read data
FG = read_excel(here("Data","PrevComp.xlsx"),sheet = "Sheet1",na =".")

#select rooms
Rm1 = subset(FG, FG$ROOM==9)  
Rm2 = subset(FG, FG$ROOM==10)
Rm3 = subset(FG, FG$ROOM==11)
Rm4 = subset(FG, FG$Farm=='EAST' & FG$ROOM==12)
Rm6 = subset(FG, FG$ROOM==3)
Rm8 = subset(FG, FG$Farm=='G&W LEIBOL' & FG$ROOM==6)
Rm10 = subset(FG, FG$Farm=='PETE SCHNEIDER' & FG$ROOM==6)

#sort rooms

Rm1<-Rm1[order(-Rm1$Prop_pos),]
Rm2<-Rm2[order(-Rm2$Prop_pos),] 
Rm3<-Rm3[order(-Rm3$Prop_pos),] 
Rm4<-Rm4[order(-Rm4$Prop_pos),] 
Rm6<-Rm6[order(-Rm6$Prop_pos),] 
Rm8<-Rm8[order(-Rm8$Prop_pos),] 
Rm10<-Rm10[order(-Rm10$Prop_pos),] 




x<-list(Rm1$Positive,Rm2$Positive,Rm3$Positive,Rm4$Positive,Rm6$Positive,Rm8$Positive,Rm10$Positive)
tot<-list(Rm1$`TOTAL PIGS`,Rm2$`TOTAL PIGS`,Rm3$`TOTAL PIGS`,Rm4$`TOTAL PIGS`,Rm6$`TOTAL PIGS`,Rm8$`TOTAL PIGS`,Rm10$`TOTAL PIGS`)

tot_tc<-NULL
posc<-NULL
pc<-NULL

for (i in 1:length(x)){
  tot_tc[i]<-sum(tot[[i]])
  pc[i]<-sum(x[[i]])/tot_tc[i]
  posc[i]<-round(pc[i]*tot_tc[i],0)
  
}




prob<-vector(mode = "list", length = length(x))
N<-vector(mode = "list", length = length(x))
z_prob<-vector(mode = "list", length = length(x))
like<-vector(mode = "list", length = length(x))
sd_d2<-vector(mode = "list", length = length(x))
like2<-NULL
sd_d<-NULL

#Define likelihood and prior, then multiply them to get the posterior

#likelihood

likelihood <- function(ck){
  
  for(i in 1:length(x)){
    
    prob[[i]][1]<-(
      (tot[[i]][1]/(tot_tc[i])
      +(1-(tot[[i]][1]/tot_tc[i]))*ck)
    )
    
    N[[i]][1]<-round(posc[i]*prob[[i]][1])
   
    if (!is.infinite(dbinom(x=x[[i]][1],size=tot[[i]][1],prob=prob[[i]][1],log=T))       
        )  
    {
      like[[i]][1]<-NA}    else{
        
        like[[i]][1]<-dbinom(x=x[[i]][1],size=tot[[i]][1],prob=prob[[i]][1],log=T) 
      } 
   
    for(j in 2:length(x[[i]])){ #Initial binomial recursive model
    
      prob[[i]][j]<-(
        (tot[[i]][j]/(tot_tc[i]-sum(tot[[i]][1:j-1])))
             +(1-(tot[[i]][j]/(tot_tc[i]-sum(tot[[i]][1:j-1]))))*ck
      )
  
      N[[i]][j]<-min(round((posc[i]-sum(N[[i]][1:j-1]))*prob[[i]][j]),tot_tc[i])
      
      if (!is.infinite(dbinom(x=x[[i]][j],size=tot[[i]][j],prob=prob[[i]][j],log=T)) 
          )  
      {
        like[[i]][j]<-NA}    else{
          
          like[[i]][j]<-dbinom(x=x[[i]][j],size=tot[[i]][j],prob=prob[[i]][j],log=T) 
        }
      
      
    } #End binomial recursive
  
    like2[i]<-sum(like[[i]],na.rm = T)
  }
  likelihood<-sum(like2)
  
  return(likelihood)   
}

# Prior distribution (beta)
prior <- function(ck){
  prior <- sum(dbeta(ck, shape1=1, shape2=1, log = TRUE))
  return (prior)
}

# Posterior is the sum of logs, or multiplication in natural scale
posterior <- function(ck){
  return (likelihood(ck) + prior(ck))
}

#Metropolis-Hasting algorithm

#initialize parameters
burnIn = 100
startvalue = c(0.5)
iterations <- 10000
chain = array(dim = c(iterations+1, 1)) #1 parameter
chain[1] = startvalue

pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                     max = iterations, # Maximum value of the progress bar
                     style = 3,    # Progress bar style (also available style = 1 and style = 2)
                     width = 50,   # Progress bar width. Defaults to getOption("width")
                     char = "=")   # Character used to create the bar

for (i in 1:iterations){

  #sampling from a proposal function one value, based on the current value
  proposal <- rnorm(1,mean = chain[i],sd=1)
  
  #Metropolis-hastings ratio (subtraction of logs)
  acceptance = exp(posterior(proposal) - #posterior value at the new parameter value sampled
                     posterior(chain[i])) #posterior value at the current parameter value sampled
  
  if (runif(1) < acceptance){ # Logic function for "acceptance"
    chain[i+1] <- proposal #Move the chain to the next value
  }else{
    chain[i+1] = chain[i] #Retain the old value
  }
  
  setTxtProgressBar(pb, i)  
  
}
close(pb)

#remove burn in
chain <- chain[-(1:burnIn)]

quantile(chain,0.025);mean(chain);quantile(chain,0.975)
summary(chain)
hist(chain)
plot(chain,type = "l")
