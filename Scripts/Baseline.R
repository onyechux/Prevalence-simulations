
# Creating folders
dir.create(here("Output","Baseline"),showWarnings = F) #Create the directory for the output
dir.create(here("Figures","Baseline"),showWarnings = F) #Create the directory for the figure

# Reading and tidy the data
FG = read_excel(here("Data","PrevComp.xlsx"),sheet = "Sheet1",na =".")

FG$Farm = as.factor(FG$Farm)
FG$ROOM = as.factor(FG$ROOM)
FG$FOFC2 = FG$FOFC
FG$FOFC=as.factor(FG$FOFC2)
FG$ROOM=as.factor(FG$ROOM)
#Categorical variable of FOF positive (1) or negative (0)

# Logistic regression

MF2 = glmer(FOFC~Prop_pos + (1|Farm), family = binomial(link = "logit"), data = FG, na.action = na.omit)
pred<-predictInterval(MF2, newdata = MF2@frame, n.sims = 999,type = c("probability"),which = "fixed")
  
 # predict(MF2,  MF2@frame, type = "response", re.form = NA, alpha = 0.1)

pred.final<-cbind.data.frame(MF2@frame,fitted=pred)
pred.final$FOFC<-as.numeric(ifelse(pred.final$FOFC==1,1,0))


capture.output(pred.final, file = here("Output","Baseline","Logistic_prediction.txt"), append=FALSE) #Saving crude outcome


#################################################################
#Creating plot to compare the observed data and predictive model#
#################################################################

ggplot(pred.final, aes(x=Prop_pos,y=FOFC))+
  theme_minimal()+
  geom_jitter(position = position_jitter(width = 0.01,height = 0.06),alpha=0.5,size=1.5,shape=21)+
  geom_line(aes(x=Prop_pos,y=fitted.fit),cex=1.5) + 
  geom_ribbon(aes(ymin = fitted.lwr, ymax = fitted.upr), alpha = 0.1)+
  labs(x = "Within-litter prevalence", 
       y = "Probability of PRRSV detection in FOF")+ 
    theme(axis.title.y = element_text(size = 14L, face = "bold"), axis.title.x = element_text(size = 12L, face = "bold"))+
  scale_x_continuous(labels = function(Prop_pos) paste0(Prop_pos*100, "%"),breaks = seq(from = 0, to = 1, by = 0.25))


ggsave(here("Figures","Baseline","Figure1.png"),width = 7,height = 5,device = "png",dpi=300,bg="white")  #save the plot



#################################
# Stochastic monte carlo model  #
#################################

#Read source

source(here("Data", "data.R"))

#Model as function
baseline <- function(nrpl,n,c){

for (w in 1:length(prevalence)){ ##Iterations for each prevalence
  cat("\n\n piglet_prev= ",paste0(prevalence[w]*100,"%",sep=""),"\n")
  tempo_inicial = proc.time()

  
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = nrpl, # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar
  
  
  for(j in 1:nrpl){ #Each iteration is a different room by chance
    
      
    # Creating objects for i (litter)  
    prob<-NULL
    N<-NULL
    
    
    # Room characteristics

    tot<-sample(sort(unique(FG$`TOTAL PIGS`)),n,replace=T,prob=prop.table(table(FG$`TOTAL PIGS`))) # Litters' size (stochastic in j) based on Almeida's data
    
    tot_t<-sum(tot) # Number of piglets in the room
    
    pos<-round(prevalence[w]*tot_t,0) #Expected number of positive piglets in the room


    #Initial values for the binomial recursive model

    prob[1]<-(tot[1]/tot_t)+(1-(tot[1]/tot_t))*c #Initial condition for the within litter prevalence

    N[1]<-min(rbinom(1,pos,prob[1]),tot[1]) #Initial condition for the number of positive piglets in the litter

    

    
    for(i in 2:n){ # Start the binomial recursive model

      prob[i]<-tot[i]/(tot_t-sum(tot[1:i-1]))+(1-(tot[i]/(tot_t-sum(tot[1:i-1]))))*c # Vector of within litter prevalence

      N[i]<-min(rbinom(1,pos-sum(N[1:i-1]),prob[i]),tot[i]) #Vector for the within litter number of positive piglets
      ###

    } #End of the recursive Binomial model


    prev<-data.frame(Prop_pos=N/tot)
    ff = predict(MF2, prev, type = "response", re.form = NA)
    
    for (rt in 1:2000){
      
      x = rbinom(length(ff),1,ff)
      
      yu[rt] = (sum(x))/n
    }
    ALP[j] = mean(yu)
    
    #pos_prop<-rbinom(length(ff),1,ff) # status of each crate in the iteration j
    
    #ALP[j]<-mean(pos_prop) # Proportion of observed positive crates in the iteration j
    
    TLP[j]<-mean(N>0) #Proportion of positive litters in the iteration j
    

    setTxtProgressBar(pb, j)
  
} # End of the stochastic simulation node j
  close(pb)
  
  final[[w]]<-data.frame(ALP=ALP,TPL=TLP)

} #End of w scenarios for piglet prevalence
  
  names(final)<-paste(prevalence*100,"%",sep="") #Names for lists
  
  
  saveRDS(final, file = here("Output","Baseline","Baseline.rds"))

  # Plot of the ALP and TLP in boxplot
  data_names<-paste(prevalence*100,"%",sep="")
  
  
  for(i in 1:length(data_names)){
    PP[[i]]<-rep(data_names[i],nrpl)
  }
  
  PP<-unlist(PP)
  
    data_plot<-do.call(rbind, final)
    data_plot<-cbind.data.frame(data_plot,PP=PP)
  
    data_plot1<-data_plot%>%
      gather(key="A_TLP", value="prev",-PP)
    
data_plot1%>%
    arrange(data_plot1) %>%
      mutate(PP = factor(PP, levels=data_names)) %>%
    subset(PP %in% c("5%","30%","50%","70%","95%" ) )%>%
    ggplot(aes(x=PP,y=prev,fill=A_TLP))+
      theme_minimal()+
      geom_boxplot()+
      xlab("Piglet prevalence")+
      ylab("Litter prevalence")+
      theme(legend.position="right")+
      scale_fill_manual(values=c(4,2),labels=c("Apparent","True"))+
    labs(fill = "Litter prevalence")
  
  ggsave(here("Figures","Baseline","Figure2.png"),width = 7,height = 5,device = "png",dpi=300,bg="white")  #save the plot   
  
  
  
  
#  data_plot1%>%
#    arrange(data_plot1) %>%
#    mutate(PP = factor(PP, levels=data_names)) %>%
#    subset(PP %in% c("5%","30%","50%","70%","95%" ) )%>%
#    ggplot(aes(x=prev,y=PP,fill=A_TLP))+
#    theme_minimal()+
#    geom_density_ridges(alpha = 0.4, scale = 1.5) +
#    ylab("Piglet prevalence")+
#    xlab("Litter prevalence")+
#    theme(legend.position="right")+
#    scale_fill_manual(values=c(4,2),labels=c("Apparent","True"))+
#    labs(fill = "Litter prevalence")
  
   
  
  #Summary of the simulation
  for (i in 1:length(prevalence)){
                                    prev1[i]=prevalence[i]
                                    
                                    ALP_mean[i]=mean(final[[i]]$ALP)
                                    TLP_mean[i]=mean(final[[i]]$TPL)
                                    
                                    ALP2.5[i]=quantile(final[[i]]$ALP,0.025)
                                    ALP_median[i]=quantile(final[[i]]$ALP,0.5)
                                    ALP97.5[i]=quantile(final[[i]]$ALP,0.975)
                                    
                                    TLP2.5[i]=quantile(final[[i]]$TPL,0.025)
                                    TLP_median[i]=quantile(final[[i]]$TPL,0.5)
                                    TLP97.5[i]=quantile(final[[i]]$TPL,0.975)
                                    
                                    
  }
  
  
  plot_data<-cbind.data.frame(prev1,
                              ALP_mean,
                              TLP_mean,
                              ALP2.5,
                              ALP_median,
                              ALP97.5,
                              TLP2.5,
                              TLP_median,
                              TLP97.5
                              )
  capture.output(plot_data, file = here("Output", "Baseline","Plot_data.txt"), append=FALSE) #Saving summary outcome
  
  
  plot_data1<-plot_data[,c(1:3)]%>%
    gather(key="true_false",value="Prevalence",-prev1) #creating dataset for plot
  
  ggplot(plot_data1,aes(x=prev1*100,y=Prevalence,group=true_false))+
    theme_minimal()+
    geom_line(aes(linetype = true_false))+
    xlab("Piglet prevalence")+
    ylab("Litter prevalence")+
    theme(legend.position="right")+
    scale_linetype_manual(values=c(1,2),labels=c("Apparent","True"))+
    labs(linetype = "Litter prevalence")
    
ggsave(here("Figures","Baseline","Figure3.png"),width = 7,height = 5,device = "png",dpi=300,bg="white")  #save the plot

} # End of the function