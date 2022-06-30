

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
  
  
  sce<-crossing(sce_n,sce_c)
  sce_n<-sce$sce_n
  sce_c<-sce$sce_c
  
  
  
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
  
  
  
  #Creating labels for crates
  crates<-list()

  lab1<-unique(sce_n)
  for (i in 1:length(lab1)){
  crates[[i]]<-c(rep(lab1[i],length(lab1)*length(prevalence)))
   }
  crates<-unlist(crates)
  
  #Creating labels for cluster
  cluster<-list()
  
  lab2<-unique(sce_c)
  
  for (i in 1:length(lab2)){
    cluster[[i]]<-c(rep(lab2[i],length(prevalence)))
  }
  
  cluster<-rep(cluster,length(lab2))
  
  cluster<-unlist(cluster)
  
  
  #Get the simulation results
  lit_prev1<-unlist(ALP_me)
  lit_prev2<-unlist(TLP_me)
  
  #Labels for APL or TLP
  LP1<-rep("ALP",length(sce_n)*length(sce_c)*length(prevalence))
  LP2<-rep("TLP",length(sce_n)*length(sce_c)*length(prevalence))
  
  #Piglet prevalence
  pig_prev<-rep(prevalence,length(sce_n)*length(sce_c))
  
  
  #creating one dataset to ALP and another to TLP
  data0<-cbind.data.frame(lit_prev=lit_prev1,crates,cluster,LP=LP1,pig_prev)
  
  data1<-cbind.data.frame(lit_prev=lit_prev2,crates,cluster,LP=LP2,pig_prev)
  
  #Combinig both
  data<-rbind.data.frame(data0,data1)
  
  ####Save the plot#####
 #plotp<-function(){
   plotp<-ggplot(data, aes(x=pig_prev,y=lit_prev,color=LP))+
    theme_minimal()+
    facet_grid(crates~cluster,)+
    geom_line()+
    ylab("Litter prevalence")+
    xlab("Piglet prevalence")+
    theme(legend.title=element_blank(),
          axis.text.x.top = element_blank(),       # do not show top / right axis labels
          axis.text.y.right = element_blank(),     # for secondary axis
          axis.title.x.top = element_blank(),      # as above, don't show axis titles for
          axis.title.y.right = element_blank())    # secondary axis either)
  
 # Labels 
 labelR = "Number of crates"
 labelT = "Clustering level"
 
 # Get the ggplot grob
 z <- ggplotGrob(plotp)
 
 # Get the positions of the strips in the gtable: t = top, l = left, ...
 posR <- subset(z$layout, grepl("strip-r", name), select = t:r)
 posT <- subset(z$layout, grepl("strip-t", name), select = t:r)
 
 # Add a new column to the right of current right strips, 
 # and a new row on top of current top strips
 width <- z$widths[max(posR$r)]    # width of current right strips
 height <- z$heights[min(posT$t)]  # height of current top strips
 
 z <- gtable_add_cols(z, width, max(posR$r))  
 z <- gtable_add_rows(z, height, min(posT$t)-1)
 
 # Construct the new strip grobs
 stripR <- gTree(name = "Strip_right", children = gList(
   rectGrob(gp = gpar(col = NA, fill = "grey85")),
   textGrob(labelR, rot = -90, gp = gpar(fontsize = 9.0, col = "grey10"))))
 
 stripT <- gTree(name = "Strip_top", children = gList(
   rectGrob(gp = gpar(col = NA, fill = "grey85")),
   textGrob(labelT, gp = gpar(fontsize = 9.0, col = "grey10"))))
 
 # Position the grobs in the gtable
 z <- gtable_add_grob(z, stripR, t = min(posR$t)+1, l = max(posR$r) + 1, b = max(posR$b)+1, name = "strip-right")
 z <- gtable_add_grob(z, stripT, t = min(posT$t), l = min(posT$l), r = max(posT$r), name = "strip-top")
 
 # Add small gaps between strips
 z <- gtable_add_cols(z, unit(1/5, "line"), max(posR$r))
 z <- gtable_add_rows(z, unit(1/5, "line"), min(posT$t))
 
 # Draw it
  grid.newpage()
 g<-grid.arrange(z)
 
   ggsave(here("Figures","Scenarios","Figure2.png"),width = 10,height = 5,device = "png",dpi=300,bg = "white",g)  #save the plot
  
  
  #Sva the dataframe
  write.table(data,file=here("Output", "Scenarios","Scenarios.txt"),sep=",")
  
} # End of the function
