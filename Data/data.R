
#Set parameters
prevalence = c(0.01,seq(0.05,1,0.05)) #List of possible prevalences
sce_n<-c(10,30,50,70,56,56,56,56)
sce_c<-c(0.63,0.63,0.63,0.63,0.01,0.2,0.7,0.9)
n = 56
c = 0.63

# Creating objects
p<-NULL
prob<-NULL
N<-NULL
TLP<-NULL
teste<-NULL
ff <- NULL
prev<-NULL
ALP<-NULL
TLP<- NULL
final<-list()

prev1<-NULL
ALP_mean<-NULL
TLP_mean<-NULL
ALP2.5<-NULL
ALP_median<-NULL
ALP97.5<-NULL
TLP2.5<-NULL
TLP_median<-NULL
TLP97.5<-NULL

ALP_me<-vector(mode = "list", length = length(sce_n))
TLP_me<-vector(mode = "list", length = length(sce_c))

