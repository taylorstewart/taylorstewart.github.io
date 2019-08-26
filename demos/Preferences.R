####
###Author: Andrew Nguyen
### Start Date: 2017-04-05
### Last modified: 2017-04-05

#Source script to include preferences for plotting and 
# fitting non-linear functions

##Ggplot theme
T<-theme_bw()+theme(text=element_text(size=30),axis.text=element_text(size=30),
                    legend.text=element_text(size=28),panel.grid.major=element_blank(),
                    legend.position="none",panel.grid.minor.x = element_blank(),
                    panel.grid = element_blank(),legend.key = element_blank())

########################################
########################################

### Gaussian function for plotting
gaus<-function(a=2,b=5,c=30,x=seq(-5,25,.05)){
  #y=a*exp(((x-b)^2)/(2c^2))
  y=a*exp(-5*(((x-b))/(c))^2)
  return(y)
}
########################################
########################################

########################################
########################################

#####fitting data with boltzmann function
Boltz<-function(data=x){
  B<-nls(gxp ~ (1+(max-1)/(1+exp((Tm-temperature)/a))),data=data, start=list(max=80,Tm=35,a=1.05), trace=TRUE,control=nls.control(warnOnly = TRUE, tol = 1e-05, maxiter=1000))
  #summary(B)
  return(summary(B)$parameters)
}

##Plotting boltzman
###creates data based off of the parameters of a function
boltzmann<-function(temperature=seq(25,70,.1),Tm=40,slope=1.8,max=50){
  y<-1+ (max-1)/(1+exp(((Tm-temperature)/slope)))
  return(y)
}

###gene expression dataset
temperature<-c(25,28,30,31.5,33,35,36.5,40,41)
gxp<-c(1.139380725,
       1.495138067,
       1.31816746,
       2.39787468,
       3.341707929,
       6.387151393,
       6.266289656,
       11.39597512,
       11.99697887)
dat<-as.data.frame(cbind(temperature,gxp))#;dat

########################################
########################################
#### Unfolding curves
UFfun<-function(data=data){
  y<-nls(unfolding ~ min+ (1-min)/(1+exp((-slope*(Tm-temperature)))),data=data, 
           start=list(slope=.5,Tm=45,min=.3),
           trace=TRUE,control=nls.control(warnOnly = TRUE, tol = 1e-05, maxiter=1000))
  #return(y)
  return(summary(y)$coefficients)
}

###Plotting unfolding
ufold<-function(temperature=seq(25,70,1),Tm=40,slope=.5,max=1,min=0){
  y<-min+ (max-min)/(1+exp((-slope*(Tm-temperature))))
  return(y)
}

########################################
########################################

###
