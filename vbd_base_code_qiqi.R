# If you have not yet installed these packages, you will need to do so, but only one time
#install.packages("ggplot2","proxy","plot.matrix","gridExtra")
rm(list=ls())
# Code for Spatial Model for T-Dependent Mosquito Population 
library(ggplot2) # required for plotting
library(proxy) # required for distance calculations 
library(plot.matrix) # this helps with plotting matrices
library(gridExtra) # this allows for setting plots in a grid 

####################################### Functions for Temperature and Traits ##########
## You should not need to modify anything between this line and line 58, except in  ###
## the exploration of control policies                                      ############

# Function to Control Temperature
temp_func <- function(t,cVec){
  tf = cVec[1] - cVec[2]*cos((2*pi/365)*(t + cVec[3]))
  
}

# Function to output temperature-dependent rates
vecDynamics_params <- function(t, cVec, aVec, jVec, eVec, dd_par, JV){
  
  # calculate temperature from the temperature function
  temp <- temp_func(t, cVec)
  
  # Calculate Density-dependent juvenile mortality fraction
  dd_muJV = 1-exp(-dd_par*JV)
  
  
  # Calculate Adult Lifespan, Juvenile Development Time, and EIP as a function of temperature
  muV0 = 1/adult_lifespan(aVec, temp) # adult mortality rate

  # Turn Rates into Exponential Probabilities
  muV = max(0,1-exp(-muV0))
  nuJV = max(0, 1-exp(-juv_development_rate(jVec, temp)))
  sigma_v = max(0, 1-exp(-extrinsic_incubation_rate(eVec,temp)))
  
  return(c(dd_muJV, muV, nuJV, sigma_v))
  
}

# adult Lifespan function
adult_lifespan <- function(aVec, temp){
  return(aVec[1]*(temp-aVec[2])*(aVec[3]-temp))
}

# juvenile development rate function 
juv_development_rate <-function(jVec, temp) {
  return(jVec[1]*temp*(temp-jVec[2])*(pmax(jVec[3]-temp,0))^(1/2))
}

# incubation rate function
extrinsic_incubation_rate <-function(eVec, temp){
  return(eVec[1]*temp*(temp-eVec[2])*(pmax(eVec[3]-temp,0))^(1/2))
}

#####################################################################################
####################################### Code for Model Implementation ###############

NPOP = 9    # Number of Populations
maxT = 365*2  # Number of Days to run the Model 


# Define initial population sizes
# Population sizes are listed as populations 1-9
#NH = c(5000, 4000, 5000, 
#       2400, 50000, 2400,
#       400, 4500, 4000)

NH = rep.int(5000, NPOP)

# (x,y) coordinates of populations
# This creates a matrix where each row gives (x,y) coordinates for a population. 
xy_mat = matrix(nrow=NPOP, ncol=2,
                data=c(1, 1, 1, 2, 1, 3,
                       2, 1, 2, 2, 2, 3, 
                       3, 1, 3, 2, 3, 3), byrow=T)

# Parameter Data Frame to make access to parameters easier
par <- data.frame(beta_vh = .09, # Probability of transmission from vectors to humans (probability of contact * probability of transmission upon contact)
                  beta_hv = .09, # Probability of transmission from humans to vectors (probability of contact * probability of transmission upon contact)
                  sigma_h = 1-exp(-1/5), # Incubation time of virus in humans, approximately 7 days
                  gamma_h = 1-exp(-1/7), # Duration of infectious period of virus in humans, approximately 7 days
                  alphaV = 8, # Number of new offspring that develop to larvae per female mosquito per day
                  dd_par = 1/50000 # Parameter controlling density dependence
)


# Introduction Information
time_of_intro = seq(1, 181, by = 30)
num_of_intro = 1
pop_of_intro = 5

# Gravity Model Parameters
grav_a = .3 # exponent for the "donor" population
grav_b = .5 # exponent for the "receiving" population
grav_c = 2 # exponent for the distance
grav_theta = 0.000001 # weight controlling contribution of neighboring populations

# Combine Population Numbers and (x,y) coordinates into one data frame
population_matrix = data.frame(NH, xy_mat)

# define names of columns of data frame
names(population_matrix) = c("NH", "x", "y")

# Plot population size to visualize
# Note we plot these population sizes on a log scale to highlight differences
# in population magnitude
plot1 <- ggplot(population_matrix, aes(x=x, y=y, fill=log10(NH)))
plot1 + geom_tile() + ggtitle("Population Size (log scale)")



# This for loop will create a matrix of distances between populations and a matrix for 
# determining movement between populations, which will be used in the model
dist_mat = matrix(0, nrow=NPOP, ncol=NPOP)
move_mat = matrix(0, nrow=NPOP, ncol=NPOP)

for(i in 1:NPOP){
  for(j in 1:NPOP){
    
    # This calculates the distance between two points given their (x,y) coordinates.  
    # we will use Euclidean distance here, but for larger regions, you may want to use
    # the Haversine formula, which takes into account the curvature of the Earth
    dist_mat[i,j] = dist(population_matrix[i,2:3], population_matrix[j,2:3], method="Euclidean")
    
    # Movement matrix from population i to population j
    # This will be driven by the gravity model with parameters a, b, c
    if(i!=j){
      move_mat[i,j] = grav_theta*population_matrix$NH[i]^grav_a*population_matrix$NH[j]^grav_b/(dist_mat[i,j]^grav_c)
    }else{
      move_mat[i,j] = 0
    }
    

  }
}



# Temperature parameter vector
# Written as c(annual average, amplitude, phase shift)
cVec0 = c(21.5, 8.5, 347) # Baseline Temperature

# Matrix for 
cMat = matrix(cVec0,nrow=NPOP, ncol=3, byrow=TRUE)

# Outer Populations all have same temperature
# Center Population is +2 degrees warmer on average
cMat[5,] = c(21.5+2, 8.5, 347)



# Parameters for Temperature-dependent functions
# written as c(rate constant, thermal min, thermal max)
jVec = c(0.00856, 14.58, 34.61) # Juvenile development rate
aVec = c(0.148, 9.16, 37.73) # adult lifespan
eVec = c(.0000665, 10.68, 45.90) # extrinsic incubation rate 
  

# Initiate time
SV_list = list()
EV_list = list()
IV_list = list()
SH_list = list()
EH_list = list()
IH_list = list()
RH_list = list()
new_cases_list = list()

for(k in 1:length(time_of_intro)){
t = 0

# Initiate Storage Matrices
SH = matrix(0, nrow=maxT+1, ncol=NPOP * length(time_of_intro))
EH = SH
IH = SH
RH = SH
JV = SH
SV = SH
EV = SH
IV = SH
lambda_vh = matrix(0, nrow=maxT+1, ncol=NPOP)
lambda_hv = lambda_vh
new_cases = matrix(0, nrow=maxT+1, ncol=NPOP)

vecPars = data.frame(dd_muJV=0, muV=0, nuJV=0, sigma_v=0)

# Initial Conditions for Each Populations - Humans
SH[1,] = NH
EH[1,] = c(0,0,0,0,0,0,0,0,0)
IH[1,] = c(0,0,0,0,0,0,0,0,0)
RH[1,] = c(0,0,0,0,0,0,0,0,0)

#  Initial Conditions for Each Populations - Vectors
maxT_burnin = 365*5
JVb = matrix(0, nrow=maxT_burnin+1, ncol=NPOP)
SVb = matrix(0, nrow=maxT_burnin+1, ncol=NPOP)
JVb[1,] = 10*NH
SVb[1,] = NH

EV[1,] = 0
IV[1,] = 0

NH_sum = sum(SH[1,]+EH[1,]+IH[1,]+RH[1,])


# Quick Burn-in Period to allow for population dynamics of mosquitoes to settle
for(t in 1:(maxT_burnin)){
  for(i in 1:NPOP){
    # Extract temperature-dependent parameters
    vecPars_output = vecDynamics_params(t, cVec0, aVec, jVec, eVec, par$dd_par, JVb[t,i])
    vecPars$dd_muJV = vecPars_output[1]
    vecPars$muV = vecPars_output[2]
    vecPars$nuJV = vecPars_output[3]
    vecPars$sigma_v = vecPars_output[4]
    
    JVb[t+1, i] = JVb[t,i] + par$alphaV*(SVb[t, i]) - (1-vecPars$dd_muJV)*(vecPars$nuJV)*JVb[t, i] - vecPars$dd_muJV*JVb[t, i]
    SVb[t+1, i] = SVb[t, i] + (1-vecPars$dd_muJV)*(vecPars$nuJV)*JVb[t, i] - vecPars$muV*SVb[t, i] 
  }
}

JV[1,] = tail(JVb,1)  
SV[1,] = tail(SVb,1)  



  # Loop for Time
  for(t in 1:maxT){
    
    move_mat_hv = matrix(0,nrow=NPOP, ncol=NPOP)
    move_mat_vh = matrix(0,nrow=NPOP, ncol=NPOP)
    
    # Calculate Transmission Matrices with movement matrix and 
    # Population numbers at time t
    # Note: there are more efficient ways to do this with matrix/vector operations
    for(i in 1:NPOP){
      for(j in 1:NPOP){
  
        # Host to vector Transmission
        # This will be driven by the move_ment matrix move_mat
        # This can be interpreted as the infection of vectors from hosts that are moving around
        # That is, a vector in population i can be infected by hosts traveling from population j
        if(i!=j){
          move_mat_hv[i,j] = move_mat[j,i]*IH[t,j]/NH[i]
        }else{
          move_mat_hv[i,j] = par$beta_hv*IH[t,i]/NH[i]
        }
        
        
        # Vector to Host Transmission
        # This will be driven by the move_ment matrix move_mat
        # This can be interpreted as the infection of hosts from vectors in visited populations
        # That is, a host in population i can be infected by moving and visiting vectors in population j
        if(i!=j){
          move_mat_vh[i,j] = move_mat[i,j]*IV[t,j]/NH[j]
        }else{
          move_mat_vh[i,j] = par$beta_vh*IV[t,i]/NH[i]
        }
        
      }
    }


    # Loop for each population 
    for(i in 1:NPOP){
      
      # Introduce a case
      if(t==time_of_intro[k] && i==pop_of_intro){
        IH[t,i] = IH[t,i] + num_of_intro
        SH[t,i] = SH[t,i] - num_of_intro
      }
      
      
      # Extract temperature-dependent parameters
      cVec = cMat[i,]
      vecPars_output = vecDynamics_params(t, cVec, aVec, jVec, eVec, par$dd_par, JV[t,i])
      vecPars$dd_muJV = vecPars_output[1]
      vecPars$muV = vecPars_output[2]
      vecPars$nuJV = vecPars_output[3]
      vecPars$sigma_v = vecPars_output[4]
      
      # Calculate vector to host transmission when a person in population i visits
      # other populations
      lambda_vh[t,i] = max(0,1-exp(-sum(move_mat_vh[i,])))
      
      # Human Dengue Dynamics 
      SH[t+1, i] = SH[t, i] - lambda_vh[t,i]*SH[t,i]
      EH[t+1, i] = EH[t, i] + lambda_vh[t,i]*SH[t,i] - par$sigma_h*EH[t, i]
      IH[t+1, i] = IH[t, i] + par$sigma_h*EH[t, i] - par$gamma_h*IH[t, i]
      RH[t+1, i] = RH[t, i] + par$gamma_h*IH[t, i]
      
      # Calculate the incidence, or new infections in humans in population i on day t
      new_cases[t+1,i] = par$sigma_h*EH[t, i]
      
      # This equation controls the juvenile population dynamics
      JV[t+1, i] = JV[t,i] + par$alphaV*(SV[t, i]+EV[t, i]+IV[t, i]) - (1-vecPars$dd_muJV)*(vecPars$nuJV)*JV[t, i] - vecPars$dd_muJV*JV[t, i]
      
      # Calculate host to vector transmission from people in other populations visiting
      # population i
      lambda_hv[t,i] = max(0,1-exp(-sum(move_mat_hv[i,])))
      
      # Vector Dengue Dynamics
      SV[t+1, i] = SV[t, i] + (1-vecPars$dd_muJV)*(vecPars$nuJV)*JV[t, i] - (1-vecPars$muV)*lambda_hv[t,i]*SV[t, i] - vecPars$muV*SV[t, i] 
      EV[t+1, i] = EV[t, i] + (1-vecPars$muV)*lambda_hv[t,i]*SV[t, i] - (1-vecPars$muV)*vecPars$sigma_v*EV[t, i] - vecPars$muV*EV[t, i] 
      IV[t+1, i] = IV[t, i] + (1-vecPars$muV)*vecPars$sigma_v*EV[t, i] - vecPars$muV*IV[t, i]
      
      }
        
      
  }
  SV_list[[k]] <- SV
  EV_list[[k]] <- EV
  IV_list[[k]] <- IV
  SH_list[[k]] <- SH
  EH_list[[k]] <- EH
  IH_list[[k]] <- IH
  RH_list[[k]] <- RH
  new_cases_list[[k]] <- new_cases
}

### Plotting to Visualize Model Output ##########################################
for (k in 1:length(time_of_intro)){
  SV_list[[k]] -> SV
  EV_list[[k]] -> EV
  IV_list[[k]] -> IV
  SH_list[[k]] -> SH
  EH_list[[k]] -> EH
  IH_list[[k]] -> IH
  RH_list[[k]] -> RH
  new_cases_list[[k]] -> new_cases
  # Combine Infection Numbers into a data frame for easier plotting
  infections = data.frame(0:maxT, IH)
  names(infections) = c('day','I1', 'I2', 'I3', 'I4', 'I5', 'I6', 'I7', 'I8', 'I9')
  
  # Combine Adult Vector Numbers into a data frame for easier plotting - Vector to Host Ratio
  vectors = data.frame(0:maxT, SV+EV+IV)
  names(vectors) = c('day','V1', 'V2', 'V3', 'V4', 'V5', 'V6', 'V7', 'V8', 'V9')
  
  # Combine incidence in a vector for plotting
  incidence = data.frame(0:maxT, new_cases)
  names(incidence) = c('day','I1', 'I2', 'I3', 'I4', 'I5', 'I6', 'I7', 'I8', 'I9')
  
  # Calculate total number of cases
  outbreak_size = colSums(new_cases)
  
  # Find Time of first case in each population and other metrics of interest
  # which() returns the index of the vector at which the argument is true
  # For example which(cases>1) would return all of the indices of cases where the value was
  # greater than 1. 
  time_first_new_case = numeric(length=NPOP)
  peak_case_number = numeric(length=NPOP)
  time_peak_cases = numeric(length=NPOP)
  for(i in 1:NPOP){
    time_first_new_case[i] = min(which(new_cases[,i]>=1))-1 # subtract 1 because we start at t=0
    peak_case_number[i] = max(new_cases[,i])
    time_peak_cases[i] = min(which(new_cases[,i]==peak_case_number[i]))
  }
  
  epidemic_metrics = data.frame(NH, 1:9, xy_mat, outbreak_size, time_first_new_case, peak_case_number, time_peak_cases)
  names(epidemic_metrics)[2]= 'popNum'
  names(epidemic_metrics)[3]= 'x'
  names(epidemic_metrics)[4]= 'y'
  
  
  # # Visualize Temperature-dependent rates
  # timevec = seq(0,365, by=1)
  # tempvec = seq(0,40,by=.1)
  # 
  # temp_plot = temp_func(timevec,cVec0)
  # alife_plot = pmax(0,adult_lifespan(aVec,tempvec))
  # jd_plot = pmax(0,juv_development_rate(jVec,tempvec))
  # eip_plot = pmax(0,extrinsic_incubation_rate(eVec,tempvec))
  # 
  # time_temp = data.frame(timevec, temp_plot)
  # temp_traits = data.frame(tempvec, alife_plot, jd_plot, eip_plot)
  # 
  # r1 <- ggplot(data=time_temp, aes(x=timevec, y=temp_plot)) + geom_line() + ylab("Temperature")
  # r2 <- ggplot(data=temp_traits, aes(x=tempvec, y=jd_plot)) + geom_line() + ylab("Juvenile Development Rate")
  # r3 <- ggplot(data=temp_traits, aes(x=tempvec, y=alife_plot)) + geom_line() + ylab("Adult Lifespan")
  # r4 <- ggplot(data=temp_traits, aes(x=tempvec, y=eip_plot)) + geom_line() + ylab("Extrinisc Incubation Rate")
  # 
  # grid.arrange(r1,r2,r3,r4, nrow=2, ncol=2)
  
  
  ###
  # temp_plot_mat = matrix(0, ncol=9, nrow=length(0:maxT))
  # 
  # for(i in 1:NPOP){
  #   temp_plot_mat[,i] = temp_func(timevec,cMat[i,])
  # }
  # 
  # tpm_df = data.frame(timevec, temp_plot_mat)
  # names(tpm_df) = c("time", "pop1", "pop2", "pop3", "pop4","pop5", "pop6", "pop7", "pop8", "pop9")
  # 
  # r1 <- ggplot(data=tpm_df, aes(x=time, y=pop1)) + geom_line() + ylab("Temp, Population 1") + ylim(c(0,40))
  # r2 <- ggplot(data=tpm_df, aes(x=time, y=pop2)) + geom_line() + ylab("Temp, Population 2") + ylim(c(0,40))
  # r3 <- ggplot(data=tpm_df, aes(x=time, y=pop3)) + geom_line() + ylab("Temp, Population 3") + ylim(c(0,40))
  # 
  # r4 <- ggplot(data=tpm_df, aes(x=time, y=pop4)) + geom_line() + ylab("Temp, Population 4") + ylim(c(0,40))
  # r5 <- ggplot(data=tpm_df, aes(x=time, y=pop5)) + geom_line() + ylab("Temp, Population 5") + ylim(c(0,40))
  # r6 <- ggplot(data=tpm_df, aes(x=time, y=pop6)) + geom_line() + ylab("Temp, Population 6") + ylim(c(0,40))
  # 
  # r7 <- ggplot(data=tpm_df, aes(x=time, y=pop7)) + geom_line() + ylab("Temp, Population 7") + ylim(c(0,40))
  # r8 <- ggplot(data=tpm_df, aes(x=time, y=pop8)) + geom_line() + ylab("Temp, Population 8") + ylim(c(0,40))
  # r9 <- ggplot(data=tpm_df, aes(x=time, y=pop9)) + geom_line() + ylab("Temp, Population 9") + ylim(c(0,40))
  # 
  # grid.arrange(r7, r8, r9, r4, r5, r6, r1, r2, r3, nrow=3, ncol=3)
  
  #################################### Visualize Other Epidemic Metrics
  # p1 <- ggplot(epidemic_metrics, aes(popNum, outbreak_size/NH*1000)) + geom_col(fill="darkblue") + 
  #   scale_x_discrete(name="Population", limits = factor(epidemic_metrics$popNum)) +theme_bw() +ylab("Total Incidence (per 1000 people)")
  # 
  # p2 <- ggplot(epidemic_metrics, aes(popNum, time_first_new_case)) + geom_col(fill="mediumblue") + 
  #   scale_x_discrete(name="Population", limits = factor(epidemic_metrics$popNum)) +theme_bw() +ylab("Time of First Case")
  # 
  # p3 <- ggplot(epidemic_metrics, aes(popNum, peak_case_number/NH*1000)) + geom_col(fill="violet") + 
  #   scale_x_discrete(name="Population", limits = factor(epidemic_metrics$popNum)) +theme_bw() +ylab("Peak Incidence (per 1000 people)")
  # 
  # p4 <- ggplot(epidemic_metrics, aes(popNum, time_peak_cases)) + geom_col(fill="darkviolet") + 
  #   scale_x_discrete(name="Population", limits = factor(epidemic_metrics$popNum)) +theme_bw() +ylab("Time of Peak")
  # 
  # grid.arrange(p1,p2,p3,p4, nrow=2, ncol=2, top= paste0("Introduction time: Day ",time_of_intro[k]))
  # 
  # 
  # #################################### Visualize epidemic metrics in a tile grid
  # p1 <- ggplot(epidemic_metrics, aes(x=x,y=y, fill=outbreak_size/NH*1000)) + geom_tile(color="black") + 
  #   ggtitle("Total Incidence (per 1000 people)") + scale_fill_gradientn(colors=hcl.colors(7,"Reds")) + 
  #   coord_fixed() + guides(fill=guide_colourbar(title='', ticks=TRUE))
  # 
  # p2 <- ggplot(epidemic_metrics, aes(x=x,y=y, fill=time_first_new_case)) + geom_tile(color="black")+ 
  #   ggtitle("Time to first case") + scale_fill_gradientn(colors=hcl.colors(7,"Hawaii")) +
  #   coord_fixed() + guides(fill=guide_colourbar(title='', ticks=TRUE))
  # 
  # p3 <- ggplot(epidemic_metrics, aes(x=x,y=y, fill=peak_case_number/NH*1000)) + geom_tile(color="black") + 
  #   ggtitle("Peak Incidence (per 1000 people)") + scale_fill_gradientn(colors=hcl.colors(7,"PurpOr")) +
  #   coord_fixed() + guides(fill=guide_colourbar(title='', ticks=TRUE))
  # 
  # p4 <- ggplot(epidemic_metrics, aes(x=x,y=y, fill=time_peak_cases)) + geom_tile(color="black") + 
  #   ggtitle("Time of Peak") + scale_fill_gradientn(colors=hcl.colors(7,"Batlow"))+ 
  #   coord_fixed() + guides(fill=guide_colourbar(title='', ticks=TRUE))
  # 
  # grid.arrange(p1,p2,p3,p4, nrow=2, ncol=2, top= paste0("Introduction time: Day ",time_of_intro[k]))
  # 
  #################################### Visual of Cases at different time points
  # incidencet01 = data.frame(cases=t(t(new_cases[50,])/NH*1000), x=xy_mat[,1], y=xy_mat[,2])
  # incidencet02 = data.frame(cases=t(t(new_cases[100,])/NH*1000), x=xy_mat[,1], y=xy_mat[,2])
  # incidencet03 = data.frame(cases=t(t(new_cases[150,])/NH*1000), x=xy_mat[,1], y=xy_mat[,2])
  # incidencet04 = data.frame(cases=t(t(new_cases[200,])/NH*1000), x=xy_mat[,1], y=xy_mat[,2])
  # incidencet05 = data.frame(cases=t(t(new_cases[250,])/NH*1000), x=xy_mat[,1], y=xy_mat[,2])
  # incidencet06 = data.frame(cases=t(t(new_cases[300,])/NH*1000), x=xy_mat[,1], y=xy_mat[,2])
  # 
  # maxINC = max(c(incidencet01$cases,incidencet02$cases,incidencet03$cases,incidencet04$cases,
  #                incidencet05$cases,incidencet06$cases))
  # 
  # p1 <- ggplot(incidencet01, aes(x=x,y=y, fill=cases)) + geom_tile(color="black") + 
  #   ggtitle("t=50") + scale_fill_gradientn(colors=hcl.colors(7,"PurpOr"),limits=c(0,maxINC)) + 
  #   coord_fixed() + guides(fill=guide_colourbar(title='', ticks=TRUE)) 
  # 
  # 
  # p2 <- ggplot(incidencet02, aes(x=x,y=y, fill=cases)) + geom_tile(color="black") + 
  #   ggtitle("t=100") + scale_fill_gradientn(colors=hcl.colors(7,"PurpOr"),limits=c(0,maxINC)) + 
  #   coord_fixed() + guides(fill=guide_colourbar(title='', ticks=TRUE)) 
  # 
  # p3 <- ggplot(incidencet03, aes(x=x,y=y, fill=cases)) + geom_tile(color="black") + 
  #   ggtitle("t=150") + scale_fill_gradientn(colors=hcl.colors(7,"PurpOr"),limits=c(0,maxINC)) + 
  #   coord_fixed() + guides(fill=guide_colourbar(title='', ticks=TRUE))
  # 
  # p4 <- ggplot(incidencet04, aes(x=x,y=y, fill=cases)) + geom_tile(color="black") + 
  #   ggtitle("t=200") + scale_fill_gradientn(colors=hcl.colors(7,"PurpOr"),limits=c(0,maxINC)) + 
  #   coord_fixed() + guides(fill=guide_colourbar(title='', ticks=TRUE)) 
  # 
  # p5 <- ggplot(incidencet05, aes(x=x,y=y, fill=cases)) + geom_tile(color="black") + 
  #   ggtitle("t=250") + scale_fill_gradientn(colors=hcl.colors(7,"PurpOr"),limits=c(0,maxINC)) + 
  #   coord_fixed() + guides(fill=guide_colourbar(title='', ticks=TRUE))
  # 
  # p6 <- ggplot(incidencet06, aes(x=x,y=y, fill=cases)) + geom_tile(color="black") + 
  #   ggtitle("t=300") + scale_fill_gradientn(colors=hcl.colors(7,"PurpOr"),limits=c(0,maxINC)) + 
  #   coord_fixed() + guides(fill=guide_colourbar(title='', ticks=TRUE)) 
  # 
  # grid.arrange(p1,p2,p3,p4,p5, p6, nrow=2, ncol=3, top=paste0("Introduction time: Day ",
  #                                                             time_of_intro[k],
  #                                                             " - Incidence (per 1000 people)"))
  # 
  
  #################################### Visualize Model Output - Daily Incidence
  r1 <- ggplot(data=incidence, aes(x=day, y=I1/NH[1]*1000)) + geom_line() + ylab('incidence') +
    ggtitle('population 1') + theme(plot.title=element_text(size=10))
  r2 <- ggplot(data=incidence, aes(x=day, y=I2/NH[2]*1000)) + geom_line() + ylab('incidence') +
    ggtitle('population 2') + theme(plot.title=element_text(size=10))
  r3 <- ggplot(data=incidence, aes(x=day, y=I3/NH[3]*1000)) + geom_line() + ylab('incidence') +
    ggtitle('population 3') + theme(plot.title=element_text(size=10))
  r4 <- ggplot(data=incidence, aes(x=day, y=I4/NH[4]*1000)) + geom_line() + ylab('incidence') +
    ggtitle('population 4') + theme(plot.title=element_text(size=10))
  r5 <- ggplot(data=incidence, aes(x=day, y=I5/NH[5]*1000)) + geom_line() + ylab('incidence') +
    ggtitle('population 5') + theme(plot.title=element_text(size=10))
  r6 <- ggplot(data=incidence, aes(x=day, y=I6/NH[6]*1000)) + geom_line() +  ylab('incidence') +
    ggtitle('population 6') + theme(plot.title=element_text(size=10))
  r7 <- ggplot(data=incidence, aes(x=day, y=I7/NH[7]*1000)) + geom_line() +  ylab('incidence') +
    ggtitle('population 7') + theme(plot.title=element_text(size=10))
  r8 <- ggplot(data=incidence, aes(x=day, y=I8/NH[8]*1000)) + geom_line() +  ylab('incidence') +
    ggtitle('population 8') + theme(plot.title=element_text(size=10))
  r9 <- ggplot(data=incidence, aes(x=day, y=I9/NH[9]*1000)) + geom_line() +  ylab('incidence') +
    ggtitle('population 9') + theme(plot.title=element_text(size=10))
  
  p <- grid.arrange(r1,r2,r3,r4,r5,r6,r7,r8,r9, nrow=3, ncol=3, top=paste0("Introduction time: Day ",
                                                                      time_of_intro[k],
                                                                      " - Incidence (per 1000 people)"))
  ggsave(paste0("inci_d",time_of_intro[k],".png"), p)
}

# reference: https://www.nagraj.net/notes/gifs-in-r/

## list file names and read in
library(magick)
imgs <- list.files("/Users/qiqiyang/Documents/GitHub/qiqi-yang-eeid/inci", full.names = TRUE)
img_list <- lapply(imgs, image_read)

## join the images together
img_joined <- image_join(img_list)

## animate at 2 frames per second
img_animated <- image_animate(img_joined, fps = 1)

## view animated image
img_animated

## save to disk
image_write(image = img_animated,
            path = "/Users/qiqiyang/Documents/GitHub/qiqi-yang-eeid/inci_intro_time.gif")


# # Visualize Model Output - Vector - Host Ratio 
# r1 <- ggplot(data=vectors, aes(x=day, y=V1/NH[1])) + geom_line() + ylim(0, 1.025*max(vectors$V1/NH[1])) + ylab("VHR, Population 1")
# r2 <- ggplot(data=vectors, aes(x=day, y=V2/NH[2])) + geom_line() + ylim(0, 1.025*max(vectors$V2/NH[2])) + ylab("VHR, Population 2")
# r3 <- ggplot(data=vectors, aes(x=day, y=V3/NH[3])) + geom_line() + ylim(0, 1.025*max(vectors$V3/NH[3])) + ylab("VHR, Population 3")
# 
# r4 <- ggplot(data=vectors, aes(x=day, y=V4/NH[4])) + geom_line() + ylim(0, 1.025*max(vectors$V4/NH[4])) + ylab("VHR, Population 4")
# r5 <- ggplot(data=vectors, aes(x=day, y=V5/NH[5])) + geom_line() + ylim(0, 1.025*max(vectors$V5/NH[5])) + ylab("VHR, Population 5")
# r6 <- ggplot(data=vectors, aes(x=day, y=V6/NH[6])) + geom_line() + ylim(0, 1.025*max(vectors$V6/NH[6])) + ylab("VHR, Population 6")
# 
# r7 <- ggplot(data=vectors, aes(x=day, y=V7/NH[7])) + geom_line() + ylim(0, 1.025*max(vectors$V7/NH[7])) + ylab("VHR, Population 7")
# r8 <- ggplot(data=vectors, aes(x=day, y=V8/NH[8])) + geom_line()+ ylim(0, 1.025*max(vectors$V8/NH[8])) + ylab("VHR, Population 8")
# r9 <- ggplot(data=vectors, aes(x=day, y=V9/NH[9])) + geom_line() + ylim(0, 1.025*max(vectors$V9/NH[9])) + ylab("VHR, Population 9")
# grid.arrange(r1,r2,r3,r4,r5,r6,r7,r8,r9, nrow=3, ncol=3)
# 



