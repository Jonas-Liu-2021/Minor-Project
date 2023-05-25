library(deSolve)


# Build the differential equation for enzyme: dE/dt = ks*mRNA-(kd+u)*E
EnzymeRates <- function(parms,starts){
  x=vector()
  num_enzyme=parms[1]
  # ################################################################################################
  # Change enzyme_params of enzyme synthesis as needed
  # The parameters are in this order: enzyme_parms(ks1,mRNA1,kd1,ks2,mRNA2,kd2,...)
  # ################################################################################################
  enzyme_parms=c(0.045,85.87,0.02,
                 0.045,126.16,0.02,
                 0.045,22.15,0.02,
                 0.045,269.49,0.02)
  for (i in 1:num_enzyme) {
    x[i] <- enzyme_parms[3*i-2]*enzyme_parms[3*i-1]-(enzyme_parms[3*i]+parms[3])*starts[i]
  }
  
  return(x)
}


# Build the reaction rate equation: mostly Michaelis-Menten v = kc*E*M/(kM+M)
ReactionRate <- function(parms,starts){
  x=vector()
  #################################################################################################
  # Change reaction_params of reaction rate as needed
  # The parameters are in this order: reaction_parms(kc1,KM1,kc2,KM2,...)
  #################################################################################################
  reaction_parms=c(4798.3,932,
                   4954,476,
                   9507,463,
                   4639,358,
                   1206,253,
                   102.7, 81,
                   1702.8, 11.3,
                   3701.7, 34,
                   88.8, 37,
                   148.1,194,
                   59.2,216,
                   29.6, 74,
                   36.9,374,
                   36.9, 26,
                   157630,3000,
                   350289,640)
  # The equations should be listed one by one since some substrates are catalyzed by several enzymes
  x[1] <- reaction_parms[1]*starts[1]*starts[5]/(reaction_parms[2]+starts[5])
  x[2] <- reaction_parms[3]*starts[1]*starts[6]/(reaction_parms[4]+starts[6])
  x[3] <- reaction_parms[5]*starts[1]*starts[7]/(reaction_parms[6]+starts[7])
  x[4] <- reaction_parms[7]*starts[1]*starts[8]/(reaction_parms[8]+starts[8])
  x[5] <- reaction_parms[9]*starts[1]*starts[9]/(reaction_parms[10]+starts[9])
  x[6] <- reaction_parms[11]*starts[1]*starts[10]/(reaction_parms[12]+starts[10])
  x[7] <- reaction_parms[13]*starts[2]*starts[6]/(reaction_parms[14]+starts[6])
  x[8] <- reaction_parms[15]*starts[2]*starts[7]/(reaction_parms[16]+starts[7])
  x[9] <- reaction_parms[17]*starts[2]*starts[8]/(reaction_parms[18]+starts[8])
  x[10] <- reaction_parms[19]*starts[2]*starts[9]/(reaction_parms[20]+starts[9])
  x[11] <- reaction_parms[21]*starts[2]*starts[10]/(reaction_parms[22]+starts[10])
  x[12] <- reaction_parms[23]*starts[2]*starts[11]/(reaction_parms[24]+starts[11])
  x[13] <- reaction_parms[25]*starts[3]*starts[10]/(reaction_parms[26]+starts[10])
  x[14] <- reaction_parms[27]*starts[3]*starts[11]/(reaction_parms[28]+starts[11])
  x[15] <- reaction_parms[29]*starts[4]*starts[5]/(reaction_parms[30]+starts[5])
  x[16] <- reaction_parms[31]*starts[4]*starts[6]/(reaction_parms[32]+starts[6])
  
  return(x)
}



# Build the function including all ODE functions.
MetNetwork <- function(t, starts, parms) {
  
  # The rate of enzyme synthesis
  dEdt= EnzymeRates(parms,starts)
  
  #################################################################################################
  # Change the stoichiometry matrix as needed
  # The stoichiometric matrix is in this format: S=rbind(c(-1, 0, 0), c(1, -1, 0),c(0, 2, -1),...)
  #################################################################################################
  S=rbind(c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 
          c( 1,-1, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 1,-1),
          c( 0, 1,-1, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 1), 
          c( 0, 0, 1,-1, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0),
          c( 0, 0, 0, 1,-1, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0),
          c( 0, 0, 0, 0, 1,-1, 0, 0, 0, 0,-1, 0,-1, 0, 0, 0),
          c( 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,-1, 0,-1, 0, 0),
          c( 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
          c( 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0),
          c( 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0),
          c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0),
          c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0),
          c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0))
  
  # The rate of metabolite synthesis
  #################################################################################################
  # The function is dM/dt = Reaction Rate - Dilution Rate - Degradation Rate (for end metabolite)
  # Change the metabolite vector as needed (Only the initial substrate does not have dilution rate. 
  # Only the end metabolite have degradation rate)
  #################################################################################################
  dMdt=S%*%ReactionRate(parms,starts)-parms[3]*c(0,starts[6:17])-parms[4]*c(rep(0,7),starts[12:17])
  
  print(list(c(dEdt,dMdt)))
  list(c(dEdt,dMdt))
}



# Build the function which receive parameters and plot the result
Concentration <- function() {
  #################################################################################################
  # Change the tspan, params and starts as needed
  # The parameters are in this order: parms(enzyme number, metabolite number, dilution rate, 
  # degradation rate for end metabolite)
  # The initial states are in this order: starts(E1,E2,...,M1,M2,...)
  #################################################################################################
  tspan <- seq(0,400,by=0.1)
  parms <- c(4,13,0.015,0.0000007)
  starts <- c(0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0)
  
  out <- ode(starts,tspan,MetNetwork,parms)
  
  # Print the output. You could delete it if you want to
  print(out)
  print(out[4000,2:5])
  print(c(out[4000,13:18],sum(out[4000,13:18])))
  # Plot the result
  num_enzyme=parms[1]
  num_metabolite=parms[2]
  lenStarts=parms[1]
  
  windows()
  
  matplot(out[,1],out[,2:(num_enzyme+1)],type="l", ylab="Relative Enzyme Concentration", 
          xlab="Time", main="Relative Enzyme Concentration vs Time", lwd = 1,lty = 1,col=1:4)   
  legend("topleft",c("MAM3","CYP79F1","CYP79F2","MAM1"),lty = 1,col = 1:4,xjust = 0,
         bty="n",x.intersp=1)
  windows()
  matplot(out[,1],out[,13:18],type ="l",ylab ="Relative Metabolite Concentration",
          xlab = "Time",main="Relative Metabolite Concentration vs Time",lwd = 2,lty = 1,col=1:13)
  legend(0,1000000,c("C3","C4","C5","C6","C7","C8"),lty = 1,col = 1:13,xjust = 0,
         bty="n",x.intersp=1)
  windows()
  matplot(out[1:30,1],out[1:30,6:12],type ="l",ylab ="Relative Metabolite Concentration",
          xlab = "Time",main="Relative Metabolite Concentration vs Time",lwd = 2,lty = 1,col=1:13)
  legend(0,11,c("MTOB","MTOP","MTOH","MTOHp","MTOO","MTON","MTOD"),lty = 1,col = 1:13,xjust = 0,
         bty="n",x.intersp=1)

  
}


# For test
Concentration()

