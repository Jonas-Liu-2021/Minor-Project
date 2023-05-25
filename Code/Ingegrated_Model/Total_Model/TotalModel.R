library(deSolve)


# Build the differential equation for enzyme: dE/dt = ks*mRNA-(kd+u)*E
EnzymeRates <- function(parms,starts){
  x=vector()
  num_enzyme=parms[1]
  # ################################################################################################
  # Change enzyme_params of enzyme synthesis as needed
  # The parameters are in this order: enzyme_parms(ks1,mRNA1,kd1,ks2,mRNA2,kd2,...)
  # ################################################################################################
  enzyme_parms=c(0.045,277.4,0.0478,
                 0.045,269.5,  0.02,
                 0.045, 85.9,  0.02,
                 0.045, 65.8,  0.02,
                 0.045, 91.4,  0.02,
                 0.045,323.0,  0.02,
                 0.045,261.9,  0.02,
                 0.045,137.8,  0.02,
                 0.045,126.2,  0.02,
                 0.045, 1.83,  0.02,
                 0.045,30.46,  0.02,
                 0.045,101.5,  0.02,
                 0.045, 68.0,  0.02,
                 0.045, 59.3,  0.02,
                 0.045, 32.7,  0.02,
                 0.045, 22.2,  0.02)
  for (i in 1:num_enzyme) {
    enzyme_parms[3*i-1] <- enzyme_parms[3*i-1]/1
  }
  
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
  reaction_parms=c(206.23, 930,
                   157630,3000,
                   4798.3, 932,
                   36000, 100,
                   36000, 100,
                   16200,1900,
                   1192.4,68.6,
                   56692, 1920,
                   1702.8,  34,
                   1500, 556,
                   36000, 100,
                   3600, 5.8,
                   87.96,  70,
                   63.01,  88,
                   402.5,  55,
                   350289, 640,
                   4954, 476,
                   3701.7,  34,
                   70,  80,
                   82.53,  65,
                   274.3,  43,
                   9507, 463,
                   88.8,  37,
                   4639, 358,
                   148.06, 194,
                   1206, 253,
                   59.23, 216,
                   36.86, 374,
                   102.7,  81,
                   29.61,  74,
                   36.86,  26,
                   62000, 1920,
                   7000, 1920,
                   100, 1920,
                   20,1920,
                   0.2,1920,
                   3.1,1920)
  # The equations should be listed one by one since some substrates are catalyzed by several enzymes
  # First part
  x[1] <- reaction_parms[1]*starts[1]*starts[parms[1]+1]/(reaction_parms[2]+starts[parms[1]+1])
  x[2] <- reaction_parms[3]*starts[2]*starts[parms[1]+2]/(reaction_parms[4]+starts[parms[1]+2])
  x[3] <- reaction_parms[5]*starts[3]*starts[parms[1]+2]/(reaction_parms[6]+starts[parms[1]+2])
  x[4] <- reaction_parms[7]*starts[4]*starts[parms[1]+3]/(reaction_parms[8]+starts[parms[1]+3])
  x[5] <- reaction_parms[9]*starts[5]*starts[parms[1]+3]/(reaction_parms[10]+starts[parms[1]+3])
  x[6] <- reaction_parms[11]*starts[6]*starts[parms[1]+3]/(reaction_parms[12]+starts[parms[1]+3])
  x[7] <- reaction_parms[13]*starts[7]*starts[parms[1]+4]/(reaction_parms[14]+starts[parms[1]+4])
  x[8] <- reaction_parms[63]*starts[8]*starts[parms[1]+5]/(reaction_parms[16]+starts[parms[1]+5])
  x[9] <- reaction_parms[17]*starts[9]*starts[parms[1]+6]/(reaction_parms[18]+starts[parms[1]+6])
  x[10] <- reaction_parms[19]*starts[10]*starts[parms[1]+7]/(reaction_parms[20]+starts[parms[1]+7])
  x[11] <- reaction_parms[21]*starts[11]*starts[parms[1]+8]/(reaction_parms[22]+starts[parms[1]+9])
  x[12] <- reaction_parms[23]*starts[12]*starts[parms[1]+9]/(reaction_parms[24]+starts[parms[1]+9])
  x[13] <- reaction_parms[25]*starts[13]*starts[parms[1]+10]/(reaction_parms[26]+starts[parms[1]+10])
  x[14] <- reaction_parms[27]*starts[14]*starts[parms[1]+10]/(reaction_parms[28]+starts[parms[1]+10])
  x[15] <- reaction_parms[29]*starts[15]*starts[parms[1]+10]/(reaction_parms[30]+starts[parms[1]+10])
  # Second part
  x[16] <- reaction_parms[31]*starts[2]*starts[parms[1]+5]/(reaction_parms[32]+starts[parms[1]+5])
  x[17] <- reaction_parms[33]*starts[3]*starts[parms[1]+5]/(reaction_parms[34]+starts[parms[1]+5])
  x[18] <- reaction_parms[7]*starts[4]*starts[parms[1]+12]/(reaction_parms[8]+starts[parms[1]+12])
  x[19] <- reaction_parms[9]*starts[5]*starts[parms[1]+12]/(reaction_parms[10]+starts[parms[1]+12])
  x[20] <- reaction_parms[11]*starts[6]*starts[parms[1]+12]/(reaction_parms[12]+starts[parms[1]+12])         
  x[21] <- reaction_parms[13]*starts[7]*starts[parms[1]+13]/(reaction_parms[14]+starts[parms[1]+13])
  x[22] <- reaction_parms[65]*starts[8]*starts[parms[1]+14]/(reaction_parms[16]+starts[parms[1]+14])
  x[23] <- reaction_parms[35]*starts[9]*starts[parms[1]+15]/(reaction_parms[36]+starts[parms[1]+15])
  x[24] <- reaction_parms[19]*starts[10]*starts[parms[1]+16]/(reaction_parms[20]+starts[parms[1]+16])
  x[25] <- reaction_parms[21]*starts[11]*starts[parms[1]+17]/(reaction_parms[22]+starts[parms[1]+17])
  x[26] <- reaction_parms[23]*starts[12]*starts[parms[1]+18]/(reaction_parms[24]+starts[parms[1]+18])
  x[27] <- reaction_parms[37]*starts[13]*starts[parms[1]+19]/(reaction_parms[38]+starts[parms[1]+19])
  x[28] <- reaction_parms[39]*starts[14]*starts[parms[1]+19]/(reaction_parms[40]+starts[parms[1]+19])
  x[29] <- reaction_parms[41]*starts[15]*starts[parms[1]+19]/(reaction_parms[42]+starts[parms[1]+19])
  # Third part
  x[30] <- reaction_parms[43]*starts[3]*starts[parms[1]+14]/(reaction_parms[44]+starts[parms[1]+14])
  x[31] <- reaction_parms[7]*starts[4]*starts[parms[1]+21]/(reaction_parms[8]+starts[parms[1]+21])
  x[32] <- reaction_parms[9]*starts[5]*starts[parms[1]+21]/(reaction_parms[10]+starts[parms[1]+21])
  x[33] <- reaction_parms[11]*starts[6]*starts[parms[1]+21]/(reaction_parms[12]+starts[parms[1]+21])
  x[34] <- reaction_parms[13]*starts[7]*starts[parms[1]+22]/(reaction_parms[14]+starts[parms[1]+22])
  x[35] <- reaction_parms[67]*starts[8]*starts[parms[1]+23]/(reaction_parms[16]+starts[parms[1]+23])
  x[36] <- reaction_parms[45]*starts[9]*starts[parms[1]+24]/(reaction_parms[46]+starts[parms[1]+24])
  x[37] <- reaction_parms[19]*starts[10]*starts[parms[1]+25]/(reaction_parms[20]+starts[parms[1]+25])
  x[38] <- reaction_parms[21]*starts[11]*starts[parms[1]+26]/(reaction_parms[22]+starts[parms[1]+26])
  x[39] <- reaction_parms[23]*starts[12]*starts[parms[1]+27]/(reaction_parms[24]+starts[parms[1]+27])
  x[40] <- reaction_parms[37]*starts[13]*starts[parms[1]+28]/(reaction_parms[38]+starts[parms[1]+28])
  x[41] <- reaction_parms[39]*starts[14]*starts[parms[1]+28]/(reaction_parms[40]+starts[parms[1]+28])
  x[42] <- reaction_parms[41]*starts[15]*starts[parms[1]+28]/(reaction_parms[42]+starts[parms[1]+28])
  # Fourth part
  x[43] <- reaction_parms[47]*starts[3]*starts[parms[1]+23]/(reaction_parms[48]+starts[parms[1]+23])
  x[44] <- reaction_parms[7]*starts[4]*starts[parms[1]+30]/(reaction_parms[8]+starts[parms[1]+30])
  x[45] <- reaction_parms[9]*starts[5]*starts[parms[1]+30]/(reaction_parms[10]+starts[parms[1]+30])
  x[46] <- reaction_parms[11]*starts[6]*starts[parms[1]+30]/(reaction_parms[12]+starts[parms[1]+30])
  x[47] <- reaction_parms[13]*starts[7]*starts[parms[1]+31]/(reaction_parms[14]+starts[parms[1]+31])
  x[48] <- reaction_parms[69]*starts[8]*starts[parms[1]+32]/(reaction_parms[16]+starts[parms[1]+32])
  x[49] <- reaction_parms[49]*starts[9]*starts[parms[1]+33]/(reaction_parms[50]+starts[parms[1]+33])
  x[50] <- reaction_parms[19]*starts[10]*starts[parms[1]+34]/(reaction_parms[20]+starts[parms[1]+34])
  x[51] <- reaction_parms[21]*starts[11]*starts[parms[1]+35]/(reaction_parms[22]+starts[parms[1]+35])
  x[52] <- reaction_parms[23]*starts[12]*starts[parms[1]+36]/(reaction_parms[24]+starts[parms[1]+36])
  x[53] <- reaction_parms[37]*starts[13]*starts[parms[1]+37]/(reaction_parms[38]+starts[parms[1]+37])
  x[54] <- reaction_parms[39]*starts[14]*starts[parms[1]+37]/(reaction_parms[40]+starts[parms[1]+37])
  x[55] <- reaction_parms[41]*starts[15]*starts[parms[1]+37]/(reaction_parms[42]+starts[parms[1]+37])
  # Fifth part
  x[56] <- reaction_parms[51]*starts[3]*starts[parms[1]+32]/(reaction_parms[52]+starts[parms[1]+32])
  x[57] <- reaction_parms[7]*starts[4]*starts[parms[1]+39]/(reaction_parms[8]+starts[parms[1]+39])
  x[58] <- reaction_parms[9]*starts[5]*starts[parms[1]+39]/(reaction_parms[10]+starts[parms[1]+39])
  x[59] <- reaction_parms[11]*starts[6]*starts[parms[1]+39]/(reaction_parms[12]+starts[parms[1]+39])
  x[60] <- reaction_parms[13]*starts[7]*starts[parms[1]+40]/(reaction_parms[14]+starts[parms[1]+40])
  x[61] <- reaction_parms[71]*starts[8]*starts[parms[1]+41]/(reaction_parms[16]+starts[parms[1]+41])
  x[62] <- reaction_parms[53]*starts[9]*starts[parms[1]+42]/(reaction_parms[54]+starts[parms[1]+42])
  x[63] <- reaction_parms[55]*starts[16]*starts[parms[1]+42]/(reaction_parms[56]+starts[parms[1]+42])
  x[64] <- reaction_parms[19]*starts[10]*starts[parms[1]+43]/(reaction_parms[20]+starts[parms[1]+43])
  x[65] <- reaction_parms[21]*starts[11]*starts[parms[1]+44]/(reaction_parms[22]+starts[parms[1]+44])
  x[66] <- reaction_parms[23]*starts[12]*starts[parms[1]+45]/(reaction_parms[24]+starts[parms[1]+45])
  x[67] <- reaction_parms[39]*starts[14]*starts[parms[1]+46]/(reaction_parms[40]+starts[parms[1]+46])
  x[68] <- reaction_parms[41]*starts[15]*starts[parms[1]+46]/(reaction_parms[42]+starts[parms[1]+46])
  # Sixth part
  x[69] <- reaction_parms[57]*starts[3]*starts[parms[1]+41]/(reaction_parms[58]+starts[parms[1]+41])
  x[70] <- reaction_parms[7]*starts[4]*starts[parms[1]+48]/(reaction_parms[8]+starts[parms[1]+48])
  x[71] <- reaction_parms[9]*starts[5]*starts[parms[1]+48]/(reaction_parms[10]+starts[parms[1]+48])
  x[72] <- reaction_parms[11]*starts[6]*starts[parms[1]+48]/(reaction_parms[12]+starts[parms[1]+48])
  x[73] <- reaction_parms[13]*starts[7]*starts[parms[1]+49]/(reaction_parms[14]+starts[parms[1]+49])
  x[74] <- reaction_parms[73]*starts[8]*starts[parms[1]+50]/(reaction_parms[16]+starts[parms[1]+50])
  x[75] <- reaction_parms[59]*starts[9]*starts[parms[1]+51]/(reaction_parms[60]+starts[parms[1]+51])
  x[76] <- reaction_parms[61]*starts[16]*starts[parms[1]+51]/(reaction_parms[62]+starts[parms[1]+51])
  x[77] <- reaction_parms[19]*starts[10]*starts[parms[1]+52]/(reaction_parms[20]+starts[parms[1]+52])
  x[78] <- reaction_parms[21]*starts[11]*starts[parms[1]+53]/(reaction_parms[22]+starts[parms[1]+53])
  x[79] <- reaction_parms[23]*starts[12]*starts[parms[1]+54]/(reaction_parms[24]+starts[parms[1]+54])
  x[80] <- reaction_parms[39]*starts[14]*starts[parms[1]+55]/(reaction_parms[40]+starts[parms[1]+55])
  x[81] <- reaction_parms[41]*starts[15]*starts[parms[1]+55]/(reaction_parms[42]+starts[parms[1]+55])
  
  
  
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
  S=as.matrix(read.csv("StoichiometryMatrix.csv",header = FALSE))
  
  # The rate of metabolite synthesis
  #################################################################################################
  # The function is dM/dt = Reaction Rate - Dilution Rate - Degradation Rate (for end metabolite)
  # Change the metabolite vector as needed (Only the initial substrate does not have dilution rate. 
  # Only the end metabolite have degradation rate)
  #################################################################################################
  dMdt=S%*%ReactionRate(parms,starts)-parms[3]*c(0,starts[(parms[1]+2):(parms[1]+parms[2])])-
    parms[4]*c(rep(0,10),starts[(parms[1]+11)],rep(0,8),starts[(parms[1]+20)],rep(0,8),starts[(parms[1]+29)],
               rep(0,8),starts[(parms[1]+38)],rep(0,8),starts[(parms[1]+47)],rep(0,8),starts[(parms[1]+56)])
  
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
  tspan <- seq(0,400,by=1)
  parms <- c(16,56,0.015,0.0000007)
  
  starts <- c(rep(0,parms[1]),785,rep(0,parms[2]-1))
  
  out <- ode(starts,tspan,MetNetwork,parms)
  
  # Print the output. You could delete it if you want to
  out <<- out
  print(out)
  print(c(out[400,c((parms[1]+12),(parms[1]+21),(parms[1]+30),(parms[1]+39),(parms[1]+48),(parms[1]+57))],sum(out[300,c((parms[1]+12),(parms[1]+21),(parms[1]+30),(parms[1]+39),(parms[1]+48),(parms[1]+57))])))
  print(out[400,18:19])
  
  # Plot the result
  num_enzyme=parms[1]
  num_metabolite=parms[2]
  windows()
  par(mfrow=c(2,2))
  matplot(out[,1],out[,2:(num_enzyme+1)],type="l", ylab="Concentration", 
          xlab="Time", main="Enzyme Concentration", lwd = 2)   
  matplot(out[,1],out[,(num_enzyme+2):(num_enzyme+num_metabolite+1)],type ="l",ylab ="Relative Concentration",
          xlab = "Time",main="Metabolite Concentration",lwd = 2)
  matplot(out[,1],out[,c((num_enzyme+2):66,68:(num_enzyme+num_metabolite+1))],type ="l",ylab ="Relative Concentration",
          xlab = "Time",main="Metabolite Concentration",lwd = 2)
  barplot(c(out[400,c((parms[1]+12),(parms[1]+21),(parms[1]+30),(parms[1]+39),(parms[1]+48),(parms[1]+57))],sum(out[300,c((parms[1]+12),(parms[1]+21),(parms[1]+30),(parms[1]+39),(parms[1]+48),(parms[1]+57))])),
          main="Metabolite level",ylab ="Concentration of Metabolites",names.arg=c("C3","C4","C5","C6","C7","C8","Total GSLs"))
  windows()
  barplot(c(out[400,c((parms[1]+12),(parms[1]+21),(parms[1]+30),(parms[1]+39),(parms[1]+48),(parms[1]+57))],sum(out[300,c((parms[1]+12),(parms[1]+21),(parms[1]+30),(parms[1]+39),(parms[1]+48),(parms[1]+57))])),
          ylab ="Relative Concentration of GSLs",names.arg=c("C3","C4","C5","C6","C7","C8","Total GSLs"))
  windows()
  matplot(out[,1],out[,c((parms[1]+12),(parms[1]+21),(parms[1]+30),(parms[1]+39),(parms[1]+48),(parms[1]+57))],type = "l",
          ylab ="Relative Concentration of GSLs")
}


# For test
Concentration()