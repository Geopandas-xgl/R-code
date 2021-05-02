# carbon cycle OMZ model
# by Kimberly Lau

library(seacarb)

#### model variables, constants, and initial conditions ####
#constants, modified from Payne and Kump
G = 3.3 #Berner (2004)
Z = 0.09 #sil weathering, Berner (2004)
Zc = 0.07 #carb weathering, Berner (2004)
frac.org = -29.5
frac.carb = 0#
oceanV = 1.4e21 #L of ocean
dCvolc = -5
dCwcarb = 2.54
dCworg = -23.19 

#initial reservoirs and variables, modified from Payne and Kump
PO4 <- vector()
DIC <- vector()
d13C <- vector()
Alk <- vector()
kcp <- vector()
fanox <- vector()
PO4.i = 3e15 #mol initial P
PO4[1] = PO4.i 
DIC[1] = 7.42e18 #mol initial DIC
d13C[1] = 2.54 #initial DIC permil
Alk[1] = 8.06e18 #inital Alk
kcp[1] = 106 #initial C:P ratio
fanox.i = 0.002
fanox[1] = fanox.i
korg = 2935384.61538462
kcarb = 2.8e13

#inital carbonate system, modified from ABJ CCaPS
omega <- vector()
pCO2 <- vector()
RCO2 <- vector()
depth=0 #meters
temp.i = 20 #initial ocean surface temp in °C
salinity = 35 #salinity in salinity units
sal = (salinity/1000)+1 #salinity in kg/L
Ca.i = 15e-3 #3.5e-3 #initial calcium conc. in mol/kg
carb.parms.i = pCa(flag=15, var1=Alk[1]/oceanV, var2=DIC[1]/oceanV, Ca=Ca.i, S=salinity, T=temp.i, P=0, Pt=PO4[1]/oceanV) #abj_co2sys(DIC[1]/oceanV, Alk[1]/oceanV, temp.i, salinity, depth, Ca.i) #calculates the initial seawater carbonate parameters using the seacarb package 
print(carb.parms.i)
omega.i = carb.parms.i$OmegaAragonite[2] #initial omega
omega[1]=omega.i
pCO2.i = carb.parms.i$pCO2[2] #initial pCO2
pCO2[1]=pCO2.i
RCO2[1] = pCO2.i/pCO2.i
pH.i = carb.parms.i$pH[2]

#Phosphate fluxes, modified from Payne and Kump
Fwp <- vector()
Fbp <- vector()
Fwp[1] = 1.3e11*RCO2[1]^(G*Z)*(1+G*Z*log(RCO2[1]))^0.65 #weathered PO4, mol/yr
Fbp[1] = Fwp[1]*(PO4.i/PO4.i) #buried PO4, mol/yr, initial balanced
  
#Carbon fluxes, modified from Payne and Kump
Fwcarb <- vector()
Fborg <- vector()
Fbcarb <- vector()
Fvolc = 4e12 #volcanic CO2, mol/yr
Fworg = 8e12 #organic C weathering, mol/yr
Fwcarb[1] = 24e12*RCO2[1]^(G*Zc)*(1+G*Zc*log(RCO2[1])) #carbonate weathering, mol/yr
Fborg.i=8e12 #burial of organic carbon, mol/yr
Fborg[1] = Fborg.i#7.2e12 #*(fanox[i]/fanox.i)^(1/3)*(Fbp[i-1]/Fbp[1])*(kcp[i-1]/kcp[1]) #organic C burial, mol/yr
Fbcarb[1] = kcarb*(omega.i/omega.i) #carbonate burial, mol/yr

#Alk fluxes, modified from Payne and Kump
Fwsil <- vector()
Falkin <- vector()
Falkout <- vector()
Fwsil[1] = 8e12*RCO2[1]^(G*Z)*(1+G*Z*log(RCO2[1]))^0.65 #silicate weathering, mol/yr # is 8e12 in Payne and Kump
Falkin[1] = 2*Fwcarb[1]+Fwsil[1]
Falkout[1] = 2*Fbcarb[1]

#### model set up ####
t.end = 5e6 #length of model run in yrs 
dt = 10000#timestep; interval at which model will be evaluated in yrs 
t = seq(0,t.end,dt) #times when model will be evaluated
Zub <- vector()
dampen <- vector()
lambda <- vector()
lambda.Tr <- vector()
for (i in 1:length(t)) {
  lambda[i]=40
  Zub[i] = -50*cos(pi*i/lambda[i])
}

DIC.ans=array(dim=c(length(t),3))
Alk.ans=array(dim=c(length(t),3))
pH.ans=array(dim=c(length(t),3))
pCO2.ans=array(dim=c(length(t),3))
d13C.ans=array(dim=c(length(t),3))
PO4.ans=array(dim=c(length(t),3))
Fborg.ans=array(dim=c(length(t),3))
Fbcarb.ans=array(dim=c(length(t),3))
kcp.ans=array(dim=c(length(t),3))
fanox.ans=array(dim=c(length(t),3))
Fwsil.ans=array(dim=c(length(t),3))
Fwcarb.ans=array(dim=c(length(t),3))
omega.ans=array(dim=c(length(t),3))
Fbp.ans=array(dim=c(length(t),3))

DIC.ans[1,] = DIC[1]
Alk.ans[1,] = Alk[1]
pH.ans[1,] = pH.i
pCO2.ans[1,] = pCO2.i
d13C.ans[1,] = d13C[1]
PO4.ans[1,] = PO4[1]
Fborg.ans[1,] = Fborg[1]
Fbcarb.ans[1,] = Fbcarb[1]
kcp.ans[1,] = kcp[1]
fanox.ans[1,] = fanox.i
Fwsil.ans[1,] = Fwsil[1]
Fwcarb.ans[1,] = Fwcarb[1]
omega.ans[1,] = omega.i
Fbp.ans[1,] = Fbp[1]

#### sea level relationships ####
fanox1 <- vector()
fanox1[1] = fanox.i
fanox2 <- vector()
fanox2[1] = fanox.i
fanox3 <- vector()
fanox3[1] = fanox.i
fanox4 <- vector()
fanox4[1] = fanox.i
for (i in 2:100) {
  fanox1[i] = (0.2/100)*(i) #shallow OMZ
  fanox2[i] = (0.01/100)*(i) #deep OMZ
  fanox3[i]=fanox.i #steady state
}

#### model equations ####
for (j in 1:3) {
  for (i in 2:length(t)) {
    dampen[i]=1
    if (j==1) { 
      fanox[i] = dampen[i]*(0.2/100)*(Zub[i]+50)+fanox.i #shallow OMZ
    } else if (j==2) {
      fanox[i] = dampen[i]*(0.01/100)*(Zub[i]+50)+fanox.i #deep OMZ
    } else if (j==3) {
      fanox[i] = fanox.i #steady state
    } 
  
    #fluxes
    Fborg[i] = Fborg.i*(Fbp[i-1]/Fbp[1])*(kcp[i-1]/kcp[1]) #organic C burial, mol/yr
    Fwcarb[i] = 24e12*RCO2[i-1]^(G*Zc)*(1+G*Zc*log(RCO2[i-1])) #carbonate weathering, mol/yr
    Fbcarb[i] = kcarb*(omega[i-1]/omega.i) #carbonate burial, mol/yr
    Fwp[i] = 1.3e11*RCO2[i-1]^(G*Z)*(1+G*Z*log(RCO2[i-1]))^0.65 #weathered PO4, mol/yr
    Fwsil[i] = 8e12*RCO2[i-1]^(G*Z)*(1+G*Z*log(RCO2[i-1]))^0.65 #silicate weathering, mol/yr
    Falkin[i] = 2*Fwcarb[i]+Fwsil[i]
    Falkout[i] = 2*Fbcarb[i]
    
    #reservoirs
    Alk[i] = Alk[i-1] + (Falkin[i] - Falkout[i])*dt
    DIC[i] = DIC[i-1] + (Fvolc+Fworg+Fwcarb[i]-Fborg[i]-Fbcarb[i])*dt
    PO4[i] = PO4[i-1] + (Fwp[i]-Fbp[i-1])*dt
    
    #carb system params
    carb.parms = pCa(flag=15, var1=Alk[i]/oceanV, var2=DIC[i]/oceanV,Ca=Ca.i,S=salinity,T=temp.i,P=0,Pt=PO4[i]/oceanV)
    pCO2[i] = carb.parms$pCO2[2]
    pH.temp = carb.parms$pH[2]
    omega[i] = carb.parms$OmegaAragonite[2]
    RCO2[i] = pCO2[i]/pCO2.i

    #other variables
    kcp[i] = kcp[1]*(fanox[i]/fanox.i)^(1/3)*PO4[i]/PO4.i #C:P ratio, as the phosphate pool gets larger, Fborg goes up, presumably O2 goes down, less iron oxides which phosphate can adsorp onto
    Fbp[i] = Fwp[1]*(PO4[i]/PO4.i) #buried PO4, mol/kyr, initial balanced
    Iso.volc = Fvolc*(dCvolc - d13C[i-1])/DIC[i]
    Iso.wcarb = Fwcarb[i]*(dCwcarb - d13C[i-1])/DIC[i]
    Iso.worg = Fworg*(dCworg - d13C[i-1])/DIC[i]
    Iso.bcarb = Fbcarb[i]*frac.carb/DIC[i]
    Iso.borg = Fborg[i]*(frac.org)/DIC[i]
    d13C[i] = d13C[i-1] + (Iso.volc+Iso.wcarb+Iso.worg-Iso.bcarb-Iso.borg)*dt
    
    DIC.ans[i,j] = DIC[i]
    Alk.ans[i,j] = Alk[i]
    pH.ans[i,j] = pH.temp
    pCO2.ans[i,j] = pCO2[i]
    d13C.ans[i,j] = d13C[i]
    PO4.ans[i,j] = PO4[i]
    Fborg.ans[i,j] = Fborg[i]
    Fbcarb.ans[i,j] = Fbcarb[i]
    kcp.ans[i,j] = kcp[i]
    fanox.ans[i,j] = fanox[i]
    Fwsil.ans[i,j] = Fwsil[i]
    Fwcarb.ans[i,j] = Fwcarb[i]
    omega.ans[i,j] = omega[i]
    Fbp.ans[i,j] = Fbp[i]
  }
}

#### plots ####
#### panel graph for figure ####
par(mfrow=c(5,1), mar=c(2,4,1,1))

plot(fanox1*100,type="l",col="dark gray",xlab="Sea level (m)",ylab="fanox (%)",lwd=2)
lines(fanox3*100,col="dark green",lwd=2,lty="dotted")
lines(fanox2*100,col="blue",lwd=2,lty="dashed")
legend("topleft",c("Shallow", "Deep"),col=c("dark gray","blue"),lty=c("solid","dashed"),lwd=2,bty="n")

plot(t/1e6, Zub+50, type="l",xlab="Model time (Myr)",ylab="Relative change in OMZ_ub (m)",lwd=2)

plot(t/1e6, fanox.ans[,1]*100,type="l",col="dark gray",xlab="Model time (Myr)",ylab="fanox (%)",lwd=2)
lines(t/1e6, fanox.ans[,3]*100, col="dark green",lwd=2,lty="dotted")
lines(t/1e6, fanox.ans[,2]*100, col="blue",lwd=2,lty="dashed")

plot(t/1e6, Fborg.ans[,1],type="l",col="dark gray",xlab="Model time (Myr)",ylab="Fborg (mol/yr)",lwd=2)
lines(t/1e6, Fborg.ans[,3], col="dark green",lwd=2,lty="dotted")
lines(t/1e6, Fborg.ans[,2], col="blue",lwd=2,lty="dashed")

plot(t/1e6, d13C.ans[,1],type="l",col="dark gray",xlab="Model time (Myr)",ylab="d13C",lwd=2)
lines(t/1e6, d13C.ans[,3],col="dark green",lwd=2,lty="dotted")
lines(t/1e6, d13C.ans[,2],col="blue",lwd=2,lty="dashed")

#### all model output ####
if(False){
par(mfrow=c(5,2), mar=c(2,4,1,1))
plot(t/1e6, Zub, type="l",xlab="Model time (Myr)",ylab="Sea level change (m)",lwd=2)

plot(t/1e6, Fborg.ans[,1]/(Fborg.ans[,1]+Fbcarb.ans[,1]),type="l",xlab="Model time (Myr)",ylab="forg",lwd=2)
lines(t/1e6, Fborg.ans[,2]/(Fborg.ans[,2]+Fbcarb.ans[,2]),col="blue",lwd=2,lty="dashed")
lines(t/1e6, Fborg.ans[,3]/(Fborg.ans[,3]+Fbcarb.ans[,3]),col="dark green",lwd=2,lty="dotted")

plot(t/1e6, Fbcarb.ans[,1],type="l",xlab="Model time (Myr)",ylab="Fbcarb",lwd=2)
lines(t/1e6, Fbcarb.ans[,2],col="blue",lwd=2,lty="dashed")
lines(t/1e6, Fbcarb.ans[,3],col="dark green",lwd=2,lty="dotted")

plot(t/1e6, DIC.ans[,1],type="l",xlab="Model time (Myr)",ylab="DIC (mol)",lwd=2)
lines(t/1e6, DIC.ans[,2], col="blue",lwd=2,lty="dashed")
lines(t/1e6, DIC.ans[,3], col="dark green",lwd=2,lty="dotted")

plot(t/1e6, Alk.ans[,1],type="l",xlab="Model time (Myr)",ylab="Alk (mol)",lwd=2)
lines(t/1e6, Alk.ans[,2], col="blue",lwd=2,lty="dashed")
lines(t/1e6, Alk.ans[,3], col="dark green",lwd=2,lty="dotted")

plot(t/1e6, pCO2.ans[,1],type="l",xlab="Model time (Myr)",ylab="pCO2 (ppm)",lwd=2)
lines(t/1e6, pCO2.ans[,2], col="blue",lwd=2,lty="dashed")
lines(t/1e6, pCO2.ans[,3], col="dark green",lwd=2,lty="dotted")

plot(t/1e6, Fborg.ans[,1],type="l",xlab="Model time (Myr)",ylab="Fborg (mol/yr)",lwd=2)
lines(t/1e6, Fborg.ans[,2], col="blue",lwd=2,lty="dashed")
lines(t/1e6, Fborg.ans[,3], col="dark green",lwd=2,lty="dotted")

plot(t/1e6, kcp.ans[,1],type="l",xlab="Model time (Myr)",ylab="kcp",lwd=2)
lines(t/1e6, kcp.ans[,2], col="blue",lwd=2,lty="dashed")
lines(t/1e6, kcp.ans[,3], col="dark green",lwd=2,lty="dotted")

plot(t/1e6, PO4.ans[,1],type="l",xlab="Model time (Myr)",ylab="PO4",lwd=2)
lines(t/1e6, PO4.ans[,2], col="blue",lwd=2,lty="dashed")
lines(t/1e6, PO4.ans[,3], col="dark green",lwd=2,lty="dotted")

plot(t/1e6, Fbp.ans[,1],type="l",xlab="Model time (Myr)",ylab="Fbp",lwd=2)
lines(t/1e6, Fbp.ans[,2], col="blue",lwd=2,lty="dashed")
lines(t/1e6, Fbp.ans[,3], col="dark green",lwd=2,lty="dotted")
}

#### Triassic matching ####
if(false){
t.end.Tr = 10e6 #length of model run in yrs
dt.Tr = 10000#timestep; interval at which model will be evaluated in yrs 
t.Tr = seq(0,t.end.Tr,dt.Tr) #times when model will evaluated
Zub.Tr <- vector()
dampen <- vector()
lambda.Tr <- vector()
DIC.ans.Tr <- vector()
Alk.ans.Tr <- vector()
pH.ans.Tr <- vector()
pCO2.ans.Tr <- vector()
d13C.ans.Tr <- vector()
PO4.ans.Tr <- vector()
Fborg.ans.Tr <- vector()
Fbcarb.ans.Tr <- vector()
kcp.ans.Tr <- vector()
fanox.ans.Tr <- vector()
Fwsil.ans.Tr <- vector()
Fwcarb.ans.Tr <- vector()
omega.ans.Tr <- vector()
Fbp.ans.Tr <- vector()

for (i in 1:length(t.Tr)) {
  lambda.Tr[i] = 35+i/20
  Zub.Tr[i] = -50*cos(pi*i/lambda.Tr[i])
}

for (i in 2:length(t.Tr)) {
  if (i <= 500) {
    dampen[i] = (length(t.Tr)-1.5*i)/(length(t.Tr)) }
  else if (i >500) {
    dampen[i]=exp(-(i-500)/100)*dampen[500]
  }
  fanox[i] = dampen[i]*(0.2/100)*(Zub.Tr[i]+50)+fanox.i #shallow
  
  #fluxes
  Fborg[i] = 7.2e12*(fanox[i]/fanox.i)^(1/3)*(Fbp[i-1]/Fbp[1])*(kcp[i-1]/kcp[1]) #organic C burial, mol/yr
  Fwcarb[i] = 24e12*RCO2[i-1]^(G*Z)*(1+G*Zc*log(RCO2[i-1])) #carbonate weathering, mol/yr
  Fbcarb[i] = kcarb*(omega[i-1]/omega.i) #carbonate burial, mol/yr
  Fwp[i] = 1.3e11*RCO2[i-1]^(G*Z)*(1+G*Z*log(RCO2[i-1]))^0.65 #weathered PO4, mol/yr
  Fwsil[i] = 9.6e12*RCO2[i-1]^(G*Z)*(1+G*Z*log(RCO2[i-1]))^0.65 #silicate weathering, mol/yr
  Falkin[i] = 2*Fwcarb[i]+Fwsil[i]
  Falkout[i] = 2*Fbcarb[i]
  
  #reservoirs
  Alk[i] = Alk[i-1] + (Falkin[i] - Falkout[i])*dt
  DIC[i] = DIC[i-1] + (Fvolc+Fworg+Fwcarb[i]-Fborg[i]-Fbcarb[i])*dt
  PO4[i] = PO4[i-1] + (Fwp[i]-Fbp[i-1])*dt
  
  #carb system params
  carb.parms = pCa(flag=15, var1=Alk[i]/oceanV, var2=DIC[i]/oceanV,Ca=Ca.i,S=salinity,T=temp.i,P=0,Pt=PO4[i]/oceanV)
  pCO2[i] = carb.parms$pCO2[2]
  pH.temp = carb.parms$pH[2]
  omega[i] = carb.parms$OmegaAragonite[2]
  RCO2[i] = pCO2[i]/pCO2.i
  
  #other variables
  kcp[i] = kcp[i-1]*PO4[i]/PO4.i #C:P ratio, as the phosphate pool gets larger, Fborg goes up, presumably O2 goes down, less iron oxides which phosphate can adsorp onto
  Fbp[i] = Fwp[i]*(PO4[i]/PO4.i) #buried PO4, mol/kyr, initial balanced
  Iso.volc = Fvolc*(dCvolc - d13C[i-1])/DIC[i]
  Iso.wcarb = Fwcarb[i]*(dCwcarb - d13C[i-1])/DIC[i]
  Iso.worg = Fworg*(dCworg - d13C[i-1])/DIC[i]
  Iso.bcarb = Fbcarb[i]*frac.carb/DIC[i]
  Iso.borg = Fborg[i]*(frac.org)/DIC[i]
  d13C[i] = d13C[i-1] + (Iso.volc+Iso.wcarb+Iso.worg-Iso.bcarb-Iso.borg)*dt
  
  DIC.ans.Tr[i] = DIC[i]
  Alk.ans.Tr[i] = Alk[i]
  pH.ans.Tr[i] = pH.temp
  pCO2.ans.Tr[i] = pCO2[i]
  d13C.ans.Tr[i] = d13C[i]
  PO4.ans.Tr[i] = PO4[i]
  Fborg.ans.Tr[i] = Fborg[i]
  Fbcarb.ans.Tr[i] = Fbcarb[i]
  kcp.ans.Tr[i] = kcp[i]
  fanox.ans.Tr[i] = fanox[i]
  Fwsil.ans.Tr[i] = Fwsil[i]
  Fwcarb.ans.Tr[i] = Fwcarb[i]
  omega.ans.Tr[i] = omega[i]
  Fbp.ans.Tr[i] = Fbp[i]
}

#### Triassic model output ####
par(mfrow=c(4,1), mar=c(2,4,1,1))
plot(fanox.ans.Tr*100,type="l",xlab="Sea level (m)",ylab="fanox (%)",lwd=2)
plot(t.Tr/1e6, Zub.Tr+50-100*(1-dampen), type="l",xlab="Model time (Myr)",ylab="Depth of OMZ_ub relative to SSB (m)",lwd=2)
plot(t.Tr/1e6, fanox.ans.Tr*100,type="l",xlab="Model time (Myr)",ylab="fanox (%)",lwd=2)
plot(t.Tr/1e6, d13C.ans.Tr,type="l",xlab="Model time (Myr)",ylab="d13C",lwd=2)
}
