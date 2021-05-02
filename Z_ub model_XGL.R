library(seacarb)
G = 3.3
Z = 0.09
Zc = 0.07
frac.org = -29.5
frac.carb = 0
oceanV = 1.4e21
dCvolc = -5
dCwcarb = 2.54
dCworg = -23.19

# initial reservoirs and variables, modified from Payne and Kump
PO4 <- vector()
DIC <- vector()
d13C <- vector()
Alk <- vector()
kcp <- vector()
fanox <- vector()
PO4.i = 3e15
PO4[1] = PO4.i
DIC[1] = 7.42e18
d13C[1] = 2.54
Alk[1] = 8.06e18
kcp[1] = 106
fanox.i = 0.002
fanox[1] = fanox.i
korg = 2935384.61538462
kcarb = 2.8e13

#initial carbonate system
omega <-  vector()
pCO2 <-  vector()
RCO2 <-  vector()
depth = 0
temp.i = 20
salinity = 35
sal = (salinity/1000)+1
Ca.i = 15e-3
carb.parms.i = pCa(flag=15, var1=Alk[1]/oceanV, var2 = DIC[1]/oceanV, Ca=Ca.i, S=salinity, T=temp.i, P=0, Pt=PO4[1]/oceanV)
print(carb.parms.i)

omega.i = carb.parms.i$OmegaAragonite[2] 
omega[1] = omega.i
pCO2.i = carb.parms.i$pCO2[2]
pCO2[1] = pCO2.i
RCO2[1] = pCO2.i/pCO2.i
pH.i = carb.parms.i$pH[2]

Fwp <- vector()
Fbp <- vector()
Fwp[1] = 1.3e11*RCO2[1]^(G*Z)*(1+G*Z*log(RCO2[1]))^0.65
Fbp[1] = Fwp[1]*(PO4.i/PO4.i)

Fwcarb <-  vector()
Fborg <- vector()
Fbcarb <-  vector()
Fvolc = 4e12
Fworg = 8e12
Fwcarb[1] = 24e12*RCO2[1]^(G*Zc)*(1+G*Zc*log(RCO2[1]))
Fborg[1] = Fborg.i
Fbcarb[1] = kcarb*(omega.i/omega.i)

Fwsil <- vector()
Falkin <-  vector()
Falkout <-  vector()
Fwsil[1] = 8e12*RCO2[1]^(G*Z)*(1+G*Z*log(RCO2[1]))^0.65
Falkin[1] = 2*Fwcarb[1]+Fwsil[1]
Falkout[1] = 2*Fbcarb[1]

# model set up
t.end = 5e6
dt = 10000
t = seq(0,t.end,dt)

Zub <- vector()
dampen <- vector()
lambda <-  vector()
lambda.Tr <- vector()

for (i in 1:length(t)) {
  lambda[i]=40
  Zub[i] = -50*cos(pi*i/lambda[i])
}

DIC.ans = array(dim=c(length(t),3))
Alk.ans = array(dim=c(length(t),3))
pH.ans = array(dim=c(length(t),3))
pCO2.ans = array(dim=c(length(t),3))
d13C.ans = array(dim=c(length(t),3))
PO4.ans = array(dim = c(length(t),3))
Fborg.ans = array(dim = c(length(t),3))
Fbcarb = array(dim = c(length(t),3))
kcp.ans = array(dim=c(length(t),3))
fanox.ans = array(dim=c(length(t),3))
Fwsil.ans = array(dim=c(length(t),3))
Fwcarb.ans = array(dim=c(length(t),3))
omega.ans =array(dim=c(length(t),3))
Fbp.ans = array(dim=c(length(t),3))

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

fanox1 <-  vector()
fanox1[1] = fanox.i
fanox2 <- vector()
fanox2[1] = fanox.i
fanox3 <- vector()
fanox3[1] = fanox.i
fanox4 <- vector()
fanox4[1] = fanox.i
for (i in 2:100) {
  fanox1[i] = (0.2/100)*(i)
  fanox2[i] = (0.01/100)*(i)
  fanox3[i] = fanox.i
}

for (j in 1:3){
  for (i in 2:length(t)){
    dampen[i] = 1
    if (j==1){
      fanox[i] = dampen[i]*(0.2/100)*(Zub[i]+50)+fanox.i
    }else if (j==2){
      fanox[i] = dampen[i]*(0.01/100)*(Zub[i]+50)+fanox.i
    }else if (j==3){
      fanox[i] = fanox.i
    }
    Fborg[i] = Fborg.i*(Fbp[i-1]/Fbp[1])*(kcp[i-1]/kcp[1])
    Fwcarb[i] = 24e12*RCO2[i-1]^(G*Zc)*(1+G*Zc*log(RCO2[i-1]))
    Fbcarb[i] = kcarb*(omega[i-1]/omega.i)
    Fwp[i] = 1.3e11*RCO2[i-1]^(G*Z)*(1+G*Z*log(RCO2[i-1]))^0.65
    Fwsil[i] = 8e12*RCO2[i-1]^(G*Z)*(1+G*Z*log(RCO2[i-1]))^0.65
    Falkin[i] = 2*Fwcarb[i] + Fwsil[i]
    Falkout[i] = 2*Fbcarb[i]
    
    Alk[i] = Alk[i-1] + (Falkin[i]-Falkout[i])*dt
    DIC[i] = DIC[i-1] + (Fvolc+Fworg+Fwcarb[i]-Fborg[i]-Fbcarb[i])*dt
    PO4[i] = PO4[i-1] + (Fwp[i]-Fbp[i-1])*dt
    
    carb.parms = pCa(flag=15, var1=Alk[i]/oceanV, var2=DIC[i]/oceanV, Ca=Ca.i, S=salinity, T= temp.i, P=0, Pt=PO4[i]/oceanV)

    pCO2[i] = carb.parms$pCO2[2]
    pH.temp = carb.parms$pH[2]
    omega[i] = carb.parms$OmegaAragonite[2]
    RCO2[i] = pCO2[i]/pCO2.i
    
    kcp[i] = kcp[1]*(fanox[i]/fanox.i)^(1/3)*PO4[i]/PO4.i
    Fbp[i] = Fwp[1]*(PO4[i]/PO4.i)
    Iso.volc = Fvolc*(dCvolc- d13C[i-1])/DIC[i]
    Iso.wcarb = Fwcarb[i]*(dCwcarb - d13C[i-1])/DIC[i]
    Iso.worg = Fworg*(dCworg - d13C[i-1])/DIC[i]
    Iso.bcarb = Fbcarb[i]*frac.carb/DIC[i]
    Iso.borg = Fborg[i]*(frac.org)/DIC[i]
    d13C[i] = d13C[i-1] + (Iso.volc + Iso.wcarb + Iso.worg - Iso.bcarb- Iso.borg)*dt
    
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
    Fwcarb.ans[i,j] =Fwcarb[i]
    omega.ans[i,j] = omega[i]
    Fbp.ans[i,j] = Fbp[i]
  }
}

###plots###
par(mfrow=c(5,1), mar=c(2,4,1,1))

plot(fanox1*100, type='l', col='dark gray', xlab='Sea level(m)',ylab='fanox(%)',lwd=2)
lines(fanox3*100,col="dark green", lwd=2, lty='dotted')
lines(fanox2*100,col='blue',lwd=2, lty='dashed')
legend('topleft',c('shallow','Deep'), col=c('dark gray', 'blue'), lty=c('solid','dashed'),lwd=2, bty='n')

plot(t/1e6, Zub+50, type='l', xlab='Model time (Myr)', ylab='Relative change in OMZ_ub (m)', lwd=2)

plot(t/1e6, fanox.ans[,1]*100, type='l', col='dark gray', xlab='Model time (Myr)', ylab='fanox(%)', lwd=2)
lines(t/1e6, fanox.ans[,3]*100, col='dark green', lwd=2, lty='dotted')
lines(t/1e6, fanox.ans[,2]*100, col='blue', lwd=2, lty='dashed')

plot(t/1e6, Fborg.ans[,1], type='l', col='dark gray', xlab='Model time year (Myr)', ylab='Fborg (mol/yr)', lwd=2)
lines(t/1e6, Fborg.ans[,3], col='dark green', lwd=2, lty='dotted')
lines(t/1e6, Fborg.ans[,2], col='blue', lwd=2, lty='dashed')

plot(t/1e6, d13C.ans[,1], type='l', col='dark gray', xlab ='Model time year (Myr)', ylab='d13C', lwd=2)
lines(t/1e6, d13C.ans[,3], col='dark green', lwd=2, lty='dotted')
lines(t/1e6, d13C.ans[,2], col='blue', lwd=2, lty='dashed')




