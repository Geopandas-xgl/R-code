#Forward box model to explore alterative forcings on U systematics
#by Kimberly Lau

library(Hmisc)

#### U mass set up ####
sw.conc <- numeric() #ppb
sw.moles <- numeric() #in moles
sinks <- numeric()
carb.conc <- numeric() #in ppm
sw.conc[1] <- 3.3 #modern: 3.3 ppb, 1.96e13 mol U (Morford and Emerson 1999)#º£Ë®UÅ¨¶È
ocean.mass <- 1.41357575757576e21 #mass of seawater in g
M.riv <- 0.4e8
extent <- 0.0021 ##Set initial extent of anoxic sink (modern = 0.0021)
K.a <- 1.45772594752187e-4 #calculated
K.o <- 1.73834440079269e-6 #calculated to reach steady state
D.u <- 1.4 # partition coefficient factor to convert from moles cc to moles sw #1.6 original
sw.moles[1] <- sw.conc[1]*ocean.mass/(1e6*238) #coverts ppb to moles U

#### U isotope set up ####
d238U <- numeric() 
d238U[1] <- -0.165 #Late Permian average
dUriv <- -0.05 #calculated steady state input to produce d238U ~ -0.165
Frac.a <- 0.6 #Andersen et al. 2014
Frac.o <- 0.0294118

#### time set up ####
duration <- 2e6 #2 million years
dt <- 500 #500 year time steps
t<- seq(1, duration, dt)
carb1 <- numeric()
carb2 <- numeric()
carb3 <- numeric()
carb4 <- numeric()
d1 <- numeric()
d2 <- numeric()
d3 <- numeric()
d4 <- numeric()
dUriv1 <- numeric()
dUriv2 <- numeric()
dUriv3 <- numeric()
dUriv4 <- numeric()
Uriv1 <- numeric()
Uriv2 <- numeric()
Uriv3 <- numeric()
Uriv4 <- numeric()
fanox1 <- numeric()
fanox2 <- numeric()
fanox3 <- numeric()
fanox4 <- numeric()

#### model equations ####
for (j in 1:4) {
for (i in 2: length(t)) {
  if (i <= 718) {
    extent=0.0021
    M.riv=0.4e8
    dUriv=-0.05
  } else if (i <= 780) {
      if (j<=1) {
        extent = 0.2 } #rapid expansion of anoxia to 20%
      else if (j<=2) {
        M.riv <- 0.4e7 #decrease riverine flux and isotopic composition
        dUriv = -0.8 }
      else if (j<=3) {
        dUriv = -0.8 } #decrease riverine isotopic composition
      else if (j<=4) {
        M.riv <- 0.4e9 #increase riverine flux and decrease isotopic composition
        dUriv = -0.8
      }
  } else {
    if (j<=1) {
    extent=0.05 } #anoxia decreases to 5% for first scenario
  } 

  M.anox <- sw.moles[i-1]*K.a*extent
  M.other <- sw.moles[i-1]*K.o*(1-extent)
  sw.moles[i] = sw.moles[i-1] + (M.riv - M.anox - M.other)*dt
  sw.conc[i] = sw.moles[i]*238*1e6/ocean.mass
  sinks[i] = M.anox + M.other
  carb.conc[i] <- sw.conc[i]/D.u
  
  Iso.riv <- M.riv*(dUriv - d238U[i-1])/sw.moles[i-1]
  Iso.other <- M.other*(Frac.o)/sw.moles[i-1]
  Iso.anox <- M.anox*(Frac.a)/sw.moles[i-1]
  d238U[i] = d238U[i-1] + (Iso.riv - Iso.other - Iso.anox)*dt 
  
  if (j<=1) {
    carb1[i] = carb.conc[i]
    d1[i] = d238U[i] 
    dUriv1[i] = dUriv
    Uriv1[i] = M.riv
    fanox1[i] = extent
  } else if (j<=2) {
    carb2[i] = carb.conc[i]
    d2[i] = d238U[i]
    dUriv2[i] = dUriv
    Uriv2[i] = M.riv
    fanox2[i] = extent
  } else if (j<=3) {
    carb3[i] = carb.conc[i]
    d3[i] = d238U[i]
    dUriv3[i] = dUriv
    Uriv3[i] = M.riv
    fanox3[i] = extent
  } else if (j<=4) {
    carb4[i] = carb.conc[i]
    d4[i] = d238U[i]
    dUriv4[i] = dUriv
    Uriv4[i] = M.riv
    fanox4[i] = extent
  }
}}

#### plots ####
#stacked figure
plot(carb1,t,ylab="Time (yr)",xlab="[U] (ppm)",cex=0.5,xaxt="n",xlim=c(0.05,30),log="x",type="l",col="black")
lines(carb2,t,lty="longdash", cex=10)
lines(carb3,t,lty="dotdash", cex=10)
lines(carb4,t,lty="dotted", cex=10)
axis(1,c(0.1,1,10))

###plot(d1,t,ylab="Time (yr)", xlab="d238U",xlim=c(-1.4,0.4),cex=0.5,type="l",col="black")
lines(d2,t,lty="longdash")
lines(d3,t,lty="dotdash")
lines(d4,t,lty="dotted")
minor.tick(nx=5)###
