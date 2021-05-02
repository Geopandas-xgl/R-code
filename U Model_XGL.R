libary(Hmisc)
###U mass set up###
sw.conc <- numeric() #ppb
sw.moles <- numeric() #in moles
sinks <- numeric()
carb.conc <- numeric() #in ppm
sw.conc[1] <- 3.3 #modern:3.3 ppb
ocean.mass <- 1.41357575757576 # mass of seawater in g
M.riv <- 0.4e8
extent <- 0.0021
K.a <- 1.45772594752187e-4 #calculated
K.o <- 1.73834440079269e-6 #calulated to reach steady state
D.u <- 1.4
sw.moles[1] <-sw.conc[1]*ocean.mass/(1e6*238)

###U isotope set up ###
d238U <- numeric()
d238U[1] <- -0.165 #Late Permian average
dUriv <- -0.05
Frac.a <- 0.6 
Frac.o <- 0.0294118

###time set up###
duration <- 2e6 #2 million years
dt <- 500 #500 year time steps
t <- seq(1,duration,dt)
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

### model equations ###
for(j in 1:4)
{
  for(i in 2:length(t))
{
    if(i <= 718){
      extent = 0.0021
      M.riv = 0.4e8
      dUriv = -0.05}
    else if(i <= 780)
    {
      if(j <= 1){
         extent = 0.2}
      else if(j <= 2){
         M.riv <- 0.4e7
         dUriv = -0.8}
      else if(j <=3){
         dUriv = -0.8}
      else if(j<=4){
         M.riv <- 0.4e9
         dUriv = -0.8}
    }
    else
    {
      if(j<=1){
      extent = 0.05}
    }
    M.anox <- sw.moles[i-1]*K.a*extent
       

























