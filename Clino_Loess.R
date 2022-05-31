rm(list = ls())

library(xlsx)

data_1 = read.xlsx("C:\\Users\\19328\\Desktop\\Xiong-2022-U isotope\\图件\\大气水带\\U isotope_Clino_NK_XK.xlsx",1) 
data_2 = read.xlsx("C:\\Users\\19328\\Desktop\\Xiong-2022-U isotope\\图件\\大气水带\\U isotope_Clino_NK_XK.xlsx",2) 
data_3 = read.xlsx("C:\\Users\\19328\\Desktop\\Xiong-2022-U isotope\\图件\\大气水带\\U isotope_Clino_NK_XK.xlsx",3)

model1 = loess(delta_U ~ depth, data=data_1, span=0.6, degree=2)
model2 = loess(delta_U ~ depth, data=data_2, span=0.9, degree=2)
model3 = loess(delta_U ~ depth, data=data_3, span=0.8, degree=2)

smoothed1 = predict(model1)
smoothed2 = predict(model2)
smoothed3 = predict(model3)

print(smoothed2)
print(smoothed3)

plot(data_1$delta_U, data_1$depth)
lines(smoothed1, y=data_1$depth, col = 'red')

plot(data_2$delta_U, data_2$depth)
lines(smoothed2, y=data_2$depth, col = 'red')

plot(data_3$delta_U, data_3$depth)
lines(smoothed3, y=data_3$depth, col = 'red')
