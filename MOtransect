
#------------------The MO transect------------------------#
library(dplyr)
library(tidyr)
library(tibble)


cline <- read.csv("/Users/sonjalecic/Downloads/t2-mo-percentages_full.csv")
cline[complete.cases(cline),]
cline_t <- as_tibble(cline)
# extract T2 only
mo <- cline[c(51:82),]
mo <- as_tibble(mo)

mo.long <- mo %>% as_tibble() %>% pivot_longer(cols = c("X2015", "X2016", "X2017", "X2018", "X2019", "X2020"), names_to ="year", values_to = "infection")

# add an imaginary counts as if the sample size for all locs/years was 25
mo.long$total <- 25
mo.long$inf <- (mo.long$infection*mo.long$total)/100
mo.long$uninf <- mo.long$total - mo.long$inf
#t2.long[is.na(t2.long)] = 0
mo.long[!is.na(mo.long$inf),]
mo.long.na <- na.omit(mo.long)
mo.long.na$infectionFreq <- mo.long.na$inf/mo.long.na$total


# remove locations with NA for lon and/or lat
mo.long.na.na <-subset(mo.long.na, !is.na(latitude) & !is.na(longitude))
dist2mo <- c()
for(i in 1:nrow(mo.long.na.na)) {
  
  row <- mo.long.na.na[i,] # extract row
  
  # approximate radius of earth in km
  R = 6373.0
  lat1 = rad(48.69315)
  lon1 = rad(16.24286)
  lat2 = rad(row$latitude)
  lon2 = rad(row$longitude)
  
  dlon = lon2 - lon1
  dlat = lat2 - lat1
  
  a = sin(dlat / 2)**2 + cos(lat1) * cos(lat2) * sin(dlon / 2)**2
  c = 2 * atan2(sqrt(a), sqrt(1 - a))
  # distance in meters
  distance = (R * c)*1000
  #print(distance)
  dist2mo  <- c(dist2mo, distance) # save as a vector
}
dist2mo

# add distance as a column
mo.long.na.na$dist <- dist2mo/1000

# transfrom years into generation from 1 to 6
mo.long.na.na.t <- mo.long.na.na %>%
  mutate(time = case_when(
    year == "X2015" ~ 1,
    year == "X2016" ~ 2,
    year == "X2017" ~ 3,
    year == "X2018" ~ 4,
    year == "X2019" ~ 5,
    year == "X2020" ~ 6
    
  ))

moData = mo.long.na.na.t[,c(10, 11, 12)]


#-----------------nls difusion model fitting-----------------------#

# the diffusion model
f <- function(x, phat, t){1/(1+(exp(x-((1-(2*phat))*t))))}

#our infection frequencies
y=mo.long.na.na.t$infectionFreq
#distance from location 1 (calculated from online distance calculaors based on coordinates)
x=mo.long.na.na.t$dist
# time
t=mo.long.na.na.t$time

#curve fitting

plot(x,y)
#curve(f(x,5,7), add = TRUE) #test the function plot with arbitrary values
curve(f(x, -2, 5), add = TRUE)

#sigma.t.fit <- nls(y~f(x,sigma,t), start = list(sigma=1,t=1))
#library(broom)
moTidy <- tibble()
for(i in 1:nrow(mo.long.na.na.t)) {
  row = mo.long.na.na.t[i,]
  phat.mo.fit <- nls(y~f(x, phat, row$time), start = list(phat=0)
                     #,lower=c(0), algorithm = "port"
  )
  #print(sigma.t.fit) 
  broomed <- as.data.frame(cbind(broom::tidy(phat.mo.fit), row$time))
  moTidy <- rbind(moTidy, broomed)
}

phat.mo.fit
summary(phat.mo.fit)
subset(moTidy, row$time==2)
moTidy[which(moTidy$`row$time` == 2),]$estimate[1]

# plot the model fitted curves for each time
tiff("/Volumes/LaCie/the_cline/MOnlsFit.tiff", width = 8, height = 6, units = "in", res = 300) #this is for saving higher quality plots only
#plotCI(x,y, main="MO transect", ylab="Infection Frequency", xlab="Distance(km)", err = "y", uiw = CIh, liw = CIl, col = "grey", xlim = c(0,80), slty = 2, cex.axis=1.5,cex.lab=1.5, cex.main=1.5)
plot(x,y, lwd = "1", ylab="Infection Frequency", xlab="Distance(km)", frame.plot = FALSE, ylim = c(0,1), pch = 19,  col= "grey60")
curve(f(x,moTidy[which(moTidy$`row$time` == 2),]$estimate[1],t=2), add = TRUE, lwd =2, col = "darkred", lty = 2 )
curve(f(x,moTidy[which(moTidy$`row$time` == 3),]$estimate[1],t=3), add = TRUE, lwd =2, col = "darkred", lty = 2 )
curve(f(x,moTidy[which(moTidy$`row$time` == 4),]$estimate[1],t=4), add = TRUE, lwd =2, col = "darkred", lty = 2 )
curve(f(x,moTidy[which(moTidy$`row$time` == 5),]$estimate[1],t=5), add = TRUE, lwd =2, col = "darkred", lty = 2 )
text(x= 15, y= 0.9, paste("t1_phat=", round(tidy[which(tidy$`row$time` == 1),]$estimate[1],2)), cex=1)
text(x= 15, y= 0.85, paste("t2_phat=", round(tidy[which(tidy$`row$time` == 2),]$estimate[1],2)), cex=1)
text(x= 15, y= 0.8, paste("t3_phat=", round(tidy[which(tidy$`row$time` == 3),]$estimate[1],2)), cex=1)
text(x= 15, y= 0.75, paste("t4_phat=", round(tidy[which(tidy$`row$time` == 4),]$estimate[1],2)), cex=1)
text(x= 15, y= 0.7, paste("t5_phat=", round(tidy[which(tidy$`row$time` == 5),]$estimate[1],2)), cex=1)
text(x= 15, y= 0.65, paste("t6_phat=", round(tidy[which(tidy$`row$time` == 6),]$estimate[1],2)), cex=1)
text(x= 15, y= 0.6, paste("sigma=", round(broom::tidy(sigma.t.fit)$estimate[1],3)), cex=1.1)
dev.off()

summary(phat.mo.fit)

#residual sum of squares and R value
RSS <- sum(residuals(phat.mo.fit)^2)
TSS <- sum((y - mean(y))^2)
R.square <- 1 - (RSS/TSS)
R.square


#-----------------Maximum Likelihood difusion model fitting-----------------------#

# make a function for the diffusion approximation formula
dtfun <- function(dis, phat, t){
  # Prediction of the model
  1/(1+(exp(dis-((1-(2*phat))*t))))
}

# plot quickly to check
plot(x=mo.long.na.na.t$dist, y=mo.long.na.na.t$infectionFreq, xlab = "Distance(km)", ylab = "Infection frequency", ylim=c(0,1))
curve(dtfun(dis=x, phat=0, t=25), mo.long.na.na.t$dist, add = TRUE)

# extract only dist, frequncy & time
modata <- mo.long.na.na.t[,c(10, 11,12)]
options(scipen = 999)
modata <- as.data.frame(modata)

NLL = function(phat, t=5, data) {
  # Values predicted by the model
  ppred = dtfun(modata$dist, phat, t=5)
  #print(ppred)
  # binomial distibution
  -sum(dbinom(x = ppred, 25, modata$infectionFreq), log=T)
  #dbinom(x = data$infectionFreq, 92, ppred)
  #-sum(dnorm(x = data$G, mean = Gpred, sd = pars[4], log = TRUE))
}
#data3 <- data[,c(2,1,3)]
#fit = optim(par = par0, fn = NLL, data = data, hessian = T)
fitmo = optim(par = c(-1.75), data=modata, NLL, lower=c(-Inf), upper=c(Inf))
fitmo$par
optimize(NLL, lower=c(-Inf), upper=c(Inf))


### Maximum likelihood fitting function ####
tidyMLmo <- data.frame()
for(i in 1:nrow(modata)){
  row = modata[i,]
  NLL = function(phat, t=row$time, modata) {
    # Values predicted by the model
    ppred = dtfun(modata$dist, phat, t=row$time)
    #print(ppred)
    # binomial distibution
    -sum(dbinom(x = ppred, 25, modata$infectionFreq), log=T)
  }
  if (row$time == 2) {
    fit = optim(par = c(-5.17), NLL, modata=modata ,lower=c(-Inf), upper=c(Inf)) 
  } else if (row$time == 3) {
    fit = optim(par = c(-3.27), NLL, modata=modata ,lower=c(-Inf), upper=c(Inf)) 
  } else if (row$time == 4) {
    fit = optim(par = c(-2.33), NLL, modata=modata ,lower=c(-Inf), upper=c(Inf)) 
  } else if (row$time == 5) {
    fit = optim(par = c(-1.75), NLL, modata=modata ,lower=c(-Inf), upper=c(Inf)) 
  }
  #fit = optim(par = c(1), NLL, data=data ,lower=c(-Inf), upper=c(Inf))
  #fit$par
  broomedMl <- as.data.frame(cbind(broom::tidy(fit), row$time))
  tidyMLmo <- rbind(tidyMLmo, broomedMl)
}


tiff("/Volumes/LaCie/the_cline/MoMLFit.tiff", width = 8, height = 6, units = "in", res = 300)
plot(x=modata$dist, y=modata$infectionFreq, xlab = "Distance(km)", ylab = "Infection frequency", frame.plot = FALSE, ylim = c(0,1), pch = 19,  col= "grey60")
curve(dtfun(x, tidyMLmo[which(tidyMLmo$`row$time` == 2),]$value[1], t=2), data$dist, add = TRUE, lwd =2, col = "darkred", lty = 2)
curve(dtfun(x, tidyMLmo[which(tidyMLmo$`row$time` == 3),]$value[1], t=3), data$dist, add = TRUE, lwd =2, col = "red3", lty = 2)
curve(dtfun(x, tidyMLmo[which(tidyMLmo$`row$time` == 4),]$value[1], t=4), data$dist, add = TRUE, lwd =2, col = "red3", lty = 2)
curve(dtfun(x, tidyMLmo[which(tidyMLmo$`row$time` == 5),]$value[1], t=5), data$dist, add = TRUE, lwd =2, col = "red3", lty = 2)
dev.off()


#---------wave speed MO stransect for year 2016-------------#
fsigma <- function(x, sigma, t){0.5 * (1 - tanh((x - sigma * (1 - 2 * 0) * sqrt(0.98) * t / 2)/(2 * sigma / sqrt(0.98))))}

#subset year 2016
mo19 <- subset(mo.long.na.na.t, time == 5)

# infection frequencies for year 2016
y=mo19$infectionFreq
#distance from location 1 (calculated from online distance calculaors based on coordinates)
x=mo19$dist

plot(x, y)
#find a best fitting sigma and t
sigma.mo.fit <- nls(y~fsigma(x,sigma,t), start = list(sigma=1,t=5))

# wave speed formula = σ*√lCI)/2
moSpeed <- (9.272*sqrt(0.98))/2

#residual sum of squares and R value
RSS <- sum(residuals(sigma.mo.fit)^2)
TSS <- sum((y - mean(y))^2)
Rnls.MO <- 1 - (RSS/TSS)
Rnls.MO


tiff("/Volumes/LaCie/the_cline/MO_sigma2019.tiff", width = 8, height = 6, units = "in", res = 300) #this is for saving higher quality plots only
#plotCI(x,y, main="Transition Zone CZ", ylab="Infection Frequency", xlab="Distance(km)", err = "y", uiw = CIh, liw = CIl, col = "grey", xlim = c(0,80), slty = 2, cex.axis=1.5,cex.lab=1.5, cex.main=1.5)
plot(x,y, lwd = "1", ylab="Infection Frequency", xlab="Distance(km)", frame.plot = FALSE, ylim = c(0,1), pch = 19,  col= "grey60")
curve(fsigma(x,broom::tidy(sigma.mo.fit)$estimate[1],t=broom::tidy(sigma.mo.fit)$estimate[2]), add = TRUE, lwd =2, col = "blue", lty = 2 )
text(x= 60, y= 0.9, paste("sigma=", round(broom::tidy(sigma.mo.fit)$estimate[1],3)), cex=1.1)
text(x= 60, y= 0.8, paste("wave speed=", round(moSpeed,3)), cex=1.1)
text(x= 60, y= 0.7, paste("R=", round(Rnls.MO,3)), cex=1.1)
dev.off()


----------#########--------------
moData

X = sqrt((2*1*(3.65554012791^2))/(9.272^2))

moData$scaledX <- sqrt((2*1*(moData$dist^2))/(9.272^2))

tidyMLmo <- data.frame()
for(i in 1:nrow(moData)){
  row = moData[i,]
  NLL = function(phat, t=row$time, moData) {
    # Values predicted by the model
    ppred = dtfun(moData$scaledX, phat, t=row$time)
    #print(ppred)
    # binomial distibution
    -sum(dbinom(x = ppred, 25, moData$infectionFreq), log=T)
  }
  if (row$time == 2) {
    fit = optim(par = c(-5.17), NLL, moData=moData ,lower=c(-Inf), upper=c(Inf)) 
  } else if (row$time == 3) {
    fit = optim(par = c(-3.27), NLL, moData=moData ,lower=c(-Inf), upper=c(Inf)) 
  } else if (row$time == 4) {
    fit = optim(par = c(-2.33), NLL, moData=moData ,lower=c(-Inf), upper=c(Inf)) 
  } else if (row$time == 5) {
    fit = optim(par = c(-1.75), NLL, moData=moData ,lower=c(-Inf), upper=c(Inf)) 
  }
  #fit = optim(par = c(1), NLL, data=data ,lower=c(-Inf), upper=c(Inf))
  #fit$par
  broomedMl <- as.data.frame(cbind(broom::tidy(fit), row$time))
  tidyMLmo <- rbind(tidyMLmo, broomedMl)
}


y=moData$infectionFreq
#distance from location 1 (calculated from online distance calculaors based on coordinates)
x=moData$scaledX
# time
t=moData$time

#curve fitting

plot(x,y)
#curve(f(x,5,7), add = TRUE) #test the function plot with arbitrary values
curve(f(x, -2, 5), add = TRUE)

#find a best fitting sigma and t
#sigma.t.fit <- nls(y~f(x,sigma,t), start = list(sigma=1,t=1))
library(broom)
moTidy <- tibble()
for(i in 1:nrow(moData)) {
  row = moData[i,]
  phat.mo.fit <- nls(y~f(x, phat, row$time), start = list(phat=0)
                     #,lower=c(0), algorithm = "port"
  )
  #print(sigma.t.fit) 
  broomed <- as.data.frame(cbind(broom::tidy(phat.mo.fit), row$time))
  moTidy <- rbind(moTidy, broomed)
}


tiff("/Volumes/LaCie/the_cline/MoNLSFit_scaledX.tiff", width = 8, height = 6, units = "in", res = 300)
plot(x=moData$scaledX, y=moData$infectionFreq, xlab = "Distance(km)", ylab = "Infection frequency", frame.plot = FALSE, ylim = c(0,1), xlim = c(0, 12), pch = 19,  col= "grey60")
curve(f(x,moTidy[which(moTidy$`row$time` == 2),]$estimate[1],t=2), add = TRUE, lwd =2, col = "darkred", lty = 2 )
curve(f(x,moTidy[which(moTidy$`row$time` == 3),]$estimate[1],t=3), add = TRUE, lwd =2, col = "darkred", lty = 2 )
curve(f(x,moTidy[which(moTidy$`row$time` == 4),]$estimate[1],t=4), add = TRUE, lwd =2, col = "darkred", lty = 2 )
curve(f(x,moTidy[which(moTidy$`row$time` == 5),]$estimate[1],t=5), add = TRUE, lwd =2, col = "darkred", lty = 2 )
#text(x= 15, y= 0.9, paste("t1_phat=", round(tidy[which(tidy$`row$time` == 1),]$estimate[1],2)), cex=1)
text(x= 9, y= 0.8, paste("t2_phat=", round(moTidy[which(moTidy$`row$time` == 2),]$estimate[1],2)), cex=1)
text(x= 9, y= 0.75, paste("t3_phat=", round(moTidy[which(moTidy$`row$time` == 3),]$estimate[1],2)), cex=1)
text(x= 9, y= 0.7, paste("t4_phat=", round(moTidy[which(moTidy$`row$time` == 4),]$estimate[1],2)), cex=1)
text(x= 9, y= 0.65, paste("t5_phat=", round(moTidy[which(moTidy$`row$time` == 5),]$estimate[1],2)), cex=1)
#text(x= 15, y= 0.65, paste("t6_phat=", round(tidy[which(tidy$`row$time` == 6),]$estimate[1],2)), cex=1)
text(x= 9, y= 0.6, paste("sigma=", round(broom::tidy(sigma.mo.fit)$estimate[1],3)), cex=1)
dev.off()


##----------------two-dimensional space---------------##
f <- function(x, phat, t, k){1/(2*(0.5 - phat))*(1 + log(x) * (exp(((2*(0.5 - phat))^2)*(t + k)-1)))}

nlc <- nls.control(maxiter= 100, tol=1e-02, warnOnly=TRUE)
moTidy <- tibble()
for(i in 1:nrow(moData)) {
  row = moData[i,]
  phat.mo.fit <- nls(y~f(x, phat, row$time, k=0.4), start = list(phat=-1),
                     #,lower=c(0), algorithm = "port"
  control = nlc)
  #print(sigma.t.fit) 
  broomed <- as.data.frame(cbind(broom::tidy(phat.mo.fit), row$time))
  moTidy <- rbind(moTidy, broomed)
}

#tiff("/Volumes/LaCie/the_cline/MoMLFit.tiff", width = 8, height = 6, units = "in", res = 300)
plot(x=moData$scaledX, y=moData$infectionFreq, xlab = "Distance(km)", ylab = "Infection frequency", frame.plot = FALSE, ylim = c(0,1), xlim = c(0,15), pch = 19,  col= "grey60")
curve(dtfun(x, 0.2988578, t=15), moData$scaledX, add = TRUE, lwd =2, col = "darkred", lty = 2)
curve(dtfun(x, tidyMLmo[which(tidyMLmo$`row$time` == 3),]$value[1], t=3), moData$dist, add = TRUE, lwd =2, col = "red3", lty = 2)
curve(dtfun(x, tidyMLmo[which(tidyMLmo$`row$time` == 4),]$value[1], t=4), moData$dist, add = TRUE, lwd =2, col = "red3", lty = 2)
curve(dtfun(x, moTidy[which(moTidy$`row$time` == 5),]$estimate[1], t=5), modata$scaledX, add = TRUE, lwd =2, col = "red3", lty = 2)
#dev.off()
