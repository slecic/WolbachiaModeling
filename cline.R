
#--- libraries
library(dplyr)
library(tidyr)
library(tibble)
library(circular) # calculate distance from x,y coordinates
library(broom) # for tidying the data

## cline data ####

cline <- read.csv("/Users/sonjalecic/Downloads/t2-mo-percentages_full.csv")
cline[complete.cases(cline),]
cline_t <- as_tibble(cline)
# extract T2 only
t2 <- cline_t[c(1:50),]
tail(t2)

t2.long <- t2 %>% as_tibble() %>% pivot_longer(cols = c("X2015", "X2016", "X2017", "X2018", "X2019", "X2020"), names_to ="year", values_to = "infection")

#---------calculate summary data for each individual trapping site using weighted means---------------------
# ---year 2015-----
y15 <- subset(t2.long.na, year == "X2015")

trapData <- tibble()
splitData <- split(y15, y15$pop)

for(i in 1:length(splitData)){
  trapData[i,1] <- splitData[[i]]$id[1] #id
  trapData[i,2] <- splitData[[i]]$pop[1] #pop
  trapData[i,3] <- splitData[[i]]$latitude[1] # lat
  trapData[i,4] <- splitData[[i]]$longitude[1] #lon
  trapData[i,5] <- splitData[[i]]$year[1] # year
  trapData[i,6] <- splitData[[i]]$infection[1] # infection
  trapData[i,7] <- splitData[[i]]$inf[1] # number of infected individuals
  trapData[i,8] <- splitData[[i]]$uninf[1] # number of uninfected individuals
  #trapData[i,5] <- length(which(splitData[[i]]$Infected == "I" )) #n infected families
  #trapData[i,6] <- length(which(splitData[[i]]$Infected == "U" )) #n uninfected families
  
  #Estimate infection freq. and 95% binomial confidence intervals
  trapData[i,9] <- binconf(x=as.numeric(trapData[i,7]), n=as.numeric(trapData[i,7] + trapData[i,8]))[1] #infection freq.
  trapData[i,10] <- binconf(x=as.numeric(trapData[i,7]), n=as.numeric(trapData[i,7] + trapData[i,8]))[2]
  trapData[i,11] <- binconf(x=as.numeric(trapData[i,7]), n=as.numeric(trapData[i,7] + trapData[i,8]))[3]

} 
rm(i)
rm(splitData)

colnames(trapData) <- c("id", "pop", "lat", "lon", 
                        "year", "infection", "nInf", "nUnInf", "infectionFreq", "infectionFreq_lowerCI", "infectionFreq_upperCI")

# round the decimals
trapData[,9:11] <- round(trapData[,9:11], digits=3)

# re-order the data frame based on spatial location of the cline
trapData.ord <- trapData[match(y15$pop, trapData$pop),]

#plot yakuba infection freq. data
p1 <- ggplot(trapData.ord, aes(x=factor(pop, level = pop), y=infectionFreq)) + 
  geom_errorbar(aes(ymin=infectionFreq_lowerCI, 
                    ymax=infectionFreq_upperCI), 
                width=.2, size=1.5) +
  geom_point(size=4) +
  theme_bw() +
  theme(axis.title.y = element_text(vjust= 1.8),
        axis.title.x = element_text(vjust= 0.5),
        axis.text.x = element_text(angle = 90, vjust= 0.5),
        axis.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(limits = c(0, 1.10), breaks=c(0,0.25,0.5,0.75,1.0)) +
  labs(y = "Infection Frequency", x = "Site") +
  ggtitle("Year 2015")

p1


# ---year 2016-----
y16 <- subset(t2.long.na, year == "X2016")

trapData <- tibble()
splitData <- split(y16, y16$pop)

for(i in 1:length(splitData)){
  trapData[i,1] <- splitData[[i]]$id[1] #id
  trapData[i,2] <- splitData[[i]]$pop[1] #pop
  trapData[i,3] <- splitData[[i]]$latitude[1] # lat
  trapData[i,4] <- splitData[[i]]$longitude[1] #lon
  trapData[i,5] <- splitData[[i]]$year[1] # year
  trapData[i,6] <- splitData[[i]]$infection[1] # infection
  trapData[i,7] <- splitData[[i]]$inf[1] # number of infected individuals
  trapData[i,8] <- splitData[[i]]$uninf[1] # number of uninfected individuals
  #trapData[i,5] <- length(which(splitData[[i]]$Infected == "I" )) #n infected families
  #trapData[i,6] <- length(which(splitData[[i]]$Infected == "U" )) #n uninfected families
  
  #Estimate infection freq. and 95% binomial confidence intervals
  trapData[i,9] <- binconf(x=as.numeric(trapData[i,7]), n=as.numeric(trapData[i,7] + trapData[i,8]))[1] #infection freq.
  trapData[i,10] <- binconf(x=as.numeric(trapData[i,7]), n=as.numeric(trapData[i,7] + trapData[i,8]))[2]
  trapData[i,11] <- binconf(x=as.numeric(trapData[i,7]), n=as.numeric(trapData[i,7] + trapData[i,8]))[3]
  
} 
rm(i)
rm(splitData)

colnames(trapData) <- c("id", "pop", "lat", "lon", 
                        "year", "infection", "nInf", "nUnInf", "infectionFreq", "infectionFreq_lowerCI", "infectionFreq_upperCI")

# round the decimals
trapData[,9:11] <- round(trapData[,9:11], digits=3)

# re-order the data frame based on spatial location of the cline
trapData.ord <- trapData[match(y16$pop, trapData$pop),]

#plot yakuba infection freq. data
p2 <- ggplot(trapData.ord, aes(x=factor(pop, level = pop), y=infectionFreq)) + 
  geom_errorbar(aes(ymin=infectionFreq_lowerCI, 
                    ymax=infectionFreq_upperCI), 
                width=.2, size=1.5) +
  geom_point(size=4) +
  theme_bw() +
  theme(axis.title.y = element_text(vjust= 1.8),
        axis.title.x = element_text(vjust= 0.5),
        axis.text.x = element_text(angle = 90, vjust= 0.5),
        axis.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(limits = c(0, 1.10), breaks=c(0,0.25,0.5,0.75,1.0)) +
  labs(y = "Infection Frequency", x = "Site") +
  ggtitle("Year 2016")

p2


# ---year 2017-----
y17 <- subset(t2.long.na, year == "X2017")

trapData <- tibble()
splitData <- split(y17, y17$pop)

for(i in 1:length(splitData)){
  trapData[i,1] <- splitData[[i]]$id[1] #id
  trapData[i,2] <- splitData[[i]]$pop[1] #pop
  trapData[i,3] <- splitData[[i]]$latitude[1] # lat
  trapData[i,4] <- splitData[[i]]$longitude[1] #lon
  trapData[i,5] <- splitData[[i]]$year[1] # year
  trapData[i,6] <- splitData[[i]]$infection[1] # infection
  trapData[i,7] <- splitData[[i]]$inf[1] # number of infected individuals
  trapData[i,8] <- splitData[[i]]$uninf[1] # number of uninfected individuals
  #trapData[i,5] <- length(which(splitData[[i]]$Infected == "I" )) #n infected families
  #trapData[i,6] <- length(which(splitData[[i]]$Infected == "U" )) #n uninfected families
  
  #Estimate infection freq. and 95% binomial confidence intervals
  trapData[i,9] <- binconf(x=as.numeric(trapData[i,7]), n=as.numeric(trapData[i,7] + trapData[i,8]))[1] #infection freq.
  trapData[i,10] <- binconf(x=as.numeric(trapData[i,7]), n=as.numeric(trapData[i,7] + trapData[i,8]))[2]
  trapData[i,11] <- binconf(x=as.numeric(trapData[i,7]), n=as.numeric(trapData[i,7] + trapData[i,8]))[3]
  
} 
rm(i)
rm(splitData)

colnames(trapData) <- c("id", "pop", "lat", "lon", 
                        "year", "infection", "nInf", "nUnInf", "infectionFreq", "infectionFreq_lowerCI", "infectionFreq_upperCI")

# round the decimals
trapData[,9:11] <- round(trapData[,9:11], digits=3)

# re-order the data frame based on spatial location of the cline
trapData.ord <- trapData[match(y17$pop, trapData$pop),]

#plot yakuba infection freq. data
p3 <- ggplot(trapData.ord, aes(x=factor(pop, level = pop), y=infectionFreq)) + 
  geom_errorbar(aes(ymin=infectionFreq_lowerCI, 
                    ymax=infectionFreq_upperCI), 
                width=.2, size=1.5) +
  geom_point(size=4) +
  theme_bw() +
  theme(axis.title.y = element_text(vjust= 1.8),
        axis.title.x = element_text(vjust= 0.5),
        axis.text.x = element_text(angle = 90, vjust= 0.5),
        axis.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(limits = c(0, 1.10), breaks=c(0,0.25,0.5,0.75,1.0)) +
  labs(y = "Infection Frequency", x = "Site") +
  ggtitle("Year 2017")

p3


# ---year 2018-----
y18 <- subset(t2.long.na, year == "X2018")

trapData <- tibble()
splitData <- split(y18, y18$pop)

for(i in 1:length(splitData)){
  trapData[i,1] <- splitData[[i]]$id[1] #id
  trapData[i,2] <- splitData[[i]]$pop[1] #pop
  trapData[i,3] <- splitData[[i]]$latitude[1] # lat
  trapData[i,4] <- splitData[[i]]$longitude[1] #lon
  trapData[i,5] <- splitData[[i]]$year[1] # year
  trapData[i,6] <- splitData[[i]]$infection[1] # infection
  trapData[i,7] <- splitData[[i]]$inf[1] # number of infected individuals
  trapData[i,8] <- splitData[[i]]$uninf[1] # number of uninfected individuals
  #trapData[i,5] <- length(which(splitData[[i]]$Infected == "I" )) #n infected families
  #trapData[i,6] <- length(which(splitData[[i]]$Infected == "U" )) #n uninfected families
  
  #Estimate infection freq. and 95% binomial confidence intervals
  trapData[i,9] <- binconf(x=as.numeric(trapData[i,7]), n=as.numeric(trapData[i,7] + trapData[i,8]))[1] #infection freq.
  trapData[i,10] <- binconf(x=as.numeric(trapData[i,7]), n=as.numeric(trapData[i,7] + trapData[i,8]))[2]
  trapData[i,11] <- binconf(x=as.numeric(trapData[i,7]), n=as.numeric(trapData[i,7] + trapData[i,8]))[3]
  
} 
rm(i)
rm(splitData)

colnames(trapData) <- c("id", "pop", "lat", "lon", 
                        "year", "infection", "nInf", "nUnInf", "infectionFreq", "infectionFreq_lowerCI", "infectionFreq_upperCI")

# round the decimals
trapData[,9:11] <- round(trapData[,9:11], digits=3)

# re-order the data frame based on spatial location of the cline
trapData.ord <- trapData[match(y18$pop, trapData$pop),]

#plot yakuba infection freq. data
p4 <- ggplot(trapData.ord, aes(x=factor(pop, level = pop), y=infectionFreq)) + 
  geom_errorbar(aes(ymin=infectionFreq_lowerCI, 
                    ymax=infectionFreq_upperCI), 
                width=.2, size=1.5) +
  geom_point(size=4) +
  theme_bw() +
  theme(axis.title.y = element_text(vjust= 1.8),
        axis.title.x = element_text(vjust= 0.5),
        axis.text.x = element_text(angle = 90, vjust= 0.5),
        axis.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(limits = c(0, 1.10), breaks=c(0,0.25,0.5,0.75,1.0)) +
  labs(y = "Infection Frequency", x = "Site") +
  ggtitle("Year 2018")

p4


# ---year 2019-----
y19 <- subset(t2.long.na, year == "X2019")

trapData <- tibble()
splitData <- split(y19, y19$pop)

for(i in 1:length(splitData)){
  trapData[i,1] <- splitData[[i]]$id[1] #id
  trapData[i,2] <- splitData[[i]]$pop[1] #pop
  trapData[i,3] <- splitData[[i]]$latitude[1] # lat
  trapData[i,4] <- splitData[[i]]$longitude[1] #lon
  trapData[i,5] <- splitData[[i]]$year[1] # year
  trapData[i,6] <- splitData[[i]]$infection[1] # infection
  trapData[i,7] <- splitData[[i]]$inf[1] # number of infected individuals
  trapData[i,8] <- splitData[[i]]$uninf[1] # number of uninfected individuals
  #trapData[i,5] <- length(which(splitData[[i]]$Infected == "I" )) #n infected families
  #trapData[i,6] <- length(which(splitData[[i]]$Infected == "U" )) #n uninfected families
  
  #Estimate infection freq. and 95% binomial confidence intervals
  trapData[i,9] <- binconf(x=as.numeric(trapData[i,7]), n=as.numeric(trapData[i,7] + trapData[i,8]))[1] #infection freq.
  trapData[i,10] <- binconf(x=as.numeric(trapData[i,7]), n=as.numeric(trapData[i,7] + trapData[i,8]))[2]
  trapData[i,11] <- binconf(x=as.numeric(trapData[i,7]), n=as.numeric(trapData[i,7] + trapData[i,8]))[3]
  
} 
rm(i)
rm(splitData)

colnames(trapData) <- c("id", "pop", "lat", "lon", 
                        "year", "infection", "nInf", "nUnInf", "infectionFreq", "infectionFreq_lowerCI", "infectionFreq_upperCI")

# round the decimals
trapData[,9:11] <- round(trapData[,9:11], digits=3)

# re-order the data frame based on spatial location of the cline
trapData.ord <- trapData[match(y19$pop, trapData$pop),]

#plot yakuba infection freq. data
p5 <- ggplot(trapData.ord, aes(x=factor(pop, level = pop), y=infectionFreq)) + 
  geom_errorbar(aes(ymin=infectionFreq_lowerCI, 
                    ymax=infectionFreq_upperCI), 
                width=.2, size=1.5) +
  geom_point(size=4) +
  theme_bw() +
  theme(axis.title.y = element_text(vjust= 1.8),
        axis.title.x = element_text(vjust= 0.5),
        axis.text.x = element_text(angle = 90, vjust= 0.5),
        axis.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(limits = c(0, 1.10), breaks=c(0,0.25,0.5,0.75,1.0)) +
  labs(y = "Infection Frequency", x = "Site") +
  ggtitle("Year 2019")

p5


# ---year 2020-----
y20 <- subset(t2.long.na, year == "X2020")

trapData <- tibble()
splitData <- split(y20, y20$pop)

for(i in 1:length(splitData)){
  trapData[i,1] <- splitData[[i]]$id[1] #id
  trapData[i,2] <- splitData[[i]]$pop[1] #pop
  trapData[i,3] <- splitData[[i]]$latitude[1] # lat
  trapData[i,4] <- splitData[[i]]$longitude[1] #lon
  trapData[i,5] <- splitData[[i]]$year[1] # year
  trapData[i,6] <- splitData[[i]]$infection[1] # infection
  trapData[i,7] <- splitData[[i]]$inf[1] # number of infected individuals
  trapData[i,8] <- splitData[[i]]$uninf[1] # number of uninfected individuals
  #trapData[i,5] <- length(which(splitData[[i]]$Infected == "I" )) #n infected families
  #trapData[i,6] <- length(which(splitData[[i]]$Infected == "U" )) #n uninfected families
  
  #Estimate infection freq. and 95% binomial confidence intervals
  trapData[i,9] <- binconf(x=as.numeric(trapData[i,7]), n=as.numeric(trapData[i,7] + trapData[i,8]))[1] #infection freq.
  trapData[i,10] <- binconf(x=as.numeric(trapData[i,7]), n=as.numeric(trapData[i,7] + trapData[i,8]))[2]
  trapData[i,11] <- binconf(x=as.numeric(trapData[i,7]), n=as.numeric(trapData[i,7] + trapData[i,8]))[3]
  
} 
rm(i)
rm(splitData)

colnames(trapData) <- c("id", "pop", "lat", "lon", 
                        "year", "infection", "nInf", "nUnInf", "infectionFreq", "infectionFreq_lowerCI", "infectionFreq_upperCI")

# round the decimals
trapData[,9:11] <- round(trapData[,9:11], digits=3)

# re-order the data frame based on spatial location of the cline
trapData.ord <- trapData[match(y20$pop, trapData$pop),]

#plot yakuba infection freq. data
p6 <- ggplot(trapData.ord, aes(x=factor(pop, level = pop), y=infectionFreq)) + 
  geom_errorbar(aes(ymin=infectionFreq_lowerCI, 
                    ymax=infectionFreq_upperCI), 
                width=.2, size=1.5) +
  geom_point(size=4) +
  theme_bw() +
  theme(axis.title.y = element_text(vjust= 1.8),
        axis.title.x = element_text(vjust= 0.5),
        axis.text.x = element_text(angle = 90, vjust= 0.5),
        axis.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(limits = c(0, 1.10), breaks=c(0,0.25,0.5,0.75,1.0)) +
  labs(y = "Infection Frequency", x = "Site") +
  ggtitle("Year 2020")

p6

p7 <- ggarrange(p1, p2, p3, p4, p5, p6)



# ---location Ivan 2-----
t2.long.na = t2.long.latlon.pasohl23hayani.na
ivan2 <- subset(t2.long.na, pop == "Ivan 2")

trapData <- tibble()
splitData <- split(ivan2, ivan2$year)

for(i in 1:length(splitData)){
  trapData[i,1] <- splitData[[i]]$id[1] #id
  trapData[i,2] <- splitData[[i]]$pop[1] #pop
  trapData[i,3] <- splitData[[i]]$latitude[1] # lat
  trapData[i,4] <- splitData[[i]]$longitude[1] #lon
  trapData[i,5] <- splitData[[i]]$year[1] # year
  trapData[i,6] <- splitData[[i]]$infection[1] # infection
  trapData[i,7] <- splitData[[i]]$inf[1] # number of infected individuals
  trapData[i,8] <- splitData[[i]]$uninf[1] # number of uninfected individuals
  #trapData[i,5] <- length(which(splitData[[i]]$Infected == "I" )) #n infected families
  #trapData[i,6] <- length(which(splitData[[i]]$Infected == "U" )) #n uninfected families
  
  #Estimate infection freq. and 95% binomial confidence intervals
  trapData[i,9] <- binconf(x=as.numeric(trapData[i,7]), n=as.numeric(trapData[i,7] + trapData[i,8]))[1] #infection freq.
  trapData[i,10] <- binconf(x=as.numeric(trapData[i,7]), n=as.numeric(trapData[i,7] + trapData[i,8]))[2]
  trapData[i,11] <- binconf(x=as.numeric(trapData[i,7]), n=as.numeric(trapData[i,7] + trapData[i,8]))[3]
  
} 
rm(i)
rm(splitData)

colnames(trapData) <- c("id", "pop", "lat", "lon", 
                        "year", "infection", "nInf", "nUnInf", "infectionFreq", "infectionFreq_lowerCI", "infectionFreq_upperCI")

# round the decimals
trapData[,9:11] <- round(trapData[,9:11], digits=3)

# re-order the data frame based on spatial location of the cline
#trapData.ord <- trapData[match(y15$pop, trapData$pop),]

#plot yakuba infection freq. data
pl1 <- ggplot(trapData, aes(x=factor(year, level = year), y=infectionFreq)) + 
  geom_errorbar(aes(ymin=infectionFreq_lowerCI, 
                    ymax=infectionFreq_upperCI), 
                width=.2, size=1.5) +
  geom_point(size=4) +
  theme_bw() +
  theme(axis.title.y = element_text(vjust= 1.8),
        axis.title.x = element_text(vjust= 0.5),
        axis.text.x = element_text(angle = 90, vjust= 0.5),
        axis.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(limits = c(0, 1.10), breaks=c(0,0.25,0.5,0.75,1.0)) +
  labs(y = "Infection Frequency", x = "Site") +
  ggtitle("location Ivan2")
pl1

# ---location Velke Nemcice 1-----
velnem1 <- subset(t2.long.na, pop == "Velke Nemcice 1")

trapData <- tibble()
splitData <- split(velnem1, velnem1$year)

for(i in 1:length(splitData)){
  trapData[i,1] <- splitData[[i]]$id[1] #id
  trapData[i,2] <- splitData[[i]]$pop[1] #pop
  trapData[i,3] <- splitData[[i]]$latitude[1] # lat
  trapData[i,4] <- splitData[[i]]$longitude[1] #lon
  trapData[i,5] <- splitData[[i]]$year[1] # year
  trapData[i,6] <- splitData[[i]]$infection[1] # infection
  trapData[i,7] <- splitData[[i]]$inf[1] # number of infected individuals
  trapData[i,8] <- splitData[[i]]$uninf[1] # number of uninfected individuals
  #trapData[i,5] <- length(which(splitData[[i]]$Infected == "I" )) #n infected families
  #trapData[i,6] <- length(which(splitData[[i]]$Infected == "U" )) #n uninfected families
  
  #Estimate infection freq. and 95% binomial confidence intervals
  trapData[i,9] <- binconf(x=as.numeric(trapData[i,7]), n=as.numeric(trapData[i,7] + trapData[i,8]))[1] #infection freq.
  trapData[i,10] <- binconf(x=as.numeric(trapData[i,7]), n=as.numeric(trapData[i,7] + trapData[i,8]))[2]
  trapData[i,11] <- binconf(x=as.numeric(trapData[i,7]), n=as.numeric(trapData[i,7] + trapData[i,8]))[3]
  
} 
rm(i)
rm(splitData)

colnames(trapData) <- c("id", "pop", "lat", "lon", 
                        "year", "infection", "nInf", "nUnInf", "infectionFreq", "infectionFreq_lowerCI", "infectionFreq_upperCI")

# round the decimals
trapData[,9:11] <- round(trapData[,9:11], digits=3)

# re-order the data frame based on spatial location of the cline
#trapData.ord <- trapData[match(y15$pop, trapData$pop),]

#plot yakuba infection freq. data
pl2 <- ggplot(trapData, aes(x=factor(year, level = year), y=infectionFreq)) + 
  geom_errorbar(aes(ymin=infectionFreq_lowerCI, 
                    ymax=infectionFreq_upperCI), 
                width=.2, size=1.5) +
  geom_point(size=4) +
  theme_bw() +
  theme(axis.title.y = element_text(vjust= 1.8),
        axis.title.x = element_text(vjust= 0.5),
        axis.text.x = element_text(angle = 90, vjust= 0.5),
        axis.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(limits = c(0, 1.10), breaks=c(0,0.25,0.5,0.75,1.0)) +
  labs(y = "Infection Frequency", x = "Site") +
  ggtitle("location Velike Nemcice1")

pl2


# ---location NosislavZidlochovice 2-----
nonzid2 <- subset(t2.long.na, pop == "NosislavZidlochovice 2")

trapData <- tibble()
splitData <- split(nonzid2, nonzid2$year)

for(i in 1:length(splitData)){
  trapData[i,1] <- splitData[[i]]$id[1] #id
  trapData[i,2] <- splitData[[i]]$pop[1] #pop
  trapData[i,3] <- splitData[[i]]$latitude[1] # lat
  trapData[i,4] <- splitData[[i]]$longitude[1] #lon
  trapData[i,5] <- splitData[[i]]$year[1] # year
  trapData[i,6] <- splitData[[i]]$infection[1] # infection
  trapData[i,7] <- splitData[[i]]$inf[1] # number of infected individuals
  trapData[i,8] <- splitData[[i]]$uninf[1] # number of uninfected individuals
  #trapData[i,5] <- length(which(splitData[[i]]$Infected == "I" )) #n infected families
  #trapData[i,6] <- length(which(splitData[[i]]$Infected == "U" )) #n uninfected families
  
  #Estimate infection freq. and 95% binomial confidence intervals
  trapData[i,9] <- binconf(x=as.numeric(trapData[i,7]), n=as.numeric(trapData[i,7] + trapData[i,8]))[1] #infection freq.
  trapData[i,10] <- binconf(x=as.numeric(trapData[i,7]), n=as.numeric(trapData[i,7] + trapData[i,8]))[2]
  trapData[i,11] <- binconf(x=as.numeric(trapData[i,7]), n=as.numeric(trapData[i,7] + trapData[i,8]))[3]
  
} 
rm(i)
rm(splitData)

colnames(trapData) <- c("id", "pop", "lat", "lon", 
                        "year", "infection", "nInf", "nUnInf", "infectionFreq", "infectionFreq_lowerCI", "infectionFreq_upperCI")

# round the decimals
trapData[,9:11] <- round(trapData[,9:11], digits=3)

# re-order the data frame based on spatial location of the cline
#trapData.ord <- trapData[match(y15$pop, trapData$pop),]

#plot yakuba infection freq. data
pl3 <- ggplot(trapData, aes(x=factor(year, level = year), y=infectionFreq)) + 
  geom_errorbar(aes(ymin=infectionFreq_lowerCI, 
                    ymax=infectionFreq_upperCI), 
                width=.2, size=1.5) +
  geom_point(size=4) +
  theme_bw() +
  theme(axis.title.y = element_text(vjust= 1.8),
        axis.title.x = element_text(vjust= 0.5),
        axis.text.x = element_text(angle = 90, vjust= 0.5),
        axis.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(limits = c(0, 1.10), breaks=c(0,0.25,0.5,0.75,1.0)) +
  labs(y = "Infection Frequency", x = "Site") +
  ggtitle("location NosislavZidlochovice 2")

pl3


##--------------Make location groupings---------------------------------------------
#install.packages('circular')
library(circular)
# remove locations with NA for lon and/or lat
t2.na <-subset(t2, !is.na(latitude) & !is.na(longitude))
dist <- c()
for(i in 1:nrow(t2.na)) {
  
  row <- t2.na[i,] # extract row
  
  # approximate radius of earth in km
  R = 6373.0
  lat1 = rad(48.79926)
  lon1 = rad(16.49265)
  lat2 = rad(row$latitude)
  lon2 = rad(row$longitude)
  
  dlon = lon2 - lon1
  dlat = lat2 - lat1
  
  a = sin(dlat / 2)**2 + cos(lat1) * cos(lat2) * sin(dlon / 2)**2
  c = 2 * atan2(sqrt(a), sqrt(1 - a))
  # distance in meters
  distance = (R * c)*1000
  #print(distance)
  dist  <- c(dist, distance) # save as a vector
}
dist

# add distance as a column
t2.na$dist <- dist
# divide the dataferame into bins based on the distance. 
# A bin is defined as a set of location within a geographic separation of 10km; therefore we will have 8 bins
t2.na.bin <- t2.na %>% mutate(locbins =
                     case_when(dist < 10000 ~ "bin1",
                               dist > 10000 & dist < 20000 ~ "bin2",
                               dist > 20000 & dist < 30000 ~ "bin3",
                               dist > 30000 & dist < 40000 ~ "bin4",
                               dist > 40000 & dist < 50000 ~ "bin5",
                               dist > 50000 & dist < 60000 ~ "bin6",
                               dist > 50000 & dist < 60000 ~ "bin7",
                               dist > 60000 & dist < 70000 ~ "bin8")
)

# trasnform to long format
t2.na.bin.long <- t2.na.bin %>% as_tibble() %>% pivot_longer(cols = c("X2015", "X2016", "X2017", "X2018", "X2019", "X2020"), names_to ="year", values_to = "infection")
t2.na.bin.long$total <- 25
t2.na.bin.long$inf <- (t2.na.bin.long$infection*t2.na.bin.long$total)/100
t2.na.bin.long$uninf <- t2.na.bin.long$total - t2.na.bin.long$inf
#t2.long[is.na(t2.long)] = 0
t2.na.bin.long[!is.na(t2.na.bin.long$inf),]
t2.na.bin.long.na <- na.omit(t2.na.bin.long)

#summarise data for each elevation
finalData <- tibble()
splitData <- split(t2.na.bin.long.na, t2.na.bin.long.na$locbins)

for(i in 1:length(splitData)){
  finalData[i,1] <- splitData[[i]]$id[1] #id
  finalData[i,2] <- splitData[[i]]$locbins[1] #pop
  finalData[i,3] <- splitData[[i]]$latitude[1] # lat
  finalData[i,4] <- splitData[[i]]$longitude[1] #lon
  finalData[i,5] <- splitData[[i]]$year[1] # year
  finalData[i,6] <- splitData[[i]]$infection[1] # infection
  finalData[i,7] <- splitData[[i]]$inf[1] # number of infected individuals
  finalData[i,8] <- splitData[[i]]$uninf[1] # number of uninfected individuals
  #trapData[i,5] <- length(which(splitData[[i]]$Infected == "I" )) #n infected families
  #trapData[i,6] <- length(which(splitData[[i]]$Infected == "U" )) #n uninfected families
  
  #Estimate infection freq. and 95% binomial confidence intervals
  finalData[i,9] <- binconf(x=as.numeric(finalData[i,7]), n=as.numeric(finalData[i,7] + finalData[i,8]))[1] #infection freq.
  finalData[i,10] <- binconf(x=as.numeric(finalData[i,7]), n=as.numeric(finalData[i,7] + finalData[i,8]))[2]
  finalData[i,11] <- binconf(x=as.numeric(finalData[i,7]), n=as.numeric(finalData[i,7] + finalData[i,8]))[3]
  
} 
rm(i)
rm(splitData)

colnames(finalData) <- c("id", "bin", "lat", "lon", 
                         "year", "infection", "nInf", "nUnInf", "infectionFreq", "infectionFreq_lowerCI", "infectionFreq_upperCI")

#finalData[,3] <- round(finalData[,3], digits=3)
finalData[,9:11] <- round(finalData[,9:11], digits=3)


#plot yakuba infection freq. data
p.bin <- ggplot(finalData, aes(x=factor(bin, level = bin), y=infectionFreq)) + 
  geom_errorbar(aes(ymin=infectionFreq_lowerCI, 
                    ymax=infectionFreq_upperCI), 
                width=.2, size=1.5) +
  geom_point(size=4) +
  theme_bw() +
  theme(axis.title.y = element_text(vjust= 1.8),
        axis.title.x = element_text(vjust= 0.5),
        axis.text.x = element_text(angle = 90, vjust= 0.5),
        axis.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(limits = c(0, 1.10), breaks=c(0,0.25,0.5,0.75,1.0)) +
  labs(y = "Infection Frequency", x = "Bins")

p.bin


# ---bin1-----
bin1 <- subset(t2.na.bin.long.na, locbins == "bin1")

#summarise data for each elevation
finalData <- tibble()
splitData <- split(bin1, bin1$year)

for(i in 1:length(splitData)){
  finalData[i,1] <- splitData[[i]]$id[1] #id
  finalData[i,2] <- splitData[[i]]$locbins[1] #pop
  finalData[i,3] <- splitData[[i]]$latitude[1] # lat
  finalData[i,4] <- splitData[[i]]$longitude[1] #lon
  finalData[i,5] <- splitData[[i]]$year[1] # year
  finalData[i,6] <- splitData[[i]]$infection[1] # infection
  finalData[i,7] <- splitData[[i]]$inf[1] # number of infected individuals
  finalData[i,8] <- splitData[[i]]$uninf[1] # number of uninfected individuals
  
  #Estimate infection freq. and 95% binomial confidence intervals
  finalData[i,9] <- binconf(x=as.numeric(finalData[i,7]), n=as.numeric(finalData[i,7] + finalData[i,8]))[1] #infection freq.
  finalData[i,10] <- binconf(x=as.numeric(finalData[i,7]), n=as.numeric(finalData[i,7] + finalData[i,8]))[2]
  finalData[i,11] <- binconf(x=as.numeric(finalData[i,7]), n=as.numeric(finalData[i,7] + finalData[i,8]))[3]
  
} 
rm(i)
rm(splitData)

colnames(finalData) <- c("id", "bin", "lat", "lon", 
                         "year", "infection", "nInf", "nUnInf", "infectionFreq", "infectionFreq_lowerCI", "infectionFreq_upperCI")

finalData[,9:11] <- round(finalData[,9:11], digits=3)

#plotinfection freq. data
p.bin1 <- ggplot(finalData, aes(x=factor(year, level = year), y=infectionFreq)) + 
  geom_errorbar(aes(ymin=infectionFreq_lowerCI, 
                    ymax=infectionFreq_upperCI), 
                width=.2, size=1.5) +
  geom_point(size=4) +
  theme_bw() +
  theme(axis.title.y = element_text(vjust= 1.8),
        axis.title.x = element_text(vjust= 0.5),
        axis.text.x = element_text(angle = 90, vjust= 0.5),
        axis.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(limits = c(0, 1.10), breaks=c(0,0.25,0.5,0.75,1.0)) +
  labs(y = "Infection Frequency", x = "year") +
  ggtitle("Binned locations bin1")

p.bin1

# ---bin2-----
bin2 <- subset(t2.na.bin.long.na, locbins == "bin2")

#summarise data for each elevation
finalData <- tibble()
splitData <- split(bin2, bin2$year)

for(i in 1:length(splitData)){
  finalData[i,1] <- splitData[[i]]$id[1] #id
  finalData[i,2] <- splitData[[i]]$locbins[1] #pop
  finalData[i,3] <- splitData[[i]]$latitude[1] # lat
  finalData[i,4] <- splitData[[i]]$longitude[1] #lon
  finalData[i,5] <- splitData[[i]]$year[1] # year
  finalData[i,6] <- splitData[[i]]$infection[1] # infection
  finalData[i,7] <- splitData[[i]]$inf[1] # number of infected individuals
  finalData[i,8] <- splitData[[i]]$uninf[1] # number of uninfected individuals
  
  #Estimate infection freq. and 95% binomial confidence intervals
  finalData[i,9] <- binconf(x=as.numeric(finalData[i,7]), n=as.numeric(finalData[i,7] + finalData[i,8]))[1] #infection freq.
  finalData[i,10] <- binconf(x=as.numeric(finalData[i,7]), n=as.numeric(finalData[i,7] + finalData[i,8]))[2]
  finalData[i,11] <- binconf(x=as.numeric(finalData[i,7]), n=as.numeric(finalData[i,7] + finalData[i,8]))[3]
  
} 
rm(i)
rm(splitData)

colnames(finalData) <- c("id", "bin", "lat", "lon", 
                         "year", "infection", "nInf", "nUnInf", "infectionFreq", "infectionFreq_lowerCI", "infectionFreq_upperCI")

finalData[,9:11] <- round(finalData[,9:11], digits=3)

#plotinfection freq. data
p.bin2 <- ggplot(finalData, aes(x=factor(year, level = year), y=infectionFreq)) + 
  geom_errorbar(aes(ymin=infectionFreq_lowerCI, 
                    ymax=infectionFreq_upperCI), 
                width=.2, size=1.5) +
  geom_point(size=4) +
  theme_bw() +
  theme(axis.title.y = element_text(vjust= 1.8),
        axis.title.x = element_text(vjust= 0.5),
        axis.text.x = element_text(angle = 90, vjust= 0.5),
        axis.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(limits = c(0, 1.10), breaks=c(0,0.25,0.5,0.75,1.0)) +
  labs(y = "Infection Frequency", x = "year") +
  ggtitle("Binned locations bin2")

p.bin2


# ---bin3-----
bin3 <- subset(t2.na.bin.long.na, locbins == "bin3")

#summarise data for each elevation
finalData <- tibble()
splitData <- split(bin3, bin3$year)

for(i in 1:length(splitData)){
  finalData[i,1] <- splitData[[i]]$id[1] #id
  finalData[i,2] <- splitData[[i]]$locbins[1] #pop
  finalData[i,3] <- splitData[[i]]$latitude[1] # lat
  finalData[i,4] <- splitData[[i]]$longitude[1] #lon
  finalData[i,5] <- splitData[[i]]$year[1] # year
  finalData[i,6] <- splitData[[i]]$infection[1] # infection
  finalData[i,7] <- splitData[[i]]$inf[1] # number of infected individuals
  finalData[i,8] <- splitData[[i]]$uninf[1] # number of uninfected individuals
  
  #Estimate infection freq. and 95% binomial confidence intervals
  finalData[i,9] <- binconf(x=as.numeric(finalData[i,7]), n=as.numeric(finalData[i,7] + finalData[i,8]))[1] #infection freq.
  finalData[i,10] <- binconf(x=as.numeric(finalData[i,7]), n=as.numeric(finalData[i,7] + finalData[i,8]))[2]
  finalData[i,11] <- binconf(x=as.numeric(finalData[i,7]), n=as.numeric(finalData[i,7] + finalData[i,8]))[3]
  
} 
rm(i)
rm(splitData)

colnames(finalData) <- c("id", "bin", "lat", "lon", 
                         "year", "infection", "nInf", "nUnInf", "infectionFreq", "infectionFreq_lowerCI", "infectionFreq_upperCI")

finalData[,9:11] <- round(finalData[,9:11], digits=3)

#plotinfection freq. data
p.bin3 <- ggplot(finalData, aes(x=factor(year, level = year), y=infectionFreq)) + 
  geom_errorbar(aes(ymin=infectionFreq_lowerCI, 
                    ymax=infectionFreq_upperCI), 
                width=.2, size=1.5) +
  geom_point(size=4) +
  theme_bw() +
  theme(axis.title.y = element_text(vjust= 1.8),
        axis.title.x = element_text(vjust= 0.5),
        axis.text.x = element_text(angle = 90, vjust= 0.5),
        axis.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(limits = c(0, 1.10), breaks=c(0,0.25,0.5,0.75,1.0)) +
  labs(y = "Infection Frequency", x = "year") +
  ggtitle("Binned locations bin3")

p.bin3


# ---bin4-----
bin4 <- subset(t2.na.bin.long.na, locbins == "bin4")

#summarise data for each elevation
finalData <- tibble()
splitData <- split(bin4, bin4$year)

for(i in 1:length(splitData)){
  finalData[i,1] <- splitData[[i]]$id[1] #id
  finalData[i,2] <- splitData[[i]]$locbins[1] #pop
  finalData[i,3] <- splitData[[i]]$latitude[1] # lat
  finalData[i,4] <- splitData[[i]]$longitude[1] #lon
  finalData[i,5] <- splitData[[i]]$year[1] # year
  finalData[i,6] <- splitData[[i]]$infection[1] # infection
  finalData[i,7] <- splitData[[i]]$inf[1] # number of infected individuals
  finalData[i,8] <- splitData[[i]]$uninf[1] # number of uninfected individuals
  
  #Estimate infection freq. and 95% binomial confidence intervals
  finalData[i,9] <- binconf(x=as.numeric(finalData[i,7]), n=as.numeric(finalData[i,7] + finalData[i,8]))[1] #infection freq.
  finalData[i,10] <- binconf(x=as.numeric(finalData[i,7]), n=as.numeric(finalData[i,7] + finalData[i,8]))[2]
  finalData[i,11] <- binconf(x=as.numeric(finalData[i,7]), n=as.numeric(finalData[i,7] + finalData[i,8]))[3]
  
} 
rm(i)
rm(splitData)

colnames(finalData) <- c("id", "bin", "lat", "lon", 
                         "year", "infection", "nInf", "nUnInf", "infectionFreq", "infectionFreq_lowerCI", "infectionFreq_upperCI")

finalData[,9:11] <- round(finalData[,9:11], digits=3)

#plotinfection freq. data
p.bin4 <- ggplot(finalData, aes(x=factor(year, level = year), y=infectionFreq)) + 
  geom_errorbar(aes(ymin=infectionFreq_lowerCI, 
                    ymax=infectionFreq_upperCI), 
                width=.2, size=1.5) +
  geom_point(size=4) +
  theme_bw() +
  theme(axis.title.y = element_text(vjust= 1.8),
        axis.title.x = element_text(vjust= 0.5),
        axis.text.x = element_text(angle = 90, vjust= 0.5),
        axis.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(limits = c(0, 1.10), breaks=c(0,0.25,0.5,0.75,1.0)) +
  labs(y = "Infection Frequency", x = "year") +
  ggtitle("Binned locations bin4")

p.bin4


# ---bin5-----
bin5 <- subset(t2.na.bin.long.na, locbins == "bin5")

#summarise data for each elevation
finalData <- tibble()
splitData <- split(bin5, bin5$year)

for(i in 1:length(splitData)){
  finalData[i,1] <- splitData[[i]]$id[1] #id
  finalData[i,2] <- splitData[[i]]$locbins[1] #pop
  finalData[i,3] <- splitData[[i]]$latitude[1] # lat
  finalData[i,4] <- splitData[[i]]$longitude[1] #lon
  finalData[i,5] <- splitData[[i]]$year[1] # year
  finalData[i,6] <- splitData[[i]]$infection[1] # infection
  finalData[i,7] <- splitData[[i]]$inf[1] # number of infected individuals
  finalData[i,8] <- splitData[[i]]$uninf[1] # number of uninfected individuals
  
  #Estimate infection freq. and 95% binomial confidence intervals
  finalData[i,9] <- binconf(x=as.numeric(finalData[i,7]), n=as.numeric(finalData[i,7] + finalData[i,8]))[1] #infection freq.
  finalData[i,10] <- binconf(x=as.numeric(finalData[i,7]), n=as.numeric(finalData[i,7] + finalData[i,8]))[2]
  finalData[i,11] <- binconf(x=as.numeric(finalData[i,7]), n=as.numeric(finalData[i,7] + finalData[i,8]))[3]
  
} 
rm(i)
rm(splitData)

colnames(finalData) <- c("id", "bin", "lat", "lon", 
                         "year", "infection", "nInf", "nUnInf", "infectionFreq", "infectionFreq_lowerCI", "infectionFreq_upperCI")

finalData[,9:11] <- round(finalData[,9:11], digits=3)

#plotinfection freq. data
p.bin5 <- ggplot(finalData, aes(x=factor(year, level = year), y=infectionFreq)) + 
  geom_errorbar(aes(ymin=infectionFreq_lowerCI, 
                    ymax=infectionFreq_upperCI), 
                width=.2, size=1.5) +
  geom_point(size=4) +
  theme_bw() +
  theme(axis.title.y = element_text(vjust= 1.8),
        axis.title.x = element_text(vjust= 0.5),
        axis.text.x = element_text(angle = 90, vjust= 0.5),
        axis.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(limits = c(0, 1.10), breaks=c(0,0.25,0.5,0.75,1.0)) +
  labs(y = "Infection Frequency", x = "year") +
  ggtitle("Binned locations bin5")

p.bin5


# ---bin6-----
bin6 <- subset(t2.na.bin.long.na, locbins == "bin6")

#summarise data for each elevation
finalData <- tibble()
splitData <- split(bin6, bin6$year)

for(i in 1:length(splitData)){
  finalData[i,1] <- splitData[[i]]$id[1] #id
  finalData[i,2] <- splitData[[i]]$locbins[1] #pop
  finalData[i,3] <- splitData[[i]]$latitude[1] # lat
  finalData[i,4] <- splitData[[i]]$longitude[1] #lon
  finalData[i,5] <- splitData[[i]]$year[1] # year
  finalData[i,6] <- splitData[[i]]$infection[1] # infection
  finalData[i,7] <- splitData[[i]]$inf[1] # number of infected individuals
  finalData[i,8] <- splitData[[i]]$uninf[1] # number of uninfected individuals
  
  #Estimate infection freq. and 95% binomial confidence intervals
  finalData[i,9] <- binconf(x=as.numeric(finalData[i,7]), n=as.numeric(finalData[i,7] + finalData[i,8]))[1] #infection freq.
  finalData[i,10] <- binconf(x=as.numeric(finalData[i,7]), n=as.numeric(finalData[i,7] + finalData[i,8]))[2]
  finalData[i,11] <- binconf(x=as.numeric(finalData[i,7]), n=as.numeric(finalData[i,7] + finalData[i,8]))[3]
  
} 
rm(i)
rm(splitData)

colnames(finalData) <- c("id", "bin", "lat", "lon", 
                         "year", "infection", "nInf", "nUnInf", "infectionFreq", "infectionFreq_lowerCI", "infectionFreq_upperCI")

finalData[,9:11] <- round(finalData[,9:11], digits=3)

#plotinfection freq. data
p.bin6 <- ggplot(finalData, aes(x=factor(year, level = year), y=infectionFreq)) + 
  geom_errorbar(aes(ymin=infectionFreq_lowerCI, 
                    ymax=infectionFreq_upperCI), 
                width=.2, size=1.5) +
  geom_point(size=4) +
  theme_bw() +
  theme(axis.title.y = element_text(vjust= 1.8),
        axis.title.x = element_text(vjust= 0.5),
        axis.text.x = element_text(angle = 90, vjust= 0.5),
        axis.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(limits = c(0, 1.10), breaks=c(0,0.25,0.5,0.75,1.0)) +
  labs(y = "Infection Frequency", x = "year") +
  ggtitle("Binned locations bin6")

p.bin6


# ---bin8-----
bin8 <- subset(t2.na.bin.long.na, locbins == "bin8")

#summarise data for each elevation
finalData <- tibble()
splitData <- split(bin8, bin8$year)

for(i in 1:length(splitData)){
  finalData[i,1] <- splitData[[i]]$id[1] #id
  finalData[i,2] <- splitData[[i]]$locbins[1] #pop
  finalData[i,3] <- splitData[[i]]$latitude[1] # lat
  finalData[i,4] <- splitData[[i]]$longitude[1] #lon
  finalData[i,5] <- splitData[[i]]$year[1] # year
  finalData[i,6] <- splitData[[i]]$infection[1] # infection
  finalData[i,7] <- splitData[[i]]$inf[1] # number of infected individuals
  finalData[i,8] <- splitData[[i]]$uninf[1] # number of uninfected individuals
  
  #Estimate infection freq. and 95% binomial confidence intervals
  finalData[i,9] <- binconf(x=as.numeric(finalData[i,7]), n=as.numeric(finalData[i,7] + finalData[i,8]))[1] #infection freq.
  finalData[i,10] <- binconf(x=as.numeric(finalData[i,7]), n=as.numeric(finalData[i,7] + finalData[i,8]))[2]
  finalData[i,11] <- binconf(x=as.numeric(finalData[i,7]), n=as.numeric(finalData[i,7] + finalData[i,8]))[3]
  
} 
rm(i)
rm(splitData)

colnames(finalData) <- c("id", "bin", "lat", "lon", 
                         "year", "infection", "nInf", "nUnInf", "infectionFreq", "infectionFreq_lowerCI", "infectionFreq_upperCI")

finalData[,9:11] <- round(finalData[,9:11], digits=3)

#plotinfection freq. data
p.bin8 <- ggplot(finalData, aes(x=factor(year, level = year), y=infectionFreq)) + 
  geom_errorbar(aes(ymin=infectionFreq_lowerCI, 
                    ymax=infectionFreq_upperCI), 
                width=.2, size=1.5) +
  geom_point(size=4) +
  theme_bw() +
  theme(axis.title.y = element_text(vjust= 1.8),
        axis.title.x = element_text(vjust= 0.5),
        axis.text.x = element_text(angle = 90, vjust= 0.5),
        axis.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(limits = c(0, 1.10), breaks=c(0,0.25,0.5,0.75,1.0)) +
  labs(y = "Infection Frequency", x = "year") +
  ggtitle("Binned locations bin8")

p.bin8

p.bins1to8 <- ggarrange(p.bin1, p.bin2, p.bin3, p.bin4, p.bin5, p.bin6, p.bin8)


#--------analysis of mu and F using determinitic models assuming no CI (equation 1)--------------------
#First, let's assume that low levels of CI can be ignored, and evaluate an equilibrium between imperfect
#trans and positive fitness effects. With F(1-mu)>1, the stable equilibrium is 

##yak_low analysis
#Lower CI
mu <- 0.003 #Lower CI
# (1 - mu) * relF > 1
relF <- 1 / (1 - mu)
relF

relF <- seq(0.5,4, by=0.01) #range of relative fitness values
phat <- 1 - ((mu*relF)/(relF-1)) #stable equilibrium frequency calculation base on mu and relF

trapData.ord$infectionFreq_lowerCI[1] #Lower CI
relF[which.min(abs(phat - trapData.ord$infectionFreq_lowerCI[1]))]
trapData.ord$infectionFreq[1] #infection freq.
relF[which.min(abs(phat - trapData.ord$infectionFreq[1]))] #F value require to explain infection freq.
trapData.ord$infectionFreq_upperCI[1] #Upper CI
relF[which.min(abs(phat - trapData.ord$infectionFreq_upperCI[1]))]


#Upper CI
mu <- 0.184 #Upper CI
# (1 - mu) * relF > 1
relF <- 1 / (1 - mu)
relF

relF <- seq(1,4, by=0.01) #range of relative fitness values
phat <- 1 - ((mu*relF)/(relF-1)) #stable equilibrium frequency calculation base on mu and relF

trapData.ord$infectionFreq_lowerCI[1] #Lower CI
relF[which.min(abs(phat - trapData.ord$infectionFreq_lowerCI[1]))]
trapData.ord$infectionFreq[1] #infection freq.
relF[which.min(abs(phat - trapData.ord$infectionFreq[1]))] #F value require to explain infection freq.
trapData.ord$infectionFreq_upperCI[1] #Upper CI
relF[which.min(abs(phat - trapData.ord$infectionFreq_upperCI[1]))]


#Upper CI
mu <- 0.1 #Upper CI
# (1 - mu) * relF > 1
relF <- 1 / (1 - mu)
relF

relF <- seq(1,4, by=0.01) #range of relative fitness values
phat <- 1 - ((mu*relF)/(relF-1)) #stable equilibrium frequency calculation base on mu and relF

trapData$infectionFreq_lowerCI[1] #Lower CI
relF[which.min(abs(phat - trapData$infectionFreq_lowerCI[1]))]
trapData$infectionFreq[1] #infection freq.
relF[which.min(abs(phat - trapData$infectionFreq[1]))] #F value require to explain infection freq.
trapData$infectionFreq_upperCI[1] #Upper CI
relF[which.min(abs(phat - trapData$infectionFreq_upperCI[1]))]


#-------------------------Plots-----------------------------------------------------------
# no effect of CI and a low mu
relF <- seq(1.0, 1.5, by=0.01)
mu <- 0.2 # maternal transmission
# different strengths of CI with the low mu
phatCI_sh1 <- 1 - ((mu*relF)/(relF-1)) 
is.na(phatCI_sh1) <- sapply(phatCI_sh1, is.infinite)

sh <- (0.98) #experimental estimate of CI
phatCI_sh2 <- (sh + 1 - relF + sqrt((sh + 1 - relF)^2 + 4* sh * ((relF - relF*mu) - 1) * (1 - relF * mu)))/(2*sh*(1-relF*mu)) 

sh<-(0.45) #moderate CI
phatCI_sh3 <- (sh + 1 - relF + sqrt((sh + 1 - relF)^2 + 4* sh * ((relF - relF*mu) - 1) * (1 - relF * mu)))/(2*sh*(1-relF*mu)) 

sh<-(0.20) #low CI
phatCI_sh4 <- (sh + 1 - relF + sqrt((sh + 1 - relF)^2 + 4* sh * ((relF - relF*mu) - 1) * (1 - relF * mu)))/(2*sh*(1-relF*mu)) 

yak_low_data <- as.data.frame(cbind(relF, phatCI_sh1
                                    , phatCI_sh2, phatCI_sh3, phatCI_sh4
                                    ))
yak_low_data <- melt(yak_low_data, id.vars="relF")

p13 <- ggplot(data=yak_low_data, aes(x=relF, y=value, col=variable)) +
  geom_line() +
  theme_classic() +
  scale_color_manual(values=c('black','blue', 'orange', 'red'), labels = c("sh = 0", "sh = 0.98", "sh = 0.45", "sh = 0.10"), 
                     guide = guide_legend(reverse=TRUE)) +
  scale_y_continuous(limits=c(0,1)) +
  geom_hline(yintercept=mean(trapData$infectionFreq), linetype="dashed", color = "black") +
  geom_hline(yintercept=mean(trapData$infectionFreq_lowerCI), color = "grey") +
  geom_hline(yintercept=mean(trapData$infectionFreq_upperCI), color = "grey") +
  labs(y = "phat", x = "relative fitness", color = "sh") +
  ggtitle(paste("low mu =", mu))
p13


# no effect of CI with a higher mu
relF <- seq(1.0, 1.5, by=0.01)
mu <- 0.1
# different strengths of CI with the higher mu
phatCI_sh1 <- 1 - ((mu*relF)/(relF-1)) 
is.na(phatCI_sh1) <- sapply(phatCI_sh1, is.infinite)

yak_low_data <- as.data.frame(cbind(relF, phatCI_sh1
                                    , phatCI_sh2, phatCI_sh3, phatCI_sh4
))
yak_low_data <- melt(yak_low_data, id.vars="relF")

p14 <- ggplot(data=yak_low_data, aes(x=relF, y=value, col=variable)) +
  geom_line() +
  theme_classic() +
  scale_color_manual(values=c('black','blue', 'orange', 'red'), labels = c("sh = 0", "sh = 0.98", "sh = 0.45", "sh = 0.10"), 
                     guide = guide_legend(reverse=TRUE)) +
  scale_y_continuous(limits=c(0,1)) +
  geom_hline(yintercept=trapData.ord$infectionFreq[1], linetype="dashed", color = "black") +
  geom_hline(yintercept=trapData.ord$infectionFreq_lowerCI[1], color = "grey") +
  geom_hline(yintercept=trapData.ord$infectionFreq_upperCI[1], color = "grey") +
  labs(y = "phat", x = "relative fitness", color = "sh") +
  ggtitle(paste("high mu =", mu))
p14


#yak_low upper CI of mu
relF <- seq(1.0, 1.5, by=0.01)
mu <- 0.23

phatCI_sh1 <- 1 - ((mu*relF)/(relF-1)) 
is.na(phatCI_sh1) <- sapply(phatCI_sh1, is.infinite)

yak_low_data <- as.data.frame(cbind(relF, phatCI_sh1
                                    , phatCI_sh2, phatCI_sh3, phatCI_sh4
))
yak_low_data <- melt(yak_low_data, id.vars="relF")

p15 <- ggplot(data=yak_low_data, aes(x=relF, y=value, col=variable)) +
  geom_line() +
  scale_color_manual(values=c('black','blue', 'orange', 'red'), labels = c("sh = 0", "sh = 0.98", "sh = 0.45", "sh = 0.10"), 
                     guide = guide_legend(reverse=TRUE)) +
  scale_y_continuous(limits=c(0,1)) +
  geom_hline(yintercept=trapData$infectionFreq[1], linetype="dashed", color = "black") +
  geom_hline(yintercept=trapData$infectionFreq_lowerCI[1], color = "grey") +
  geom_hline(yintercept=trapData$infectionFreq_upperCI[1], color = "grey") +
  labs(y = "phat", x = "relative fitness", color = "sh") +
  ggtitle(paste("low mu =", mu))
p15


pmu <- ggarrange(p13, p14, p15)


#------------------ fit the T2 cline and estimate phat for each year with an ML model ------ 

#--- prepare the data for model fitting -----#

##### caclulate distance from the x,y coordinates

library(circular)
# remove locations with NA for lon and/or lat
t2.long.na.na <-subset(t2.long.na, !is.na(latitude) & !is.na(longitude))
dist2 <- c()
for(i in 1:nrow(t2.long.na.na)) {
  
  row <- t2.long.na.na[i,] # extract row
  
  # approximate radius of earth in km
  R = 6373.0
  lat1 = rad(48.79926)
  lon1 = rad(16.49265)
  lat2 = rad(row$latitude)
  lon2 = rad(row$longitude)
  
  dlon = lon2 - lon1
  dlat = lat2 - lat1
  
  a = sin(dlat / 2)**2 + cos(lat1) * cos(lat2) * sin(dlon / 2)**2
  c = 2 * atan2(sqrt(a), sqrt(1 - a))
  # distance in meters
  distance = (R * c)*1000
  #print(distance)
  dist2  <- c(dist2, distance) # save as a vector
}
dist2

# add distance as a column
t2.long.na.na$dist <- dist2/1000

# transfrom years into generation from 1 to 6
t2.long.na.na.t <- t2.long.na.na %>%
  mutate(time = case_when(
    year == "X2015" ~ 1,
    year == "X2016" ~ 2,
    year == "X2017" ~ 3,
    year == "X2018" ~ 4,
    year == "X2019" ~ 5,
    year == "X2020" ~ 6
    
  ))

# subset only 
data = t2.long.na.na.t[,c(10, 11, 12)]


# the diffusion model
dist = data$dist
time = 6
phat = sf/sh

f <- 1/(1+exp(dist-((1-(2*(phat)))*time)))
#plot(x=data.nout$dist, y=data.nout$infectionFreq, xlab = "Distance", ylab = "Infection frequency", ylim=c(0,1))
plot(f, type = "l", ylim = c(0,1), xlim = c(0,70), lwd = 3, ylab = "Infection frequency", xlab = "Distance (km)")

# make a function for the diffusion approximation formula
dtfun <- function(dis, phat, t){
  # Prediction of the model
  1/(1+(exp(dis-((1-(2*phat))*t))))
}

# plot quickly to check
plot(x=data$dist, y=data$infectionFreq, xlab = "Distance", ylab = "Infection frequency", ylim=c(0,1))
curve(dtfun(dis=x, phat=0, t=6), data$dist, add = TRUE)



NLL = function(phat, t=6, data20) {
  # Values predicted by the model
  ppred = dtfun(data20$dist, phat, t=6)
  #print(ppred)
  # binomial distibution
  -sum(dbinom(x = ppred, 25, data20$infectionFreq), log=T)
  #dbinom(x = data$infectionFreq, 92, ppred)
  #-sum(dnorm(x = data$G, mean = Gpred, sd = pars[4], log = TRUE))
}
#data3 <- data[,c(2,1,3)]
#fit = optim(par = par0, fn = NLL, data = data, hessian = T)
fit = optim(par = c(-1.6), NLL, data=data20, lower=c(-Inf), upper=c(Inf))
fit$par

plot(x=data$dist, y=data$infectionFreq, xlab = "Distance", ylab = "Infection frequency")
curve(dtfun(x, -10.8, t=1), data15$dist, add = TRUE, lwd =2, col = "yellow", lty = 2 )
curve(dtfun(x, -5.8, t=2), data16$dist, add = TRUE, lwd =2, col = "orange", lty = 2 )
curve(dtfun(x, -3.42, t=3), data17$dist, add = TRUE, lwd =2, col = "red", lty = 2 )
curve(dtfun(x, -2.43, t=4), data18$dist, add = TRUE, lwd =2, col = "red2", lty = 2 )
curve(dtfun(x, -2.16, t=5), data19$dist, add = TRUE, lwd =2, col = "brown", lty = 2 )
curve(dtfun(x, -1.44, t=6), data20$dist, add = TRUE, lwd =2, col = "grey10", lty = 2 )


### the Mximum likelihood function ####
tidyML <- data.frame()
for(i in 1:nrow(data)){
  row = data[i,]
  NLL = function(phat, t=row$time, data) {
    # Values predicted by the model
    ppred = dtfun(data$dist, phat, t=row$time)
    #print(ppred)
    # binomial distibution
    -sum(dbinom(x = ppred, 25, data$infectionFreq), log=T)
  }
  if (row$time == 1) {
    fit = optim(par = c(-12), NLL, data=data ,lower=c(-Inf), upper=c(Inf)) 
  } else if (row$time == 2) {
    fit = optim(par = c(-5.8), NLL, data=data ,lower=c(-Inf), upper=c(Inf)) 
  } else if (row$time == 3) {
    fit = optim(par = c(-3.8), NLL, data=data ,lower=c(-Inf), upper=c(Inf)) 
  } else if (row$time == 4) {
    fit = optim(par = c(-2.7), NLL, data=data ,lower=c(-Inf), upper=c(Inf)) 
  } else if (row$time == 5) {
    fit = optim(par = c(-2.1), NLL, data=data ,lower=c(-Inf), upper=c(Inf)) 
  } else if (row$time == 6) {
    fit = optim(par = c(-1.6), NLL, data=data ,lower=c(-Inf), upper=c(Inf)) 
  }
  #fit = optim(par = c(1), NLL, data=data ,lower=c(-Inf), upper=c(Inf))
  #fit$par
  broomedMl <- as.data.frame(cbind(broom::tidy(fit), row$time))
  tidyML <- rbind(tidyML, broomedMl)
}

tiff("/Volumes/LaCie/the_cline/T2_MLfit.tiff", width = 8, height = 6, units = "in", res = 300)
plot(x=data$dist, y=data$infectionFreq, xlab = "Distance (km)", ylab = "Infection frequency", frame.plot = FALSE, ylim = c(0,1), xlim=c(0,50), pch = 19,  col= "grey60")
curve(dtfun(x, -10.80, t=1), data$dist, add = TRUE, lwd =2, col = "red3", lty = 2)
curve(dtfun(x, -5.22, t=2), data$dist, add = TRUE, lwd =2, col = "red3", lty = 2)
curve(dtfun(x, -3.42, t=3), data$dist, add = TRUE, lwd =2, col = "red3", lty = 2)
curve(dtfun(x, -2.43, t=4), data$dist, add = TRUE, lwd =2, col = "red3", lty = 2)
curve(dtfun(x, -1.89, t=5), data$dist, add = TRUE, lwd =2, col = "red3", lty = 2)
curve(dtfun(x, -1.44, t=6), data$dist, add = TRUE, lwd =2, col = "red3", lty = 2)
dev.off()

#----------------the NSL model fitting -------------------------#

# subset only infection frequency, distnce and generation time
data = t2.long.na.na.t[,c(10, 11, 12)]

data.nout <-t2.long.na.na.t[-c(5, 7,9),] 
dat <- subset(t2.long.na.na.t, year == "X2015")
dat20 <- subset(t2.long.na.na.t, year == "X2020")
dat19 <- subset(t2.long.na.na.t, year == "X2019")
dat18 <- subset(t2.long.na.na.t, year == "X2018")
dat17 <- subset(t2.long.na.na.t, year == "X2017")
dat15 <- subset(t2.long.na.na.t, year == "X2015")

#-----------------sigma fitting with an nls model---------------#
f <- function(x, sigma, t){0.5 * (1 - tanh((x - sigma * (1 - 2 * 0) * sqrt(0.98) * t / 2)/(2 * sigma / sqrt(0.98))))}

#our infection frequencies
y=t2.long.na.na.t$infectionFreq
#distance from location 1 (calculated from online distance calculaors based on coordinates)
x=t2.long.na.na.t$dist
# time
t=t2.long.na.na.t$time

#confidence interval calculation see bellow *after CI calculation, infection frequency is subtracted by the low and high CI because this function just adds or substracts from the dot (for graphing purposes)
for(i in 1:nrow(t2.long.na.na.t)){
  
  #Estimate infection freq. and 95% binomial confidence intervals
  #t2.long.na.na.t[i,13] <- binconf(x=as.numeric(t2.long.na.na.t[i,8]), n=as.numeric(t2.long.na.na.t[i,8] + t2.long.na.na.t[i,9]))[1] #infection freq.
  t2.long.na.na.t[i,13] <- binconf(x=as.numeric(t2.long.na.na.t[i,8]), n=as.numeric(t2.long.na.na.t[i,8] + t2.long.na.na.t[i,9]))[2]
  t2.long.na.na.t[i,14] <- binconf(x=as.numeric(t2.long.na.na.t[i,8]), n=as.numeric(t2.long.na.na.t[i,8] + t2.long.na.na.t[i,9]))[3]
  
}

colnames(t2.long.na.na.t) <- c("id", "pop", "latitude", "longitude", 
                         "year", "infection", "total", "inf", "uninf", "infectionFreq", "dist", "time", "infectionFreq_lowerCI", "infectionFreq_upperCI")

# plot quickly to just check
plot(x,y)
curve(f(x, 3, 11), add = TRUE)

#find a best fitting sigma and t
sigma.t.fit <- nls(y~f(x,sigma,t), start = list(sigma=1, t=3))

# wave speed formula = σ*√lCI)/2
t2Speed <- (broom::tidy(sigma.t.fit)$estimate[1]*sqrt(0.98))/2

#residual sum of squares and R value
RSS <- sum(residuals(sigma.t.fit)^2)
TSS <- sum((y - mean(y))^2)
Rnls.t2 <- 1 - (RSS/TSS)
Rnls.t2


# plot
tiff("/Volumes/LaCie/the_cline/T2_sigmaAll.tiff", width = 8, height = 6, units = "in", res = 300)
#plot(t2.long.na.na.t$dist, t2.long.na.na.t$infectionFreq, lwd = "1", ylab="Infection Frequency", xlab="Distance(km)", frame.plot = FALSE, ylim = c(0,1), xlim=c(0,50), pch = 19,  col= "grey60")
plotCI(x = t2.long.na.na.t$dist,               # plotrix plot with confidence intervals
       y = t2.long.na.na.t$infectionFreq,
       li = t2.long.na.na.t$infectionFreq_lowerCI,
       ui = t2.long.na.na.t$infectionFreq_upperCI, 
       lwd = "1", ylab="Infection Frequency", xlab="Distance(km)", frame.plot = FALSE, ylim = c(0,1), xlim=c(0,50), pch = 19,  col= "grey60")
curve(f(x,broom::tidy(sigma.t.fit)$estimate[1],t=broom::tidy(sigma.t.fit)$estimate[2]), add = TRUE, lwd =2, col = "blue", lty = 2 )
text(x= 40, y= 0.9, paste("sigma=", round(broom::tidy(sigma.t.fit)$estimate[1],3)), cex=1.1)
text(x= 40, y= 0.8, paste("wave speed=", round(t2Speed,3)), cex=1.1)
text(x= 40, y= 0.7, paste("R=", round(Rnls.t2,3)), cex=1.1)
dev.off()


#------------Nls fitting of the phat for each t---------------#
f <- function(x, phat, t){1/(1+(exp(x-((1-(2*phat))*t))))}


y=t2.long.na.na.t$infectionFreq
x=t2.long.na.na.t$dist
t=t2.long.na.na.t$time

# plot quickly to just check
plot(x,y)
curve(f(x, -2, 5), add = TRUE)


#library(broom)
tidy <- tibble()
for(i in 1:nrow(t2.long.na.na.t)) {
  row = t2.long.na.na.t[i,]
  sigma.t2.fit <- nls(y~f(x, phat, row$time), start = list(phat=0)
                      #,lower=c(0), algorithm = "port"
  )
  #print(sigma.t.fit) 
  broomed <- as.data.frame(cbind(broom::tidy(sigma.t2.fit), row$time))
  tidy <- rbind(tidy, broomed)
}

sigma.t2.fit
summary(sigma.t2.fit)

#residual sum of squares and R value
RSS <- sum(residuals(sigma.t2.fit)^2)
TSS <- sum((y - mean(y))^2)
Rnls.t2 <- 1 - (RSS/TSS)
Rnls.t2


### nls fitted phat 
tiff("/Volumes/LaCie/the_cline/T2_nlsfit.tiff", width = 8, height = 6, units = "in", res = 300)
plot(x,y, lwd = "1", ylab="Infection Frequency", xlab="Distance(km)", ylim = c(0,1), xlim=c(0,50), pch = 19,  col= "grey60")
curve(dtfun(x, tidy[which(tidy$`row$time` == 1),]$estimate[1], t=1), t2.long.na.na.t$dist, add = TRUE, lwd =2, col = "red3", lty = 2)
curve(dtfun(x, tidy[which(tidy$`row$time` == 2),]$estimate[1], t=2), t2.long.na.na.t$dist, add = TRUE, lwd =2, col = "red3", lty = 2)
curve(dtfun(x, tidy[which(tidy$`row$time` == 3),]$estimate[1], t=3), t2.long.na.na.t$dist, add = TRUE, lwd =2, col = "red3", lty = 2)
curve(dtfun(x, tidy[which(tidy$`row$time` == 4),]$estimate[1], t=4), t2.long.na.na.t$dist, add = TRUE, lwd =2, col = "red3", lty = 2)
curve(dtfun(x, tidy[which(tidy$`row$time` == 5),]$estimate[1], t=5), t2.long.na.na.t$dist, add = TRUE, lwd =2, col = "red3", lty = 2)
curve(dtfun(x, tidy[which(tidy$`row$time` == 6),]$estimate[1], t=6), t2.long.na.na.t$dist, add = TRUE, lwd =2, col = "red3", lty = 2)
text(x= 40, y= 0.9, paste("t1_phat=", round(tidy[which(tidy$`row$time` == 1),]$estimate[1],2)), cex=1)
text(x= 40, y= 0.85, paste("t2_phat=", round(tidy[which(tidy$`row$time` == 2),]$estimate[1],2)), cex=1)
text(x= 40, y= 0.8, paste("t3_phat=", round(tidy[which(tidy$`row$time` == 3),]$estimate[1],2)), cex=1)
text(x= 40, y= 0.75, paste("t4_phat=", round(tidy[which(tidy$`row$time` == 4),]$estimate[1],2)), cex=1)
text(x= 40, y= 0.7, paste("t5_phat=", round(tidy[which(tidy$`row$time` == 5),]$estimate[1],2)), cex=1)
text(x= 40, y= 0.65, paste("t6_phat=", round(tidy[which(tidy$`row$time` == 6),]$estimate[1],2)), cex=1)
#text(x= 15, y= 0.6, paste("sigma=", round(broom::tidy(sigma.t.fit)$estimate[1],3)), cex=1.1)
dev.off()



#------------------------------nls model fitting with a scaled X----------------------------#
# subset only infection frequncy, distnce and generation time
data = t2.long.na.na.t[,c(10, 11, 12)]
# scaled X
data <- data %>% add_column(scaledX = sqrt((2*1*(data$dist^2))/2.927^2))
### the Mximum likelihood function ####
tidyML <- data.frame()
for(i in 1:nrow(data)){
  row = data[i,]
  NLL = function(phat, t=row$time, data) {
    # Values predicted by the model
    ppred = dtfun(data$scaledX, phat, t=row$time)
    #print(ppred)
    # binomial distibution
    -sum(dbinom(x = ppred, 25, data$infectionFreq), log=T)
  }
  if (row$time == 1) {
    fit = optimize(interval = c(-3, 1), NLL, lower=c(-3), upper=c(1)) 
  } else if (row$time == 2) {
    fit = optimize(interval = c(-3, 1), NLL, lower=c(-3), upper=c(1)) 
  } else if (row$time == 3) {
    fit = optimize(interval = c(-3, 1), NLL, lower=c(-3), upper=c(1)) 
  } else if (row$time == 4) {
    fit = optimize(interval = c(-3, 1), NLL, lower=c(-3), upper=c(1)) 
  } else if (row$time == 5) {
    fit = optimize(interval = c(-3, 1), NLL, lower=c(-3), upper=c(1)) 
  } else if (row$time == 6) {
    fit = optimize(interval = c(-3, 1), NLL, lower=c(-3), upper=c(1)) 
  }
  #fit = optim(par = c(1), NLL, data=data ,lower=c(-Inf), upper=c(Inf))
  #fit$par
  broomedMl <- as.data.frame(cbind(broom::tidy(fit), row$time))
  tidyML <- rbind(tidyML, broomedMl)
}


#tiff("/Volumes/LaCie/the_cline/T2_MLfit.tiff", width = 8, height = 6, units = "in", res = 300)
plot(x=bla$scaledX, y=bla$infectionFreq, xlab = "Distance (km)", ylab = "Infection frequency", frame.plot = FALSE, ylim = c(0,1), xlim=c(0,20), pch = 19,  col= "grey60")
curve(dtfun(x, tidyML[which(tidyML$`row$time` == 1),]$value[1], t=1), bla$scaledX, add = TRUE, lwd =2, col = "red3", lty = 2)
curve(dtfun(x, tidyML[which(tidyML$`row$time` == 2),]$value[1], t=2), bla$scaledX, add = TRUE, lwd =2, col = "red3", lty = 2)
curve(dtfun(x, tidyML[which(tidyML$`row$time` == 3),]$value[1], t=3), bla$scaledX, add = TRUE, lwd =2, col = "red3", lty = 2)
curve(dtfun(x, tidyML[which(tidyML$`row$time` == 4),]$value[1], t=4), bla$scaledX, add = TRUE, lwd =2, col = "red3", lty = 2)
curve(dtfun(x, tidyML[which(tidyML$`row$time` == 5),]$value[1], t=5), bla$scaledX, add = TRUE, lwd =2, col = "red3", lty = 2)
curve(dtfun(x, tidyML[which(tidyML$`row$time` == 6),]$value[1], t=6), bla$scaledX, add = TRUE, lwd =2, col = "red3", lty = 2)
#dev.off()

#---- nsl fitting ----##

y=data$infectionFreq
x=data$scaledX
t=data$time

tidy <- tibble()
for(i in 1:nrow(data)) {
  row = data[i,]
  sigma.t2.fit <- nls(y~f(x, phat, row$time), start = list(phat=0)
                      #,lower=c(0), algorithm = "port"
  )
  #print(sigma.t.fit) 
  broomed <- as.data.frame(cbind(broom::tidy(sigma.t2.fit), row$time))
  tidy <- rbind(tidy, broomed)
}

sigma.t2.fit
summary(sigma.t2.fit)
#library(plotrix)
tiff("/Volumes/LaCie/the_cline/T2_NLSfit_scaledX.tiff", width = 8, height = 6, units = "in", res = 300)
plotCI(x = data$scaledX,               # plotrix plot with confidence intervals
     y = data$infectionFreq,
     li = t2.long.na.na.t$infectionFreq_lowerCI,
     ui = t2.long.na.na.t$infectionFreq_upperCI, xlab = "Distance (km)", ylab = "Infection frequency", frame.plot = FALSE, ylim = c(0,1), xlim=c(0,35), pch = 19,  col= "grey60")
curve(dtfun(x, tidy[which(tidy$`row$time` == 1),]$estimate[1], t=1), data$scaledX, add = TRUE, lwd =2, col = "red3", lty = 2)
curve(dtfun(x, tidy[which(tidy$`row$time` == 2),]$estimate[1], t=2), data$scaledX, add = TRUE, lwd =2, col = "red3", lty = 2)
curve(dtfun(x, tidy[which(tidy$`row$time` == 3),]$estimate[1], t=3), data$scaledX, add = TRUE, lwd =2, col = "red3", lty = 2)
curve(dtfun(x, tidy[which(tidy$`row$time` == 4),]$estimate[1], t=4), data$scaledX, add = TRUE, lwd =2, col = "red3", lty = 2)
curve(dtfun(x, tidy[which(tidy$`row$time` == 5),]$estimate[1], t=5), data$scaledX, add = TRUE, lwd =2, col = "red3", lty = 2)
curve(dtfun(x, tidy[which(tidy$`row$time` == 6),]$estimate[1], t=6), data$scaledX, add = TRUE, lwd =2, col = "red3", lty = 2)
text(x= 30, y= 0.9, paste("t1_phat=", round(tidy[which(tidy$`row$time` == 1),]$estimate[1],2)), cex=1)
text(x= 30, y= 0.85, paste("t2_phat=", round(tidy[which(tidy$`row$time` == 2),]$estimate[1],2)), cex=1)
text(x= 30, y= 0.8, paste("t3_phat=", round(tidy[which(tidy$`row$time` == 3),]$estimate[1],2)), cex=1)
text(x= 30, y= 0.75, paste("t4_phat=", round(tidy[which(tidy$`row$time` == 4),]$estimate[1],2)), cex=1)
text(x= 30, y= 0.7, paste("t5_phat=", round(tidy[which(tidy$`row$time` == 5),]$estimate[1],2)), cex=1)
text(x= 30, y= 0.65, paste("t6_phat=", round(tidy[which(tidy$`row$time` == 6),]$estimate[1],2)), cex=1)
text(x= 30, y= 0.6, paste("sigma=", round(broom::tidy(sigma.t.fit)$estimate[1],3)), cex=1.1)
#text(x= 40, y= 0.7, paste("R=", round(Rnls.t2,3)), cex=1.1)
dev.off()


##----------------two-dimensional space---------------##
f <- function(x, phat, t, k){1/(2*(0.5 - phat))*(1 + log(x) * (exp(((2*(0.5 - phat))^2)*(t + k)-1)))}


# add the MO transect to this
all <- cline
all <- as_tibble(all)

all.long <- all %>% as_tibble() %>% pivot_longer(cols = c("X2015", "X2016", "X2017", "X2018", "X2019", "X2020"), names_to ="year", values_to = "infection")

# add an imaginary counts as if the sample size for all locs/years was 25
all.long$total <- 25
all.long$inf <- (all.long$infection*all.long$total)/100
all.long$uninf <- all.long$total - all.long$inf
#t2.long[is.na(t2.long)] = 0
all.long[!is.na(all.long$inf),]
all.long.na <- na.omit(all.long)
all.long.na$infectionFreq <- all.long.na$inf/all.long.na$total


# remove locations with NA for lon and/or lat
all.long.na.na <-subset(all.long.na, !is.na(latitude) & !is.na(longitude))
dist2all <- c()
for(i in 1:nrow(all.long.na.na)) {
  
  row <- all.long.na.na[i,] # extract row
  
  # approximate radius of earth in km
  R = 6373.0
  lat1 = rad(48.79926)
  lon1 = rad(16.49265)
  lat2 = rad(row$latitude)
  lon2 = rad(row$longitude)
  
  dlon = lon2 - lon1
  dlat = lat2 - lat1
  
  a = sin(dlat / 2)**2 + cos(lat1) * cos(lat2) * sin(dlon / 2)**2
  c = 2 * atan2(sqrt(a), sqrt(1 - a))
  # distance in meters
  distance = (R * c)*1000
  #print(distance)
  dist2all  <- c(dist2all, distance) # save as a vector
}
dist2all

# add distance as a column
all.long.na.na$dist <- dist2all/1000

# plot the radius of the 2D cline (T2 + MO)
library(leaflet)
central_point_coordinates <- c(lng = 16.49265, lat = 48.79926)
distance_max <- max(all.long.na.na$dist)
map_background <- leaflet() %>% 
  addTiles() %>% 
  setView(lng = central_point_coordinates[1], lat = central_point_coordinates[2], zoom = 10) %>%
  addCircles(lng = as.numeric(central_point_coordinates[1]), lat = as.numeric(central_point_coordinates[2]), radius = distance_max * 1000, color = "yellow") 

map_background


# transfrom years into generation from 1 to 6
all.long.na.na.t <- all.long.na.na %>%
  mutate(time = case_when(
    year == "X2015" ~ 1,
    year == "X2016" ~ 2,
    year == "X2017" ~ 3,
    year == "X2018" ~ 4,
    year == "X2019" ~ 5,
    year == "X2020" ~ 6
    
  ))

# extract infection freqiency, distance and time columns
allData = all.long.na.na.t[,c(10, 11, 12)]

# calculate scaled distance 
funsig <- function(x, sigma, t){0.5 * (1 - tanh((x - sigma * (1 - 2 * 0) * sqrt(0.98) * t / 2)/(2 * sigma / sqrt(0.98))))}
y=allData$infectionFreq
x=allData$dist
t=allData$time
sigma.t2mo.fit <- nls(y~funsig(x,sigma,t), start = list(sigma=1, t=1))
# wave speed formula = σ*√lCI)/2
t2moSpeed <- (broom::tidy(sigma.t2mo.fit)$estimate[1]*sqrt(0.98))/2

# caluclate scaled X with the calculated sigma
allData <- allData %>% add_column(scaledX = sqrt((2*1*(allData$dist^2))/broom::tidy(sigma.t2mo.fit)$estimate[1]^2))
#allData <- allData %>% add_column(scaledR = sqrt((2*1*(max(all.long.na.na$dist))^2/broom::tidy(sigma.t2mo.fit)$estimate[1]^2)))

y=allData$infectionFreq
x=allData$scaledX
t=allData$time

nlc <- nls.control(maxiter= 1000, warnOnly=TRUE)
tidyk <- tibble()
for(i in 1:nrow(allData)) {
  row = allData[i,]
  sigma.t2.fit <- nls(y~f(x, phat=-2, row$time, k), start = list(k=0)
                      ,control = nlc
                     )
  #print(sigma.t.fit) 
  broomed <- as.data.frame(cbind(broom::tidy(sigma.t2.fit), row$time))
  tidyk <- rbind(tidyk, broomed)
}

# numerical solution for k for the large time examined (t=6)
k6 <- tidyk[which(tidyk$`row$time` == 6),]$estimate[1]

nlc <- nls.control(maxiter= 1000, warnOnly=TRUE)
tidy <- tibble()
for(i in 1:nrow(allData)) {
  row = allData[i,]
  if (row$time == 1) {
    sigma.t2.fit <- nls(y~f(x, phat, row$time, k=k6), start = list(phat=0),control = nlc)
  } else if (row$time == 2) {
    sigma.t2.fit <- nls(y~f(x, phat, row$time, k=k6), start = list(phat=0),control = nlc)
  } else if (row$time == 3) {
    sigma.t2.fit <- nls(y~f(x, phat, row$time, k=k6), start = list(phat=0),control = nlc)
  } else if (row$time == 4) {
    sigma.t2.fit <- nls(y~f(x, phat, row$time, k=k6), start = list(phat=0),control = nlc)
  } else if (row$time == 5) {
    sigma.t2.fit <- nls(y~f(x, phat, row$time, k=k6), start = list(phat=0),control = nlc)
  } else if (row$time == 6) {
    sigma.t2.fit <- nls(y~f(x, phat, row$time, k=k6), start = list(phat=0),control = nlc)
  }
 # sigma.t2.fit <- nls(y~f(x, phat, row$time, k=-1.1), start = list(phat=-1.5)
  #                    ,control = nlc)
  #print(sigma.t.fit) 
  broomed <- as.data.frame(cbind(broom::tidy(sigma.t2.fit), row$time))
  tidy <- rbind(tidy, broomed)
}

sigma.t2.fit
#residual sum of squares and R value
RSS <- sum(residuals(sigma.t2.fit)^2)
TSS <- sum((y - mean(y))^2)
Rnls.t2 <- 1 - (RSS/TSS)
Rnls.t2

#-----speed wave for p(R,T)=0.5 ----#
phats <- unique(tidy)
cs <- data.frame()
ts <- data.frame()
for (i in 1:nrow(phats)) {
  row = phats[i,]
  c <- (2*(0.5 - row$estimate)) - 1/0.5
  t <- row$`row$time`
  ts <- rbind(ts, t)
  cs <- rbind(cs, c)
  csts <- cbind(cs, ts)
}
colnames(csts) <- c("c", "t")

tiff("/Volumes/LaCie/the_cline/t2_2D_scaledX_CI.tiff", width = 8, height = 6, units = "in", res = 300)
plot(x=allData$scaledX, y=allData$infectionFreq, xlab = "Distance (km)", ylab = "Infection frequency", frame.plot = FALSE, ylim = c(0,1), xlim = c(0, 50), pch = 19, col= "grey60", main = "two dimensions")
curve(dtfun(x, tidy[which(tidy$`row$time` == 1),]$estimate[1], t=1), x, add = TRUE, lwd =2, col = "red3", lty = 2)
curve(dtfun(x, tidy[which(tidy$`row$time` == 2),]$estimate[1], t=2), x, add = TRUE, lwd =2, col = "red3", lty = 2)
curve(dtfun(x, tidy[which(tidy$`row$time` == 3),]$estimate[1], t=3), x, add = TRUE, lwd =2, col = "red3", lty = 2)
curve(dtfun(x, tidy[which(tidy$`row$time` == 4),]$estimate[1], t=4), x, add = TRUE, lwd =2, col = "red3", lty = 2)
curve(dtfun(x, tidy[which(tidy$`row$time` == 5),]$estimate[1], t=5), x, add = TRUE, lwd =2, col = "red3", lty = 2)
curve(dtfun(x, tidy[which(tidy$`row$time` == 6),]$estimate[1], t=6), x, add = TRUE, lwd =2, col = "red3", lty = 2)
text(x= 40, y= 0.9, paste("c_t1=", round(csts[which(csts$t == 1),]$c,2)), cex=1)
text(x= 40, y= 0.85, paste("c_t2=", round(csts[which(csts$t == 2),]$c,2)), cex=1)
text(x= 40, y= 0.8, paste("c_t3=", round(csts[which(csts$t == 3),]$c,2)), cex=1)
text(x= 40, y= 0.75, paste("c_t4=", round(csts[which(csts$t == 4),]$c,2)), cex=1)
text(x= 40, y= 0.7, paste("c_t5=", round(csts[which(csts$t == 5),]$c,2)), cex=1)
text(x= 40, y= 0.65, paste("c_t6=", round(csts[which(csts$t == 6),]$c,2)), cex=1)
text(x= 40, y= 0.6, paste("sigma=", round(broom::tidy(sigma.t2mo.fit)$estimate[1],2)), cex=1.1)
dev.off()


sf = seq(-1,1, by=0.1)
sh = 1
sv = -1 #it feels that only after the CI effect picks up, the negative Wol effects come into play
phat = (sf + sv - sf*sv)/sh 

d1 = 1
p = 0.5
fp = (sh*d1*p*(1- p)*(p - phat))/(1 - sf*p - sh*p*(1 - p))
fp



#----------decomposing phat in 1D model-----------#
f <- function(x, sr, sh, t){1/(1+(exp(x-((1-(2*((sr)/sh)))*t))))}

y=data$infectionFreq
x=data$scaledX
t=data$time
#sv=0
sh=0.98
#sf=0

tidy <- tibble()
for(i in 1:nrow(data)) {
  row = data[i,]
  sigma.t2.fit <- nls(y~f(x, sr, sh, row$time), start = list(sr=0)
                      #,lower=c(0), algorithm = "port"
  )
  #print(sigma.t.fit) 
  broomed <- as.data.frame(cbind(broom::tidy(sigma.t2.fit), row$time))
  tidy <- rbind(tidy, broomed)
}

sigma.t2.fit
summary(sigma.t2.fit)

RSS <- sum(residuals(sigma.t2.fit)^2)
TSS <- sum((y - mean(y))^2)
Rnls.t2 <- 1 - (RSS/TSS)
Rnls.t2

tiff("/Volumes/LaCie/the_cline/t2_1D_scaledX_sr.tiff", width = 8, height = 6, units = "in", res = 300)
plotCI(x=data$scaledX, y=data$infectionFreq,
       li = t2.long.na.na.t$infectionFreq_lowerCI,
       ui = t2.long.na.na.t$infectionFreq_upperCI, xlab = "Distance (km)", ylab = "Infection frequency", frame.plot = FALSE, ylim = c(0,1), xlim=c(0,35), pch = 19,  col= "grey60")
curve(dtfun(x, tidy[which(tidy$`row$time` == 1),]$estimate[1], t=1), data$scaledX, add = TRUE, lwd =2, col = "red3", lty = 2)
curve(dtfun(x, tidy[which(tidy$`row$time` == 2),]$estimate[1], t=2), data$scaledX, add = TRUE, lwd =2, col = "red3", lty = 2)
curve(dtfun(x, tidy[which(tidy$`row$time` == 3),]$estimate[1], t=3), data$scaledX, add = TRUE, lwd =2, col = "red3", lty = 2)
curve(dtfun(x, tidy[which(tidy$`row$time` == 4),]$estimate[1], t=4), data$scaledX, add = TRUE, lwd =2, col = "red3", lty = 2)
curve(dtfun(x, tidy[which(tidy$`row$time` == 5),]$estimate[1], t=5), data$scaledX, add = TRUE, lwd =2, col = "red3", lty = 2)
curve(dtfun(x, tidy[which(tidy$`row$time` == 6),]$estimate[1], t=6), data$scaledX, add = TRUE, lwd =2, col = "red3", lty = 2)
text(x= 30, y= 0.9, paste("sr_t1=", round(tidy[which(tidy$`row$time` == 1),]$estimate[1],2)), cex=1)
text(x= 30, y= 0.85, paste("sr_t2=", round(tidy[which(tidy$`row$time` == 2),]$estimate[1],2)), cex=1)
text(x= 30, y= 0.8, paste("sr_t3=", round(tidy[which(tidy$`row$time` == 3),]$estimate[1],2)), cex=1)
text(x= 30, y= 0.75, paste("sr_t4=", round(tidy[which(tidy$`row$time` == 4),]$estimate[1],2)), cex=1)
text(x= 30, y= 0.7, paste("sr_t5=", round(tidy[which(tidy$`row$time` == 5),]$estimate[1],2)), cex=1)
text(x= 30, y= 0.65, paste("sr_t6=", round(tidy[which(tidy$`row$time` == 6),]$estimate[1],2)), cex=1)
text(x= 30, y= 0.6, paste("sigma=", round(broom::tidy(sigma.t.fit)$estimate[1],3)), cex=1.1)
dev.off()



########-------------------------Modeling fitness effect of Wolbachia--------------------------------########
t2.long.na.na.t$transect <- "t2"
all.long.na.na.t$transect <- c(rep("t2", 92), rep("mo", 34))

t2.long.na.na.t.mid <- subset(t2.long.na.na.t, dist > 10)
t2.long.na.na.t.mid <- subset(t2.long.na.na.t, pop == "NosislavZidlochovice 2")
t2.long.na.na.t.mid <- subset(t2.long.na.na.t, pop == "Ivan 2")
table(t2.long.na.na.t$pop)
t2.long.na.na.t.mid <- subset(t2.long.na.na.t, pop == "Zidlochovice")
t2.long.na.na.t.mid <- subset(t2.long.na.na.t, pop == "Pasohlavky 3")
t2.long.na.na.t.mid <- subset(t2.long.na.na.t, pop == "Rajhrad")

confData <- tibble()
splitData <- split(t2.long.na.na.t, t2.long.na.na.t$transect)

for(i in 1:length(splitData)){
  confData[i,1] <- splitData[[i]]$id[1] #id
  confData[i,2] <- splitData[[i]]$pop[1] #pop
  confData[i,3] <- splitData[[i]]$latitude[1] # lat
  confData[i,4] <- splitData[[i]]$longitude[1] #lon
  confData[i,5] <- splitData[[i]]$time[1] # year
  confData[i,6] <- splitData[[i]]$infection[1] # infection
  confData[i,7] <- splitData[[i]]$inf[1] # number of infected individuals
  confData[i,8] <- splitData[[i]]$uninf[1] # number of uninfected individuals
  #trapData[i,5] <- length(which(splitData[[i]]$Infected == "I" )) #n infected families
  #trapData[i,6] <- length(which(splitData[[i]]$Infected == "U" )) #n uninfected families
  
  #Estimate infection freq. and 95% binomial confidence intervals
  confData[i,9] <- binconf(x=as.numeric(confData[i,7]), n=as.numeric(confData[i,7] + confData[i,8]))[1] #infection freq.
  confData[i,10] <- binconf(x=as.numeric(confData[i,7]), n=as.numeric(confData[i,7] + confData[i,8]))[2]
  confData[i,11] <- binconf(x=as.numeric(confData[i,7]), n=as.numeric(confData[i,7] + confData[i,8]))[3]
  
} 
rm(i)
rm(splitData)

colnames(confData) <- c("id", "pop", "lat", "lon", 
                        "year", "infection", "nInf", "nUnInf", "infectionFreq", "infectionFreq_lowerCI", "infectionFreq_upperCI")

# round the decimals
confData[,9:11] <- round(confData[,9:11], digits=3)

# re-order the data frame based on spatial location of the cline
#confData.ord <- confData[match(unique(y15$pop), confData$pop),]

#plot yakuba infection freq. data
p1 <- ggplot(confData, aes(x=factor(year, level = year), y=infectionFreq)) + 
  geom_errorbar(aes(ymin=infectionFreq_lowerCI, 
                    ymax=infectionFreq_upperCI), 
                width=.2, size=1.5) +
  geom_point(size=4) +
  theme_bw() +
  theme(axis.title.y = element_text(vjust= 1.8),
        axis.title.x = element_text(vjust= 0.5),
        axis.text.x = element_text(angle = 90, vjust= 0.5),
        axis.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(limits = c(0, 1.10), breaks=c(0,0.25,0.5,0.75,1.0)) +
  labs(y = "Infection Frequency", x = "Site") +
  ggtitle("Year 2015")
p1

ggplot(confData, aes(x = year, y = infectionFreq)) + 
  geom_line(col='red') + 
  geom_ribbon(aes(ymin = infectionFreq_lowerCI, ymax = infectionFreq_upperCI), alpha = 0.1) +
  theme_classic() +
  scale_x_continuous(limits = c(1, 6)) +
  scale_y_continuous(limits = c(0, 1)) 
  #+geom_line(data=yak_low_data, aes(x = y=value, col=variable))

# no effect of CI and a low mu
relF <- seq(0.9, 4, by=0.01)
mu <- 0.05 # maternal transmission
fiteffect <- relF*(1- mu) # fitness effect
# different strengths of CI with the low mu
phatCI_sh1 <- 1 - ((mu*relF)/(relF-1)) 
is.na(phatCI_sh1) <- sapply(phatCI_sh1, is.infinite)

sh <- (0.1) #experimental estimate of CI
phatCI_sh2 <- (sh + 1 - relF + sqrt((sh + 1 - relF)^2 + 4* sh * ((relF - relF*mu) - 1) * (1 - relF * mu)))/(2*sh*(1-relF*mu)) 

sh<-(0.25) #moderate CI
phatCI_sh3 <- (sh + 1 - relF + sqrt((sh + 1 - relF)^2 + 4* sh * ((relF - relF*mu) - 1) * (1 - relF * mu)))/(2*sh*(1-relF*mu)) 

sh<-(0.45) #low CI
phatCI_sh4 <- (sh + 1 - relF + sqrt((sh + 1 - relF)^2 + 4* sh * ((relF - relF*mu) - 1) * (1 - relF * mu)))/(2*sh*(1-relF*mu)) 

sh<-(0.98) #low CI
phatCI_sh5 <- (sh + 1 - relF + sqrt((sh + 1 - relF)^2 + 4* sh * ((relF - relF*mu) - 1) * (1 - relF * mu)))/(2*sh*(1-relF*mu)) 

yak_low_data <- as.data.frame(cbind(relF, phatCI_sh1
                                    , phatCI_sh2, phatCI_sh3, phatCI_sh4, phatCI_sh5
))
yak_low_data <- melt(yak_low_data, id.vars="relF")
yak_low_data$fiteffect <- rep(fiteffect, 5)


mu0003 <- ggplot(data=yak_low_data, aes(x=relF, y=value, col=variable)) +
  geom_line(data=subset(yak_low_data, fiteffect<=1), linetype=2) +
  geom_line(data=subset(yak_low_data, fiteffect>=1), linetype=1) +
  theme_classic() +
  scale_color_manual(values=c('black','blue', 'orange', 'forestgreen', 'red'), labels = c("sh = 0", "sh = 0.1", "sh = 0.25", "sh = 0.45", "sh=0.98"), 
                     guide = guide_legend(reverse=TRUE)) +
  scale_y_continuous(limits=c(0,1)) +
  geom_hline(yintercept=confData$infectionFreq, linetype="dashed", color = "black") +
  geom_hline(yintercept=confData$infectionFreq_lowerCI, color = "grey") +
  geom_hline(yintercept=confData$infectionFreq_upperCI, color = "grey") +
  labs(y = "equlibrium frequency (phat)", x = "relative fitness (F)", color = "sh") +
  ggtitle(paste("mu =", mu))
mu0003

mu001 <- ggplot(data=yak_low_data, aes(x=relF, y=value, col=variable)) +
  geom_line(data=subset(yak_low_data, fiteffect<=1), linetype=2) +
  geom_line(data=subset(yak_low_data, fiteffect>=1), linetype=1) +
  theme_classic() +
  scale_color_manual(values=c('black','blue', 'orange', 'forestgreen', 'red'), labels = c("sh = 0", "sh = 0.1", "sh = 0.25", "sh = 0.45", "sh=0.98"), 
                     guide = guide_legend(reverse=TRUE)) +
  scale_y_continuous(limits=c(0,1)) +
  geom_hline(yintercept=confData$infectionFreq, linetype="dashed", color = "black") +
  geom_hline(yintercept=confData$infectionFreq_lowerCI, color = "grey") +
  geom_hline(yintercept=confData$infectionFreq_upperCI, color = "grey") +
  labs(y = "equlibrium frequency (phat)", x = "relative fitness (F)", color = "sh") +
  ggtitle(paste("mu =", mu))
mu001

mu0025 <- ggplot(data=yak_low_data, aes(x=relF, y=value, col=variable)) +
  geom_line(data=subset(yak_low_data, fiteffect<=1), linetype=2) +
  geom_line(data=subset(yak_low_data, fiteffect>=1), linetype=1) +
  theme_classic() +
  scale_color_manual(values=c('black','blue', 'orange', 'forestgreen', 'red'), labels = c("sh = 0", "sh = 0.1", "sh = 0.25", "sh = 0.45", "sh=0.98"), 
                     guide = guide_legend(reverse=TRUE)) +
  scale_y_continuous(limits=c(0,1)) +
  geom_hline(yintercept=confData$infectionFreq, linetype="dashed", color = "black") +
  geom_hline(yintercept=confData$infectionFreq_lowerCI, color = "grey") +
  geom_hline(yintercept=confData$infectionFreq_upperCI, color = "grey") +
  labs(y = "equlibrium frequency (phat)", x = "relative fitness (F)", color = "sh") +
  ggtitle(paste("mu =", mu)) 
  # + theme(legend.position = "none")
mu0025

mu005 <- ggplot(data=yak_low_data, aes(x=relF, y=value, col=variable)) +
  geom_line(data=subset(yak_low_data, fiteffect<=1), linetype=2) +
  geom_line(data=subset(yak_low_data, fiteffect>=1), linetype=1) +
  theme_classic() +
  scale_color_manual(values=c('black','blue', 'orange', 'forestgreen', 'red'), labels = c("sh = 0", "sh = 0.1", "sh = 0.25", "sh = 0.45", "sh=0.98"), 
                     guide = guide_legend(reverse=TRUE)) +
  scale_y_continuous(limits=c(0,1)) +
  geom_hline(yintercept=confData$infectionFreq, linetype="dashed", color = "black") +
  geom_hline(yintercept=confData$infectionFreq_lowerCI, color = "grey") +
  geom_hline(yintercept=confData$infectionFreq_upperCI, color = "grey") +
  labs(y = "equlibrium frequency (phat)", x = "relative fitness (F)", color = "sh") +
  ggtitle(paste("mu =", mu))
  # + theme(legend.position = "none")
mu005

mus <- ggarrange(mu0003, mu001, mu0025, mu005, nrow = 1)
ggsave("/Volumes/LaCie/the_cline/mus_t2all_equilibrium.png", plot = mus, width = 16, height = 4, dpi = 300, units = "in", device='png')


##--------Interpretation of the Fisherian cline modeling where positive finess effect are modeled------##
# Observed frequency range is from 0 to about 0.96
# Stable equilibrium at very low frequncies can be produced when CI is non-existant and mu values higher (0.05).
# In contrast, to produce equilibria close to 0.96 strong selection comparable to F= 1.2 
# with CI levels higher than 0.1 and a low mu of 0.003 to 0.01.
# our frequency etimate of p=1 (1, 0.87) can be explained by both a very strong levels of 
# CI (lab estimate sh = 0.98), but also weak levels of CI when mu is low (0.003 or 0.01). Here, 
# F = 1.25 and F = 1.5 respectively, are required to explain the upper interval for p (1).
# If the mu is higher (mu > 0.025), sh > 0.25 and F = 4 are required to explain the upper interval for p (1)


#-----------logistic regression of the distance on infection frequency across all years together---------------#
mylogit <- glm(infectionFreq ~ dist, family=binomial(link="logit"), data = t2.long.na.na.t)
summary(mylogit)
res <- broom::tidy(mylogit)

#plot logistic regression curve
ggall <- ggplot(data=t2.long.na.na.t, aes(x=dist, y=infectionFreq)) + 
  theme_classic() +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14)) +
  geom_point(mapping = aes(size = total), col = "grey50", alpha=.5) +
  scale_size(name = expression("sample size"), breaks = c(7, 48)) +
  #geom_errorbar(aes(ymin=infectionFreq_lowerCI, 
  #                  ymax=infectionFreq_upperCI), 
  #              width=.02, size=0.3, col="grey60") +
  stat_smooth(method="glm", method.args = list(family=binomial(link="logit")),lty=1, se=T,
              formula = y ~ x) +
  scale_y_continuous(limits = c(0,1)) +
  annotate("text", x=60, y=0.80, label = paste("b = ", round(res$estimate[2], 3), "(±",(round(res$std.error[2], 3)),")"), size = 5) +
  annotate("text", x=60, y=0.75, label = paste("p < ", round(res$p.value[2], 5)), size = 5) +
  labs(x = "Distance (km)", y = "Infection frequency", color = "sh" , size = 10)

ggsave("/Volumes/LaCie/the_cline/logisticReg_t2all.png", plot = ggall, width = 10, height = 10, dpi = 300, units = "in", device='png')

#--accross years---#
yearslogit <- glm(infectionFreq ~ dist + year, family=binomial(link="logit"), data = t2.long.na.na.t)
summary(yearslogit)

ggbyyear <- ggplot(data=t2.long.na.na.t, aes(x=dist, y=infectionFreq, group = year, col=year)) + 
  theme_classic() +
  #geom_point(alpha=.5) +
  stat_smooth(method="glm", method.args = list(family=binomial(link="logit")),lty=1, se=T,
              formula = y ~ x, aes(fill = year)) +
  scale_y_continuous(limits = c(0,1)) +
  scale_colour_viridis_d(option = "plasma") +
  scale_fill_viridis_d(option = "plasma") +
  labs(x = "Distance (km)", y = "Infection frequency")
ggsave("/Volumes/LaCie/the_cline/logisticReg_t2yearsplit.png", plot = ggbyyear, width = 10, height = 10, dpi = 300, units = "in", device='png')

# also save as facets
ggbyyear.facet <- ggplot(data=t2.long.na.na.t, aes(x=dist, y=infectionFreq, group = year, col=year)) + 
  theme_classic() +
  #geom_point(alpha=.5) +
  stat_smooth(method="glm", method.args = list(family=binomial(link="logit")),lty=1, se=T,
              formula = y ~ x, aes(fill = year)) +
  scale_y_continuous(limits = c(0,1)) +
  scale_colour_viridis_d(option = "plasma") +
  scale_fill_viridis_d(option = "plasma") +
  labs(x = "Distance (km)", y = "Infection frequency") +
  facet_wrap(~year, nrow=2)
ggsave("/Volumes/LaCie/the_cline/logisticReg_t2yearsplit.facet.png", plot = ggbyyear.facet, width = 10, height = 10, dpi = 300, units = "in", device='png')


#----GLME for each year with pop as a mixed effect--------
# none of the years explains the infection frequency because this model does not allow for correcting for
# mixed effect of the locations. Therefore, I will use Generalized Linear Mixed Model with for 
# binomial distibution (logit) to account for this effect.
mixed <- glmer(infectionFreq ~ dist + year + (1|pop), data = t2.long.na.na.t, family = binomial)
isSingular(mixed, tol = 1e-4)
summary(mixed)

library(ggeffects) # for the function ggpredict

t2.glme <- plot(ggpredict(mixed, terms=c("dist[all]", "year"))) +
  theme_classic() +
  #theme_bw(panel.grid = element_blank()) +
  scale_y_continuous(limits = c(0,1)) +
  scale_colour_viridis_d(option = "plasma") +
  scale_fill_viridis_d(option = "plasma") +
  labs(x = "Distance (km)", y = "Infection frequency")
ggsave("/Volumes/LaCie/the_cline/glme_t2.png", plot = t2.glme, width = 8, height = 8, dpi = 300, units = "in", device='png')

# also save as facets
t2.glme.facet <- plot(facets = T, ggpredict(mixed, terms=c("dist[all]", "year"))) +
  theme_classic() +
  #theme_bw(panel.grid = element_blank()) +
  scale_y_continuous(limits = c(0,1)) +
  scale_colour_viridis_d(option = "plasma") +
  scale_fill_viridis_d(option = "plasma") +
  labs(x = "Distance (km)", y = "Infection frequency")
ggsave("/Volumes/LaCie/the_cline/glme_t2.facet.png", plot = t2.glme.facet, width = 8, height = 8, dpi = 300, units = "in", device='png')

# save logit distance, logit distance+year and glme model predictions in a table #
library(stargazer)
stargazer(mylogit, yearslogit, mixed, type = "html",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "",
          dep.var.labels=c("Infection frequency"),
          covariate.labels=c("Distance"),
          out="/Volumes/LaCie/the_cline/models.html")

#---logistic regression of the distance on infection frequency for each year separately now---#
y15 <- subset(t2.long.na.na.t, year == "X2015")
y15logit <- glm(infectionFreq ~ dist, family=binomial(link="logit"), data = y15)
summary(y15logit)
res15 <- broom::tidy(y15logit)

#plot logistic regression curve
gg15 <- ggplot(data=y15, aes(x=dist, y=infectionFreq)) + 
  theme_classic() +
  geom_point(mapping = aes(size = total), alpha=.5, show.legend = FALSE) +
  #stat_smooth(method="glm", method.args = list(family=binomial(link="logit")),lty=1, se=T,
  #            formula = y ~ x) +
  scale_y_continuous(limits = c(0,1)) +
  scale_x_continuous(limits = c(0,70)) +
  ggtitle(paste("year", y15$year)) +
  #annotate("text", x=60, y=0.80, label = paste("b = ", round(res15$estimate[2], 3), "(±",(round(res15$std.error[2], 3)),")")) +
  #annotate("text", x=60, y=0.75, label = paste("p* = ", round(res15$p.value[2], 5)), color = "red") +
  labs(x = "Distance (km)", y = "Infection frequency", color = "sh")

y16 <- subset(t2.long.na.na.t, year == "X2016")
y16logit <- glm(infectionFreq ~ dist, family=binomial(link="logit"), data = y16)
summary(y16logit)
res16 <- broom::tidy(y16logit)
#plot logistic regression curve
gg16 <- ggplot(data=y16, aes(x=dist, y=infectionFreq)) + 
  theme_classic() +
  geom_point(mapping = aes(size = total), alpha=.5, show.legend = FALSE) +
  #stat_smooth(method="glm", method.args = list(family=binomial(link="logit")),lty=1, se=T, formula = y ~ x) +
  scale_y_continuous(limits = c(0,1)) +
  scale_x_continuous(limits = c(0,70)) +
  ggtitle(paste("year", y16$year)) +
  #annotate("text", x=60, y=0.80, label = paste("b = ", round(res16$estimate[2], 3), "(±",(round(res16$std.error[2], 3)),")")) +
  #annotate("text", x=60, y=0.75, label = paste("p = ", round(res16$p.value[2], 5))) +
  labs(x = "Distance (km)", y = "Infection frequency", color = "sh")

y17 <- subset(t2.long.na.na.t, year == "X2017")
y17logit <- glm(infectionFreq ~ dist, family=binomial(link="logit"), data = y17)
summary(y17logit)
res17 <- broom::tidy(y17logit)
#plot logistic regression curve
gg17 <- ggplot(data=y17, aes(x=dist, y=infectionFreq)) + 
  theme_classic() +
  geom_point(mapping = aes(size = total), alpha=.5, show.legend = FALSE) +
  #stat_smooth(method="glm", method.args = list(family=binomial(link="logit")),lty=1, se=T,
              #formula = y ~ x) +
  scale_y_continuous(limits = c(0,1)) +
  ggtitle(paste("year", y17$year)) +
  scale_x_continuous(limits = c(0,70)) +
  #annotate("text", x=60, y=0.80, label = paste("b = ", round(res17$estimate[2], 3), "(±",(round(res17$std.error[2], 3)),")")) +
  #annotate("text", x=60, y=0.75, label = paste("p = ", round(res17$p.value[2], 5))) +
  labs(x = "Distance (km)", y = "Infection frequency", color = "sh")


y18 <- subset(t2.long.na.na.t, year == "X2018")
y18logit <- glm(infectionFreq ~ dist, family=binomial(link="logit"), data = y18)
summary(y18logit)
res18 <- broom::tidy(y18logit)
#plot logistic regression curve
gg18 <- ggplot(data=y18, aes(x=dist, y=infectionFreq)) + 
  theme_classic() +
  geom_point(mapping = aes(size = total), alpha=.5, show.legend = FALSE) +
  #stat_smooth(method="glm", method.args = list(family=binomial(link="logit")),lty=1, se=T,
              #formula = y ~ x) +
  scale_y_continuous(limits = c(0,1)) +
  ggtitle(paste("year", y18$year)) +
  scale_x_continuous(limits = c(0,70)) +
  #annotate("text", x=60, y=0.80, label = paste("b = ", round(res18$estimate[2], 3), "(±",(round(res18$std.error[2], 3)),")")) +
  #annotate("text", x=60, y=0.75, label = paste("p* = ", round(res18$p.value[2], 5)), color="red") +
  labs(x = "Distance (km)", y = "Infection frequency", color = "sh")


y19 <- subset(t2.long.na.na.t, year == "X2019")
y19logit <- glm(infectionFreq ~ dist, family=binomial(link="logit"), data = y19)
summary(y19logit)
res19 <- broom::tidy(y19logit)
#plot logistic regression curve
gg19 <- ggplot(data=y19, aes(x=dist, y=infectionFreq)) + 
  theme_classic() +
  geom_point(mapping = aes(size = total), alpha=.5, show.legend = FALSE) +
  #stat_smooth(method="glm", method.args = list(family=binomial(link="logit")),lty=1, se=T,
              #formula = y ~ x) +
  scale_y_continuous(limits = c(0,1)) +
  scale_x_continuous(limits = c(0,70)) +
  ggtitle(paste("year", y19$year)) +
  #annotate("text", x=60, y=0.80, label = paste("b = ", round(res19$estimate[2], 3), "(±",(round(res19$std.error[2], 3)),")")) +
  #annotate("text", x=60, y=0.75, label = paste("p = ", round(res19$p.value[2], 5))) +
  labs(x = "Distance (km)", y = "Infection frequency", color = "sh")


y20 <- subset(t2.long.na.na.t, year == "X2020")
y20logit <- glm(infectionFreq ~ dist, family=binomial(link="logit"), data = y20)
summary(y20logit)
res20 <- broom::tidy(y20logit)
#plot logistic regression curve
gg20 <- ggplot(data=y20, aes(x=dist, y=infectionFreq)) + 
  theme_classic() +
  geom_point(mapping = aes(size = total), alpha=.5, show.legend = FALSE) +
  #stat_smooth(method="glm", method.args = list(family=binomial(link="logit")),lty=1, se=T,
              #formula = y ~ x) +
  scale_y_continuous(limits = c(0,1)) +
  scale_x_continuous(limits = c(0,70)) +
  ggtitle(paste("year", y20$year)) +
  #annotate("text", x=60, y=0.80, label = paste("b = ", round(res20$estimate[2], 3), "(±",(round(res20$std.error[2], 3)),")")) +
  #annotate("text", x=60, y=0.75, label = paste("p* = ", round(res20$p.value[2], 5)), color="red") +
  labs(x = "Distance (km)", y = "Infection frequency", color = "sh")

logreg <- ggarrange(gg15, gg16, gg17, gg18, gg19, gg20, ncol = 2, nrow = 3)
ggsave("/Volumes/LaCie/the_cline/ifreq_t2years.png", plot = logreg, width = 10, height = 10, dpi = 300, units = "in", device='png')

# predictions from the model
pred = predict(mylogit, y15, se.fit = TRUE)
df_pred = data.frame(y15, pred=pred$fit, se=pred$se)
ggplot(df_pred, aes(x=dist,y=infectionFreq, group=year, col=year)) +
  geom_point(alpha = .5) + 
  geom_ribbon(aes(ymin=pred-1.96*se,ymax=pred+1.96*se, fill=year),alpha=0.1) #+
  #scale_y_continuous(limits = c(0,1))



#---------------investigation of wCer2 0.2<p<0.8 accros years with the largest sampling size----------#
t2.long.na = t2.long.latlon.pasohl23hayani.na
#confidence interval calculation see bellow *after CI calculation, infection frequency is subtracted by the low and high CI because this function just adds or substracts from the dot (for graphing purposes)
for(i in 1:nrow(t2.long.na)){
  
  #Estimate infection freq. and 95% binomial confidence intervals
  #t2.long.na.na.t[i,13] <- binconf(x=as.numeric(t2.long.na.na.t[i,8]), n=as.numeric(t2.long.na.na.t[i,8] + t2.long.na.na.t[i,9]))[1] #infection freq.
  t2.long.na[i,11] <- binconf(x=as.numeric(t2.long.na[i,8]), n=as.numeric(t2.long.na[i,8] + t2.long.na[i,9]))[2]
  t2.long.na[i,12] <- binconf(x=as.numeric(t2.long.na[i,8]), n=as.numeric(t2.long.na[i,8] + t2.long.na[i,9]))[3]
}
colnames(t2.long.na) <- c("id", "pop", "latitude", "longitude", 
                               "year", "infection", "total", "inf", "uninf", "infectionFreq", "infectionFreq_lowerCI", "infectionFreq_upperCI")

dist2 <- c()
for(i in 1:nrow(t2.long.na)) {
  
  row <- t2.long.na[i,] # extract row
  
  # approximate radius of earth in km
  R = 6373.0
  lat1 = rad(48.79926)
  lon1 = rad(16.49265)
  lat2 = rad(row$latitude)
  lon2 = rad(row$longitude)
  
  dlon = lon2 - lon1
  dlat = lat2 - lat1
  
  a = sin(dlat / 2)**2 + cos(lat1) * cos(lat2) * sin(dlon / 2)**2
  c = 2 * atan2(sqrt(a), sqrt(1 - a))
  # distance in meters
  distance = (R * c)*1000
  #print(distance)
  dist2  <- c(dist2, distance) # save as a vector
}
dist2

# add distance as a column
t2.long.na$dist <- dist2/1000

#p0208 <- subset(t2.long.na, infectionFreq >= 0.2 & infectionFreq <= 0.8)

intrapop <- c("Pasohlavky 3", 
              "Ivan 2", 
              "Velke Nemcice 1", 
              "NosislavZidlochovice 2",
              "Zidlochovice",
              "Rajhrad", 
              "Obchodní Centrum Futurum", 
              "Brno  Spilberk Castle")


t2.long.na.intra <- subset(t2.long.na, pop %in% intrapop)
t2.long.na.intra$pop <- factor(t2.long.na.intra$pop, levels=c("Pasohlavky 3", 
                                                              "Ivan 2", 
                                                              "Velke Nemcice 1", 
                                                              "NosislavZidlochovice 2",
                                                              "Zidlochovice",
                                                              "Rajhrad", 
                                                              "Obchodní Centrum Futurum", 
                                                              "Brno  Spilberk Castle"))

instapopggy <- ggplot(t2.long.na.intra, aes(x=dist, y=infectionFreq, col = pop, group = pop)) + 
  geom_errorbar(aes(ymin=infectionFreq_lowerCI, 
                    ymax=infectionFreq_upperCI), 
                width=.2, size=1.5) +
  geom_point(size=2) +
  #geom_line() +
  theme_bw() +
  theme(axis.title.y = element_text(vjust= 1.8),
        axis.title.x = element_text(vjust= 0.5),
        axis.text.x = element_text(angle = 90, vjust= 0.5),
        axis.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  #scale_y_continuous(limits = c(0, 1.10), breaks=c(0,0.25,0.5,0.75,1.0)) +
  labs(y = "Infection Frequency", x = "Distance (km)") +
  facet_wrap(~year, nrow=1) +
  ggtitle("intrapopulation wCre2 freq")
ggsave("/Volumes/LaCie/the_cline/intrapopT2y.png", plot = instapopggy, width = 15, height = 10, dpi = 300, units = "in", device='png')

# ---location NosislavZidlochovice 2-----
#nonzid2 <- subset(t2.long.na, pop == "NosislavZidlochovice 2")

#-- insta dist----#
trapData <- tibble()
splitData <- split(t2.long.na.intra, t2.long.na.intra$pop)

for(i in 1:length(splitData)){
  trapData[i,1] <- splitData[[i]]$id[1] #id
  trapData[i,2] <- splitData[[i]]$pop[1] #pop
  trapData[i,3] <- splitData[[i]]$latitude[1] # lat
  trapData[i,4] <- splitData[[i]]$longitude[1] #lon
  trapData[i,5] <- splitData[[i]]$year[1] # year
  trapData[i,6] <- splitData[[i]]$infection[1] # infection
  trapData[i,7] <- splitData[[i]]$total[1] # infection
  trapData[i,8] <- splitData[[i]]$inf[1] # number of infected individuals
  trapData[i,9] <- splitData[[i]]$uninf[1]
  trapData[i,10] <- splitData[[i]]$dist[1] # number of uninfected individuals
  #trapData[i,5] <- length(which(splitData[[i]]$Infected == "I" )) #n infected families
  #trapData[i,6] <- length(which(splitData[[i]]$Infected == "U" )) #n uninfected families
  
  #Estimate infection freq. and 95% binomial confidence intervals
  trapData[i,11] <- binconf(x=as.numeric(trapData[i,8]), n=as.numeric(trapData[i,8] + trapData[i,9]))[1] #infection freq.
  trapData[i,12] <- binconf(x=as.numeric(trapData[i,8]), n=as.numeric(trapData[i,8] + trapData[i,9]))[2]
  trapData[i,13] <- binconf(x=as.numeric(trapData[i,8]), n=as.numeric(trapData[i,8] + trapData[i,9]))[3]
  
} 
rm(i)
rm(splitData)

colnames(trapData) <- c("id", "pop", "lat", "lon", 
                        "year", "infection", "total", "inf", "uninf", "dist", "infFreq", "infectionFreq_lowerCI", "infectionFreq_upperCI")

# round the decimals
trapData[,9:11] <- round(trapData[,9:11], digits=3)
# re-order the data frame based on spatial location of the cline
trapData$pop <- factor(trapData$pop, levels=c("Pasohlavky 3", 
                                                              "Ivan 2", 
                                                              "Velke Nemcice 1", 
                                                              "NosislavZidlochovice 2",
                                                              "Zidlochovice",
                                                              "Rajhrad", 
                                                              "Obchodní Centrum Futurum", 
                                                              "Brno  Spilberk Castle"))

#plot yakuba infection freq. data
instadist <- ggplot(trapData, aes(x=dist, y=infFreq, group = pop, col=pop)) + 
  geom_errorbar(aes(ymin=infectionFreq_lowerCI, 
                    ymax=infectionFreq_upperCI), 
                width=.2, size=1.5) +
  geom_point(size=4) +
  theme_classic() +
  labs(y = "Infection Frequency", x = "Distance (km)") +
  ggtitle("")
ggsave("/Volumes/LaCie/the_cline/instadistT2.png", plot = instadist, width = 10, height = 10, dpi = 300, units = "in", device='png')



### > 3 years ####

intrapop3 <- c("Pasohlavky 3", 
              "Ivan 2",
              "Ivan-Vranovice 1",
              "Ivan-Vranovice 2",
              "Vranovice 2",
              "Prisnotice",
              "Velke Nemcice 1", 
              "NosislavZidlochovice 2",
              "Zidlochovice",
              "Holasice",
              "Rajhrad", 
              "Obchodní Centrum Futurum", 
              "Brno  Spilberk Castle",
              "Brno  Erbenova",
              "Soběšice 1",
              "Utechov")


t2.long.na.intra3 <- subset(t2.long.na, pop %in% intrapop3)
t2.long.na.intra3$pop <- factor(t2.long.na.intra3$pop, levels=c("Pasohlavky 3", 
                                                                "Ivan 2",
                                                                "Ivan-Vranovice 1",
                                                                "Ivan-Vranovice 2",
                                                                "Vranovice 2",
                                                                "Prisnotice",
                                                                "Velke Nemcice 1", 
                                                                "NosislavZidlochovice 2",
                                                                "Zidlochovice",
                                                                "Holasice",
                                                                "Rajhrad", 
                                                                "Obchodní Centrum Futurum", 
                                                                "Brno  Spilberk Castle",
                                                                "Brno  Erbenova",
                                                                "Soběšice 1",
                                                                "Utechov"))

instapopggy3 <- ggplot(t2.long.na.intra3, aes(x=dist, y=infectionFreq, col = pop, group = pop)) + 
  geom_errorbar(aes(ymin=infectionFreq_lowerCI, 
                    ymax=infectionFreq_upperCI), 
                width=.2, size=1.5) +
  geom_point(size=2) +
  #geom_line() +
  theme_bw() +
  theme(axis.title.y = element_text(vjust= 1.8),
        axis.title.x = element_text(vjust= 0.5),
        axis.text.x = element_text(angle = 90, vjust= 0.5),
        axis.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_x_continuous(limits = c(10, 60), breaks=c(10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60)) +
  labs(y = "Infection Frequency", x = "Distance (km)") +
  facet_wrap(~year, nrow=1) +
  ggtitle("intrapopulation wCre2 freq")
ggsave("/Volumes/LaCie/the_cline/intrapopT23y.png", plot = instapopggy3, width = 20, height = 10, dpi = 300, units = "in", device='png')


ggplot(t2.long.na.intra3, aes(x=year, y=infectionFreq, col = pop, group = pop)) + 
  geom_errorbar(aes(ymin=infectionFreq_lowerCI, 
                    ymax=infectionFreq_upperCI), 
                width=.2, size=1.5) +
  geom_point(size=2) +
  #geom_line() +
  theme_bw() +
  theme(axis.title.y = element_text(vjust= 1.8),
        axis.title.x = element_text(vjust= 0.5),
        axis.text.x = element_text(angle = 90, vjust= 0.5),
        axis.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  #scale_y_continuous(limits = c(0, 1.10), breaks=c(0,0.25,0.5,0.75,1.0)) +
  labs(y = "Infection Frequency", x = "Distance (km)") +
  facet_wrap(~pop, nrow=3) +
  ggtitle("intrapopulation wCer2 freq")


### fisher extact test ####
pas3 <- subset(t2.long.na.intra3, pop == "Pasohlavky 3")
#2020 vs 2017
pas3inf <- matrix(c(pas3$inf[4], pas3$uninf[4],
                          pas3$inf[1], pas3$uninf[1]),
                        ncol=2)
fisher.test(pas3inf)

#2020 vs 2018
pas3inf <- matrix(c(pas3$inf[4], pas3$uninf[4],
                    pas3$inf[2], pas3$uninf[2]),
                  ncol=2)
fisher.test(pas3inf)

#2020 vs 2018
pas3inf <- matrix(c(pas3$inf[4], pas3$uninf[4],
                    pas3$inf[2], pas3$uninf[2]),
                  ncol=2)
fisher.test(pas3inf)

ivan2 <- subset(t2.long.na.intra3, pop == "Prisnotice ")


# ---location NosislavZidlochovice 2-----
#nonzid2 <- subset(t2.long.na, pop == "NosislavZidlochovice 2")

#-- insta dist----#
trapData <- tibble()
splitData <- split(t2.long.na.intra3, t2.long.na.intra3$pop)

for(i in 1:length(splitData)){
  trapData[i,1] <- splitData[[i]]$id[1] #id
  trapData[i,2] <- splitData[[i]]$pop[1] #pop
  trapData[i,3] <- splitData[[i]]$latitude[1] # lat
  trapData[i,4] <- splitData[[i]]$longitude[1] #lon
  trapData[i,5] <- splitData[[i]]$year[1] # year
  trapData[i,6] <- splitData[[i]]$infection[1] # infection
  trapData[i,7] <- splitData[[i]]$total[1] # infection
  trapData[i,8] <- splitData[[i]]$inf[1] # number of infected individuals
  trapData[i,9] <- splitData[[i]]$uninf[1]
  trapData[i,10] <- splitData[[i]]$dist[1] # number of uninfected individuals
  #trapData[i,5] <- length(which(splitData[[i]]$Infected == "I" )) #n infected families
  #trapData[i,6] <- length(which(splitData[[i]]$Infected == "U" )) #n uninfected families
  
  #Estimate infection freq. and 95% binomial confidence intervals
  trapData[i,11] <- binconf(x=as.numeric(trapData[i,8]), n=as.numeric(trapData[i,8] + trapData[i,9]))[1] #infection freq.
  trapData[i,12] <- binconf(x=as.numeric(trapData[i,8]), n=as.numeric(trapData[i,8] + trapData[i,9]))[2]
  trapData[i,13] <- binconf(x=as.numeric(trapData[i,8]), n=as.numeric(trapData[i,8] + trapData[i,9]))[3]
  
} 
rm(i)
rm(splitData)

colnames(trapData) <- c("id", "pop", "lat", "lon", 
                        "year", "infection", "total", "inf", "uninf", "dist", "infFreq", "infectionFreq_lowerCI", "infectionFreq_upperCI")

# round the decimals
trapData[,9:11] <- round(trapData[,9:11], digits=3)
# re-order the data frame based on spatial location of the cline
trapData$pop <- factor(trapData$pop, levels=c("Pasohlavky 3", 
                                              "Ivan 2",
                                              "Ivan-Vranovice 1",
                                              "Ivan-Vranovice 2",
                                              "Vranovice 2",
                                              "Prisnotice",
                                              "Velke Nemcice 1", 
                                              "NosislavZidlochovice 2",
                                              "Zidlochovice",
                                              "Holasice",
                                              "Rajhrad", 
                                              "Obchodní Centrum Futurum", 
                                              "Brno  Spilberk Castle",
                                              "Brno  Erbenova",
                                              "Soběšice 1",
                                              "Utechov"))

instadist <- ggplot(trapData, aes(x=dist, y=infFreq, group = pop, col=pop)) + 
  geom_errorbar(aes(ymin=infectionFreq_lowerCI, 
                    ymax=infectionFreq_upperCI), 
                width=.2, size=1.5) +
  geom_point(size=4) +
  theme_classic() +
  labs(y = "Infection Frequency", x = "Distance (km)") +
  ggtitle("")
ggsave("/Volumes/LaCie/the_cline/instadistT2.png", plot = instadist, width = 10, height = 10, dpi = 300, units = "in", device='png')





