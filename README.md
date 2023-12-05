# WolbachiaModeling

---
title: "Intrapopulation_wCer2_frequency_dynamics"
author: "Sonja Lecic"
date: "16/12/2022"
output: 
   cleanrmd::html_document_clean:
     theme: minicss
     toc: true
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = F)
```

```{css}
.columns {display: flex;}
h1 {color: purple;}
h2 {color: darkorange;}
```

# Author

author: Sonja Lecic 

emmail: slecic8@gmail.com

date: "16/12/2022"

# Investigation of intrapolulation frequency dynamics of the wCer2

```{R, echo=T, include=T, message=FALSE, warning=FALSE}
# load the libraries
library(ggplot2) # plotting
library(tibble) # data manipulation
library(Hmisc) # statistics
library(circular) # distances calculation
library(ggpubr)
library(geosphere)
```


I first start by calculatting the 95% binomial confidence intervals for the whole dataset. Assuming a binomial distribution I will estimate the exact 95% confidence intervals of p for each population/year using 'binconf' function in R package Hmisc.
```{R, echo=T, warning=FALSE}

t2 <- read.csv("/Volumes/LaCie/the_cline/T2transect.csv")
t2.long <- as_tibble(t2)
t2.long.na <- na.omit(t2.long)
for(i in 1:nrow(t2.long.na)){
  
#Estimate infection freq. and 95% binomial confidence intervals
t2.long.na[i,11] <- binconf(x=as.numeric(t2.long.na[i,8]), n=as.numeric(t2.long.na[i,8] + t2.long.na[i,9]))[2]
t2.long.na[i,12] <- binconf(x=as.numeric(t2.long.na[i,8]), n=as.numeric(t2.long.na[i,8] + t2.long.na[i,9]))[3]
}
colnames(t2.long.na) <- c("id", "pop", "latitude", "longitude", 
                               "year", "infection", "total", "inf", "uninf", "infectionFreq", "infectionFreq_lowerCI", "infectionFreq_upperCI")

```

I then cacluate the distance by taking the souther-most location (Alt prerau) as the starting point. I then calculate the distance in meters by taking the radius of the Earth into account and convert the distances to km. The distance between the southern-most and the northern-most location is ~ 69km.
```{R, echo=T, warning=FALSE}
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
round(dist2/1000,3)

# add distance as a column
t2.long.na$dist <- round(dist2/1000,3)
```

Now, we focus on populations with appreciable data; here I subset the populations that have wCer2 infection frequencies recorded for at least 3 sampling years. 16 out of the total of 43 populations pass this filtering step.
```{R, echo=T, warning=FALSE}

#table(t2.long.na$pop)
intrapop3 <- c("Pasohlavky 3", 
              "Ivan 2",
              "Ivan-Vranovice 1",
              "Ivan-Vranovice 2 ",
              "Vranovice 2",
              "Prisnotice ",
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
t2.long.na.intra3$pop <- factor(t2.long.na.intra3$pop, levels=intrapop3)
```

Because there are 16 populations I plot the frequency change (y axis) over time (x axis) for each population separately. I then inspect which populations follow the expected pattern of frequency change over time. The order of populations from left to right (also in the legend from top to bottom) reflects the position of each population along the cline in the South-North direction.
```{R, echo=T, warning=FALSE}
# pairwise comparisons of means

#prm <- compare_means(infectionFreq ~ year, group.by="pop", data = t2.long.na.intra3)
#prm <- as.data.frame(prm)

intrapop2 <- ggplot(t2.long.na.intra3, aes(x=year, y=infectionFreq, group = pop)) + 
  geom_errorbar(aes(ymin=infectionFreq_lowerCI, 
                    ymax=infectionFreq_upperCI), 
                width=.2, size=1) +
  geom_point(size=1) +
  geom_line() +
  theme_bw() +
  theme(axis.title.y = element_text(vjust= 1.8),
        axis.title.x = element_text(vjust= 0.5),
        axis.text.x = element_text(angle = 90, vjust= 0.5),
        axis.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  labs(y = "Infection Frequency", x = "Collection year") +
  facet_wrap(~pop, nrow=3) +
  ggtitle("Intrapopulation wCer2 infection frequency")
  #stat_compare_means(label.y = 50)
intrapop2
```

## Populations following the expected pattern of frequency change under strong CI

### Populations starting at high wCer2 frequency (> 0.2)
OK, so it looks like populations "Ivan 2", "Ivan-Vranovice 2 " and "Velke Nemcice 1" follow the expected pattern of frequency change by quickly going to fixation due to the high levels of CI in our system. However, populations "Pasohlavky 3", "Vranovice 2" show a less clear pattern of wCer2 frequency change under strong CI. "Pasohlavky 3" fluctuates between 1 and 0.8 over the four years of sampling, while "Ivan-Vranovice 2" seems to be stable at 0.9 frequency. Could this less clear pattern be explained with imperfect maternal transmission (mu>0)?
Next, I plot the change in frequency over time of these five populations together in one plot here.
```{R, echo=T, warning=FALSE}

popexp <- c("Pasohlavky 3", 
              "Ivan 2",
              "Ivan-Vranovice 1",
              "Ivan-Vranovice 2 ",
              "Velke Nemcice 1")
t2.long.na.intra3.exp <- subset(t2.long.na.intra3, pop %in% popexp)
intrapop1 <- ggplot(t2.long.na.intra3.exp, aes(x=year, y=infectionFreq, group = pop, col=pop)) + 
  geom_errorbar(aes(ymin=infectionFreq_lowerCI, 
                    ymax=infectionFreq_upperCI), 
                width=.01, size=0.5) +
  geom_point(size=2, aes(shape=pop)) +
  geom_line(aes(linetype=pop)) +
  theme_bw() +
  theme(axis.title.y = element_text(vjust= 1.8),
        axis.title.x = element_text(vjust= 0.5),
        axis.text.x = element_text(angle = 90, vjust= 0.5),
        axis.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(y = "Infection Frequency", x = "Collection year") +
  scale_color_brewer(palette = "Dark2") +
  #facet_wrap(~pop, nrow=3) +
  ggtitle("Intrapopulation wCer2 infection frequency")
intrapop1
```

### Populations starting at low frequency (p < 0.2)
Populations "Zidlochovice", "Holasice", "Rajhrad", "Brno  Spilberk Castle", "Brno  Erbenova", "Soběšice 1", "Utechov" with wCer2 infections frequency below 20% fluctuate at low frequencies across the colletions years.
I will plot these populations here.
```{R, echo=T, warning=FALSE}

popexp <- c("Zidlochovice", 
              "Holasice",
              "Rajhrad",
              "Brno  Spilberk Castle",
              "Brno  Erbenova",
            "Soběšice 1",
            "Utechov")
t2.long.na.intra3.exp <- subset(t2.long.na.intra3, pop %in% popexp)
intrapop1 <- ggplot(t2.long.na.intra3.exp, aes(x=year, y=infectionFreq, group = pop, col=pop)) + 
  geom_errorbar(aes(ymin=infectionFreq_lowerCI, 
                    ymax=infectionFreq_upperCI), 
                width=.01, size=0.5) +
  geom_point(size=2, aes(shape=pop)) +
  geom_line(aes(linetype=pop)) +
  theme_bw() +
  theme(axis.title.y = element_text(vjust= 1.8),
        axis.title.x = element_text(vjust= 0.5),
        axis.text.x = element_text(angle = 90, vjust= 0.5),
        axis.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(y = "Infection Frequency", x = "Collection year") +
  scale_color_brewer(palette = "Dark2") +
  #facet_wrap(~pop, nrow=3) +
  ggtitle("Intrapopulation wCer2 infection frequency")
intrapop1
```

## Populations following an unexpected pattern of frequency change under strong CI

### Populations starting at high wCer2 frequency (p > 0.2)
Populations "Vranovice 2", "Prisnotice" show an unexpected pattern of wCer2 frequency change under strong CI by going down in frequency from high initial frequencies. Could this again be explained by imperfect maternal transmission (mu >0)? Population "NosislavZidlochovice 2" also shows an unexpected pattern with the frequency in the first collection year going down in frequency from 0.25 (2017) to 0.04 in 2018 and slowly back up again in the subsequent collection years. Could this slow increase in frequency be explained by fintess effect rather than the effect of strong CI? And perhaps indicate that the unstable equlibirium frequncy of around 0.5 for the CI effects to pick up?
I will plot this populations here.
```{R, echo=T, warning=FALSE}

popexp <- c("Vranovice 2", 
              "Prisnotice ",
              "NosislavZidlochovice 2")
t2.long.na.intra3.exp <- subset(t2.long.na.intra3, pop %in% popexp)
intrapop1 <- ggplot(t2.long.na.intra3.exp, aes(x=year, y=infectionFreq, group = pop, col=pop)) + 
  geom_errorbar(aes(ymin=infectionFreq_lowerCI, 
                    ymax=infectionFreq_upperCI), 
                width=.01, size=0.5) +
  geom_point(size=2, aes(shape=pop)) +
  geom_line(aes(linetype=pop)) +
  theme_bw() +
  theme(axis.title.y = element_text(vjust= 1.8),
        axis.title.x = element_text(vjust= 0.5),
        axis.text.x = element_text(angle = 90, vjust= 0.5),
        axis.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(y = "Infection Frequency", x = "Collection year") +
  scale_color_brewer(palette = "Dark2") +
  #facet_wrap(~pop, nrow=3) +
  ggtitle("Intrapopulation wCer2 infection frequency")
intrapop1
```

### Populations starting at low wCer2 frequency (p < 0.2)
Populations "Obchondi Centrum Futurum" stays fixed at 0% infection accross all four collection years. Is this an isolated population?
I will plot this populations here.
```{R, echo=T, warning=FALSE}

popexp <- c("Obchodní Centrum Futurum")
t2.long.na.intra3.exp <- subset(t2.long.na.intra3, pop %in% popexp)
intrapop1 <- ggplot(t2.long.na.intra3.exp, aes(x=year, y=infectionFreq, group = pop, col=pop)) + 
  geom_errorbar(aes(ymin=infectionFreq_lowerCI, 
                    ymax=infectionFreq_upperCI), 
                width=.01, size=0.5) +
  geom_point(size=2, aes(shape=pop)) +
  geom_line(aes(linetype=pop)) +
  theme_bw() +
  theme(axis.title.y = element_text(vjust= 1.8),
        axis.title.x = element_text(vjust= 0.5),
        axis.text.x = element_text(angle = 90, vjust= 0.5),
        axis.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(y = "Infection Frequency", x = "Collection year") +
  scale_color_brewer(palette = "Dark2") +
  #facet_wrap(~pop, nrow=3) +
  ggtitle("Intrapopulation wCer2 infection frequency")
intrapop1
```

## Pairwise distances between populations
Here, I cacluate pairwise distances in km between focal populations using geospehre R package for easier estimation of distances between these focal populations.
```{R, echo=T, warning=FALSE}
# extract only populations longitude and latitude
un <- unique(t2.long.na.intra3[, c(2, 3, 4)])
un <- un[-9,]
lonlat <- unique(un[, c(2,3)])
# and create a matrix with pop column as row names
lonlat_mat <- matrix(c(un$longitude, un$latitude), nrow=16, ncol=2,  byrow = T, dimnames = list(intrapop3, c("longitude", "latitude")))

# compute distance matrix using library(geosphere) package
dist <- (distm(lonlat_mat, fun = distGeo))/1000 
# add ids 
colnames(dist) <- intrapop3
rownames(dist) <- intrapop3

dist
```
#### minimum distance between any two populations
```{R, echo=F, warning=FALSE}
min(dist)
```
#### maximum distance between any two populations
```{R, echo=F, warning=FALSE}
max(dist) 
```
