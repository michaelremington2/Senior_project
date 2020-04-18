library(dplyr) 
library(ggplot2) 
library(GGally)
library(car)
library(plotly)
library(cowplot)
library(alr3) 
library(KernSmooth)
library(leaps) 
library(tree)
library(rpart)
library(randomForest)
library(ggdendro)
library(instaR)
library(mclust)

#park Lake, Jones pond,Doney lake, Gypsey Lake
Spotted_Frog_Data <- read.csv("~/Senior Project/SpottedFrog.csv")
Spotted_Frog_Data <- na.omit(Spotted_Frog_Data)

spotted<-Spotted_Frog_Data %>% filter(Bd.total.copy..<5*10^7)

Spot<-spotted %>% rename(Bd=Bd.total.copy..,
PC1=Principle.Coordinates.Axis..1,
PC2=Principle.Coordinates.Axis..2,
PC3=Principle.Coordinates.Axis..3,
PC4=Principle.Coordinates.Axis..4,
PC5=Principle.Coordinates.Axis..5)
#rename(mtcars, c("disp" = "displacement"))
############ EDA################

ggpairs(Spotted_Frog_Data, c(3:6))
ggpairs(Spot, c(3,7:11))

lake<- spotted %>% group_by(Lake) %>% summarise(Bdtotal=sum(Bd.total.copy..)) 
ggplot(lake,aes(x=Lake,y=Bdtotal))+geom_bar(stat='identity')+ggtitle('Bd Total')

ggplot(spotted,aes(x=Principle.Coordinates.Axis..1,y=Bd.total.copy..))+geom_point()+geom_smooth(method = 'lm')
spotted
#################### multiple regression 1 of principle coordinates########
######## principle coordinates are the representation of important genes sequences
m1<-lm(Bd.total.copy..~Principle.Coordinates.Axis..1+
         Principle.Coordinates.Axis..2+
         Principle.Coordinates.Axis..3+
         Principle.Coordinates.Axis..4+
         Principle.Coordinates.Axis..5,
       data=spotted)
summary(m1)
spotted

vif(m1)
stres<-rstandard(m1)
p1<- ggplot(spotted, aes(x = Principle.Coordinates.Axis..1,y = stres)) + geom_point()+
  ylab("Standardized Residuals")
p2<- ggplot(spotted, aes(x = Principle.Coordinates.Axis..2,y = stres)) +   geom_point()+
  ylab("Standardized Residuals")
p3<- ggplot(spotted, aes(x = Principle.Coordinates.Axis..3,y = stres)) +   geom_point()+
  ylab("Standardized Residuals")
p4<- ggplot(spotted, aes(x = Principle.Coordinates.Axis..4,y = stres)) +   geom_point()+
  ylab("Standardized Residuals")
p<- ggplot(spotted, aes(x = Principle.Coordinates.Axis..5,y = stres)) +   geom_point()+
  ylab("Standardized Residuals")
p5<- ggplot(spotted, aes(x = m1$fitted.values,y = stres)) +geom_point()+
  ylab("Standardized Residuals") + 
  xlab("Fitted Values")
p6<-ggplot(Spotted_Frog_Data, aes(x = lm.influence(m1)$hat, y = stres)) + geom_point()+
  ylab("Standardized Residuals") + 
  xlab("Influence")+
  geom_vline(xintercept = 2*(5+1)/nrow(spotted))+
  geom_hline(yintercept = 2)+
  geom_hline(yintercept = -2)
p1
p6
plot_grid(p1,p2,p3,p4,p,p5,p6)
influence<-lm.influence(m1)$hat
spotted<-spotted %>% mutate(influence=influence)


################ single linear models to avoid colinearity#######


#### ph
mph<-lm(Bd.total.copy..~
         Mean.pH,
       data=spotted)
summary(mph)

###### DO
mMean.DO.......mg.L.<-lm(Bd.total.copy..~
          Mean.DO.......mg.L.,
        data=spotted)
summary(mMean.DO.......mg.L.)

###### Mean.Water.Temp..C.
mMean.Water.Temp..C.<-lm(Bd.total.copy..~
                           Mean.Water.Temp..C.,
                         data=spotted)
summary(mMean.Water.Temp..C.)


# corelation water temp p value of .01 r^2 of .12


###### Mean.Total.K..mg.L.     
mMean.Total.K..mg.L.<-lm(Bd.total.copy..~
                           Mean.Total.K..mg.L.,
                         data=spotted)
summary(mMean.Total.K..mg.L.)




###### Mean.Total.Ca..mg.L.
mMean.Total.Ca..mg.L.<-lm(Bd.total.copy..~
                            Mean.Total.Ca..mg.L.,
                         data=spotted)
summary(mMean.Total.Ca..mg.L.)

###### Mean.Total.Na..mg.L.
mMean.Total.Na..mg.L.<-lm(Bd.total.copy..~
                            Mean.Total.Na..mg.L.,
                          data=spotted)
summary(mMean.Total.Na..mg.L.)
### almost

###### Mean.Total.Mg..mg.L.
mMean.Total.Mg..mg.L.<-lm(Bd.total.copy..~
                            Mean.Total.Mg..mg.L.,
                          data=spotted)
summary(mMean.Total.Mg..mg.L.)



###### Mean.Total.Fe..mg.L.
mMean.Total.Fe..mg.L.<-lm(Bd.total.copy..~
                            Mean.Total.Fe..mg.L.,
                          data=spotted)
summary(mMean.Total.Fe..mg.L.)

###### Mean.Total.S..mg.L.
mMean.Total.S..mg.L.<-lm(Bd.total.copy..~
                           Mean.Total.S..mg.L.,
                          data=spotted)
summary(mMean.Total.S..mg.L.)

###### Mean.Total.P..mg.L.
mMean.Total.P..mg.L.<-lm(Bd.total.copy..~
                           Mean.Total.P..mg.L.,
                         data=spotted)
summary(mMean.Total.P..mg.L.)


# p value of .01 and r^2 of .13

###### Mean.Total.N..mg.L.
mMean.Total.N..mg.L.<-lm(Bd.total.copy..~
                           Mean.Total.N..mg.L.,
                         data=spotted)
summary(mMean.Total.N..mg.L.)

# p value of .01 r^2 of .14

###### Mean.Total.N..mg.L.
mElevation<-lm(Bd.total.copy..~
                           Elevation..feet.,
                         data=lakenew)
summary(mElevation)
 
# 1.park lake
# 2. Gypsy lake
# 3. jones pond
# 4. Downey lake

############ Regression plot matrix##############

#### PC matrix
pc1<-ggplot(spotted,aes(x=Principle.Coordinates.Axis..1,y=Bd.total.copy..))+geom_point()+geom_smooth(method = 'lm')
pc2<-ggplot(spotted,aes(x=Principle.Coordinates.Axis..2,y=Bd.total.copy..))+geom_point()+geom_smooth(method = 'lm')
pc4<-ggplot(spotted,aes(x=Principle.Coordinates.Axis..4,y=Bd.total.copy..))+geom_point()+geom_smooth(method = 'lm')
plot_grid(pc1,pc2,pc4,nrow=1)


################################ full regression model###########
frog.regression<-lm(Bd.total.copy..~Principle.Coordinates.Axis..1+
                      Principle.Coordinates.Axis..2+
                      Mean.Water.Temp..C., data = spotted)
summary(frog.regression)
vif(frog.regression)
# Mean.Total.P..mg.L.+
#Mean.Water.Temp..C.Principle.Coordinates.Axis..4+
#Mean.Total.N..mg.L.+
 # Mean.Total.P..mg.L.+
b<-regsubsets(Bd.total.copy..~Principle.Coordinates.Axis..1+
                Principle.Coordinates.Axis..2+
                Principle.Coordinates.Axis..4+
                Mean.Total.N..mg.L.+
                Mean.Total.P..mg.L.+
                Mean.Water.Temp..C., data = spotted,method = 'exhaustive')
bs<-summary(b)
bs$adjr

#####adjr2############ regression tree #################

frog.tree<-tree(Bd.total.copy..~Principle.Coordinates.Axis..1+
                  Principle.Coordinates.Axis..2+
                  Mean.Water.Temp..C.,model =T, data = spotted)
frog.tree
plot(frog.tree, pch = 2)
text(frog.tree)
title("Regression Tree Predicting BD")

summary(frog.tree)

rf<-randomForest(Bd.total.copy..~Principle.Coordinates.Axis..1+
                   Principle.Coordinates.Axis..2+
                   Principle.Coordinates.Axis..4+
                   Mean.Total.N..mg.L.+
                   Mean.Total.P..mg.L.+
                   Mean.Total.N..mg.L.+
                   Mean.Total.Na..mg.L., data = spotted)

print(rf)
############################### clustering ##############################
frogcluster<-hclust(dist(spotted[,c(3,7)]),method="single")
froghc <- as.dendrogram(frogcluster)

ddata <- dendro_data(froghc, type = "rectangle")
Rows<-row.names(spotted)
spotted<-spotted %>% mutate(Rows = Rows)
tempdata<-ddata$labels %>% tbl_df()
tempdata<-left_join(tempdata,spotted, by = c('label'='Rows'))
ddata$labels[,3]<-tempdata$V3


#### Other correlations####
wt<-ggplot(spotted,aes(x=Mean.Water.Temp..C.,y=Bd.total.copy..))+geom_point()+geom_smooth(method = 'lm')+ggtitle('Water Tempurature')
pho<-ggplot(spotted,aes(x=Mean.Total.P..mg.L.,y=Bd.total.copy..))+geom_point()+geom_smooth(method = 'lm')+ggtitle('Phosphorus')
N<-ggplot(spotted,aes(x= Mean.Total.N..mg.L.,y=Bd.total.copy..))+geom_point()+geom_smooth(method = 'lm')+ggtitle('Nitrogen')
plot_grid(wt,pho,N)
####################################################
ggplot(lake,aes(x=Lake,y=sum))+geom_bar(stat='identity')+ggtitle('Bd total sum')

# conductivity, solar exposure, graphs vs lake, elevation 
#park lake
#### eig vec val
sp <- read.csv("~/Senior Project/matrixsetup.csv")


########### Clustering ###########
spottednew<-spotted %>% select(-Lake,-Identifier,-Elevation..feet.)
kmfrog<-kmeans(spottednew,3)
kmclusters<-as.factor(kmfrog$cluster)  

p2<-ggplot(spotted, aes(x = Principle.Coordinates.Axis..2, y =Bd.total.copy.., color = kmclusters))+ geom_point()

p1<-ggplot(spotted, aes(x = Principle.Coordinates.Axis..1, y =Bd.total.copy.., color = kmclusters))+ geom_point()

p3<-ggplot(spotted, aes(x = Mean.Water.Temp..C., y =Bd.total.copy.., color = kmclusters))+ geom_point()

plot_grid(p1,p2,p3,nrow=1)
