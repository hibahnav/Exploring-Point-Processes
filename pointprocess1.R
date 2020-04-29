########################### LIBRARIES ##################################
library(INLA)
library(inlabru)
library(spatstat)
library(Distance)
library(DSpat)
library(raster)
library(sp)
library(geosphere)
library(rgeos)
library(sf)

########################################################################

#Simulate homogenous poisson process (hpp)
set.seed(12994)
W <- owin(c(0,100), c(0,50)) #study window
hpp <- rpoispp(0.13, win=W)
hpp$n
plot(hpp)

#Let's make these spatial objects
hpp$locations<-cbind(hpp$x,hpp$y) 
hpp$locations <- data.frame(hpp$locations)
colnames(hpp$locations) <- c("X","Y")
hpp$points <- SpatialPoints(hpp$locations)
hpp$points = SpatialPointsDataFrame(hpp$points,as.data.frame(hpp$locations),proj4string = CRS(as.character("+proj=utm +zone=32N +datum=WGS84 +units=km +ellps=WGS84 +towgs84=0,0,0")))
proj4string(hpp$points) = CRS("+proj=utm +zone=32N +datum=WGS84 +units=km +ellps=WGS84 +towgs84=0,0,0")

#boundary - make it SpatialPolygon
hpp$bound <- spoly(data.frame(x=c(0, 100,100, 0), 
                              y=c(0, 0, 50, 50)))

########################## QUESTION 2 & 3 ###############################

##### THINNING~~~~~~~~~~

### STEP 1: create 10 transects

#10 VERTICAL EVENLY SPACED TRANSECTS
w = 2 #half width
transects = create.lines(W,nlines=10,width=w,angle=90)
transects <- transects[2:5]

#10 RANDOM TRANSECTS

randcoords <- function(index){
  
  x_0 <- runif(1,0,90)
  x_1 <- runif(1,x_0,100)
  #xs <- c(x_0,x_1)
  
  y_0 <- runif(1,0,40)
  y_1 <- runif(1,y_0,50)
  #ys <- c(y_0,y_1)
  return(cbind(x_0,x_1,y_0,y_1))
}

set.seed(49907)
rand.coords <- lapply(1:10,function(x) randcoords(x))
rand.df <-do.call('rbind',rand.coords)
rand.df <- as.data.frame(rand.df)

xLen <- apply(rand.df,1, function(x) abs(x[1] - x[2]))
xLen
##NOTE: if any value in xLen is < 4, then call those lines out and change x coordinates
#      to make line transect longer using code below:

#rand.df[6,] <- c(8.27050, 33.01201, 8.31372,17.59716)
#rand.df[10,] <- c(8.27050, 39.01201, 22.31372,42.59716)
#rand.df[2,] <- c(71.806333, 85.50646, 22.046188, 38.86261)
#rand.df[7,] <- c(35.787855, 55.12804, 6.265916, 11.18146)

##NOTE: There may be some overlapping lines....


### STEP 2: Convert to SpatialLinesDataFrame 
transects <- as.matrix(transects)  #vertical lines
rand.df <- as.matrix(rand.df)  #random lines

##Create Lines
createVertLines <- sapply(1:10,FUN=function(i){Lines(list(Line(matrix(transects[i,],ncol=2))),ID=as.character(i))})
createRandLines <- sapply(1:10,FUN=function(i){Lines(list(Line(matrix(rand.df[i,],ncol=2))),ID=as.character(i))})

##Create Spatial Lines
VertSPlines <- SpatialLines(createVertLines)
RandomSPlines <- SpatialLines(createRandLines)


##Create SpatialLinesDataFrame

#Vertical Lines
VertDF <- data.frame(id=c(1:10),filler_dat =rep(NA,times=10))
SPDF <- SpatialLinesDataFrame(VertSPlines,VertDF)
hpp$SPDF <- SPDF

#Random Lines
RandDF <- data.frame(id=c(1:10),filler_dat =rep(NA,times=10))
rSPDF <- SpatialLinesDataFrame(RandomSPlines,RandDF)
hpp$rSPDF <- rSPDF

##Assign Projection
#Vertical:
proj4string(hpp$SPDF) <- CRS("+proj=utm +zone=32N +datum=WGS84 +units=km +ellps=WGS84 +towgs84=0,0,0")
#Random:
proj4string(hpp$rSPDF) <- CRS("+proj=utm +zone=32N +datum=WGS84 +units=km +ellps=WGS84 +towgs84=0,0,0")


### Step 3: Calculate perpendicular distance from points to nearest transect
calc.dist=function(object.ppp, transmat){
  x.bar = object.ppp$x
  y.bar = object.ppp$y
  x0 = as.double(transmat$x_0)
  y0 = as.double(transmat$y_0)
  x1 = as.double(transmat$x_1)
  y1 = as.double(transmat$y_1)
  d = sum((x0-x1)^2 + (y0-y1)^2)
  u.bar = ((x.bar-x0)*(x1-x0) + (y.bar-y0)*(y1-y0))/d
  p.x = x0+u.bar*(x1-x0)
  p.y = y0+u.bar*(y1-y0)
  distVals = sqrt((x.bar-p.x)^2 + (y.bar-p.y)^2)
  return(list(distance=distVals))
}

#Vertical transects:
transects = as.data.frame(transects)
colnames(transects) <- c("x_0","x_1","y_0","y_1") 
VertDist <- lapply(1:10,function(x) calc.dist(hpp,transects[x,]))

#Random transects:
rand.df <- as.data.frame(rand.df)
RandDist <- lapply(1:10,function(x) calc.dist(hpp,rand.df[x,]))

####Now that we have the distance matrix for each transect, we need to compare each pair
####to extract the smallest distances b/c smaller distance mean nearer transects
####so this is a way to identify the distance of each pt from its nearest transect

#Vertical Transects:
for(i in 1:9){
  vdist1 <- VertDist[[i]]$distance[which(VertDist[[i]]$distance < VertDist[[i+1]]$distance)]
  vdist2 <- VertDist[[i+1]]$distance[which(VertDist[[i+1]]$distance < VertDist[[i]]$distance)]
  distMatVert <- c(vdist1,vdist2)}


#Random Transects:
for(i in 1:9){
  rdist1 <- RandDist[[i]]$distance[which(RandDist[[i]]$distance < RandDist[[i+1]]$distance)]
  rdist2 <- RandDist[[i+1]]$distance[which(RandDist[[i+1]]$distance < RandDist[[i]]$distance)]
  distMatRand <- c(rdist1,rdist2)}


#distMatVert and distMatRand are our final distance matrices

#save to hpp:
hpp$VertDistMat <- distMatVert
hpp$RandDistMat <- distMatRand


##### THINNING ####

#Half Normal probability function:
hn = function(distance, lsig){ 
  exp(-0.5*(distance/exp(lsig))^2)}


##Vertical Transects:
set.seed(12)
vert.seen <- rbinom(hpp$n,1,hn(hpp$VertDistMat,2))
length(vert.seen[which(vert.seen == 1)]) 

thinnedVert <- rbinom(hpp$n,1,vert.seen) * hpp$locations
length(thinnedVert$X[which(thinnedVert$X > 0.00001)]) 

hpp$thinnedVert <- thinnedVert

##### Random Transects:
set.seed(473)
random.seen <- rbinom(hpp$n,1,hn(hpp$RandDistMat,2))
length(random.seen[which(random.seen == 1)])  

thinnedRandom <- rbinom(hpp$n,1,random.seen) * hpp$locations
length(thinnedRandom$X[which(thinnedRandom$X > 0.00001)]) 

hpp$thinnedRandom <- thinnedRandom

#Plotting everything

#Vertical Transects - before and after thinning:
par(mfrow=c(2,2))
ggplot() + gg(hpp$SPDF) + geom_point(data=hpp$locations,aes(X,Y))
ggplot() + gg(hpp$SPDF) + geom_point(data=hpp$thinnedVert,aes(X,Y))

#Random Transects - before and after thinning:
par(mfrow=c(2,2))
ggplot() + gg(hpp$rSPDF) + geom_point(data=hpp$locations,aes(X,Y))
ggplot() + gg(hpp$rSPDF) + geom_point(data=hpp$thinnedRandom,aes(X,Y))


