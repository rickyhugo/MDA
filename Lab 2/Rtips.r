# R tips for HA2

# tips 1: importing data

data <- read.table(gzfile("zip.train.gz"))
data <- as.matrix(data)

## Choose my lucky number and save data into 'Four' 
Four <- data[which(data[,1]==4),2:257]

#------------------------------------------------------------------------------------#

# tips 2: displaying an image
## Color 
colors <- c('white','black'); cus_col <- colorRampPalette(colors=colors)
## image
z <- matrix(Four[8,256:1],16,16,byrow=T)[,16:1]
image(t(z),col=cus_col(256))

#------------------------------------------------------------------------------------#

# tips 3: multiple plots in R
par(mfrow = c(1,2)) # after running this line, R will put you next two plots in to a 1 row and 2 columns frame.
par(mfrow = c(1,1)) # change the frame back

#------------------------------------------------------------------------------------#

# tips 4: Randomly split a data set into training set and testing set.
n = dim(data)[1] # sample size
id = sample(1:n, floor(0.8*n)) # 'id' will contains the ids of all observations that will be included in training set.
tr_d = data[id, ] # training set
te_d = data[-id, ] # testing set