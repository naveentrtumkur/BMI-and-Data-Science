## variable assignment
a <- 'apple'
a

b <- 'banana'
b

c <- paste(b,a)
c

## matrix
x <- seq(1, 9)
m <- matrix(x, nrow = 3, ncol = 3)
dim(m)
class(m)

m[2, ] # select a row
m[, 1] # select a column
m[2, 3] # select an element
summary(m) # generate summaries of the matrix

## list
l <- list(x = 1:5, y = c('A', 'B'))
l
class(l)
l$x
l$y

## data frame
df <- data.frame(x = 1:3, y = c('A', 'B', 'C'))
df
class(df)
dim(df) 
nrow(df)
ncol(df)
df[2, ] # select a row
df[, 2] # select a column

## reading and writing data
write.csv(df, 'file.csv')

df1 <- read.csv('file.csv')
df2 <- read.csv('file.csv', header = TRUE)
df3 <- read.csv('file.csv')[, -1]
df4 <- read.csv('file.csv')[-1, ]

df
df1
df2
df3
df4



