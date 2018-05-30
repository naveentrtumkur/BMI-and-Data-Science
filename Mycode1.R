#usage of for loop
for(i in 1:5){
  j <- i + 10
  print(j)
}

#usage of while loop
while(i < 5){
  print(i)
  i <- i + 1
}

#if statement
if(i>3)
{
  print("yes")
  
}else
{
  print("no")
}

#usage of functions
square <- function(x){
  squared <- x*x
  return(squared)
}
