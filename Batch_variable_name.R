for(i in 1:3){
  assign(paste("p", i, sep=""), i)
  tmp <- get(paste("p", i, sep=""))
  print(tmp)
}
