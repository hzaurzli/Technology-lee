#apply function(x),x means each row's variables

data = data.frame(a = c(1,3,6,7),
                  b = c(4,1,2,3),
                  c = c(8,7,1,2),
                  d = c(5,4,3,1))

apply(data, 1, function(x){
  fisher.test(matrix(c(x[1],x[2],x[3],x[4]), nrow = 2, ncol = 2))$p.value
})
