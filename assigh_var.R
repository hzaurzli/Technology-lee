b = vector()
for (i in 1:10) {
  b[i] = paste0("a", i)
  assign(b[i],cbind(S@x.values[[i]],S@y.values[[i]]))
  
}

b1 = cbind(P@x.values[[1]],P@y.values[[1]])
uu = function(x){
  for (i in 1:10) {
    c = as.data.frame(x)
    c[,1] = as.character(c[,1])
    c = c[!duplicated(c[,1]),]
    c = c[!duplicated(c[,2]),]
  }
  return(c)
}

write.xlsx(a1, file="/home/lirz/111/file1.xlsx", sheetName="sheet1", row.names=FALSE)
write.xlsx(a2, file="/home/lirz/111/file1.xlsx", sheetName="sheet2",append = TRUE, row.names=FALSE)
write.xlsx(a3, file="/home/lirz/111/file1.xlsx", sheetName="sheet3",append = TRUE, row.names=FALSE)
write.xlsx(a4, file="/home/lirz/111/file1.xlsx", sheetName="sheet4",append = TRUE, row.names=FALSE)
write.xlsx(a5, file="/home/lirz/111/file1.xlsx", sheetName="sheet5",append = TRUE, row.names=FALSE)
write.xlsx(a6, file="/home/lirz/111/file1.xlsx", sheetName="sheet6",append = TRUE, row.names=FALSE)
write.xlsx(a7, file="/home/lirz/111/file1.xlsx", sheetName="sheet7",append = TRUE, row.names=FALSE)
write.xlsx(a8, file="/home/lirz/111/file1.xlsx", sheetName="sheet8",append = TRUE, row.names=FALSE)
write.xlsx(a9, file="/home/lirz/111/file1.xlsx", sheetName="sheet9",append = TRUE, row.names=FALSE)
write.xlsx(a10, file="/home/lirz/111/file1.xlsx", sheetName="sheet10",append = TRUE, row.names=FALSE)
