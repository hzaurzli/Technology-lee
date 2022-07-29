diff_cor = function(x,z,w){
  apply(x,2,function(y){
    a = cor(y,z) - cor(y,w)
    return(a)
  })
} #z,w作为第二个函数内部的变量，在外部函数进行传参

a8_1 = diff_cor(dmso_8,dmso8_f_md,dmso8_m_md)
a8_2 = diff_cor(e2_8,dmso8_f_md,dmso8_m_md)
a8_3 = diff_cor(em_8,dmso8_f_md,dmso8_m_md)
a8_4 = diff_cor(em_e2_8,dmso8_f_md,dmso8_m_md)


