#include <stdio.h>
#include <stdlib.h>


int countchar(char *str,char a){
    int n=0;
    int i = 0;
    while(*(str+i)!='\0'){
        if(*(str+i) == a)
            n++;
        i++;
    }
    return n;
}

int main(){
    char str[20] = "AANNNTTAA";
    char a = 'A';
    int n;
    n = countchar(str,a);
    printf("字符：%c 在字符串：%s 中出现了%d次\n",a,str,n);
    
    return 0;
}
