#include <stdio.h>
#include <stdlib.h>

#define MAX_LINE_LENGTH 1000

float countchar(char *str,char a){
    float n = 0;
    int i = 0;
    while(*(str+i)!='\0'){
        if(*(str+i) == a)
            n++;
        i++;
    }
    return n;
}

int main(int argc, char *argv[])
{
        FILE *fi;
        FILE *fo;

        if ( argc == 2) {
            fi = fopen(argv[1], "r");
        } else if (argc == 3){
            fo = fopen(argv[2], "w");
        } else {
                exit(EXIT_FAILURE);
        }

        char line[MAX_LINE_LENGTH];
        char a = 'A';
        float n;
        int count=0;
        int m = 2;
        while ( fgets(line, MAX_LINE_LENGTH, fi) != NULL ){
	       count = count + 1;	
               if ( count = m ){
                   n = countchar(line,a)/150;
                   printf("Reads:%d Proportion:%.2f%%\n",count-1,n*100);
        }
               m = m + 4;
      }

        return 0;
}
