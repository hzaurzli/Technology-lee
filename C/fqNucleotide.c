#include <stdio.h>
#include <stdlib.h>

#define MAX_LINE_LENGTH 1000

float countword(char *str,char w){
    float n = 0;
    int i = 0;
    while(*(str+i)!='\0'){
        if(*(str+i) == w)
            n++;
        i++;
    }
    return n;
}

int main(int argc, char *argv[])
{
        FILE *fi;
        FILE *fo;
        char *w;

        fi = fopen(argv[1], "r");
        w = (char)*argv[2];

        char line[MAX_LINE_LENGTH];
        float n;
        int count=0;
        int m = 2;
        //char w = 'A';
        while ( fgets(line, MAX_LINE_LENGTH, fi) != NULL ){
	       count = count + 1;	
               if ( count = m ){
                   n = countword(line,w)/150;
                   printf("Reads:%d Proportion:%.2f%%\n",count-1,n*100);
        }
               m = m + 4;
      }

        return 0;
}
