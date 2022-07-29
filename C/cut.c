#include <stdio.h>
#include <string.h>


int main(int argc, char *argv[])
{
        printf("%d\n",argc);
        char a[221];
        sscanf(argv[1],"%s",a);
        printf("%s\n",a); 

}



