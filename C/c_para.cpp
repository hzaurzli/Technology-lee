#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

int debug;

void show_version(char* name)
{
    printf("%s by Late Lee, version: 1.0\n", name);
}

void usage(char* name)
{
    show_version(name);

    printf("         -h,  --help           short help\n");
    printf("         -v,  --version        show version\n");
    printf("         -c,  --version        show char\n");
    printf("         -p,  --version        show plus\n");
}

int plus(int m)
{
   int a;
   a = m + 1;
   return a;
}

int main(int argc, char *argv[])
{
    int i = 0;
    char* aa;
    int bb;
    int c;

    /* early check for debug and config parameter */

    for (i = 1; i < argc; i++)
    {
         if ((strcmp(argv[i], "-h")==0) || (strcmp(argv[i], "--help")==0))
         {
             usage(argv[i]);
         }
         if ((strcmp(argv[i], "-c")==0) || (strcmp(argv[i], "--char")==0))
         {
             aa=argv[i+1];
             printf("Used interface: %s\n", aa);
         }
         if ((strcmp(argv[i], "-p")==0) || (strcmp(argv[i], "--plus")==0))
         {
             sscanf(argv[i+1],"%d",&c);
             bb = plus(c);
             printf("Used interface: %d\n", bb);
         }
    }
  return 1;
}
