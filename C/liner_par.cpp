#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <iostream>
#include <algorithm>
#include <valarray>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_LINE_LENGTH 1000
using namespace std;


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
    printf("         -l,  --version        show liner\n");

}

int flus(int m)
{
   int a;
   a = m + 1;
   return a;
}


int split_str(const char *psrc,int psrc_len,char sign,char **pdest,int pdest_len)
{
        if(NULL == psrc || NULL == pdest)
        {
                return 0;
        }

        int i = 0;
        int result = 0;
        int start = 0,end = 0;
        int len = 0;

        for(i = 0;i < psrc_len;i++)
        {
                if(psrc[i] == sign)
                {
                        end = i;
                        len = end - start;
                        pdest[result] = (char *)malloc(len+1);
                        if(pdest[result] != NULL)
                        {
                                memcpy(pdest[result],psrc+start,len);
                                pdest[result][len] = '\0';
                                result++;
                                start = end+1;
                        }else
                        {
                                return result;
                        }
                }

        }

        if(start != psrc_len-1)
        {
                len = psrc_len-start;
                pdest[result] = (char *)malloc(len+1);
                if(pdest[result] != NULL)
                {
                        memcpy(pdest[result],psrc+start,len);
                        pdest[result][len] = '\0';
                        result++;
                        start = end+1;
                }else
                {
                        return result;
                }
        }
        return result;
}


//释放内存空间
void free_space(char **pdest,int pdest_len)
{
	if(pdest == NULL)
		return;
	int i = 0;
	for(i = 0;i < pdest_len;i++)
	{
		if(pdest[i] != NULL)
		{
			free(pdest[i]);
		}
	}
	return;
}

template<class T>
int length(T& arr)
{
   return sizeof(arr) / sizeof(arr[0]);
}

// 类模板,函数可以输入array
template<class T>
int cal(T& x,T& y)
{
        valarray<double> data_x(x,length(x));
        valarray<double> data_y(y,length(y));

        float A = 0.0;
        float B = 0.0;
        float C = 0.0;
        float D = 0.0;
        A = (data_x*data_x).sum();
        B = data_x.sum();
        C = (data_x*data_y).sum();
        D = data_y.sum();
        float tmp = A*data_x.size() - B*B;
        float k,b;
        k = (C*data_x.size() - B*D) /tmp;
        b = (A*D -C*B) /tmp;
        cout << "y=" << k << "x+" << b << endl;
}



int main(int argc, char *argv[])
{
    int i = 0;
    char* aa;
    int bb;
    double c;

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
             c = strtod(argv[i+1],NULL);
             bb = flus(c);
             printf("Used interface: %d\n", bb);
         }
         if((strcmp(argv[i], "-l")==0) || (strcmp(argv[i], "--liner")==0))
         {
            double x[10];
            double y[10];
            double z[10];

            int count=0;
            char *dst[100];

            FILE *fi;
            fi = fopen(argv[i+1], "r");

            char line[MAX_LINE_LENGTH];
            while ( fgets(line, MAX_LINE_LENGTH, fi) != NULL ){
                split_str(line, strlen(line),',',dst,sizeof(dst)/sizeof(char *));
                count +=1;

            if(count % 2 == 0){
              for(int i=0;i<strlen(dst[i]);i++){
                double a = strtod(dst[i],NULL);
                y[i] = a;
              }
            }
            else{
              for(int i=0;i<strlen(dst[i]);i++){
                double a = strtod(dst[i],NULL);
                x[i] = a;
              } 
             }
            if(count % 2 == 0){
               cal(x,y);
            }
            //释放dst的内存空间
            //free_space(dst,sizeof(dst)/sizeof(char *));
          }
            // system("pause");
            return 0;
        }
    }
  return 1;
}
