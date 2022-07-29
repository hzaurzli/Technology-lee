#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_LINE_LENGTH 1000

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

int main(int argc, char *argv[])
{
        FILE *fi;
        FILE *fo;
        char *dst[1000];
        int cnt,ccnt;
        char rna[] = "mRNA";
        char *ddst[1000];
            
   
        fi = fopen(argv[1], "r");

        char line[MAX_LINE_LENGTH];
       
        //char w = 'A';
        while ( fgets(line, MAX_LINE_LENGTH, fi) != NULL ){
                cnt = split_str(line, strlen(line),'\t',dst,sizeof(dst)/sizeof(char *));
                if(strcmp(dst[2], rna) == 0){
                   ccnt = split_str(dst[8], strlen(dst[8]),';',ddst,sizeof(ddst)/sizeof(char *));
                   printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s_mRNA;%s;\n",dst[0],dst[1],dst[2],dst[3],dst[4],dst[5],dst[6],dst[7],ddst[0],ddst[1]);
                }
                else{
                   printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s",dst[0],dst[1],dst[2],dst[3],dst[4],dst[5],dst[6],dst[7],dst[8]);
                }
		    //释放dst的内存空间
           free_space(dst,sizeof(dst)/sizeof(char *));
  }
      return 0;
}
