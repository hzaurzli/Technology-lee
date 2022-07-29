#include <stdio.h>
#include <stdlib.h>

#define MAX_LINE_LENGTH 1000

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
	int count=0;

	while ( fgets(line, MAX_LINE_LENGTH, fi) != NULL ){
		if (count == 0){
		    line[0] = '>';
		}
		if (count < 2){
			if ( argc == 2){
			    fprintf(stdout, "%s", line);
			} else{
			    fprintf(fo, "%s", line);
			}
		}
		if (++count > 3){
			count = 0;
		}
	}

	fclose(fi);
	if (argc == 3) fclose(fo);

	return 0;

}

