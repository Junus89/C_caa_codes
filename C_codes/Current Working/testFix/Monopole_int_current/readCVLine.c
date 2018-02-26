#include <stdio.h>
#include <string.h>
#include <stdbool.h>

#define datafile "ConstantV.txt"

//function prototypes
void GetNextLine (FILE *fp, char *line, int lineLen);

int main(int argc, char **argv){
		
	FILE *fp;
	char str[1000]; // character array to hold lines read from the file
	
	printf("Demonstrating file I/O in Standard C\n");
	printf("Opening file...\n");
	// Open the file for reading, checking to make sure it was successfully opened
	if((fp = fopen(datafile, "r"))==NULL)
	{
		// fopen() returns NULL if the file could not be found or if for some other reason the open failed.
		printf("Unable to open file %s\nProgram terminating...\n");
		return 0;
	}
	// file was successfully oopened so call the C function to read
	// all lines from the file.
	printf("Printing lines in file...\n\n");
	while(true)
	{
		GetNextLine(fp, str, 1000);
			if(strlen(str)==0)   // end of file reached
				break; 
		else
		{
			printf("%s\n",str);  // print the line
		}
	}// end while loop
	fclose(fp);
	printf("\n\ndone...\n\n");
	
	return 0;
}

void GetNextLine(FILE *fp, char *line, int lineLen)
{
	int done = false;
	while(!done)
	{
		if((fgets(line, lineLen, fp))== NULL) // read a line from the file
		{
			strcpy(line,""); // return an empty string to flag EOF
			return;
			if(strlen(line)>0)
				line[strlen(line)-1] = '\0';
			
			if(strlen(line)==0) // skip any blank lines
				continue;
			else if(line[0] == '#') // skip any comment lines
				continue;
			else done = true;  // Got a valid data line so return with this line
		}
	}
	
	
	
	
}