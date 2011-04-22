/*
  Dent Earl, dearl (a) soe ucsc edu
  June 2010 
  Pipe in a maf file, blocks that are preceeded
  with a line that reads '# Paralog=1' will be
  commented out. I.e., this

  # Paralog=1
  a score=0 pctid=0.0
  s evoHg19.chr20            795 1 + 73767698 A
  s evoHg19-evoPanTro2.chr20 794 1 + 73522843 G

  becomes this

  # # Paralog=1
  # a score=0 pctid=0.0
  # s evoHg19.chr20            795 1 + 73767698 A
  # s evoHg19-evoPanTro2.chr20 794 1 + 73522843 G

##################################################
Copyright (C) 2009-2011 by
Dent Earl (dearl@soe.ucsc.edu, dentearl@gmail.com)
Benedict Paten (benedict@soe.ucsc.edu, benedictpaten@gmail.com)
... and other members of the Reconstruction Team of David Haussler's
lab (BME Dept. UCSC)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
##################################################

*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define BUFF_SIZE 32768

void usage(void){
   fprintf(stderr, "USAGE: \n");
   exit(2);
}

void processLine(char *sptr, int *isParalog, int wasNewLine){
   if (strncmp(sptr, "# Paralog=1", 11) == 0){
      *isParalog = 1;
   }
   if(*isParalog){
      if((strcmp(sptr, "\n") == 0) && wasNewLine){
         *isParalog = 0;
         printf("\n");
      }else{
         if (wasNewLine)
            printf("# %s", sptr);
         else
            printf("%s", sptr);
      }
   }else{
      printf("%s", sptr);
   }
}

void checkLine(char *sptr, int *isNewLine){
   if (sptr[strlen(sptr)-1] == '\n' )
      *isNewLine = 1;
   else
      *isNewLine = 0;
}

int main(int argc, char **argv){
   char *sptr;
   char buff[BUFF_SIZE];
   int isParalog = 0;
   int isNewLine = 0;
   int wasNewLine = 0;
   sptr = fgets(buff, BUFF_SIZE, stdin);
   while (sptr != NULL){
      checkLine(sptr, &isNewLine);
      processLine(sptr, &isParalog, wasNewLine);
      sptr = fgets(buff, BUFF_SIZE, stdin);
      wasNewLine=isNewLine;
   }
   return 0;
}
