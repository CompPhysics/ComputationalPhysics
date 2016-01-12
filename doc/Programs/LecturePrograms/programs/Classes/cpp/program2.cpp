using namespace std;
#include <iostream>   

   int main(int argc, char *argv[])
   {
      FILE *in_file, *out_file;
 
      if( argc < 3)  {
          printf("The programs has the following structure :\n");
          printf("write in the name of the input and output files \n");
          exit(0);
      }
      in_file = fopen( argv[1], "r");    // returns pointer to the  input file
      if( in_file == NULL )  {           // NULL means that the file is missing
         printf("Can't find the input file %s\n", argv[1]);
         exit(0);
      }
      out_file = fopen( argv[2], "w");     // returns a pointer to the output file
      if( out_file == NULL )  {            // can't find the file
          printf("Can't find the output file%s\n", argv[2]);
          exit(0);
       }
       fclose(in_file);
       fclose(out_file);
       return 0;
  } 
