#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../include/transform.h"

char testpath[256] = "../test/test_signals/";
double* mat;
double* oracle;
void argument_error(){
    printf("Usage: ./run -n1 0 -n2 10000 -f1 5 -f2 10000 -Ts 0.01 -transform S -file full_path_to_file \n");
    printf("n1,n2,f1,f2, transform optional, defaults: -n1 0 -n2 10000 -f1 5 -f2 10000 -transform S\n");
    printf("Please specify file and sampleing period \n");
}
    //defaults
int n_1 = 0;
int n_2 = 100000;
double f_1 = 0;
double f_2 = 1000;
char transform_type= 'S'; //S for stft, W for wavelet
double sampleing_period = 0.01; //idk

char* filename = NULL;


int set_input(int argc, char* argv[]){

    char* eptr;
//parse arguments
    for(int i = 1; i < argc; i++){
        if (i + 1 != argc){
            if (strcmp(argv[i], "-n1") == 0){
                n_1 = atoi(argv[i + 1]);
                i++;       // Move to the next flag
                continue;
                }
            else if (strcmp(argv[i], "-n2") == 0){
                n_2 = atoi(argv[i + 1]);
                i++;
                continue;
            }
            else if (strcmp(argv[i], "-f1") == 0){
                f_1 = strtod(argv[i + 1],&eptr);
                if (f_1 == 0){
                    printf("f1 cannot be 0, check if out of range\n");
                    return(1);
                }
                    
                i++;    // Move to the next flag
                continue;
                }
            else if (strcmp(argv[i], "-f2") == 0){
                f_2 = strtod(argv[i + 1],&eptr);
                if (f_2 == 0){
                 printf("f2 cannot be 0, check if out of range\n");
                 return(1);
                }
                i++;
                continue;
            }
            else if (strcmp(argv[i], "-Ts") == 0){
                sampleing_period = strtod(argv[i + 1], &eptr);
                if (sampleing_period  == 0){
                    printf("sampling period cannot be 0, check if out of range\n");
                    return(1);
                }
                i++;
                continue;
            }

            else if (strcmp(argv[i], "-transform") == 0){
                transform_type = *argv[i + 1];
                if (transform_type == 'S'){
                    i++;
                    continue;
                }
                else{
                    argument_error();
                    return(1);
                }
            }
            else if (strcmp(argv[i], "-file") == 0){
                filename = argv[i + 1];
                i++;
                continue;
            }
            else{
                argument_error();
                return(1);
            }
        }
        else{
            argument_error();
            return(1);
        }
    

    }
    if (filename != NULL) { //filename is not null
     printf("File: %s Trasnform: %c \n", filename, transform_type);
        return(0);
    }
    else{//filename is null
        printf("filename not specified. Use: -file full_path_to_file\n");
        return(1);
    }
}



int main(int argc, char* argv[]){
    set_input(argc, argv);
    //for testing
    //strcat(testpath, filename); // use if you just want to type the filename instead of the whole file path, but check the test path at the top of the document
  
    mat = comp_transform(testpath, n_1, n_2, f_1, f_2, transform_type, sampleing_period );

    return(0);
}
