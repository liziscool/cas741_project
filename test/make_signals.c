
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define N 100000
#define SCALE 1

void write_signal(int freq){
    double x;
    char filename[64];
    sprintf(filename,"test_signals/%d.txt", freq);
    FILE *file;
    file = fopen(filename, "w");
    for (int n = 0; n < N; n++){
        x = SCALE * sin((2*M_PI*freq * n)/(100000));
        fprintf(file, "%f\n", x);
        
    }
}
void write_piecewise_signal(int freq1, int freq2){
    double x;
    char filename[64];
    sprintf(filename,"test_signals/%d_%d.txt", freq1, freq2);
    FILE *file;
    file = fopen(filename, "w");
    for (int n = 0; n < N/4; n++){
        x = sin((2*M_PI*freq1 * n)/(100000));
        fprintf(file, "%f\n", x);
        
    }
    for (int n = N/4; n < N; n++){
        x = sin((2*M_PI*freq2 * n)/(100000));
        fprintf(file, "%f\n", x);
        
    }
}

void write_combo_signal(int freq1, int freq2){
    double x;
    char filename[64];
    sprintf(filename,"test_signals/%d+%d.txt", freq1, freq2);
    FILE *file;
    file = fopen(filename, "w");
    for (int n = 0; n < N; n++){
        x = sin((2*M_PI*freq1 * n)/(100000)) + sin((2*M_PI*freq2 * n)/(100000));
        //printf("%f\n", x);
        fprintf(file, "%f\n", x);
        
    }
}

int main(void){
    
    write_signal(100);
    write_signal(500);
    write_signal(1000);
    write_signal(2000);
    write_signal(3000);
    write_signal(4000);
    write_signal(5000);
    write_signal(10000);
     
     
    write_signal(1000);
    write_signal(2000);
    write_signal(3000);
    write_signal(4000);
    write_signal(5000);
    write_signal(10000);
    write_signal(11000);
    write_signal(12000);
    write_signal(13000);
    write_signal(30000);
    write_signal(50000);
    write_signal(100000);
    
    write_piecewise_signal(3000, 10000);
    write_piecewise_signal(4000, 5000);
    write_piecewise_signal(3000, 5000);
    write_piecewise_signal(1000, 5000);
    write_piecewise_signal(4000, 8000);
    write_piecewise_signal(4000, 6000);
    write_piecewise_signal(10000, 5000);
    write_combo_signal(2000, 10000);
    write_combo_signal(3333, 10000);
    write_combo_signal(4000, 8000);
    return(0);
}
