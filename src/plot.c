#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../include/plot.h"

//I borrowed this from my supervisor (to help with consistancy between implimentaions)
//and made some small changes


float min(float a,float b)
{
	if(a>b) return(b);
    else return(a);
}

float max(float a,float b)
{
	if(a>b) return(a);
    else return(b);
}


int get_R(double valor)
{
    if (valor <= 0.25)
        return 0; 
    else if (valor <= 0.5)
        return(int)(1020*valor - 255);
    else return 255; 
}

int get_G(double valor)
{


    if (valor <= 0.25)
        return(int)(1020*valor)/255;
    else if (valor <= 0.5)
        return 1;
    else
        return(int)(-510*valor + 510)/255;

}

int get_B(double valor)
{
    if (valor <= 0.25)
        return (int)(-512*valor + 128);
    else return 0;

}

double normalize(double * data,int num_x,int num_y)
{
	double max=0.;
	double min=0.;
	for(int i = 0; i < num_x;i++) { 
		for(int j=0;j<num_y;j++){
			if(data[j*num_x+i]>max){
				max=data[j*num_x+i];
			}
			if(data[j*num_x+i]<min){
                                min=data[j*num_x+i];
                        }
		}
	}
	//printf("MAX IS %f \n",max);
	for(int i = 0; i < num_x;i++) { 
		for(int j=0;j<num_y;j++){
			data[j*num_x+i]/=max;
			if(data[j*num_x+i]>1.1){
				printf("WOW NOT WORING %f \n",data[j*num_y+i]);
				exit(0);
			}
		}
	}

	printf("MAX Is %f MIN is %f \n",max,min);
	return(max);

}


/*------ LOG SCALES -------------------------------- */


double scale_log(double val)
{
	return(log10(val+1.)/0.30102999566391);
	//return((log10(val+.001)+2.)/2.);
//	return(log1p(10.*val)/2.39);
	return(val);

}



void find_log_map(void)
{
	double a,b;
	double x1,x2,y1,y2;

	x1=0.001; x2=1;
	y1=0.001;y2=1;

	b=log10(y1/y2)/(x1-x2);
	a=y1/exp(b*x1);

	printf("A=%f  B=%f \n",a,b);


}

/*------------------------------------------------- */

#define PLOT_OX  40
#define PLOT_OY  40

void plot_array(double* data,int num_x,int num_y,const char* name)
{
        normalize(data,num_x,num_y);

        FILE *f = fopen(name, "wb");
        fprintf(f, "P6\n%i %i 255\n", num_x, num_y);
        for(int i = 0; i < num_x;i++) {
                for(int j=0;j<num_y;j++){
                        double value= scale_log( (double) (data[i*num_y+j]));

			int r = 255.* min(max(0, 1.5-fabs(1-4*(value-0.5))),1);
			int g = 255.* min(max(0, 1.5-fabs(1-4*(value-0.25))),1);
			int b = 255.* min(max(0, 1.5-fabs(1-4*value)),1);

                        fputc(r, f);   // 0 .. 255
                        fputc(g, f); // 0 .. 255
                        fputc(b, f);  // 0 .. 255
                }
        }
        fclose(f);


}


