#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>

#include "../include/plot.h"

//#define N 100000 //length of signal


double *in_signal;
double* padded_signal;

double *window;
double *mor_wavelet;
double *out_transform;

fftw_complex *fft_in;

int size_0pad = 100;
int hop_size =  100;
int time_div; //=50;
int fft_buffer_size;

int N;

typedef struct complex_t{
    double r;
    double i; //real and imaginary, apprarently 'im' is already reserved for something in c
} complex_t;

complex_t* init_complex_wavelet(int size, complex_t* wavelet){
    wavelet = (complex_t*)malloc(size*sizeof(complex_t));
    //wavelet->r = (double*)malloc(size*sizeof(double));
    //wavelet->i = (double*)malloc(size*sizeof(double));
    return(wavelet);
}

void free_complex_wavelet(complex_t* wavelet){
    free(wavelet);
}

double* init_mat(int X, int Y){
    double* mat = (double*)malloc(X*Y*sizeof(double));
    return(mat);
}

void free_mat(double* mat){
    free(mat);
}

double* read_signal(char* infilename, double* signal, int n_1, int n_2){
    FILE *infile = fopen(infilename, "r");
    double temp;
    if (! infile ){ //in_file == NULL
        printf("file does not exists or can't be read\n");
        return(NULL);
    }
    //skip over entries 0 to n_1
    for(int j = 0; j < n_1; j++){
        fscanf(infile, "%lf", &temp);
    }
    for(int j = 0; j < N; j++){
        if (fscanf(infile, "%lf", &temp) != 1 ){
            printf("error reading val = %f\n",temp);
            break;
        }
        signal[j]= temp;
        //printf("We just read value = %f \n",temp);
    }
    fclose(infile);
    return(signal);
}

int zero_pad_signal(double* sig, double* padded_sig){
    //is this a weird way to do this? - yes.
    //but at least its readable
    for(int j = 0; j < size_0pad; j++){
        padded_sig[j] = 0.0;
    }
    for(int j = size_0pad; j < N; j++){
        padded_sig[j] = in_signal[j];
    }
    for(int j = N; j < N+size_0pad; j++){
        padded_sig[j] = 0.0;
    }
    return(0);
}

void make_window(int size){
    //Hann Window: https://en.wikipedia.org/wiki/Window_function
    //Size size
    for(int m = 0; m < size; m++){
        window[m] = sin((M_PI * m)/size) * sin((M_PI * m)/size);
        //printf("window at %d = %f\n", m,window[m] );
    }
}
double* comp_stft(double * sig){
    
    fftw_complex *fft_buffer, *fft_out;
    fft_buffer = (fftw_complex*)fftw_malloc(fft_buffer_size*sizeof(fftw_complex));
    fft_out = (fftw_complex*)fftw_malloc(fft_buffer_size*sizeof(fftw_complex));
    fftw_plan plan = fftw_plan_dft_1d(fft_buffer_size, fft_buffer, fft_out, FFTW_FORWARD, FFTW_ESTIMATE);
    window = (double*)malloc(fft_buffer_size*sizeof(double));
    
    out_transform = init_mat(time_div,fft_buffer_size);
    make_window(fft_buffer_size);
    double re;
    double im;
    int window_start;
   
    for (int t = 0; t < time_div; t++ ){
    //for (int t = 0; t < N/hop_size; t++ ){
        //copy data into buffer
        window_start = t*hop_size;
        for(int i = 0; i < fft_buffer_size; i++){
            fft_buffer[i][0] = sig[i + window_start] * window[i];//real
            fft_buffer[i][1] = 0; //imaginary
        }
    
        fftw_execute(plan);
        
        for(int j = 0; j <fft_buffer_size; j++){
        re = fft_out[j][0] ;
        im = fft_out[j][1] ;
        out_transform[t+ time_div*j] = sqrt(re*re + im*im);
        
        if(out_transform[t+ time_div*j]>30){
            printf("t = %d j = %d re = %f im = %f  X = %f \n", t, j, re, im, out_transform[t+ time_div*j]);
        
            }
        }
   
    }
    free(window);
    fftw_destroy_plan(plan);
    return(out_transform);
}

complex_t* make_mor_wavelet(int size, complex_t* wavelet){
    //more information: https://paos.colorado.edu/research/wavelets/wavelet2.html
    /*
     //complicated version
    double norm_constant;
    int addmis_criterion;
    norm_constant = pow(M_PI,-0.25)* 1/sqrt(1 + exp(- sigma*sigma) - 2* exp(-0.75*sigma*sigma));
    addmis_criterion = exp(-0.5 *sigma *sigma);
    for(int m = 0; m < size; m++){
        mor_wavelet[m] = norm_constant * exp(-0.5*m*m)*(cos(sigma);
     */
    //simple version
    //double w_0 = 2*M_PI /(size); //honestly this is just a guess
    double w_0 = 2*M_PI;
    double c1 = pow(M_PI,-0.25);
    double eta;
    double real; double im;
    for(int m = 0; m < size; m++){
        eta = m/(size*1.0) * 8;
        //printf("m = %d eta = %f\n", m, eta);
        real =  c1 * cos(w_0*eta) * exp(-0.5*eta*eta);
        im = c1 * sin(w_0*eta) * exp(-0.5*eta*eta);
        //printf("here1.1\n");
        wavelet[m].r = real;
        //printf("here1.2\n");
        wavelet[m].i = im;
        //printf("m = %d r = %f i = %f\n", m, wavelet[m].r, wavelet[m].i);
    }
    return(wavelet);
    
}


double* comp_wavelet_transform(double * sig){
    /*
    fftw_complex *fft_buffer, *fft_out;
    fft_buffer = (fftw_complex*)fftw_malloc(fft_buffer_size*sizeof(fftw_complex));
    fft_out = (fftw_complex*)fftw_malloc(fft_buffer_size*sizeof(fftw_complex));
    fftw_plan plan = fftw_plan_dft_1d(fft_buffer_size, fft_buffer, fft_out, FFTW_FORWARD, FFTW_ESTIMATE);
    window = (double*)malloc(fft_buffer_size*sizeof(double));
    */
    out_transform = init_mat(time_div,fft_buffer_size);
    complex_t* wavelet = init_complex_wavelet(fft_buffer_size, wavelet);
    wavelet = make_mor_wavelet(fft_buffer_size, wavelet);
    double re;
    double im;
    int window_start;
    for (int t = 0; t < time_div; t++ ){
        window_start = t*hop_size;
        re = 0;
        im = 0;
        for(int s = 0; s < fft_buffer_size; s++){//for each scale
            
            //WE NEED TO SCALE THE WAVELET
            for(int i = 0; i < fft_buffer_size; i++){
                re += sig[i + window_start] * wavelet[i].r;//real
                im += sig[i + window_start] * wavelet[i].i;//imaginary
            }
        out_transform[t+ time_div*s] = sqrt(re*re + im*im);
        if(out_transform[t+ time_div*s]>30){
            printf("t = %d s = %d re = %f im = %f  X = %f \n", t, s, re, im, out_transform[t+ time_div*s]);
            }
        }
   
    }
    free_complex_wavelet(wavelet);
    return(out_transform);
}

double dominant_frequency(double* mat){
    int max = 0;
    double freq = -1;
    for (int t = 0; t < time_div; t++ ){
        for(int j = 0; j <fft_buffer_size; j++){
            if (mat[t+ time_div*j] >= max){
                max = mat[t+ time_div*j];
                freq = (fft_buffer_size-j)*250;//(fft_buffer_size*1.0*0.01*2*M_PI);
            }
        }
    }
    return(freq);
}

double* comp_transform(char sig_filename[256], int n_1, int n_2, double f_1, double f_2, char transform_type, double Ts){
    
    
    N = n_2 - n_1;
    time_div =  N/hop_size;
    fft_buffer_size = 400; // freq_div;
    

    in_signal = (double*)malloc(N*sizeof(double));
    padded_signal = (double*)malloc((N + 2*size_0pad)*sizeof(double));
    //fft_in = (fftw_complex*)fftw_malloc(N*sizeof(fftw_complex) + 2*size_0pad*sizeof(fftw_complex));

    double* out_mat;
    
    read_signal(sig_filename, in_signal, n_1, n_2);
    zero_pad_signal(in_signal,padded_signal);

    //out_mat = comp_stft(padded_signal);
    out_mat = comp_wavelet_transform(padded_signal);
    printf("Matrix Computed\n");
    double df = dominant_frequency(out_mat);
    printf("dominant frequency detected: %f Hz\n", df);
    plot_array(out_mat,time_div,fft_buffer_size/2,"out.png");
    free(in_signal);
    free(padded_signal);
    return(out_mat);
 }

/*-------------------------------------------------------*/
//PSEDO ORACLE
//fix these later

/*
#include "psudo_oracle/wavelet.h"


double* comp_psudo_oracle(char sig_filename[256], int n_1, int n_2, double f_1, double f_2, char transform_type, double Ts){

    fft_buffer_size = 400;

    N = n_2 - n_1;
    time_div =  N/hop_size;
    
    double* data = (double*)malloc(N*sizeof(double));
    double* frequencys = (double*)malloc(N*sizeof(double));
    double* scales = (double*)malloc(N*sizeof(double));
    double* res = init_mat(N, fft_buffer_size);
    int n = fft_buffer_size; //?

    read_signal(sig_filename, in_signal, n_1, n_2);
    printf("signals read\n");
    //generate_scales(f_1,f_2 , N,frequencys,scales,W_MORLET);
    generate_scales(0,fft_buffer_size ,N,frequencys,scales,W_MORLET);
    printf("scales generated\n");
    res = ctm_trans(data, scales, 1/Ts,n,N,res,&mor_file_coeff);
    printf("ctm transform\n");
    plot_array(res,N,fft_buffer_size/2,"oracle_out.png");
    free(data);free(frequencys); free(scales);
    return(res);
}

void normalize_01(double* mat, int size){
    double max = -100;
    double min = 100;
    for (int i = 0; i < size; i++){
        if (mat[i] > max) max = mat[i];
        if (mat[i] < min) min = mat[i];
    }
    //printf("max: %.3f  min: %.3f\n   mean: %.3f\n", max, min, mean);
    if (max == min || max == -1|| min == 1){
        printf("error normalizing matrix\n");
        //return;
    }
    double range = max - min;
        
    for (int i = 0; i < size; i++){
        mat[i] = (mat[i] - min)/(range);
    }
    return;
}

double compare(double* expected, double* calc){
    double error = 0;
    int size = time_div * fft_buffer_size;
    //normalize both mat between 0 and 1?
    normalize_01(calc, size);
    normalize_01(expected, size);
    //accuacy = (expected-calculated)/expected
    for (int t = 0; t < time_div; t++ ){
         for(int j = 0; j <fft_buffer_size; j++){
             error += fabs((expected[t+ time_div*j] - calc[t+ time_div*j])/expected[t+ time_div*j]);
         }
    }
    error = error/size;
    return(error);
}
 */
