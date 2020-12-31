double get_x(double n);
double* comp_transform(char sig_filename[256],int n_1, int n_2, double f_1, double f_2, char transform_type, double Ts);

//the following are only here for testing
//otherwise i think its bad because its not information hiding....
//but I have to test these somehow
double dominant_frequency(double* mat);

double* comp_stft(double * sig);

double* read_signal(char* infilename, double* signal, int n_1, int n_2);

int zero_pad_signal(double* sig, double* padded_sig);


//double* comp_psudo_oracle(char sig_filename[256], int n_1, int n_2, double f_1, double f_2, char transform_type, double Ts);

//double compare(double* expected, double* calc);
