#ifndef __FFT_H__
#define __FFT_H__

#ifdef __cplusplus
extern "C"{
#endif


#ifndef M_PI
#define M_PI        3.14159265358979323846
#endif

typedef float fftfloat;

typedef struct
{
    fftfloat **costbl;
    fftfloat **negsintbl;
    unsigned short **reordertbl;
} FFT_Tables;


void fft_initialize	( FFT_Tables *fft_tables );
void fft_terminate	( FFT_Tables *fft_tables );
void rfft			( FFT_Tables *fft_tables, float *x, int logm );
void fft			( FFT_Tables *fft_tables, float *xr, float *xi, int logm );
void ffti			( FFT_Tables *fft_tables, float *xr, float *xi, int logm );

#ifdef __cplusplus
}
#endif

#endif
