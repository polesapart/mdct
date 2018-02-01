#ifndef __MDCT_H__
#define __MDCT_H__
#include <string.h>


#ifdef __cplusplus
extern "C"{
#endif

#define M_PI        3.14159265358979323846
#define TWOPI       2*M_PI
#define FRAME_LEN   1024
#define BLOCK_LEN_LONG 1024
#define BLOCK_LEN_SHORT 128


void		CalculateKBDWindow	( float* win, float alpha, int length );
float	    Izero				( float x);
void		MDCT				( FFT_Tables *fft_tables, float *data, int N );
void		IMDCT				( FFT_Tables *fft_tables, float *data, int N );
void		MDCTKbdWindow		( FFT_Tables *fft_tables, float *p_overlap,float *p_in_data,float* p_kbd_window,float* p_out_mdct);
void		IMDCTKbdWindow		( FFT_Tables *fft_tables, float *p_overlap,float *p_in_data,float* p_kbd_window,float* p_out_data);

#ifdef __cplusplus
};
#endif

#endif