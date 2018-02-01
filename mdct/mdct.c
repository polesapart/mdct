#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "fft.h"
#include "mdct.h"


float Izero(float x)
{
	const float IzeroEPSILON = 1E-41;  /* Max error acceptable in Izero */
	float sum, u, halfx, temp;
	int n;

	sum = u = n = 1;
	halfx = x/2.0;
	do {
		temp = halfx/(float)n;
		n += 1;
		temp *= temp;
		u *= temp;
		sum += u;
	} while (u >= IzeroEPSILON*sum);

	return(sum);
}

/******************************************************************************************
  *  @brief      ����KBD��������ֵ  ֻ��Ҫ����ǰ�벿�� ��벿����ǰ�벿�ֶԳ�.
  *  @param[in]  win���洢������ֵ  alpha =4��һ�������  length=N/2=1024 ֻ����ǰ�벿��
  *  @param[out] win��Ϊ����Ĵ�����ֵ
  *  @Return:    
  *  @note       
  *  @date       2016.09.22
  *  @see        
  *****************************************************************************************/
void CalculateKBDWindow(float* win, float alpha, int length)
{
	int i;
	float IBeta;
	float tmp;
	float sum = 0.0;

	alpha *= M_PI;
	IBeta = 1.0/Izero(alpha);

	/* calculate lower half of Kaiser Bessel window */
	for(i=0; i<(length>>1); i++) {
		tmp = 4.0*(float)i/(float)length - 1.0;
		win[i] = Izero(alpha*sqrt(1.0-tmp*tmp))*IBeta;
		sum += win[i];
	}

	sum = 1.0/sum;
	tmp = 0.0;

	/* calculate lower half of window */
	for(i=0; i<(length>>1); i++) {
		tmp += win[i];
		win[i] = sqrt(tmp*sum);
	}
}

/******************************************************************************************
  *  @brief      MDCT�任.
  *  @param[in]  data��ʱ�����ݣ�N�����ã�N=2048
			     N:�任���� N=2048
  *  @param[out] data:mdct�任���Ƶ������ N�㣨ǰN/2�����N/2�����Գƣ�
  *  @Return:    
  *  @note       
  *  @date       2016.09.26
  *  @see        
  *****************************************************************************************/

void MDCT( FFT_Tables *fft_tables, float *data, int N )//N=2048
{
	float *xi, *xr;
	float tempr, tempi, c, s, cold, cfreq, sfreq; /* temps for pre and post twiddle */
	float freq = TWOPI / N;
	float cosfreq8, sinfreq8;
	int i, n;

	xi = (float*)malloc((N >> 2)*sizeof(float)); //xi N/4��
	xr = (float*)malloc((N >> 2)*sizeof(float)); //xr N/4��

	/* prepare for recurrence�����ƣ� relation in pre-twiddle */
	cfreq = cos (freq);
	sfreq = sin (freq);
	cosfreq8 = cos (freq * 0.125);
	sinfreq8 = sin (freq * 0.125);
	c = cosfreq8;
	s = sinfreq8;

	for (i = 0; i < (N >> 2); i++) {//i=0 1 2 .... 511
		/* calculate real and imaginary parts of g(n) or G(p)   ����ʵ�����鲿*/
		n = (N >> 1) - 1 - 2 * i;//n=1023 1021 ...  3 1
		
		//����N/4��512���㸴��FFT�任�� ʵ�����ݵĻ�ȡ
		if (i < (N >> 3)) //i=0 1 2 ....254 255
			//ǰ256��ʵ�����ݵĻ�ȡ
			//(N >> 2) + n=0.75*N-1-2*i=1535 1533 1531 ... 1027 1025 ;   N + (N >> 2) - 1 - n=0.75*N+2*i=1536 1538 ...2044  2046
			tempr = data [(N >> 2) + n] + data [N + (N >> 2) - 1 - n]; /* use second form of e(n) for n = N / 2 - 1 - 2i     */
		else  //i=256 257....510 511
			//��256��ʵ�����ݵĻ�ȡ
			//(N >> 2) + n=0.75*N-1-2*i=1023 1021 ....515 513 ;   (N >> 2) - 1 - n= -0.25*N+2*i= 0 2 4... 508 510
			tempr = data [(N >> 2) + n] - data [(N >> 2) - 1 - n]; /* use first form of e(n) for n = N / 2 - 1 - 2i */


		n = 2 * i;//n=0 2 4 ....1020 1022

		//����N/4��512���㸴��FFT�任�� �鲿���ݵĻ�ȡ
		if (i < (N >> 3)) //i=0 1 2 ....254 255
			//ǰ256���鲿���ݵĻ�ȡ
			//(N >> 2) + n=0.25*N+2*i=512 514  ... 1020 1022  ;   (N >> 2) - 1 - n= 0.25*N-1-2*i= 511 509 ....3 1
			tempi = data [(N >> 2) + n] - data [(N >> 2) - 1 - n]; /* use first form of e(n) for n=2i */
		else  //i=256 257....510 511
			//��256���鲿���ݵĻ�ȡ
			//(N >> 2) + n=0.25*N+2*i=1024 1026...1532 1534  ;    N + (N >> 2) - 1 - n=1.25*N-1-2*i=2047 2045 ...1539 1537
			tempi = data [(N >> 2) + n] + data [N + (N >> 2) - 1 - n]; /* use second form of e(n) for n=2i*/


		/* calculate pre-twiddled FFT input */
		xr[i] = tempr * c + tempi * s;
		xi[i] = tempi * c - tempr * s;

		/* use recurrence to prepare cosine and sine for next value of i */
		cold = c;
		c = c * cfreq - s * sfreq;
		s = s * cfreq + cold * sfreq;
	}


	/* Perform in-place complex FFT of length N/4 */
	switch (N) {
	case BLOCK_LEN_SHORT * 2: //128*2==256==N ��MDCT�任  ���Ӧ��fft�任Ϊ256/4=64�㣨N/4�㣩
		fft( fft_tables, xr, xi, 6);
		break;
	case BLOCK_LEN_LONG * 2: //1024*2==2048==N ��MDCT�任  ���Ӧ��fft�任Ϊ2048/4=512�㣨N/4�㣩
		fft( fft_tables, xr, xi, 9);
	}

	/* prepare for recurrence relations in post-twiddle */
	c = cosfreq8;
	s = sinfreq8;

	/* post-twiddle FFT output and then get output data */
	for (i = 0; i < (N >> 2); i++) {
		/* get post-twiddled FFT output  */
		tempr = 2. * (xr[i] * c + xi[i] * s);
		tempi = 2. * (xi[i] * c - xr[i] * s);

		//�任���ϵ�������ﻹ�Ǵ����data���棩 ����ԳƵģ�ʱ������N�㣬��Ƶ��ǰN/2�����N/2����Գƣ���ֻ��ҪǰN/2������ݱ������ȷmdct�任��
		/* fill in output values */
		data [2 * i] = -tempr;   /* first half even�� 0 2 4 6 ....1018 1020 1022 */  
		data [(N >> 1) - 1 - 2 * i] = tempi;  /* first half odd:  1023 1021 1019 .... 5 3 1 */


		data [(N >> 1) + 2 * i] = -tempi;  /* second half even :  1024 1026 1028 ... 2042 2044 2046*/
		data [N - 1 - 2 * i] = tempr;  /* second half odd :  2047 2045 2043 ... 1029 1027 1025*/

		/* use recurrence to prepare cosine and sine for next value of i */
		cold = c;
		c = c * cfreq - s * sfreq;
		s = s * cfreq + cold * sfreq;
	}

	if (xr) free(xr);
	if (xi) free(xi);

}
/******************************************************************************************
  *  @brief      IMDCT�任.
  *  @param[in]  data��Ƶ��MDCT���ݣ�N/2�����ã�
			     N:�任���� N=2048
  *  @param[out] data:imdct�任���ʱ������ N��
  *  @Return:    
  *  @note       
  *  @date       2016.09.26
  *  @see        
  *****************************************************************************************/
void IMDCT( FFT_Tables *fft_tables, float *data, int N)
{
	float *xi, *xr;
	float tempr, tempi, c, s, cold, cfreq, sfreq; /* temps for pre and post twiddle */
	float freq = 2.0 * M_PI / N;
	float fac, cosfreq8, sinfreq8;
	int i;

	xi = (float*)malloc((N >> 2)*sizeof(float));
	xr = (float*)malloc((N >> 2)*sizeof(float));

	/* Choosing to allocate 2/N factor to Inverse Xform! */
	//////////-----------------wjs --------------------------///
	//fac = 2. / N; /* remaining 2/N from 4/N IFFT factor */
	fac = 0.5;

	/* prepare for recurrence relation in pre-twiddle */
	cfreq = cos (freq);
	sfreq = sin (freq);
	cosfreq8 = cos (freq * 0.125);
	sinfreq8 = sin (freq * 0.125);
	c = cosfreq8;
	s = sinfreq8;

	for (i = 0; i < (N >> 2); i++) {
		/* calculate real and imaginary parts of g(n) or G(p) */
		tempr = -data[2 * i];
		tempi = data[(N >> 1) - 1 - 2 * i];

		/* calculate pre-twiddled FFT input */
		xr[i] = tempr * c - tempi * s;
		xi[i] = tempi * c + tempr * s;

		/* use recurrence to prepare cosine and sine for next value of i */
		cold = c;
		c = c * cfreq - s * sfreq;
		s = s * cfreq + cold * sfreq;
	}

	/* Perform in-place complex IFFT of length N/4 */
	switch (N) {
	case BLOCK_LEN_SHORT * 2:
		ffti( fft_tables, xr, xi, 6);
		break;
	case BLOCK_LEN_LONG * 2:
		ffti( fft_tables, xr, xi, 9);
	}

	/* prepare for recurrence relations in post-twiddle */
	c = cosfreq8;
	s = sinfreq8;

	/* post-twiddle FFT output and then get output data */
	for (i = 0; i < (N >> 2); i++) {

		/* get post-twiddled FFT output  */

		tempr = fac*(xr[i] * c - xi[i] * s);
		tempi = fac*(xi[i] * c + xr[i] * s);

		/* fill in output values */
		data [(N >> 1) + (N >> 2) - 1 - 2 * i] = tempr;
		if (i < (N >> 3))
			data [(N >> 1) + (N >> 2) + 2 * i] = tempr;
		else
			data [2 * i - (N >> 2)] = -tempr;

		data [(N >> 2) + 2 * i] = tempi;
		if (i < (N >> 3))
			data [(N >> 2) - 1 - 2 * i] = -tempi;
		else
			data [(N >> 2) + N - 1 - 2*i] = tempi;

		/* use recurrence to prepare cosine and sine for next value of i */
		cold = c;
		c = c * cfreq - s * sfreq;
		s = s * cfreq + cold * sfreq;
	}

	if (xr) free(xr);
	if (xi) free(xi);


}



/******************************************************************************************
  *  @brief      ��50%�ص�KBD����MDCT�任.
  *  @param[in]  p_in_data    1024 ��ǰ֡ʱ������ ��Ϊ��ǰ�任��ĺ�1024������ 
				 p_overlap    1024 ��һ֡ʱ������ ��Ϊ��ǰ�任���ǰ1024������ 
			     p_kbd_window 1024 kbd����ǰ1024����  ��1024������ǰ1024����Գ� 	
				 
  *  @param[out] p_out_mdct   2048 mdctƵ������ ����Գ� ��IMDCT�任��ֻ��ǰ1024������
  *  @Return:    
  *  @note       
  *  @date       2016.09.23
  *  @see        
  *****************************************************************************************/

void MDCTKbdWindow(FFT_Tables *fft_tables,float *p_overlap,float *p_in_data,float* p_kbd_window,float* p_out_mdct)
{

	int i;
	//float *trans_buf=(float *)malloc(2*BLOCK_LEN_LONG*sizeof(float));
	memcpy(p_out_mdct,p_overlap,BLOCK_LEN_LONG*sizeof(float));
	memcpy(p_out_mdct+BLOCK_LEN_LONG,p_in_data,BLOCK_LEN_LONG*sizeof(float));
	//���浱ǰ֡����p_in_data��p_overlap
	memcpy(p_overlap,p_in_data,BLOCK_LEN_LONG*sizeof(float));

	for(i=0;i<BLOCK_LEN_LONG;i++){
		p_out_mdct[i]*=p_kbd_window[i];
		p_out_mdct[i+BLOCK_LEN_LONG]*=p_kbd_window[BLOCK_LEN_LONG-1-i];
	}

	MDCT(fft_tables,p_out_mdct,2*BLOCK_LEN_LONG);	

	
	
}

/******************************************************************************************
  *  @brief      IMDCT�任��ȥ50%��KBD�ص���.
  *  @param[in]  p_in_data    1024 ��ǰ֡MDCT���� ��Ϊ��ǰ�任��ĺ�1024������ 
				 p_overlap    1024 ��һ֡��IMDCTȥ50%�ص�����ĺ�1024��ʱ������ 
				 �뵱ǰ֡ȥ�����ǰ1024��ʱ�����ݵ�����Ϊ��ǰ֡ʱ���������
			     p_kbd_window 1024 kbd����ǰ1024����  ��1024������ǰ1024����Գ� 	
				 
  *  @param[out] p_out_data   1024 IMDCT�任ȥ50%�ص������ʱ������
  *  @Return:    
  *  @note       
  *  @date       2016.09.23
  *  @see        
  *****************************************************************************************/

void IMDCTKbdWindow(FFT_Tables *fft_tables,float *p_overlap,float *p_in_data,float* p_kbd_window,float* p_out_data)
{
	int i,j;
	float *o_buf, *transf_buf, *overlap_buf;

	transf_buf = (float*)malloc(2*BLOCK_LEN_LONG*sizeof(float)); //2048
	overlap_buf = (float*)malloc(2*BLOCK_LEN_LONG*sizeof(float));//2048

	memcpy(overlap_buf,p_overlap,BLOCK_LEN_LONG*sizeof(float));//1024  ǰ1024������������һ֡�ĺ�벿�ֵļӴ�ʱ������
	o_buf = overlap_buf;//2048  ǰ1024�����ݼ�Ϊ��һ֡�ĺ�벿�ֵļӴ�ʱ������

	memcpy(transf_buf,p_in_data,BLOCK_LEN_LONG*sizeof(float));//1024  ��ǰ֡��mdct�����ݱ任��ʱ��
	IMDCT(fft_tables,transf_buf,2*BLOCK_LEN_LONG);//transf_buf�д洢�� ʱ������ 2048��

	for ( i = 0 ; i < BLOCK_LEN_LONG ; i++){// ����transf_buf�д洢�� ��ǰ֡�Ӵ����ʱ������ 2048��
		transf_buf[i] *= p_kbd_window[i];
		transf_buf[i+BLOCK_LEN_LONG]*=p_kbd_window[BLOCK_LEN_LONG-1-i];
	}

	for ( i = 0 ; i < BLOCK_LEN_LONG ; i++){
		o_buf[i]+=transf_buf[i];//��ǰһ֡�ĺ�벿�ּӴ�ʱ�����ݵ��� 1024 ���õ���ǰ֡���������1024��  
		o_buf[i+BLOCK_LEN_LONG] = transf_buf[i+BLOCK_LEN_LONG];//o_buf��1024������ ���浱ǰ֡�ļӴ���ĺ�1024������ ������һ֡��Ϊ�ص���Ӳ���
	}

	memcpy(p_out_data,o_buf,BLOCK_LEN_LONG*sizeof(float));//p_out_data �������1024�� 
	/* save unused output data */
	memcpy(p_overlap,o_buf+BLOCK_LEN_LONG,BLOCK_LEN_LONG*sizeof(float));//��ǰ֡ʱ��Ӵ����� ��벿��


	if (overlap_buf) free(overlap_buf);
	if (transf_buf) free(transf_buf);

}