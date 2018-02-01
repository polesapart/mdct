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
  *  @brief      计算KBD窗函数的值  只需要计算前半部分 后半部分与前半部分对称.
  *  @param[in]  win：存储窗函数值  alpha =4（一般情况）  length=N/2=1024 只计算前半部分
  *  @param[out] win即为输出的窗函数值
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
  *  @brief      MDCT变换.
  *  @param[in]  data：时域数据（N点有用）N=2048
			     N:变换点数 N=2048
  *  @param[out] data:mdct变换后的频域数据 N点（前N/2点与后N/2点成奇对称）
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

	xi = (float*)malloc((N >> 2)*sizeof(float)); //xi N/4点
	xr = (float*)malloc((N >> 2)*sizeof(float)); //xr N/4点

	/* prepare for recurrence（递推） relation in pre-twiddle */
	cfreq = cos (freq);
	sfreq = sin (freq);
	cosfreq8 = cos (freq * 0.125);
	sinfreq8 = sin (freq * 0.125);
	c = cosfreq8;
	s = sinfreq8;

	for (i = 0; i < (N >> 2); i++) {//i=0 1 2 .... 511
		/* calculate real and imaginary parts of g(n) or G(p)   计算实部和虚部*/
		n = (N >> 1) - 1 - 2 * i;//n=1023 1021 ...  3 1
		
		//用于N/4（512）点复数FFT变换的 实部数据的获取
		if (i < (N >> 3)) //i=0 1 2 ....254 255
			//前256个实部数据的获取
			//(N >> 2) + n=0.75*N-1-2*i=1535 1533 1531 ... 1027 1025 ;   N + (N >> 2) - 1 - n=0.75*N+2*i=1536 1538 ...2044  2046
			tempr = data [(N >> 2) + n] + data [N + (N >> 2) - 1 - n]; /* use second form of e(n) for n = N / 2 - 1 - 2i     */
		else  //i=256 257....510 511
			//后256个实部数据的获取
			//(N >> 2) + n=0.75*N-1-2*i=1023 1021 ....515 513 ;   (N >> 2) - 1 - n= -0.25*N+2*i= 0 2 4... 508 510
			tempr = data [(N >> 2) + n] - data [(N >> 2) - 1 - n]; /* use first form of e(n) for n = N / 2 - 1 - 2i */


		n = 2 * i;//n=0 2 4 ....1020 1022

		//用于N/4（512）点复数FFT变换的 虚部数据的获取
		if (i < (N >> 3)) //i=0 1 2 ....254 255
			//前256个虚部数据的获取
			//(N >> 2) + n=0.25*N+2*i=512 514  ... 1020 1022  ;   (N >> 2) - 1 - n= 0.25*N-1-2*i= 511 509 ....3 1
			tempi = data [(N >> 2) + n] - data [(N >> 2) - 1 - n]; /* use first form of e(n) for n=2i */
		else  //i=256 257....510 511
			//后256个虚部数据的获取
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
	case BLOCK_LEN_SHORT * 2: //128*2==256==N 点MDCT变换  则对应的fft变换为256/4=64点（N/4点）
		fft( fft_tables, xr, xi, 6);
		break;
	case BLOCK_LEN_LONG * 2: //1024*2==2048==N 点MDCT变换  则对应的fft变换为2048/4=512点（N/4点）
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

		//变换后的系数（这里还是存放在data里面） 是奇对称的（时域输入N点，则频域前N/2点与后N/2点奇对称，故只需要前N/2点的数据便可以明确mdct变换）
		/* fill in output values */
		data [2 * i] = -tempr;   /* first half even： 0 2 4 6 ....1018 1020 1022 */  
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
  *  @brief      IMDCT变换.
  *  @param[in]  data：频域MDCT数据（N/2点有用）
			     N:变换点数 N=2048
  *  @param[out] data:imdct变换后的时域数据 N点
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
  *  @brief      加50%重叠KBD窗的MDCT变换.
  *  @param[in]  p_in_data    1024 当前帧时域数据 作为当前变换块的后1024个数据 
				 p_overlap    1024 上一帧时域数据 作为当前变换块的前1024个数据 
			     p_kbd_window 1024 kbd窗的前1024个点  后1024个点与前1024个点对称 	
				 
  *  @param[out] p_out_mdct   2048 mdct频域数据 呈奇对称 在IMDCT变换中只用前1024个数据
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
	//缓存当前帧数据p_in_data到p_overlap
	memcpy(p_overlap,p_in_data,BLOCK_LEN_LONG*sizeof(float));

	for(i=0;i<BLOCK_LEN_LONG;i++){
		p_out_mdct[i]*=p_kbd_window[i];
		p_out_mdct[i+BLOCK_LEN_LONG]*=p_kbd_window[BLOCK_LEN_LONG-1-i];
	}

	MDCT(fft_tables,p_out_mdct,2*BLOCK_LEN_LONG);	

	
	
}

/******************************************************************************************
  *  @brief      IMDCT变换后去50%的KBD重叠窗.
  *  @param[in]  p_in_data    1024 当前帧MDCT数据 作为当前变换块的后1024个数据 
				 p_overlap    1024 上一帧经IMDCT去50%重叠窗后的后1024点时域数据 
				 与当前帧去窗后的前1024个时域数据叠加作为当前帧时域数据输出
			     p_kbd_window 1024 kbd窗的前1024个点  后1024个点与前1024个点对称 	
				 
  *  @param[out] p_out_data   1024 IMDCT变换去50%重叠窗后的时域数据
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

	memcpy(overlap_buf,p_overlap,BLOCK_LEN_LONG*sizeof(float));//1024  前1024个数据来自上一帧的后半部分的加窗时域数据
	o_buf = overlap_buf;//2048  前1024个数据即为上一帧的后半部分的加窗时域数据

	memcpy(transf_buf,p_in_data,BLOCK_LEN_LONG*sizeof(float));//1024  当前帧的mdct域数据变换到时域
	IMDCT(fft_tables,transf_buf,2*BLOCK_LEN_LONG);//transf_buf中存储着 时域数据 2048点

	for ( i = 0 ; i < BLOCK_LEN_LONG ; i++){// 现在transf_buf中存储着 当前帧加窗后的时域数据 2048点
		transf_buf[i] *= p_kbd_window[i];
		transf_buf[i+BLOCK_LEN_LONG]*=p_kbd_window[BLOCK_LEN_LONG-1-i];
	}

	for ( i = 0 ; i < BLOCK_LEN_LONG ; i++){
		o_buf[i]+=transf_buf[i];//与前一帧的后半部分加窗时域数据叠加 1024 ，得到当前帧的输出数据1024点  
		o_buf[i+BLOCK_LEN_LONG] = transf_buf[i+BLOCK_LEN_LONG];//o_buf后1024个数据 保存当前帧的加窗后的后1024个数据 留给下一帧作为重叠相加部分
	}

	memcpy(p_out_data,o_buf,BLOCK_LEN_LONG*sizeof(float));//p_out_data 输出数据1024点 
	/* save unused output data */
	memcpy(p_overlap,o_buf+BLOCK_LEN_LONG,BLOCK_LEN_LONG*sizeof(float));//当前帧时域加窗数据 后半部分


	if (overlap_buf) free(overlap_buf);
	if (transf_buf) free(transf_buf);

}