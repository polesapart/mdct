#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "fft.h"
#define MAXLOGM 9
#define MAXLOGR 8

/******************************************************************************************
  *  @brief      FFT ��ʼ��.
  *  @param[in]  
  *  @param[out] 
  *  @Return:    
  *  @note       
  *  @date       2016.09.29
  *  @see        
  *****************************************************************************************/
void fft_initialize( FFT_Tables *fft_tables )     //FFT��ʼ����fft_tables
{
	int i;    // FFT_Table�ǽṹ��
	//(MAXLOGM+1)��fft�����㼶�� costbl  negsintbl reordertbl �൱�ڶ�ά�����׵�ַ
	// costbl[0]  negsintbl[0] reordertbl[0] �൱��һά�����׵�ַ 
	//sizeof( fft_tables->costbl[0])�õ�����һά����������Ԫ���ܹ�ռ���ֽ���
	//costbl[i][j]  negsintbl[i][j]�����ǵ�i�������е�j����ת���ӵ����ҷ��������ҷ���
	//fft_tables->costbl ���Ƕ�λ�����׵�ַ  ��ΪNULL
	fft_tables->costbl		= (fftfloat **)malloc( (MAXLOGM+1) * sizeof( fft_tables->costbl[0] ) );
	fft_tables->negsintbl	= (fftfloat **)malloc( (MAXLOGM+1) * sizeof( fft_tables->negsintbl[0] ) );
	fft_tables->reordertbl	= (unsigned short **)malloc( (MAXLOGM+1) * sizeof( fft_tables->reordertbl[0] ) );
	
	//һά�����׵�ַΪnull
	for( i = 0; i< MAXLOGM+1; i++ )
	{
		fft_tables->costbl[i]		= NULL;
		fft_tables->negsintbl[i]	= NULL;
		fft_tables->reordertbl[i]	= NULL;
	}
}

/******************************************************************************************
  *  @brief      FFT�任������ �����ͷ�����ڴ�ȣ�.
  *  @param[in]  
  *  @param[out] 
  *  @Return:    
  *  @note       
  *  @date       2016.09.29
  *  @see        
  *****************************************************************************************/

void fft_terminate( FFT_Tables *fft_tables )  //FFT��ֹ
{ 
	int i;
	for( i = 0; i< MAXLOGM+1; i++ )
	{
		if( fft_tables->costbl[i] != NULL )
			free( fft_tables->costbl[i] );
		
		if( fft_tables->negsintbl[i] != NULL )
			free( fft_tables->negsintbl[i] );
			
		if( fft_tables->reordertbl[i] != NULL )
			free( fft_tables->reordertbl[i] );
	}

	free( fft_tables->costbl );
	free( fft_tables->negsintbl );
	free( fft_tables->reordertbl );

	fft_tables->costbl		= NULL;
	fft_tables->negsintbl	= NULL;
	fft_tables->reordertbl	= NULL;
}

/*************************************************************************
*Function����λ���� �����������ݵ��˳��
*Paras��fft_tables�� x��x[]�д�ŵ��������ݵ� logm��fft���㼶��
*Return��void
*************************************************************************/
static void reorder( FFT_Tables *fft_tables, float *x, int logm)   
{
	int i;
	int size = 1 << logm;//size��fft�ĵ���
	unsigned short *r;	//size

	//��fft_initialize( FFT_Tables *fft_tables )���ú�reordertbl[logm]��reordertbl[i]����Ϊ�� 
	// ����һ��λ��ת�� create bit reversing table
	if ( fft_tables->reordertbl[logm] == NULL ) 
	{
		fft_tables->reordertbl[logm] =(unsigned short *) malloc(size * sizeof(*(fft_tables->reordertbl[0])));
		//reordertbl[logm] ����Ϊ��
		for (i = 0; i < size; i++)//i��ԭλ���±� 
		{
			int reversed = 0;//��λ��ʶ
			int b0;
			int tmp = i;

			for (b0 = 0; b0 < logm; b0++)//
			{
				reversed = (reversed << 1) | (tmp & 1);//(tmp & 1)��ֵΪ1
				tmp >>= 1;//  tmp/2
			}
			fft_tables->reordertbl[logm][i] = reversed;
		}
	}

	r = fft_tables->reordertbl[logm];

	for (i = 0; i < size; i++)
	{
		int j = r[i];
		float tmp;

		if (j <= i)
			continue;
		//����x[i]��x[j] ���������ݵ����x[]��
		tmp = x[i];
		x[i] = x[j];
		x[j] = tmp;
	}
}

/*************************************************************************
*Function��fft�任����
*Paras��xr���������ݵ�ʵ�������׵�ַ xi���������ݵ��鲿�����׵�ַ
*            refac��ʵ�����ϵ��  imfac���鲿���ϵ��  size:fft����
*Return��void
*************************************************************************/
static void fft_proc(
		float *xr, 
		float *xi,
		fftfloat *refac, 
		fftfloat *imfac, 
		int size)	
{
	int step, shift, pos;
	int exp, estep;

	estep = size;
	for (step = 1; step < size; step *= 2)
	{
		int x1;
		int x2 = 0;
		estep >>= 1;
		for (pos = 0; pos < size; pos += (2 * step))
		{
			x1 = x2;
			x2 += step;
			exp = 0;
			for (shift = 0; shift < step; shift++)
			{
				float v2r, v2i;

				v2r = xr[x2] * refac[exp] - xi[x2] * imfac[exp];
				v2i = xr[x2] * imfac[exp] + xi[x2] * refac[exp];

				xr[x2] = xr[x1] - v2r;
				xr[x1] += v2r;

				xi[x2] = xi[x1] - v2i;

				xi[x1] += v2i;

				exp += estep;

				x1++;
				x2++;
			}
		}
	}
}

static void check_tables( FFT_Tables *fft_tables, int logm)
{
	if( fft_tables->costbl[logm] == NULL )
	{
		int i;
		int size = 1 << logm;

		if( fft_tables->negsintbl[logm] != NULL )
			free( fft_tables->negsintbl[logm] );

		fft_tables->costbl[logm]	= malloc((size / 2) * sizeof(*(fft_tables->costbl[0])));
		fft_tables->negsintbl[logm]	= malloc((size / 2) * sizeof(*(fft_tables->negsintbl[0])));

		for (i = 0; i < (size >> 1); i++)
		{
			float theta = 2.0 * M_PI * ((float) i) / (float) size;
			fft_tables->costbl[logm][i]		= cos(theta);
			fft_tables->negsintbl[logm][i]	= -sin(theta);
		}
	}
}


/*************************************************************************
*Function��fft�任���ú���
*Paras��ffttables ��xr���������ݵ�ʵ�������׵�ַ xi���������ݵ��鲿�����׵�ַ
*                  logm�����㼶��
*Return��void
*************************************************************************/
void fft( FFT_Tables *fft_tables, float *xr, float *xi, int logm)			//FFT�任
{
	if (logm > MAXLOGM)
	{
		fprintf(stderr, "fft size too big\n");
		exit(1);
	}

	if (logm < 1)
	{
		//printf("logm < 1\n");
		return;
	}

	check_tables( fft_tables, logm);       //���FFT����
	//��������ʵ��xr�� �鲿xi����
	reorder( fft_tables, xr, logm);
	reorder( fft_tables, xi, logm);

	fft_proc( xr, xi, fft_tables->costbl[logm], fft_tables->negsintbl[logm], 1 << logm );   //FFT�Ĺ���
}

//��������ȫ��Ϊʵ���� fft�任
void rfft( FFT_Tables *fft_tables, float *x, int logm)     
{
	float xi[1 << MAXLOGR];

	if (logm > MAXLOGR)
	{
		fprintf(stderr, "rfft size too big\n");
		exit(1);
	}

	memset(xi, 0, (1 << logm) * sizeof(xi[0]));

	fft( fft_tables, x, xi, logm);

	memcpy(x + (1 << (logm - 1)), xi, (1 << (logm - 1)) * sizeof(*x));
}

//fft���任
void ffti( FFT_Tables *fft_tables, float *xr, float *xi, int logm)
{
	int i, size;
	float fac;
	float *xrp, *xip;

	fft( fft_tables, xi, xr, logm);

	size = 1 << logm;
	fac = 1.0 / size;
	xrp = xr;
	xip = xi;

	for (i = 0; i < size; i++)
	{
		*xrp++ *= fac;
		*xip++ *= fac;
	}
}

