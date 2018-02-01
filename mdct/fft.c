#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "fft.h"
#define MAXLOGM 9
#define MAXLOGR 8

/******************************************************************************************
  *  @brief      FFT 初始化.
  *  @param[in]  
  *  @param[out] 
  *  @Return:    
  *  @note       
  *  @date       2016.09.29
  *  @see        
  *****************************************************************************************/
void fft_initialize( FFT_Tables *fft_tables )     //FFT初始化，fft_tables
{
	int i;    // FFT_Table是结构体
	//(MAXLOGM+1)是fft的运算级数 costbl  negsintbl reordertbl 相当于二维数组首地址
	// costbl[0]  negsintbl[0] reordertbl[0] 相当于一维数组首地址 
	//sizeof( fft_tables->costbl[0])得到的是一维数组中所有元素总共占的字节数
	//costbl[i][j]  negsintbl[i][j]中则是第i级运算中第j个旋转因子的余弦分量和正弦分量
	//fft_tables->costbl 则是二位数组首地址  不为NULL
	fft_tables->costbl		= (fftfloat **)malloc( (MAXLOGM+1) * sizeof( fft_tables->costbl[0] ) );
	fft_tables->negsintbl	= (fftfloat **)malloc( (MAXLOGM+1) * sizeof( fft_tables->negsintbl[0] ) );
	fft_tables->reordertbl	= (unsigned short **)malloc( (MAXLOGM+1) * sizeof( fft_tables->reordertbl[0] ) );
	
	//一维数组首地址为null
	for( i = 0; i< MAXLOGM+1; i++ )
	{
		fft_tables->costbl[i]		= NULL;
		fft_tables->negsintbl[i]	= NULL;
		fft_tables->reordertbl[i]	= NULL;
	}
}

/******************************************************************************************
  *  @brief      FFT变换结束后 后处理（释放相关内存等）.
  *  @param[in]  
  *  @param[out] 
  *  @Return:    
  *  @note       
  *  @date       2016.09.29
  *  @see        
  *****************************************************************************************/

void fft_terminate( FFT_Tables *fft_tables )  //FFT终止
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
*Function：倒位函数 调整输入数据点的顺序
*Paras：fft_tables： x：x[]中存放的输入数据点 logm：fft运算级数
*Return：void
*************************************************************************/
static void reorder( FFT_Tables *fft_tables, float *x, int logm)   
{
	int i;
	int size = 1 << logm;//size是fft的点数
	unsigned short *r;	//size

	//在fft_initialize( FFT_Tables *fft_tables )调用后reordertbl[logm]（reordertbl[i]）即为空 
	// 创建一个位翻转表 create bit reversing table
	if ( fft_tables->reordertbl[logm] == NULL ) 
	{
		fft_tables->reordertbl[logm] =(unsigned short *) malloc(size * sizeof(*(fft_tables->reordertbl[0])));
		//reordertbl[logm] 不再为空
		for (i = 0; i < size; i++)//i是原位序下标 
		{
			int reversed = 0;//倒位标识
			int b0;
			int tmp = i;

			for (b0 = 0; b0 < logm; b0++)//
			{
				reversed = (reversed << 1) | (tmp & 1);//(tmp & 1)的值为1
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
		//交换x[i]和x[j] 真正的数据点放在x[]中
		tmp = x[i];
		x[i] = x[j];
		x[j] = tmp;
	}
}

/*************************************************************************
*Function：fft变换函数
*Paras：xr：输入数据点实部数组首地址 xi：输入数据点虚部数组首地址
*            refac：实部相乘系数  imfac：虚部相乘系数  size:fft点数
*Return：void
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
*Function：fft变换调用函数
*Paras：ffttables ：xr：输入数据点实部数组首地址 xi：输入数据点虚部数组首地址
*                  logm：运算级数
*Return：void
*************************************************************************/
void fft( FFT_Tables *fft_tables, float *xr, float *xi, int logm)			//FFT变换
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

	check_tables( fft_tables, logm);       //检查FFT长度
	//输入数据实部xr， 虚部xi倒序
	reorder( fft_tables, xr, logm);
	reorder( fft_tables, xi, logm);

	fft_proc( xr, xi, fft_tables->costbl[logm], fft_tables->negsintbl[logm], 1 << logm );   //FFT的过程
}

//输入数据全部为实数的 fft变换
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

//fft反变换
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

