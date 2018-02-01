#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include "wavread.h"
#include "fft.h"
#include "mdct.h"
#include "wavwrite.h"



int main(int argc, char *argv[])
{
	char *input_wav_name="test.wav";
	char *output_wav_name="test1.wav";
	int samplesInput=FRAME_LEN;
	int samplesRead=0;
	int map[1]={0};
	int frame_num=0;

	PcmFile *infile;
	WavFile *outfile;
	
	FFT_Tables* fft_tables;
	float pcm_overlap_buf[FRAME_LEN];
	float pcm_buf[FRAME_LEN];
	float *kbd_window_long;
	float *p_out_mdct;

	float mdct_freq_buf[FRAME_LEN];
	float imdct_overlap_buf[FRAME_LEN];
	float imdct_time_buf[FRAME_LEN];

	float diff[1024];
	float max_diff=0;
	infile= openWavRead(input_wav_name, 0);
	outfile=openWavWrite(output_wav_name,infile->samplerate,infile->channels,WAV_FMT_16BIT,OUTPUT_WAV,0);

	fft_tables=(FFT_Tables*)malloc(sizeof(FFT_Tables));
	fft_initialize(fft_tables);

	memset(pcm_overlap_buf,0,FRAME_LEN*sizeof(float));
	memset(pcm_buf,0,FRAME_LEN*sizeof(float));

	kbd_window_long=(float*)malloc(FRAME_LEN*sizeof(float));
	CalculateKBDWindow(kbd_window_long,4,2*FRAME_LEN);

	p_out_mdct=(float *)malloc(2*FRAME_LEN*sizeof(float));
	memset(p_out_mdct,0,2*FRAME_LEN*sizeof(float));

	memset(mdct_freq_buf,0,FRAME_LEN*sizeof(float));
	memset(imdct_overlap_buf,0,FRAME_LEN*sizeof(float));
	memset(imdct_time_buf,0,FRAME_LEN*sizeof(float));



	while(1){
		int i=0,j=0,k=0;
		short mixPcmData[FRAME_LEN]={0};
		samplesRead = wavReadFloat32(infile, pcm_buf, samplesInput,map);
		frame_num++;
		
		if(samplesRead<1024) 
			break; //丢掉最后一帧 		
		/*对每帧时域信号进行2048点的mdct变换*/
		MDCTKbdWindow(fft_tables,pcm_overlap_buf,pcm_buf,kbd_window_long,p_out_mdct);
		memcpy(mdct_freq_buf, p_out_mdct, FRAME_LEN*sizeof(float));
		IMDCTKbdWindow(fft_tables,imdct_overlap_buf,mdct_freq_buf,kbd_window_long,imdct_time_buf);
		
		//for(j=0;j<FRAME_LEN;j++){
		//	diff[j]=pcm_buf[j]-imdct_time_buf[j];
		//	if(diff[j]>max_diff){
		//		max_diff=diff[j];
		//	}
		//}
		//fprintf(stdout,"%d:  %f\n",frame_num,max_diff);

		for(j=0;j<FRAME_LEN;j++){
			mixPcmData[j]=imdct_time_buf[j]*32767.0f;
		}	
		
		writeWavData(outfile,mixPcmData,FRAME_LEN,0);
		
	}

	if(infile){
		closeWavRead(infile);
		infile=NULL;
	}

	if(outfile){
		closeWavWrite(outfile);
		outfile=NULL;
	}

	fft_terminate(fft_tables);

	if(kbd_window_long){
		free(kbd_window_long);
		kbd_window_long=NULL;
	}
	if(p_out_mdct){
		free(p_out_mdct);
		p_out_mdct=NULL;
	}



	system("pause");
	return 0;
}
