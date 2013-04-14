/**
 * @file colortransfer.c
 * @brief Global color transfer
 *
 * Usage: 
 *    #./colortransfer source.png target.png output.png
 *
 * Parameters:
 *
 * source image (format: png)
 * target image (format: png)
 * output image (format: png)
 *
 * Read/write operations (png format) make use
 * of io_png.c and io_png.h

 * @author Shan An <anshan.tju@gmail.com>
 * 
 * Erik Reinhard, et al., "Color Transfer between Images", 2001.
 */

/**
 * @mainpage Global color transfer
 *
 * README.txt:
 * @verbinclude README.txt
 */


/**
 * @file   colortransfer.c
 * @brief  Main executable file
 *
 *
 * @author Shan An <anshan.tju@gmail.com>
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "io_png/io_png.h"

#define NEWA(type,n) (type*)malloc(sizeof(type)*(n))
#define NEW(type) NEWA(type,1)

typedef unsigned char uchar;

/*Transformation Matrix*/
static float rgb2lms[3][3] = {
	{ 0.3811, 0.5783, 0.0402 },
	{ 0.1967, 0.7244, 0.0782 },
	{ 0.0241, 0.1288, 0.8444 }
};

static float lms2rgb[3][3] = {
	{  4.4679, -3.5873,  0.1193 },
	{ -1.2186,  2.3809, -0.1624 },
	{  0.0497, -0.2439,  1.2045 }
};

static float lms2lab1[3][3] = {
	{ 1.0 / 1.7321, 0.0, 0.0 },
	{ 0.0, 1.0 / 2.4495, 0.0 },
	{ 0.0, 0.0, 1.0 / 1.4142 }
};

static float lms2lab2[3][3] = {
	{ 1.0,  1.0,  1.0 },
	{ 1.0,  1.0, -2.0 },
	{ 1.0, -1.0,  0.0 }
};

/* Matrix operation */
float *mat_mul(float *A, float *B) 
{
	int i, j, k;
	float *C = NEWA(float, 9);
	for(i = 0; i < 3; i++) {
		for(j = 0; j < 3; j++) {
			*(C + i*3 + j) = 0.0;
			for(k = 0; k < 3; k++) {
				*(C + i*3 + j)  += *(A + i*3 + k) * (*(B + k*3 + j));
			}
		}
	}
	return C;
}

float *mat_mul_vec(float *A, float *b) 
{
	int i, j;
	float *C = NEWA(float, 3);
	for(i = 0; i < 3; i++) {
		*(C + i)= 0.0;
		for(j = 0; j < 3; j++) {
			*(C + i) += *(A + i*3 + j) * (*(b + j));
		}
	}
	return C;
}

float *mat_trans(float *A) 
{
	int i, j;
	float *T = NEWA(float, 9);
	for(i = 0; i < 3; i++) {
		for(j = 0; j < 3; j++) {
			*(T + i*3 + j) = *(A + j*3 + i);
		}
	}
	return T;
}

/*Calculate the statistics in Lab color space.*/
float *getStatLab(float *lab, int w, int h) 
{
	int i, j;
	float *stat = NEWA(float, 6);
	for (i = 0; i < 6; i++){
		*(stat + i) = 0.0;
	}

	for(j = 0; j < h; j++) {
		for(i = 0; i < w; i++) {
			*(stat + 0) += lab[j*w*3 + i*3 + 0];
			*(stat + 1) += lab[j*w*3 + i*3 + 1];
			*(stat + 2) += lab[j*w*3 + i*3 + 2];
			*(stat + 3) += lab[j*w*3 + i*3 + 0] * lab[j*w*3 + i*3 + 0];
			*(stat + 4) += lab[j*w*3 + i*3 + 1] * lab[j*w*3 + i*3 + 1];
			*(stat + 5) += lab[j*w*3 + i*3 + 2] * lab[j*w*3 + i*3 + 2];
		}
	}

	for(i = 0; i < 6; i++) {
		*(stat + i) = stat[i] / (float)(w * h);
	}

	for(i = 3; i < 6; i++) {
		*(stat + i) = sqrt(stat[i] - stat[i-3] * stat[i-3]);
	}

	return stat;
}

/*Trans RGB to Lab.*/
void transRGB2Lab(float *lab, float *lms2lab, uchar *R, uchar *G, uchar *B, int w, int h) 
{
	int i, j;
	for(j = 0; j < h; j++) {
		for(i = 0; i < w; i++) {
			float c[3] = {(float)R[j*w + i], (float)G[j*w + i], (float)B[j*w + i]};
			float *pLMS = mat_mul_vec(*rgb2lms, c);
			pLMS[0] = (pLMS[0] > 1.0e-5) ? log10(pLMS[0]) : -5.0;
			pLMS[1] = (pLMS[1] > 1.0e-5) ? log10(pLMS[1]) : -5.0;
			pLMS[2] = (pLMS[2] > 1.0e-5) ? log10(pLMS[2]) : -5.0;
			float *pLab = mat_mul_vec(lms2lab, pLMS);

			lab[j*w*3 + i*3 + 0] = pLab[0];
			lab[j*w*3 + i*3 + 1] = pLab[1];
			lab[j*w*3 + i*3 + 2] = pLab[2];
			free(pLMS);
			free(pLab);
		}
	}
}

/*Trans Lab to RGB.*/
void transLab2RGB(float *lab, float *lab2lms, uchar *R, uchar *G, uchar *B, int w, int h) 
{
	int i, j;
	for(j = 0; j < h; j++) {
		for(i =0 ; i < w; i++) {
			float pLab[3];
			pLab[0] = lab[j*w*3 + i*3 + 0];
			pLab[1] = lab[j*w*3 + i*3 + 1];
			pLab[2] = lab[j*w*3 + i*3 + 2];

			float *pLMS = mat_mul_vec(lab2lms, pLab);
			pLMS[0] = (pLMS[0] > -5.0) ? pow(10.0, pLMS[0]) : 1.0e-5;
			pLMS[1] = (pLMS[1] > -5.0) ? pow(10.0, pLMS[1]) : 1.0e-5;
			pLMS[2] = (pLMS[2] > -5.0) ? pow(10.0, pLMS[2]) : 1.0e-5;

			float *c = mat_mul_vec(*lms2rgb, pLMS);
			int r = (int)c[0];
			int g = (int)c[1];
			int b = (int)c[2];

			r = (r >= 0) ? (r < 256) ? r : 255 : 0;
			g = (g >= 0) ? (g < 256) ? g : 255 : 0;
			b = (b >= 0) ? (b < 256) ? b : 255 : 0;

			R[j*w + i] = (uchar)r;
			G[j*w + i] = (uchar)g;
			B[j*w + i] = (uchar)b;
			free(pLMS);
			free(c);
		}
	}
}

/* 
Converts 'input' array of color data from 
format RRRR...GGGG...BBBB to 3 independent arrays
'RR', 'GG', 'BB' each one of them storing the data
corresponding to Red, Green and Blue channels, respectively.
'size' is the length of the 'input' array. 
*/
void input2RGB(uchar *input,
               uchar **RR, uchar **GG, uchar **BB,
               int size) 
{
	uchar *R, *G, *B;
	int n;

	R= (uchar*) malloc(size*sizeof(uchar));
	G= (uchar*) malloc(size*sizeof(uchar));
	B= (uchar*) malloc(size*sizeof(uchar));

    for (n=0; n < size; n++) {
		R[n]=input[n];
		G[n]=input[size + n];
        B[n]=input[2*size + n];
    }

    	*RR= R;
   	*GG= G;
   	*BB= B;
}


/* 
Converts data stored in 3 independent arrays 'R', 'G' and 'B'
to a single array 'output', with format RRR...GGG...BBB.
Data are converted to uchar (range [0, 255]).
'size' is the length of the 'output' array. 
*/
void RGB2output(uchar *R, uchar *G, uchar *B,
                uchar *output, int size) 
{
	int n;

	for (n=0; n < size; n++) {
       	output[n] = R[n];
        output[size + n] = G[n];
        output[2*size + n] = B[n];
    }
}

/**
 * @brief main function call
 */
int main(int argc, const char **argv) 
{
	const char *namesrc, *nametar, *nameout;
	uchar *imgsrc, *imgtar, *output;
	int src_w, src_h, tar_w, tar_h;
    size_t w1, h1, w2, h2;
	int nc, i, j;
	uchar *Rsrc, *Gsrc, *Bsrc, *Rtar, *Gtar, *Btar, *Rout, *Gout, *Bout;

	namesrc = argv[1];
	nametar = argv[2];
	nameout = argv[3];

 	if (argc < 4) {
        printf("Usage: colortransfer source.png target.png output.png \n");
        return EXIT_FAILURE;
  	}
    /* read the PNG source and target image */
  	imgsrc = read_png_u8_rgb(namesrc, &w1, &h1);
	imgtar = read_png_u8_rgb(nametar, &w2, &h2);
	src_w = (int) w1;
	src_h = (int) h1;
	tar_w = (int) w2;
	tar_h = (int) h2;

  	input2RGB(imgsrc, &Rsrc, &Gsrc, &Bsrc, src_w*src_h);
	input2RGB(imgtar, &Rtar, &Gtar, &Btar, tar_w*tar_h);

	/* transform pixels' value from RGB to Lab */
    float *lab_src = NEWA(float, src_w*src_h*3*sizeof(float));
	float *lab_tar = NEWA(float, tar_w*tar_h*3*sizeof(float));

	float *lms2lab = mat_mul(*lms2lab1, *lms2lab2);
	float *lab2lms = mat_mul(mat_trans(*lms2lab2), *lms2lab1);

	transRGB2Lab(lab_src, lms2lab, Rsrc, Gsrc, Bsrc, src_w, src_h);
	transRGB2Lab(lab_tar, lms2lab, Rtar, Gtar, Btar, tar_w, tar_h);

	/* calculate the statistics */
	float *stat_lab_src = getStatLab(lab_src, src_w, src_h);
	float *stat_lab_tar = getStatLab(lab_tar, tar_w, tar_h);
	printf("source: mean = (%1.4f, %1.4f, %1.4f), std dev = (%1.4f, %1.4f, %1.4f)\n", stat_lab_src[0], stat_lab_src[1], stat_lab_src[2], stat_lab_src[3], stat_lab_src[4], stat_lab_src[5]);
	printf("target: mean = (%1.4f, %1.4f, %1.4f), std dev = (%1.4f, %1.4f, %1.4f)\n", stat_lab_tar[0], stat_lab_tar[1], stat_lab_tar[2], stat_lab_tar[3], stat_lab_tar[4], stat_lab_tar[5]);

	/* color transfer */
	for(j = 0; j < src_h; j++) {
		for(i = 0; i < src_w; i++) {
			lab_src[j*src_w*3 + i*3 + 0] = (lab_src[j*src_w*3 + i*3 + 0] - stat_lab_src[0]) * stat_lab_tar[3] / stat_lab_src[3] + stat_lab_tar[0];
			lab_src[j*src_w*3 + i*3 + 1] = (lab_src[j*src_w*3 + i*3 + 1] - stat_lab_src[1]) * stat_lab_tar[4] / stat_lab_src[4] + stat_lab_tar[1];
			lab_src[j*src_w*3 + i*3 + 2] = (lab_src[j*src_w*3 + i*3 + 2] - stat_lab_src[2]) * stat_lab_tar[5] / stat_lab_src[5] + stat_lab_tar[2];
		}
	}
	
	/* output */
   	Rout = NEWA(uchar, src_w*src_h*sizeof(uchar));
   	Gout = NEWA(uchar, src_w*src_h*sizeof(uchar));
   	Bout = NEWA(uchar, src_w*src_h*sizeof(uchar));

	transLab2RGB(lab_src, lab2lms, Rout, Gout, Bout, src_w, src_h);

    	output=NEWA(uchar, 3*src_w*src_h);
   	RGB2output(Rout, Gout, Bout, output, src_w*src_h);

    /* write the PNG output image */
   	nc=3; /*color image*/
    write_png_u8(nameout, output, src_w, src_h, nc);

	free(lab_src);
	free(lab_tar);
	free(lms2lab);
	free(lab2lms);
	free(stat_lab_src);
	free(stat_lab_tar);
	free(Rout);
	free(Gout);
	free(Bout);
	free(output);

    return EXIT_SUCCESS;
}
