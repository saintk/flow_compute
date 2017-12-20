#include "stdafx.h"
#include <iostream>  
#include "opencv2/opencv.hpp"  
#include "opencv2/imgproc.hpp"
#include "opencv2/highgui.hpp"
#include "cmath"
#include "flow.h"
using namespace std;
#define FRAME_SIZE	64
#define SEARCH_SIZE	4 // maximum offset to search: 4 + 1/2 pixels
#define TILE_SIZE	8               						// x & y tile size
#define NUM_BLOCKS	5 // x & y number of tiles to check


uint8_t compute_flow(uint8_t *image1, uint8_t *image2, float *pixel_flow_x, float *pixel_flow_y)
{

	/* constants */
	const int16_t winmin = -SEARCH_SIZE; //search_size == 4
	const int16_t winmax = SEARCH_SIZE;
	const uint16_t hist_size = 2 * (winmax - winmin + 1) + 1;                  //19

																			   /* variables */
	uint16_t pixLo = SEARCH_SIZE + 1;                                               //pixLo=5
	uint16_t pixHi = FRAME_SIZE - (SEARCH_SIZE + 1) - TILE_SIZE; //frame_size == 64 //pixHi=51
	uint16_t pixStep = (pixHi - pixLo) / NUM_BLOCKS + 1;         //10
	uint16_t i, j;
	uint32_t acc[8]; // subpixels
	uint16_t histx[hist_size]; // counter for x shift
	uint16_t histy[hist_size]; // counter for y shift
	int8_t  dirsx[64]; // shift directions in x
	int8_t  dirsy[64]; // shift directions in y
	uint8_t  subdirs[64]; // shift directions of best subpixels
	float meanflowx = 0.0f;
	float meanflowy = 0.0f;
	uint16_t meancount = 0;
	float histflowx = 0.0f;
	float histflowy = 0.0f;

	/* initialize with 0 */
	for (j = 0; j < hist_size; j++) { histx[j] = 0; histy[j] = 0; }

	/* iterate over all patterns
	*/
	for (j = pixLo; j < pixHi; j += pixStep)
	{
		for (i = pixLo; i < pixHi; i += pixStep)
		{
			/* test pixel if it is suitable for flow tracking */
			uint32_t diff = compute_diff(image1, i, j, FRAME_SIZE);
			if (diff < 30)  //threshold == 30
			{
				continue;
			} //找到一个梯度变化大的图像，即为适合进行光流计算的图像

			uint32_t dist = 0xFFFFFFFF; // set initial distance to "infinity"
			int8_t sumx = 0;
			int8_t sumy = 0;
			int8_t ii, jj;

			for (jj = winmin; jj <= winmax; jj++)
			{
				for (ii = winmin; ii <= winmax; ii++)
				{
					uint32_t temp_dist = compute_sad_8x8(image1, image2, i, j, i + ii, j + jj, 64);
					if (temp_dist < dist)
					{
						sumx = ii;
						sumy = jj;
						dist = temp_dist;
					}
				}
			}

			/* acceptance SAD distance threshhold */
			if (dist < 5000) //threshold==5000;即相关性可接受范围
			{
				meanflowx += (float)sumx;
				meanflowy += (float)sumy;

				compute_subpixel(image1, image2, i, j, i + sumx, j + sumy, acc, 64);
				uint32_t mindist = dist; // best SAD until now
				uint8_t mindir = 8; // direction 8 for no direction
				for (uint8_t k = 0; k < 8; k++)
				{
					if (acc[k] < mindist)
					{
						// SAD becomes better in direction k
						mindist = acc[k];
						mindir = k;
					}
				}
				dirsx[meancount] = sumx;
				dirsy[meancount] = sumy;
				subdirs[meancount] = mindir;
				meancount++;

				/* feed histogram filter*/
				uint8_t hist_index_x = 2 * sumx + (winmax - winmin + 1);
				if (subdirs[i] == 0 || subdirs[i] == 1 || subdirs[i] == 7) hist_index_x += 1;
				if (subdirs[i] == 3 || subdirs[i] == 4 || subdirs[i] == 5) hist_index_x += -1;
				uint8_t hist_index_y = 2 * sumy + (winmax - winmin + 1);
				if (subdirs[i] == 5 || subdirs[i] == 6 || subdirs[i] == 7) hist_index_y += -1;
				if (subdirs[i] == 1 || subdirs[i] == 2 || subdirs[i] == 3) hist_index_y += 1;

				histx[hist_index_x]++;
				histy[hist_index_y]++;

			}
		}
	}


	/* evaluate flow calculation */
	if (meancount > 10)
	{
		meanflowx /= meancount;
		meanflowy /= meancount;

		int16_t maxpositionx = 0;
		int16_t maxpositiony = 0;
		uint16_t maxvaluex = 0;
		uint16_t maxvaluey = 0;

		/* position of maximal histogram peek */
		for (j = 0; j < hist_size; j++)
		{
			if (histx[j] > maxvaluex)
			{
				maxvaluex = histx[j];
				maxpositionx = j;
			}
			if (histy[j] > maxvaluey)
			{
				maxvaluey = histy[j];
				maxpositiony = j;
			}
		}

		if (1)
		{

			/* use histogram filter peek value */
			uint16_t hist_x_min = maxpositionx;
			uint16_t hist_x_max = maxpositionx;
			uint16_t hist_y_min = maxpositiony;
			uint16_t hist_y_max = maxpositiony;

			/* x direction */
			if (maxpositionx > 1 && maxpositionx < hist_size - 2)
			{
				hist_x_min = maxpositionx - 2;
				hist_x_max = maxpositionx + 2;
			}
			else if (maxpositionx == 0)
			{
				hist_x_min = maxpositionx;
				hist_x_max = maxpositionx + 2;
			}
			else  if (maxpositionx == hist_size - 1)
			{
				hist_x_min = maxpositionx - 2;
				hist_x_max = maxpositionx;
			}
			else if (maxpositionx == 1)
			{
				hist_x_min = maxpositionx - 1;
				hist_x_max = maxpositionx + 2;
			}
			else  if (maxpositionx == hist_size - 2)
			{
				hist_x_min = maxpositionx - 2;
				hist_x_max = maxpositionx + 1;
			}

			/* y direction */
			if (maxpositiony > 1 && maxpositiony < hist_size - 2)
			{
				hist_y_min = maxpositiony - 2;
				hist_y_max = maxpositiony + 2;
			}
			else if (maxpositiony == 0)
			{
				hist_y_min = maxpositiony;
				hist_y_max = maxpositiony + 2;
			}
			else if (maxpositiony == hist_size - 1)
			{
				hist_y_min = maxpositiony - 2;
				hist_y_max = maxpositiony;
			}
			else if (maxpositiony == 1)
			{
				hist_y_min = maxpositiony - 1;
				hist_y_max = maxpositiony + 2;
			}
			else if (maxpositiony == hist_size - 2)
			{
				hist_y_min = maxpositiony - 2;
				hist_y_max = maxpositiony + 1;
			}

			float hist_x_value = 0.0f;
			float hist_x_weight = 0.0f;

			float hist_y_value = 0.0f;
			float hist_y_weight = 0.0f;

			for (uint8_t h = hist_x_min; h < hist_x_max + 1; h++)
			{
				hist_x_value += (float)(h*histx[h]);
				hist_x_weight += (float)histx[h];
			}

			for (uint8_t h = hist_y_min; h < hist_y_max + 1; h++)
			{
				hist_y_value += (float)(h*histy[h]);
				hist_y_weight += (float)histy[h];
			}

			histflowx = (hist_x_value / hist_x_weight - (winmax - winmin + 1)) / 2.0f;
			histflowy = (hist_y_value / hist_y_weight - (winmax - winmin + 1)) / 2.0f;

		}
		else {
			/* use average of accepted flow values */
			uint32_t meancount_x = 0;
			uint32_t meancount_y = 0;

			for (uint8_t h = 0; h < meancount; h++)
			{
				float subdirx = 0.0f;
				if (subdirs[h] == 0 || subdirs[h] == 1 || subdirs[h] == 7) subdirx = 0.5f;
				if (subdirs[h] == 3 || subdirs[h] == 4 || subdirs[h] == 5) subdirx = -0.5f;
				histflowx += (float)dirsx[h] + subdirx;
				meancount_x++;

				float subdiry = 0.0f;
				if (subdirs[h] == 5 || subdirs[h] == 6 || subdirs[h] == 7) subdiry = -0.5f;
				if (subdirs[h] == 1 || subdirs[h] == 2 || subdirs[h] == 3) subdiry = 0.5f;
				histflowy += (float)dirsy[h] + subdiry;
				meancount_y++;
			}

			histflowx /= meancount_x;
			histflowy /= meancount_y;
		}
		*pixel_flow_x = histflowx;
		*pixel_flow_y = histflowy;
	}
	else
	{
		*pixel_flow_x = 0.0f;
		*pixel_flow_y = 0.0f;
		return 0;
	}

	/* calc quality */
	uint8_t qual = (uint8_t)(meancount * 255 / (5 * 5));

	return qual; //品质由0~255.
}




static inline uint32_t compute_sad_8x8(uint8_t *image1, uint8_t *image2, uint16_t off1X, uint16_t off1Y, uint16_t off2X, uint16_t off2Y, uint16_t row_size)
{
	/* calculate position in image buffer */
	uint16_t off1 = off1Y * row_size + off1X; // image1
	uint16_t off2 = off2Y * row_size + off2X; // image2

	uint32_t acc = 0;
	int i, j;
	for (i = 0; i<8; i++)
	{
		for (j = 0; j<8; j++)
		{
			acc += abs(image1[off1 + i + j*row_size] - image2[off2 + i + j*row_size]);
		}
	}
	return acc;
}

static inline uint32_t compute_diff(uint8_t *image, uint16_t offX, uint16_t offY, uint16_t row_size)
{
	/* calculate position in image buffer */
	uint16_t off = (offY + 2) * row_size + (offX + 2); // we calc only the 4x4 pattern
	uint32_t acc = 0;
	int i, j;
	for (i = 0; i<3; i++)
	{
		for (j = 0; j < 4; j++)
		{
			acc += abs(image[off + i*row_size + j] - image[off + (i + 1)*row_size + j]);
		}
	}
	for (j = 0; j<3; j++)
	{
		for (i = 0; i<4; i++)
		{
			acc += abs(image[off + i*row_size + j] - image[off + i*row_size + j + 1]);
		}
	}
	return acc;

}

static inline uint32_t compute_subpixel(uint8_t *image1, uint8_t *image2, uint16_t off1X, uint16_t off1Y, uint16_t off2X, uint16_t off2Y, uint32_t *acc, uint16_t row_size)
{
	/* calculate position in image buffer */
	uint16_t off1 = off1Y * row_size + off1X; // image1
	uint16_t off2 = off2Y * row_size + off2X; // image2

	uint32_t s0, s1, s2, s3, s4, s5, s6, s7, t1, t3, t5, t7;

	int i, j;
	for (i = 0; i < 8; i++)
	{
		acc[i] = 0;
	}

	for (i = 0; i<8; i++)
	{
		for (j = 0; j<8; j++)
		{
			s0 = (image2[off2 + 0 + (i + 0) * row_size] + image2[off2 + 1 + (i + 0) * row_size]) / 2;
			s1 = (image2[off2 + 0 + (i + 1) * row_size] + image2[off2 + 1 + (i + 1) * row_size]) / 2;
			s2 = (image2[off2 + 0 + (i + 0) * row_size] + image2[off2 + 0 + (i + 1) * row_size]) / 2;
			s3 = (image2[off2 + 0 + (i + 1) * row_size] + image2[off2 - 1 + (i + 1) * row_size]) / 2;
			s4 = (image2[off2 + 0 + (i + 0) * row_size] + image2[off2 - 1 + (i + 0) * row_size]) / 2;
			s5 = (image2[off2 + 0 + (i - 1) * row_size] + image2[off2 - 1 + (i - 1) * row_size]) / 2;
			s6 = (image2[off2 + 0 + (i + 0) * row_size] + image2[off2 + 0 + (i - 1) * row_size]) / 2;
			s7 = (image2[off2 + 0 + (i - 1) * row_size] + image2[off2 + 1 + (i - 1) * row_size]) / 2;
			t1 = (s0 + s1) / 2;
			t3 = (s3 + s4) / 2;
			t5 = (s4 + s5) / 2;
			t7 = (s7 + s0) / 2;

			acc[0] += abs(uint8_t(image1[off1 + j + i * row_size] - s0));
			acc[1] += abs(uint8_t(image1[off1 + j + i * row_size] - t1));
			acc[2] += abs(uint8_t(image1[off1 + j + i * row_size] - s2));
			acc[3] += abs(uint8_t(image1[off1 + j + i * row_size] - t3));
			acc[4] += abs(uint8_t(image1[off1 + j + i * row_size] - s4));
			acc[5] += abs(uint8_t(image1[off1 + j + i * row_size] - t5));
			acc[6] += abs(uint8_t(image1[off1 + j + i * row_size] - s6));
			acc[7] += abs(uint8_t(image1[off1 + j + i * row_size] - t7));
		}
	}
	return 0;
}
