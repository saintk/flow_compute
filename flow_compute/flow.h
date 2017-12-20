#include "stdafx.h"
#include <iostream>  
#include "opencv2/opencv.hpp"  
#include "opencv2/imgproc.hpp"
#include "opencv2/highgui.hpp"
#include "cmath"

using namespace cv;
using namespace std;


uint8_t compute_flow(uint8_t *image1, uint8_t *image2, float *pixel_flow_x, float *pixel_flow_y);
static inline uint32_t compute_diff(uint8_t *image, uint16_t offX, uint16_t offY, uint16_t row_size);
static inline uint32_t compute_subpixel(uint8_t *image1, uint8_t *image2, uint16_t off1X, uint16_t off1Y, uint16_t off2X, uint16_t off2Y, uint32_t *acc, uint16_t row_size);
static inline uint32_t compute_sad_8x8(uint8_t *image1, uint8_t *image2, uint16_t off1X, uint16_t off1Y, uint16_t off2X, uint16_t off2Y, uint16_t row_size);
