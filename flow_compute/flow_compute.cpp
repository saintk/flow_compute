// flow_compute.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include <iostream>  
#include "opencv2/opencv.hpp"  
#include "opencv2/imgproc.hpp"
#include "opencv2/highgui.hpp"
#include "cmath"
#include "flow.h"

using namespace cv;
using namespace std;


int main()
{
	VideoCapture cap(0);
	if (!cap.isOpened())
	{
		return -1;
	}
	Mat frame;
	Mat edge;
	Mat A;
	Mat frame_now(64,64,CV_8UC1), frame_last(64, 64, CV_8UC1);
	uint8_t*p_now = frame_now.data;
	uint8_t*p_last = frame_last.data;

	float pixel_flow_x = 0.0f;
	float pixel_flow_y = 0.0f;
	float pixel_flow_x_sum = 0.0f;
	float pixel_flow_y_sum = 0.0f;
	float velocity_x_sum = 0.0f;
	float velocity_y_sum = 0.0f;
	float velocity_x_lp = 0.0f;
	float velocity_y_lp = 0.0f;
	int valid_frame_count = 0;
	int pixel_flow_count = 0;

	static float accumulated_flow_x = 0;
	static float accumulated_flow_y = 0;
	static float accumulated_gyro_x = 0;
	static float accumulated_gyro_y = 0;
	static float accumulated_gyro_z = 0;
	static uint16_t accumulated_framecount = 0;
	static uint16_t accumulated_quality = 0;
	static uint32_t integration_timespan = 0;
	static uint32_t lasttime = 0;
	uint32_t time_since_last_sonar_update = 0;

	int num = 1, qual;
	uint32_t counter = 0;

	bool stop = false;

	while (!stop)
	{
		const float focal_length_px = 24.0/(20+4)*1000; //original focal (16.0mm is focal_len)
		double t = (double)getTickCount();
		frame_now.copyTo(frame_last);//取上一帧图像
		cap >> frame;
		cvtColor(frame, edge, CV_BGR2GRAY);
		resize(edge, A, Size(64, 64));
		A.copyTo(frame_now);;
		if (num == 1) { frame_now.copyTo(frame_last); num = 0; }//第一帧等于相等
		qual = compute_flow(p_last, p_now, &pixel_flow_x, &pixel_flow_y);
		cout << "pixel_flow_x = " << pixel_flow_x << ' ' << "pixel_flow_y = " << pixel_flow_y << endl;
		imshow("now", frame_now);
		t = ((double)getTickCount() - t) / getTickFrequency();

		float flow_compx = pixel_flow_x / focal_length_px / t ;
		float flow_compy = pixel_flow_y / focal_length_px / t ;

		float new_velocity_x = -flow_compx * 1.0;
		float new_velocity_y = -flow_compy * 1.0; //假定距离为1m,distance.

		if (qual > 0)
		{
			velocity_x_sum += new_velocity_x;
			velocity_y_sum += new_velocity_y;
			valid_frame_count++;

			uint32_t deltatime = t;
			integration_timespan += deltatime;
			accumulated_flow_x += pixel_flow_y / focal_length_px * 1.0f; //rad axis swapped to align x flow around y axis
			accumulated_flow_y += pixel_flow_x / focal_length_px * -1.0f;//rad
/*			accumulated_gyro_x += x_rate * deltatime / 1000000.0f;	//rad
			accumulated_gyro_y += y_rate * deltatime / 1000000.0f;	//rad
			accumulated_gyro_z += z_rate * deltatime / 1000000.0f;	//rad*/
			accumulated_framecount++;
			accumulated_quality += qual;

			/* lowpass velocity output */
			velocity_x_lp = 0.3 * new_velocity_x +
				(1.0f - 0.3) * velocity_x_lp;
			velocity_y_lp = 0.3 * new_velocity_y +
				(1.0f - 0.3) * velocity_y_lp;
		}
		else
		{
			/* taking flow as zero */
			velocity_x_lp = (1.0f - 0.3) * velocity_x_lp;
			velocity_y_lp = (1.0f - 0.3) * velocity_y_lp;
		}

		pixel_flow_x_sum += pixel_flow_x;
		pixel_flow_y_sum += pixel_flow_y;
		pixel_flow_count++;

		counter++;

//		cout<<"Vx = "<<new_velocity_x << "  Vy = " << new_velocity_y<<endl;
		if (waitKey(30)>=0)
		{
			stop = true;
		}
	}
}
/*	uint8_t*image1 = inputImage.data;
	uint8_t*image2 = transImage.data;
	float pixel_flow_x = 0.0;
	float pixel_flow_y = 0.0;
	uchar qual;
	qual = compute_flow(image1, image2, 0, 0, 0, &pixel_flow_x, &pixel_flow_y);
	cout << "pixel_flow_x = " << pixel_flow_x << endl << "pixel_flow_y = " << pixel_flow_y << endl;
	cout << "quality is " << (int)qual;
	waitKey();
}

*/

