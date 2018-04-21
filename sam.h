#ifndef SAM_H
#define SAM_H

#include <iostream>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include "Eigen3/Eigen/Dense"


#define PI 3.14159265

using namespace std;
using namespace Eigen;

typedef int LandMarkType;

struct Pose2D
{
	double x;
	double y;
	double theta;
};

struct Point
{
	double x;
	double y;
};
struct Line
{
	float A;
	float B;
	float C;
};
struct Segment
{

	bool IsSegment;
	Line Param;
	int StartIndex;
	int EndIndex;
	float Angle;
	double Start_x, Start_y;
	double End_x, End_y;
	void Angle2Radian() {
		Angle = Angle*PI / 180;
	}
	void Radian2Angle() {
		Angle = Angle * 180 / PI;
	}
};
struct LandMark {

	bool HaveMark;
	float Theta, Theta_W_B;
	Point O_B_W, O_W_B;
	int Serial;
};
class SAM
{
public:
	SAM(double angle_min, double angle_increment, int nrays);
	~SAM();

	bool   landmark_found;
	double result_B_W_x, result_B_W_y, result_B_W_theta;
	double result_W_B_x, result_W_B_y, result_W_B_theta;
	double nrays_public;

	//method
	//void doSam(double *laserdata, double L1, double L2);
	void doSam(const vector<float>& laserdata, int RaysAmount, 
				float L1, float L2, float angle_mark, 
				double AngleIncrease, double angle_start, LandMarkType LMTYPE,
				float sam_thresh, float dist_thresh, float segment_thresh);
	void creatCache(double angle_min, double angle_increment,int nrays_arg);
	
	Pose2D get_Result_W_B();
	Pose2D get_Result_B_W();
	int get_start_index();
	int get_end_index();

private:

	//int RAYSAMOUNT;
	int NUM_SCAN;
	int RAYSTART;
	double DEG_INCREASE;
	int GROUP;
	double DIST_THRESHOLD;
	double LASER_THRESHOLD;
	int SEG_THRESHOLD;
	double ANGLE_START;
	LandMark Mark;

	double nrays_private;
	vector<double> a_sin_, a_cos_;
	double EDGE1, EDGE2;
	double ANGLE_LANDMARK;
	LandMarkType LMtype;
	int start_index_;
	int end_index_;

	double angle_min_;
	double angle_increment_;

	//method
	LandMark SamforIcp(const vector<float>& ScanData, Point *PointList, vector<Segment>& Segs, int RaysAmount = 541, double AngleIncrease = 0.5, double AngleStart = -135.0);
	void DistanceFilter(Point *Points, int PointNum, vector<int> &Index, int &IndexLength);
	void FindSegment(Point *PointList, int PointLength, vector<Point> &ResultList, int ResultLength, vector<int> &SegIndex);
	void sam(Point *TempList, int TempLength, vector<Point> &ResultList);
	LandMark FindMark(vector<Segment> &Segs, int SegLength, Point *PointList);
	double norm(Point PointA, Point PointB);
	Line LineFit(Point PointA, Point PointB);
	double P2LDist(Point P, Line L);
	Line LSFit(Point *Points, vector<int> &Index, int i);
	float AngleBtwSeg(Segment Seg1, Segment Seg2);
	Point Intersect(Segment Seg1, Segment Seg2);

	template <class T>
	int length(T& Array) {
		return (sizeof(Array) / sizeof(Array[0]));
	}

	template <class T>
	void DeleteVector(vector<T> &v) {
		vector<T> Tempv;
		Tempv.swap(v);
	}

};

#endif