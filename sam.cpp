#include "sam.h"

SAM::SAM(double angle_min, double angle_increment, int nrays)
{
	creatCache(angle_min, angle_increment, nrays);
}
SAM::~SAM() {}


void SAM::doSam(const vector<float> &laserdata, int RaysAmount, 
				float L1, float L2, float angle_mark, 
				double AngleIncrease, double angle_start, LandMarkType LMTYPE,
				float sam_thresh = 0.05, float dist_thresh = 0.2, float segment_thresh = 5)
{
	
	static int i=0;
	EDGE1 = L1;
	EDGE2 = L2;
	ANGLE_LANDMARK = angle_mark;
	LMtype = LMTYPE;
	LASER_THRESHOLD = sam_thresh;
	DIST_THRESHOLD = dist_thresh;
	SEG_THRESHOLD = segment_thresh;
	//int nrays = (sizeof(laserdata) / sizeof(laserdata[0]));
	//    sam for icp
	Point *PointList = new Point[RaysAmount];
	vector<Segment> Segs;
	
	Mark = SamforIcp(laserdata, PointList, Segs, RaysAmount, AngleIncrease, angle_start);
	delete[]PointList;
	//ROS_INFO("after samfroICP");
	if (Mark.HaveMark)
	{
		//i++;
		//std::cout << i << std::endl;
		landmark_found = true;
		result_B_W_x = Mark.O_B_W.x;
		result_B_W_y = Mark.O_B_W.y;
		result_B_W_theta = Mark.Theta;
		result_W_B_x = Mark.O_W_B.x;
		result_W_B_y = Mark.O_W_B.y;
		result_W_B_theta = Mark.Theta_W_B;

		if(Mark.Theta > 50)
		{
			std::cout <<"Mark.Theta=" << Mark.Theta << std::endl;
			std::cout <<"Mark.O_B_W.x=" << Mark.O_B_W.x <<std::endl;
			std::cout <<"Mark.Serial = " <<Mark.Serial << std::endl;
		}

	}
	else
	{
		
		landmark_found = false;
		result_B_W_x = 0;
		result_B_W_y = 0;
		result_B_W_theta = 0;
		result_W_B_x = 0;
		result_W_B_y = 0;
		result_W_B_theta = 0;
	}

	DeleteVector(Segs);
}

LandMark SAM::SamforIcp(const vector<float>& ScanData, Point *PointList, vector<Segment> &Segs, int RaysAmount, double Angle_Increase, double Angle_Start) {
	//    read data and generate PointList
	//std::cout << RaysAmount << std::endl;
	//ROS_INFO("into SamforICP");
	double *Theta = new double[RaysAmount];
	for (int i = 0; i < RaysAmount; i++) {
		Theta[i] = (i*Angle_Increase + Angle_Start) / 180 * PI;

		PointList[i].x = a_cos_[i]*ScanData[i];
		PointList[i].y = a_sin_[i]*ScanData[i] ;
	}

	//    preprocessing, clustering according to d between points
	int DistIndexLength = 1;
	vector<int> DistIndex;
	DistIndex.push_back(0);
	
	//ROS_INFO("before DistanceFilter");
	DistanceFilter(PointList, RaysAmount, DistIndex, DistIndexLength);

	//    do sam, find SegIndex
	int Start, Last;
	vector<Point> ResultList; 
	for (int i = 0; i < (DistIndexLength - 1); i++) 
	{
		Start = DistIndex[i];
		Last = DistIndex[i + 1];
		Point *ScanPartial = new Point[Last - Start + 1];
		for (int j = Start; j <= Last; j++) 
		{
			*(ScanPartial + j - Start) = PointList[j];
	    }
		vector<Point> TempList;
		sam(ScanPartial, Last - Start + 1, TempList);
		ResultList.insert(ResultList.end(), TempList.begin(), TempList.end());
		DeleteVector(TempList);
		delete[]ScanPartial;
	}
	DeleteVector(DistIndex);
	vector<int> SegmentIndex;
	SegmentIndex.push_back(0);
	FindSegment(PointList, RaysAmount, ResultList, ResultList.size(), SegmentIndex);
	DeleteVector(ResultList);


	//find segment
	for (int i = 0; i < (SegmentIndex.size() - 1); i++) {
		Segment tempseg;
		if (abs(SegmentIndex[i] - SegmentIndex[i + 1]) > SEG_THRESHOLD)
		{
			// least square fitting
			tempseg.IsSegment = true;
			tempseg.Param = LSFit(PointList, SegmentIndex, i);
			tempseg.StartIndex = SegmentIndex[i];
			tempseg.EndIndex = SegmentIndex[i + 1];
			tempseg.Start_x = PointList[SegmentIndex[i]].x;
			tempseg.Start_y = PointList[SegmentIndex[i]].y;
			tempseg.End_x = PointList[SegmentIndex[i + 1]].x;
			tempseg.End_y = PointList[SegmentIndex[i + 1]].y;
			if (tempseg.Param.B == 0)
			{
				tempseg.Angle = 90;
			}
			else {
				float k = -tempseg.Param.A / tempseg.Param.B;
				tempseg.Angle = atan(k);
				tempseg.Radian2Angle();
				//if(k < 0) {
				//    tempseg.Angle += 180;
				//}
			}
		}
		else {
			tempseg.IsSegment = false;
			tempseg.StartIndex = SegmentIndex[i];
			tempseg.EndIndex = SegmentIndex[i + 1];
			tempseg.Angle = 0xFFFFFFFF;
		}
		Segs.push_back(tempseg);
	}	

	LandMark Mark;
	//std::cout << "start FindLandMark" << std::endl;
	Mark = FindMark(Segs, SegmentIndex.size() - 1, PointList);
	DeleteVector(SegmentIndex);
	if (Mark.HaveMark == false) {
		//cout << "Cannot find landmark!" << endl;
	}
	else {


		// Mark.Theta is in degree
		double the = Mark.Theta / 180 * M_PI; 
		Mark.O_W_B.x = -cos(the)*Mark.O_B_W.x - sin(the)*Mark.O_B_W.y;
		Mark.O_W_B.y = sin(the)*Mark.O_B_W.x - cos(the)*Mark.O_B_W.y;

		// Change Mark.Theta into the angle of O_W_B
		Mark.Theta_W_B = -Mark.Theta;
		start_index_ = Segs[Mark.Serial].StartIndex;
		end_index_ = Segs[Mark.Serial+1].EndIndex;

		//cout << "Find out landmark!" << endl;
		//cout << "Line" << Mark.Serial << " & " << "Line" << Mark.Serial + 1 << endl;
		//cout << "StartIndex: " << Segs[Mark.Serial].StartIndex << endl;
		//cout << "MiddleIndex" << Segs[Mark.Serial].EndIndex << endl;
		//cout << "EndIndex: " << Segs[Mark.Serial + 1].EndIndex << endl;

		//cout << "===============" << endl;
		//cout << "Theta_rad = " << the << endl;
		//cout << "O_B_W.x,y = " << Mark.O_B_W.x << '\t' << Mark.O_B_W.y << endl;
		//cout << "O_W_B.x,y = " << Mark.O_W_B.x << '\t' << Mark.O_W_B.y << endl;
		//cout << "===============" << endl;

	}

	delete[]Theta;
	return Mark;
}

void SAM::DistanceFilter(Point *Points, int PointNum, vector<int> &Index, int &IndexLength) {
	for (int i = 0; i < (PointNum - 1); i++) {
		if (norm(Points[i], Points[i + 1]) > DIST_THRESHOLD) {
			if (i == Index[IndexLength - 1]) {
				IndexLength++;
				Index.push_back(i + 1);
			}
			else {
				IndexLength += 2;
				Index.push_back(i);
				Index.push_back(i + 1);
			}
		}
	}
	IndexLength++;
	Index.push_back(PointNum - 1);
}

void SAM::sam(Point *TempList, int TempLength, vector<Point> &ResultList) {
	double DistMax = 0;
	double Dist;
	int Index = 0;
	int Last = TempLength;
	//    Find the point with the maximum distance
	for (int i = 1; i < (Last - 1); i++) {
		Dist = P2LDist(*(TempList + i), LineFit(*TempList, *(TempList + Last - 1)));
		if (Dist > DistMax) {
			Index = i;
			DistMax = Dist;
		}
	}
	//    If max distance is greater than epsilon, recursively simplify
	if (DistMax > LASER_THRESHOLD) {
		vector<Point> RecResult1;
		vector<Point> RecResult2;
		Point *TempList1 = new Point[Index + 1];
		Point *TempList2 = new Point[Last - Index];
		for (int i = 0; i <= Index; i++) {
			*(TempList1 + i) = *(TempList + i);
		}
		for (int i = Index; i < Last; i++) {
			*(TempList2 + i - Index) = *(TempList + i);
		}
		sam(TempList1, Index + 1, RecResult1);
		sam(TempList2, Last - Index, RecResult2);
		vector<Point>(RecResult1).swap(RecResult1);  //minimize the space use of vector
		vector<Point>(RecResult2).swap(RecResult2);
		ResultList.clear();
		ResultList.insert(ResultList.begin(), RecResult2.begin(), RecResult2.end());
		ResultList.insert(ResultList.begin(), RecResult1.begin(), RecResult1.end());
		vector<Point>(ResultList).swap(ResultList);
		DeleteVector(RecResult1);  //release the memory
		DeleteVector(RecResult2);
		delete[]TempList1;
		delete[]TempList2;
	}
	else {
		ResultList.clear();
		ResultList.reserve(2);
		ResultList.push_back(*TempList);
		ResultList.push_back(*(TempList + Last - 1));
	}
}

void SAM::FindSegment(Point *PointList, int PointLength, vector<Point> &ResultList, int ResultLength, vector<int> &SegIndex)
{
	vector<Point> Value;
	for (int i = 0; i < ((ResultLength)-1); i++) {
		if (ResultList[i].x - ResultList[i + 1].x == 0) {
			Value.push_back(ResultList[i]);
		}
	}
	
	int v = 0;
	for (int i = 0; i < Value.size(); i++) 
	{
		for(int k=v; k<PointLength; k++)
		{
			if((PointList+k)->x == Value[i].x && (PointList+k)->y == Value[i].y)
			{
				SegIndex.push_back(k);
				v = k+1;
				break;
			}
		}
	}

}

LandMark SAM::FindMark(vector<Segment> &Segs, int SegLength, Point *PointList)
{
	//std::cout << "size of LINE = " << Segs.size()<<std::endl;
	LandMark Mark;
	vector<LandMark> MarkList;
	for (int i = 0; i < (SegLength - 1); i++)
	{
		if (Segs[i].IsSegment && Segs[i + 1].IsSegment)
		{
			if (abs(AngleBtwSeg(Segs[i], Segs[i + 1]) - (180 - ANGLE_LANDMARK)) < 10)
			{
				//cout << "angle ok, check length" << endl;
				int start1 = Segs[i].StartIndex;
				int start2 = Segs[i + 1].StartIndex;
				int end1 = Segs[i].EndIndex;
				int end2 = Segs[i + 1].EndIndex;
				double length1 = norm(PointList[start1 + 2], PointList[end1]);
				double length2 = norm(PointList[start2], PointList[end2 - 2]);
				double error1 = abs(length1 - EDGE1);
				double error2 = abs(length2 - EDGE2);

				if (error1<0.05 && error2<0.05)
				{
					LandMark Mark_temp;
					Mark_temp.O_B_W = Intersect(Segs[i], Segs[i + 1]);
					double delta_angle = (180 - AngleBtwSeg(Segs[i], Segs[i + 1])) / 2;
					double theta_rad, theta_angle;
					double k;
					
					switch (LMtype)
					{
					case 0:
						delta_angle = (180 - AngleBtwSeg(Segs[i], Segs[i + 1])) / 2;
						cout << "CONVEX" << endl;
						break;
					case 1:
						delta_angle = 90 + AngleBtwSeg(Segs[i], Segs[i + 1]) / 2;
						//cout << "CONCAVE" << endl;
						break;
					}
					if (Segs[i].Param.B != 0)
					{
						k = -Segs[i].Param.A / Segs[i].Param.B;
						theta_rad = atan(k);
						theta_angle = theta_rad*180.0 / M_PI;

						if (theta_angle >= 0 && Segs[i].Start_x < Segs[i].End_x)
							theta_angle = theta_angle - 180;
						else if (theta_angle<0 && Segs[i].Start_x < Segs[i].End_x)
							theta_angle = theta_angle + 180;
					}
					else
					{
						if (Segs[i].Start_y > Segs[i].End_y)
							theta_angle = 90;
						else
							theta_angle = -90;
					}

					Mark_temp.Theta = theta_angle + delta_angle;
					Mark_temp.Serial = i;
					Mark_temp.HaveMark = true;
					MarkList.push_back(Mark_temp);

/*
					cout << "++++++++++FindMark++++++++++" << endl;
					cout << "i = " << i << endl;
					cout << "Segs[i].StartIndex = " << Segs[i].StartIndex << endl;
					cout << "Segs[i].EndIndex = " << Segs[i].EndIndex << endl;
					cout << "Segs[i].Param =" << Segs[i].Param.A << '\t' << Segs[i].Param.B << '\t' << Segs[i].Param.C << endl;
					cout << "MarkTemp.O_B_W.x = " << Mark_temp.O_B_W.x <<endl;
					cout << "MarkTemp.O_B_W.y = " << Mark_temp.O_B_W.y <<endl;
					cout << "MarkTemp.Theta = " << Mark_temp.Theta <<endl;
					//cout << "k = " << k << endl;
					cout << "++++++++++++++++++++++++++++" << endl;
*/
				}
			}
		}
	}
	
	int MarkListSize = MarkList.size();
	//cout << "+++++++++++FindMark+++++++++" << endl;
	//cout << "MarkListSize  " << MarkListSize << endl;
	if (!MarkListSize)
	{
		Mark.HaveMark = false;
		Mark.Theta = 0;
		Mark.Theta_W_B = 0;
		Mark.Serial = 0;
		Mark.O_B_W.x = 0;
		Mark.O_B_W.y = 0;
		Mark.O_W_B.x = 0;
		Mark.O_W_B.y = 0;

		return Mark;
	}
	else
	{
		//find out true landmark
		//Method in use: return the mark with smallest abs(theta); 
		double smallest_dist = MarkList[0].O_B_W.y*MarkList[0].O_B_W.y + MarkList[0].O_B_W.x*MarkList[0].O_B_W.x;
		double smallest_y = fabs(MarkList[0].O_B_W.y);
		int smallest_index = 0;
		for (int i = 0; i < MarkListSize; i++)
		{
			double dist = MarkList[i].O_B_W.y*MarkList[i].O_B_W.y + MarkList[i].O_B_W.x*MarkList[i].O_B_W.x;
			if (dist < smallest_dist)
			{	
				smallest_index = i;
				smallest_dist = dist;
			}
		} 
		//std::cout << "smallest_inMarkList_index=" << smallest_index << endl;

		return MarkList[smallest_index];
	}
}

double SAM::norm(Point PointA, Point PointB) {
	return sqrt((PointA.x - PointB.x)*(PointA.x - PointB.x) + (PointA.y - PointB.y)*(PointA.y - PointB.y));
}

Line SAM::LineFit(Point PointA, Point PointB) {
	double X = PointB.x - PointA.x;
	double Y = PointB.y - PointA.y;
	Line L;
	if (X == 0) {
		L.A = 1;
		L.B = 0;
		L.C = -PointA.x;
	}
	else {
		L.A = Y;
		L.B = -X;
		L.C = PointA.y*X - PointA.x*Y;
	}
	return L;
}

double SAM::P2LDist(Point P, Line L) {
	return (abs(L.A*P.x + L.B*P.y + L.C) / sqrt(L.A*L.A + L.B*L.B));
}

Line SAM::LSFit(Point *Points, vector<int> &Index, int i) {
	Line L;
	L.A = 0;
	L.B = 0;
	L.C = 0;
	float SumX = 0, SumY = 0, SumXY = 0, SumX2 = 0, b;
	int N = Index[i + 1] - Index[i] + 1;
	for (int j = Index[i]; j <= Index[i + 1]; j++) {
		SumX += Points[j].x;
		SumY += Points[j].y;
		SumXY += Points[j].x*Points[j].y;
		SumX2 += Points[j].x*Points[j].x;
	}
	b = (SumXY / N - SumX*SumY / (N*N)) / (SumX2 / N - (SumX / N)*(SumX / N));
	if (abs(b) >= 5000) {
		L.A = 1;
		L.B = 0;
		L.C = -SumX / N;
	}
	else {
		L.A = b;
		L.B = -1;
		L.C = SumY / N - b*SumX / N;
	}
	return L;
}

float SAM::AngleBtwSeg(Segment Seg1, Segment Seg2) {
	float A1 = Seg1.Param.A;
	float B1 = Seg1.Param.B;
	float A2 = Seg2.Param.A;
	float B2 = Seg2.Param.B;
	float CosA = abs(A1*A2 + B1*B2) / (sqrt(A1*A1 + B1*B1)) / (sqrt(A2*A2 + B2*B2));
	float AlphaRad = acos(CosA);
	float Alpha = AlphaRad * 180 / M_PI;
	return Alpha;
}

Point SAM::Intersect(Segment Seg1, Segment Seg2) {
	MatrixXd M1(2, 2);
	MatrixXd M2(2, 1);
	M1(0, 0) = Seg1.Param.A;
	M1(0, 1) = Seg1.Param.B;
	M1(1, 0) = Seg2.Param.A;
	M1(1, 1) = Seg2.Param.B;
	M2(0, 0) = -Seg1.Param.C;
	M2(1, 0) = -Seg2.Param.C;
	M2 = M1.inverse()*M2;
	Point P;
	P.x = M2(0, 0);
	P.y = M2(1, 0);
	return P;
}

void SAM::creatCache(double angle_min, double angle_increment, int nrays_arg)
{
	a_cos_.clear();
  	a_sin_.clear();

  	for (unsigned int i = 0; i < nrays_arg; ++i)
  {
    double angle = angle_min + i * angle_increment;
    a_cos_.push_back(cos(angle));
    a_sin_.push_back(sin(angle));
  }

}

Pose2D SAM::get_Result_W_B()
{
	Pose2D result_W_B;
	result_W_B.x = result_W_B_x;
    result_W_B.y = result_W_B_y;
    result_W_B.theta = result_W_B_theta/180.0*M_PI;

    return result_W_B;
}

Pose2D SAM::get_Result_B_W()
{
	Pose2D result_B_W;
	result_B_W.x = result_B_W_x;
    result_B_W.y = result_B_W_y;
    result_B_W.theta = result_B_W_theta;

    return result_B_W;
}

int SAM::get_start_index(){return start_index_;}
int SAM::get_end_index(){return end_index_;}

//========================================================

