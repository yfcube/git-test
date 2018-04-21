# 实现路标识别与定位的类——SAM

#### 0. 程序文件说明

类的声明在sam.h文件中， 类的实现在sam.cpp文件中。

#### 1.初始化

调用构造函数 `SAM::SAM(double, double, int)`

```c++
/* 参数说明
angle_min:激光数据开始的角度(即第一条扫描线的角度) (单位：弧度)
angle_increment: 激光数据的角分辨率 (单位：弧度)
nrays: 一帧激光数据的扫描点的个数
*/
SAM sam(double angle_min, double angle_increment, int nrays);
```

#### 2.利用当前帧的激光数据进行路标识别与定位

调用成员函数` SAM::doSam( )`

```c++
/*参数说明
laserdata		激光数据
nrays			激光数据的个数
L1,L2			三角路标两条边长
angle_mark		三角路标的角度(单位：度)
angle_increment	激光数据角分辨率(单位：弧度)
angle_min		激光数据开始的角度(单位：弧度)
LMTYPE			三角路标的类型，0代表凸形， 1代表凹形

sam_thresh, dist_thresh, segment_thresh 为路标识别算法用到的阈值，已给出默认值。调用函数时可以不用给出这三个参数
*/
void SAM::doSam(const vector<float>& laserdata, int	nrays,
				float L1, float L2, float angle_mark,
				double angle_increment, double angle_min, LandMarkType LMTYPE,
                 float sam_thresh=0.05, float dist_thres=0.2, float segment_thresh=5
				);
```

#### 3. 获取路标识别结果

执行`doSam()`	之后，可从SAM类的成员变量获取结果

```c++
//example
SAM sam(angle_min, angle_increment, nrays);
sam.doSam(...);

if(sam.landmark_found == 1)// 识别出路标
{   
    //路标相对于激光传感器的位姿
  	sam.result_B_W_x;
  	sam.result_B_W_y;
  	sam.result_B_W_theta; // rad
  
  	//激光传感器相对于路标的位姿
  	sam.result_W_B_x;
    sam.result_W_B_y;
    sam.result_W_B_theta; //rad
}
else //没有找到路标
{
  	// do something else.
}
```



