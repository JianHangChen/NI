#include "plane.h"
#include <cmath>
#include <iostream>
using namespace std;


cPlaneEKF::cPlaneEKF() 
{  	 

    setDim(3, 1, 3, 3, 3);  //設定matrix dimention  n state = 4 (X = A X   X vector 的大小)   m (measurement vector  的大小)
	/*設定matrix dimention
		setDim setDim  (  
		K_UINT_32  n_,    # of state      
		K_UINT_32  nu_,   # of input 
		K_UINT_32  nw_,   # of standrd process 
		K_UINT_32  m_,    # of measure  
		K_UINT_32  nv_    # of measurement   noise 

    K_UINT_32 n;        //!< Size of the state vector.
    K_UINT_32 nu;       //!< Size of the input vector.
    K_UINT_32 nw;       //!< Size of the process noise vector.
    K_UINT_32  m;       //!< Size of the measurement vector.
    K_UINT_32 nv;       //!< Size of the measurement noise vector.


*/
	   COG_height = 0.666 ; // COG 離地高度 (單位mm)
	   samplingtime = 0.005;    //step time 5ms 
	   ww = sqrt( 9.8/ COG_height) ;
	   mass = 730; //機器人質量 單位 kg  from adams z-force
	  

}



void cPlaneEKF::makeBaseA()
{
    A(1,1) = 1.0;
	A(1,2) = samplingtime;
	A(1,3) = 0.0;

	A(2,1) = ww*ww*samplingtime;
	A(2,2) = 1;
	A(2,3) = -ww*ww*samplingtime;
	
	A(3,1) = 0.0;
	A(3,2) = 0.0;
	A(3,3) = 1.0;


	}


void cPlaneEKF::makeA()
{
	// A(1,1) = 1.0;
	//A(1,2) = samplingtime;
	//// A(1,3) = 0.0;
	//// A(1,4) = 0.0;

	//// A(2,1) = 0.0;
	//A(2,2) = 1;
	//// A(2,3) = 0.0;
	//// A(2,4) = 0.0;

	//// A(3,1) = 0.0;
	//A(3,2) = 0;
	// A(3,3) = 1.0;
	// A(3,4) = Period;
	

}



void cPlaneEKF::makeBaseW()   //根據f 做partial w 為process noise 每個state 都有process noise 
{
	W(1,1) = 1.0; 
	W(1,2) = 0.0;
	W(1,3) = 0.0;
	
	W(2,1) = 0.0;
	W(2,2) = 1.0;
	W(2,3) = 0.0;
		
	W(3,1) = 0.0;
	W(3,2) = 0.0;
	W(3,3) = 1.0;





}



void cPlaneEKF::makeBaseQ()  //standard  process covariance matrix 為 4X4 描述各個variable間的covariance 
{     
	//自己跟自己 0.001 
	Q(1,1) = 0.001;          //將x(4) external force 調大  讓其相信measurement
	Q(1,2) = 0.001/50;
	Q(1,3) = 0.001/50;
	
	Q(2,1) = 0.001/50;
	Q(2,2) = 0.001;
	Q(2,3) = 0.001/50;
	
	Q(3,1) = 0.001/50;
	Q(3,2) = 0.001/50;
	Q(3,3) = 0.001;



}





void cPlaneEKF::makeBaseH()  ///matrix C 
{
	H(1,1) = 1.0;
	H(1,2) = 0.0;
	H(1,3) = 0.0;
	
	H(2,1) = 0.0;
	H(2,2) = 0.0;
	H(2,3) = 1.0;
	
	H(3,1) = ww*ww;
	H(3,2) = 0.0;
	H(3,3) = -ww*ww;
	


}







void cPlaneEKF::makeH()
{
	/*H(1,1) = 1.0;
	H(1,2) = 0.0;
	H(1,3) = 0.0;
	
	H(2,1) = 0.0;
	H(2,2) = 0.0;
	H(2,3) = 1.0;
	
	H(3,1) = ww*ww;
	H(3,2) = 0.0;
	H(3,3) = -ww*ww;*/
}



void cPlaneEKF::makeBaseV()
{
	V(1,1) = 1.0;
	V(1,2) = 0.0;
	V(1,3) = 0.0;
	
	
	V(2,1) = 0.0;
	V(2,2) = 1.0;
    V(2,3) = 0.0;
	

	V(3,1) = 0.0;
	V(3,2) = 0.0;
    V(3,3) = 1.0;
	
}

void cPlaneEKF::makeBaseR()  //measurement noise  (依據每個sensor會有不同) 測量方式:將sensor平放，獲得covariance
{

	////EXP MODE 參數
	R(1,1) = 0.0001;  // 由joint encoder 回算的COG 位置
	R(2,2) = 0.001223;  // 由力規算出的 ZMP
	R(3,3) = 0.00022;  //IMU 量測的 COG 加速度 
	

	//process noise = 0.001
	//ADAMS MODE 參數
	//R(1,1) = 0.01;  // 由joint encoder 回算的COG 位置
	//R(2,2) = 0.0001;  // 由力規算出的 ZMP
	//R(3,3) = 0.0001;  //IMU 量測的 COG 加速度 
	


}


void cPlaneEKF::makeProcess() //  doing model  描述state間的關係 依據matrix A   
{
	Vector x_(x.size());
	
	x_(1) = x(1) +  samplingtime*x(2);
	x_(2) = ww*ww*samplingtime*x(1)+x(2)-ww*ww*samplingtime*x(3)+(samplingtime/mass)* x(4);
	x_(3) = x(3); 
	x.swap(x_);
}


void cPlaneEKF::makeMeasure()    //描述measurement與state間的關係  依據 matrix C  x(4) 為 (1/kg)*N   1/mass*x(4)*1000 要乘1000 
{
	z(1)=x(1);
	z(2)=x(3);
	z(3)= ww*ww*x(1) - ww*ww*x(3) ;
	
}



