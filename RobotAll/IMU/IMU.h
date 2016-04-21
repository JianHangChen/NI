#include <stdio.h>
#include <windows.h>
#include <sys/timeb.h>
#include <time.h>
#include <fstream>
#include <iostream>



#define PI 3.1415926535828 
#define IMUbuffersize 10000
using namespace std;
extern "C"       //include c library 特殊寫法

{
#include "i3dmgx3_Errors.h"
#include "i3dmgx3_Cont.h"
#include "i3dmgx3_Utils.h"
#include "i3dmgx3_Utils_CM.h"
#include "i3dmgx3_readWriteDR.h"
#include "i3dmgx3_Serial.h"
}

#ifndef _IMU_H_
#define _IMU_H_
#pragma once



class IMU {

public:
	
	IMU (); // 建構子
	~IMU (); // 解構子


public:	

	
	int absangle ; //  flag 絕對角度
	int count   ;
	int datacount  ;    //存下的data數
	int calicount  ;   //cali的筆數
	int adamssimcount ; // adams counts start from 1000
	double Ts ;
	float  ang1,ang2; 
	float microstraincalix , microstraincaliy  ;  //cali值  與 microstrain MIP 比較
	





	float anglex[IMUbuffersize];        // vel 積分
	float angley[IMUbuffersize];
	float absanglex [IMUbuffersize] ;  //abs 絕對角度
	float absangley [IMUbuffersize] ;
	
	float finalanglex[IMUbuffersize];   //機器人pitch方向
	float finalangley[IMUbuffersize];   //機器人roll方向
	float finalanglez[IMUbuffersize];
	






	float anglefilterx[IMUbuffersize];
	float anglefiltery[IMUbuffersize];
	float anglefilterz[IMUbuffersize];//20140921 dora add
	float anglez[IMUbuffersize];
	float filterx [IMUbuffersize];
	float filtery [IMUbuffersize];
	float filterabsx [IMUbuffersize];
	float filterabsy [IMUbuffersize];
		
	float finalanglex_test[IMUbuffersize];

	float velx [IMUbuffersize];	
	float vely [IMUbuffersize];	
	float velz [IMUbuffersize];	
	
	float accelx[IMUbuffersize];
	float accely[IMUbuffersize];
	float accelz[IMUbuffersize];

	float angvelx[IMUbuffersize];//20140921 DORA add
	float angvely[IMUbuffersize];//20140921 DORA add
	float angvelz[IMUbuffersize];//20140921 DORA add
	float posx[IMUbuffersize];//20140921 DORA add
	float posy[IMUbuffersize];//20140921 DORA add
	float posz[IMUbuffersize];//20140921 DORA add
	
	float filteraccelx[IMUbuffersize];
	float filteraccely[IMUbuffersize];
	float filteraccelz[IMUbuffersize];
   


	double waisttheta[IMUbuffersize];
	double waistthetadot[IMUbuffersize];


	//------------------------variables---------------------------//
	float  velbiasx , velbiasx1 , velbiasx2 ;
	float  velbiasy , velbiasy1 , velbiasy2 ;
	float  velbiasz , velbiasz1 , velbiasz2 ;// vel cali量


	float  anglebiasx ; 
	float  anglebiasy ;
	float  anglebiasz ; 

	float anglexbias  ;
	float angleybias  ; 
	float absxbias    ;
	float absybias    ;

	float abscalix    ;
	float abscaliy    ;
	float anglecalix  ;  
	float anglecaliy  ;  
	float anglecaliz  ;  

    ////// FOR IMU rotation matrix 
	float IMU_roll;
	float IMU_pitch;
	float IMU_yaw;


//------------------------------------------------ 相關function-----------------------------------------------------------//
 void ReadContinuousData(int portNum); 
 void integation(float *velx ,float *anglex , float *vely ,float *angley ,float *velz ,float *anglez );
 void complementaryfilter(float* ang ,  float * abs , float * filterang , float *filterabs , float *filter , float bias  );
 void adamssim(int gFlagSimulation,int adamssimcount);
 bool IMU_Lock;
 //------------------------------------------------ 相關function-----------------------------------------------------------//




}  ;

#endif

