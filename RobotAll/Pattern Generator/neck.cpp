/**************************************************************************************************
Copyright, 2010-2012, Robotics Lab., Dept. of M.E., National Taiwan University
File Name: neck.h

Author: che-hsuan chang
Version: 1
Date: 2013/03/19 
Functions:
	    setvalue()
Classes: 
        neck
Description:
     本程式主要用在計算給予脖子機構兩顆馬達的位置軌跡
	 各函式與變數之說明請詳見下方宣告與定義處之說明

Note: None

//脖子COM57
//臉 COM56

***************************************************************************************************/

#include "stdafx.h"
#include "serial_port.h"

#ifndef NECK_H
#define NECK_H
#include "neck.h"
#endif

//extern SerialPort *gpPortneck ; 


neck :: neck (void) {
      theta = 0    ;
      thetastart = 0 ; 
      thetafinal = 0  ;   
      neck_point = 0 ;   // thetastart  到  thetafinal  內差點數 
	  neck_omega = 0  ;  //轉動速度  300~ 500  一般速度300  
	  neck_motor = 0    ;//上馬達1  下馬達0
      extern SerialPort *gneck; 
      
}




void neck :: setvalue(int  thetastart  , int thetafinal , int neck_point  ,int  neck_motor , int neck_omega){
	/******************************************************************
	input:  (int  thetastart  , int thetafinal , int neck_point  ,int  neck_motor , int neck_omega)
	output: void

	Note:
	
	******************************************************************/
	
    extern SerialPort *gneck; 	
	
	for (int yy = 0 ;  yy<= neck_point ; yy++)
	{   
		
		theta  = thetastart + (((thetafinal - thetastart)/ neck_point) * yy)     ;
	    
	   
	char s[40];
	sprintf(s,"# %d P %d  S %d   \r  ", neck_motor , theta , neck_omega );
	
	
	printf("%d %d  \n "  , theta , neck_motor);
	
	unsigned char buffer[40]; 
	for( int i=0;i<strlen(s);i++)
		buffer[i]= s[i];
	gneck ->write(buffer,strlen(s));
	
		
}
}

neck ::~neck (void) {
}