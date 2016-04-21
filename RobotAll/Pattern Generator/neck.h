/**************************************************************************************************
Copyright, 2010-2012, Robotics Lab., Dept. of M.E., National Taiwan University
File Name: neck.h

Author: che-hsuan chang
Version: 1.0
Date: 2013/03/19 
Functions:
	    setvalue()
Classes: 
        neck
Description:
     本程式主要用在計算給予脖子機構兩顆馬達的位置軌跡
	 各函式與變數之說明請詳見下方宣告與定義處之說明

Note: None
***************************************************************************************************/
#include "stdafx.h"
#include <stdio.h>
#include <iostream>
using namespace std;

class  neck 
{
public  :
	neck(void) ;  // 建構子
	~neck(void) ; //解構子
	void setvalue(int thetastart , int thetafinal  , int neck_point , int neck_motor , int neck_omega);
	int  theta;   // 角度的軌跡值
    int  thetastart; // 初始角度
    int  thetafinal; // 最終角度 
    int  neck_point;  // 初始角度與最終角度間的內差點數
	int  neck_motor  ;  // 控制脖子機構的馬達index  上馬達1  下馬達0
	int  neck_omega ;   //轉動速度  300~ 500  一般速度300  
	};

