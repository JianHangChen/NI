/**************************************************************************************************
Copyright, 2013-2015, Robotics Lab., Dept. of M.E., National Taiwan University
File Name: TwinCAT_COM.h

Author: Kenneth Yi-wen Chao, Dora
Version: 0.85
Date: 2013/07/13 by Slongz

Functions:
     請詳見 Class-TwinCAT_COM

Classes: TwinCAT_COM

Description:
     本函式庫包含了對Maxon Epos3馬達控制器進行控制/資料存取的函式
	 透過EtherCAT和TwinCAT進行連結並利用TwinCAT® I/O 和EPOS3進行資料的傳輸或存取

Note: 
	 修改建議:日後需要進行大量修改或次批次的傳值取值時,可先建立指標陣列將相關I/O的關係於建構子指定
	 可使敘述更整潔
	 Class下方有預列出可能可以增加的Func名稱
***************************************************************************************************/

#include "stdafx.h"  // 註解順序可能會影響compile成功與否 請注意
#include "windows.h" 
#include "stdio.h" 
#include <iostream>
#include "fstream"

#if TwinCAT_Mode        // Defined in Mainlop.h
#include <TCatIoApi.h> 	// header file shipped with TwinCAT® I/O 

#ifndef TASK_RLEG_H 
#define TASK_RLEG_H 
#include "Task_RLeg.h"  //右腳I/O定義檔,由TwinCAT匯出而得 
#endif    

#ifndef TASK_LLEG_H
#define TASK_LLEG_H 
#include "Task_LLeg.h"  //左腳I/O定義檔,由TwinCAT匯出而得 
#endif    

using namespace std;


#pragma once

class TwinCAT_COM
{
public:
	
	TwinCAT_COM(void);
	~TwinCAT_COM(void);
	
	//TwinCAT R3IO 範例程式設定值
	#define IOERR_IOSTATEBUSY 			0x2075 
	#define TASK1MS_DELAY 			1 // ms
	#define TASK10MS_DELAY 			10 // ms
	#define DELAY_IN_100NS(x) 			(x*10000)

	PTask_RLeg_Outputs 		pTRL1msOut; //Task_RLeg.h內的資料結構
	PTask_RLeg_Inputs 		pTRL1msIn; 

	PTask_LLeg_Outputs 		pTLL1msOut; //Task_LLeg.h內的資料結構
	PTask_LLeg_Inputs 		pTLL1msIn; 

	//Funcs of Motor Control and Data acquisition
	void EtherCATInit(void);					// TwinCAT初始化

	void EtherCATHomingRLeg(void);				// 右腳homing
	void EtherCATSetCSPRLeg(void);				// Set右腳為CSP Mode
	void EtherCATHomingLLeg(void);				// 左腳homing
	void EtherCATSetCSPLLeg(void);			  	// Set左腳為CSP Mode
	void EtherCATReadEncoder(long* buf);		// 讀兩腳Encoder值
	 
	void EtherCATSetEncoder(long* buf);			// 設定兩腳要追的Encoder Command (CSP)
	void EtherCATReadVel(long* buf);			// 讀兩腳轉速值(RPM)
	void EtherCATSetVelOffset(long* buf);		// 設定Velocity offset (CSP)
	void EtherCATSetTorqueOffset(short* buf);	// 設定Torqueoffset (CSP)
	void EtherCATReadTorque(short* buf);		// 讀兩腳扭矩值(千分之Rated torque)
	void EtherCATReadEncDiff(long* buf);		// 讀兩腳Encoder error值
	void EtherCATEmergentStop(void);			// 緊急停止-Disable任何模式的控制


	unsigned short* pCtrlWords[12];				// 兩腳Ctrlword 的指標陣列

	//Signal Processing
	void MovingAve(int numb, long* buf, long* result);
												//Moving Average
	static const int MovingIndex=20;			//The size of window for moving average
	long MovingBuffer[12*MovingIndex];			//The buffer of window for moving average

	//Error Handling/Check
	void EtherCATFaultReset(int index);			// Reset Epos3 (可手動reset的錯誤限定)
	bool BitCheck(unsigned short statword,bool value, int bitnumb);
												// 檢查特定bit位址之二近位值
	bool DisabledCheck(unsigned short statword);// 確認Epos3狀態是否為Disabled狀態(若是disabled則可接受initialization)
	unsigned char EPOS_ErrorCheck(unsigned short statword);
												// 確認Epos3使否產生Error及其種類
	
	//預定
	
	//void EtherCATResetInit 
	//EtherCATSetCSTLLeg
	//EtherCATSetCSTRLeg
};

#endif
