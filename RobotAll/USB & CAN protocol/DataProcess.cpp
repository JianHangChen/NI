/**************************************************************************************************
Copyright, 2010-2012, Robotics Lab., Dept. of M.E., National Taiwan University
File Name: DataProcess.cpp

Author: Jiu-Lou Yan
Version: 1.0
Date: 2010/06/20

Functions:
	DataProcess() ~DataProcess() DPInit() TransENC()

Classes: DataProcess

Description:
     本程式主要資料傳輸物件
	 目前只有encoder轉碼功能 主要轉換電腦與local端的encoder值
	 各函式與變數之說明請詳見下方宣告與定義處之說明

Note: None
***************************************************************************************************/

#include "StdAfx.h"
#include "DataProcess.h"
#include "ConstIO.h"
#include <iostream>
#include "math.h"
using namespace std;
DataProcess::DataProcess(void)
{
	/******************************************************************
	input: void
	output: void

	Note:
	// Class constructor 
	******************************************************************/
}

DataProcess::~DataProcess(void)
{
	/******************************************************************
	input: void
	output: void

	Note:
	// Class destructor 
	******************************************************************/
}
void DataProcess::DPInit(void)
{
	/******************************************************************
	input: void
	output: void

	Note:
	// 初始化變數 State = [255 0 0 0 255]
	// 這是用來與C32溝通用的5 byte 變數
	******************************************************************/
		State[0]=0xff;
		State[1]=0x00;
		State[2]=0x00;
		State[3]=0x00;
		State[4]=0xff;
}

void DataProcess::TransENC(long double TarENC,unsigned int TxID)
{
	/******************************************************************
	input: TarENC 目標軸的encoder命令  TxID 目標軸的ID
	output: void

	Note:
	// 設定單軸encoder位置 控制馬達位置用
	******************************************************************/

		TxData[0]=TxID;
		if (TarENC > 0)
		{
			// 32bit資料被拆成4 bytes 才相容於CAN-Bus格式
			TxData[1] = (signed long)TarENC;
            TxData[2] = (signed long)TarENC>>8;
			TxData[3] = (signed long)TarENC>>16;
			TxData[4] = (signed long)TarENC>>24;
			TxData[5] = TxData[6] = TxData[7] = 0;	
			TxData[8] = 1;
		}
		else if (TarENC == 0)
		{
			TxData[1] = TxData[2] = TxData[3] = TxData[4] = TxData[5] = TxData[6] = TxData[7] = 0;
			TxData[8] = 1; 
		}
		else
		{
			// 32bit資料被拆成4 bytes 才相容於CAN-Bus格式
			TxData[1] = (signed long)(-TarENC);
			TxData[2] = (signed long)(-TarENC)>>8;
			TxData[3] = (signed long)(-TarENC)>>16;
			TxData[4] = (signed long)(-TarENC)>>24;
			TxData[5] = TxData[6] = TxData[7] = 0;	
			TxData[8] = 1; 	

			TxData[4] = 255-TxData[4]; // 告知對方是負的
		}

}
