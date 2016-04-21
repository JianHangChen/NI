/**************************************************************************************************
Copyright, 2010-2012, Robotics Lab., Dept. of M.E., National Taiwan University
File Name: DataProcess.h

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

#pragma once

#include <iostream>
using namespace std;
class DataProcess
{
public:
	DataProcess(void); // 建構子
public:
	~DataProcess(void); // 解構子

public:
	void DPInit(void); // 初始化物件
	unsigned char State[5];
	unsigned char TxData[9];
	void TransENC(long double TarENC,unsigned int TXID); // 傳輸encoder值 以及發送對象ID

};
