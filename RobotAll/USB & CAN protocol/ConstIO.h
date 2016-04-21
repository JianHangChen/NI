/**************************************************************************************************
Copyright, 2010-2012, Robotics Lab., Dept. of M.E., National Taiwan University
File Name: ConstIO.h

Author: Jiu-Lou Yan
Version: 1.0
Date: 2012/02/01

Functions:None

Classes: None

Description:
     本標頭檔定義了所有CAN-Bus通訊的指令定義以及通訊協定
	 詳細內容與定義請參照內文註解

Note: None
***************************************************************************************************/

// ID parameters
#define MID 0 // 個軸的ID 
#define CMD_SET 16 // 每軸

// 說明請參照下方舊版
enum ID_Sets{ID_InitMech,ID_InitPara,ID_SetENC,ID_SetPID,ID_ModeIdle,ID_ModePID,ID_AllReturnENC,ReservedID1,ID_SetOneAxis,
	         ID_SetIdle,ID_SetCaliENC,ID_ReadEncoder,ReservedID2,ID_SetPWMLimit,ID_EndOfContTraj,ID_SetTrajPacket};

#define ID_SendReadENC		MID*CMD_SET+1
#define ID_ReadENCcount		MID*CMD_SET+9
#define ID_SetpotPID		MID*CMD_SET+11

//// 舊版
//// Master Output ID Command Definitions
//#define ID_InitMech		MID*CMD_SET+0 // 機構初始化歸零
//#define ID_InitPara		MID*CMD_SET+1 // 初始化參數 目前不需要使用
//#define ID_SetENC		MID*CMD_SET+2 // 設定單軸目標ENCDOR位置
//#define ID_SetPID		MID*CMD_SET+3 // 設定單軸PID參數
//#define ID_ModeIdle		MID*CMD_SET+4 // 設定單軸為IDLE MODE
//#define ID_ModePID		MID*CMD_SET+5 // 設定單軸為PID MODE
//#define ID_AllReturnENC MID*CMD_SET+6 // all motor nodes return their ENC counts
//#define ID_ReadEncoder  MID*CMD_SET+11 //讀取單軸encoder
////////////////////////  MID*CMD_SET+7 7號已經被定義為回傳 //////////////////////
//#define ID_SetOneAxis	MID*CMD_SET+8 // 設定單軸命令 -- 跟C32溝通專用
//#define ID_SetIdle		MID*CMD_SET+9 // 設定IDLE MODE -- 跟C32溝通專用
//#define ID_SetCaliENC 	MID*CMD_SET+10 // 設定機構誤差補償之encoder值 
//#define ID_SetPWMLimit		13 // 設PWM的上下限值
//#define ID_EndOfContTraj	14 // 告訴C32連續軌跡已經傳完
//#define ID_SetTrajPacket	15 // 告訴C32即將開始連續傳送大量軌跡

// 說明請看下方舊版
enum SP_Sets{SP_IDLE,SP_PID,SP_INIT_MECH,SP_WATCH_NODES,SP_SetOneAxis,SP_SetIdle};
//// State Machine -- Primary State Definition
//#define SP_IDLE			0 // main state is idle
//#define SP_PID			1 // main state is pid
//#define SP_INIT_MECH	2 // main state is initial mechanism
//#define SP_WATCH_NODES  3 // the master node watches the slave nodes -- not used now
//#define SP_SetOneAxis	4 // set the position of one specified axis 
//#define SP_SetIdle		5 // set C32 IDLE

//#define PI 3.1415926  // 在Global Constant 設定

