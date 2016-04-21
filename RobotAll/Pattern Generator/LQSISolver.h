/**************************************************************************************************
Copyright, 2010-2012, Robotics Lab., Dept. of M.E., National Taiwan University
File Name: LQSISolver.h

Author: Jiu-Lou Yan
Version: 1.1
Date: 2010/07/20

Functions:
     DiffEqn() DiffEqnShift() LowPassFilter() GetMTSAcceleration()
     MovingAvergeSmooth() LQSISolver()  ~LQSISolver() Initval()
	 BackRiccati() C2DmHumanoid() DummyControl()
	 ZMPFeedbackControl() tic()toc()

Classes: LQSISolver

Description:
     本程式主要用在解出基於倒單擺模型的inverted pendulum 問題
	 輸入ZMP以及COG高度軌跡 就可以解出水平方向ZMP擺動之軌跡
	 各函式與變數之說明請詳見下方宣告與定義處之說明

Note: None
***************************************************************************************************/

#include "Kine.h"
extern Kine gKineAll;//20140921DORA
// Optimzal control: weightings in the performance index
#define Q 100//100 // 10^6
#define R 1 // 10^-3
#define P 100//100 // 10^6

void DiffEqn(double* A, int* LenA, double* result, double* d_t); // 一階微分 回傳少了左二右二值的陣列 (微分運算捨棄最左最右兩個元素)
void DiffEqnShift(double* A, int* LenA, double* result, double* d_t); // 一階微分 並且回傳跟原陣列一樣長度的陣列
void LowPassFilter(double* data, int dataLen, double* result, double gain); // 數位低通濾波器
void GetMTSAcceleration(double* sig, double* h, double* q, int InputTrajLen, double* Boundary, double* ddQ); // 計算MTS參數
void MovingAvergeSmooth(double* Traj, int LengthOfTraj, int FrameSize, double* result); // moving average 平均法 平滑曲線用

class LQSISolver
{
public:

	LQSISolver(void); // 建構子
	~LQSISolver(void); // 解構子

	void Initval(double* InputCOGz); // 輸入重心高度軌跡 並且做初始運算
	void BackRiccati(double* ZMPx, double* ZMPy); // 輸入ZMP軌跡 並且解 riccati equation
	void BackRiccati2(double* ZMPx, double* ZMPy, int period, int index_now); // 輸入ZMP軌跡 並且解 riccati equation 
	void BackRiccati3(double* ZMPx, double* ZMPy, int period, int index_now); // 20141117修正後的推導版本，輸入ZMP軌跡 並且解 riccati equation 
	void C2DmHumanoid(int i); // 將倒單擺模型的 state matrices 做 zero order hold
	void DummyControl(void); // 將預估值直接當作回授值 預先疊代出一系列水平COG軌跡 將來用以當作機器人重心reference trajectory
	void DummyControl2(int index_now); // 將預估值直接當作回授值 預先疊代出一系列水平COG軌跡 將來用以當作機器人重心reference trajectory
	void DummyControl3(int index_now); // 將預估值直接當作回授值 預先疊代出一系列水平COG軌跡 將來用以當作機器人重心reference trajectory //20141117修正後的推導版本
	void ZMPFeedbackControl(double* ZMPx, double* ZMPy,double* COGx ,double* COGy , int index);// 計算每個iteration的重心軌跡
	void tic(void); // 計時用函數 計時開始
	void toc(void); // 計時用函數 計時結束並且在command window 顯示此次計時時間
	void toc2(void); // 計時用函數 計時結束並且在command window 顯示最大計時時間(重複計時的時候)
	void toc3(void); // 計時用函數 計時結束寫入txt 紀錄每個Step 的時間
	void ticPMS(void); // PMS計時用函數 計時開始
	void toc3PMS(void); // PMS計時用函數 計時結束寫入txt 紀錄每個Step 的時間


	YMatLite* Buffer0; // buffer for computation

	// LQSI 使用運算矩陣與狀態矩陣定義
	YMatLite* Qx;
	double Dc, Dd;

	YMatLite* Cd;
	YMatLite* CdT;

	YMatLite AdStk[LQSIBufferSize];
	YMatLite AIStk[LQSIBufferSize];
	YMatLite BdStk[LQSIBufferSize];
	YMatLite BinvMBT_Stk[LQSIBufferSize];
	YMatLite WkStk[LQSIBufferSize];
	YMatLite WkTStk[LQSIBufferSize];
	YMatLite invWkStk[LQSIBufferSize];
	YMatLite NkStk[LQSIBufferSize];
	YMatLite EkStk[LQSIBufferSize];
	YMatLite NkTStk[LQSIBufferSize];
	YMatLite SsStk[LQSIBufferSize+1];
	YMatLite vsStkX[LQSIBufferSize+1];
	YMatLite vsStkY[LQSIBufferSize+1];
	YMatLite invDkStk[LQSIBufferSize];
	YMatLite XState[LQSIBufferSize+1];
	YMatLite YState[LQSIBufferSize+1];
	YMatLite Z_Schur[LQSIBufferSize];
	YMatLite Z_Schur_U[LQSIBufferSize];
	double invMkStk[LQSIBufferSize];
	YMatLite DkStk[LQSIBufferSize/*2001*/];//max preview length

	//20140823 now state
	YMatLite AdStk2[LQSIBufferSize];//for now state
	YMatLite BdStk2[LQSIBufferSize];
	YMatLite XState2[LQSIBufferSize+1];
	YMatLite YState2[LQSIBufferSize+1];
	YMatLite Sk_test[10000][600];//20141118 test for Sk 太大會爆
	double AngAccX[LQSIBufferSize];
	double AngAccY[LQSIBufferSize];
	double AngAccZ[LQSIBufferSize];
	double AngVelX[LQSIBufferSize];
	double AngVelY[LQSIBufferSize];
	double AngVelZ[LQSIBufferSize];
	//
	// timer 計時使用變數
	LARGE_INTEGER gStartTime;
	LARGE_INTEGER gCurrentTime;
	LARGE_INTEGER gnFreq;
	float gFreqT;

	LARGE_INTEGER gStartTime_PMS;
	LARGE_INTEGER gCurrentTime_PMS;
	LARGE_INTEGER gnFreq_PMS;
	float gFreqT_PMS;
	// timer


	double Wp, sh, ch; // sh = sinh ch = cosh
	double* COGz; // 重心高度軌跡
	double delCOG[LQSIBufferSize]; // 重心高度速度
	double ddelCOG[LQSIBufferSize]; // 重心高度加速度

	double* CompuTemp; // 矩陣乘法等需要重複利用的function需要兩個暫存區 不然會覆寫到正在用的資料
	double* CompuTemp2; // 矩陣乘法等需要重複利用的function需要兩個暫存區 不然會覆寫到正在用的資料
	double* CompuTemp3; // 矩陣乘法等需要重複利用的function需要兩個暫存區 不然會覆寫到正在用的資料
	double* CompuTemp4; // 矩陣乘法等需要重複利用的function需要兩個暫存區 不然會覆寫到正在用的資料

	unsigned int LQDataLen; // 儲存LQ所使用資料長度

	fstream TimerRecord;

	static const int TimeDataTotal = 1000000;
	int TimeCount;
	int TimeCountPMS;
	double LogSysTime[TimeDataTotal];
	double LogPMSTime[TimeDataTotal];

	//20140916 dora test
	double test01[3][3000];
	double test02[3][3000];
	double test03[3][3000];
};
