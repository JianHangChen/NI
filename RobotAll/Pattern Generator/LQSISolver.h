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
     ���{���D�n�Φb�ѥX���˳��\�ҫ���inverted pendulum ���D
	 ��JZMP�H��COG���׭y�� �N�i�H�ѥX������VZMP�\�ʤ��y��
	 �U�禡�P�ܼƤ������иԨ��U��ŧi�P�w�q�B������

Note: None
***************************************************************************************************/

#include "Kine.h"
extern Kine gKineAll;//20140921DORA
// Optimzal control: weightings in the performance index
#define Q 100//100 // 10^6
#define R 1 // 10^-3
#define P 100//100 // 10^6

void DiffEqn(double* A, int* LenA, double* result, double* d_t); // �@���L�� �^�Ǥ֤F���G�k�G�Ȫ��}�C (�L���B��˱�̥��̥k��Ӥ���)
void DiffEqnShift(double* A, int* LenA, double* result, double* d_t); // �@���L�� �åB�^�Ǹ��}�C�@�˪��ת��}�C
void LowPassFilter(double* data, int dataLen, double* result, double gain); // �Ʀ�C�q�o�i��
void GetMTSAcceleration(double* sig, double* h, double* q, int InputTrajLen, double* Boundary, double* ddQ); // �p��MTS�Ѽ�
void MovingAvergeSmooth(double* Traj, int LengthOfTraj, int FrameSize, double* result); // moving average �����k ���Ʀ��u��

class LQSISolver
{
public:

	LQSISolver(void); // �غc�l
	~LQSISolver(void); // �Ѻc�l

	void Initval(double* InputCOGz); // ��J���߰��׭y�� �åB����l�B��
	void BackRiccati(double* ZMPx, double* ZMPy); // ��JZMP�y�� �åB�� riccati equation
	void BackRiccati2(double* ZMPx, double* ZMPy, int period, int index_now); // ��JZMP�y�� �åB�� riccati equation 
	void BackRiccati3(double* ZMPx, double* ZMPy, int period, int index_now); // 20141117�ץ��᪺���ɪ����A��JZMP�y�� �åB�� riccati equation 
	void C2DmHumanoid(int i); // �N�˳��\�ҫ��� state matrices �� zero order hold
	void DummyControl(void); // �N�w���Ȫ�����@�^�­� �w���|�N�X�@�t�C����COG�y�� �N�ӥΥH��@�����H����reference trajectory
	void DummyControl2(int index_now); // �N�w���Ȫ�����@�^�­� �w���|�N�X�@�t�C����COG�y�� �N�ӥΥH��@�����H����reference trajectory
	void DummyControl3(int index_now); // �N�w���Ȫ�����@�^�­� �w���|�N�X�@�t�C����COG�y�� �N�ӥΥH��@�����H����reference trajectory //20141117�ץ��᪺���ɪ���
	void ZMPFeedbackControl(double* ZMPx, double* ZMPy,double* COGx ,double* COGy , int index);// �p��C��iteration�����߭y��
	void tic(void); // �p�ɥΨ�� �p�ɶ}�l
	void toc(void); // �p�ɥΨ�� �p�ɵ����åB�bcommand window ��ܦ����p�ɮɶ�
	void toc2(void); // �p�ɥΨ�� �p�ɵ����åB�bcommand window ��̤ܳj�p�ɮɶ�(���ƭp�ɪ��ɭ�)
	void toc3(void); // �p�ɥΨ�� �p�ɵ����g�Jtxt �����C��Step ���ɶ�
	void ticPMS(void); // PMS�p�ɥΨ�� �p�ɶ}�l
	void toc3PMS(void); // PMS�p�ɥΨ�� �p�ɵ����g�Jtxt �����C��Step ���ɶ�


	YMatLite* Buffer0; // buffer for computation

	// LQSI �ϥιB��x�}�P���A�x�}�w�q
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
	YMatLite Sk_test[10000][600];//20141118 test for Sk �Ӥj�|�z
	double AngAccX[LQSIBufferSize];
	double AngAccY[LQSIBufferSize];
	double AngAccZ[LQSIBufferSize];
	double AngVelX[LQSIBufferSize];
	double AngVelY[LQSIBufferSize];
	double AngVelZ[LQSIBufferSize];
	//
	// timer �p�ɨϥ��ܼ�
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
	double* COGz; // ���߰��׭y��
	double delCOG[LQSIBufferSize]; // ���߰��׳t��
	double ddelCOG[LQSIBufferSize]; // ���߰��ץ[�t��

	double* CompuTemp; // �x�}���k���ݭn���ƧQ�Ϊ�function�ݭn��ӼȦs�� ���M�|�мg�쥿�b�Ϊ����
	double* CompuTemp2; // �x�}���k���ݭn���ƧQ�Ϊ�function�ݭn��ӼȦs�� ���M�|�мg�쥿�b�Ϊ����
	double* CompuTemp3; // �x�}���k���ݭn���ƧQ�Ϊ�function�ݭn��ӼȦs�� ���M�|�мg�쥿�b�Ϊ����
	double* CompuTemp4; // �x�}���k���ݭn���ƧQ�Ϊ�function�ݭn��ӼȦs�� ���M�|�мg�쥿�b�Ϊ����

	unsigned int LQDataLen; // �x�sLQ�Ҩϥθ�ƪ���

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
