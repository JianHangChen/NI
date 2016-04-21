/**************************************************************************************************
Copyright, 2010-2012, Robotics Lab., Dept. of M.E., National Taiwan University
File Name: LQSISolver.cpp

Author: Jiu-Lou Yan
Version: 1.0
Date: 2010/07/20

Functions:
     DiffEqn() DiffEqnShift() LowPassFilter() GetMTSAcceleration()
     MovingAvergeSmooth() LQSISolver()  ~LQSISolver() Initval()
	 BackRiccati() C2DmHumanoid() DummyControl()
	 ZMPFeedbackControl() tic()toc()toc2()

Classes: LQSISolver

Description:
     ���{���D�n�Φb�ѥX���˳��\�ҫ���inverted pendulum ���D
	 ��JZMP�H��COG���׭y�� �N�i�H�ѥX������VZMP�\�ʤ��y��
	 �U�禡�P�ܼƤ������иԨ��U��ŧi�P�w�q�B������

Note: None
***************************************************************************************************/

#include "stdafx.h"
#include <fstream>
#include "LQSISolver.h"
#include "IMU.h"

extern double MaxIterationTime; //Parameter for Computational Speed Calculation
extern int gFlagSimulation;//20140921 DORA
extern IMU IMU1;//20140921 DORA

LQSISolver::LQSISolver()
{
	/******************************************************************
	input: void
	output: void

	Note:
	// Class constructor  ��l�ƩҦ��ݭn�Ψ쪺�ܼƻP�x�}
	******************************************************************/
	for (int i = 0; i<LQSIBufferSize ; i++)
	{
		AdStk[i].InitPara(3,3);
		AIStk[i].InitPara(3,3);
		BdStk[i].InitPara(3,1);
		BinvMBT_Stk[i].InitPara(3,3);
		WkStk[i].InitPara(3,3);
		WkTStk[i].InitPara(3,3);
		invWkStk[i].InitPara(3,3);
		NkStk[i].InitPara(3,3);
		EkStk[i].InitPara(3,3);
		NkTStk[i].InitPara(3,3);
		SsStk[i].InitPara(3,3);
		vsStkX[i].InitPara(3,1);
		vsStkY[i].InitPara(3,1);
		invDkStk[i].InitPara(3,3);
		XState[i].InitPara(3,1);
		YState[i].InitPara(3,1);
		Z_Schur[i].InitPara(6,6);
		Z_Schur_U[i].InitPara(6,6);
		//20140825 for state now
		AdStk2[i].InitPara(3,3);
		BdStk2[i].InitPara(3,1);
		XState2[i].InitPara(3,1);
		YState2[i].InitPara(3,1);
		
	}
	for(int k=0;k<LQSIBufferSize/*2001*/;k++)
		DkStk[k].InitPara(3,3);

	for (int i=0;i<10000;i++)
		for (int j=0;j<600;j++)
			Sk_test[i][j].InitPara(3,3);

	SsStk[int(LQSIBufferSize)].InitPara(3,3); 
	vsStkX[int(LQSIBufferSize)].InitPara(3,1);
	vsStkY[int(LQSIBufferSize)].InitPara(3,1);
	XState[int(LQSIBufferSize)].InitPara(3,1);
	YState[int(LQSIBufferSize)].InitPara(3,1);

	CompuTemp = new double[400]; // �x�}���k���ݭn���ƧQ�Ϊ�function�ݭn��ӼȦs�� ���M�|�мg�쥿�b�Ϊ����
	CompuTemp2 = new double[400];
	CompuTemp3 = new double[400];
	CompuTemp4 = new double[400];

	Qx = new YMatLite[1];
	Qx->InitPara(3,3);

	Qx->data = new double[Qx->MSize];
	Qx->data[0] = 10;//400000;//10
	Qx->data[1] = 0;
	Qx->data[2] = 0;
	Qx->data[3] = 0;
	Qx->data[4] = 0;//400000;
	Qx->data[5] = 0;
	Qx->data[6] = 0;
	Qx->data[7] = 0;
	Qx->data[8] = 0;//40000;

	//Qx = 10, 0, 0,
	//	    0, 0, 0,
	//	    0, 0, 0;


	Cd = new YMatLite[1];
	Cd->InitPara(1,3);

	Cd->data[0] = 0;
	Cd->data[1] = 0;
	Cd->data[2] = 1;

	CdT = new YMatLite[1];
	CdT->InitPara(3,1);
	CdT->data[0] = 0;
	CdT->data[1] = 0;
	CdT->data[2] = 1;

	Dc = 0;
	Dd = 0;

	LQDataLen = LQSIBufferSize; // Back Riccati Length

	// initial ���|�ܪ���
	for (int i = 0; i< LQDataLen ; i++)
	{
		AdStk[i].data[6] = 0;
		AdStk[i].data[7] = 0;
		AdStk[i].data[8] = 1;
		BdStk[i].data[2] = dt;
	}

	if (QueryPerformanceFrequency(&gnFreq))
		gFreqT = (float)(gnFreq.QuadPart);
	
	if (QueryPerformanceFrequency(&gnFreq_PMS))
		gFreqT_PMS = (float)(gnFreq_PMS.QuadPart);

	TimeCount=0;
	TimeCountPMS = 0;
}

LQSISolver::~LQSISolver(void)
{
	/******************************************************************
	input: void
	output: void

	Note:
	// Class destructor 
	******************************************************************/

	// �M���ʺA�O����
	delete[] Cd;
	delete[] CdT;
	delete[] Qx;

	delete[] CompuTemp;
	delete[] CompuTemp2;
	delete[] CompuTemp3;
	delete[] CompuTemp4;

}

void LQSISolver::Initval(double* InputCOGz)
{
	/******************************************************************
	input: InputCOGz
	output: void

	Note:
	// �q�~����JCOG���׭y�� �o�ˤ~���X�Ҧ��˳��\�ҫ��b�C�@���������A�x�}
	******************************************************************/

	int LenTemp = LQDataLen+8;
	double SampTime = dt;
								
	COGz = InputCOGz;

	// �ﭫ�߰��׭y��L��
	// �ɶ����ܵu ��10^-5
	DiffEqn(COGz,&LenTemp,delCOG,&SampTime);
	LenTemp = LQDataLen+4;
	DiffEqn(delCOG,&LenTemp,ddelCOG,&SampTime);
	// �ɶ����ܵu ��10^-5

	int cnt = 0;

	for (int i = 0 ; i < LQDataLen ; i++)
	{
		if (ddelCOG[i] <= (-GravityConst))
		{
			cnt += 1;
			ddelCOG[i] = (-GravityConst)+1;
			//printf("index = %d \n",i);
		}
	}

	if (cnt > 0)
	{
		printf("�Фp�� �� %d ��COG�[�t�׶W�L 9810mm/s^2 �Ъ`�N!! �ȳQ�j��令-9809�F",cnt);
	}
	//*****test print
	fstream dcogz;
	dcogz.open("dInpCOGz.txt",ios::out);
	dcogz.precision(10);
	for (int i=0 ; i<= LQDataLen ; i++)
	{
		dcogz << delCOG[i]<< endl;
	}
	dcogz.close();
	fstream ddcogz;
	ddcogz.open("ddInpCOGz.txt",ios::out);
	ddcogz.precision(10);
	for (int i=0 ; i<= LQDataLen ; i++)
	{
		ddcogz << ddelCOG[i]<< endl;
	}
	ddcogz.close();

}


void LQSISolver::BackRiccati(double* ZMPx, double* ZMPy)
{
	/******************************************************************
	input: ZMPx ZMPy, ZMP�b������V���y�� �o�OLQSI controller�� reference input
	output: void

	Note:
	// ��JMZPx ZMPy��LQSI controller��@reference input
	******************************************************************/

	YMatLite EYE3(3,3);
	EYE3.data[0] = 1;
	EYE3.data[1] = 0;
	EYE3.data[2] = 0;
	EYE3.data[3] = 0;
	EYE3.data[4] = 1;
	EYE3.data[5] = 0;
	EYE3.data[6] = 0;
	EYE3.data[7] = 0;
	EYE3.data[8] = 1;

	// ���Xzero order hold �� ���A�x�}(state-space matrices)
	// 0.24ms --> 0.13ms
	for (int i = 0 ; i < LQDataLen ; i++)
	{
		C2DmHumanoid(i);
	}
	// 0.24ms


	//tic();

	// LQSI�����ܼ� AI = Ad - (identity matrix 3-by-3)
	// 0.2ms --> 0.018ms
	for (int i = 0; i < LQDataLen ; i++)
	{
		AIStk[i].data[0] = AdStk[i].data[0]-1;
		AIStk[i].data[1] = AdStk[i].data[1];
		AIStk[i].data[2] = AdStk[i].data[2];
		AIStk[i].data[3] = AdStk[i].data[3];
		AIStk[i].data[4] = AdStk[i].data[4]-1;
		AIStk[i].data[5] = AdStk[i].data[5];
		AIStk[i].data[6] = AdStk[i].data[6];
		AIStk[i].data[7] = AdStk[i].data[7];
		AIStk[i].data[8] = AdStk[i].data[8]-1;
	}
	// 0.2ms

	// �p��LQSI Control Law~ 
	// 0.08ms
	for (int i = 0; i < LQDataLen ; i++)
	{
		// ��l equation
		//gTempM = BdT_stk[i]*Qx*BdStk[i]; 
		//invMkStk[i] = 1.0/(gTempM.data()[0]+R); // invMkStk(1,i) = 1/(Bd'*Qx*Bd+R);

		MatMulAtB(BdStk[i].data,3,1,Qx->data,3,3,CompuTemp);
		MatMulAB(CompuTemp,1,3,BdStk[i].data,3,1,CompuTemp2);
		invMkStk[i] = 1.0/(CompuTemp2[0]+R); // invMkStk(1,i) = 1/(Bd'*Qx*Bd+R);

		//temp11 = (BdStk[i]>Qx[0])*BdStk[i];
		//invMkStk[i] = 1.0/(temp11.data[0]+R); // invMkStk(1,i) = 1/(Bd'*Qx*Bd+R);

	}
	//0.08ms


	// �p��LQSI Control Law~ 
	for (int i = 0; i < LQDataLen ; i++)
	{
		// 0.1ms
		MatScalarMul(BdStk[i].data,3,invMkStk+i,CompuTemp);
		MatMulABt(CompuTemp,3,1,BdStk[i].data,3,1,BinvMBT_Stk[i].data);
		// 0.1ms
	}

	// �p��LQSI Control Law~ 
	for (int i = 0; i < LQDataLen ; i++)
	{
		// 0.31ms
		MatMulAB(BinvMBT_Stk[i].data,3,3,Qx->data,3,3,CompuTemp);
		MatMulAB(CompuTemp,3,3,AIStk[i].data,3,3,CompuTemp2);
		MatMiuAB(AdStk[i].data,CompuTemp2,WkStk[i].data,9);
	}


	// �p��LQSI Control Law~ 
	// 0.24ms
	for (int i = 0; i < LQDataLen ; i++)
	{
	   // NkStk[i] = WkStk[i] - EYE3;
		NkStk[i].data[0] = WkStk[i].data[0]-1;
		NkStk[i].data[1] = WkStk[i].data[1];
		NkStk[i].data[2] = WkStk[i].data[2];
		NkStk[i].data[3] = WkStk[i].data[3];
		NkStk[i].data[4] = WkStk[i].data[4]-1;
		NkStk[i].data[5] = WkStk[i].data[5];
		NkStk[i].data[6] = WkStk[i].data[6];
		NkStk[i].data[7] = WkStk[i].data[7];
		NkStk[i].data[8] = WkStk[i].data[8]-1;

	}
	// 0.24ms 0.188

	double tempScale = 0;

	// �p��LQSI Control Law~ 
	/////// �o�̥i�H�� �w����n���̲׭Ȩ��N
	tempScale = P;

	vsStkX[(int)(LQDataLen)].data[0] = CdT->data[0]*(P*ZMPx[(int)(LQDataLen-1)]); // ���w�̫᪺��(backward riccati���}�l)
	vsStkX[(int)(LQDataLen)].data[1] = CdT->data[1]*(P*ZMPx[(int)(LQDataLen-1)]); // ���w�̫᪺��(backward riccati���}�l)
	vsStkX[(int)(LQDataLen)].data[2] = CdT->data[2]*(P*ZMPx[(int)(LQDataLen-1)]); // ���w�̫᪺��(backward riccati���}�l) vn final value

	vsStkY[(int)(LQDataLen)].data[0] = CdT->data[0]*(P*ZMPy[(int)(LQDataLen-1)]); // ���w�̫᪺��(backward riccati���}�l)
	vsStkY[(int)(LQDataLen)].data[1] = CdT->data[1]*(P*ZMPy[(int)(LQDataLen-1)]); // ���w�̫᪺��(backward riccati���}�l)
	vsStkY[(int)(LQDataLen)].data[2] = CdT->data[2]*(P*ZMPy[(int)(LQDataLen-1)]); // ���w�̫᪺��(backward riccati���}�l)

	MatScalarMul(Cd->data,3,&tempScale,CompuTemp);
	MatMulAB(CdT->data,3,1,CompuTemp,1,3,SsStk[(int)(LQDataLen)].data);//Sn final value

	/////// �o�̥i�H�� �w����n���̲׭Ȩ��N

	// �p��LQSI Control Law~ 
	// 88ms original
	for (int i = (int)(LQDataLen-1) ; i >=0 ; i--)
	{
		// MATLAB code
		////invDkStk(:,:,i) = inv(eye(3)+BinvMBStk(:,:,i)*Ss(:,:,i+1));
		////Ss(:,:,i) = Nk'*Qx*Nk+(invMkStk(1,i)*invMkStk(1,i)*R)*(AI'*Qx*Bd*Bd'*Qx*AI)+CQC+Ad'*Ss(:,:,i+1)*invDkStk(:,:,i)*Wk;
		////vs(:,i) = Cd'*Q*ZMPr(i)+Ad'*(eye(3)-Ss(:,:,i+1)*invDkStk(:,:,i)*BinvMBStk(:,:,i))*vs(:,i+1);
		// MATLAB code

		//invDkStk[i] = EYE3+BinvMBT_Stk[i]*SsStk[i+1]; // �U����檺�즡		Dk
		MatMulAB(BinvMBT_Stk[i].data,3,3,SsStk[i+1].data,3,3,CompuTemp);
		MatAddAB(EYE3.data,CompuTemp,invDkStk[i].data,9);

		InvSqMat(invDkStk[i].data,3); // Take Inverse: call by reference ���ݭn�t�~�[�t inv Dk

		// Ss �즡
		//SsStk[i] = NkTStk[i]*Qx*NkStk[i]+((invMkStk[i]*invMkStk[i]*R)*(AIT_stk[i]*Qx*BdStk[i]))*(BdT_stk[i]*Qx*AIStk[i]) + CdT*Q*Cd + AdT_stk[i]*SsStk[i+1]*invDkStk[i]*WkStk[i];  Sk

		// ��� Ss, �٤U 50ms
		MatMulAtB(NkStk[i].data,3,3,Qx->data,3,3,CompuTemp);
		MatMulAB(CompuTemp,3,3,NkStk[i].data,3,3,CompuTemp2); // �⧹ N'QxN �s�b CompuTemp2

		MatMulAB(Qx->data,3,3,BdStk[i].data,3,1,CompuTemp);
		MatMulAtB(AIStk[i].data,3,3,CompuTemp,3,1,CompuTemp3);

		tempScale = invMkStk[i]*invMkStk[i]*R;
		MatScalarMul(CompuTemp3,3,&tempScale,CompuTemp4);
		MatMulABt(CompuTemp4,3,1,CompuTemp3,3,1,CompuTemp); // �⧹���q�A�s�b CompuTemp

		//MatScalarMul(CompuTemp3,3,&tempScale,CompuTemp3);
		//MatMulAtB(BdStk[i].data(),3,1,Qx.data(),3,3,CompuTemp);
		//MatMulAB(CompuTemp,1,3,AIStk[i].data(),3,3,CompuTemp4);
		//MatMulAB(CompuTemp3,3,1,CompuTemp4,1,3,CompuTemp); // �⧹���q�A�s�b CompuTemp

		MatAddAB(CompuTemp,CompuTemp2,CompuTemp,9); // N'QN + ���q�A�s�b CompuTemp

		tempScale = Q;
		MatScalarMul(Cd->data,3,&tempScale,CompuTemp2);
		MatMulAtB(Cd->data,1,3,CompuTemp2,1,3,CompuTemp3); // �⧹ C'QC

		MatAddAB(CompuTemp,CompuTemp3,CompuTemp,9); // N'QN + ���q + C'QC�A�s�b CompuTemp

		MatMulAtB(AdStk[i].data,3,3,SsStk[i+1].data,3,3,CompuTemp2);
		MatMulAB(CompuTemp2,3,3,invDkStk[i].data,3,3,CompuTemp3);
		MatMulAB(CompuTemp3,3,3,WkStk[i].data,3,3,CompuTemp2);

		MatAddAB(CompuTemp,CompuTemp2,SsStk[i].data,9); // �⧹ Ss�A�s�b SsStk[i]
		// ��� Ss

		//// ��� vs �٤U30ms
		////vs_Stk[i] = CdT*Q*Input_ZMP[i]+AdT_stk[i]*(EYE3-SsStk[i+1]*invDkStk[i]*BinvMBT_Stk[i])*vs_Stk[i+1];

		//tempScale = Q*Input_ZMP[i];
		//MatScalarMul(Cd->data,3,&tempScale,CompuTemp); // C'Qr

		//MatMulAB(BinvMBT_Stk[i].data,3,3,vs_Stk[i+1].data,3,1,CompuTemp2);
		//MatMulAB(invDkStk[i].data,3,3,CompuTemp2,3,1,CompuTemp3);
		//MatMulAB(SsStk[i+1].data,3,3,CompuTemp3,3,1,CompuTemp2);

		//MatMiuAB(vs_Stk[i+1].data,CompuTemp2,CompuTemp3,3);
		//MatMulAtB(AdStk[i].data,3,3,CompuTemp3,3,1,CompuTemp2); // �⧹��b

		//MatAddAB(CompuTemp,CompuTemp2,vs_Stk[i].data,3);
		//// ��� vs

	}

	// for ZMPx
	for (int i = (int)(LQDataLen-1) ; i >=0 ; i--)
	{
		// ��� vs �٤U30ms
		//vs_Stk[i] = CdT*Q*Input_ZMP[i]+AdT_stk[i]*(EYE3-SsStk[i+1]*invDkStk[i]*BinvMBT_Stk[i])*vs_Stk[i+1]; vkx

		tempScale = Q*ZMPx[i];
		MatScalarMul(Cd->data,3,&tempScale,CompuTemp); // C'Qr

		MatMulAB(BinvMBT_Stk[i].data,3,3,vsStkX[i+1].data,3,1,CompuTemp2);
		MatMulAB(invDkStk[i].data,3,3,CompuTemp2,3,1,CompuTemp3);
		MatMulAB(SsStk[i+1].data,3,3,CompuTemp3,3,1,CompuTemp2);

		MatMiuAB(vsStkX[i+1].data,CompuTemp2,CompuTemp3,3);
		MatMulAtB(AdStk[i].data,3,3,CompuTemp3,3,1,CompuTemp2); // �⧹��b

		MatAddAB(CompuTemp,CompuTemp2,vsStkX[i].data,3);
		// ��� vs
	}

	// for ZMPy
	for (int i = (int)(LQDataLen-1) ; i >=0 ; i--)
	{
		// ��� vs �٤U30ms
		//vs_Stk[i] = CdT*Q*Input_ZMP[i]+AdT_stk[i]*(EYE3-SsStk[i+1]*invDkStk[i]*BinvMBT_Stk[i])*vs_Stk[i+1]; vky

		tempScale = Q*ZMPy[i];
		MatScalarMul(Cd->data,3,&tempScale,CompuTemp); // C'Qr

		MatMulAB(BinvMBT_Stk[i].data,3,3,vsStkY[i+1].data,3,1,CompuTemp2);
		MatMulAB(invDkStk[i].data,3,3,CompuTemp2,3,1,CompuTemp3);
		MatMulAB(SsStk[i+1].data,3,3,CompuTemp3,3,1,CompuTemp2);

		MatMiuAB(vsStkY[i+1].data,CompuTemp2,CompuTemp3,3);
		MatMulAtB(AdStk[i].data,3,3,CompuTemp3,3,1,CompuTemp2); // �⧹��b

		MatAddAB(CompuTemp,CompuTemp2,vsStkY[i].data,3);
		// ��� vs

	}

}

void LQSISolver::BackRiccati2(double* ZMPx, double* ZMPy ,int period, int index_now)
{
	/******************************************************************
	input: ZMPx ZMPy, ZMP�b������V���y�� �o�OLQSI controller�� reference input
	output: void

	Note:
	// ��JMZPx ZMPy��LQSI controller��@reference input
	******************************************************************/

	YMatLite EYE3(3,3);
	EYE3.data[0] = 1;
	EYE3.data[1] = 0;
	EYE3.data[2] = 0;
	EYE3.data[3] = 0;
	EYE3.data[4] = 1;
	EYE3.data[5] = 0;
	EYE3.data[6] = 0;
	EYE3.data[7] = 0;
	EYE3.data[8] = 1;

	// ���Xzero order hold �� ���A�x�}(state-space matrices)
	// 0.24ms --> 0.13ms original
	//tic();//release avg 0.094ms 600��/loop --> 0.002ms 1��/loop
	if (index_now==0)
	{
		for (int i = index_now/*0*/ ; i < index_now + period ; i++)
		{
			C2DmHumanoid(i);
		}
	} 
	else
	{
		C2DmHumanoid(index_now + period-1);
	}
	//toc();
	// 0.24ms


	//tic();

	// LQSI�����ܼ� AI = Ad - (identity matrix 3-by-3)
	// 0.2ms --> 0.018ms original
	//tic();//release avg 0.02ms 600��/loop --> 0.00016ms 1��/loop
	if (index_now==0)
	{
		for (int i = index_now/*0*/; i < index_now + period ; i++)
		{
			AIStk[i].data[0] = AdStk[i].data[0]-1;
			AIStk[i].data[1] = AdStk[i].data[1];
			AIStk[i].data[2] = AdStk[i].data[2];
			AIStk[i].data[3] = AdStk[i].data[3];
			AIStk[i].data[4] = AdStk[i].data[4]-1;
			AIStk[i].data[5] = AdStk[i].data[5];
			AIStk[i].data[6] = AdStk[i].data[6];
			AIStk[i].data[7] = AdStk[i].data[7];
			AIStk[i].data[8] = AdStk[i].data[8]-1;
		}
	} 
	else
	{
		AIStk[index_now + period-1].data[0] = AdStk[index_now + period-1].data[0]-1;
		AIStk[index_now + period-1].data[1] = AdStk[index_now + period-1].data[1];
		AIStk[index_now + period-1].data[2] = AdStk[index_now + period-1].data[2];
		AIStk[index_now + period-1].data[3] = AdStk[index_now + period-1].data[3];
		AIStk[index_now + period-1].data[4] = AdStk[index_now + period-1].data[4]-1;
		AIStk[index_now + period-1].data[5] = AdStk[index_now + period-1].data[5];
		AIStk[index_now + period-1].data[6] = AdStk[index_now + period-1].data[6];
		AIStk[index_now + period-1].data[7] = AdStk[index_now + period-1].data[7];
		AIStk[index_now + period-1].data[8] = AdStk[index_now + period-1].data[8]-1;
	}
	//toc();
	// 0.2ms

	// �p��LQSI Control Law~ 
	// 0.08ms original
	//tic();//release avg 0.048ms 600��/loop --> 0.00067ms 1��/loop
	if (index_now==0)
	{
		for (int i = index_now/*0*/; i < index_now + period ; i++)
		{
			// ��l equation
			//gTempM = BdT_stk[i]*Qx*BdStk[i]; 
			//invMkStk[i] = 1.0/(gTempM.data()[0]+R); // invMkStk(1,i) = 1/(Bd'*Qx*Bd+R);

			MatMulAtB(BdStk[i].data,3,1,Qx->data,3,3,CompuTemp);
			MatMulAB(CompuTemp,1,3,BdStk[i].data,3,1,CompuTemp2);
			invMkStk[i] = 1.0/(CompuTemp2[0]+R); // invMkStk(1,i) = 1/(Bd'*Qx*Bd+R);

			//temp11 = (BdStk[i]>Qx[0])*BdStk[i];
			//invMkStk[i] = 1.0/(temp11.data[0]+R); // invMkStk(1,i) = 1/(Bd'*Qx*Bd+R);
		}
	} 
	else
	{
		MatMulAtB(BdStk[index_now + period-1].data,3,1,Qx->data,3,3,CompuTemp);
		MatMulAB(CompuTemp,1,3,BdStk[index_now + period-1].data,3,1,CompuTemp2);
		invMkStk[index_now + period-1] = 1.0/(CompuTemp2[0]+R); // invMkStk(1,i) = 1/(Bd'*Qx*Bd+R);
	}
	//toc();
	//0.08ms


	// �p��LQSI Control Law~ 
	//�H�U��̲׭ȫe���i���@�p�ɤF*******20140918 Start
	//tic();//release avg  0.00074ms 1��/loop
	if (index_now==0)
	{
		for (int i = index_now/*0*/; i < index_now + period ; i++)
		{
			// 0.1ms
			MatScalarMul(BdStk[i].data,3,invMkStk+i,CompuTemp);
			MatMulABt(CompuTemp,3,1,BdStk[i].data,3,1,BinvMBT_Stk[i].data);
			// 0.1ms
		}
	} 
	else
	{
		MatScalarMul(BdStk[index_now + period-1].data,3,invMkStk+(index_now + period-1),CompuTemp);
		MatMulABt(CompuTemp,3,1,BdStk[index_now + period-1].data,3,1,BinvMBT_Stk[index_now + period-1].data);
	}
	

	// �p��LQSI Control Law~ 
	if (index_now==0)
	{
		for (int i = index_now/*0*/; i < index_now + period ; i++)
		{
			// 0.31ms
			MatMulAB(BinvMBT_Stk[i].data,3,3,Qx->data,3,3,CompuTemp);
			MatMulAB(CompuTemp,3,3,AIStk[i].data,3,3,CompuTemp2);
			MatMiuAB(AdStk[i].data,CompuTemp2,WkStk[i].data,9);
		}
	} 
	else
	{
		MatMulAB(BinvMBT_Stk[index_now + period-1].data,3,3,Qx->data,3,3,CompuTemp);
		MatMulAB(CompuTemp,3,3,AIStk[index_now + period-1].data,3,3,CompuTemp2);
		MatMiuAB(AdStk[index_now + period-1].data,CompuTemp2,WkStk[index_now + period-1].data,9);
	}

	// �p��LQSI Control Law~ 
	// 0.24ms
	if (index_now==0)
	{
		for (int i = index_now/*0*/; i < index_now + period ; i++)
		{
		   // NkStk[i] = WkStk[i] - EYE3;
			NkStk[i].data[0] = WkStk[i].data[0]-1;
			NkStk[i].data[1] = WkStk[i].data[1];
			NkStk[i].data[2] = WkStk[i].data[2];
			NkStk[i].data[3] = WkStk[i].data[3];
			NkStk[i].data[4] = WkStk[i].data[4]-1;
			NkStk[i].data[5] = WkStk[i].data[5];
			NkStk[i].data[6] = WkStk[i].data[6];
			NkStk[i].data[7] = WkStk[i].data[7];
			NkStk[i].data[8] = WkStk[i].data[8]-1;
		}
	} 
	else
	{
		NkStk[index_now + period-1].data[0] = WkStk[index_now + period-1].data[0]-1;
		NkStk[index_now + period-1].data[1] = WkStk[index_now + period-1].data[1];
		NkStk[index_now + period-1].data[2] = WkStk[index_now + period-1].data[2];
		NkStk[index_now + period-1].data[3] = WkStk[index_now + period-1].data[3];
		NkStk[index_now + period-1].data[4] = WkStk[index_now + period-1].data[4]-1;
		NkStk[index_now + period-1].data[5] = WkStk[index_now + period-1].data[5];
		NkStk[index_now + period-1].data[6] = WkStk[index_now + period-1].data[6];
		NkStk[index_now + period-1].data[7] = WkStk[index_now + period-1].data[7];
		NkStk[index_now + period-1].data[8] = WkStk[index_now + period-1].data[8]-1;
	}
	//toc();
	// 0.24ms 0.188
	//*******************20140918 end
	double tempScale = 0;

	// �p��LQSI Control Law~ 
	/////// �o�̥i�H�� �w����n���̲׭Ȩ��N ++period
	tempScale = P;

	vsStkX[(int)(index_now + period+1)].data[0] = CdT->data[0]*(P*ZMPx[(int)(index_now + period/*-1*/)]); // ���w�̫᪺��(backward riccati���}�l)//20141020 vsStkX[(int)(index_now + period)].data[0] = CdT->data[0]*(P*ZMPx[(int)(index_now + period-1)]); // ���w�̫᪺��(backward riccati���}�l)
	vsStkX[(int)(index_now + period+1)].data[1] = CdT->data[1]*(P*ZMPx[(int)(index_now + period/*-1*/)]); // ���w�̫᪺��(backward riccati���}�l)
	vsStkX[(int)(index_now + period+1)].data[2] = CdT->data[2]*(P*ZMPx[(int)(index_now + period/*-1*/)]); // ���w�̫᪺��(backward riccati���}�l) vn regulation final value

	vsStkY[(int)(index_now + period+1)].data[0] = CdT->data[0]*(P*ZMPy[(int)(index_now + period/*-1*/)]); // ���w�̫᪺��(backward riccati���}�l)
	vsStkY[(int)(index_now + period+1)].data[1] = CdT->data[1]*(P*ZMPy[(int)(index_now + period/*-1*/)]); // ���w�̫᪺��(backward riccati���}�l)
	vsStkY[(int)(index_now + period+1)].data[2] = CdT->data[2]*(P*ZMPy[(int)(index_now + period/*-1*/)]); // ���w�̫᪺��(backward riccati���}�l)

	MatScalarMul(Cd->data,3,&tempScale,CompuTemp);
	MatMulAB(CdT->data,3,1,CompuTemp,1,3,SsStk[(int)(index_now + period+1)].data);//Sn final value  //20141020 MatMulAB(CdT->data,3,1,CompuTemp,1,3,SsStk[(int)(index_now + period)].data);//Sn final value

	/////// �o�̥i�H�� �w����n���̲׭Ȩ��N

	// �p��LQSI Control Law~ 
	// 88ms original
	for (int i = (int)(index_now + period) ; i >= index_now + 0 ; i--)//back period  20140823 ����i��n��state�@�_��s//20141020 for (int i = (int)(index_now + period-1) ; i >=index_now + 0 ; i--)
	{
		// MATLAB code
		////invDkStk(:,:,i) = inv(eye(3)+BinvMBStk(:,:,i)*Ss(:,:,i+1));
		////Ss(:,:,i) = Nk'*Qx*Nk+(invMkStk(1,i)*invMkStk(1,i)*R)*(AI'*Qx*Bd*Bd'*Qx*AI)+CQC+Ad'*Ss(:,:,i+1)*invDkStk(:,:,i)*Wk;
		////vs(:,i) = Cd'*Q*ZMPr(i)+Ad'*(eye(3)-Ss(:,:,i+1)*invDkStk(:,:,i)*BinvMBStk(:,:,i))*vs(:,i+1);
		// MATLAB code

		//invDkStk[i] = EYE3+BinvMBT_Stk[i]*SsStk[i+1]; // �U����檺�즡		Dk
		MatMulAB(BinvMBT_Stk[i].data,3,3,SsStk[i+1].data,3,3,CompuTemp);
		MatAddAB(EYE3.data,CompuTemp,invDkStk[i].data,9);

		InvSqMat(invDkStk[i].data,3); // Take Inverse: call by reference ���ݭn�t�~�[�t inv Dk

		// Ss �즡
		//SsStk[i] = NkTStk[i]*Qx*NkStk[i]+((invMkStk[i]*invMkStk[i]*R)*(AIT_stk[i]*Qx*BdStk[i]))*(BdT_stk[i]*Qx*AIStk[i]) + CdT*Q*Cd + AdT_stk[i]*SsStk[i+1]*invDkStk[i]*WkStk[i];  Sk

		// ��� Ss, �٤U 50ms
		MatMulAtB(NkStk[i].data,3,3,Qx->data,3,3,CompuTemp);
		MatMulAB(CompuTemp,3,3,NkStk[i].data,3,3,CompuTemp2); // �⧹ N'QxN �s�b CompuTemp2

		MatMulAB(Qx->data,3,3,BdStk[i].data,3,1,CompuTemp);
		MatMulAtB(AIStk[i].data,3,3,CompuTemp,3,1,CompuTemp3);

		tempScale = invMkStk[i]*invMkStk[i]*R;
		MatScalarMul(CompuTemp3,3,&tempScale,CompuTemp4);
		MatMulABt(CompuTemp4,3,1,CompuTemp3,3,1,CompuTemp); // �⧹���q�A�s�b CompuTemp

		//MatScalarMul(CompuTemp3,3,&tempScale,CompuTemp3);
		//MatMulAtB(BdStk[i].data(),3,1,Qx.data(),3,3,CompuTemp);
		//MatMulAB(CompuTemp,1,3,AIStk[i].data(),3,3,CompuTemp4);
		//MatMulAB(CompuTemp3,3,1,CompuTemp4,1,3,CompuTemp); // �⧹���q�A�s�b CompuTemp

		MatAddAB(CompuTemp,CompuTemp2,CompuTemp,9); // N'QN + ���q�A�s�b CompuTemp

		tempScale = Q;
		MatScalarMul(Cd->data,3,&tempScale,CompuTemp2);
		MatMulAtB(Cd->data,1,3,CompuTemp2,1,3,CompuTemp3); // �⧹ C'QC

		MatAddAB(CompuTemp,CompuTemp3,CompuTemp,9); // N'QN + ���q + C'QC�A�s�b CompuTemp

		MatMulAtB(AdStk[i].data,3,3,SsStk[i+1].data,3,3,CompuTemp2);
		MatMulAB(CompuTemp2,3,3,invDkStk[i].data,3,3,CompuTemp3);
		MatMulAB(CompuTemp3,3,3,WkStk[i].data,3,3,CompuTemp2);

		MatAddAB(CompuTemp,CompuTemp2,SsStk[i].data,9); // �⧹ Ss�A�s�b SsStk[i]
		// ��� Ss

		//// ��� vs �٤U30ms
		////vs_Stk[i] = CdT*Q*Input_ZMP[i]+AdT_stk[i]*(EYE3-SsStk[i+1]*invDkStk[i]*BinvMBT_Stk[i])*vs_Stk[i+1];

		//tempScale = Q*Input_ZMP[i];
		//MatScalarMul(Cd->data,3,&tempScale,CompuTemp); // C'Qr

		//MatMulAB(BinvMBT_Stk[i].data,3,3,vs_Stk[i+1].data,3,1,CompuTemp2);
		//MatMulAB(invDkStk[i].data,3,3,CompuTemp2,3,1,CompuTemp3);
		//MatMulAB(SsStk[i+1].data,3,3,CompuTemp3,3,1,CompuTemp2);

		//MatMiuAB(vs_Stk[i+1].data,CompuTemp2,CompuTemp3,3);
		//MatMulAtB(AdStk[i].data,3,3,CompuTemp3,3,1,CompuTemp2); // �⧹��b

		//MatAddAB(CompuTemp,CompuTemp2,vs_Stk[i].data,3);
		//// ��� vs

	}

	// for ZMPx
	for (int i = (int)(index_now + period) ; i >=index_now + 0 ; i--)//back period  20140823 ����i��n��state�@�_��s
	{
		// ��� vs �٤U30ms
		//vs_Stk[i] = CdT*Q*Input_ZMP[i]+AdT_stk[i]*(EYE3-SsStk[i+1]*invDkStk[i]*BinvMBT_Stk[i])*vs_Stk[i+1]; vkx

		tempScale = Q*ZMPx[i];//regulation goal
		MatScalarMul(Cd->data,3,&tempScale,CompuTemp); // C'Qr

		MatMulAB(BinvMBT_Stk[i].data,3,3,vsStkX[i+1].data,3,1,CompuTemp2);
		MatMulAB(invDkStk[i].data,3,3,CompuTemp2,3,1,CompuTemp3);
		MatMulAB(SsStk[i+1].data,3,3,CompuTemp3,3,1,CompuTemp2);

		MatMiuAB(vsStkX[i+1].data,CompuTemp2,CompuTemp3,3);
		MatMulAtB(AdStk[i].data,3,3,CompuTemp3,3,1,CompuTemp2); // �⧹��b

		MatAddAB(CompuTemp,CompuTemp2,vsStkX[i].data,3);
		// ��� vs
	}

	// for ZMPy
	for (int i = (int)(index_now + period) ; i >=index_now + 0 ; i--)
	{
		// ��� vs �٤U30ms
		//vs_Stk[i] = CdT*Q*Input_ZMP[i]+AdT_stk[i]*(EYE3-SsStk[i+1]*invDkStk[i]*BinvMBT_Stk[i])*vs_Stk[i+1]; vky

		tempScale = Q*ZMPy[i];
		MatScalarMul(Cd->data,3,&tempScale,CompuTemp); // C'Qr

		MatMulAB(BinvMBT_Stk[i].data,3,3,vsStkY[i+1].data,3,1,CompuTemp2);
		MatMulAB(invDkStk[i].data,3,3,CompuTemp2,3,1,CompuTemp3);
		MatMulAB(SsStk[i+1].data,3,3,CompuTemp3,3,1,CompuTemp2);

		MatMiuAB(vsStkY[i+1].data,CompuTemp2,CompuTemp3,3);
		MatMulAtB(AdStk[i].data,3,3,CompuTemp3,3,1,CompuTemp2); // �⧹��b

		MatAddAB(CompuTemp,CompuTemp2,vsStkY[i].data,3);
		// ��� vs

	}

}
void LQSISolver::BackRiccati3(double* ZMPx, double* ZMPy ,int period, int index_now)
{
	/******************************************************************
	input: ZMPx ZMPy, ZMP�b������V���y�� �o�OLQSI controller�� reference input
	output: void

	Note:
	// ��JMZPx ZMPy��LQSI controller��@reference input
	******************************************************************/

	YMatLite EYE3(3,3);
	EYE3.data[0] = 1;
	EYE3.data[1] = 0;
	EYE3.data[2] = 0;
	EYE3.data[3] = 0;
	EYE3.data[4] = 1;
	EYE3.data[5] = 0;
	EYE3.data[6] = 0;
	EYE3.data[7] = 0;
	EYE3.data[8] = 1;

	/*double AA[9]={1,2,3,4,5,6,7,8,9};
	double BB[9]={1,3,5,7,9,11,13,15,17};
	double CC[9]={0,0,0,0,0,0,0,0,0};
	MatMulAtB(AA,3,3,BB,3,3,CC);*/
	// ���Xzero order hold �� ���A�x�}(state-space matrices)
	// 0.24ms --> 0.13ms original
	//tic();//release avg 0.094ms 600��/loop --> 0.002ms 1��/loop
	if (index_now==0)
	{
		for (int i = index_now/*0*/ ; i <= index_now + period ; i++)
		{
			C2DmHumanoid(i);
		}
	} 
	else
	{
		C2DmHumanoid(index_now + period);
	}
	//toc();
	// 0.24ms


	//tic();

	// LQSI�����ܼ� AI = Ad - (identity matrix 3-by-3)
	// 0.2ms --> 0.018ms original
	//tic();//release avg 0.02ms 600��/loop --> 0.00016ms 1��/loop
	if (index_now==0)
	{
		for (int i = index_now/*0*/; i <= index_now + period ; i++)
		{
			AIStk[i].data[0] = AdStk[i].data[0]-1;
			AIStk[i].data[1] = AdStk[i].data[1];
			AIStk[i].data[2] = AdStk[i].data[2];
			AIStk[i].data[3] = AdStk[i].data[3];
			AIStk[i].data[4] = AdStk[i].data[4]-1;
			AIStk[i].data[5] = AdStk[i].data[5];
			AIStk[i].data[6] = AdStk[i].data[6];
			AIStk[i].data[7] = AdStk[i].data[7];
			AIStk[i].data[8] = AdStk[i].data[8]-1;
		}
	} 
	else
	{
		AIStk[index_now + period].data[0] = AdStk[index_now + period].data[0]-1;
		AIStk[index_now + period].data[1] = AdStk[index_now + period].data[1];
		AIStk[index_now + period].data[2] = AdStk[index_now + period].data[2];
		AIStk[index_now + period].data[3] = AdStk[index_now + period].data[3];
		AIStk[index_now + period].data[4] = AdStk[index_now + period].data[4]-1;
		AIStk[index_now + period].data[5] = AdStk[index_now + period].data[5];
		AIStk[index_now + period].data[6] = AdStk[index_now + period].data[6];
		AIStk[index_now + period].data[7] = AdStk[index_now + period].data[7];
		AIStk[index_now + period].data[8] = AdStk[index_now + period].data[8]-1;
	}
	//toc();
	// 0.2ms

	// �p��LQSI Control Law~ 
	// 0.08ms original
	//tic();//release avg 0.048ms 600��/loop --> 0.00067ms 1��/loop
	if (index_now==0)//20141121 checkOK
	{
		for (int i = index_now/*0*/; i <= index_now + period ; i++)
		{
			// ��l equation
			//gTempM = BdT_stk[i]*Qx*BdStk[i]; 
			//invMkStk[i] = 1.0/(gTempM.data()[0]+R); // invMkStk(1,i) = 1/(Bd'*Qx*Bd+R);

			MatMulAtB(BdStk[i].data,3,1,Qx->data,3,3,CompuTemp);
			MatMulAB(CompuTemp,1,3,BdStk[i].data,3,1,CompuTemp2);
			invMkStk[i] = 1.0/(CompuTemp2[0]+R); // invMkStk(1,i) = 1/(Bd'*Qx*Bd+R);

			//temp11 = (BdStk[i]>Qx[0])*BdStk[i];
			//invMkStk[i] = 1.0/(temp11.data[0]+R); // invMkStk(1,i) = 1/(Bd'*Qx*Bd+R);
		}
	} 
	else
	{
		MatMulAtB(BdStk[index_now + period].data,3,1,Qx->data,3,3,CompuTemp);
		MatMulAB(CompuTemp,1,3,BdStk[index_now + period].data,3,1,CompuTemp2);
		invMkStk[index_now + period] = 1.0/(CompuTemp2[0]+R); // invMkStk(1,i) = 1/(Bd'*Qx*Bd+R);
	}
	//toc();
	//0.08ms


	// �p��LQSI Control Law~ 
	//�H�U��̲׭ȫe���i���@�p�ɤF*******20140918 Start
	//tic();//release avg  0.00074ms 1��/loop
	if (index_now==0)//20141121 check OK
	{
		for (int i = index_now/*0*/; i <= index_now + period ; i++)
		{
			// 0.1ms
			MatScalarMul(BdStk[i].data,3,invMkStk+i,CompuTemp);
			MatMulABt(CompuTemp,3,1,BdStk[i].data,3,1,BinvMBT_Stk[i].data);
			// 0.1ms
		}
	} 
	else
	{
		MatScalarMul(BdStk[index_now + period].data,3,invMkStk+(index_now + period),CompuTemp);
		MatMulABt(CompuTemp,3,1,BdStk[index_now + period].data,3,1,BinvMBT_Stk[index_now + period].data);
	}
	

	// �p��LQSI Control Law~ 
	double tempScale = 0;
	if (index_now==0)//Wk 20141121 check OK Fk and invFk
	{
		for (int i = index_now/*0*/; i <= index_now + period ; i++)
		{
			// 0.31ms
			MatMulAB(BinvMBT_Stk[i].data,3,3,Qx->data,3,3,CompuTemp);
			MatMulAB(CompuTemp,3,3,AIStk[i].data,3,3,CompuTemp2);
			MatMiuAB(AdStk[i].data,CompuTemp2,WkStk[i].data,9);

			for (int j=0;j<9;j++)
			{
				invWkStk[i].data[j] = WkStk[i].data[j];
				 
			}
			InvSqMat(invWkStk[i].data,3);

			WkTStk[i].data[0] = WkStk[i].data[0];
			WkTStk[i].data[1] = WkStk[i].data[3];
			WkTStk[i].data[2] = WkStk[i].data[6];
			WkTStk[i].data[3] = WkStk[i].data[1];
			WkTStk[i].data[4] = WkStk[i].data[4];
			WkTStk[i].data[5] = WkStk[i].data[7];
			WkTStk[i].data[6] = WkStk[i].data[2];
			WkTStk[i].data[7] = WkStk[i].data[5];
			WkTStk[i].data[8] = WkStk[i].data[8];
		}
	} 
	else
	{
		MatMulAB(BinvMBT_Stk[index_now + period].data,3,3,Qx->data,3,3,CompuTemp);
		MatMulAB(CompuTemp,3,3,AIStk[index_now + period].data,3,3,CompuTemp2);
		MatMiuAB(AdStk[index_now + period].data,CompuTemp2,WkStk[index_now + period].data,9);//Wk = Fk

		for (int j=0;j<9;j++)
		{
			invWkStk[index_now + period].data[j] = WkStk[index_now + period].data[j];
		}
		InvSqMat(invWkStk[index_now + period].data,3);

		WkTStk[index_now + period].data[0] = WkStk[index_now + period].data[0];
		WkTStk[index_now + period].data[1] = WkStk[index_now + period].data[3];
		WkTStk[index_now + period].data[2] = WkStk[index_now + period].data[6];
		WkTStk[index_now + period].data[3] = WkStk[index_now + period].data[1];
		WkTStk[index_now + period].data[4] = WkStk[index_now + period].data[4];
		WkTStk[index_now + period].data[5] = WkStk[index_now + period].data[7];
		WkTStk[index_now + period].data[6] = WkStk[index_now + period].data[2];
		WkTStk[index_now + period].data[7] = WkStk[index_now + period].data[5];
		WkTStk[index_now + period].data[8] = WkStk[index_now + period].data[8];
	}

	if (index_now==0)//Ek 
	{
		for (int i = index_now/*0*/; i <= index_now + period ; i++)
		{
			// Ss �즡
			//SsStk[i] = NkTStk[i]*Qx*NkStk[i]+((invMkStk[i]*invMkStk[i]*R)*(AIT_stk[i]*Qx*BdStk[i]))*(BdT_stk[i]*Qx*AIStk[i]) + CdT*Q*Cd + AdT_stk[i]*SsStk[i+1]*invDkStk[i]*WkStk[i];  Sk

			// ��� Ss, �٤U 50ms
			MatMulAtB(AIStk[i].data,3,3,Qx->data,3,3,CompuTemp);
			MatMulAB(CompuTemp,3,3,AIStk[i].data,3,3,CompuTemp2); // �⧹ N'QxN �s�b CompuTemp2 //20141117 Nk��AI  Aik Qx AikT

			MatMulAB(Qx->data,3,3,BdStk[i].data,3,1,CompuTemp);
			MatMulAtB(AIStk[i].data,3,3,CompuTemp,3,1,CompuTemp3);

			//tempScale = invMkStk[i]*invMkStk[i]*R;
			tempScale = -1*invMkStk[i];
			MatScalarMul(CompuTemp3,3,&tempScale,CompuTemp4);
			MatMulABt(CompuTemp4,3,1,CompuTemp3,3,1,CompuTemp); // �⧹���q�A�s�b CompuTemp

			//MatScalarMul(CompuTemp3,3,&tempScale,CompuTemp3);
			//MatMulAtB(BdStk[i].data(),3,1,Qx.data(),3,3,CompuTemp);
			//MatMulAB(CompuTemp,1,3,AIStk[i].data(),3,3,CompuTemp4);
			//MatMulAB(CompuTemp3,3,1,CompuTemp4,1,3,CompuTemp); // �⧹���q�A�s�b CompuTemp

			MatAddAB(CompuTemp,CompuTemp2,CompuTemp,9); // N'QN + ���q�A�s�b CompuTemp

			tempScale = Q;
			MatScalarMul(Cd->data,3,&tempScale,CompuTemp2);
			MatMulAtB(Cd->data,1,3,CompuTemp2,1,3,CompuTemp3); // �⧹ C'QC

			MatAddAB(CompuTemp,CompuTemp3,EkStk[i].data,9); // N'QN + ���q + C'QC�A�s�b CompuTemp  20141117�즹�⧹Ek
		}
	} 
	else
	{
		// Ss �즡
		//SsStk[i] = NkTStk[i]*Qx*NkStk[i]+((invMkStk[i]*invMkStk[i]*R)*(AIT_stk[i]*Qx*BdStk[i]))*(BdT_stk[i]*Qx*AIStk[i]) + CdT*Q*Cd + AdT_stk[i]*SsStk[i+1]*invDkStk[i]*WkStk[i];  Sk

		// ��� Ss, �٤U 50ms
		MatMulAtB(AIStk[index_now + period].data,3,3,Qx->data,3,3,CompuTemp);
		MatMulAB(CompuTemp,3,3,AIStk[index_now + period].data,3,3,CompuTemp2); // �⧹ N'QxN �s�b CompuTemp2 //20141117 Nk��AI  Aik Qx AikT

		MatMulAB(Qx->data,3,3,BdStk[index_now + period].data,3,1,CompuTemp);
		MatMulAtB(AIStk[index_now + period].data,3,3,CompuTemp,3,1,CompuTemp3);

		//tempScale = invMkStk[i]*invMkStk[i]*R;
		tempScale = -1*invMkStk[index_now + period];
		MatScalarMul(CompuTemp3,3,&tempScale,CompuTemp4);
		MatMulABt(CompuTemp4,3,1,CompuTemp3,3,1,CompuTemp); // �⧹���q�A�s�b CompuTemp

		//MatScalarMul(CompuTemp3,3,&tempScale,CompuTemp3);
		//MatMulAtB(BdStk[i].data(),3,1,Qx.data(),3,3,CompuTemp);
		//MatMulAB(CompuTemp,1,3,AIStk[i].data(),3,3,CompuTemp4);
		//MatMulAB(CompuTemp3,3,1,CompuTemp4,1,3,CompuTemp); // �⧹���q�A�s�b CompuTemp

		MatAddAB(CompuTemp,CompuTemp2,CompuTemp,9); // N'QN + ���q�A�s�b CompuTemp

		tempScale = Q;
		MatScalarMul(Cd->data,3,&tempScale,CompuTemp2);
		MatMulAtB(Cd->data,1,3,CompuTemp2,1,3,CompuTemp3); // �⧹ C'QC

		MatAddAB(CompuTemp,CompuTemp3,EkStk[index_now + period].data,9); // N'QN + ���q + C'QC�A�s�b CompuTemp  20141117�즹�⧹Ek
	}

	if (index_now==0)//�զXschur decomposition �ݭn���x�}
	{
		for (int i = index_now; i <= index_now + period ; i++)
		{
			//BinvMBT_Stk[i].data;
			MatMulAB(invWkStk[i].data,3,3,BinvMBT_Stk[i].data,3,3,CompuTemp);//inv(Fk)*Bk*invM*Bk'
			MatMulAB(EkStk[i].data,3,3,invWkStk[i].data,3,3,CompuTemp2);//Ek*inv(Fk)
			MatMulAB(EkStk[i].data,3,3,CompuTemp,3,3,CompuTemp3);
			MatAddAB(WkTStk[i].data,CompuTemp3,CompuTemp4,9);//Fk'+Ek*inv(Fk)*Bk*invM*Bk'
			//Z transpose
			Z_Schur[i].data[0] = invWkStk[i].data[0];//fill left up 3*3 and transpose
			Z_Schur[i].data[6] = invWkStk[i].data[1];
			Z_Schur[i].data[12] = invWkStk[i].data[2];
			Z_Schur[i].data[1] = invWkStk[i].data[3];
			Z_Schur[i].data[7] = invWkStk[i].data[4];
			Z_Schur[i].data[13] = invWkStk[i].data[5];
			Z_Schur[i].data[2] = invWkStk[i].data[6];
			Z_Schur[i].data[8] = invWkStk[i].data[7];
			Z_Schur[i].data[14] = invWkStk[i].data[8];
			Z_Schur[i].data[18] = CompuTemp[0];//fill right up 3*3
			Z_Schur[i].data[24] = CompuTemp[1];
			Z_Schur[i].data[30] = CompuTemp[2];
			Z_Schur[i].data[19] = CompuTemp[3];
			Z_Schur[i].data[25] = CompuTemp[4];
			Z_Schur[i].data[31] = CompuTemp[5];
			Z_Schur[i].data[20] = CompuTemp[6];
			Z_Schur[i].data[26] = CompuTemp[7];
			Z_Schur[i].data[32] = CompuTemp[8];
			Z_Schur[i].data[3] = CompuTemp2[0];//fill left down 3*3
			Z_Schur[i].data[9] = CompuTemp2[1];
			Z_Schur[i].data[15] = CompuTemp2[2];
			Z_Schur[i].data[4] = CompuTemp2[3];
			Z_Schur[i].data[10] = CompuTemp2[4];
			Z_Schur[i].data[16] = CompuTemp2[5];
			Z_Schur[i].data[5] = CompuTemp2[6];
			Z_Schur[i].data[11] = CompuTemp2[7];
			Z_Schur[i].data[17] = CompuTemp2[8];
			Z_Schur[i].data[21] = CompuTemp4[0];//fill right down 3*3
			Z_Schur[i].data[27] = CompuTemp4[1];
			Z_Schur[i].data[33] = CompuTemp4[2];
			Z_Schur[i].data[22] = CompuTemp4[3];
			Z_Schur[i].data[28] = CompuTemp4[4];
			Z_Schur[i].data[34] = CompuTemp4[5];
			Z_Schur[i].data[23] = CompuTemp4[6];
			Z_Schur[i].data[29] = CompuTemp4[7];
			Z_Schur[i].data[35] = CompuTemp4[8];
			//************
			dare_sqrmat_matlab(Z_Schur[i].data,Z_Schur_U[i].data,6);//U�O��m�L��
			CompuTemp[0] = Z_Schur_U[i].data[0];//U11
			CompuTemp[1] = Z_Schur_U[i].data[6];
			CompuTemp[2] = Z_Schur_U[i].data[12];
			CompuTemp[3] = Z_Schur_U[i].data[1];
			CompuTemp[4] = Z_Schur_U[i].data[7];
			CompuTemp[5] = Z_Schur_U[i].data[13];
			CompuTemp[6] = Z_Schur_U[i].data[2];
			CompuTemp[7] = Z_Schur_U[i].data[8];
			CompuTemp[8] = Z_Schur_U[i].data[14];

			CompuTemp2[0] = Z_Schur_U[i].data[3];//U21
			CompuTemp2[1] = Z_Schur_U[i].data[9];
			CompuTemp2[2] = Z_Schur_U[i].data[15];
			CompuTemp2[3] = Z_Schur_U[i].data[4];
			CompuTemp2[4] = Z_Schur_U[i].data[10];
			CompuTemp2[5] = Z_Schur_U[i].data[16];
			CompuTemp2[6] = Z_Schur_U[i].data[5];
			CompuTemp2[7] = Z_Schur_U[i].data[11];
			CompuTemp2[8] = Z_Schur_U[i].data[17];

			InvSqMat(CompuTemp,3);//inv U11
			MatMulAB(CompuTemp2,3,3,CompuTemp,3,3,SsStk[period/*i*/].data);//Sn of DARE  need?? 20150204 SS�ϥ�PERIOD���׼Ȧs
			//system("pause");
		}
	} 
	else
	{
		//BinvMBT_Stk[i].data;
			MatMulAB(invWkStk[index_now + period].data,3,3,BinvMBT_Stk[index_now + period].data,3,3,CompuTemp);//inv(Fk)*Bk*invM*Bk'
			MatMulAB(EkStk[index_now + period].data,3,3,invWkStk[index_now + period].data,3,3,CompuTemp2);//Ek*inv(Fk)
			MatMulAB(EkStk[index_now + period].data,3,3,CompuTemp,3,3,CompuTemp3);
			MatAddAB(WkTStk[index_now + period].data,CompuTemp3,CompuTemp4,9);//Fk'+Ek*inv(Fk)*Bk*invM*Bk'
			if (index_now==3000)
			{
				index_now=index_now;
			}
			//Z transpose
			Z_Schur[index_now + period].data[0] = invWkStk[index_now + period].data[0];//fill left up 3*3 and transpose
			Z_Schur[index_now + period].data[6] = invWkStk[index_now + period].data[1];
			Z_Schur[index_now + period].data[12] = invWkStk[index_now + period].data[2];
			Z_Schur[index_now + period].data[1] = invWkStk[index_now + period].data[3];
			Z_Schur[index_now + period].data[7] = invWkStk[index_now + period].data[4];
			Z_Schur[index_now + period].data[13] = invWkStk[index_now + period].data[5];
			Z_Schur[index_now + period].data[2] = invWkStk[index_now + period].data[6];
			Z_Schur[index_now + period].data[8] = invWkStk[index_now + period].data[7];
			Z_Schur[index_now + period].data[14] = invWkStk[index_now + period].data[8];
			Z_Schur[index_now + period].data[18] = CompuTemp[0];//fill right up 3*3
			Z_Schur[index_now + period].data[24] = CompuTemp[1];
			Z_Schur[index_now + period].data[30] = CompuTemp[2];
			Z_Schur[index_now + period].data[19] = CompuTemp[3];
			Z_Schur[index_now + period].data[25] = CompuTemp[4];
			Z_Schur[index_now + period].data[31] = CompuTemp[5];
			Z_Schur[index_now + period].data[20] = CompuTemp[6];
			Z_Schur[index_now + period].data[26] = CompuTemp[7];
			Z_Schur[index_now + period].data[32] = CompuTemp[8];
			Z_Schur[index_now + period].data[3] = CompuTemp2[0];//fill left down 3*3
			Z_Schur[index_now + period].data[9] = CompuTemp2[1];
			Z_Schur[index_now + period].data[15] = CompuTemp2[2];
			Z_Schur[index_now + period].data[4] = CompuTemp2[3];
			Z_Schur[index_now + period].data[10] = CompuTemp2[4];
			Z_Schur[index_now + period].data[16] = CompuTemp2[5];
			Z_Schur[index_now + period].data[5] = CompuTemp2[6];
			Z_Schur[index_now + period].data[11] = CompuTemp2[7];
			Z_Schur[index_now + period].data[17] = CompuTemp2[8];
			Z_Schur[index_now + period].data[21] = CompuTemp4[0];//fill right down 3*3
			Z_Schur[index_now + period].data[27] = CompuTemp4[1];
			Z_Schur[index_now + period].data[33] = CompuTemp4[2];
			Z_Schur[index_now + period].data[22] = CompuTemp4[3];
			Z_Schur[index_now + period].data[28] = CompuTemp4[4];
			Z_Schur[index_now + period].data[34] = CompuTemp4[5];
			Z_Schur[index_now + period].data[23] = CompuTemp4[6];
			Z_Schur[index_now + period].data[29] = CompuTemp4[7];
			Z_Schur[index_now + period].data[35] = CompuTemp4[8];
			if (index_now==3000)
			{
				index_now=index_now;
			}
			//************
			dare_sqrmat_matlab(Z_Schur[index_now + period].data,Z_Schur_U[index_now + period].data,6);//U�O��m�L��
			CompuTemp[0] = Z_Schur_U[index_now + period].data[0];//U11
			CompuTemp[1] = Z_Schur_U[index_now + period].data[6];
			CompuTemp[2] = Z_Schur_U[index_now + period].data[12];
			CompuTemp[3] = Z_Schur_U[index_now + period].data[1];
			CompuTemp[4] = Z_Schur_U[index_now + period].data[7];
			CompuTemp[5] = Z_Schur_U[index_now + period].data[13];
			CompuTemp[6] = Z_Schur_U[index_now + period].data[2];
			CompuTemp[7] = Z_Schur_U[index_now + period].data[8];
			CompuTemp[8] = Z_Schur_U[index_now + period].data[14];

			CompuTemp2[0] = Z_Schur_U[index_now + period].data[3];//U21
			CompuTemp2[1] = Z_Schur_U[index_now + period].data[9];
			CompuTemp2[2] = Z_Schur_U[index_now + period].data[15];
			CompuTemp2[3] = Z_Schur_U[index_now + period].data[4];
			CompuTemp2[4] = Z_Schur_U[index_now + period].data[10];
			CompuTemp2[5] = Z_Schur_U[index_now + period].data[16];
			CompuTemp2[6] = Z_Schur_U[index_now + period].data[5];
			CompuTemp2[7] = Z_Schur_U[index_now + period].data[11];
			CompuTemp2[8] = Z_Schur_U[index_now + period].data[17];
			/*if (index_now==3000)
			{
				index_now=index_now;
			}*/
			InvSqMat(CompuTemp,3);//inv U11
			MatMulAB(CompuTemp2,3,3,CompuTemp,3,3,SsStk[/*index_now + */period].data);//Sn of DARE

			/*if (index_now==3000)
			{
			index_now=index_now;
			}*/
	}
	
	//****Vk(:,period+1)=inv(eye(3)-Fk(:,:,i+period)'*(eye(3)-S(:,:,period+1)*invDk(:,:,period+1)*Bk(:,i+period)*invM(i+period)*Bk(:,i+period)'))*Ck(i+period,:)'*(Q*InpZMPx(i+period));%V301
	//****Dk(:,:,period+1)=eye(3)+Bk(:,i+period)*invM(i+period)*Bk(:,i+period)'*S(:,:,period+1);%D301
		//�p��Vn 20150205
		MatMulAB(BinvMBT_Stk[index_now + period].data,3,3,SsStk[/*index_now + */period].data,3,3,CompuTemp);//??
		MatAddAB(EYE3.data,CompuTemp,invDkStk[/*index_now + */period].data,9);
		InvSqMat(invDkStk[/*index_now + */period].data,3); //inv Dk 20141121 check OK
		MatMulAB(SsStk[/*index_now + */period].data,3,3,invDkStk[/*index_now + */period].data,3,3,CompuTemp);
		MatMulAB(CompuTemp,3,3,BinvMBT_Stk[index_now + period].data,3,3,CompuTemp2);
		MatMiuAB(EYE3.data,CompuTemp2,CompuTemp,9);
		MatMulAB(WkTStk[index_now + period].data,3,3,CompuTemp,3,3,CompuTemp3);//?transpose
		MatMiuAB(EYE3.data,CompuTemp3,CompuTemp,9);
		InvSqMat(CompuTemp,3);//inv(eye(3)-Fk(:,:,i+period)'*(eye(3)-S(:,:,period+1)*invDk(:,:,period+1)*Bk(:,i+period)*invM(i+period)*Bk(:,i+period)'))
		MatScalarMul(CdT->data,3,ZMPx[index_now + period],CompuTemp4);
		MatScalarMul(CompuTemp4,3,P,CompuTemp4);//Ck(i+period,:)'*(Q*InpZMPx(i+period))    
		MatMulAB(CompuTemp,3,3,CompuTemp4,3,1,vsStkX[/*index_now + */period].data);//Vx
		MatScalarMul(CdT->data,3,ZMPy[index_now + period],CompuTemp4);
		MatScalarMul(CompuTemp4,3,P,CompuTemp4);//Ck(i+period,:)'*(Q*InpZMPy(i+period))    
		MatMulAB(CompuTemp,3,3,CompuTemp4,3,1,vsStkY[/*index_now + */period].data);//Vy
		if (index_now==500)
			index_now=index_now;
		//vsStkX[(int)(index_now + period)].data[0] = CdT->data[0]*(P*ZMPx[(int)(index_now + period/*-1*/)]); // ���w�̫᪺��(backward riccati���}�l)//20141020 vsStkX[(int)(index_now + period)].data[0] = CdT->data[0]*(P*ZMPx[(int)(index_now + period-1)]); // ���w�̫᪺��(backward riccati���}�l)
		//vsStkX[(int)(index_now + period)].data[1] = CdT->data[1]*(P*ZMPx[(int)(index_now + period/*-1*/)]); // ���w�̫᪺��(backward riccati���}�l)
		//vsStkX[(int)(index_now + period)].data[2] = CdT->data[2]*(P*ZMPx[(int)(index_now + period/*-1*/)]); // ���w�̫᪺��(backward riccati���}�l) vn regulation final value
	//vsStkY[(int)(index_now + period)].data[0] = CdT->data[0]*(P*ZMPy[(int)(index_now + period/*-1*/)]); // ���w�̫᪺��(backward riccati���}�l)
	//vsStkY[(int)(index_now + period)].data[1] = CdT->data[1]*(P*ZMPy[(int)(index_now + period/*-1*/)]); // ���w�̫᪺��(backward riccati���}�l)
	//vsStkY[(int)(index_now + period)].data[2] = CdT->data[2]*(P*ZMPy[(int)(index_now + period/*-1*/)]); // ���w�̫᪺��(backward riccati���}�l)

	/*MatMulAB(BinvMBT_Stk[index_now + period].data,3,3,SsStk[index_now + period].data,3,3,DkStk[period+1].data);
	MatAddAB(EYE3.data,DkStk[period+1].data,9);*/


	// �p��LQSI Control Law~ 
	// 0.24ms
	//if (index_now==0)//20141117����Nk
	//{
	//	for (int i = index_now/*0*/; i < index_now + period ; i++)
	//	{
	//	   // NkStk[i] = WkStk[i] - EYE3;
	//		NkStk[i].data[0] = WkStk[i].data[0]-1;
	//		NkStk[i].data[1] = WkStk[i].data[1];
	//		NkStk[i].data[2] = WkStk[i].data[2];
	//		NkStk[i].data[3] = WkStk[i].data[3];
	//		NkStk[i].data[4] = WkStk[i].data[4]-1;
	//		NkStk[i].data[5] = WkStk[i].data[5];
	//		NkStk[i].data[6] = WkStk[i].data[6];
	//		NkStk[i].data[7] = WkStk[i].data[7];
	//		NkStk[i].data[8] = WkStk[i].data[8]-1;
	//	}
	//} 
	//else
	//{
	//	NkStk[index_now + period-1].data[0] = WkStk[index_now + period-1].data[0]-1;
	//	NkStk[index_now + period-1].data[1] = WkStk[index_now + period-1].data[1];
	//	NkStk[index_now + period-1].data[2] = WkStk[index_now + period-1].data[2];
	//	NkStk[index_now + period-1].data[3] = WkStk[index_now + period-1].data[3];
	//	NkStk[index_now + period-1].data[4] = WkStk[index_now + period-1].data[4]-1;
	//	NkStk[index_now + period-1].data[5] = WkStk[index_now + period-1].data[5];
	//	NkStk[index_now + period-1].data[6] = WkStk[index_now + period-1].data[6];
	//	NkStk[index_now + period-1].data[7] = WkStk[index_now + period-1].data[7];
	//	NkStk[index_now + period-1].data[8] = WkStk[index_now + period-1].data[8]-1;
	//}
	//toc();
	// 0.24ms 0.188
	//*******************20140918 end
	

	// �p��LQSI Control Law~ 20141121 check OK
	/////// �o�̥i�H�� �w����n���̲׭Ȩ��N ++period
	//tempScale = P;

	//vsStkX[(int)(index_now + period)].data[0] = CdT->data[0]*(P*ZMPx[(int)(index_now + period/*-1*/)]); // ���w�̫᪺��(backward riccati���}�l)//20141020 vsStkX[(int)(index_now + period)].data[0] = CdT->data[0]*(P*ZMPx[(int)(index_now + period-1)]); // ���w�̫᪺��(backward riccati���}�l)
	//vsStkX[(int)(index_now + period)].data[1] = CdT->data[1]*(P*ZMPx[(int)(index_now + period/*-1*/)]); // ���w�̫᪺��(backward riccati���}�l)
	//vsStkX[(int)(index_now + period)].data[2] = CdT->data[2]*(P*ZMPx[(int)(index_now + period/*-1*/)]); // ���w�̫᪺��(backward riccati���}�l) vn regulation final value

	//vsStkY[(int)(index_now + period)].data[0] = CdT->data[0]*(P*ZMPy[(int)(index_now + period/*-1*/)]); // ���w�̫᪺��(backward riccati���}�l)
	//vsStkY[(int)(index_now + period)].data[1] = CdT->data[1]*(P*ZMPy[(int)(index_now + period/*-1*/)]); // ���w�̫᪺��(backward riccati���}�l)
	//vsStkY[(int)(index_now + period)].data[2] = CdT->data[2]*(P*ZMPy[(int)(index_now + period/*-1*/)]); // ���w�̫᪺��(backward riccati���}�l)

	//MatScalarMul(Cd->data,3,&tempScale,CompuTemp);
	//MatMulAB(CdT->data,3,1,CompuTemp,1,3,SsStk[(int)(index_now + period)].data);//Sn final value  //20141020 MatMulAB(CdT->data,3,1,CompuTemp,1,3,SsStk[(int)(index_now + period)].data);//Sn final value

	/////// �o�̥i�H�� �w����n���̲׭Ȩ��N

	// �p��LQSI Control Law~ 
	// 88ms original
	for (int i = (int)(index_now + period-1) ; i >= index_now + 0 ; i--)//back period  20140823 ����i��n��state�@�_��s//20141020 for (int i = (int)(index_now + period-1) ; i >=index_now + 0 ; i--)
	{
		if (index_now==3000)
			index_now=index_now;

			//�B�zDARE
		
		// MATLAB code
		////invDkStk(:,:,i) = inv(eye(3)+BinvMBStk(:,:,i)*Ss(:,:,i+1));
		////Ss(:,:,i) = Nk'*Qx*Nk+(invMkStk(1,i)*invMkStk(1,i)*R)*(AI'*Qx*Bd*Bd'*Qx*AI)+CQC+Ad'*Ss(:,:,i+1)*invDkStk(:,:,i)*Wk;
		////vs(:,i) = Cd'*Q*ZMPr(i)+Ad'*(eye(3)-Ss(:,:,i+1)*invDkStk(:,:,i)*BinvMBStk(:,:,i))*vs(:,i+1);
		// MATLAB code

		//invDkStk[i] = EYE3+BinvMBT_Stk[i]*SsStk[i+1]; // �U����檺�즡		Dk
		MatMulAB(BinvMBT_Stk[i].data,3,3,SsStk[i+1-index_now].data,3,3,CompuTemp);
		MatAddAB(EYE3.data,CompuTemp,invDkStk[i-index_now].data,9);
		InvSqMat(invDkStk[i-index_now].data,3); // Take Inverse: call by reference ���ݭn�t�~�[�t inv Dk 20141121 check OK
		
		

		////// Ss �즡
		//////SsStk[i] = NkTStk[i]*Qx*NkStk[i]+((invMkStk[i]*invMkStk[i]*R)*(AIT_stk[i]*Qx*BdStk[i]))*(BdT_stk[i]*Qx*AIStk[i]) + CdT*Q*Cd + AdT_stk[i]*SsStk[i+1]*invDkStk[i]*WkStk[i];  Sk

		////// ��� Ss, �٤U 50ms
		////MatMulAtB(AIStk[i].data,3,3,Qx->data,3,3,CompuTemp);
		////MatMulAB(CompuTemp,3,3,AIStk[i].data,3,3,CompuTemp2); // �⧹ N'QxN �s�b CompuTemp2 //20141117 Nk��AI  Aik Qx AikT

		////MatMulAB(Qx->data,3,3,BdStk[i].data,3,1,CompuTemp);
		////MatMulAtB(AIStk[i].data,3,3,CompuTemp,3,1,CompuTemp3);

		//////tempScale = invMkStk[i]*invMkStk[i]*R;
		////tempScale = -1*invMkStk[i];
		////MatScalarMul(CompuTemp3,3,&tempScale,CompuTemp4);
		////MatMulABt(CompuTemp4,3,1,CompuTemp3,3,1,CompuTemp); // �⧹���q�A�s�b CompuTemp

		//////MatScalarMul(CompuTemp3,3,&tempScale,CompuTemp3);
		//////MatMulAtB(BdStk[i].data(),3,1,Qx.data(),3,3,CompuTemp);
		//////MatMulAB(CompuTemp,1,3,AIStk[i].data(),3,3,CompuTemp4);
		//////MatMulAB(CompuTemp3,3,1,CompuTemp4,1,3,CompuTemp); // �⧹���q�A�s�b CompuTemp

		////MatAddAB(CompuTemp,CompuTemp2,CompuTemp,9); // N'QN + ���q�A�s�b CompuTemp

		////tempScale = Q;
		////MatScalarMul(Cd->data,3,&tempScale,CompuTemp2);
		////MatMulAtB(Cd->data,1,3,CompuTemp2,1,3,CompuTemp3); // �⧹ C'QC

		////MatAddAB(CompuTemp,CompuTemp3,CompuTemp,9); // N'QN + ���q + C'QC�A�s�b CompuTemp  20141117�즹�⧹Ek

		MatMulAtB(WkStk[i].data,3,3,SsStk[i+1-index_now].data,3,3,CompuTemp2);//20141117 Ak to Wk
		MatMulAB(CompuTemp2,3,3,invDkStk[i-index_now].data,3,3,CompuTemp3);
		MatMulAB(CompuTemp3,3,3,WkStk[i].data,3,3,CompuTemp2);

		MatAddAB(EkStk[i].data/*CompuTemp*/,CompuTemp2,SsStk[i-index_now].data,9); // �⧹ Ss�A�s�b SsStk[i] 20141121 check OK
		// ��� Ss

		//// ��� vs �٤U30ms
		////vs_Stk[i] = CdT*Q*Input_ZMP[i]+AdT_stk[i]*(EYE3-SsStk[i+1]*invDkStk[i]*BinvMBT_Stk[i])*vs_Stk[i+1];

		//tempScale = Q*Input_ZMP[i];
		//MatScalarMul(Cd->data,3,&tempScale,CompuTemp); // C'Qr

		//MatMulAB(BinvMBT_Stk[i].data,3,3,vs_Stk[i+1].data,3,1,CompuTemp2);
		//MatMulAB(invDkStk[i].data,3,3,CompuTemp2,3,1,CompuTemp3);
		//MatMulAB(SsStk[i+1].data,3,3,CompuTemp3,3,1,CompuTemp2);

		//MatMiuAB(vs_Stk[i+1].data,CompuTemp2,CompuTemp3,3);
		//MatMulAtB(AdStk[i].data,3,3,CompuTemp3,3,1,CompuTemp2); // �⧹��b

		//MatAddAB(CompuTemp,CompuTemp2,vs_Stk[i].data,3);
		//// ��� vs

	}

	// for ZMPx
	for (int i = (int)(index_now + period-1) ; i >=index_now + 0 ; i--)//back period  20140823 ����i��n��state�@�_��s
	{
		// ��� vs �٤U30ms
		//vs_Stk[i] = CdT*Q*Input_ZMP[i]+AdT_stk[i]*(EYE3-SsStk[i+1]*invDkStk[i]*BinvMBT_Stk[i])*vs_Stk[i+1]; vkx
		if (index_now==3000 /*|| i==index_now*/)
			index_now=index_now;
		tempScale = Q*ZMPx[i];//regulation goal
		MatScalarMul(Cd->data,3,&tempScale,CompuTemp); // C'Qr

		MatMulAB(BinvMBT_Stk[i].data,3,3,vsStkX[i+1-index_now].data,3,1,CompuTemp2);
		MatMulAB(invDkStk[i-index_now].data,3,3,CompuTemp2,3,1,CompuTemp3);
		MatMulAB(SsStk[i+1-index_now].data,3,3,CompuTemp3,3,1,CompuTemp2);

		MatMiuAB(vsStkX[i+1-index_now].data,CompuTemp2,CompuTemp3,3);
		MatMulAtB(WkStk[i].data,3,3,CompuTemp3,3,1,CompuTemp2); // �⧹��b 20141117 Ak to Wk

		MatAddAB(CompuTemp,CompuTemp2,vsStkX[i-index_now].data,3);//20141121 check OK
		// ��� vs
		/*if(i==2000&&index_now==2000)
			system("pause");*/
	}

	// for ZMPy
	for (int i = (int)(index_now + period-1) ; i >=index_now + 0 ; i--)
	{
		// ��� vs �٤U30ms
		//vs_Stk[i] = CdT*Q*Input_ZMP[i]+AdT_stk[i]*(EYE3-SsStk[i+1]*invDkStk[i]*BinvMBT_Stk[i])*vs_Stk[i+1]; vky

		tempScale = Q*ZMPy[i];
		MatScalarMul(Cd->data,3,&tempScale,CompuTemp); // C'Qr

		MatMulAB(BinvMBT_Stk[i].data,3,3,vsStkY[i+1-index_now].data,3,1,CompuTemp2);
		MatMulAB(invDkStk[i-index_now].data,3,3,CompuTemp2,3,1,CompuTemp3);
		MatMulAB(SsStk[i+1-index_now].data,3,3,CompuTemp3,3,1,CompuTemp2);

		MatMiuAB(vsStkY[i+1-index_now].data,CompuTemp2,CompuTemp3,3);
		MatMulAtB(WkStk[i].data,3,3,CompuTemp3,3,1,CompuTemp2); // �⧹��b 20141117 Ak to Wk

		MatAddAB(CompuTemp,CompuTemp2,vsStkY[i-index_now].data,3);//20141121 check OK
		// ��� vs
		/*if(i==2000&&index_now==2000)
			system("pause");*/
	}
	//**********Sk test
	/*for (int i=0;i<period;i++)
	{
		for(int j=0;j<9;j++)
			Sk_test[index_now][i].data[j]=SsStk[i].data[j];
	}*/
	//if (index_now==2000)
	//{
	//	fstream Sk;
	//	/*Sk.open("Sk.txt",ios::out);
	//	Sk.precision(10);
	//	for (int i=index_now+1 ; i<= index_now+1200 ; i++)
	//	{
	//		Sk << SsStk[i].data[0] << " " << SsStk[i].data[1] << " " << SsStk[i].data[2] << " " << SsStk[i].data[3] <<" " << SsStk[i].data[4] <<" " << SsStk[i].data[5] <<" " << SsStk[i].data[6] <<" " << SsStk[i].data[7] <<" " << SsStk[i].data[8] << endl;
	//	}
	//	Sk.close();*/
	//}

}


void LQSISolver::DummyControl(void)
{
	/******************************************************************
	input: void
	output: void

	Note:
	// DummyControl �N�Ҧ��w�g��n��state matrices and control laws �Ӻ�X�������߭y��
	// �Ҧ���feedback���ϥ��|�N�����G �����i�H��X����reference ������COG�y��
	******************************************************************/

	YMatLite EYE3(3,3);
	EYE3.data[0] = 1;
	EYE3.data[1] = 0;
	EYE3.data[2] = 0;
	EYE3.data[3] = 0;
	EYE3.data[4] = 1;
	EYE3.data[5] = 0;
	EYE3.data[6] = 0;
	EYE3.data[7] = 0;
	EYE3.data[8] = 1;

	// state variable
	XState[0].data[0] = 0;
	XState[0].data[1] = 0;
	XState[0].data[2] = 0;

	YMatLite Bkuk(3,1);

	for (int i = 0 ; i<LQDataLen ; i++)
	{
		// sample matlab code
		// Bkuk = BinvMBStk(:,:,i)*(-(Qx*AI+Ss(:,:,i+1)*invDkStk(:,:,i)*WkStk(:,:,i))*x+(eye(3)-Ss(:,:,i+1)*invDkStk(:,:,i)*BinvMBStk(:,:,i))*vs(:,i+1)); 	
		// x = AkStk(:,:,i)*x + Bkuk; 

		//Bkuk = BinvMBT_Stk[i]*(-(Qx*AIStk[i]+SsStk[i+1]*invDkStk[i]*WkStk[i])*XState[i]+(EYE3-SsStk[i+1]*invDkStk[i]*BinvMBT_Stk[i])*vs_Stk[i+1]);
		//XState[i+1] = AdStk[i]*XState[i] + Bkuk;

		MatMulAB(WkStk[i].data,3,3,XState[i].data,3,1,CompuTemp); // Wx
		MatMulAB(invDkStk[i].data,3,3,CompuTemp,3,1,CompuTemp2); // inv(D)Wx
		MatMulAB(SsStk[i+1].data,3,3,CompuTemp2,3,1,CompuTemp); // Sinv(D)Wx

		MatMulAB(AIStk[i].data,3,3,XState[i].data,3,1,CompuTemp2); // AIx
		MatMulAB(Qx->data,3,3,CompuTemp2,3,1,CompuTemp3); // QxAIx

		MatAddAB(CompuTemp,CompuTemp3,CompuTemp2,3); // CompuTemp2 = Sinv(D)Wx+QxAIx

		// get (EYE3-SsStk[i+1]*invDkStk[i]*BinvMBT_Stk[i])*vs_Stk[i+1]
		MatMulAB(BinvMBT_Stk[i].data,3,3,vsStkX[i+1].data,3,1,CompuTemp); // BMBv
		MatMulAB(invDkStk[i].data,3,3,CompuTemp,3,1,CompuTemp3); // inv(D)BMBv
		MatMulAB(SsStk[i+1].data,3,3,CompuTemp3,3,1,CompuTemp); // Sinv(D)BMBv

		MatMiuAB(vsStkX[i+1].data,CompuTemp,CompuTemp,3); // v-Sinv(D)BMBv
		MatMiuAB(CompuTemp,CompuTemp2,CompuTemp,3); // (v-Sinv(D)BMBv) - (Sinv(D)Wx+QxAIx)
		
		MatMulAB(BinvMBT_Stk[i].data,3,3,CompuTemp,3,1,Bkuk.data); // Bkuk found

		MatMulAB(AdStk[i].data,3,3,XState[i].data,3,1,CompuTemp);
		//MatAddAB(CompuTemp,Bkuk.data(),XState[i+1].data(),3);
		MatAddAB(CompuTemp,Bkuk.data,CompuTemp,3);

		// xk+1 = Axk+Bkuk
		MatMulAB(AdStk[i].data,3,3,XState[i].data,3,1,CompuTemp);
		MatAddAB(CompuTemp,Bkuk.data,XState[i+1].data,3);
	}


	for (int i = 0 ; i<LQDataLen ; i++)
	{
		// sample matlab code
		// Bkuk = BinvMBStk(:,:,i)*(-(Qx*AI+Ss(:,:,i+1)*invDkStk(:,:,i)*WkStk(:,:,i))*x+(eye(3)-Ss(:,:,i+1)*invDkStk(:,:,i)*BinvMBStk(:,:,i))*vs(:,i+1)); 	
		// x = AkStk(:,:,i)*x + Bkuk; 

		//Bkuk = BinvMBT_Stk[i]*(-(Qx*AIStk[i]+SsStk[i+1]*invDkStk[i]*WkStk[i])*XState[i]+(EYE3-SsStk[i+1]*invDkStk[i]*BinvMBT_Stk[i])*vs_Stk[i+1]);
		//XState[i+1] = AdStk[i]*XState[i] + Bkuk;

		MatMulAB(WkStk[i].data,3,3,YState[i].data,3,1,CompuTemp); // Wx
		MatMulAB(invDkStk[i].data,3,3,CompuTemp,3,1,CompuTemp2); // inv(D)Wx
		MatMulAB(SsStk[i+1].data,3,3,CompuTemp2,3,1,CompuTemp); // Sinv(D)Wx

		MatMulAB(AIStk[i].data,3,3,YState[i].data,3,1,CompuTemp2); // AIx
		MatMulAB(Qx->data,3,3,CompuTemp2,3,1,CompuTemp3); // QxAIx

		MatAddAB(CompuTemp,CompuTemp3,CompuTemp2,3); // CompuTemp2 = Sinv(D)Wx+QxAIx

		// get (EYE3-SsStk[i+1]*invDkStk[i]*BinvMBT_Stk[i])*vs_Stk[i+1]
		MatMulAB(BinvMBT_Stk[i].data,3,3,vsStkY[i+1].data,3,1,CompuTemp); // BMBv
		MatMulAB(invDkStk[i].data,3,3,CompuTemp,3,1,CompuTemp3); // inv(D)BMBv
		MatMulAB(SsStk[i+1].data,3,3,CompuTemp3,3,1,CompuTemp); // Sinv(D)BMBv

		MatMiuAB(vsStkY[i+1].data,CompuTemp,CompuTemp,3); // v-Sinv(D)BMBv
		MatMiuAB(CompuTemp,CompuTemp2,CompuTemp,3); // (v-Sinv(D)BMBv) - (Sinv(D)Wx+QxAIx)
		
		MatMulAB(BinvMBT_Stk[i].data,3,3,CompuTemp,3,1,Bkuk.data); // Bkuk found

		MatMulAB(AdStk[i].data,3,3,YState[i].data,3,1,CompuTemp);
		//MatAddAB(CompuTemp,Bkuk.data(),XState[i+1].data(),3);
		MatAddAB(CompuTemp,Bkuk.data,CompuTemp,3);

		// xk+1 = Axk+Bkuk
		MatMulAB(AdStk[i].data,3,3,YState[i].data,3,1,CompuTemp);
		MatAddAB(CompuTemp,Bkuk.data,YState[i+1].data,3);

	}

	#if SaveLQSIStates
		// �nprint���G�X�Ӫ��ܡA�Х��}�U��

		fstream Fx;

		Fx.open("statesx.txt",ios::out);
		Fx.precision(10);
		for (int i=0 ; i<= LQDataLen ; i++)
		{
			Fx << XState[i].data[0] << " " << XState[i].data[1] << " " << XState[i].data[2] <<  endl;
		}

		Fx.close();

		fstream Fy;

		Fy.open("statesy.txt",ios::out);
		Fy.precision(10);
		for (int i=0 ; i<= LQDataLen ; i++)
		{
			Fy << YState[i].data[0] << " " << YState[i].data[1] << " " << YState[i].data[2] <<  endl;
		}

		Fy.close();

		fstream Ak5;

		Ak5.open("Ak5.txt",ios::out);
		Ak5.precision(10);
		for (int i=0 ; i<= LQDataLen ; i++)
		{
			Ak5 << AdStk[i].data[0]<< " " << AdStk[i].data[1]<< " " << AdStk[i].data[2]<< " " << AdStk[i].data[3] << " "<< AdStk[i].data[4] << " "<< AdStk[i].data[5]<< " " << AdStk[i].data[6] << " "<< AdStk[i].data[7] << " "<< AdStk[i].data[8] << endl;
		}

		Ak5.close();

		fstream vk;

		vk.open("vkx.txt",ios::out);
		vk.precision(10);
		for (int i=0 ; i<= LQDataLen ; i++)
		{
			vk << vsStkX[i].data[0] << " " << vsStkX[i].data[1] << " " << vsStkX[i].data[2] <<  endl;
		}

		vk.close();
		// �nprint���G�X�Ӫ��ܡA�Х��}�W��
	#endif
}

void LQSISolver::DummyControl2( int index_now)
{
	/******************************************************************
	input: void
	output: void

	Note:
	// DummyControl �N�Ҧ��w�g��n��state matrices and control laws �Ӻ�X�������߭y��
	// �Ҧ���feedback���ϥ��|�N�����G �����i�H��X����reference ������COG�y��
	******************************************************************/

	YMatLite EYE3(3,3);
	EYE3.data[0] = 1;
	EYE3.data[1] = 0;
	EYE3.data[2] = 0;
	EYE3.data[3] = 0;
	EYE3.data[4] = 1;
	EYE3.data[5] = 0;
	EYE3.data[6] = 0;
	EYE3.data[7] = 0;
	EYE3.data[8] = 1;
	YMatLite Bkuk(3,1);
	// state variable
	if (index_now==0)
	{
		XState[0].data[0] = 0;//original ideal state
		XState[0].data[1] = 0;
		XState[0].data[2] = 0;
		YState[0].data[0] = 0;
		YState[0].data[1] = 0;
		YState[0].data[2] = 0;
		XState2[0].data[0] = 0;//now state �i�����~���ܤ�
		XState2[0].data[1] = 0;
		XState2[0].data[2] = 0;
		YState2[0].data[0] = 0;
		YState2[0].data[1] = 0;
		YState2[0].data[2] = 0;
	} 
	else
	{
		if(gFlagSimulation == ADAMSSimu)
		{
		//Kine.IMUSensorData(index_now,);
			XState2[index_now].data[0] = XState[index_now].data[0];//*(IMU1.posx + index_now+1000-2);//
			XState2[index_now].data[1] = XState[index_now].data[1];
			XState2[index_now].data[2] = gKineAll.FS_ZMP[(index_now+1000-1)*2+1];//XState[index_now].data[2];//
			YState2[index_now].data[0] = YState[index_now].data[0];//*(IMU1.posy + index_now+1000-2);//
			YState2[index_now].data[1] = YState[index_now].data[1];
			YState2[index_now].data[2] = gKineAll.FS_ZMP[(index_now+1000-1)*2];//YState[index_now].data[2];//
		}
		else
		{
			XState2[index_now].data[0] = XState[index_now].data[0];//�������ӭn���_�^�Ƕi�ӧ�sstate2�A���B���N�W�@���z�Q�ȷ������{���
			XState2[index_now].data[1] = XState[index_now].data[1];
			XState2[index_now].data[2] = XState[index_now].data[2];
			YState2[index_now].data[0] = YState[index_now].data[0];
			YState2[index_now].data[1] = YState[index_now].data[1];
			YState2[index_now].data[2] = YState[index_now].data[2];
			
		}
		
		/*if (index_now%200==10)
		{
			XState2[index_now].data[2]+=20;
		}else if (index_now%200==70)
		{
			XState2[index_now].data[2]-=20;
		}*/

	}
	//for (int i = 0 ; i<index_now/*LQDataLen*/ ; i++)
	//{
		// sample matlab code
		// Bkuk = BinvMBStk(:,:,i)*(-(Qx*AI+Ss(:,:,i+1)*invDkStk(:,:,i)*WkStk(:,:,i))*x+(eye(3)-Ss(:,:,i+1)*invDkStk(:,:,i)*BinvMBStk(:,:,i))*vs(:,i+1)); 	
		// x = AkStk(:,:,i)*x + Bkuk; 

		//Bkuk = BinvMBT_Stk[i]*(-(Qx*AIStk[i]+SsStk[i+1]*invDkStk[i]*WkStk[i])*XState[i]+(EYE3-SsStk[i+1]*invDkStk[i]*BinvMBT_Stk[i])*vs_Stk[i+1]);
		//XState[i+1] = AdStk[i]*XState[i] + Bkuk;

		MatMulAB(WkStk[index_now].data,3,3,XState2[index_now].data,3,1,CompuTemp); // Wx  nowX
		MatMulAB(invDkStk[index_now].data,3,3,CompuTemp,3,1,CompuTemp2); // inv(D)Wx
		MatMulAB(SsStk[index_now+1].data,3,3,CompuTemp2,3,1,CompuTemp); // Sinv(D)Wx

		MatMulAB(AIStk[index_now].data,3,3,XState2[index_now].data,3,1,CompuTemp2); // AIx  nowX
		MatMulAB(Qx->data,3,3,CompuTemp2,3,1,CompuTemp3); // QxAIx

		MatAddAB(CompuTemp,CompuTemp3,CompuTemp2,3); // CompuTemp2 = Sinv(D)Wx+QxAIx

		// get (EYE3-SsStk[i+1]*invDkStk[i]*BinvMBT_Stk[i])*vs_Stk[i+1]
		MatMulAB(BinvMBT_Stk[index_now].data,3,3,vsStkX[index_now+1].data,3,1,CompuTemp); // BMBv
		MatMulAB(invDkStk[index_now].data,3,3,CompuTemp,3,1,CompuTemp3); // inv(D)BMBv
		MatMulAB(SsStk[index_now+1].data,3,3,CompuTemp3,3,1,CompuTemp); // Sinv(D)BMBv

		MatMiuAB(vsStkX[index_now+1].data,CompuTemp,CompuTemp,3); // v-Sinv(D)BMBv  ??�P�@�Ȧs��
		MatMiuAB(CompuTemp,CompuTemp2,CompuTemp,3); // (v-Sinv(D)BMBv) - (Sinv(D)Wx+QxAIx)
		
		MatMulAB(BinvMBT_Stk[index_now].data,3,3,CompuTemp,3,1,Bkuk.data); // Bkuk found
		

		//MatMulAB(AdStk[index_now].data,3,3,XState[index_now].data,3,1,CompuTemp);//angel ���ƤF?
		//MatAddAB(CompuTemp,Bkuk.data(),XState[i+1].data(),3);
		//MatAddAB(CompuTemp,Bkuk.data,CompuTemp,3);//angel ���ƤF?

		// xk+1 = Axk+Bkuk
		MatMulAB(AdStk[index_now].data,3,3,XState2[index_now].data,3,1,CompuTemp);
		MatAddAB(CompuTemp,Bkuk.data,XState[index_now+1].data,3);//estimate next ideal state
	//}
		//********20140916�t��í�w���� DORA
		/*for (int i=0;i<3;i++)
		{
		test01[i][index_now] = CompuTemp[i];
		test02[i][index_now] = Bkuk.data[i];
		test03[i][index_now] = CompuTemp[i] + Bkuk.data[i];
		}

		if (index_now%100==0)
		{
		index_now=index_now;
		}*/
		//********20140916�t��í�w���� DORA END


	//for (int i = 0 ; i<index_now/*LQDataLen*/ ; i++)
	//{
		// sample matlab code
		// Bkuk = BinvMBStk(:,:,i)*(-(Qx*AI+Ss(:,:,i+1)*invDkStk(:,:,i)*WkStk(:,:,i))*x+(eye(3)-Ss(:,:,i+1)*invDkStk(:,:,i)*BinvMBStk(:,:,i))*vs(:,i+1)); 	
		// x = AkStk(:,:,i)*x + Bkuk; 

		//Bkuk = BinvMBT_Stk[i]*(-(Qx*AIStk[i]+SsStk[i+1]*invDkStk[i]*WkStk[i])*XState[i]+(EYE3-SsStk[i+1]*invDkStk[i]*BinvMBT_Stk[i])*vs_Stk[i+1]);
		//XState[i+1] = AdStk[i]*XState[i] + Bkuk;

		MatMulAB(WkStk[index_now].data,3,3,YState2[index_now].data,3,1,CompuTemp); // Wx
		MatMulAB(invDkStk[index_now].data,3,3,CompuTemp,3,1,CompuTemp2); // inv(D)Wx
		MatMulAB(SsStk[index_now+1].data,3,3,CompuTemp2,3,1,CompuTemp); // Sinv(D)Wx

		MatMulAB(AIStk[index_now].data,3,3,YState2[index_now].data,3,1,CompuTemp2); // AIx
		MatMulAB(Qx->data,3,3,CompuTemp2,3,1,CompuTemp3); // QxAIx

		MatAddAB(CompuTemp,CompuTemp3,CompuTemp2,3); // CompuTemp2 = Sinv(D)Wx+QxAIx

		// get (EYE3-SsStk[i+1]*invDkStk[i]*BinvMBT_Stk[i])*vs_Stk[i+1]
		MatMulAB(BinvMBT_Stk[index_now].data,3,3,vsStkY[index_now+1].data,3,1,CompuTemp); // BMBv
		MatMulAB(invDkStk[index_now].data,3,3,CompuTemp,3,1,CompuTemp3); // inv(D)BMBv
		MatMulAB(SsStk[index_now+1].data,3,3,CompuTemp3,3,1,CompuTemp); // Sinv(D)BMBv

		MatMiuAB(vsStkY[index_now+1].data,CompuTemp,CompuTemp,3); // v-Sinv(D)BMBv
		MatMiuAB(CompuTemp,CompuTemp2,CompuTemp,3); // (v-Sinv(D)BMBv) - (Sinv(D)Wx+QxAIx)
		
		MatMulAB(BinvMBT_Stk[index_now].data,3,3,CompuTemp,3,1,Bkuk.data); // Bkuk found

		//MatMulAB(AdStk[index_now].data,3,3,YState[index_now].data,3,1,CompuTemp);//angel ���ƤF?
		//MatAddAB(CompuTemp,Bkuk.data(),XState[i+1].data(),3);
		//MatAddAB(CompuTemp,Bkuk.data,CompuTemp,3);//angel ���ƤF?

		// xk+1 = Axk+Bkuk
		MatMulAB(AdStk[index_now].data,3,3,YState2[index_now].data,3,1,CompuTemp);
		MatAddAB(CompuTemp,Bkuk.data,YState[index_now+1].data,3);

	//}
		//******test for Sk
		/*if (index_now==2000)
		{
			fstream Sk;
			Sk.open("Sk.txt",ios::out);
			Sk.precision(10);
			for (int i=index_now ; i<= index_now+300 ; i++)
			{
				Sk << SsStk[i].data[0] << " " << SsStk[i].data[1] << " " << SsStk[i].data[2] << " " << SsStk[i].data[3] <<" " << SsStk[i].data[4] <<" " << SsStk[i].data[5] <<" " << SsStk[i].data[6] <<" " << SsStk[i].data[7] <<" " << SsStk[i].data[8] << endl;
			}
			Sk.close();
		}*/
	//#if SaveLQSIStates
	//	// �nprint���G�X�Ӫ��ܡA�Х��}�U��
		/*if (index_now==7215)
		{

			fstream Fx;

			Fx.open("statesx_ideal.txt",ios::out);
			Fx.precision(10);
			for (int i=0 ; i<= LQDataLen ; i++)
			{
				Fx << XState[i].data[0] << " " << XState[i].data[1] << " " << XState[i].data[2] <<  endl;
			}

			Fx.close();
			Fx.open("statesx_adams.txt",ios::out);
			Fx.precision(10);
			for (int i=0 ; i<= LQDataLen ; i++)
			{
				Fx << XState2[i].data[0] << " " << XState2[i].data[1] << " " << XState2[i].data[2] <<  endl;
			}
			Fx.close();

			fstream Fy;

			Fy.open("statesy_ideal.txt",ios::out);
			Fy.precision(10);
			for (int i=0 ; i<= LQDataLen ; i++)
			{
				Fy << YState[i].data[0] << " " << YState[i].data[1] << " " << YState[i].data[2] <<  endl;
			}
			Fy.close();
			Fy.open("statesy_adams.txt",ios::out);
			Fy.precision(10);
			for (int i=0 ; i<= LQDataLen ; i++)
			{
				Fy << YState2[i].data[0] << " " << YState2[i].data[1] << " " << YState2[i].data[2] <<  endl;
			}
		}*/
	//	// �nprint���G�X�Ӫ��ܡA�Х��}�W��
	//#endif
}
void LQSISolver::DummyControl3( int index_now)//20141117 �PDummyControl2
{
	/******************************************************************
	input: void
	output: void

	Note:
	// DummyControl �N�Ҧ��w�g��n��state matrices and control laws �Ӻ�X�������߭y��
	// �Ҧ���feedback���ϥ��|�N�����G �����i�H��X����reference ������COG�y��
	******************************************************************/

	YMatLite EYE3(3,3);
	EYE3.data[0] = 1;
	EYE3.data[1] = 0;
	EYE3.data[2] = 0;
	EYE3.data[3] = 0;
	EYE3.data[4] = 1;
	EYE3.data[5] = 0;
	EYE3.data[6] = 0;
	EYE3.data[7] = 0;
	EYE3.data[8] = 1;
	YMatLite Bkuk(3,1);
	// state variable
	if (index_now==0)
	{
		XState[0].data[0] = 0;//original ideal state
		XState[0].data[1] = 0;
		XState[0].data[2] = 0;
		YState[0].data[0] = 0;
		YState[0].data[1] = 0;
		YState[0].data[2] = 0;
		XState2[0].data[0] = 0;//now state �i�����~���ܤ�
		XState2[0].data[1] = 0;
		XState2[0].data[2] = 0;
		YState2[0].data[0] = 0;
		YState2[0].data[1] = 0;
		YState2[0].data[2] = 0;
	} 
	else
	{
		if(gFlagSimulation == ADAMSSimu)
		{
		//Kine.IMUSensorData(index_now,);
			XState2[index_now].data[0] = XState[index_now].data[0];//*(IMU1.posx + index_now+1000-2);//
			XState2[index_now].data[1] = XState[index_now].data[1];
			XState2[index_now].data[2] = gKineAll.FS_ZMP[(index_now+1000-1)*2+1];//XState[index_now].data[2];//�쥻index_now���_�ơA�ݦA�T�{20150301
			/*if (index_now<400)
			{
				XState2[index_now].data[2]+=1.165;
			}*/
 			YState2[index_now].data[0] = YState[index_now].data[0];//*(IMU1.posy + index_now+1000-2);//
			YState2[index_now].data[1] = YState[index_now].data[1];
			YState2[index_now].data[2] = YState[index_now].data[2];//gKineAll.FS_ZMP[(index_now+1000-1)*2];//
		}
		else
		{
			XState2[index_now].data[0] = XState[index_now].data[0];//�������ӭn���_�^�Ƕi�ӧ�sstate2�A���B���N�W�@���z�Q�ȷ������{���
			XState2[index_now].data[1] = XState[index_now].data[1];
			XState2[index_now].data[2] = XState[index_now].data[2];
			YState2[index_now].data[0] = YState[index_now].data[0];
			YState2[index_now].data[1] = YState[index_now].data[1];
			YState2[index_now].data[2] = YState[index_now].data[2];
			
		}
		
		/*if (index_now%50==20)
		{
			XState2[index_now].data[2]+=30;
		}else if (index_now%50==7)
		{
			XState2[index_now].data[2]-=30;
		}*/

	}
	//for (int i = 0 ; i<index_now/*LQDataLen*/ ; i++)
	//{
		// sample matlab code
		// Bkuk = BinvMBStk(:,:,i)*(-(Qx*AI+Ss(:,:,i+1)*invDkStk(:,:,i)*WkStk(:,:,i))*x+(eye(3)-Ss(:,:,i+1)*invDkStk(:,:,i)*BinvMBStk(:,:,i))*vs(:,i+1)); 	
		// x = AkStk(:,:,i)*x + Bkuk; 

		//Bkuk = BinvMBT_Stk[i]*(-(Qx*AIStk[i]+SsStk[i+1]*invDkStk[i]*WkStk[i])*XState[i]+(EYE3-SsStk[i+1]*invDkStk[i]*BinvMBT_Stk[i])*vs_Stk[i+1]);
		//XState[i+1] = AdStk[i]*XState[i] + Bkuk;

		MatMulAB(WkStk[index_now].data,3,3,XState2[index_now].data,3,1,CompuTemp); // Wx  nowX
		MatMulAB(invDkStk[index_now-index_now].data,3,3,CompuTemp,3,1,CompuTemp2); // inv(D)Wx
		MatMulAB(SsStk[index_now+1-index_now].data,3,3,CompuTemp2,3,1,CompuTemp); // Sinv(D)Wx

		MatMulAB(AIStk[index_now].data,3,3,XState2[index_now].data,3,1,CompuTemp2); // AIx  nowX
		MatMulAB(Qx->data,3,3,CompuTemp2,3,1,CompuTemp3); // QxAIx

		MatAddAB(CompuTemp,CompuTemp3,CompuTemp2,3); // CompuTemp2 = Sinv(D)Wx+QxAIx

		// get (EYE3-SsStk[i+1]*invDkStk[i]*BinvMBT_Stk[i])*vs_Stk[i+1]
		MatMulAB(BinvMBT_Stk[index_now].data,3,3,vsStkX[index_now+1-index_now].data,3,1,CompuTemp); // BMBv
		MatMulAB(invDkStk[index_now-index_now].data,3,3,CompuTemp,3,1,CompuTemp3); // inv(D)BMBv
		MatMulAB(SsStk[index_now+1-index_now].data,3,3,CompuTemp3,3,1,CompuTemp); // Sinv(D)BMBv

		MatMiuAB(vsStkX[index_now+1-index_now].data,CompuTemp,CompuTemp,3); // v-Sinv(D)BMBv  ??�P�@�Ȧs��
		MatMiuAB(CompuTemp,CompuTemp2,CompuTemp,3); // (v-Sinv(D)BMBv) - (Sinv(D)Wx+QxAIx)
		
		MatMulAB(BinvMBT_Stk[index_now].data,3,3,CompuTemp,3,1,Bkuk.data); // Bkuk found
		
		/*if (index_now==3000)
		{
			index_now=index_now;
		}*/
		//MatMulAB(AdStk[index_now].data,3,3,XState[index_now].data,3,1,CompuTemp);//angel ���ƤF?
		//MatAddAB(CompuTemp,Bkuk.data(),XState[i+1].data(),3);
		//MatAddAB(CompuTemp,Bkuk.data,CompuTemp,3);//angel ���ƤF?

		// xk+1 = Axk+Bkuk
		MatMulAB(AdStk[index_now].data,3,3,XState2[index_now].data,3,1,CompuTemp);
		MatAddAB(CompuTemp,Bkuk.data,XState[index_now+1].data,3);//estimate next ideal state
	//}
		//********20140916�t��í�w���� DORA
		/*for (int i=0;i<3;i++)
		{
		test01[i][index_now] = CompuTemp[i];
		test02[i][index_now] = Bkuk.data[i];
		test03[i][index_now] = CompuTemp[i] + Bkuk.data[i];
		}

		if (index_now%100==0)
		{
		index_now=index_now;
		}*/
		//********20140916�t��í�w���� DORA END


	//for (int i = 0 ; i<index_now/*LQDataLen*/ ; i++)
	//{
		// sample matlab code
		// Bkuk = BinvMBStk(:,:,i)*(-(Qx*AI+Ss(:,:,i+1)*invDkStk(:,:,i)*WkStk(:,:,i))*x+(eye(3)-Ss(:,:,i+1)*invDkStk(:,:,i)*BinvMBStk(:,:,i))*vs(:,i+1)); 	
		// x = AkStk(:,:,i)*x + Bkuk; 

		//Bkuk = BinvMBT_Stk[i]*(-(Qx*AIStk[i]+SsStk[i+1]*invDkStk[i]*WkStk[i])*XState[i]+(EYE3-SsStk[i+1]*invDkStk[i]*BinvMBT_Stk[i])*vs_Stk[i+1]);
		//XState[i+1] = AdStk[i]*XState[i] + Bkuk;

		MatMulAB(WkStk[index_now].data,3,3,YState2[index_now].data,3,1,CompuTemp); // Wx
		MatMulAB(invDkStk[index_now-index_now].data,3,3,CompuTemp,3,1,CompuTemp2); // inv(D)Wx
		MatMulAB(SsStk[index_now+1-index_now].data,3,3,CompuTemp2,3,1,CompuTemp); // Sinv(D)Wx

		MatMulAB(AIStk[index_now].data,3,3,YState2[index_now].data,3,1,CompuTemp2); // AIx
		MatMulAB(Qx->data,3,3,CompuTemp2,3,1,CompuTemp3); // QxAIx

		MatAddAB(CompuTemp,CompuTemp3,CompuTemp2,3); // CompuTemp2 = Sinv(D)Wx+QxAIx

		// get (EYE3-SsStk[i+1]*invDkStk[i]*BinvMBT_Stk[i])*vs_Stk[i+1]
		MatMulAB(BinvMBT_Stk[index_now].data,3,3,vsStkY[index_now+1-index_now].data,3,1,CompuTemp); // BMBv
		MatMulAB(invDkStk[index_now-index_now].data,3,3,CompuTemp,3,1,CompuTemp3); // inv(D)BMBv
		MatMulAB(SsStk[index_now+1-index_now].data,3,3,CompuTemp3,3,1,CompuTemp); // Sinv(D)BMBv

		MatMiuAB(vsStkY[index_now+1-index_now].data,CompuTemp,CompuTemp,3); // v-Sinv(D)BMBv
		MatMiuAB(CompuTemp,CompuTemp2,CompuTemp,3); // (v-Sinv(D)BMBv) - (Sinv(D)Wx+QxAIx)
		
		MatMulAB(BinvMBT_Stk[index_now].data,3,3,CompuTemp,3,1,Bkuk.data); // Bkuk found

		//MatMulAB(AdStk[index_now].data,3,3,YState[index_now].data,3,1,CompuTemp);//angel ���ƤF?
		//MatAddAB(CompuTemp,Bkuk.data(),XState[i+1].data(),3);
		//MatAddAB(CompuTemp,Bkuk.data,CompuTemp,3);//angel ���ƤF?

		// xk+1 = Axk+Bkuk
		MatMulAB(AdStk[index_now].data,3,3,YState2[index_now].data,3,1,CompuTemp);
		MatAddAB(CompuTemp,Bkuk.data,YState[index_now+1].data,3);

	//}

	//#if SaveLQSIStates
	//	// �nprint���G�X�Ӫ��ܡA�Х��}�U��
		/*if (index_now==7215)
		{
		
		fstream Fx;

		Fx.open("statesx_ideal.txt",ios::out);
		Fx.precision(10);
		for (int i=0 ; i<= LQDataLen ; i++)
		{
			Fx << XState[i].data[0] << " " << XState[i].data[1] << " " << XState[i].data[2] <<  endl;
		}

		Fx.close();
		Fx.open("statesx_adams.txt",ios::out);
		Fx.precision(10);
		for (int i=0 ; i<= LQDataLen ; i++)
		{
			Fx << XState2[i].data[0] << " " << XState2[i].data[1] << " " << XState2[i].data[2] <<  endl;
		}
		Fx.close();

		fstream Fy;

		Fy.open("statesy_ideal.txt",ios::out);
		Fy.precision(10);
		for (int i=0 ; i<= LQDataLen ; i++)
		{
			Fy << YState[i].data[0] << " " << YState[i].data[1] << " " << YState[i].data[2] <<  endl;
		}
		Fy.close();
		Fy.open("statesy_adams.txt",ios::out);
		Fy.precision(10);
		for (int i=0 ; i<= LQDataLen ; i++)
		{
			Fy << YState2[i].data[0] << " " << YState2[i].data[1] << " " << YState2[i].data[2] <<  endl;
		}
		}*/
		//******test for Sk
		/*if (index_now==2000)
		{
			fstream Sk;
			Sk.open("Sk.txt",ios::out);
			Sk.precision(10);
			for (int i=index_now ; i<= index_now+1200 ; i++)
			{
				Sk << SsStk[i].data[0] << " " << SsStk[i].data[1] << " " << SsStk[i].data[2] << " " << SsStk[i].data[3] <<" " << SsStk[i].data[4] <<" " << SsStk[i].data[5] <<" " << SsStk[i].data[6] <<" " << SsStk[i].data[7] <<" " << SsStk[i].data[8] << endl;
			}
			Sk.close();
		}*/
		
		
	//	}// �nprint���G�X�Ӫ��ܡA�Х��}�W��
	//#endif
}

void LQSISolver::C2DmHumanoid(int i)
{
	/******************************************************************
	input: i
	output: void

	Note:
	// ��Ji ���禡�N�|�p���i�� sampling time�� state-space �Azero order hold ���᪺��
	******************************************************************/
	//extern Matrix<double,3,3>* AdStk;

	Wp = sqrt((ddelCOG[i]+GravityConst)/(COGz[i+4])); // �o�䪺wp = sqrt(wp) in paper
	sh = sinh(Wp*dt);
	ch = cosh(Wp*dt);

	AdStk[i].data[0] = ch;
	AdStk[i].data[1] = sh/Wp;
	AdStk[i].data[2] = 1-ch;
	AdStk[i].data[3] = Wp*sh;
	AdStk[i].data[4] = ch;
	AdStk[i].data[5] = -Wp*sh;

	BdStk[i].data[0] = dt-(sh/Wp);
	BdStk[i].data[1] = -ch+1;

	// sample matlab code
	////// From Matlab code LQSI_Planner2.m ///////////
    //sw = sqrt(w);
    //sh = sinh(sw*T);
    //ch = cosh(sw*T);
    //Ad = [ch sh/sw 1-ch ; sw*sh ch -sw*sh ; 0 0 1];
    //Bd = [T-sh/sw ; -ch+1 ; T];
    //Cd = C;
    //Dd = D;
	////// From Matlab code LQSI_Planner2.m ///////////

}
 
void LQSISolver::tic(void)
{
	/******************************************************************
	input: void
	output: void

	Note:
	// �p�PMATLAB�� tic(); function �p�ɶ}�l
	******************************************************************/
	QueryPerformanceCounter(&gStartTime);
}

void LQSISolver::ticPMS(void)
{
	/******************************************************************
	input: void
	output: void

	Note:
	// �p�PMATLAB�� tic(); function �p�ɶ}�l
	******************************************************************/
	QueryPerformanceCounter(&gStartTime_PMS);
}


void LQSISolver::toc(void)
{
	/******************************************************************
	input: void
	output: void

	Note:
	// �p�PMATLAB�� toc(); function �p�ɵ��� �åB�bcommand window �L�X�g�L�ɶ�
	******************************************************************/
	QueryPerformanceCounter(&gCurrentTime);	
	double gSysTime = (gCurrentTime.QuadPart - gStartTime.QuadPart)/gFreqT;
	//cout << gSysTime << endl;
			TimerRecord.open("timerdata.txt",ios::app);
			TimerRecord << gSysTime << "\t";
			TimerRecord.close();
}

void LQSISolver::toc2(void)
{
	/******************************************************************
	input: void
	output: void

	Note:
	// ���MATLAB�� toc(); function �p�ɵ��� �|�bcommand window �L�X�̤j�B��ɶ�
	******************************************************************/
	QueryPerformanceCounter(&gCurrentTime);	
	double gSysTime = (gCurrentTime.QuadPart - gStartTime.QuadPart)/gFreqT;
	if (gSysTime>MaxIterationTime)
	MaxIterationTime=gSysTime;
	cout << MaxIterationTime << endl;
		// save encoder
			TimerRecord.open("timerdata.txt",ios::app);
			TimerRecord << MaxIterationTime << "\t";
			TimerRecord.close();
		// save encoder
}

void LQSISolver::toc3(void)
{
	/******************************************************************
	input: void
	output: void

	Note:
	// �p�PMATLAB�� toc(); function �p�ɵ��� �åB�bcommand window �L�X�g�L�ɶ�
	******************************************************************/
	QueryPerformanceCounter(&gCurrentTime);	
	double gSysTime = (gCurrentTime.QuadPart - gStartTime.QuadPart)/gFreqT;

	LogSysTime[TimeCount]=gSysTime;
	TimeCount++;
}

void LQSISolver::toc3PMS(void)
{
	/******************************************************************
	input: void
	output: void

	Note:
	// �p�PMATLAB�� toc(); function �p�ɵ��� �åB�bcommand window �L�X�g�L�ɶ�
	******************************************************************/
	QueryPerformanceCounter(&gCurrentTime_PMS);	
	double gSysTime = (gCurrentTime_PMS.QuadPart - gStartTime_PMS.QuadPart)/gFreqT_PMS;

	LogPMSTime[TimeCountPMS]=gSysTime;
	TimeCountPMS++;
}

void DiffEqn(double* A, int* LenA, double* result, double* d_t)
{
	/******************************************************************
	input:  A ��J�n�Q�L�����}�C�A LenA �}�CA�����סA ���G�s�b result�A d_t�O�L����delta t 
	output: void

	Note:
	// ���L����k�|�ϱo�̥��̥k��Ӥ����Q�˱�A�ҥH�}�C���׷|��4
	******************************************************************/
	double temp = 1.0/12.0/(*d_t);
	for (int i=2 ; i< (*LenA)-2; i++)
	{
		result[i-2] = temp*(A[i-2]-8*A[i-1]+8*A[i+1]-A[i+2]);
	}
}

void DiffEqnShift(double* A, int* LenA, double* result, double* d_t)
{
	/******************************************************************
	input:  A ��J�n�Q�L�����}�C�A LenA �}�CA�����סA ���G�s�b result�A d_t�O�L����delta t 
	output: void

	Note:
	// ���L����k�|�ϱo�̥��̥k��Ӥ����Q�˱�A�ҥH�}�C���׷|��4 
	// ���O�ѩ�b�{�����N�s���Y������ƻs�⦸ �ҥH�i�H���������
	******************************************************************/
	double temp = 1.0/12.0/(*d_t);
	for (int i=2 ; i< (*LenA)-2; i++)
	{
		result[i] = temp*(A[i-2]-8*A[i-1]+8*A[i+1]-A[i+2]);
	}
	result[0] = result[2];
	result[1] = result[2];
	result[*LenA-1] = result[*LenA-3];
	result[*LenA-2] = result[*LenA-3];

}

void LowPassFilter(double* data, int dataLen, double* result, double gain)
{
	/******************************************************************
	input:  data ��J�n�Qfilter���}�C�A dataLen �}�CA�����סA ���G�s�b result�A gain�Ofilter��filter gain 
	output: void

	Note:
	// �Q�ήt����{���F���Ʀ��o�i��
	******************************************************************/
	result[0] = data[0];

	for (int i=1;i<dataLen;i++)
	{
		result[i] = 1.0/(gain+1.0)*(data[i]+data[i-1])+(gain-1.0)/(gain+1.0)*result[i-1];
	}
	// Matlab sample codes

	//L = length(data);
	//result = zeros(L,1);
	//result(1) = data(1);
	//for i = 2:L
	//    result(i) = 1/(k+1)*(data(i)+data(i-1))+(k-1)/(k+1)*result(i-1);
	//end
}

void GetMTSAcceleration(double* sig, double* h, double* q, int InputTrajLen, double* Boundary, double* ddQ)
{
	/******************************************************************
	input: sig,  h,  q,  InputTrajLen,  Boundary,  ddQ �����p�U
	output: void

	Note:
	// �p��MTS��k���y�񤺴�
	// sig = sigma
	// h = time intervals �w�g���J��q�۰ʥͦ��I���ɶ��q�}�C ����=TrajLen+1�q
	// q = the traj knots ��y��A�|�����J���I
	// boundary are the boundry conditions boundary = [v_1 v_end a_1 a_end]
	// TrajLen = ��l�y���
	// ddQ = double derivative of Q
	******************************************************************/

	// �ЫؼȦs��
	double* a;
	a = new double[InputTrajLen+1];
	double* b;
	b = new double[InputTrajLen+1];

	double* A;
	A = new double[InputTrajLen*InputTrajLen];
	double* B;
	B = new double[InputTrajLen];

	double temp1, temp2;
	double q_hat_2, q_hat_nm1; // q2 and qn-1

	// get all a and b
	for (int i = 0 ; i < InputTrajLen ; i++)
	{
		temp1 = 1.0/sig[i]/sig[i]/h[i];
		temp2 = sig[i]*h[i];
		a[i] = temp1-1.0/sig[i]/sinh(temp2);
		b[i] = 1.0/sig[i]/tanh(temp2)-temp1;
	}


	// Find A
		// clear A
	for (int i = 0 ; i < InputTrajLen*InputTrajLen ; i++)
	{
		A[i] = 0;
	}

	A[0] = b[0]+b[1]+h[0]*h[0]/6.0/h[1]+h[0]/6.0; // first one
	A[InputTrajLen*InputTrajLen-1] = b[InputTrajLen-2]+b[InputTrajLen-1]+h[InputTrajLen-1]*h[InputTrajLen-1]/6.0/h[InputTrajLen-2]+h[InputTrajLen-1]/6.0; // final one

	for ( int i = 1 ; i < InputTrajLen ; i++)
	{
		A[InputTrajLen*(i-1)+i] = a[i];
		A[InputTrajLen*i+i-1] = a[i];
	}
	for ( int i = 1 ; i < InputTrajLen-1 ; i++)
	{
		A[InputTrajLen*i+i] = b[i]+b[i+1];
	}
	// ��K�g�{�� �̫�b�ɤW�t��
	A[InputTrajLen] -= h[0]*h[0]/6.0/h[1];
	A[InputTrajLen*(InputTrajLen-1)-1] -= h[InputTrajLen-2]*h[InputTrajLen-2]/6.0/h[InputTrajLen-3];


	// Find B
		// q2 and qn-1
		// boundary = [v_1 v_end a_1 a_end]
	q_hat_2 = q[0]+h[0]*Boundary[0]+h[0]*h[0]*Boundary[2]/3.0;
	q_hat_nm1 = q[InputTrajLen-1]-h[InputTrajLen-1]*Boundary[1]+h[InputTrajLen-1]*h[InputTrajLen-1]*Boundary[3]/3.0; 


	B[0] = (q[1]-q_hat_2)/h[1]-(q_hat_2-q[0])/h[0]-a[0]*Boundary[2];
	B[1] = (q[2]-q[1])/h[2]-(q[1]-q_hat_2)/h[1];
	B[InputTrajLen-1] = (q[InputTrajLen-1]-q_hat_nm1)/h[InputTrajLen-1]-(q_hat_nm1-q[InputTrajLen-2])/h[InputTrajLen-2]-a[InputTrajLen-1]*Boundary[3];
	B[InputTrajLen-2] = (q_hat_nm1-q[InputTrajLen-2])/h[InputTrajLen-2]-(q[InputTrajLen-2]-q[InputTrajLen-3])/h[InputTrajLen-3];

	for ( int i = 2 ; i < InputTrajLen-2 ; i++)
	{
		B[i] = (q[i+1]-q[i])/h[i+1]-(q[i]-q[i-1])/h[i];
	}

	InvSqMat(A,InputTrajLen);

	// �p��[�t��
	MatMulAB(A,InputTrajLen,InputTrajLen,B,InputTrajLen,1,ddQ);


	for (int pp = 0 ; pp<InputTrajLen ; pp++)
	{
		printf(" %f\n",ddQ[pp]);
	}

	// �M���ʺA�O����
	delete[] a;
	delete[] b;
	delete[] A;
	delete[] B;

}

void MovingAvergeSmooth(double* Traj, int LengthOfTraj, int FrameSize, double* result)
{
	/******************************************************************
	input: Traj ��l��J�y��,  LengthOfTraj ��l�y�����,  FrameSize �n�������]�t�y��d��,  ���G�s�bresult
	output: void

	Note:
	// �Q��Moving Average��k�D�o������
	// ����k���u�I�O�S��phase lag �i�H�������u 
	// �ʽ��low pass filter ���I���P

	******************************************************************/

	int a,b;
	int Size;

	if (FrameSize%2 == 0)
		Size = FrameSize-1;
	else
		Size = FrameSize;

	a = (Size-1)/2;

	for (int i = 0 ; i < LengthOfTraj ; i++)
	{
		result[i] = 0;
		b = LengthOfTraj-1-i;

		if (i < a)
		{
			for (int j = 0 ; j <= 2*i ; j++)
			{
				result[i] += Traj[j];
			}
			result[i] /= (i*2+1);
		}
		else if (i+a > LengthOfTraj-1)
		{
			for (int j = i-b ; j <= LengthOfTraj-1 ; j++)
			{
				result[i] += Traj[j];
			}
			result[i] /= (b*2+1);
		}
		else
		{
			for (int j = i-a ; j <= i+a ; j++)
			{
				result[i] += Traj[j];
			}
			result[i] /= (Size);
		}
	}
}

void LQSISolver::ZMPFeedbackControl(double* ZMPx, double* ZMPy,double* COGx ,double* COGy , int index)
{
	/******************************************************************
	input: ZMPx, ZMPy, COGx, COGy, �}�C����(�W���w�w)
	output: void
	Note:
	// ��DummyControl ���ܦӨ�, 
	   �N�Ҧ��w�g��n��state matrices and control laws �Ӻ�X�������߭y��
	// �MDummyControl���P�B��,feedback���A�ϥ��|�N�����G 
	   �ӥ�ManinLoops.cpp  C2MWrite2Txt() ����Case 5 ��update�s��ZMP��m
	   �]���ݨC��iteration(5ms)�N�u����ZMP��^state space����������m
	   �M��ΦA�M�w����n��control laws �Ӻ�X�C��5ms��COG���������߭y��
	******************************************************************/
//	double CompuTemp[400]; // �x�}���k���ݭn���ƧQ�Ϊ�function�ݭn��ӼȦs�� ���M�|�мg�쥿�b�Ϊ����
//	double CompuTemp2[400];
//	double CompuTemp3[400];
//	double CompuTemp4[400];

	YMatLite EYE3(3,3);
	EYE3.data[0] = 1;
	EYE3.data[1] = 0;
	EYE3.data[2] = 0;
	EYE3.data[3] = 0;
	EYE3.data[4] = 1;
	EYE3.data[5] = 0;
	EYE3.data[6] = 0;
	EYE3.data[7] = 0;
	EYE3.data[8] = 1;

	//XState[0].data[0] = 0;
	//XState[0].data[1] = 0;
	//XState[0].data[2] = 0;

	//Matrix<double,3,1> Bkuk;

	YMatLite Bkuk(3,1);
	double Aka=0;
	//tic();

	//for (int i = 0 ; i<LQDataLen ; i++)
	//{
		// Bkuk = BinvMBStk(:,:,i)*(-(Qx*AI+Ss(:,:,i+1)*invDkStk(:,:,i)*WkStk(:,:,i))*x+(eye(3)-Ss(:,:,i+1)*invDkStk(:,:,i)*BinvMBStk(:,:,i))*vs(:,i+1)); 	
		// x = AkStk(:,:,i)*x + Bkuk; 

		//Bkuk = BinvMBT_Stk[i]*(-(Qx*AIStk[i]+SsStk[i+1]*invDkStk[i]*WkStk[i])*XState[i]+(EYE3-SsStk[i+1]*invDkStk[i]*BinvMBT_Stk[i])*vs_Stk[i+1]);
		//XState[i+1] = AdStk[i]*XState[i] + Bkuk;


		MatMulAB(WkStk[index].data,3,3,XState[index].data,3,1,CompuTemp); // Wx
		MatMulAB(invDkStk[index].data,3,3,CompuTemp,3,1,CompuTemp2); // inv(D)Wx
		MatMulAB(SsStk[index+1].data,3,3,CompuTemp2,3,1,CompuTemp); // Sinv(D)Wx

		MatMulAB(AIStk[index].data,3,3,XState[index].data,3,1,CompuTemp2); // AIx
		MatMulAB(Qx->data,3,3,CompuTemp2,3,1,CompuTemp3); // QxAIx

		MatAddAB(CompuTemp,CompuTemp3,CompuTemp2,3); // CompuTemp2 = Sinv(D)Wx+QxAIx

		// get (EYE3-SsStk[i+1]*invDkStk[i]*BinvMBT_Stk[i])*vs_Stk[i+1]
		MatMulAB(BinvMBT_Stk[index].data,3,3,vsStkX[index+1].data,3,1,CompuTemp); // BMBv
		MatMulAB(invDkStk[index].data,3,3,CompuTemp,3,1,CompuTemp3); // inv(D)BMBv
		MatMulAB(SsStk[index+1].data,3,3,CompuTemp3,3,1,CompuTemp); // Sinv(D)BMBv

		MatMiuAB(vsStkX[index+1].data,CompuTemp,CompuTemp,3); // v-Sinv(D)BMBv
		MatMiuAB(CompuTemp,CompuTemp2,CompuTemp,3); // (v-Sinv(D)BMBv) - (Sinv(D)Wx+QxAIx)
		
		MatMulAB(BinvMBT_Stk[index].data,3,3,CompuTemp,3,1,Bkuk.data); // Bkuk found

		MatMulAB(AdStk[index].data,3,3,XState[index].data,3,1,CompuTemp);
		//MatAddAB(CompuTemp,Bkuk.data(),XState[i+1].data(),3);
		MatAddAB(CompuTemp,Bkuk.data,CompuTemp,3);

		// xk+1 = Axk+Bkuk
		MatMulAB(AdStk[index].data,3,3,XState[index].data,3,1,CompuTemp);
		MatAddAB(CompuTemp,Bkuk.data,XState[index+1].data,3);


	//}

		double QQ=0;
	//for (int i = 0 ; i<LQDataLen ; i++)
	//{
		// Bkuk = BinvMBStk(:,:,i)*(-(Qx*AI+Ss(:,:,i+1)*invDkStk(:,:,i)*WkStk(:,:,i))*x+(eye(3)-Ss(:,:,i+1)*invDkStk(:,:,i)*BinvMBStk(:,:,i))*vs(:,i+1)); 	
		// x = AkStk(:,:,i)*x + Bkuk; 

		//Bkuk = BinvMBT_Stk[i]*(-(Qx*AIStk[i]+SsStk[i+1]*invDkStk[i]*WkStk[i])*XState[i]+(EYE3-SsStk[i+1]*invDkStk[i]*BinvMBT_Stk[i])*vs_Stk[i+1]);
		//XState[i+1] = AdStk[i]*XState[i] + Bkuk;


		MatMulAB(WkStk[index].data,3,3,YState[index].data,3,1,CompuTemp); // Wx
		MatMulAB(invDkStk[index].data,3,3,CompuTemp,3,1,CompuTemp2); // inv(D)Wx
		MatMulAB(SsStk[index+1].data,3,3,CompuTemp2,3,1,CompuTemp); // Sinv(D)Wx

		MatMulAB(AIStk[index].data,3,3,YState[index].data,3,1,CompuTemp2); // AIx
		MatMulAB(Qx->data,3,3,CompuTemp2,3,1,CompuTemp3); // QxAIx

		MatAddAB(CompuTemp,CompuTemp3,CompuTemp2,3); // CompuTemp2 = Sinv(D)Wx+QxAIx

		// get (EYE3-SsStk[i+1]*invDkStk[i]*BinvMBT_Stk[i])*vs_Stk[i+1]
		MatMulAB(BinvMBT_Stk[index].data,3,3,vsStkY[index+1].data,3,1,CompuTemp); // BMBv
		MatMulAB(invDkStk[index].data,3,3,CompuTemp,3,1,CompuTemp3); // inv(D)BMBv
		MatMulAB(SsStk[index+1].data,3,3,CompuTemp3,3,1,CompuTemp); // Sinv(D)BMBv

		MatMiuAB(vsStkY[index+1].data,CompuTemp,CompuTemp,3); // v-Sinv(D)BMBv
		MatMiuAB(CompuTemp,CompuTemp2,CompuTemp,3); // (v-Sinv(D)BMBv) - (Sinv(D)Wx+QxAIx)
		
		MatMulAB(BinvMBT_Stk[index].data,3,3,CompuTemp,3,1,Bkuk.data); // Bkuk found

		MatMulAB(AdStk[index].data,3,3,YState[index].data,3,1,CompuTemp);
		//MatAddAB(CompuTemp,Bkuk.data(),XState[i+1].data(),3);
		MatAddAB(CompuTemp,Bkuk.data,CompuTemp,3);
		
		// xk+1 = Axk+Bkuk
		MatMulAB(AdStk[index].data,3,3,YState[index].data,3,1,CompuTemp);
		MatAddAB(CompuTemp,Bkuk.data,YState[index+1].data,3);
		
	//}
		QQ=AdStk[index].data[0];
		QQ=AdStk[index].data[1];
		QQ=AdStk[index].data[2];
		QQ=AdStk[index].data[3];
		QQ=AdStk[index].data[4];
		QQ=AdStk[index].data[5];
		QQ=AdStk[index].data[6];
		QQ=AdStk[index].data[7];
		QQ=AdStk[index].data[8];

	//// �nprint�X�Ӫ��ܡA�Х��}�U��

	//fstream Fx;

	//Fx.open("statesx.txt",ios::out);
	//Fx.precision(10);
	//for (int i=0 ; i<= LQDataLen ; i++)
	//{
	//	Fx << XState[i].data[0] << " " << XState[i].data[1] << " " << XState[i].data[2] <<  endl;
	//}

	//Fx.close();

	//fstream Fy;

	//Fy.open("statesy.txt",ios::out);
	//Fy.precision(10);
	//for (int i=0 ; i<= LQDataLen ; i++)
	//{
	//	Fy << YState[i].data[0] << " " << YState[i].data[1] << " " << YState[i].data[2] <<  endl;
	//}

	//Fy.close();

	//// �nprint�X�Ӫ��ܡA�Х��}�W��
}