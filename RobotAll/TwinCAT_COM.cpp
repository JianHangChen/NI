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
#include "StdAfx.h"
#include "MainLoops.h"

#if TwinCAT_Mode        // Defined in Mainlop.h
#include "TwinCAT_COM.h"



TwinCAT_COM::TwinCAT_COM(void)
{
	/******************************************************************
	input: void
	output: void

	Note:
	// Class constructor  初始化+兩腳Homing+Set SCP Mode
	******************************************************************/

	EtherCATInit();

	cout<<"Press Any Key to Start Right Leg Homing!"<<endl;
	system("pause");
	EtherCATHomingRLeg();
	EtherCATSetCSPRLeg();

	cout<<"Press Any Key to Start Left Leg Homing!"<<endl;
	system("pause");
	EtherCATHomingLLeg();
	EtherCATSetCSPLLeg();

	//Buffer initilization for moving avergae	
	for (int i=0;i<MovingIndex*12;i++)
		MovingBuffer[i]=0;

	//Pointer sequence for numerous variables
	pCtrlWords[0]=& pTLL1msOut->LL_O_CtrlWord_01;
	pCtrlWords[1]=& pTLL1msOut->LL_O_CtrlWord_02;
	pCtrlWords[2]=& pTLL1msOut->LL_O_CtrlWord_03;
	pCtrlWords[3]=& pTLL1msOut->LL_O_CtrlWord_04;
	pCtrlWords[4]=& pTLL1msOut->LL_O_CtrlWord_05;
	pCtrlWords[5]=& pTLL1msOut->LL_O_CtrlWord_06;

	pCtrlWords[6]=& pTRL1msOut->RL_O_CtrlWord_01;
	pCtrlWords[7]=& pTRL1msOut->RL_O_CtrlWord_02;
	pCtrlWords[8]=& pTRL1msOut->RL_O_CtrlWord_03;
	pCtrlWords[9]=& pTRL1msOut->RL_O_CtrlWord_04;
	pCtrlWords[10]=& pTRL1msOut->RL_O_CtrlWord_05;
	pCtrlWords[11]=& pTRL1msOut->RL_O_CtrlWord_06;
}
void TwinCAT_COM::EtherCATInit(void)
{
	/******************************************************************
	input: void
	output: void

	Note:
	// TwinCAT I/O 初始化
	   由12個獨立的while loop組成,按序列式依序進行初始化的動作
	   條件判斷式詳見TwinCAT R3IO範例程式
	******************************************************************/
int chklimit=10; //初始化嘗試連接次數 上限值
	if ( TCatIoOpen() == 0 ) 
	{ 
		// get the process image pointer 
		if ( TCatIoGetOutputPtr(TASK_RLEG_PORTNUMBER, (void**)&pTRL1msOut,sizeof(Task_RLeg_Outputs) ) == 0 && 
				TCatIoGetInputPtr( TASK_RLEG_PORTNUMBER, (void**)&pTRL1msIn, sizeof(Task_RLeg_Inputs) ) == 0 &&
			 TCatIoGetOutputPtr(TASK_LLEG_PORTNUMBER, (void**)&pTLL1msOut,sizeof(Task_LLeg_Outputs) ) == 0 &&
				TCatIoGetInputPtr( TASK_LLEG_PORTNUMBER, (void**)&pTLL1msIn, sizeof(Task_LLeg_Inputs) ) == 0
				) 
		{ 
			int counterchk=0;
			long nError;
			long nError2;
			int target=0;
			//Copy
			while(1)
			{
					pTRL1msOut->RL_O_CtrlWord_01=6;
					if ( (nError = TCatIoInputUpdate( TASK_RLEG_PORTNUMBER )) == 0   ) 
					TCatIoOutputUpdate( TASK_RLEG_PORTNUMBER );
					Sleep(20);
					TCatIoInputUpdate( TASK_RLEG_PORTNUMBER ) ;
					if( DisabledCheck(pTRL1msIn->RL_I_StatWord_01) )
					{
						cout<<"Idle Mode Activated: Device R01 Shoutdown"<<endl;
						cout<<endl;
						break;
					}
					else
					{
						//cout<<"Idle: RL 01 Failed"<<endl;
						counterchk++;
						Sleep(20);
						if(counterchk>chklimit)
						{
							cout<<"R01 EtherCAT Error!"<<endl;
							system("pause");
						}
					}
					//cout<<endl;
			}
			//Copy
			//Copy
			while(1)
			{
					pTRL1msOut->RL_O_CtrlWord_02=6;
					if ( (nError = TCatIoInputUpdate( TASK_RLEG_PORTNUMBER )) == 0) 
					TCatIoOutputUpdate( TASK_RLEG_PORTNUMBER );
					Sleep(20);
					TCatIoInputUpdate( TASK_RLEG_PORTNUMBER ) ;
					if( DisabledCheck(pTRL1msIn->RL_I_StatWord_02) )
					{
						cout<<"Idle Mode Activated: Device R02 Shoutdown"<<endl;
						cout<<endl;
						break;
					}
					else
					{
						//cout<<"Idle: RL 02 Failed"<<endl;
						counterchk++;
						Sleep(20);
						if(counterchk>chklimit)
						{
							cout<<"R02 EtherCAT Error!"<<endl;
							system("pause");
						}
					}
					//cout<<endl;
			}
			//Copy
			//Copy
			while(1)
			{
					pTRL1msOut->RL_O_CtrlWord_03=6;
					if ( (nError = TCatIoInputUpdate( TASK_RLEG_PORTNUMBER )) == 0  ) 
					TCatIoOutputUpdate( TASK_RLEG_PORTNUMBER );
					Sleep(20);
					TCatIoInputUpdate( TASK_RLEG_PORTNUMBER ) ;
					if( DisabledCheck(pTRL1msIn->RL_I_StatWord_03) )
					{
						cout<<"Idle Mode Activated: Device R03 Shoutdown"<<endl;
						cout<<endl;
						break;
					}
					else
					{
						//cout<<"Idle: RL 03 Failed"<<endl;
						counterchk++;
						Sleep(20);
						if(counterchk>chklimit)
						{
							cout<<"R03 EtherCAT Error!"<<endl;
							system("pause");
						}
					}
					//cout<<endl;
			}
			//Copy
			//Copy
			while(1)
			{
					pTRL1msOut->RL_O_CtrlWord_04=6;
					if ( (nError = TCatIoInputUpdate( TASK_RLEG_PORTNUMBER )) == 0  ) 
					TCatIoOutputUpdate( TASK_RLEG_PORTNUMBER );
					Sleep(20);
					TCatIoInputUpdate( TASK_RLEG_PORTNUMBER ) ;
					if( DisabledCheck(pTRL1msIn->RL_I_StatWord_04) )
					{
						cout<<"Idle Mode Activated: Device R04 Shoutdown"<<endl;
						cout<<endl;
						break;
					}
					else
					{
						//cout<<"Idle: RL 04 Failed"<<endl;
						counterchk++;
						Sleep(20);
						if(counterchk>chklimit)
						{
							cout<<"R04 EtherCAT Error!"<<endl;
							system("pause");
						}
					}
					//cout<<endl;
			}
			//Copy
			//Copy
			while(1)
			{
					pTRL1msOut->RL_O_CtrlWord_05=6;
					if ( (nError = TCatIoInputUpdate( TASK_RLEG_PORTNUMBER )) == 0 ) 
					TCatIoOutputUpdate( TASK_RLEG_PORTNUMBER );
					Sleep(20);
					TCatIoInputUpdate( TASK_RLEG_PORTNUMBER ) ;
					if( DisabledCheck(pTRL1msIn->RL_I_StatWord_05) )
					{
						cout<<"Idle Mode Activated: Device R05 Shoutdown"<<endl;
						cout<<endl;
						break;
					}
					else
					{
						//cout<<"Idle: RL 05 Failed"<<endl;
						counterchk++;
						Sleep(20);
						if(counterchk>chklimit)
						{
							cout<<"R05 EtherCAT Error!"<<endl;
							system("pause");
						}
					}
					//cout<<endl;
			}
			//Copy
			//Copy
			while(1)
			{
					pTRL1msOut->RL_O_CtrlWord_06=6;
					if ( (nError = TCatIoInputUpdate( TASK_RLEG_PORTNUMBER )) == 0  ) 
					TCatIoOutputUpdate( TASK_RLEG_PORTNUMBER );
					Sleep(20);
					TCatIoInputUpdate( TASK_RLEG_PORTNUMBER ) ;
					if( DisabledCheck(pTRL1msIn->RL_I_StatWord_06) )
					{
						cout<<"Idle Mode Activated: Device R06 Shoutdown"<<endl;
						cout<<endl;
						break;
					}
					else
					{
						//cout<<"Idle: RL 06 Failed"<<endl;
						counterchk++;
						Sleep(20);
						if(counterchk>chklimit)
						{
							cout<<"R06 EtherCAT Error!"<<endl;
							system("pause");
						}
					}
					//cout<<endl;
			}
			//Copy
					Sleep(100);
			//Copy
			while(1)
			{
					pTLL1msOut->LL_O_CtrlWord_01=6;
					if ( (nError = TCatIoInputUpdate( TASK_LLEG_PORTNUMBER )) == 0    ) 
					TCatIoOutputUpdate( TASK_LLEG_PORTNUMBER );
					Sleep(20);
					TCatIoInputUpdate( TASK_LLEG_PORTNUMBER ) ;
					if( DisabledCheck(pTLL1msIn->LL_I_StatWord_01) )
					{
						cout<<"Idle Mode Activated: Device L01 Shoutdown"<<endl;
						cout<<endl;
						break;
					}
					else
					{
						//cout<<"Idle: RL 01 Failed"<<endl;
						counterchk++;
						Sleep(20);
						if(counterchk>chklimit)
						{
							cout<<"L01 EtherCAT Error!"<<endl;
							system("pause");
						}
					}
					//cout<<endl;
			}
			//Copy
			//Copy
			while(1)
			{
					pTLL1msOut->LL_O_CtrlWord_02=6;
					if ( (nError = TCatIoInputUpdate( TASK_LLEG_PORTNUMBER )) == 0 ) 
					TCatIoOutputUpdate( TASK_LLEG_PORTNUMBER );
					Sleep(20);
					TCatIoInputUpdate( TASK_LLEG_PORTNUMBER ) ;
					if( DisabledCheck(pTLL1msIn->LL_I_StatWord_02) )
					{
						cout<<"Idle Mode Activated: Device L02 Shoutdown"<<endl;
						cout<<endl;
						break;
					}
					else
					{
						//cout<<"Idle: LL 02 Failed"<<endl;
						counterchk++;
						Sleep(20);
						if(counterchk>chklimit)
						{
							cout<<"L02 EtherCAT Error!"<<endl;
							system("pause");
						}
					}
					//cout<<endl;
			}
			//Copy
			//Copy
			while(1)
			{
					pTLL1msOut->LL_O_CtrlWord_03=6;
					if ( (nError = TCatIoInputUpdate( TASK_LLEG_PORTNUMBER )) == 0 ) 
					TCatIoOutputUpdate( TASK_LLEG_PORTNUMBER );
					Sleep(20);
					TCatIoInputUpdate( TASK_LLEG_PORTNUMBER ) ;
					if( DisabledCheck(pTLL1msIn->LL_I_StatWord_03) )
					{
						cout<<"Idle Mode Activated: Device L03 Shoutdown"<<endl;
						cout<<endl;
						break;
					}
					else
					{
						//cout<<"Idle: LL 03 Failed"<<endl;
						counterchk++;
						Sleep(20);
						if(counterchk>chklimit)
						{
							cout<<"L03 EtherCAT Error!"<<endl;
							system("pause");
						}
					}
					//cout<<endl;
			}
			//Copy
			//Copy
			while(1)
			{
					pTLL1msOut->LL_O_CtrlWord_04=6;
					if ( (nError = TCatIoInputUpdate( TASK_LLEG_PORTNUMBER )) == 0 ) 
					TCatIoOutputUpdate( TASK_LLEG_PORTNUMBER );
					Sleep(20);
					TCatIoInputUpdate( TASK_LLEG_PORTNUMBER ) ;
					if( DisabledCheck(pTLL1msIn->LL_I_StatWord_04) )
					{
						cout<<"Idle Mode Activated: Device L04 Shoutdown"<<endl;
						cout<<endl;
						break;
					}
					else
					{
						//cout<<"Idle: LL 04 Failed"<<endl;
						counterchk++;
						Sleep(20);
						if(counterchk>chklimit)
						{
							cout<<"L04 EtherCAT Error!"<<endl;
							system("pause");
						}
					}
					//cout<<endl;
			}
			//Copy
			//Copy
			while(1)
			{
					pTLL1msOut->LL_O_CtrlWord_05=6;
					if ( (nError = TCatIoInputUpdate( TASK_LLEG_PORTNUMBER )) == 0  ) 
					TCatIoOutputUpdate( TASK_LLEG_PORTNUMBER );
					Sleep(20);
					TCatIoInputUpdate( TASK_LLEG_PORTNUMBER ) ;
					if( DisabledCheck(pTLL1msIn->LL_I_StatWord_05) )
					{
						cout<<"Idle Mode Activated: Device L05 Shoutdown"<<endl;
						cout<<endl;
						break;
					}
					else
					{
						//cout<<"Idle: LL 05 Failed"<<endl;
						counterchk++;
						Sleep(20);
						if(counterchk>chklimit)
						{
							cout<<"L05 EtherCAT Error!"<<endl;
							system("pause");
						}
					}
					//cout<<endl;
			}
			//Copy
			//Copy
			while(1)
			{
					pTLL1msOut->LL_O_CtrlWord_06=6;
					if ( (nError = TCatIoInputUpdate( TASK_LLEG_PORTNUMBER )) == 0 ) 
					TCatIoOutputUpdate( TASK_LLEG_PORTNUMBER );
					Sleep(20);
					TCatIoInputUpdate( TASK_LLEG_PORTNUMBER ) ;
					if( DisabledCheck(pTLL1msIn->LL_I_StatWord_06) )
					{
						cout<<"Idle Mode Activated: Device L06 Shoutdown"<<endl;
						cout<<endl;
						break;
					}
					else
					{
						//cout<<"Idle: LL 06 Failed"<<endl;
						counterchk++;
						Sleep(20);
						if(counterchk>chklimit)
						{
							cout<<"L06 EtherCAT Error!"<<endl;
							system("pause");
						}
					}
					//cout<<endl;
			}
			//Copy

		} //If  TACTPARAMETER
	}//If   TACTOPEN
	else
	{
		cout<<"Warning!TwinCAT is not opened"<<endl;
		system("pause");
	}
}	
void TwinCAT_COM::EtherCATHomingRLeg(void)
{
	/******************************************************************
	input: void
	output: void

	Note:
	// 透過TwinCAT I/O 進行Homing動作
	   組成方式為六組依序(由腳踝開始)的homing function
	   //Copy
				func ....
	   //Copy
	   請先確認數位電有開啟
	   條件判斷式詳見EPOS3-Firmware Specification & Applicaion Note
	******************************************************************/

	//Copy
		pTRL1msOut->RL_O_CtrlWord_06=6;

 		TCatIoOutputUpdate( TASK_RLEG_PORTNUMBER );
		Sleep(50);
		//change mode
		pTRL1msOut->RL_O_SetMode_06=6;
		TCatIoOutputUpdate( TASK_RLEG_PORTNUMBER );
		Sleep(50);
		cout<<"06 HomingMode Set"<<endl;
		//switch on
		pTRL1msOut->RL_O_CtrlWord_06=15;
		TCatIoOutputUpdate( TASK_RLEG_PORTNUMBER );
		Sleep(50);
		//start homing
		//cout<<"Press Anykey to Start 06 Homing"<<endl;
		cout<<"Start 06 Homing ..."<<endl;
		Sleep(50);
		//system("pause");	
		pTRL1msOut->RL_O_CtrlWord_06=31;
		TCatIoOutputUpdate( TASK_RLEG_PORTNUMBER );
		Sleep(50);

		while(1)
		{
			TCatIoInputUpdate( TASK_RLEG_PORTNUMBER );
			if((pTRL1msIn->RL_I_StatWord_06 >>12)%2==1)
			{
				cout<<"06 Homing Posi Attained"<<endl;			
				break;
			}
		}

		//switch on
		pTRL1msOut->RL_O_CtrlWord_06=15;
		TCatIoOutputUpdate( TASK_RLEG_PORTNUMBER );
		Sleep(50);
		//system("pause");
	//Copy
	//Copy
		pTRL1msOut->RL_O_CtrlWord_05=6;
		TCatIoOutputUpdate( TASK_RLEG_PORTNUMBER );
		Sleep(50);
		//change mode
		pTRL1msOut->RL_O_SetMode_05=6;
		TCatIoOutputUpdate( TASK_RLEG_PORTNUMBER );
		Sleep(50);
		cout<<"05 HomingMode Set"<<endl;
		//switch on
		pTRL1msOut->RL_O_CtrlWord_05=15;
		TCatIoOutputUpdate( TASK_RLEG_PORTNUMBER );
		Sleep(50);
		//start homing
		//cout<<"Press Anykey to Start 04 Homing"<<endl;
		cout<<"Start 05 Homing ..."<<endl;
		Sleep(50);
		//system("pause");	
		pTRL1msOut->RL_O_CtrlWord_05=31;
		TCatIoOutputUpdate( TASK_RLEG_PORTNUMBER );
		//Sleep(50);

		while(1)
		{
			TCatIoInputUpdate( TASK_RLEG_PORTNUMBER );
			if((pTRL1msIn->RL_I_StatWord_05 >>12)%2==1)
			{
				cout<<"05 Homing Posi Attained"<<endl;			
				break;
			}
		}

		//switch on
		pTRL1msOut->RL_O_CtrlWord_05=15;
		TCatIoOutputUpdate( TASK_RLEG_PORTNUMBER );
		Sleep(50);
		//system("pause");
	//Copy
	//Copy
		pTRL1msOut->RL_O_CtrlWord_04=6;
		TCatIoOutputUpdate( TASK_RLEG_PORTNUMBER );
		Sleep(50);
		//change mode
		pTRL1msOut->RL_O_SetMode_04=6;
		TCatIoOutputUpdate( TASK_RLEG_PORTNUMBER );
		Sleep(50);
		cout<<"04 HomingMode Set"<<endl;
		//switch on
		pTRL1msOut->RL_O_CtrlWord_04=15;
		TCatIoOutputUpdate( TASK_RLEG_PORTNUMBER );
		Sleep(50);
		//start homing
		//cout<<"Press Anykey to Start 04 Homing"<<endl;
		cout<<"Start 04 Homing ..."<<endl;
		Sleep(50);
		//system("pause");	
		pTRL1msOut->RL_O_CtrlWord_04=31;
		TCatIoOutputUpdate( TASK_RLEG_PORTNUMBER );
		Sleep(50);

		while(1)
		{
			TCatIoInputUpdate( TASK_RLEG_PORTNUMBER );
			if((pTRL1msIn->RL_I_StatWord_04 >>12)%2==1)
			{
				cout<<"04 Homing Posi Attained"<<endl;			
				break;
			}
		}

		//switch on
		pTRL1msOut->RL_O_CtrlWord_04=15;
		TCatIoOutputUpdate( TASK_RLEG_PORTNUMBER );
		Sleep(50);
		//system("pause");
	//Copy
	//Copy
		pTRL1msOut->RL_O_CtrlWord_03=6;
		TCatIoOutputUpdate( TASK_RLEG_PORTNUMBER );
		Sleep(50);
		//change mode
		pTRL1msOut->RL_O_SetMode_03=6;
		TCatIoOutputUpdate( TASK_RLEG_PORTNUMBER );
		Sleep(50);
		cout<<"03 HomingMode Set"<<endl;
		//switch on
		pTRL1msOut->RL_O_CtrlWord_03=15;
		TCatIoOutputUpdate( TASK_RLEG_PORTNUMBER );
		Sleep(50);
		//start homing
		//cout<<"Press Anykey to Start 03 Homing"<<endl;
		cout<<"Start 03 Homing ..."<<endl;
		Sleep(50);
		//system("pause");	
		pTRL1msOut->RL_O_CtrlWord_03=31;
		TCatIoOutputUpdate( TASK_RLEG_PORTNUMBER );
		Sleep(50);

		while(1)
		{
			TCatIoInputUpdate( TASK_RLEG_PORTNUMBER );
			if((pTRL1msIn->RL_I_StatWord_03 >>12)%2==1)
			{
				cout<<"03 Homing Posi Attained"<<endl;			
				break;
			}
		}

		//switch on
		pTRL1msOut->RL_O_CtrlWord_03=15;
		TCatIoOutputUpdate( TASK_RLEG_PORTNUMBER );
		Sleep(50);
		//system("pause");
	//Copy
	//Copy
		pTRL1msOut->RL_O_CtrlWord_02=6;
		TCatIoOutputUpdate( TASK_RLEG_PORTNUMBER );
		Sleep(50);
		//change mode
		pTRL1msOut->RL_O_SetMode_02=6;
		TCatIoOutputUpdate( TASK_RLEG_PORTNUMBER );
		Sleep(50);
		cout<<"02 HomingMode Set"<<endl;
		//switch on
		pTRL1msOut->RL_O_CtrlWord_02=15;
		TCatIoOutputUpdate( TASK_RLEG_PORTNUMBER );
		Sleep(50);
		//start homing
		//cout<<"Press Anykey to Start 02 Homing"<<endl;
		cout<<"Start 02 Homing ..."<<endl;
		Sleep(50);
		//system("pause");	
		pTRL1msOut->RL_O_CtrlWord_02=31;
		TCatIoOutputUpdate( TASK_RLEG_PORTNUMBER );
		Sleep(50);

		while(1)
		{
			TCatIoInputUpdate( TASK_RLEG_PORTNUMBER );
			if((pTRL1msIn->RL_I_StatWord_02 >>12)%2==1)
			{
				cout<<"02 Homing Posi Attained"<<endl;			
				break;
			}
		}

		//switch on
		pTRL1msOut->RL_O_CtrlWord_02=15;
		TCatIoOutputUpdate( TASK_RLEG_PORTNUMBER );
		Sleep(50);
		//system("pause");
	//Copy
	//Copy
		pTRL1msOut->RL_O_CtrlWord_01=6;
		TCatIoOutputUpdate( TASK_RLEG_PORTNUMBER );
		Sleep(50);
		//change mode
		pTRL1msOut->RL_O_SetMode_01=6;
		TCatIoOutputUpdate( TASK_RLEG_PORTNUMBER );
		Sleep(50);
		cout<<"01 HomingMode Set"<<endl;
		//switch on
		pTRL1msOut->RL_O_CtrlWord_01=15;
		TCatIoOutputUpdate( TASK_RLEG_PORTNUMBER );
		Sleep(50);
		//start homing
		//cout<<"Press Anykey to Start 01 Homing"<<endl;
		cout<<"Start 01 Homing ..."<<endl;
		Sleep(50);
		//system("pause");	
		pTRL1msOut->RL_O_CtrlWord_01=31;
		TCatIoOutputUpdate( TASK_RLEG_PORTNUMBER );
		Sleep(50);

		while(1)
		{
			TCatIoInputUpdate( TASK_RLEG_PORTNUMBER );
			if((pTRL1msIn->RL_I_StatWord_01 >>12)%2==1)
			{
				cout<<"01 Homing Posi Attained"<<endl;			
				break;
			}
		}

		//switch on
		pTRL1msOut->RL_O_CtrlWord_01=15;
		TCatIoOutputUpdate( TASK_RLEG_PORTNUMBER );
		Sleep(50);
	//Copy

		cout<<"Right Leg Homing Finished"<<endl;
		cout<<endl;

}
void TwinCAT_COM::EtherCATHomingLLeg(void)
{
	/******************************************************************
	input: void
	output: void

	Note:
	// 透過TwinCAT I/O 進行Homing動作
	   組成方式為六組依序(由腳踝開始)的homing function
	   //Copy
				func ....
	   //Copy
	   請先確認數位電有開啟
	   條件判斷式詳見EPOS3-Firmware Specification & Applicaion Note
	******************************************************************/
	//Copy
		pTLL1msOut->LL_O_CtrlWord_06=6;
		TCatIoOutputUpdate( TASK_LLEG_PORTNUMBER );
		Sleep(50);
		//change mode
		pTLL1msOut->LL_O_SetMode_06=6;
		TCatIoOutputUpdate( TASK_LLEG_PORTNUMBER );
		Sleep(50);
		cout<<"06 HomingMode Set"<<endl;
		//switch on
		pTLL1msOut->LL_O_CtrlWord_06=15;
		TCatIoOutputUpdate( TASK_LLEG_PORTNUMBER );
		Sleep(50);
		//start homing
		//cout<<"Press Anykey to Start 06 Homing"<<endl;
		cout<<"Start 06 Homing ..."<<endl;
		Sleep(50);
		//system("pause");	
		pTLL1msOut->LL_O_CtrlWord_06=31;
		TCatIoOutputUpdate( TASK_LLEG_PORTNUMBER );
		Sleep(50);

		while(1)
		{
			TCatIoInputUpdate( TASK_LLEG_PORTNUMBER );
			if((pTLL1msIn->LL_I_StatWord_06 >>12)%2==1)
			{
				cout<<"06 Homing Posi Attained"<<endl;			
				break;
			}
		}

		//switch on
		pTLL1msOut->LL_O_CtrlWord_06=15;
		TCatIoOutputUpdate( TASK_LLEG_PORTNUMBER );
		Sleep(50);
	//Copy
	//Copy
		pTLL1msOut->LL_O_CtrlWord_05=6;
		TCatIoOutputUpdate( TASK_LLEG_PORTNUMBER );
		Sleep(50);
		//change mode
		pTLL1msOut->LL_O_SetMode_05=6;
		TCatIoOutputUpdate( TASK_LLEG_PORTNUMBER );
		Sleep(50);
		cout<<"05 HomingMode Set"<<endl;
		//switch on
		pTLL1msOut->LL_O_CtrlWord_05=15;
		TCatIoOutputUpdate( TASK_LLEG_PORTNUMBER );
		Sleep(50);
		//start homing
		//cout<<"Press Anykey to Start 04 Homing"<<endl;
		cout<<"Start 05 Homing ..."<<endl;
		Sleep(50);
		//system("pause");	
		pTLL1msOut->LL_O_CtrlWord_05=31;
		TCatIoOutputUpdate( TASK_LLEG_PORTNUMBER );
		Sleep(50);

		while(1)
		{
			TCatIoInputUpdate( TASK_LLEG_PORTNUMBER );
			if((pTLL1msIn->LL_I_StatWord_05 >>12)%2==1)
			{
				cout<<"05 Homing Posi Attained"<<endl;			
				break;
			}
		}

		//switch on
		pTLL1msOut->LL_O_CtrlWord_05=15;
		TCatIoOutputUpdate( TASK_LLEG_PORTNUMBER );
		Sleep(50);
	//Copy
	//Copy
		pTLL1msOut->LL_O_CtrlWord_04=6;
		TCatIoOutputUpdate( TASK_LLEG_PORTNUMBER );
		Sleep(50);
		//change mode
		pTLL1msOut->LL_O_SetMode_04=6;
		TCatIoOutputUpdate( TASK_LLEG_PORTNUMBER );
		Sleep(50);
		cout<<"04 HomingMode Set"<<endl;
		//switch on
		pTLL1msOut->LL_O_CtrlWord_04=15;
		TCatIoOutputUpdate( TASK_LLEG_PORTNUMBER );
		Sleep(50);
		//start homing
		//cout<<"Press Anykey to Start 04 Homing"<<endl;
		cout<<"Start 04 Homing ..."<<endl;
		Sleep(50);
		//system("pause");	
		pTLL1msOut->LL_O_CtrlWord_04=31;
		TCatIoOutputUpdate( TASK_LLEG_PORTNUMBER );
		Sleep(50);

		while(1)
		{
			TCatIoInputUpdate( TASK_LLEG_PORTNUMBER );
			if((pTLL1msIn->LL_I_StatWord_04 >>12)%2==1)
			{
				cout<<"04 Homing Posi Attained"<<endl;			
				break;
			}
		}

		//switch on
		pTLL1msOut->LL_O_CtrlWord_04=15;
		TCatIoOutputUpdate( TASK_LLEG_PORTNUMBER );
		Sleep(50);
	//Copy
	//Copy
		pTLL1msOut->LL_O_CtrlWord_03=6;
		TCatIoOutputUpdate( TASK_LLEG_PORTNUMBER );
		Sleep(50);
		//change mode
		pTLL1msOut->LL_O_SetMode_03=6;
		TCatIoOutputUpdate( TASK_LLEG_PORTNUMBER );
		Sleep(50);
		cout<<"03 HomingMode Set"<<endl;
		//switch on
		pTLL1msOut->LL_O_CtrlWord_03=15;
		TCatIoOutputUpdate( TASK_LLEG_PORTNUMBER );
		Sleep(50);
		//start homing
		//cout<<"Press Anykey to Start 03 Homing"<<endl;
		cout<<"Start 03 Homing ..."<<endl;
		Sleep(50);
		//system("pause");	
		pTLL1msOut->LL_O_CtrlWord_03=31;
		TCatIoOutputUpdate( TASK_LLEG_PORTNUMBER );
		Sleep(50);

		while(1)
		{
			TCatIoInputUpdate( TASK_LLEG_PORTNUMBER );
			if((pTLL1msIn->LL_I_StatWord_03 >>12)%2==1)
			{
				cout<<"03 Homing Posi Attained"<<endl;			
				break;
			}
		}

		//switch on
		pTLL1msOut->LL_O_CtrlWord_03=15;
		TCatIoOutputUpdate( TASK_LLEG_PORTNUMBER );
		Sleep(50);
	//Copy
	//Copy
		pTLL1msOut->LL_O_CtrlWord_02=6;
		TCatIoOutputUpdate( TASK_LLEG_PORTNUMBER );
		Sleep(50);
		//change mode
		pTLL1msOut->LL_O_SetMode_02=6;
		TCatIoOutputUpdate( TASK_LLEG_PORTNUMBER );
		Sleep(50);
		cout<<"02 HomingMode Set"<<endl;
		//switch on
		pTLL1msOut->LL_O_CtrlWord_02=15;
		TCatIoOutputUpdate( TASK_LLEG_PORTNUMBER );
		Sleep(50);
		//start homing
		//cout<<"Press Anykey to Start 02 Homing"<<endl;
		cout<<"Start 02 Homing ..."<<endl;
		Sleep(50);
		//system("pause");	
		pTLL1msOut->LL_O_CtrlWord_02=31;
		TCatIoOutputUpdate( TASK_LLEG_PORTNUMBER );
		Sleep(50);

		while(1)
		{
			TCatIoInputUpdate( TASK_LLEG_PORTNUMBER );
			if((pTLL1msIn->LL_I_StatWord_02 >>12)%2==1)
			{
				cout<<"02 Homing Posi Attained"<<endl;			
				break;
			}
		}

		//switch on
		pTLL1msOut->LL_O_CtrlWord_02=15;
		TCatIoOutputUpdate( TASK_LLEG_PORTNUMBER );
		Sleep(50);
	//Copy
	//Copy
		pTLL1msOut->LL_O_CtrlWord_01=6;
		TCatIoOutputUpdate( TASK_LLEG_PORTNUMBER );
		Sleep(50);
		//change mode
		pTLL1msOut->LL_O_SetMode_01=6;
		TCatIoOutputUpdate( TASK_LLEG_PORTNUMBER );
		Sleep(50);
		cout<<"01 HomingMode Set"<<endl;
		//switch on
		pTLL1msOut->LL_O_CtrlWord_01=15;
		TCatIoOutputUpdate( TASK_LLEG_PORTNUMBER );
		Sleep(50);
		//start homing
		//cout<<"Press Anykey to Start 01 Homing"<<endl;
		cout<<"Start 01 Homing ..."<<endl;
		Sleep(50);
		//system("pause");	
		pTLL1msOut->LL_O_CtrlWord_01=31;
		TCatIoOutputUpdate( TASK_LLEG_PORTNUMBER );
		Sleep(50);

		while(1)
		{
			TCatIoInputUpdate( TASK_LLEG_PORTNUMBER );
			if((pTLL1msIn->LL_I_StatWord_01 >>12)%2==1)
			{
				cout<<"01 Homing Posi Attained"<<endl;			
				break;
			}
		}

		//switch on
		pTLL1msOut->LL_O_CtrlWord_01=15;
		TCatIoOutputUpdate( TASK_LLEG_PORTNUMBER );
		Sleep(50);
	//Copy

		cout<<"Left Leg Homing Finished"<<endl;
		cout<<endl;

}
void TwinCAT_COM::EtherCATSetCSPRLeg(void)
{
	/******************************************************************
	input: void
	output: void

	Note:
	// 透過TwinCAT I/O 進行Set CSP Mode動作
	   組成方式為六組依序(由腳踝開始)的homing function
	   //Copy
				func ....
	   //Copy
	   模式啟動流程詳見EPOS3-Firmware Specification & Applicaion Note
	******************************************************************/
// Traj
	//Copy
		//shoutdown
		//pTRL1msOut->RL_O_CtrlWord_01=6;
		//TCatIoOutputUpdate( TASK_RLEG_PORTNUMBER );
		//Sleep(50);

		//change mode
		pTRL1msOut->RL_O_SetMode_01=8;
		TCatIoOutputUpdate( TASK_RLEG_PORTNUMBER );
		Sleep(50);
		//enable
		pTRL1msOut->RL_O_CtrlWord_01=15;
		TCatIoOutputUpdate( TASK_RLEG_PORTNUMBER );
		Sleep(50);


		TCatIoInputUpdate( TASK_RLEG_PORTNUMBER );
		Sleep(50);
		if(	pTRL1msIn->RL_I_Mode_01==8 	)
		cout<<"R01 CSP seted succeed"<<endl;
		else
		cout<<"R01 CSP seted failed"<<endl;
		cout<<endl;
	//Copy
	//Copy
		//shoutdown
		//pTRL1msOut->RL_O_CtrlWord_02=6;
		//TCatIoOutputUpdate( TASK_RLEG_PORTNUMBER );
		//Sleep(50);

		//change mode
		pTRL1msOut->RL_O_SetMode_02=8;
		TCatIoOutputUpdate( TASK_RLEG_PORTNUMBER );
		Sleep(50);
		//enable
		pTRL1msOut->RL_O_CtrlWord_02=15;
		TCatIoOutputUpdate( TASK_RLEG_PORTNUMBER );
		Sleep(50);


		TCatIoInputUpdate( TASK_RLEG_PORTNUMBER );
		Sleep(50);
		if(	pTRL1msIn->RL_I_Mode_02==8 	)
		cout<<"R02 CSP seted succeed"<<endl;
		else
		cout<<"R02 CSP seted failed"<<endl;
		cout<<endl;
	//Copy
	//Copy
		//shoutdown
		//pTRL1msOut->RL_O_CtrlWord_03=6;
		//TCatIoOutputUpdate( TASK_RLEG_PORTNUMBER );
		//Sleep(50);

		//change mode
		pTRL1msOut->RL_O_SetMode_03=8;
		TCatIoOutputUpdate( TASK_RLEG_PORTNUMBER );
		Sleep(50);
		//enable
		pTRL1msOut->RL_O_CtrlWord_03=15;
		TCatIoOutputUpdate( TASK_RLEG_PORTNUMBER );
		Sleep(50);


		TCatIoInputUpdate( TASK_RLEG_PORTNUMBER );
		Sleep(50);
		if(	pTRL1msIn->RL_I_Mode_03==8 	)
		cout<<"R03 CSP seted succeed"<<endl;
		else
		cout<<"R03 CSP seted failed"<<endl;
		cout<<endl;
	//Copy
	//Copy
		//shoutdown
		//pTRL1msOut->RL_O_CtrlWord_04=6;
		//TCatIoOutputUpdate( TASK_RLEG_PORTNUMBER );
		//Sleep(50);

		//change mode
		pTRL1msOut->RL_O_SetMode_04=8;
		TCatIoOutputUpdate( TASK_RLEG_PORTNUMBER );
		Sleep(50);
		//enable
		pTRL1msOut->RL_O_CtrlWord_04=15;
		TCatIoOutputUpdate( TASK_RLEG_PORTNUMBER );
		Sleep(50);


		TCatIoInputUpdate( TASK_RLEG_PORTNUMBER );
		Sleep(50);
		if(	pTRL1msIn->RL_I_Mode_04==8 	)
		cout<<"R04 CSP seted succeed"<<endl;
		else
		cout<<"R04 CSP seted failed"<<endl;
		cout<<endl;
	//Copy
	//Copy
		//shoutdown
		//pTRL1msOut->RL_O_CtrlWord_05=6;
		//TCatIoOutputUpdate( TASK_RLEG_PORTNUMBER );
		//Sleep(50);

		//change mode
		pTRL1msOut->RL_O_SetMode_05=8;
		TCatIoOutputUpdate( TASK_RLEG_PORTNUMBER );
		Sleep(50);
		//enable
		pTRL1msOut->RL_O_CtrlWord_05=15;
		TCatIoOutputUpdate( TASK_RLEG_PORTNUMBER );
		Sleep(50);


		TCatIoInputUpdate( TASK_RLEG_PORTNUMBER );
		Sleep(50);
		if(	pTRL1msIn->RL_I_Mode_05==8 	)
		cout<<"R05 CSP seted succeed"<<endl;
		else
		cout<<"R05 CSP seted failed"<<endl;
		cout<<endl;
	//Copy
	//Copy
		//shoutdown
		//pTRL1msOut->RL_O_CtrlWord_06=6;
		//TCatIoOutputUpdate( TASK_RLEG_PORTNUMBER );
		//Sleep(50);

		//change mode
		pTRL1msOut->RL_O_SetMode_06=8;
		TCatIoOutputUpdate( TASK_RLEG_PORTNUMBER );
		Sleep(50);
		//enable
		pTRL1msOut->RL_O_CtrlWord_06=15;
		TCatIoOutputUpdate( TASK_RLEG_PORTNUMBER );
		Sleep(50);


		TCatIoInputUpdate( TASK_RLEG_PORTNUMBER );
		Sleep(50);
		if(	pTRL1msIn->RL_I_Mode_06=8 	)
		cout<<"R06 CSP seted succeed"<<endl;
		else
		cout<<"R06 CSP seted failed"<<endl;
		cout<<endl;
	//Copy
}
void TwinCAT_COM::EtherCATSetCSPLLeg(void)
{
	/******************************************************************
	input: void
	output: void

	Note:
	// 透過TwinCAT I/O 進行Set CSP Mode動作
	   組成方式為六組依序(由腳踝開始)的homing function
	   //Copy
				func ....
	   //Copy
	   模式啟動流程詳見EPOS3-Firmware Specification & Applicaion Note
	******************************************************************/
// Traj
	//Copy
		//shoutdown
		//pTLL1msOut->LL_O_CtrlWord_01=6;
		//TCatIoOutputUpdate( TASK_LLEG_PORTNUMBER );
		//Sleep(50);

		//change mode
		pTLL1msOut->LL_O_SetMode_01=8;
		TCatIoOutputUpdate( TASK_LLEG_PORTNUMBER );
		Sleep(50);
		//enable
		pTLL1msOut->LL_O_CtrlWord_01=15;
		TCatIoOutputUpdate( TASK_LLEG_PORTNUMBER );
		Sleep(50);


		TCatIoInputUpdate( TASK_LLEG_PORTNUMBER );
		Sleep(50);
		if(	pTLL1msIn->LL_I_Mode_01==8 	)
		cout<<"L01 CSP seted succeed"<<endl;
		else
		cout<<"L01 CSP seted failed"<<endl;
		cout<<endl;
	//Copy
	//Copy
		//shoutdown
		//pTLL1msOut->LL_O_CtrlWord_02=6;
		//TCatIoOutputUpdate( TASK_LLEG_PORTNUMBER );
		//Sleep(50);

		//change mode
		pTLL1msOut->LL_O_SetMode_02=8;
		TCatIoOutputUpdate( TASK_LLEG_PORTNUMBER );
		Sleep(50);
		//enable
		pTLL1msOut->LL_O_CtrlWord_02=15;
		TCatIoOutputUpdate( TASK_LLEG_PORTNUMBER );
		Sleep(50);


		TCatIoInputUpdate( TASK_LLEG_PORTNUMBER );
		Sleep(50);
		if(	pTLL1msIn->LL_I_Mode_02==8 	)
		cout<<"L02 CSP seted succeed"<<endl;
		else
		cout<<"L02 CSP seted failed"<<endl;
		cout<<endl;
	//Copy
	//Copy
		//shoutdown
		//pTLL1msOut->LL_O_CtrlWord_03=6;
		//TCatIoOutputUpdate( TASK_LLEG_PORTNUMBER );
		//Sleep(50);

		//change mode
		pTLL1msOut->LL_O_SetMode_03=8;
		TCatIoOutputUpdate( TASK_LLEG_PORTNUMBER );
		Sleep(50);
		//enable
		pTLL1msOut->LL_O_CtrlWord_03=15;
		TCatIoOutputUpdate( TASK_LLEG_PORTNUMBER );
		Sleep(50);


		TCatIoInputUpdate( TASK_LLEG_PORTNUMBER );
		Sleep(50);
		if(	pTLL1msIn->LL_I_Mode_03==8 	)
		cout<<"L03 CSP seted succeed"<<endl;
		else
		cout<<"L03 CSP seted failed"<<endl;
		cout<<endl;
	//Copy
	//Copy
		//shoutdown
		//pTLL1msOut->LL_O_CtrlWord_04=6;
		//TCatIoOutputUpdate( TASK_LLEG_PORTNUMBER );
		//Sleep(50);

		//change mode
		pTLL1msOut->LL_O_SetMode_04=8;
		TCatIoOutputUpdate( TASK_LLEG_PORTNUMBER );
		Sleep(50);
		//enable
		pTLL1msOut->LL_O_CtrlWord_04=15;
		TCatIoOutputUpdate( TASK_LLEG_PORTNUMBER );
		Sleep(50);


		TCatIoInputUpdate( TASK_LLEG_PORTNUMBER );
		Sleep(50);
		if(	pTLL1msIn->LL_I_Mode_04==8 	)
		cout<<"L04 CSP seted succeed"<<endl;
		else
		cout<<"L04 CSP seted failed"<<endl;
		cout<<endl;
	//Copy
	//Copy
		//shoutdown
		//pTLL1msOut->LL_O_CtrlWord_05=6;
		//TCatIoOutputUpdate( TASK_LLEG_PORTNUMBER );
		//Sleep(50);

		//change mode
		pTLL1msOut->LL_O_SetMode_05=8;
		TCatIoOutputUpdate( TASK_LLEG_PORTNUMBER );
		Sleep(50);
		//enable
		pTLL1msOut->LL_O_CtrlWord_05=15;
		TCatIoOutputUpdate( TASK_LLEG_PORTNUMBER );
		Sleep(50);


		TCatIoInputUpdate( TASK_LLEG_PORTNUMBER );
		Sleep(50);
		if(	pTLL1msIn->LL_I_Mode_05==8 	)
		cout<<"L05 CSP seted succeed"<<endl;
		else
		cout<<"L05 CSP seted failed"<<endl;
		cout<<endl;
	//Copy
	//Copy
		//shoutdown
		//pTLL1msOut->LL_O_CtrlWord_06=6;
		//TCatIoOutputUpdate( TASK_LLEG_PORTNUMBER );
		//Sleep(50);

		//change mode
		pTLL1msOut->LL_O_SetMode_06=8;
		TCatIoOutputUpdate( TASK_LLEG_PORTNUMBER );
		Sleep(50);
		//enable
		pTLL1msOut->LL_O_CtrlWord_06=15;
		TCatIoOutputUpdate( TASK_LLEG_PORTNUMBER );
		Sleep(50);


		TCatIoInputUpdate( TASK_LLEG_PORTNUMBER );
		Sleep(50);
		if(	pTLL1msIn->LL_I_Mode_06=8 	)
		cout<<"L06 CSP seted succeed"<<endl;
		else
		cout<<"L06 CSP seted failed"<<endl;
		cout<<endl;
	//Copy
}
void TwinCAT_COM::EtherCATEmergentStop(void)
{
	/******************************************************************
	input: void
	output: void

	Note:
	// 透過TwinCAT I/O 進行Mode disalbe的動作
	   函式執行後機器人會整個鬆開軟腳(輸出終止),請留意執行時機
	   模式啟動流程詳見EPOS3-Firmware Specification & Applicaion Note
	******************************************************************/


	long nError;
	long nError2;
	static int iii=0; 
	if ( ((nError = TCatIoInputUpdate( TASK_RLEG_PORTNUMBER )) == 0 )&& ((nError2 = TCatIoOutputUpdate( TASK_RLEG_PORTNUMBER )) == 0 )&&
				(( nError = TCatIoInputUpdate( TASK_LLEG_PORTNUMBER )) == 0 )&& ( (nError2 = TCatIoOutputUpdate( TASK_LLEG_PORTNUMBER )) == 0 )
			) 
	{ 
		pTRL1msOut->RL_O_CtrlWord_01=6;
		pTRL1msOut->RL_O_CtrlWord_02=6;
		pTRL1msOut->RL_O_CtrlWord_03=6;
		pTRL1msOut->RL_O_CtrlWord_04=6;
		pTRL1msOut->RL_O_CtrlWord_05=6;
		pTRL1msOut->RL_O_CtrlWord_06=6;
		TCatIoOutputUpdate( TASK_RLEG_PORTNUMBER );
		Sleep(50);
		pTLL1msOut->LL_O_CtrlWord_01=6;
		pTLL1msOut->LL_O_CtrlWord_02=6;
		pTLL1msOut->LL_O_CtrlWord_03=6;
		pTLL1msOut->LL_O_CtrlWord_04=6;
		pTLL1msOut->LL_O_CtrlWord_05=6;
		pTLL1msOut->LL_O_CtrlWord_06=6;
		TCatIoOutputUpdate( TASK_LLEG_PORTNUMBER );

		cout<<"Emergent Stop!!!"<<endl;

	}
	else 
		printf( "TCatInputUpdate(%d) %d failed with 0x%x !\n",
			TASK_RLEG_PORTNUMBER, iii++, nError );

}
void TwinCAT_COM::EtherCATSetEncoder(long* buf)
{
	/******************************************************************
	input: long* buf
	output: void

	Note:
	// 將指標陣列儲存的數值透過TwinCAT I/O 進行指定Target Encoder動作
	   須將TwinCAT I/O Output 部分進行update才會開始動作
	******************************************************************/
	long LLegAxisEnable[6]={1,1,1,1,1,1};
	long RLegAxisEnable[6]={1,1,1,1,1,1};

	pTLL1msOut->LL_O_TarEnc_01=buf[0]*LLegAxisEnable[0];
	pTLL1msOut->LL_O_TarEnc_02=buf[1]*LLegAxisEnable[1];
	pTLL1msOut->LL_O_TarEnc_03=buf[2]*LLegAxisEnable[2];
	pTLL1msOut->LL_O_TarEnc_04=buf[3]*LLegAxisEnable[3];
	pTLL1msOut->LL_O_TarEnc_05=buf[4]*LLegAxisEnable[4];
	pTLL1msOut->LL_O_TarEnc_06=buf[5]*LLegAxisEnable[5]; 

	pTRL1msOut->RL_O_TarEnc_01=buf[6]*RLegAxisEnable[0];
	pTRL1msOut->RL_O_TarEnc_02=buf[7]*RLegAxisEnable[1];
	pTRL1msOut->RL_O_TarEnc_03=buf[8]*RLegAxisEnable[2];
	pTRL1msOut->RL_O_TarEnc_04=buf[9]*RLegAxisEnable[3];
	pTRL1msOut->RL_O_TarEnc_05=buf[10]*RLegAxisEnable[4];
	pTRL1msOut->RL_O_TarEnc_06=buf[11]*RLegAxisEnable[5];
}
void TwinCAT_COM::EtherCATSetVelOffset(long* buf)
{
	/******************************************************************
	input: long* buf
	output: void

	Note:
	// 將指標陣列儲存的數值透過TwinCAT I/O 進行指定Velocity offset動作
	   須將TwinCAT I/O Output 部分進行update才會開始動作
	******************************************************************/
	long LLegVelOffset[6]={1,1,1,1,1,1};
	long RLegVelOffset[6]={1,1,1,1,1,1};

	pTLL1msOut->LL_O_VelocityOffset_01=buf[0]*LLegVelOffset[0];
	pTLL1msOut->LL_O_VelocityOffset_02=buf[1]*LLegVelOffset[1];
	pTLL1msOut->LL_O_VelocityOffset_03=buf[2]*LLegVelOffset[2];
	pTLL1msOut->LL_O_VelocityOffset_04=buf[3]*LLegVelOffset[3];
	pTLL1msOut->LL_O_VelocityOffset_05=buf[4]*LLegVelOffset[4];
	pTLL1msOut->LL_O_VelocityOffset_06=buf[5]*LLegVelOffset[5]; 

	pTRL1msOut->RL_O_VelocityOffset_01=buf[6]*RLegVelOffset[0];
	pTRL1msOut->RL_O_VelocityOffset_02=buf[7]*RLegVelOffset[1];
	pTRL1msOut->RL_O_VelocityOffset_03=buf[8]*RLegVelOffset[2];
	pTRL1msOut->RL_O_VelocityOffset_04=buf[9]*RLegVelOffset[3];
	pTRL1msOut->RL_O_VelocityOffset_05=buf[10]*RLegVelOffset[4];
	pTRL1msOut->RL_O_VelocityOffset_06=buf[11]*RLegVelOffset[5];
}
void TwinCAT_COM::EtherCATSetTorqueOffset(short* buf)
{
	/******************************************************************
	input: long* buf
	output: void

	Note:
	// 將指標陣列儲存的數值透過TwinCAT I/O 進行指定Torque offset動作
	   須將TwinCAT I/O Output 部分進行update才會開始動作
	******************************************************************/
	short LLegTorqueOffset[6]={1,1,1,1,1,1};
	short RLegTorqueOffset[6]={1,1,1,1,1,1};

	pTLL1msOut->LL_O_TorqueOffset_01=buf[0]*LLegTorqueOffset[0];
	pTLL1msOut->LL_O_TorqueOffset_02=buf[1]*LLegTorqueOffset[1];
	pTLL1msOut->LL_O_TorqueOffset_03=buf[2]*LLegTorqueOffset[2];
	pTLL1msOut->LL_O_TorqueOffset_04=buf[3]*LLegTorqueOffset[3];
	pTLL1msOut->LL_O_TorqueOffset_05=buf[4]*LLegTorqueOffset[4];
	pTLL1msOut->LL_O_TorqueOffset_06=buf[5]*LLegTorqueOffset[5]; 

	pTRL1msOut->RL_O_TorqueOffset_01=buf[6]*RLegTorqueOffset[0];
	pTRL1msOut->RL_O_TorqueOffset_02=buf[7]*RLegTorqueOffset[1];
	pTRL1msOut->RL_O_TorqueOffset_03=buf[8]*RLegTorqueOffset[2];
	pTRL1msOut->RL_O_TorqueOffset_04=buf[9]*RLegTorqueOffset[3];
	pTRL1msOut->RL_O_TorqueOffset_05=buf[10]*RLegTorqueOffset[4];
	pTRL1msOut->RL_O_TorqueOffset_06=buf[11]*RLegTorqueOffset[5];
}
void TwinCAT_COM::EtherCATReadEncoder(long* buf)
{
	/******************************************************************
	input: long* buf
	output: void

	Note:
	// 透過TwinCAT I/O 進行讀取Encoder動作
	   須將TwinCAT I/O Input 部分進行update才會開始動作
	******************************************************************/
	buf[0]=pTLL1msIn->LL_I_Enc_01;
	buf[1]=pTLL1msIn->LL_I_Enc_02;
	buf[2]=pTLL1msIn->LL_I_Enc_03;
	buf[3]=pTLL1msIn->LL_I_Enc_04;
	buf[4]=pTLL1msIn->LL_I_Enc_05;
	buf[5]=pTLL1msIn->LL_I_Enc_06;

	buf[6]=pTRL1msIn->RL_I_Enc_01;
	buf[7]=pTRL1msIn->RL_I_Enc_02;
	buf[8]=pTRL1msIn->RL_I_Enc_03;
	buf[9]=pTRL1msIn->RL_I_Enc_04;
	buf[10]=pTRL1msIn->RL_I_Enc_05;
	buf[11]=pTRL1msIn->RL_I_Enc_06;
}
void TwinCAT_COM::EtherCATReadVel(long* buf)
{
	/******************************************************************
	input: long* buf
	output: void

	Note:
	// 透過TwinCAT I/O 進行讀取Velocity(RPM)動作
	   須將TwinCAT I/O Input 部分進行update才會開始動作	     
	******************************************************************/
	buf[0]=pTLL1msIn->LL_I_Vel_01;
	buf[1]=pTLL1msIn->LL_I_Vel_02;
	buf[2]=pTLL1msIn->LL_I_Vel_03;
	buf[3]=pTLL1msIn->LL_I_Vel_04;
	buf[4]=pTLL1msIn->LL_I_Vel_05;
	buf[5]=pTLL1msIn->LL_I_Vel_06;

	buf[6]=pTRL1msIn->RL_I_Vel_01;
	buf[7]=pTRL1msIn->RL_I_Vel_02;
	buf[8]=pTRL1msIn->RL_I_Vel_03;
	buf[9]=pTRL1msIn->RL_I_Vel_04;
	buf[10]=pTRL1msIn->RL_I_Vel_05;
	buf[11]=pTRL1msIn->RL_I_Vel_06;
}
void TwinCAT_COM::EtherCATReadTorque(short* buf)
{
	/******************************************************************
	input: long* buf
	output: void

	Note:
	// 透過TwinCAT I/O 進行讀取Torque(千分之rated torque)動作
	   須將TwinCAT I/O Input 部分進行update才會開始動作	     
	******************************************************************/
	buf[0]=pTLL1msIn->LL_I_Torque01;
	buf[1]=pTLL1msIn->LL_I_Torque02;
	buf[2]=pTLL1msIn->LL_I_Torque03;
	buf[3]=pTLL1msIn->LL_I_Torque04;
	buf[4]=pTLL1msIn->LL_I_Torque05;
	buf[5]=pTLL1msIn->LL_I_Torque06;

	buf[6]=pTRL1msIn->RL_I_Torque01;
	buf[7]=pTRL1msIn->RL_I_Torque02;
	buf[8]=pTRL1msIn->RL_I_Torque03;
	buf[9]=pTRL1msIn->RL_I_Torque04;
	buf[10]=pTRL1msIn->RL_I_Torque05;
	buf[11]=pTRL1msIn->RL_I_Torque06;
}
void TwinCAT_COM::EtherCATReadEncDiff(long* buf)
{
	/******************************************************************
	input: long* buf
	output: void

	Note:
	// 透過TwinCAT I/O 進行讀取Encoder Error動作
	   !!!CAUTION!!!須注意和Set coder放置的位置關係	     
	******************************************************************/
	buf[0]=pTLL1msOut->LL_O_TarEnc_01-pTLL1msIn->LL_I_Enc_01;
	buf[1]=pTLL1msOut->LL_O_TarEnc_02-pTLL1msIn->LL_I_Enc_02;
	buf[2]=pTLL1msOut->LL_O_TarEnc_03-pTLL1msIn->LL_I_Enc_03;
	buf[3]=pTLL1msOut->LL_O_TarEnc_04-pTLL1msIn->LL_I_Enc_04;
	buf[4]=pTLL1msOut->LL_O_TarEnc_05-pTLL1msIn->LL_I_Enc_05;
	buf[5]=pTLL1msOut->LL_O_TarEnc_06-pTLL1msIn->LL_I_Enc_06;

	buf[6]=pTRL1msOut->RL_O_TarEnc_01-pTRL1msIn->RL_I_Enc_01;
	buf[7]=pTRL1msOut->RL_O_TarEnc_02-pTRL1msIn->RL_I_Enc_02;
	buf[8]=pTRL1msOut->RL_O_TarEnc_03-pTRL1msIn->RL_I_Enc_03;
	buf[9]=pTRL1msOut->RL_O_TarEnc_04-pTRL1msIn->RL_I_Enc_04;
	buf[10]=pTRL1msOut->RL_O_TarEnc_05-pTRL1msIn->RL_I_Enc_05;
	buf[11]=pTRL1msOut->RL_O_TarEnc_06-pTRL1msIn->RL_I_Enc_06;
}
bool TwinCAT_COM::BitCheck(unsigned short statword,bool value, int bitnumb) //bitnumb count from 0
{
	/******************************************************************
	input: unsigned short statword,bool value, int bitnumb
	output: ture/false

	Note:
	// 進行特定位址之二近位數值大小判斷	     
	******************************************************************/
	if(   (statword>>bitnumb)%2==value        )
		return true;
	else
		return false;
}
bool TwinCAT_COM::DisabledCheck(unsigned short statword)
{
	/******************************************************************
	input: unsigned short statword
	output: ture/false

	Note:
	// 配合EtherCATInit使用
	   進行Statusword Disabled與否的狀態判讀
	   若是結果為disalbed(true),則可進行Initialization
	   狀態判別詳見Epos3 Firmware specifications & Application note
	******************************************************************/
/*	if(    BitCheck(statword,false,3)  &&  BitCheck(statword,false,2)  &&  BitCheck(statword,false,1)  &&  BitCheck(statword,false,0)  && 	
		   BitCheck(statword,true,6)  &&  BitCheck(statword,false,5)  &&  BitCheck(statword,false,4)  &&  
		   BitCheck(statword,true,8)  &&
		   BitCheck(statword,false,14)
	 )
		return true;
	else */if(    BitCheck(statword,false,3)  &&  BitCheck(statword,false,2)  &&  BitCheck(statword,false,1)  &&  BitCheck(statword,true,0)  && 	
		   BitCheck(statword,false,6)  &&  BitCheck(statword,true,5)  &&  BitCheck(statword,false,4)  &&  
		   BitCheck(statword,true,8)  &&
		   BitCheck(statword,false,14)
	 )
	 	return true;
	else if(    BitCheck(statword,false,3)  &&  BitCheck(statword,false,2)  &&  BitCheck(statword,true,1)  &&  BitCheck(statword,true,0)  && 	
		   BitCheck(statword,false,6)  &&  BitCheck(statword,true,5)  &&  BitCheck(statword,false,4)  &&  
		   BitCheck(statword,true,8)  &&
		   BitCheck(statword,false,14)
	 )
	 	return true;
	else
		return false;
}
unsigned char TwinCAT_COM::EPOS_ErrorCheck(unsigned short statword)
{
	/******************************************************************
	input: unsigned short statword
	output: ture/false

	Note:
	// 進行Statusword ERROR的狀態(i.e.是否有error)與種類判讀
	   狀態判別詳見Epos3 Firmware specifications & Application note
	******************************************************************/
	if(    BitCheck(statword,true,13))//following Fault
	 	return 1;
	else if(BitCheck(statword,true,3)  &&  BitCheck(statword,false,2)  &&  BitCheck(statword,false,1)  &&  BitCheck(statword,false,0)  && 	
		   BitCheck(statword,false,6)  &&  BitCheck(statword,false,5)  &&  BitCheck(statword,false,4)  &&  
		   BitCheck(statword,true,8)  &&
		   BitCheck(statword,false,14 )	 )//Fault 
	 	return 2;
	else if(	BitCheck(statword,true,3)  &&  BitCheck(statword,true,2)  &&  BitCheck(statword,true,1)  &&  BitCheck(statword,true,0)  && 	
			   BitCheck(statword,false,6)  &&  BitCheck(statword,false,5)  &&  BitCheck(statword,true,4)  &&  
			   BitCheck(statword,true,8)  &&
			   BitCheck(statword,false,14)//Fault Reaction Active (enable)
			   )
		return 3;
	else if(	BitCheck(statword,true,3)  &&  BitCheck(statword,true,2)  &&  BitCheck(statword,true,1)  &&  BitCheck(statword,true,0)  && 	
			   BitCheck(statword,false,6)  &&  BitCheck(statword,false,5)  &&  BitCheck(statword,false,4)  &&  
			   BitCheck(statword,true,8)  &&
			   BitCheck(statword,false,14)//Fault Reaction Active (disable)
			   )
		return 4;
	else
		return 0;//no error
}
void TwinCAT_COM:: MovingAve(int numb, long* buf, long* result)
{
	/******************************************************************
	input: int numb, long* buf,
	output: long* result

	Note:
	// 進行Moving Average的Filtering,對於random noise有良好的過濾效果 其實此一功能和EtherCAT運行無關
	   可移出此Class之外
	******************************************************************/
	for(int i =0;i<12;i++)
	MovingBuffer[(numb-1)*12+i]=buf[i];
	
	for(int j=0; j<12;j++)
	{
		for(int i =0;i<numb;i++)
		{
			result[j]+=MovingBuffer[i*12+j];
			
			if(i>0)
				MovingBuffer[(i-1)*12+j]=MovingBuffer[(i)*12+j];
		}
		result[j]=result[j]/numb;
	}
}
void TwinCAT_COM::EtherCATFaultReset(int index)
{
	/******************************************************************
	input: unsigned short statword
	output: ture/false

	Note:
	// 進行Fault Reset的動作 (如果錯誤可在不需重開電源排除且重置的話)
	   相關分類與Reset流程詳見Epos3 Firmware specifications & Application note
	******************************************************************/
	if(BitCheck(*pCtrlWords[index],false,7) )
	{
		*pCtrlWords[index]=*pCtrlWords[index]+128;
		if(index>=6)
			TCatIoOutputUpdate( TASK_RLEG_PORTNUMBER );
		else
			TCatIoOutputUpdate( TASK_LLEG_PORTNUMBER );
		cout<<"Axis  "<< index+1 <<"th Fault Reseted"<<endl;
	}
	else
		cout<<"Axis  "<< index+1 <<"th command aborted"<<endl;

}
TwinCAT_COM::~TwinCAT_COM(void)
{
	/******************************************************************
	input: void
	output: void

	Note:
	// Class destructor  關閉TwinCAT I/O
	******************************************************************/
	TCatIoClose(); 
}

#endif