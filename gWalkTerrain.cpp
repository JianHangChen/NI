void gWalkSlopeDora(void)
{
	/******************************************************************
	input: StepInput �`�B�� �]��preview�H�Υh�Y�h�� �ҥH�n� �i�H�o���`��B��, StepLength �B�Z�A�i���i�t
	output: void

	Note:
	// ��l�ƾ����H�欰������
	// �}�B�]�w �a�γ]�w �����H�B��Ѽ� �]�w��
	// �ۦ�����Ʀ� gInitTurnLeft gInitWalkStraight gInitStepHere gInitStair gInitSquat
	// ���שY �b�a�W��2�B�~��W�שY ���O�Ĥ@�B�N��W�שY
	// 20121214 doratom
	******************************************************************/
	gNumOfStep = 9; // �]�t��l�B�ഫ�Ppreivew���`�B��
	gCOGDown = 30; // �U�ֽ��\�U�� ������H �]����٤O

	slopeangle = -5.15*3.1415926/180; //�שY���פp��!!!�n�O�o�[�W�t��
	check_slopeangle = 1;
	rotate_pitch_time_ratio = 0.94;

	gKineAll.FlagSumoMode = 0;

	double x_val_zmp = 80; // ZMP ���k��V��m
	bool PNx = 0;
	double y_step_zmp =  200;
	double TurnRadius = 300;
	double distance2L = 160.0;

	gKineAll.StepHeight[0] = 0;
	
	for (int i = 1 ; i < gNumOfStep+30 ; i++)
		gKineAll.StepHeight[i] = 11.5;

	gKineAll.selSupport[0] = 2;
	for (int i = 1 ; i < gNumOfStep+30 ; i+=2)
	{
		gKineAll.selSupport[i] = 1;
		gKineAll.selSupport[i+1] = 0;
	}
	for (int i = gNumOfStep-4 ; i < 4000 ; i++)
		gKineAll.selSupport[i] = 2; // double support and prepare to stop

	double AngChange = 0.0;

	for (int i = 0 ; i < gNumOfStep+30 ; i++)
	{
		gLRotAngZ[i*2] = AngChange*i + gLAngZWorld;
		gLRotAngZ[i*2+1] = AngChange*i + gLAngZWorld;
	}

	gRRotAngZ[0] = gRAngZWorld;
	for (int i = 0 ; i < gNumOfStep+30 ; i++)
	{
		gRRotAngZ[i*2+1] = AngChange*i + gRAngZWorld;
		gRRotAngZ[i*2+2] = AngChange*i + gRAngZWorld;
	}

	gLLInitZMPFB = -80*sin(gRRotAngZ[0]); // ZMP Initial 
	gRLInitZMPFB = 80*sin(gRRotAngZ[0]);
	gLLInitZMPLR = 80*cos(gRRotAngZ[0]); // ZMP Initial Swing �y��|�Ψ� ���}�n��
	gRLInitZMPLR = -80*cos(gRRotAngZ[0]);

	double BodyDir = (gLAngZWorld+gRAngZWorld)/2.0;
	double StrideX = y_step_zmp*sin(BodyDir);
	double StrideY = y_step_zmp*cos(BodyDir);

	//R
     gRRotAngPitch[0] = 0.0;
	 gRRotAngPitch[1] = 0.0;
	 gRRotAngPitch[2] = 0.0;

	 for(int i = 3;i< gNumOfStep-4;i++)
      gRRotAngPitch[i] = slopeangle;

	 for(int i = gNumOfStep-4;i< gNumOfStep+30;i++)
		 gRRotAngPitch[i] = 0.0;

	 //L
	 gLRotAngPitch[0] = 0.0;
	 gLRotAngPitch[1] = 0.0;
	 gLRotAngPitch[2] = slopeangle;
	 gLRotAngPitch[3] = slopeangle;
	 
	 for(int i = 4;i< gNumOfStep-5;i++)
      gLRotAngPitch[i] = slopeangle;

	 for(int i = gNumOfStep-5;i< gNumOfStep+30;i++)
		 gLRotAngPitch[i] = 0.0;


	// ���k 
	gFstpX[0] = 0;
	gFstpX[1] = -x_val_zmp*cos(gRRotAngZ[1]); // ZMP Initial Swing �y��|�Ψ� ���}�n��
	gFstpX[2] = x_val_zmp*cos(gRRotAngZ[0])+StrideX/2.0; // ZMP Initial Swing �y��|�Ψ� ���}�n��
	for (int i = 2 ; i < gNumOfStep ; i++)
	{
		gFstpX[i*2-1] = gFstpX[(i-1)*2-1] + StrideX;
		gFstpX[i*2] = gFstpX[(i-1)*2] + StrideX;
	}

	if (gKineAll.selSupport[gNumOfStep-4] == 0) // right support
	{
		gFstpX[gNumOfStep-4] = gFstpX[gNumOfStep-5] + 160*cos((gLRotAngZ[gNumOfStep-4]+gRRotAngZ[gNumOfStep-4])/2.0);
	}
	else
	{
		gFstpX[gNumOfStep-4] = gFstpX[gNumOfStep-5] - 160*cos((gLRotAngZ[gNumOfStep-4]+gRRotAngZ[gNumOfStep-4])/2.0);
	}

	gFstpX[gNumOfStep-4] = (gFstpX[gNumOfStep-4]+gFstpX[gNumOfStep-5])/2.0;

	for (int i = gNumOfStep - 3 ; i < gNumOfStep+4 ; i++)
	{
		gFstpX[i] = gFstpX[gNumOfStep-4];
	}

	//gFstpY[0] = 0;
	//gFstpY[1] = distance2L/2.0*sin(gRRotAngZ[2]);
	//gFstpY[2] = -distance2L/2.0*sin(gLRotAngZ[0])+StrideY/2.0;
	//for (int i = 1 ; i < gNumOfStep ; i++)
	//{
	//	gFstpY[i*2+1] = gFstpY[i*2-1]+StrideY;
	//	gFstpY[i*2+2] = gFstpY[i*2]+StrideY;
	//}

	//if (gKineAll.selSupport[gNumOfStep-4] == 0) // right support
	//{
	//	gFstpY[gNumOfStep-4] = gFstpY[gNumOfStep-5] - 160*sin((gLRotAngZ[gNumOfStep-4]+gRRotAngZ[gNumOfStep-4])/2.0);
	//}
	//else
	//{
	//	gFstpY[gNumOfStep-4] = gFstpY[gNumOfStep-5] + 160*sin((gLRotAngZ[gNumOfStep-4]+gRRotAngZ[gNumOfStep-4])/2.0);
	//}

	//gFstpY[gNumOfStep-4] = (gFstpY[gNumOfStep-4]+gFstpY[gNumOfStep-5])/2.0;
	
	//20121217 �B�z���Ĥ@�B�M�̫�@�B
	gFstpY[0] = 0;
	gFstpY[1] = 0;
	gFstpY[2] = 210;
	gFstpY[3] = 300;
	gFstpY[4] = 400;
	gFstpY[5] = 500;
	gFstpY[6] = 500;
	gFstpY[7] = 500;
	
	
	//gFstpY[6] = 495+210;
	//gFstpY[7] = 495+210;


	for (int i = gNumOfStep - 3 ; i < gNumOfStep+5 ; i++)
	{
		gFstpY[i] = gFstpY[gNumOfStep-4];
	}

	//�H�W���򪽨��@�Ҥ@��//
	  
	//�]�w�n�񰪪�����//
	gGroundHeight[0] = 0.0;
	gGroundHeight[1] = 0.0;
	gGroundHeight[2] = (gFstpY[2]-114)*tan(-slopeangle);
	gGroundHeight[3] = (gFstpY[3]-114)*tan(-slopeangle);
	gGroundHeight[4] = (gFstpY[4]-114)*tan(-slopeangle);
	gGroundHeight[5] = (gFstpY[5]-114)*tan(-slopeangle);
	gGroundHeight[6] = (gFstpY[5])*tan(-slopeangle);
	
	/*for (int i = 3 ; i < gNumOfStep-4 ; i++)
		gGroundHeight[i] = gGroundHeight[i-1]+ (StrideY/2)*tan(-slopeangle);*/

	for (int i = gNumOfStep-4 ; i < gNumOfStep+25 ; i++)
		gGroundHeight[i] = gGroundHeight[i-1];

}


oid gWalkTerrain(int StepInput, double StepLength)
{
	/******************************************************************
	input: StepInput �`�B�� �]��preview�H�Υh�Y�h�� �ҥH�n� �i�H�o���`��B��, StepLength �B�Z�A�i���i�t
	output: void

	Note:

	// ��l�ƾ����H�欰��walk straight, the first 
	// �}�B�]�w �a�γ]�w �����H�B��Ѽ� �]�w��
	// �ۦ�����Ʀ� gInitTurnLeft gInitWalkStraight gInitStepHere gInitStair gInitSquat gInitDownStair

	******************************************************************/
	gNumOfStep = StepInput; // �]�t��l�B�ഫ�Ppreivew���`�B��
	gCOGDown = 25;//9;//36; // �U�ֽ��\�U�� ������H �]����٤O
	if(gNumOfStep==7)//Specific Ver. Warn.
		checkonestep = 1;
	gKineAll.FlagSumoMode = 0;

	double x_val_zmp = 80; // ZMP ���k��V��m
	bool PNx = 0;
	double y_step_zmp =  StepLength;;
	double TurnRadius = 300;


	double distance2L = 160.0;

	gKineAll.StepHeight[0] = 0;
	for (int i = 1 ; i < gNumOfStep+30 ; i++){
		if (AKHSObstacle == true){
			gKineAll.StepHeight[i] = 30;	// 0516 WZ �� �쥻20
		} 
		else{
			gKineAll.StepHeight[i] = 20;
		}	
	}


	gKineAll.selSupport[0] = 2;
	for (int i = 1 ; i < gNumOfStep+30 ; i+=2)
	{
		gKineAll.selSupport[i] = 1;
		gKineAll.selSupport[i+1] = 0;
	}
	for (int i = gNumOfStep-4 ; i < 4000 ; i++)
		gKineAll.selSupport[i] = 2; // double support and prepare to stop


	double AngChange = 0.0;

	for (int i = 0 ; i < gNumOfStep+30 ; i++)
	{
		gLRotAngZ[i*2] = AngChange*i + gLAngZWorld;
		gLRotAngZ[i*2+1] = AngChange*i + gLAngZWorld;
	}

	gRRotAngZ[0] = gRAngZWorld;
	for (int i = 0 ; i < gNumOfStep+30 ; i++)
	{
		gRRotAngZ[i*2+1] = AngChange*i + gRAngZWorld;
		gRRotAngZ[i*2+2] = AngChange*i + gRAngZWorld;
	}

	gLLInitZMPFB = -80*sin(gRRotAngZ[0]); // ZMP Initial 
	gRLInitZMPFB = 80*sin(gRRotAngZ[0]);
	gLLInitZMPLR = 80*cos(gRRotAngZ[0]); // ZMP Initial Swing �y��|�Ψ� ���}�n��
	gRLInitZMPLR = -80*cos(gRRotAngZ[0]);

	double BodyDir = (gLAngZWorld+gRAngZWorld)/2.0;
	double StrideX = y_step_zmp*sin(BodyDir);
	double StrideY = y_step_zmp*cos(BodyDir);


	// ���k 
	gFstpX[0] = 0;
	gFstpX[1] = -x_val_zmp*cos(gRRotAngZ[1]); // ZMP Initial Swing �y��|�Ψ� ���}�n��
	gFstpX[2] = x_val_zmp*cos(gRRotAngZ[0])+StrideX/2.0; // ZMP Initial Swing �y��|�Ψ� ���}�n��
	for (int i = 2 ; i < gNumOfStep ; i++)
	{
		gFstpX[i*2-1] = gFstpX[(i-1)*2-1] + StrideX;
		gFstpX[i*2] = gFstpX[(i-1)*2] + StrideX;
	}

	if (gKineAll.selSupport[gNumOfStep-4] == 0) // right support
	{
		gFstpX[gNumOfStep-4] = gFstpX[gNumOfStep-5] + 160*cos((gLRotAngZ[gNumOfStep-4]+gRRotAngZ[gNumOfStep-4])/2.0);
	}
	else
	{
		gFstpX[gNumOfStep-4] = gFstpX[gNumOfStep-5] - 160*cos((gLRotAngZ[gNumOfStep-4]+gRRotAngZ[gNumOfStep-4])/2.0);
	}

	gFstpX[gNumOfStep-4] = (gFstpX[gNumOfStep-4]+gFstpX[gNumOfStep-5])/2.0;

	for (int i = gNumOfStep - 3 ; i < gNumOfStep+4 ; i++)
	{
		gFstpX[i] = gFstpX[gNumOfStep-4];
	}

	if(gNumOfStep ==7)
	{
      for (int i = gNumOfStep-4; i < 4000 ; i++)
		gFstpX[i] = 0.0;

	    gFstpX[3] =gFstpX[1];
	}
	//else
	//{
	//	for (int i = gNumOfStep-4; i < 4000 ; i++)
	//	gFstpX[i] = 0.0;
	//}



	gFstpY[0] = 0;
	gFstpY[1] = distance2L/2.0*sin(gRRotAngZ[2]);
	gFstpY[2] = -distance2L/2.0*sin(gLRotAngZ[0])+StrideY/2.0;
	//gFstpY[0] = 0;
	//gFstpY[1] = 0;
	//gFstpY[2] = 50;
	//gFstpY[3] = 260;
	//	gFstpY[4] = 310;
	//		gFstpY[5] = 520;
	//			gFstpY[6] = 570;
	//				gFstpY[7] = 570;
	for (int i = 1 ; i < gNumOfStep ; i++)
	{
		gFstpY[i*2+1] = gFstpY[i*2-1]+StrideY;
		gFstpY[i*2+2] = gFstpY[i*2]+StrideY;
	}

	if (gKineAll.selSupport[gNumOfStep-4] == 0) // right support
	{
		gFstpY[gNumOfStep-4] = gFstpY[gNumOfStep-5] - 160*sin((gLRotAngZ[gNumOfStep-4]+gRRotAngZ[gNumOfStep-4])/2.0);
	}
	else
	{
		gFstpY[gNumOfStep-4] = gFstpY[gNumOfStep-5] + 160*sin((gLRotAngZ[gNumOfStep-4]+gRRotAngZ[gNumOfStep-4])/2.0);
	}

	gFstpY[gNumOfStep-4] = (gFstpY[gNumOfStep-4]+gFstpY[gNumOfStep-5])/2.0;

	for (int i = gNumOfStep - 3 ; i < gNumOfStep+5 ; i++)
	{
		gFstpY[i] = gFstpY[gNumOfStep-4];
	}

	for (int i = 0 ; i < gNumOfStep+25 ; i++)
		gGroundHeight[i] = 0.0;

	//***************20141022�ܨB�Z²�����
	gFstpY[3]=220;
}