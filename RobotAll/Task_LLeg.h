///////////////////////////////////////////////////////////////////////////////
// 
// Beckhoff Automation GmbH
// 
// TwinCAT IO HeaderFile
// 
///////////////////////////////////////////////////////////////////////////////
// 
// C:\Users\Robotlab\Desktop\Task_LLeg.h

#define	TASK_LLEG_PORTNUMBER	302

#define	TASK_LLEG_INPUTSIZE	78
#define	TASK_LLEG_OUTPUTSIZE	90

#pragma pack(push, 1)

typedef struct
{
	unsigned short	LL_I_StatWord_01;
	unsigned short	LL_I_StatWord_02;
	unsigned short	LL_I_StatWord_03;
	unsigned short	LL_I_StatWord_04;
	unsigned short	LL_I_StatWord_05;
	unsigned short	LL_I_StatWord_06;
	long	LL_I_Enc_01;
	long	LL_I_Enc_02;
	long	LL_I_Enc_03;
	long	LL_I_Enc_04;
	long	LL_I_Enc_05;
	long	LL_I_Enc_06;
	char	LL_I_Mode_01;
	char	LL_I_Mode_02;
	char	LL_I_Mode_03;
	char	LL_I_Mode_04;
	char	LL_I_Mode_05;
	char	LL_I_Mode_06;
	short	LL_I_Torque01;
	short	LL_I_Torque02;
	short	LL_I_Torque03;
	short	LL_I_Torque04;
	short	LL_I_Torque05;
	short	LL_I_Torque06;
	long	LL_I_Vel_01;
	long	LL_I_Vel_02;
	long	LL_I_Vel_03;
	long	LL_I_Vel_04;
	long	LL_I_Vel_05;
	long	LL_I_Vel_06;
} Task_LLeg_Inputs, *PTask_LLeg_Inputs;

typedef struct
{
	unsigned char	reserved1[12];
	long	LL_O_TarEnc_01;
	long	LL_O_TarEnc_02;
	long	LL_O_TarEnc_03;
	long	LL_O_TarEnc_04;
	long	LL_O_TarEnc_05;
	long	LL_O_TarEnc_06;
	char	LL_O_SetMode_01;
	char	LL_O_SetMode_02;
	char	LL_O_SetMode_03;
	char	LL_O_SetMode_04;
	char	LL_O_SetMode_05;
	char	LL_O_SetMode_06;
	unsigned short	LL_O_CtrlWord_01;
	unsigned short	LL_O_CtrlWord_02;
	unsigned short	LL_O_CtrlWord_03;
	unsigned short	LL_O_CtrlWord_04;
	unsigned short	LL_O_CtrlWord_05;
	unsigned short	LL_O_CtrlWord_06;
	short	LL_O_TorqueOffset_01;
	short	LL_O_TorqueOffset_02;
	short	LL_O_TorqueOffset_03;
	short	LL_O_TorqueOffset_04;
	short	LL_O_TorqueOffset_05;
	short	LL_O_TorqueOffset_06;
	long	LL_O_VelocityOffset_01;
	long	LL_O_VelocityOffset_02;
	long	LL_O_VelocityOffset_03;
	long	LL_O_VelocityOffset_04;
	long	LL_O_VelocityOffset_05;
	long	LL_O_VelocityOffset_06;
} Task_LLeg_Outputs, *PTask_LLeg_Outputs;

#pragma pack(pop)
