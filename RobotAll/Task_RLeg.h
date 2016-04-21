///////////////////////////////////////////////////////////////////////////////
// 
// Beckhoff Automation GmbH
// 
// TwinCAT IO HeaderFile
// 
///////////////////////////////////////////////////////////////////////////////
// 
// C:\Users\Robotlab\Desktop\Task_RLeg.h

#define	TASK_RLEG_PORTNUMBER	301

#define	TASK_RLEG_INPUTSIZE	78
#define	TASK_RLEG_OUTPUTSIZE	78

#pragma pack(push, 1)

typedef struct
{
	unsigned short	RL_I_StatWord_01;
	unsigned short	RL_I_StatWord_02;
	unsigned short	RL_I_StatWord_03;
	unsigned short	RL_I_StatWord_04;
	unsigned short	RL_I_StatWord_05;
	unsigned short	RL_I_StatWord_06;
	long	RL_I_Enc_01;
	long	RL_I_Enc_02;
	long	RL_I_Enc_03;
	long	RL_I_Enc_04;
	long	RL_I_Enc_05;
	long	RL_I_Enc_06;
	char	RL_I_Mode_01;
	char	RL_I_Mode_02;
	char	RL_I_Mode_03;
	char	RL_I_Mode_04;
	char	RL_I_Mode_05;
	char	RL_I_Mode_06;
	short	RL_I_Torque01;
	short	RL_I_Torque02;
	short	RL_I_Torque03;
	short	RL_I_Torque04;
	short	RL_I_Torque05;
	short	RL_I_Torque06;
	long	RL_I_Vel_01;
	long	RL_I_Vel_02;
	long	RL_I_Vel_03;
	long	RL_I_Vel_04;
	long	RL_I_Vel_05;
	long	RL_I_Vel_06;
} Task_RLeg_Inputs, *PTask_RLeg_Inputs;

typedef struct
{
	unsigned short	RL_O_CtrlWord_01;
	unsigned short	RL_O_CtrlWord_02;
	unsigned short	RL_O_CtrlWord_03;
	unsigned short	RL_O_CtrlWord_04;
	unsigned short	RL_O_CtrlWord_05;
	unsigned short	RL_O_CtrlWord_06;
	long	RL_O_TarEnc_01;
	long	RL_O_TarEnc_02;
	long	RL_O_TarEnc_03;
	long	RL_O_TarEnc_04;
	long	RL_O_TarEnc_05;
	long	RL_O_TarEnc_06;
	char	RL_O_SetMode_01;
	char	RL_O_SetMode_02;
	char	RL_O_SetMode_03;
	char	RL_O_SetMode_04;
	char	RL_O_SetMode_05;
	char	RL_O_SetMode_06;
	short	RL_O_TorqueOffset_01;
	short	RL_O_TorqueOffset_02;
	short	RL_O_TorqueOffset_03;
	short	RL_O_TorqueOffset_04;
	short	RL_O_TorqueOffset_05;
	short	RL_O_TorqueOffset_06;
	long	RL_O_VelocityOffset_01;
	long	RL_O_VelocityOffset_02;
	long	RL_O_VelocityOffset_03;
	long	RL_O_VelocityOffset_04;
	long	RL_O_VelocityOffset_05;
	long	RL_O_VelocityOffset_06;
} Task_RLeg_Outputs, *PTask_RLeg_Outputs;

#pragma pack(pop)
