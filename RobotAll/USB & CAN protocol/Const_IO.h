

// ID parameters
#define MID 0
#define CMD_SET 16
#define Motor_Nodes 2 // the number of the motor nodes on the bus

// Master Output ID Command Definitions
#define ID_InitMech		MID*CMD_SET+0
#define ID_InitPara		MID*CMD_SET+1
#define ID_SetENC		MID*CMD_SET+2
#define ID_SetPID		MID*CMD_SET+3
#define ID_ModeIdle		MID*CMD_SET+4
#define ID_ModePID		MID*CMD_SET+5
#define ID_AllReturnENC MID*CMD_SET+6 // all motor nodes return their ENC counts
#define ID_ReadEncoder  MID*CMD_SET+11
//////////////////////  MID*CMD_SET+7 7號已經被定義為回傳 //////////////////////
#define ID_SetOneAxis	MID*CMD_SET+8
#define ID_SetIdle		MID*CMD_SET+9
#define ID_SetCaliENC 	MID*CMD_SET+10 // 設定機構誤差補償之encoder值

#define ID_SetPWMLimit		13

#define ID_EndOfContTraj	14
#define ID_SetTrajPacket	15


// ID Return to Master
#define IDR_Return_ENC	 0
#define IDR_LastCMD_Succ MID*CMD_SET+7 // 上一個指令完成，回傳告訴C32

// State Machine -- Primary State Definition
#define SP_IDLE			0
#define SP_PID			1
#define SP_INIT_MECH	2
#define SP_WATCH_NODES  3 // the master node watches the slave nodes
#define SP_SetOneAxis	4 // set the position of one specified axis
#define SP_SetIdle		5

// in C32
//#define SP_SetPIDVal	6 // set the PID value of an sepcified axis
//#define SP_FeedTraj		7 // feed trajectory
//#define SP_SetPWMLimit	8

// State Machine -- Secondary State -- Reply or Execute Definition
#define ShiftNum		20 // prevent the SPs and the SSs have the same numbers

#define SS_IDLE			0 + ShiftNum
#define SS_SET_ENC		1 + ShiftNum
#define SS_SET_PID		2 + ShiftNum
#define SS_INIT_PARA	3 + ShiftNum
//#define SS_READ_SEN1	4 + ShiftNum // example for new sensor reading state
//#define SS_SET_PARA1	5 + ShiftNum // example for parameter setting state

// constants
//#define FCY_Osc 		8000000 * 8   
//#define FCY_Mech 		8000000 * 2   
#define PI 3.1415926

#define INPUT			1
#define OUTPUT			0

// IO ports
//#define	PWM1L	LATEbits.LATE0
//#define	PWM1H	LATEbits.LATE1
//#define	PWM2L	LATEbits.LATE2
//#define	PWM2H	LATEbits.LATE3
//#define DirPWM1L	TRISEbits.TRISE0
//#define DirPWM1H	TRISEbits.TRISE1
//#define DirPWM2L	TRISEbits.TRISE2
//#define DirPWM2H	TRISEbits.TRISE3