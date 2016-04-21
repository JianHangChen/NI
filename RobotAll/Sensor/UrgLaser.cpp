#include "stdafx.h"
#define _USE_MATH_DEFINES

#include<windows.h>
#include<cmath>
#include"UrgLaser.h"
using namespace std;

BOOL UrgLaser::CommOpen(LPCTSTR ComPortName)
{
	DCB	dcb;
	BOOL retVal;

//	g_CommOverlapped.hEvent = CreateEvent(NULL, TRUE, FALSE, NULL);
	//g_ReadOverlapped.hEvent = CreateEvent(NULL, TRUE, FALSE, NULL);
	//g_WriteOverlapped.hEvent = CreateEvent(NULL, TRUE, FALSE, NULL);

	// OPEN THE SERIAL OR USB COM PORT(CONNECTED TO THE SENSOR)
	// REPLACE "COM2" WITH A STRING OR "COM1", "COM3", ETC. TO OPEN ANOTHER PORT. 
	
	//g_commHandle = CreateFile("COM3",GENERIC_READ|GENERIC_WRITE,0,NULL,OPEN_EXISTING,FILE_FLAG_OVERLAPPED,NULL);
	g_commHandle = CreateFile(ComPortName,GENERIC_READ|GENERIC_WRITE,0,NULL,OPEN_EXISTING,0,NULL);

	// Communication Settings
	dcb.DCBlength = sizeof(DCB);
	dcb.BaudRate = CBR_115200;
	//dcb.BaudRate=CBR_256000;
	dcb.fBinary = TRUE;
	dcb.fParity = FALSE;
	dcb.fOutxCtsFlow = FALSE;
	dcb.fOutxDsrFlow = FALSE;
	dcb.fDtrControl = DTR_CONTROL_DISABLE;
	dcb.fDsrSensitivity = FALSE;
	dcb.fTXContinueOnXoff = FALSE;
	dcb.fOutX = FALSE;
	dcb.fInX = FALSE;
	dcb.fErrorChar = FALSE;
	dcb.fNull = FALSE;
	dcb.fRtsControl = RTS_CONTROL_DISABLE;
	dcb.fAbortOnError = FALSE;
	dcb.fDummy2 = 0;
	dcb.wReserved = 0;
	dcb.XonLim = 0;
	dcb.XoffLim = 0;
	dcb.ByteSize = 8;
	dcb.Parity = NOPARITY;
	dcb.StopBits = ONESTOPBIT;
	dcb.XonChar = 0;
	dcb.XoffChar = 0;
	dcb.ErrorChar = 0;
	dcb.EofChar = 0;
	dcb.EvtChar = 0;
	dcb.wReserved1 =0;

	retVal=true;
	retVal &= SetCommState(g_commHandle,&dcb);		// Settings for DCB structure details for serial port
	retVal &= SetCommMask(g_commHandle,EV_RXCHAR);	// Event settings (Read event) 
	retVal &= SetupComm(g_commHandle,0x1000,0x1000);	// Settings for input/ouput buffer size 

	return retVal;
}

BOOL UrgLaser::CommClose()
{
	return CloseHandle(g_commHandle);
}

BOOL UrgLaser::Scan769(int *ranges)
{
	char buffer[0x1000];
	strcpy(buffer,"G00076801\x0a\x00");
	DWORD d;
	int len=0;
	BOOL WriteOK=FALSE;
	int WriteNum=0;

	// Write Command
	while(WriteNum<strlen(buffer))
	{
		if(!(WriteOK=WriteFile(g_commHandle,buffer,strlen(buffer),&d,NULL)))
			break;
		WriteNum+=d;
	}

	if(!WriteOK)
	{
		return FALSE;
	}

	DWORD Event;
	bool Fail=false;
	while(true)
	{
		DWORD ReadCount=0;
		Sleep(0);
		WaitCommEvent(g_commHandle,&Event,NULL);//&g_ReadOverlapped);
		if(Event & EV_RXCHAR)
		{
			if(!ReadFile(g_commHandle,buffer+len,sizeof(buffer)-len,&ReadCount,NULL))
			{
				Fail=true;
				break;
			}
			len+=ReadCount;						
			if(buffer[len-1]=='\n' && buffer[len-2]=='\n')
				break;
		}
	}
	if(Fail) return FALSE;
	
	int ProcessingIndex=0,CurrentIndex=0;

	while(buffer[ProcessingIndex]!='\n') ProcessingIndex++;
	ProcessingIndex++;
	while(buffer[ProcessingIndex]!='\n') ProcessingIndex++;
	ProcessingIndex++;

	while(ProcessingIndex<len)
	{
		if(buffer[ProcessingIndex]!='\n')
			buffer[CurrentIndex++]=buffer[ProcessingIndex];
		ProcessingIndex++;
	}
	len=CurrentIndex;
	for(int i=0;i<len;i+=2)
	{
		int value=buffer[i]-'0';
		value=(value<<6)|((buffer[i+1]-'0'));
		ranges[i>>1]=value;
	}
	return TRUE;
}