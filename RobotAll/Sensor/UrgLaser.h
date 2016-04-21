#ifndef _URGLASER_H_
#define _URGLASER_H_
#define _USE_MATH_DEFINES
#include<windows.h>
#include<cmath>
using namespace std;

class UrgLaser
{
	HANDLE g_commHandle;
	BOOL CommOpen(LPCTSTR ComPortName);

	BOOL CommClose();

public:
	UrgLaser():g_commHandle(INVALID_HANDLE_VALUE)
	{
	}
	BOOL Open(LPCTSTR port)
	{
		return CommOpen(port);
	}
	BOOL Close()
	{
		return CommClose();
	}
	BOOL Scan769(int *ranges);
};
#endif