/*
	Serial port class designed specifically for SickLMS2xx.

	Platform: Linux, Windows (Win32?)
	License: LGPL

	Written by Jim (r94922159@ntu.edu.tw), 2005-12-10. (tab=4)
	Modified by Jim, 2005-12-14: make it MSVC-compliant.
	Modified by Jim, 2005-12-15: make it restartable...
	Modified by Jim, 2006-08-06: (minor) display Tx/Rx queue sizes
	Modified by Jim, 2006-08-18: (minor?) changed some communication settings
	Modified by Jim, 2006-08-28: (minor) added (commented) ClearCommBreak() & ClearCommError()
	Modified by Jim, 2006-08-28: separate into .cpp and .h files
	Modified by Jimmy Yen, 2007-04-13: (minor) hacked _set_baudrate() and _write() to be public

	Blah..blah...please refer to sicklms2xx.cpp.

	References:
	- Serial Communications in Win32, MSDN (http://tinyurl.com/8un6)
	- SICK's MSTDemo: UserPort.*
*/

#ifndef __SERIAL_PORT_H
#define __SERIAL_PORT_H

#include "util.h"
#include <fcntl.h>

#ifdef LINUX
	#include <unistd.h>
	#include <termios.h>
	#include <sys/ioctl.h>
	#include <sys/utsname.h>
	#include <linux/serial.h>

	typedef int DevHandle;
	#define BADFD (-1)
#else
	#define SERIAL_IN_BUFFER		(30*1024)
	//#define SERIAL_OUT_BUFFER		(1*1024)
	#define SERIAL_OUT_BUFFER		(4104)

	typedef HANDLE DevHandle;
	#define BADFD INVALID_HANDLE_VALUE
#endif

//----------------------------------------------------------------------------

// timeouts/delays in ms, except for TIMEOUT_CHANGE_BAUDRATE and WRITE_BYTE_MIN_DELAY
// The TIMEOUTs are obeyed by the SICK device. Others are obeyed by the PC.
//#define TIMEOUT_READ_BYTE			(14999+2)		// max delay b/t two received bytes
//#define TIMEOUT_AVG_READ_BYTE       4999			// according to MSTDemo
#define TIMEOUT_READ_BYTE			(14+2)		// max delay b/t two received bytes
#define TIMEOUT_AVG_READ_BYTE       4			// according to MSTDemo
//TODO #define WRITE_BYTE_MAX_DELAY	6			// max delay b/t two transmitted bytes -- thus you cannot write slower than ~1667 baud -- after the delay, the telegram is ignored
#define WRITE_BYTE_MIN_DELAY		(55+10)		// us -- min delay b/t two transmitted bytes
#define TIMEOUT_ACK_NAK				(60+5)
#define NAK_MIN_DELAY				(30+2)		// min delay after received NAK & before transmitting
#define TIMEOUT_CHANGE_BAUDRATE		3.0			// s
//TODO #define TIMEOUT_READ_MEASUREMENTS	60	// what is this..?
#define TIMEOUT_RECEIVE_SCAN		(60+3)		// max response time for a request of scan of 0.25deg resolution
#define TIMEOUT_MISC_RESPONSE		60			// guessed

//----------------------------------------------------------------------------

extern const int baudrate_table[][4];
enum SerialCardType { NORMAL_RS232=0, PCMCIA_RS422=1, PCI_RS422=2 };


class SerialPort {
	DevHandle fd;
	int _current_baudrate;

	#ifdef LINUX
	timespec write_delay;		// extra delay required after writing each byte
	#endif
	int shortnap_loops;
	void shortnap() const { for(int j=shortnap_loops;--j>=0;) for(int i=100;--i;); }		// wastes about 130,000 CPU cycles on a 2GHz CPU
public:
	bool _set_baudrate(int);
	bool _write (const uchar *data, int len) const;

	SerialPort();
	~SerialPort()            { close(); }
	bool opened() const      { return fd!=BADFD; }
	bool open(const char *, SerialCardType);
	bool close();
	bool set_baudrate(SerialCardType,int);
	bool cancel_io();
	int current_baudrate() const { return _current_baudrate; }

	bool read (uchar *buf, int len, int timeout_ms, const char *reason) const;
	bool write  (const uchar *data, int len) const;
	bool write  (uchar data) const { return _write(&data,1); }

protected:
	const char* m_dev_file;
};

#endif
