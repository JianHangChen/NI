//#include "../stdafx.h"
#include "stdafx.h"
#include "serial_port.h"
//#include <boost/thread/mutex.hpp>

//static boost::mutex io_mutex;
//----------------------------------------------------------------------------

const int baudrate_table[][4] = {
	{ 9600, 19200, 38400, 500000 },		// Normal RS-232
	{ 2150,  4301,  8602, 115200 },		// PCMCIA RS-422
	{ 1075,  2216,  4300,  56000 },		// PCI    RS-422 (ex: Quatech DSC-200, with 8x multiplier)
};


SerialPort::SerialPort() : fd(BADFD), _current_baudrate(0)
{
	#ifdef LINUX
	write_delay.tv_sec=0;
	write_delay.tv_nsec=0;
	#endif
	
	// calibrate shortnap_loops, in 1 second
	info("Calibrating write delay (1 sec)...");
	shortnap_loops=0;
	for_n_seconds(1.0) {
		for(int i=100000;--i;);
		shortnap_loops++;
	}
	shortnap_loops = WRITE_BYTE_MIN_DELAY*shortnap_loops/1000;
}


// settings: 8 N 1
bool SerialPort::open (const char *devfile, SerialCardType cardtype)
{
	if(opened()) { warn("already opened"); close(); }
	m_dev_file = devfile;
#ifdef LINUX
	fd = ::open (devfile, O_RDWR | O_NOCTTY | O_SYNC |	// O_SYNC is for the 55us requirement
				O_ASYNC | O_NONBLOCK);					// O_ASYNC and O_NONBLOCK are for the 14ms requirement
	if(fd<0) { infoe("unable to open file (%s)",devfile); fd=BADFD; return false; }
#else
	fd = CreateFile (devfile, GENERIC_READ|GENERIC_WRITE, 0/*FILE_SHARE_READ|FILE_SHARE_WRITE|FILE_SHARE_DELETE*/, 0, OPEN_EXISTING, 0 /*FILE_FLAG_OVERLAPPED*/, 0);
	if(fd==BADFD) { warnE("CreateFile"); return false; }
	////// DWORD err=0; ClearCommError(fd, &err, 0) || warnE("ClearCommError");	info("err=%X",err); 	// probably make things better, but untested
	////// ClearCommBreak(fd) || warnE("ClearCommBreak");												// probably make things better, but untested
	if (!SetupComm (fd, SERIAL_IN_BUFFER, SERIAL_OUT_BUFFER)) { warnE("SetupComm"); return false; }
	COMMPROP cp;   GetCommProperties(fd, &cp) || warnE("GetCommProperties");
	info("Current queue sizes: TX=%d(max %d), RX=%d(max %d)", cp.dwCurrentTxQueue, cp.dwMaxTxQueue, cp.dwCurrentRxQueue, cp.dwMaxRxQueue);
	// TODO: remember current state, so that it can be restored at close().
#endif
	set_baudrate(cardtype,9600);
	info("serial port opened (%s)",devfile);
	return true;
}


bool SerialPort::close()
{
	if(!opened()) return true;
	cancel_io();
	bool ret;
#ifdef LINUX
	(ret=(::close(fd)==0)) || warne("close");
#else
	(ret=CloseHandle(fd)!=0) || warnE("CloseHandle");
#endif
	fd=BADFD;
	_current_baudrate=0;
	return ret;
}


#ifdef LINUX
static bool kernel_minimum_version (int a, int b, int c)
{
	int ca, cb, cc;
	utsname uts; uname(&uts);
	sscanf (uts.release, "%d.%d.%d", &ca, &cb, &cc);
	return (ca*0x10000+cb*0x100+cc >= a*0x10000+b*0x100+c);
}
#endif


bool SerialPort::set_baudrate (SerialCardType ct, int brate)
{
	if(!opened()) { warn("device not opened yet"); return false; }

	// baudrate conversion (for those RS-422 cards)
	int brate2=0;			// converted baudrate
	for(int i=NELEMS(baudrate_table[ct]); --i>=0;)
		if(brate==baudrate_table[NORMAL_RS232][i]) brate2=baudrate_table[ct][i];
	if(brate2==0) { warn("invalid baudrate (ct=%d, brate=%d)",ct,brate); return false; }

	cancel_io();
	if(!_set_baudrate(brate2)) return false;
	_current_baudrate=brate;
	verbose( info( "baudrate set to %d (internally coded as %d baud)", brate, brate2 ));

	#ifdef LINUX
	write_delay.tv_nsec = (WRITE_BYTE_MIN_DELAY - 8000000/brate) * 1000;
	#endif
	msleep(10);	// (HACK) waiting for the port to stablize?
	return true;
}


bool SerialPort::_set_baudrate (int brate2)
{
#ifdef LINUX
	// set custom speed
	assert (kernel_minimum_version(2,4,20));
	serial_struct serinfo;
	serinfo.reserved_char[0] = 0;
	ioctl(fd, TIOCGSERIAL, &serinfo)>=0 || die("ioctl(TIOCGSERIAL)");
	serinfo.flags &= ~ASYNC_SPD_MASK;
	serinfo.flags |= ASYNC_SPD_CUST;
	serinfo.custom_divisor =  serinfo.baud_base / brate2;
	ioctl(fd, TIOCSSERIAL, &serinfo)>=0 || die("ioctl(TIOCSSERIAL)");

	// set terminal parameters
	termios term;
	if(tcgetattr (fd, &term)!=0) { warne("tcgetattr"); return false; }
	cfmakeraw (&term);
	term.c_iflag &= ~(INPCK|IUCLC|IXOFF|IXANY|IGNPAR);		// raw => iflag &= ~(IGNBRK|BRKINT|PARMRK|ISTRIP|INLCR|IGNCR|ICRNL|IXON);
	term.c_oflag = 0;										// raw => oflag &= ~OPOST;
	term.c_lflag &= ~(XCASE|ECHOE|ECHOK|ECHOCTL|ECHOPRT|ECHOKE|TOSTOP|PENDIN);	// raw => lflag &= ~(ECHO|ECHONL|ICANON|ISIG|IEXTEN);
	term.c_lflag = 0;
	term.c_cflag &= ~(CSTOPB|CRTSCTS|HUPCL|CRTSCTS);		// raw => cflag &= ~(CSIZE|PARENB);
	term.c_cflag |= (CREAD|CLOCAL);							// raw => cflag |= CS8;         NOTE: what is CREAD?
	term.c_cc[VMIN] = 0;
	term.c_cc[VTIME] = 0;
	cfsetispeed (&term, B38400);
	cfsetospeed (&term, B38400);
	tcflush (fd, TCIOFLUSH);
	if(tcsetattr (fd, TCSAFLUSH, &term)!=0) { warne("unable to set baudrate (%d): tcsetattr",brate2); return false; }
	tcflush (fd, TCIOFLUSH);
#else
	DCB dcb;
	if(!GetCommState (fd, &dcb)) { warnE("GetCommState"); return false; }
	dcb.DCBlength = sizeof(dcb);
	dcb.fBinary = true;
	dcb.fOutxCtsFlow = false;
	dcb.fOutxDsrFlow = false;
	dcb.fDsrSensitivity = false;
	dcb.fDtrControl = DTR_CONTROL_ENABLE;
	dcb.fRtsControl = RTS_CONTROL_ENABLE;
	dcb.fTXContinueOnXoff = true;
	dcb.fInX = false;
	dcb.fOutX = false;
	dcb.fNull = false;
	dcb.fErrorChar = false;
	dcb.fAbortOnError = false;		// NOTE: If you enable this, any error will abort current operations and forbid future operations, until ClearCommError() is called.
	dcb.BaudRate = brate2;
	dcb.ByteSize = 8;
	dcb.fParity = false;
	dcb.Parity = NOPARITY;
	dcb.StopBits = ONESTOPBIT;
	if(!SetCommState (fd, &dcb)) { warnE("SetCommState"); return false; }

	COMMTIMEOUTS cto;
	cto.ReadIntervalTimeout = TIMEOUT_READ_BYTE;
	cto.ReadTotalTimeoutMultiplier = 0;
	cto.ReadTotalTimeoutConstant = 0;
	cto.WriteTotalTimeoutMultiplier = 0;
	cto.WriteTotalTimeoutConstant = 0;
	if (!SetCommTimeouts(fd,&cto)) { warnE("SetTimeoutProperties"); return false; }
#endif
	return true;
}

// NOTE: An interval of up to 14 ms can elapse between two bytes sent from the LMS 2xx within a telegram.
// NOTE: This (select method) could possibly lead to context-switch for every byte read, degrading system performance.
// later note: my bad.. a large interval is actually easier on other programs' performance (especially when this has real-time priority).
bool SerialPort::read (uchar *buf, int len, int timeout_ms, const char *reason) const
{
#ifdef LINUX
	static timeval timeout_general = { TIMEOUT_READ_BYTE/1000, (TIMEOUT_READ_BYTE%1000)*1000 };
	timeval timeout = { timeout_ms/1000, (timeout_ms%1000)*1000 };
	bool dirty=false;

	fd_set fds;
	FD_ZERO (&fds);
	FD_SET (fd+1, &fds);

	while(len>0) {
		int sret = select (fd+1, &fds, 0, 0, &timeout);
		if(sret==-1) { warne("select"); return false; }
		if(sret==0) { warn("read (%s): timeout expires with %d bytes unread",reason,len); return false; }

		int bytes = ::read (fd, buf, len);
		if(bytes<0) { warne("read"); return false; }
		if(bytes>0 && !dirty) { timeout = timeout_general; dirty=true; }
		len-=bytes;
		buf+=bytes;
	}
#else
	COMMTIMEOUTS cto;
	cto.ReadIntervalTimeout = TIMEOUT_READ_BYTE;
	cto.ReadTotalTimeoutMultiplier = TIMEOUT_AVG_READ_BYTE;
	cto.ReadTotalTimeoutConstant = timeout_ms;
	cto.WriteTotalTimeoutMultiplier = 0;
	cto.WriteTotalTimeoutConstant = 0;
	if (!SetCommTimeouts(fd,&cto)) { /*warnE("SetTimeoutProperties")*/; return false; }
	bool virgin=true;

	while(len>0) {
		DWORD got=0;
//		debug_print_buf (len);
		if(!ReadFile (fd, buf, len, &got, 0)) { warnE("ReadFile"); return false; }
		if(got==0) { 
			//warn("[%s]ReadFile: timeout expires with %d bytes unread", m_dev_file, len); 
			return false; }
//		debug_print_buf (buf, got);
		len-=got;
		buf+=got;
		
		if(virgin) {
			virgin=false;
			COMMTIMEOUTS cto;
			cto.ReadIntervalTimeout = TIMEOUT_READ_BYTE;
			cto.ReadTotalTimeoutMultiplier = TIMEOUT_AVG_READ_BYTE;
			cto.ReadTotalTimeoutConstant = TIMEOUT_AVG_READ_BYTE*5;
			cto.WriteTotalTimeoutMultiplier = 0;
			cto.WriteTotalTimeoutConstant = 0;
			if (!SetCommTimeouts(fd,&cto)) { /*warnE("SetTimeoutProperties")*/; return false; }
		}
	}
#endif
	assert(len==0);
	return true;
}

bool SerialPort::_write (const uchar *data, int len) const {
#ifdef LINUX
	ssize_t bytes = ::write (fd, data, len);
	if(bytes==len) return true;
	warne("write (%d out of %d bytes written)", bytes, len);
#else
	DWORD written=0;
//	debug_print_buf (data,len);
	if(WriteFile (fd, data, len, &written, 0) && len==written) return true;
	warnE("WriteFile (%d out of %d bytes written)", written, len);
#endif
	return false;
}

/*	The SICK laser scanners require at least 55us between each byte received.
	1byte/55us = 145kbps, thus only the 500kbps mode needs attention.
	SICK also requires the delay to be <6ms (>167Hz) to prevent timeout.

	About nanosleep():
	On Linux, the delay is restricted to time-slicing accuracy,
	which is >10ms for kernel 2.4 and >1~4ms for kernel 2.6.
	With kernel 2.4, nanosleep() can pause "by busy waiting with microsecond
	precision when called from a process scheduled under a real-time policy
	like SCHED_FIFO or SCHED_RR. Unfortunately this requires superuser privileges.

	"On Windows Me, nanosleep() has millisecond accuracy.
	On Windows NT/2000/XP/2003, waitable timers with 100-nanosecond resolution are available,
	so nanosleep() has 100-nanosecond accuracy."
*/
bool SerialPort::write (const uchar *data, int len) const       // TODO...
{
#ifdef LINUX
	if(write_delay.tv_nsec<0) return _write(data,len);
	for(int i=0; i<len; ++i) {
		if(!write(data[i])) { warne("tried to write %d bytes, but only written %d", len, i); return false; }
		nanosleep(&write_delay,0);
	}
	return true;
#else
	for(int i=0; i<len; ++i) {
		if(!write(data[i])) { 
			
			                                              //modified by linsm 2008 0323
			//CString str; 
			//str.Format("%d (%s)", i, m_dev_file); 
			//AfxMessageBox(str);                          modified by linsm  2008 0323
//	  boost::mutex::scoped_lock
//      lock(io_mutex);
			warne("tried to write %d bytes, but only written %d", len, i); return false; }
		shortnap();
	}
	return true;
#endif
}


// Cancel pending I/O operations. Actually operating Read/Write functions will finish with an error!(?)
bool SerialPort::cancel_io()
{
	if(!opened()) return true;
/*#ifdef LINUX
	...
#else
	return CancelIo(fd);
#endif*/
	return true;
}

