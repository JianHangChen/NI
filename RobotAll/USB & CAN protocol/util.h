/*
	Utility functions & macros

	Platform: Linux, Windows (Win32?)
	License: LGPL

	Written by Jim (r94922159@ntu.edu.tw), 2005-12-10. (tab=4)
	Modified by Jim, 2005-12-14: make it MSVC-compliant.
	Modified by Jim, 2006-2-20: time_str2(), class CircularQueue.
	Modified by Jim, 2006-3-20: ...
	Modified by Jim, 2006-8-11: #define VERBOSE, fixed bug (errors wouldn't show if BE_QUIET && !LINUX)
	Modified by Jim, 2006-8-14: time_str3(), system_pause()
	Modified by Jim, 2006-8-23: getch()...NOT FINISHED
	Modified by Jim, 2006-9-11: renamed CC1.. to _QCC1..
	Modified by Jimmy Yen, 2007-4-17: killed some compiler warnings

	Known Bugs:
	- Some functions defined in this file are not inline. Therefore you might only be able to include util.h only once in your project.
	- Some compilers (VC?) does not accept multiple definitions of inlined functions in a single project.
*/

#ifndef __INQ_UTIL2_H
#define __INQ_UTIL2_H

//#define VERBOSE			// print extra debug info
//#define BE_QUIET			// suppress info() and warn() outputs
#pragma warning(disable:4996)
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <errno.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#ifdef LINUX
	#include <unistd.h>

	// from ncurses -- have to add -lncurses
/*	namespace NCURSES {
		struct WINDOW {
			NCURSES_SIZE_T _cury, _curx, _maxy, _maxx, _begy, _begx;
			short _flags;
			attr_t _attrs;
			chtype _bkgd;
			bool _notimeout, _clear, _leaveok, _scroll, _idlok, _idcok, _immed, _sync, _use_keypad;
			int _delay;
			struct ldat *_line;
			NCURSES_SIZE_T _regtop, _regbottom;
			int _parx, _pary;
			WINDOW *_parent;
			struct pdat { NCURSES_SIZE_T _pad_y, _pad_x, _pad_top, _pad_left, _pad_bottom, _pad_right; } _pad;
			NCURSES_SIZE_T _yoffset;
			#ifdef _XOPEN_SOURCE_EXTENDED
			cchar_t _bkgrnd;
			#endif
		};
		extern int cbreak();
		extern int nocbreak();
		extern int wgetch(WINDOW*);
		extern NCURSES_EXPORT_VAR(WINDOW*) stdscr;

		// not from curses
		int getch(bool echo=false, bool wait=true)		// waits for a key, no echo
		{
			// assuming currently in normal modes
			cbreak();
			if(!echo) noecho();
			if(!wait) nodelay(stdscr,true);
			int k=wgetch(stdscr);
			if(!wait) nodelay(stdscr,false);
			if(!echo) echo();
			nocbreak();
			return k;
		}
	};
	using NCURSES::getch;
*/
#else
	#include <windows.h>
	#include <mmsystem.h>	// Any's suggestion
	#pragma comment( lib, "winmm.lib" )
	#define USE_WINMM		// use winmm.lib's timeGetTime() instead of time()
#endif

//---------------------------------------------------------------------------------------------------

typedef unsigned char byte;
typedef unsigned char uchar;
typedef unsigned short ushort;
typedef unsigned int uint;
typedef unsigned long ulong;

#define NELEMS(x) ((int)(sizeof(x)/sizeof(*(x))))

//---------------------------------------------------------------------------------------------------

inline static void system_pause()
{
	#ifdef LINUX
	printf("Press any key to continue . . . "); fflush(stdout);
/*	while(getch(false,false)!=ERR);
	getch();	*/
	getchar();
	#else
	system("PAUSE");
	#endif
}

//---------------------------------------------------------------------------------------------------
#ifdef __GNUC__
#include <sys/time.h>

// limitation: doesn't work with array of strings (why?)
// limitation: is not single statement
// limitation: array[] is not in local scope
#define foreach(array,i,arg1,args...)			\
	typeof(arg1) array[]={arg1,##args};			\
	for (int (i)=0; (i)<NELEMS(array); ++(i))
	
#define foreach2(array,i)						\
	for (int (i)=0; (i)<NELEMS(array); ++(i))

// limitation: N is evaluated many times.
#define foreach3(array,i,N)						\
	for (int (i)=0; (i)<(N); ++(i))

// limitation: end is evaluated many times.
#define foreach4(array,i,start,end)				\
	for (int (i)=(start); (i)<(end); ++(i))

#endif
//---------------------------------------------------------------------------------------------------

// Get current time in seconds (since 1970), with up to microsecond precision.
inline static double current_time() {
	#ifdef LINUX
		timeval tv;
		gettimeofday (&tv,0);
		return tv.tv_sec + (double)tv.tv_usec/1e6;	// up to us precision
	#else
		#ifdef USE_WINMM
		DWORD tm = timeGetTime();
		return tm/1000.0;				// up to ms precision
		#else
		return time(0);					// second precision only!
		#endif
	#endif
}

#if defined(USE_WINMM) || defined(LINUX)
	#define for_n_seconds(t)	for(double _fns_alarm_=current_time()+t; current_time()<_fns_alarm_;)
#else
	#define for_n_seconds(t)	for(double _fns_alarm_=current_time()+t+1; current_time()<_fns_alarm_;)
#endif

// output a constant-sized time string
inline static const char *time_str() {
	static char buf[32];					// not thread-safe
	time_t the_time_t = time(0);
	strftime (buf, sizeof(buf), "%Y-%m-%d %H:%M:%S", localtime (&the_time_t));
	return buf;
}

// output a constant-sized time string, for use in file names
inline static const char *time_str2() {
	static char buf[32];					// not thread-safe
	time_t the_time_t = time(0);
	strftime (buf, sizeof(buf), "%Y%m%d_%H%M%S", localtime (&the_time_t));
	return buf;
}

// output a constant-sized time string, be short and sufficient
inline static const char *time_str3() {
	static char buf[32];					// not thread-safe
	time_t the_time_t = time(0);
	strftime (buf, sizeof(buf), "%m%d %H%M%S", localtime (&the_time_t));
	return buf;
}

inline static void msleep (const uint ms) {
	#ifdef LINUX
		usleep (ms*1000);
	#else
		Sleep(ms);
	#endif
}

//---------------------------------------------------------------------------------------------------
// Message/error-handling (eg. outputting messages..)
//

#ifdef VERBOSE
	#define verbose(x) {x;}
#else
	#define verbose(x) {(void)0;}
#endif

#ifdef BE_QUIET
	inline static int _printf (const char *fmt, ...) { return rand(); }
	inline static int _vprintf (const char *fmt, va_list ap) { return rand(); }
#else
	#define _printf printf
	#define _vprintf vprintf
#endif

#ifdef LINUX
	#define _QCC0		"\033[0m"
	#define _QCC1		"\033[0;44;37m"
	#define _QCC2i		"\033[1;40;37m"
	#define _QCC2w		"\033[1;40;31m"
	#define _QCC2e		"\033[1;40;31m"
	#define _QCC3		"\033[0;40;37m"
	#define _QCC4		"\033[1;40;33m"

#else
	#define _QCC0
	#define _QCC1
	#define _QCC2i
	#define _QCC2w
	#define _QCC2e
	#define _QCC3
	#define _QCC4

	inline static char *StrError() {
		static char msg[10240];
		assert (FormatMessage (FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS,
			NULL, GetLastError(), 0, (LPTSTR)&msg, sizeof(msg), NULL)!=0);
		return msg;
	}
#endif

#ifdef __GNUC__
	// INFO: fmt must be a string literal
	inline static int _die(int n) { exit(-1); return n; }
	#define info(fmt, args...)	( _printf(_QCC1"[%s]"_QCC3" %-33.33s "_QCC2i" INFO: "_QCC4 fmt _QCC0"\n", time_str3(), __PRETTY_FUNCTION__, ## args) )
	#define infoe(fmt, args...)	( _printf(_QCC1"[%s]"_QCC3" %-33.33s "_QCC2i" INFO: "_QCC4 fmt ": %s"_QCC0"\n", time_str3(), __PRETTY_FUNCTION__, ## args, strerror(errno)) )
	#define infoE(fmt, args...)	( _printf(_QCC1"[%s]"_QCC3" %-33.33s "_QCC2i" INFO: "_QCC4 fmt ": %s"_QCC0"\n", time_str3(), __PRETTY_FUNCTION__, ## args, StrError()) )
	#define warn(fmt, args...)	( _printf(_QCC1"[%s]"_QCC3" %-33.33s "_QCC2w" WARN: "_QCC4 fmt _QCC0"\n", time_str3(), __PRETTY_FUNCTION__, ## args) )
	#define warne(fmt, args...)	( _printf(_QCC1"[%s]"_QCC3" %-33.33s "_QCC2w" WARN: "_QCC4 fmt ": %s"_QCC0"\n", time_str3(), __PRETTY_FUNCTION__, ## args, strerror(errno)) )
	#define warnE(fmt, args...)	( _printf(_QCC1"[%s]"_QCC3" %-33.33s "_QCC2w" WARN: "_QCC4 fmt ": %s"_QCC0"\n", time_str3(), __PRETTY_FUNCTION__, ## args, StrError()) )
	#define die(fmt, args...)	( _die (printf(_QCC1"[%s]"_QCC2e"  ERROR "_QCC3"(in %s, line %d, %s): "_QCC4 fmt _QCC0"\n", time_str3(), __FILE__, __LINE__, __PRETTY_FUNCTION__, ## args)) )
	#define diee(fmt, args...)	( _die (printf(_QCC1"[%s]"_QCC2e"  ERROR "_QCC3"(in %s, line %d, %s): "_QCC4 fmt ": %s"_QCC0"\n", time_str3(), __FILE__, __LINE__, __PRETTY_FUNCTION__, ## args, strerror(errno))) )
	#define dieE(fmt, args...)	( _die (printf(_QCC1"[%s]"_QCC2e"  ERROR "_QCC3"(in %s, line %d, %s): "_QCC4 fmt ": %s"_QCC0"\n", time_str3(), __FILE__, __LINE__, __PRETTY_FUNCTION__, ## args, StrError())) )
#else
	inline static int _die(int n) { exit(-1); return n; }
	#define _ooxx_(ss,SS)	inline static int ss(const char *fmt, ...) {           \
								va_list ap;                           \
								va_start(ap,fmt);                     \
								_printf("[%s]  "SS": ", time_str3()); \
								_vprintf(fmt,ap);                     \
								va_end(ap);
	#define __xxoo_      		_printf("\n")
	#define __xxooe_			_printf(": %s\n", strerror(errno))
	#define __xxooE_         	_printf(": %s\n", StrError())
	#define _xxoo_      		printf("\n")
	#define _xxooe_				printf(": %s\n", strerror(errno))
	#define _xxooE_         	printf(": %s\n", StrError())
	_ooxx_(info, "INFO") return __xxoo_; }
	_ooxx_(infoe,"INFO") return __xxooe_; }
	_ooxx_(infoE,"INFO") return __xxooE_; }
	_ooxx_(warn, "WARN") return __xxoo_; }
	_ooxx_(warne,"WARN") return __xxooe_; }
	_ooxx_(warnE,"WARN") return __xxooE_; }
	_ooxx_(die, "ERROR") return _die(_xxoo_); }
	_ooxx_(diee,"ERROR") return _die(_xxooe_); }
	_ooxx_(dieE,"ERROR") return _die(_xxooE_); }
#endif

//---------------------------------------------------------------------------------------------------
/*
template <class T> inline static T sqr(T x) { return x*x; }

template <class D, class T>
inline void limit_range (T& i, const D low, const D high) {
	if (i<low) i=low; else
	if (i>high) i=high;
}

inline static float  rand_01f () { return (float) random()/RAND_MAX; }
inline static double rand_01d () { return (double)random()/RAND_MAX; }
template <class T> inline static T rand_01() { return (T)random()/RAND_MAX; }

inline static float  rand_11f () { return (float) random()/((RAND_MAX+1)/2)-1; }
inline static double rand_11d () { return (double)random()/((RAND_MAX+1)/2)-1; }
template <class T> inline static T rand_11() { return (T)random()/((RAND_MAX+1)/2)-1; }

template <class T> inline static int serial_lt (const T& a, const T& b, const T& c)             { return (a<b && b<c); }
template <class T> inline static int serial_lt (const T& a, const T& b, const T& c, const T& d) { return (a<b && b<c && c<d); }

template <class T> inline static T deg2rad (T deg) { return deg*(M_PI/180.0); }
template <class T> inline static T rad2deg (T rad) { return rad*(180.0/M_PI); }
*/

//---------------------------------------------------------------------------------------------------

/*
inline void debug_print_buf(int len) {
	printf ("(DEB=%d) ", len);
	fflush(stdout);
}
inline void debug_print_buf(const unsigned char *data, int len) {
	printf ("(DEB");
	for(int i=0; i<len; ++i)
		printf (" %02X", data[i]);
	printf (") ");
	fflush(stdout);
}
*/

//---------------------------------------------------------------------------------------------------
// Fixed-capacity circular queue
// Example usage:
//		CircularQueue<int,10> Q;
//		Q.push(1), Q.push(2), Q.push(3);
//		printf ("%d %d %d", Q.pop(), Q.pop(), Q.pop());
//
template <class T, int CAP>
class CircularQueue {
	T buf[CAP];
	int nw, nr, _count;		// next write index, next read index, current elements count
							// Indexes iterates backwards in buf[].
	public:
	CircularQueue() : nw(0), nr(0), _count(0) {}
	~CircularQueue() {}
	int count() const { return _count; }
	int capacity() const { return CAP; }
	void clear() { nw=nr=0; count=0; }		// empty the queue
	void push (T& x) {
		if(_count==CAP) { warn("buffer overflow! (new data discarded)"); return; }
		buf[nw]=x;
		if(--nw<0) nw=CAP-1;
		++_count;
	}
	T& pop() {		// NOTE: Keeping the returned (dangling) pointer is not thread-safe.
		if(_count==0) { warn("buffer underflow! (bogus data returned)"); return buf[nr]; }
		T &ret = buf[nr];
		if(--nr<0) nr=CAP-1;
		--_count;
		return ret;
	}
};


//---------------------------------------------------------------------------------------------------
// Error message system
// (not thread-safe)
//
/*
#define EMS_ENTER(x)	EMS_Class _ems_var_(x);


class EMS_Class {
	static vector<const char *> names;
	const char *name;
	int ID;

	public:
	EMS_Class(const char *_name) {
		names.push_back(name=_name);
		ID=names.size();
	}
	~EMS_Class() {
		while(ID<names.size()) { fprintf("EMS ERROR: popping unpopped name -- %s\n",names.back()); names.pop_back(); }
		names.pop_back();
	}
};




*/
//---------------------------------------------------------------------------------------------------
#endif
