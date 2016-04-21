#pragma once

#include <api/inc/fmod.hpp>
#include <api/inc/fmod_errors.h>

#include <vector>
#include "string"
#include <stdio.h>
#include <conio.h>
#include "atlstr.h"
#include <windows.h>
#include "fstream"
#include <atlbase.h>

extern CComModule _Module;

#include <atlcom.h>

#include "sapi.h"
#include <sphelper.h>


using namespace std;

typedef struct Voice_v
{
	char* word_v;
	int v;
}Voice_v;


class CSoundPlay
{
public:
	CSoundPlay(void);
	~CSoundPlay(void);
	FMOD::System* system;
	std::vector<FMOD::Sound*> m_sound_closet;
	//COM元件
	CComPtr	<ISpVoice>	cpVoice;

	FMOD::Sound *Sound1;
	std::vector <FMOD::Sound*> SignLanq;
	//FMOD::Sound *Sound2;
	//FMOD::Sound *Sound3;
	//FMOD::Sound *Sound4;

protected:
	void ERRCHECK(FMOD_RESULT result)
	{
		if (result != FMOD_OK)
		{
			printf("FMOD error! (%d) %s\n", result, FMOD_ErrorString(result));
			exit(-1);
		}
	}
public:
	bool Initialize(void);
	bool LoadSoundFile(const char* FileName);
	bool PlaySound(int id);
	bool PlaySoundByFile(FMOD::Sound *S1);
	FMOD::Channel* channel;
	bool ReleaseSoundFile(FMOD::Sound *S1);
	bool speak(char* YourInputText);
	bool speak(char* YourInputText, int iSpeed);
public:
	void VolumnLower(void);
public:
	void VolumnHigher(void);

	//工研院說話核心引擎函式
	int TTSEngineTest(wchar_t *pszInputString, wchar_t *szVoiceEngine, int par_speed, int iMode = 0, wchar_t *szOutputFileName = NULL, SPSTREAMFORMAT iOutputFileFormat = SPSF_16kHz16BitMono);
	
	void ReadSignLanq(int numb);
};
