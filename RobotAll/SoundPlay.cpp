#include "StdAfx.h"
#include "SoundPlay.h"
#include "iostream"
#include "stdlib.h"
using namespace std;


CSoundPlay::CSoundPlay(void)
: system(NULL)
, channel(NULL)
{
	Initialize();
}

CSoundPlay::~CSoundPlay(void)
{
	//ReleaseSoundFile();
}

bool CSoundPlay::Initialize(void)
{
	FMOD_RESULT result;
	unsigned int      version;

	/*
        Create a System object and initialize.
    */
    result = FMOD::System_Create(&system);
    ERRCHECK(result);

    result = system->getVersion(&version);
    ERRCHECK(result);

   /* if (version < FMOD_VERSION)
    {
        printf("Error!  You are using an old version of FMOD %08x.  This program requires %08x\n", version, FMOD_VERSION);
        return 0;
    }*/

    result = system->init(32, FMOD_INIT_NORMAL, 0);
    ERRCHECK(result);

	// channel
	channel = 0;

	// load sound file

	//sprintf(FileName,"./Output4.wav");

	//LoadSoundFile(FileName);
	//LoadSoundFile("./Output4.wav",speaker.Sound1);

	return true;
}

bool CSoundPlay::LoadSoundFile(const char* FileName)
{
	// loading the sound files

	FMOD_RESULT       result;

	result = system->createSound(FileName, FMOD_HARDWARE, 0, &Sound1);
	//result = system->createSound(FileName, FMOD_HARDWARE, 0, &Sound1);

    ERRCHECK(result);

	return true;
}



bool CSoundPlay::PlaySoundByFile(FMOD::Sound *S1)
{
	FMOD_RESULT result;
	result = system->playSound(FMOD_CHANNEL_FREE, S1, false, &channel);
    ERRCHECK(result);
	return false;
}


bool CSoundPlay::ReleaseSoundFile(FMOD::Sound *S1)
{
	FMOD_RESULT result;
	
		result = S1->release();
		 ERRCHECK(result);


	return true;
}

//20111018工研院的程式//

int CSoundPlay::TTSEngineTest(wchar_t *pszInputString, wchar_t *szVoiceEngine, int par_speed, int iMode, wchar_t *szOutputFileName, SPSTREAMFORMAT iOutputFileFormat)
{
	HRESULT hr;

	CComPtr<IEnumSpObjectTokens>        cpEnum;
	//設定 VoiceEngine
	CComPtr <ISpObjectToken> cpVoiceToken;
	CComPtr <ISpObjectToken> cpAudioToken;

	SpEnumTokens( SPCAT_AUDIOOUT, NULL, NULL, &cpEnum );


	hr = SpFindBestToken(SPCAT_VOICES, szVoiceEngine, NULL, &cpVoiceToken);									
	if (!SUCCEEDED(hr)) //設定 VoiceEngine 失敗
	{
		printf("Error(SpFindBestToken): VoiceEngine Setting Failed...\n");
		return FALSE;
	}
	cpVoice->SetVoice(cpVoiceToken);

	//設定 WaveFormat
	CSpStreamFormat cAudioFormat;
	hr = cAudioFormat.AssignFormat(iOutputFileFormat);
	if(!SUCCEEDED(hr))
	{
		printf("Error(AssignFormat): AssignFormat Failed...\n");
		return FALSE;
	}
	
	//兩種聲音輸出方式
	if (iMode==0) //iMode==0 將聲音輸出至喇叭
	{
		//set the output to the audio device
		if( SUCCEEDED( hr ) )
		{
			cpEnum->Item(0, &cpAudioToken);
			cpVoice->SetOutput(cpAudioToken, TRUE);
			//cpVoice->Speak(pszInputString, SPF_ASYNC | SPF_IS_XML, NULL);//從SPF_ASYNC改成SPF_DEFAULT就可以用了，原因不明

			long speed;
			cpVoice->GetRate(&speed);
			speed += par_speed;
			cpVoice->SetRate(speed);

			cpVoice->Speak(pszInputString, SPF_DEFAULT | SPF_IS_XML, NULL);
		}
		else
		{
			printf("ERROR(SetFormat)\n"); 
			return FALSE;
		}
	}
	else //iMode!=0 將聲音輸出至Wave檔案, 檔案格式為iOutputFileFormat
	{	
		CComPtr	<ISpStream>	cpWaveStream;
		
		hr = SPBindToFile(szOutputFileName, SPFM_CREATE_ALWAYS, &cpWaveStream, &cAudioFormat.FormatId(), cAudioFormat.WaveFormatExPtr() );
		
		if( SUCCEEDED( hr ) )  //假如SPBindToFile is OK.
		{
			//set the output to cpWaveStream so that the output audio data will be stored in cpWaveStream
			cpVoice->SetOutput(cpWaveStream, TRUE);	
	
			long speed;
			cpVoice->GetRate(&speed);
			speed += par_speed;
			cpVoice->SetRate(speed);

			cpVoice->Speak(pszInputString, SPF_DEFAULT | SPF_IS_XML, NULL);				
			if ( cpWaveStream )
			{
				cpWaveStream->Close();
				cpWaveStream.Release();
			}													
		}		
		else
		{
			printf("ERROR(SPBindToFile)\n"); 
			return FALSE;
		}// end of if( SUCCEEDED( hr ) )  //假如SPBindToFile is OK.		
	}//end of if (iMode==0) //iMode==0 將聲音輸出至喇叭	

	//system("Pause");
	//cpVoiceToken.Release();
	//system("pause");
	return TRUE;
}


bool CSoundPlay::speak(char* YourInputText)
{
	/*為了轉換 char* 做 System 內使用的 wchar (for converting the C-Style char to System Char)***/

	CString YourInputText_CString;

	YourInputText_CString.Format("%s",YourInputText);

	/*文字轉換作語音的程序 Text-To-Speech procedures*/

	//ISpVoice * pVoice = NULL;  /*1. 準備一支 NULL了的 pointercreate a pointer*/

	if (FAILED(CoInitialize(NULL))) /* 2. 嘗試初始化 COM Check if we can Initialize the COM*/

	{

	printf("Error to intiliaze COM"); return false;

	}

	 /*3. 初始化 Text-To-Speech object (Initialize Text-To-Speech object)*/

	//HRESULT hr = CoCreateInstance(CLSID_SpVoice, NULL, CLSCTX_ALL, IID_ISpVoice, (void **)&pVoice);
	HRESULT hr;
		// Initialize SAPI 
	printf("TTS SAPI Initializing ...... ");
	if (FAILED(::CoInitialize(NULL)))
	{ 
		printf("ERROR to initialize SAPI\n"); 
		return FALSE;
	}		
	printf("Initialized OK!\n\n");
	//Create Instance for cpVoice
	hr = cpVoice.CoCreateInstance( CLSID_SpVoice );	
	if (!SUCCEEDED(hr)) //假如CoCreateInstance is not OK.
	{ 
		printf("ERROR(CoCreateInstance)\n"); 
		return FALSE;
	}

	if( SUCCEEDED( hr ) ) /*4. 若初始化 TTS engine 成功 (If initialization of TTS succeeded)*/

	{

		/*5. 讀出 "Your Input is" (speaks the words "Your Input is") */

		//hr = pVoice->Speak(L"Your Input is", 0, NULL);

		/* 5. 讀出 YourInputText_CString 的內容 (speaks the content of ) YourInputText_CString */
		
		//debug use
		 //cout<<YourInputText_CString<<"\n\n";
		//hr = pVoice->Speak(YourInputText_CString.AllocSysString(),SPF_ASYNC | SPF_IS_XML,NULL);
			//設定pTTS合成參數
			//wchar_t szVoiceEngine[] = L"Name=ITRI Theresa pTTS (PMC)"; // ITRI Female Voice Engine
			wchar_t szVoiceEngine[] = L"Name=ITRI Bruce pTTS (PMC)"; // ITRI Male Voice Engine: "Name=ITRI Bruce pTTS (PMC)"
																	 // ITRI Female Voice Engine: "Name=ITRI Theresa pTTS (PMC)"
			//pTTS 準位 - 合理範圍:-10到10之間的整數值, 值越大則音高越高；反之則音高越低
			int iLevel = 3;
			//pTTS 方向 - 合理範圍:0到2之間的整數值 (0：正常, 1：像機器人, 2：像外國人說中文)
			int iDirection = 0;
			//pTTS 方向 - 合理範圍:0到20之間的整數值, 值越大則抑揚頓挫越明顯；反之則越趨平版
			int iVariance = 8;
			//將音源輸出至喇叭(iMode = 0)或輸出聲音檔案(iMode = 1)
			int iMode = 0;	
			//輸出Wave檔案格式: 可查詢 Microsoft Speech SDK 文件說明 (help SPSTREAMFORMAT)
			SPSTREAMFORMAT iOutputFileFormat = SPSF_16kHz16BitMono;
			//SPSTREAMFORMAT iOutputFileFormat = SPSF_8kHz16BitMono;
			//輸出Wave檔案名稱
			//wchar_t szOutputFileName[] = L"自我介紹.wav";	
			//wchar_t szOutputFileName[] = L"最開始的歡迎詞.wav";
			//wchar_t szOutputFileName[] = L"伸展身體時候的自我介紹.wav";
			wchar_t szOutputFileName[] = L"機器人出場.wav";
		wchar_t *pszInputString = YourInputText_CString.AllocSysString();	// 輸入的文字串
		wchar_t *pszTagOutputString;				// 輸出的插入標記語言的文字串
		int		iTagAdditionLen = 64;				// 插入標記語言所新增的字串長度

		pszTagOutputString = new wchar_t [wcslen(pszInputString)+iTagAdditionLen];

		swprintf(pszTagOutputString, L"<PitchModify para=\"%d;%d;%d\"/>%s", iLevel, iDirection, iVariance, pszInputString);

		TTSEngineTest(pszTagOutputString, szVoiceEngine,3, iMode, szOutputFileName, iOutputFileFormat);

		delete [] pszTagOutputString;
		

		 // pVoice->Release(); //6. 釋放TTS Engine (release the Text2Speech Engine after using)

		//pVoice = NULL;//7. 在釋放前 set NULL(Set to NULL before Uninitialize the COM)

	}

	CoUninitialize(); // 8. 解除 COM (Uninitialize the COM)

	return true;
}
bool CSoundPlay::speak(char* YourInputText,int iSpeed)
{
	/*為了轉換 char* 做 System 內使用的 wchar (for converting the C-Style char to System Char)***/

	CString YourInputText_CString;

	YourInputText_CString.Format("%s",YourInputText);

	/*文字轉換作語音的程序 Text-To-Speech procedures*/

	//ISpVoice * pVoice = NULL;  /*1. 準備一支 NULL了的 pointercreate a pointer*/

	if (FAILED(CoInitialize(NULL))) /* 2. 嘗試初始化 COM Check if we can Initialize the COM*/

	{

	printf("Error to intiliaze COM"); return false;

	}

	 /*3. 初始化 Text-To-Speech object (Initialize Text-To-Speech object)*/

	//HRESULT hr = CoCreateInstance(CLSID_SpVoice, NULL, CLSCTX_ALL, IID_ISpVoice, (void **)&pVoice);
	HRESULT hr;
		// Initialize SAPI 
	printf("TTS SAPI Initializing ...... ");
	if (FAILED(::CoInitialize(NULL)))
	{ 
		printf("ERROR to initialize SAPI\n"); 
		return FALSE;
	}		
	printf("Initialized OK!\n\n");
	//Create Instance for cpVoice
	hr = cpVoice.CoCreateInstance( CLSID_SpVoice );	
	if (!SUCCEEDED(hr)) //假如CoCreateInstance is not OK.
	{ 
		printf("ERROR(CoCreateInstance)\n"); 
		return FALSE;
	}

	if( SUCCEEDED( hr ) ) /*4. 若初始化 TTS engine 成功 (If initialization of TTS succeeded)*/

	{

		/*5. 讀出 "Your Input is" (speaks the words "Your Input is") */

		//hr = pVoice->Speak(L"Your Input is", 0, NULL);

		/* 5. 讀出 YourInputText_CString 的內容 (speaks the content of ) YourInputText_CString */
		
		//debug use
		 //cout<<YourInputText_CString<<"\n\n";
		//hr = pVoice->Speak(YourInputText_CString.AllocSysString(),SPF_ASYNC | SPF_IS_XML,NULL);
			//設定pTTS合成參數
			//wchar_t szVoiceEngine[] = L"Name=ITRI Theresa pTTS (PMC)"; // ITRI Female Voice Engine
			wchar_t szVoiceEngine[] = L"Name=ITRI Bruce pTTS (PMC)"; // ITRI Male Voice Engine: "Name=ITRI Bruce pTTS (PMC)"
																	 // ITRI Female Voice Engine: "Name=ITRI Theresa pTTS (PMC)"
			//pTTS 準位 - 合理範圍:-10到10之間的整數值, 值越大則音高越高；反之則音高越低
			int iLevel = 3;
			//pTTS 方向 - 合理範圍:0到2之間的整數值 (0：正常, 1：像機器人, 2：像外國人說中文)
			int iDirection = 0;
			//pTTS 方向 - 合理範圍:0到20之間的整數值, 值越大則抑揚頓挫越明顯；反之則越趨平版
			int iVariance = 8;
			//將音源輸出至喇叭(iMode = 0)或輸出聲音檔案(iMode = 1)
			int iMode = 1;	
			//輸出Wave檔案格式: 可查詢 Microsoft Speech SDK 文件說明 (help SPSTREAMFORMAT)
			SPSTREAMFORMAT iOutputFileFormat = SPSF_16kHz16BitMono;
			//SPSTREAMFORMAT iOutputFileFormat = SPSF_8kHz16BitMono;
			//輸出Wave檔案名稱
			//wchar_t szOutputFileName[] = L"自我介紹.wav";	
			//wchar_t szOutputFileName[] = L"最開始的歡迎詞.wav";
			//wchar_t szOutputFileName[] = L"伸展身體時候的自我介紹.wav";
			wchar_t szOutputFileName[] = L"N07.wav";
		wchar_t *pszInputString = YourInputText_CString.AllocSysString();	// 輸入的文字串
		wchar_t *pszTagOutputString;				// 輸出的插入標記語言的文字串
		int		iTagAdditionLen = 64;				// 插入標記語言所新增的字串長度

		pszTagOutputString = new wchar_t [wcslen(pszInputString)+iTagAdditionLen];

		swprintf(pszTagOutputString, L"<PitchModify para=\"%d;%d;%d\"/>%s", iLevel, iDirection, iVariance, pszInputString);

		TTSEngineTest(pszTagOutputString, szVoiceEngine,iSpeed, iMode, szOutputFileName, iOutputFileFormat);

		delete [] pszTagOutputString;
		

		 // pVoice->Release(); //6. 釋放TTS Engine (release the Text2Speech Engine after using)

		//pVoice = NULL;//7. 在釋放前 set NULL(Set to NULL before Uninitialize the COM)

	}

	CoUninitialize(); // 8. 解除 COM (Uninitialize the COM)

	return true;
}
void CSoundPlay::VolumnLower(void)
{
	float volumn;
	channel->getVolume(&volumn);

	volumn -= 0.8f;
	channel->setVolume(volumn);	
}

void CSoundPlay::VolumnHigher(void)
{
	float volumn;
	channel->getVolume(&volumn);

	volumn += 0.8f;
	channel->setVolume(volumn);
}
void CSoundPlay::ReadSignLanq(int numb)
{
char CharBuf [50];

for(int i=1;i<numb+1;i++)
	{
	  if(i<10)
	  sprintf (CharBuf, "SoundFiles/H0%d.wav", i);
	  else
	  sprintf (CharBuf, "SoundFiles/H%d.wav", i);
	
	  LoadSoundFile(CharBuf);

	  SignLanq.push_back(Sound1);
	}

}