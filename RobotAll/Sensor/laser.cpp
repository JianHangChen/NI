#include "StdAfx.h"
#include "laser.h"


laser::laser(void)
{
	j = 0;
	ya = new int [769];
	x_data =  new float [681];
	y_data =  new float [681];
	pi = 3.1415926;
	sample = (360.0/1024.0)/180.0*pi;
	urg = new UrgLaser();
	std::string _portStr= std::string( "\\\\.\\COM83" );
	urg->Open(TEXT( _portStr.c_str() ));
	for(int i=0;i<681;i++)
	{
	x_data[i]=0;
	y_data[i]=0;
	}
}

void laser::transfer(void)
{
		urg->Scan769(ya) ;  //代表我把掃描到的769筆data都存在一個ya的矩陣中
		//Sleep(0);
	   for(n = 43; n <724;n++)
	   {
		   if (ya[n]<=20)
			   ya[n] = 20;
		   if (ya[n]>=4096)
			   ya[n] = 4096;

		   x_data[n-43] = ya[n] * cos((n-384)*sample);
		   y_data[n-43] = ya[n] * sin((n-384)*sample);
		   //y_data[340] = 0.0;
	   }
	

}

laser::~laser(void)
{
	delete ya;
	delete x_data;
	delete y_data;
	delete urg;
}
