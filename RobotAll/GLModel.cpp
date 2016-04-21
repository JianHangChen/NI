#include "StdAfx.h"
#include "GLModel.h"
#include "Kine.h"

#include "model_leg\center_100.h"
#include "model_leg\LLeg01_100.h"
#include "model_leg\LLeg02_100.h"
#include "model_leg\LLeg03_100.h"
#include "model_leg\LLeg04_100.h"
#include "model_leg\LLeg05_100.h"
#include "model_leg\LLeg06_100.h"
#include "model_leg\RLeg01_100.h"
#include "model_leg\RLeg02_100.h"
#include "model_leg\RLeg03_100.h"
#include "model_leg\RLeg04_100.h"
#include "model_leg\RLeg05_100.h"
#include "model_leg\RLeg06_100.h"
#include "model_ARM\waist01_100.h"
#include "model_ARM\waist02_100.h"
#include "model_ARM\RArm01_100.h"
#include "model_ARM\RArm02_100.h"
#include "model_ARM\RArm03_100.h"
#include "model_ARM\RArm04_100.h"
#include "model_ARM\RArm05_100.h"
#include "model_ARM\RArm06_100.h"
#include "model_ARM\LArm01_100.h"
#include "model_ARM\LArm02_100.h"
#include "model_ARM\LArm03_100.h"
#include "model_ARM\LArm04_100.h"
#include "model_ARM\LArm05_100.h"
#include "model_ARM\LArm06_100.h"
extern Kine gKineAll; // 機器人全身 kinematics trains 宣告
extern double gLAngZWorld ; // Angle of the LLeg in  world coordinates
extern double gRAngZWorld ; // Angle of the RLeg in  world coordinates
#define PI 3.1415926

GLModel::GLModel(void)
{
	//gRadRatio = 180/3.1415926;
	//theta_home_LLeg[6]/*{ -90.0 , 0.0 , -90.0 , 0.0 , 0.0 , 0.0 }*/;//1-6
	//theta_home_RLeg[6]={ -90.0 , 0.0 , -90.0 , 0.0 , 0.0 , 0.0 };
	/*theta_home_LArm[0]=-90.0/gRadRatio;
	theta_home_LArm[1]=90.0/gRadRatio;
	theta_home_LArm[2]=90.0/gRadRatio;
	theta_home_LArm[3]=90.0/gRadRatio;
	theta_home_LArm[4]=0.0/gRadRatio;
	theta_home_LArm[5]=0.0/gRadRatio;
	theta_home_LArm[6]=90.0/gRadRatio;
	theta_home_RArm[0]=-90.0/gRadRatio;
	theta_home_RArm[1]=90.0/gRadRatio;
	theta_home_RArm[2]=90.0/gRadRatio;
	theta_home_RArm[3]=90.0/gRadRatio;
	theta_home_RArm[4]=0.0/gRadRatio;
	theta_home_RArm[5]=0.0/gRadRatio;
	theta_home_RArm[6]=90.0/gRadRatio;
	theta_home_waist[0]=180.0/gRadRatio;
	theta_home_waist[1]=0.0/gRadRatio;
	theta_home_waist[2]=90.0/gRadRatio;*/
	//theta_home_LArm[7]={ -90.0 , 90.0 , 90.0 , 90.0 , 0.0 , 0.0 , 90.0 };//4-9 10=endeffector
	//theta_home_RArm[7]={ -90.0 , 90.0 , 90.0 , 90.0 , 0.0 , 0.0 , 90.0};
	//theta_home_waist[3]={ 180.0 , 0.0 , 90.0 };//1-3

	//double d_LLeg[6]={ -99.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 };//1-6
	//double d_RLeg[6]={ -99.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 };
	/*d_LArm[0]=0.0;
	d_LArm[1]=245.5;
	d_LArm[2]=0.0;
	d_LArm[3]=267.0;
	d_LArm[4]=0.0;
	d_LArm[5]=208.0;
	d_LArm[6]=0.0;
	d_RArm[0]=0.0;
	d_RArm[1]=-245.5;
	d_RArm[2]=0.0;
	d_RArm[3]=-267.0;
	d_RArm[4]=0.0;
	d_RArm[5]=-208.0;
	d_RArm[6]=0.0;
	d_waist[0]=99.0;
	d_waist[1]=0.0;
	d_waist[2]=0.0;*/
	//d_LArm[7]={ 0.0 , 245.5 , 0.0 , 267.0 , 0.0 , 208.0 , 0.0 };//4-9 10=endeffector
	//d_RArm[7]={ 0.0 , -245.5 , 0.0, -267.0 , 0.0 , -208.0 , 0.0 };
	//d_waist[3]={ 99.0 , 0.0 , 0.0 };//1-3

	//double a_LLeg[6]={ -80.0 , 0.0 , 0.0 , -240.0 , -240.0 , 0.0 };//1-6
	//double a_RLeg[6]={ 80.0 , 0.0 , 0.0 , -240.0 , -240.0 , 0.0 };
	/*a_LArm[0]=0.0;
	a_LArm[1]=0.0;
	a_LArm[2]=0.0;
	a_LArm[3]=-10.0;
	a_LArm[4]=10.0;
	a_LArm[5]=0.0;
	a_LArm[6]=156.0;
	a_RArm[0]=0.0;
	a_RArm[1]=0.0;
	a_RArm[2]=0.0;
	a_RArm[3]=-10.0;
	a_RArm[4]=10.0;
	a_RArm[5]=0.0;
	a_RArm[6]=156.0;
	a_waist[0]=0.0;
	a_waist[1]=0.0;
	a_waist[2]=271.75;*/
	//a_LArm[7]={ 0.0 , 0.0 , 0.0 , -10.0 , 10.0 , 0.0 , 156.0 };//4-9 10=endeffector
	//a_RArm[7]={ 0.0 , 0.0 , 0.0 , -10.0 , 10.0 , 0.0 , 156.0 };
	//a_waist[3]={ 0.0 , 0.0 , 271.75 };//1-3

	//double alpha_LLeg[6]={ 0.0 , -90.0 , 90.0 , 0.0 , 0.0 , -90.0 };//1-6
	//double alpha_RLeg[6]={ 0.0 , -90.0 , -90.0 , 0.0 , 0.0 , 90.0 };
	/*alpha_LArm[0]=0.0/gRadRatio;
	alpha_LArm[1]=-90.0/gRadRatio;
	alpha_LArm[2]=-90.0/gRadRatio;
	alpha_LArm[3]=90.0/gRadRatio;
	alpha_LArm[4]=-90.0/gRadRatio;
	alpha_LArm[5]=90.0/gRadRatio;
	alpha_LArm[6]=0.0/gRadRatio;
	alpha_RArm[0]=0.0/gRadRatio;
	alpha_RArm[1]=90.0/gRadRatio;
	alpha_RArm[2]=90.0/gRadRatio;
	alpha_RArm[3]=-90.0/gRadRatio;
	alpha_RArm[4]=90.0/gRadRatio;
	alpha_RArm[5]=-90.0/gRadRatio;
	alpha_RArm[6]=0.0/gRadRatio;
	alpha_waist[0]=0.0/gRadRatio;
	alpha_waist[1]=90.0/gRadRatio;
	alpha_waist[2]=0.0/gRadRatio;*/
	//alpha_LArm[7]={ 0.0 , -90.0 , -90.0 , 90.0 , -90.0 , 90.0 , 0.0 };//4-9 10=endeffector
	//alpha_RArm[7]={ 0.0 , 90.0 , 90.0 , -90.0 , 90.0 , -90.0 , 0.0 };
	//alpha_waist[3]={ 0.0 , 90.0 , 0.0 };//1-3
}


GLModel::~GLModel(void)
{
}
void GLModel::InitGLModel(void)
{
	link_center=center::Gen3DObjectList();
	link_LLeg[0]=LLeg01::Gen3DObjectList();
	link_LLeg[1]=LLeg02::Gen3DObjectList();
	link_LLeg[2]=LLeg03::Gen3DObjectList();
	link_LLeg[3]=LLeg04::Gen3DObjectList();
	link_LLeg[4]=LLeg05::Gen3DObjectList();
	link_LLeg[5]=LLeg06::Gen3DObjectList();
	link_RLeg[0]=RLeg01::Gen3DObjectList();
	link_RLeg[1]=RLeg02::Gen3DObjectList();
	link_RLeg[2]=RLeg03::Gen3DObjectList();
	link_RLeg[3]=RLeg04::Gen3DObjectList();
	link_RLeg[4]=RLeg05::Gen3DObjectList();
	link_RLeg[5]=RLeg06::Gen3DObjectList();
	link_waist[0]=waist01::Gen3DObjectList();
	link_waist[1]=waist02::Gen3DObjectList();
	link_waist[2]=0;
	
	link_LArm[0]=LArm01::Gen3DObjectList();
	link_LArm[1]=LArm02::Gen3DObjectList();
	link_LArm[2]=LArm03::Gen3DObjectList();
	link_LArm[3]=LArm04::Gen3DObjectList();
	link_LArm[4]=LArm05::Gen3DObjectList();
	link_LArm[5]=LArm06::Gen3DObjectList();
	link_LArm[6]=0;
	
	link_RArm[0]=RArm01::Gen3DObjectList();
	link_RArm[1]=RArm02::Gen3DObjectList();
	link_RArm[2]=RArm03::Gen3DObjectList();
	link_RArm[3]=RArm04::Gen3DObjectList();
	link_RArm[4]=RArm05::Gen3DObjectList();
	link_RArm[5]=RArm06::Gen3DObjectList();
	link_RArm[6]=0;
}
void GLModel::gDHCallLink(double theta_home,double d,double a,double alpha,double theta,GLuint link)
{
	glRotatef(theta_home*180/PI,0,0,1);//angle home
	glTranslatef( a , 0.0 , d);	
	glRotatef(alpha*180/PI,1,0,0);
	glRotatef(theta*180/PI,0,0,1);//angle 視需要做rad轉換
	glCallList(link);
}


void GLModel::gGLWholeBody(void)
{
	
glPushMatrix();//世界原點

	//glRotatef((gRAngZWorld+gLAngZWorld)/PI*90,0,1,0);

	glRotatef(-90,1,0,0);
	glRotatef(-90,0,0,1);
	glTranslatef(gKineAll.DHOrigin[0] , gKineAll.DHOrigin[1] , gKineAll.DHOrigin[2]);	
	glRotatef((gRAngZWorld+gLAngZWorld)/PI*90,0,0,1);
	//glutSolidCube(100);//test
	gDHCallLink(0,0,0,0,0,link_center);
	
	glPushMatrix();//機器人原點
	for(int i=0;i<=5;i++){
		gDHCallLink(*(gKineAll.FKLLeg->theta_home+i),*(gKineAll.FKLLeg->d+i),*(gKineAll.FKLLeg->a+i),*(gKineAll.FKLLeg->alpha+i)/*theta_home_LLeg[i],d_LLeg[i],a_LLeg[i],alpha_LLeg[i]*/,*(gKineAll.FKLLeg->theta+i+1),link_LLeg[i]);
	}
	glPopMatrix();
	glPushMatrix();//機器人原點
	for(int i=0;i<=5;i++){
		gDHCallLink(*(gKineAll.FKRLeg->theta_home+i),*(gKineAll.FKRLeg->d+i),*(gKineAll.FKRLeg->a+i),*(gKineAll.FKRLeg->alpha+i)/*theta_home_RLeg[i],d_RLeg[i],a_RLeg[i],alpha_RLeg[i]*/,*(gKineAll.FKRLeg->theta+i+1),link_RLeg[i]);
	}
	glPopMatrix();
	//機器人原點向上移動
	for(int i=0;i<=2;i++){
		//gDHCallLink(theta_home_waist[i],d_waist[i],a_waist[i],alpha_waist[i],0/**(gKineAll.FKLArm->theta+i+1)*/,link_waist[i]);//self
		gDHCallLink(*(gKineAll.FKLArm->theta_home+i),*(gKineAll.FKLArm->d+i),*(gKineAll.FKLArm->a+i),*(gKineAll.FKLArm->alpha+i),*(gKineAll.FKLArm->theta+i+1),link_waist[i]);
	}
	
	glPushMatrix();//機器人上方原點
	for(int i=0;i<=6;i++){
		//gDHCallLink(theta_home_LArm[i],d_LArm[i],a_LArm[i],alpha_LArm[i],0/**(gKineAll.FKLArm->theta+(i+4))*/,link_LArm[i]);//self
		gDHCallLink(*(gKineAll.FKLArm->theta_home+i+3),*(gKineAll.FKLArm->d+i+3),*(gKineAll.FKLArm->a+i+3),*(gKineAll.FKLArm->alpha+i+3),*(gKineAll.FKLArm->theta+i+4),link_LArm[i]);
	}
	glPopMatrix();
	
	glPushMatrix();//機器人上方原點
	for(int i=0;i<=6;i++){
		//gDHCallLink(theta_home_RArm[i],d_RArm[i],a_RArm[i],alpha_RArm[i],0/**(gKineAll.FKRArm->theta+(i+4))*/,link_RArm[i]);//self
		gDHCallLink(*(gKineAll.FKRArm->theta_home+i+3),*(gKineAll.FKRArm->d+i+3),*(gKineAll.FKRArm->a+i+3),*(gKineAll.FKRArm->alpha+i+3),*(gKineAll.FKRArm->theta+i+4),link_RArm[i]);
	}
	glPopMatrix();
glPopMatrix();
}