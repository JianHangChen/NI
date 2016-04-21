#pragma once
#ifndef _GLMODEL_H_
#define _GLMODEL_H_

#include <GL\gl.h>



class GLModel
{
public:
	GLModel(void);
	~GLModel(void);
	
	void InitGLModel(void); // ªì©l¤Æ
	void gDHCallLink(double theta_home,double d,double a,double alpha,double theta,GLuint link);
	void gGLWholeBody(void);

	//double gRadRatio /*= 180/3.1415926*/;
	
	
	//double theta_home_LLeg[6];//1-6
	//double theta_home_RLeg[6];
	//double theta_home_LArm[7];//4-9 10=endeffector
	//double theta_home_RArm[7];
	//double theta_home_waist[3];//1-3

	//double d_LLeg[6]={ -99.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 };//1-6
	//double d_RLeg[6]={ -99.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 };
	//double d_LArm[7]/*={ 0.0 , 245.5 , 0.0 , 267.0 , 0.0 , 208.0 , 0.0 }*/;//4-9 10=endeffector
	//double d_RArm[7]/*={ 0.0 , -245.5 , 0.0, -267.0 , 0.0 , -208.0 , 0.0 }*/;
	//double d_waist[3]/*={ 99.0 , 0.0 , 0.0 }*/;//1-3

	//double a_LLeg[6]={ -80.0 , 0.0 , 0.0 , -240.0 , -240.0 , 0.0 };//1-6
	//double a_RLeg[6]={ 80.0 , 0.0 , 0.0 , -240.0 , -240.0 , 0.0 };
	//double a_LArm[7]/*={ 0.0 , 0.0 , 0.0 , -10.0 , 10.0 , 0.0 , 156.0 }*/;//4-9 10=endeffector
	//double a_RArm[7]/*={ 0.0 , 0.0 , 0.0 , -10.0 , 10.0 , 0.0 , 156.0 }*/;
	//double a_waist[3]/*={ 0.0 , 0.0 , 271.75 }*/;//1-3

	//double alpha_LLeg[6]={ 0.0 , -90.0 , 90.0 , 0.0 , 0.0 , -90.0 };//1-6
	//double alpha_RLeg[6]={ 0.0 , -90.0 , -90.0 , 0.0 , 0.0 , 90.0 };
	//double alpha_LArm[7]/*={ 0.0 , -90.0 , -90.0 , 90.0 , -90.0 , 90.0 , 0.0 }*/;//4-9 10=endeffector
	//double alpha_RArm[7]/*={ 0.0 , 90.0 , 90.0 , -90.0 , 90.0 , -90.0 , 0.0 }*/;
	//double alpha_waist[3]/*={ 0.0 , 90.0 , 0.0 }*/;//1-3

	//double theta_LLeg[6]={ 0.0 , 0.222145 , -0.609947 , 0.668909 , -0.0589622 , -0.222145 };//1-6
	//double theta_RLeg[6]={0.0};
	//double theta_LArm[7]/*={0.0}*/;//4-9 10=endeffector
	//double theta_RArm[7]/*={0.0}*/;
	//double theta_waist[3]/*={0.0}*/;//1-3

	GLuint link_center;
	GLuint link_LLeg[6];//1-6
	GLuint link_RLeg[6];
	GLuint link_LArm[7];//4-9 10=endeffector 4 empty
	GLuint link_RArm[7];
	GLuint link_waist[3];//1-2 3 empty
};
#endif