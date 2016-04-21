// OpenGLControl.cpp: Implementierungsdatei
//

#include "stdafx.h"
#include ".\openglcontrol.h"
#include "math.h"
//#include "cmodel.h"


#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

#ifndef pi
#define pi 3.1415926
#endif


/////////////////////////////////////////////////////////////////////////////
// COpenGLControl


COpenGLControl::COpenGLControl()
: r_mouse_button(false)
, m_mouse_button(false)
, l_mouse_button(false)
, r_mouse_ini_pos(0)
, m_mouse_ini_pos(0)
, L_Button_Pos(0)
, PerspectiveView(false)
{
	dc = NULL;
	rotation = 0.0f;
}

COpenGLControl::~COpenGLControl()
{
	if (dc)
	{
		delete dc;
	}
}


BEGIN_MESSAGE_MAP(COpenGLControl, CWnd)
	//{{AFX_MSG_MAP(COpenGLControl)
	ON_WM_PAINT()
	ON_WM_SIZE()
	ON_WM_CREATE()
	ON_WM_ERASEBKGND()
	//}}AFX_MSG_MAP
	ON_WM_TIMER()
	ON_WM_RBUTTONDOWN()
	ON_WM_MOUSEMOVE()
	ON_WM_RBUTTONUP()
	ON_WM_MBUTTONDOWN()
	ON_WM_MBUTTONUP()
//	ON_WM_LBUTTONDBLCLK()
	ON_WM_LBUTTONUP()
	ON_WM_LBUTTONDOWN()
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// Behandlungsroutinen f�r Nachrichten COpenGLControl 


void COpenGLControl::InitGL()
{
///////////////////////////////////new addition by slongz
	drawswitch=false;
	ballsize=0.012;
	linewidth=2.5;
	drawscale=400;
///////////////////////////////////new addition by slongz
	//glShadeModel(GL_SMOOTH);
	//glClearDepth(1.0f);							
	//glEnable(GL_DEPTH_TEST);
	//glDepthFunc(GL_LEQUAL);	
	//glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);

	////--------------Initial Camera----------------------------
	//eye[0]    = 1.0;	eye[1]    = 0;	eye[2]    = 0.0;
	//center[0] = 0;	center[1] = 0;	center[2] = 0;
	//up[0]     = 0;  up[1]     = 0;	up[2]     = 1;

	//
	////--------initial 場景轉換的參數--------------------------------------------------
	//rot_m [0] = 0 ;
	//rot_m [1] = 0 ;
	//trans_m[0] = 0;
	//trans_m[1] = 0;
	//scale_m=1.0;

	//// setup lighting
	////===================================================
	//GLfloat glfLightAmbient[]  = { 0.3f, 0.3f, 0.3f, 1.0f };
 //   GLfloat glfLightDiffuse[]  = { 0.5f, 0.5f, 0.5f, 1.0f };
 //   GLfloat glfLightSpecular[] = { 0.4f, 0.4f, 0.4f, 1.0f };
	//GLfloat glfLightPosition[] = {0,0,30,0.0f};

 //   // Add a light to the scene.   //
 //   glLightfv (GL_LIGHT0, GL_AMBIENT, glfLightAmbient);
 //   glLightfv (GL_LIGHT0, GL_DIFFUSE, glfLightDiffuse);
 //   glLightfv (GL_LIGHT0, GL_SPECULAR, glfLightSpecular);
	//glLightfv (GL_LIGHT0, GL_POSITION,glfLightPosition);
	//
 //  // glEnable (GL_LIGHTING);
 //  // glEnable (GL_LIGHT0);

	//// Spot light setting
	//GLfloat glfLightDiffuse1[]  = { 0.7f, 0.7f, 0.7f, 1.0f };
	//GLfloat glfLightSpecular1[] = { 0.5f, 0.5f, 0.5f, 1.0f };
	//GLfloat glfLightPosition1[] = {0,0.5,12,0.0f};
	//glLightfv(GL_LIGHT1, GL_DIFFUSE ,glfLightDiffuse1);
	//glLightfv(GL_LIGHT1, GL_POSITION,glfLightPosition1);
 //   glLightfv(GL_LIGHT1, GL_SPECULAR,glfLightSpecular1);
	//glLightf (GL_LIGHT1, GL_SPOT_CUTOFF   , 60.0);
	//glLightf (GL_LIGHT1, GL_SPOT_EXPONENT ,100.0);
	////glEnable (GL_LIGHTING);
	////glEnable (GL_LIGHT1);


 //   //
 //   // Enable depth testing and backface culling.
 //   //

 //   //glEnable (GL_DEPTH_TEST);    
 //   //glDisable (GL_CULL_FACE);
	////glLightModelfv(GL_LIGHT_MODEL_AMBIENT ,glfLightAmbient);
 //   glLightModeli(GL_LIGHT_MODEL_TWO_SIDE ,1);

	////glEnable(GL_COLOR_MATERIAL);
	//
	//glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	SetTimer(0,100,NULL);//開啟計時器
}

void COpenGLControl::DrawGLScene()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	
	glLoadIdentity();

	//-----------------------場景旋轉縮放---------------------------------------------	
	//-----------------------場景縮放位移----------------------
	if(PerspectiveView)
	{
		//glMatrixMode(GL_PROJECTION);						
		//glLoadIdentity();
		//gluPerspective(60.0f,1,0.1f,10.0f);		
		//glMatrixMode(GL_MODELVIEW);

		center[0] = -cos(rot_m[0])*sin(rot_m[1])*trans_m[1] + sin(rot_m[0])*trans_m[0];
		center[1] =	-sin(rot_m[0])*sin(rot_m[1])*trans_m[1] - cos(rot_m[0])*trans_m[0];
		center[2] =  cos(rot_m[1])*trans_m[1];
	}
	else
	{
		glMatrixMode(GL_PROJECTION);						
		glLoadIdentity();
		glOrtho((-scale_m-trans_m[0]),(scale_m-trans_m[0]),(-scale_m+trans_m[1]),(scale_m+trans_m[1]),-4*(scale_m+1),4*(scale_m+1));	
		glMatrixMode(GL_MODELVIEW);
	}

	//-----------------------攝影機位置------------------------
	eye[0] = scale_m*cos(rot_m[0])*cos(rot_m[1])+center[0];
	eye[1] = scale_m*sin(rot_m[0])*cos(rot_m[1])+center[1];
	eye[2] = scale_m*sin(rot_m[1])+center[2];
	up[0] = -cos(rot_m[0])*sin(rot_m[1]);
	up[1] = -sin(rot_m[0])*sin(rot_m[1]);
	up[2] =  cos(rot_m[1]);	
	gluLookAt(eye[0],eye[1],eye[2],center[0],center[1],center[2],up[0],up[1],up[2]);
	//---------------------------------------------------------------------------------
	glGetDoublev(GL_MODELVIEW_MATRIX,ModelView);
	glGetDoublev(GL_PROJECTION_MATRIX,Projection);
	glGetIntegerv(GL_VIEWPORT,ViewPort);

	glPushMatrix();


	


	


	////glLineWidth(1.0f);

	//glColor3f(1.0f,1.0f,0.0f);

	//glBegin(GL_LINES);
	//	glVertex3f(0.0f,0.0f,0.0f);
	//	glVertex3f(0.0f,1.0f,0.0f);
	//glEnd();

	////glLineWidth(1.0f);



	glPopMatrix();


	//glLineWidth(10.0f);

	//glBegin(GL_LINE);
	//	glVertex3f(0.0f,0.0f,0.0f);
	//	glVertex3f(100.0f,100.0f,0.0f);
	//glEnd();

	//glBegin(GL_LINE);
	//	glVertex3f(0.0f,0.0f,0.0f);
	//	glVertex3f(100.0f,100.0f,100.0f);
	//glEnd();

	//glLineWidth(1.0f);

	SwapBuffers(dc->m_hDC);
}

void COpenGLControl::DrawPanelText()
{
	char xyz[80];//time[80], position[120], imu[60], encoder[60], mouse_pos[60], overlap[80];
	char *c ;

	//sprintf_s(time,"Time = %6.3f, Coverage_Rate = %6.2f, Coverage_Rate2 = %6.2f ", g_current_t, g_Map.m_CoverageRate,g_Map.m_CoverageRate2);
	//sprintf_s(position,"x = %6.2f  y = %6.2f  theta = %6.2f  x2 = %6.2f y2 = %6.2f  theta2 = %6.2f",
	//	g_Robot_Pos[0], g_Robot_Pos[1], g_Robot_Pos[2]*180/PI,g_Robot_Pos2[0], g_Robot_Pos2[1], g_Robot_Pos2[2]*180/PI);
	//sprintf_s(imu," yaw = %6.2f  roll = %6.2f  pitch = %6.2f", g_IMU[0], g_IMU[1], g_IMU[2]);	
	//sprintf_s(encoder," encoder R = %d     L = %d", int(g_encoder[0]), int(g_encoder[1]));	
	//sprintf_s(mouse_pos,"Pos x = %6.2f, y = %6.2f", m_mouse_pos[0],m_mouse_pos[1]);	
	//sprintf_s(overlap,"OverlapRate = %6.2f, TotalCoverageRate = %6.2f ,toc=%6i", g_Map.m_TotalOverlapRate,g_Map.m_TotalCoverageRate,g_Map.nowCoverage);	

	sprintf_s(xyz,"X = %6.2f, Y = %6.2f ,Z=%6i", 1,22223,3,44445,6);	
			
	glPushMatrix();
		//glTranslatef(-100.0,-3.0,0);		
		glColor3f(0,0,0);

		glRasterPos2f(-0,-0);
		for ( c = xyz ; *c!= '\0';c++)
			glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24 ,*c);	

		//glRasterPos2f(-600,-360);
		//for ( c = time ; *c!= '\0';c++)
		//	glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_10 ,*c);
		//
		//glRasterPos2f(-600,-400);
		//for ( c = position ; *c!= '\0';c++)
		//	glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_10 ,*c);

		//glRasterPos3f(-600,-440,10);
		//for ( c = imu ; *c!= '\0';c++)
		//	glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_10 ,*c);

		//glRasterPos3f(-600,-480,10);
		//for ( c = encoder ; *c!= '\0';c++)
		//	glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_10 ,*c);	

		

	glPopMatrix();

	// Mouse pos information
	//glPushMatrix();
	//glRasterPos2f((float(m_mouse_pos_pixel[0])-m_WindowSize[0]/2.0)*1.6,((m_WindowSize[1]-m_mouse_pos_pixel[1])-m_WindowSize[1]/2.0)*1.7);
	//
	//for ( c = mouse_pos ; *c!= '\0';c++)
	//		glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_10 ,*c);
	glPopMatrix();	
}

void COpenGLControl::Create(CRect rect, CWnd *parent)
{
	CString className = AfxRegisterWndClass(
		CS_HREDRAW | CS_VREDRAW | CS_OWNDC,
		NULL,
		(HBRUSH)GetStockObject(BLACK_BRUSH),
		NULL);

	CString a;
	a ="OpenGL";
	CreateEx(
		0,
		className,
		a,
		WS_CHILD | WS_VISIBLE | WS_CLIPSIBLINGS | WS_CLIPCHILDREN,
		rect,
		parent,
		0);

}

void COpenGLControl::OnPaint() 
{
	/** OpenGL section **/

	openGLDevice.makeCurrent();

	DrawGLScene();
	
	ValidateRect(NULL);
}

void COpenGLControl::OnSize(UINT nType, int cx, int cy) 
{
	CWnd::OnSize(nType, cx, cy);	
	if (cy == 0)								
	{
		cy = 1;						
	}

	glViewport(0,0,cx,cy);	

	glMatrixMode(GL_PROJECTION);						
	glLoadIdentity();


//	glOrtho(-50.0f,50.0f,-50.0f,50.0f,50.0f,-50.0f);
	gluPerspective(60.0f,cx/cy,1.0f,1000.0f);
  //glFrustum(-50,50,-50,50,10,0);

	glMatrixMode(GL_MODELVIEW);						
	glLoadIdentity();
}


int COpenGLControl::OnCreate(LPCREATESTRUCT lpCreateStruct) 
{
	if (CWnd::OnCreate(lpCreateStruct) == -1)
		return -1;
	
	dc = new CClientDC(this);

	openGLDevice.create(dc->m_hDC);
	InitGL();
	//SetTimer(0,5,NULL);

	return 0;
}

BOOL COpenGLControl::OnEraseBkgnd(CDC* pDC) 
{
	return TRUE;
}

void COpenGLControl::OnTimer(UINT nIDEvent)
{
	// TODO: Add your message handler code here and/or call default

/*
	rotation += 1.0f;
	if (rotation >= 360.0f)
	{
		rotation -= 360.0f;
	}
*/

	InvalidateRect(NULL,FALSE);

	CWnd::OnTimer(nIDEvent);
}


void COpenGLControl::OnRButtonDown(UINT nFlags, CPoint point)
{
	// TODO: Add your message handler code here and/or call default
	r_mouse_button = true;
	r_mouse_ini_pos = point;

	CWnd::OnRButtonDown(nFlags, point);
}

void COpenGLControl::OnMouseMove(UINT nFlags, CPoint point)
{
	// TODO: Add your message handler code here and/or call default
	float scale_sensitivity = 0.01;
	float trans_sensitivity = 0.02;
	float rot_sensitivity   = 0.005;

	if(r_mouse_button == true && m_mouse_button == false)//只按右鍵 旋轉
	{
		float rotx=(point.x-r_mouse_ini_pos.x)*rot_sensitivity;
		float roty=(point.y-r_mouse_ini_pos.y)*rot_sensitivity;
		rot_m[0]-=rotx;
		rot_m[1]+=roty;		
		r_mouse_ini_pos=point;
	}

	if(r_mouse_button && m_mouse_button)//按右鍵和中鍵 縮放
	{
		float scale_temp =- (point.y-m_mouse_ini_pos.y)*scale_sensitivity;
		scale_m += scale_temp;
		m_mouse_ini_pos = point;
	}

	if(m_mouse_button == true && r_mouse_button == false)//指案中鍵 位移
	{
		float transx=(point.x-m_mouse_ini_pos.x)*trans_sensitivity;
		float transy=(point.y-m_mouse_ini_pos.y)*trans_sensitivity;
		trans_m[0] += transx;
		trans_m[1] += transy;
		m_mouse_ini_pos = point;
	}

	//if(l_mouse_button == true)//只按左鍵(移動Goal的部分)
	//{
	//	//===============轉換視角(視窗轉OpenGL)==============
	//	float winPoint[] = {(float)point.x , (float)(ViewPort[3]-point.y) , 0};
	//	double objPoint[3] = {0};
	//	//glReadPixels( int(winPoint[0]), int(winPoint[1]), 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &winPoint[2] );
	//	if(gluUnProject(winPoint[0],winPoint[1],winPoint[2],ModelView,Projection,ViewPort,&objPoint[0],&objPoint[1],&objPoint[2]) != GL_TRUE)
	//		printf("UnProject is ERROR\n");
	//	//====================================================
	//	if(MovingOBS)
	//	{
	//		StaticOBS[5]->Trans[0]+=(objPoint[0]-l_mouse_Pos[0])*0.5;
	//		StaticOBS[5]->Trans[1]+=(objPoint[1]-l_mouse_Pos[1])*0.5;
	//		StaticOBS[5]->Trans[2]+=(objPoint[2]-l_mouse_Pos[2])*0.5;
	//	}
	//	else
	//	{
	//		goal[0]+=(objPoint[0]-l_mouse_Pos[0])*2.5;
	//		goal[1]+=(objPoint[1]-l_mouse_Pos[1])*2.5;
	//		goal[2]+=(objPoint[2]-l_mouse_Pos[2])*2.5;
	//	}

	//	for(int i = 0 ; i < 3 ; i++)
	//		l_mouse_Pos[i] = objPoint[i];
	//}
	CWnd::OnMouseMove(nFlags, point);
}


void COpenGLControl::OnRButtonUp(UINT nFlags, CPoint point)
{
	// TODO: Add your message handler code here and/or call default
	r_mouse_button = false;

	CWnd::OnRButtonUp(nFlags, point);
}
void COpenGLControl::OnMButtonDown(UINT nFlags, CPoint point)
{
	// TODO: 在此加入您的訊息處理常式程式碼和 (或) 呼叫預設值
	m_mouse_button  = true;
	m_mouse_ini_pos = point;

	CWnd::OnMButtonDown(nFlags, point);
}

void COpenGLControl::OnMButtonUp(UINT nFlags, CPoint point)
{
	// TODO: 在此加入您的訊息處理常式程式碼和 (或) 呼叫預設值

	m_mouse_button=false;

	CWnd::OnMButtonUp(nFlags, point);
}

void COpenGLControl::OnLButtonDown(UINT nFlags, CPoint point)
{
	// TODO: Add your message handler code here and/or call default
	l_mouse_button  = true;
	l_mouse_ini_pos = point;
	
	//===============轉換視角(視窗轉OpenGL)==============
	float winPoint[] = {(float)point.x , (float)(ViewPort[3]-point.y) , 0};
	double objPoint[3] = {0};
	//glReadPixels( int(winPoint[0]), int(winPoint[1]), 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &winPoint[2] );
	if(gluUnProject(winPoint[0],winPoint[1],winPoint[2],ModelView,Projection,ViewPort,&objPoint[0],&objPoint[1],&objPoint[2]) != GL_TRUE)
		printf("UnProject is ERROR\n");
	for(int i = 0 ; i < 3 ; i++)
		l_mouse_Pos[i] = objPoint[i];

	L_Button_Pos = point;
	L_Button_Pos.x=( (long)(L_Button_Pos.x)-300)/5.3/scale_m;
	L_Button_Pos.y=(-(long)(L_Button_Pos.y)+300)/5.3/scale_m;

	CWnd::OnLButtonDown(nFlags, point);
}

void COpenGLControl::OnLButtonUp(UINT nFlags, CPoint point)
{
	// TODO: Add your message handler code here and/or call default
	l_mouse_button = false;

	CWnd::OnLButtonUp(nFlags, point);
}

