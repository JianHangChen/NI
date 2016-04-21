#if !defined(AFX_OPENGLCONTROL_H__52A6B63B_01A2_449D_8691_1FF59EECAB71__INCLUDED_)
#define AFX_OPENGLCONTROL_H__52A6B63B_01A2_449D_8691_1FF59EECAB71__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
// OpenGLControl.h : Header-Datei
//

#include "OpenGLDevice.h"
//#include <GL\gl.h>
//#include <GL\glu.h>
#include "glut.h"
//#include <GL\glaux.h>
//#include "ROBOTDlg.h"
//#include "ROBOT.h"
#include "atltypes.h"
#include "vector"



/////////////////////////////////////////////////////////////////////////////
// Fenster COpenGLControl 

class COpenGLControl : public CWnd
{
// Konstruktion
public:
	COpenGLControl();

// Attribute
public:

// Operationen
public:

// 鈁erschreibungen
	// Vom Klassen-Assistenten generierte virtuelle Funktions�berschreibungen
	//{{AFX_VIRTUAL(COpenGLControl)
	public:
	//}}AFX_VIRTUAL

// Implementierung
public:
	void Create(CRect rect,CWnd* parent);
	virtual ~COpenGLControl();




	// Generierte Nachrichtenzuordnungsfunktionen
protected:
	void InitGL();
	void DrawGLScene();
	OpenGLDevice openGLDevice;
	CClientDC* dc;
	float rotation;

	//{{AFX_MSG(COpenGLControl)
	afx_msg void OnPaint();
	afx_msg void OnSize(UINT nType, int cx, int cy);
	afx_msg int OnCreate(LPCREATESTRUCT lpCreateStruct);
	afx_msg BOOL OnEraseBkgnd(CDC* pDC);
	//}}AFX_MSG

	DECLARE_MESSAGE_MAP()
public:	
	afx_msg void OnTimer(UINT nIDEvent);
	bool Set3DViewPort(void);
public:
	afx_msg void OnRButtonUp(UINT nFlags, CPoint point);
	afx_msg void OnRButtonDown(UINT nFlags, CPoint point);	
	afx_msg void OnMButtonDown(UINT nFlags, CPoint point);
	afx_msg void OnMButtonUp(UINT nFlags, CPoint point);	
	afx_msg void OnLButtonUp(UINT nFlags, CPoint point);	
	afx_msg void OnMouseMove(UINT nFlags, CPoint point);
	afx_msg void OnLButtonDown(UINT nFlags, CPoint point);

public:
	//--------------------使用按鍵的參數--------------------------------
	bool r_mouse_button;// 看滑鼠右鍵是否被按住
	bool m_mouse_button;
	bool l_mouse_button;
	CPoint r_mouse_ini_pos; //按右鍵滑鼠座標起始位置
	CPoint m_mouse_ini_pos; //按中建滑鼠座標起始位置
	CPoint l_mouse_ini_pos;
	CPoint L_Button_Pos;
	double l_mouse_Pos[3];

	//--------------------宣告場景轉換的參數-------------------------
	float rot_m[2];
	float trans_m[2];
	float scale_m;
	bool PerspectiveView;

	//--------Camera--------------
	GLdouble eye[3];
	GLdouble center[3];
	GLdouble up[3];

	double ModelView [16];
	double Projection[16];
	int    ViewPort[4];


	struct GL_3Dmodel
	{
		GLint   ID;
		GLfloat Pos[3];
	}m3D_model;


	/////////////////////////New addition by Slongz
	void DrawPanelText();


	bool drawswitch;
	float ballsize;
	float linewidth;
	int drawscale;
	float xLLeg[13],yLLeg[13],zLLeg[13],xRLeg[13],yRLeg[13],zRLeg[13];
	float xLArm[13],yLArm[13],zLArm[13],xRArm[13],yRArm[13],zRArm[13];
	float xpv[15],ypv[15],zpv[15];
	float COGx,COGy,COGz;
};

/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ f�gt unmittelbar vor der vorhergehenden Zeile zus酹zliche Deklarationen ein.

#endif // AFX_OPENGLCONTROL_H__52A6B63B_01A2_449D_8691_1FF59EECAB71__INCLUDED_
