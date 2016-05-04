/**************************************************************************************************
Copyright, 2010-2012, Robotics Lab., Dept. of M.E., National Taiwan University
File Name: MatOptInterface.cpp

Author: Jiu-Lou Yan
Version: 1.0
Date: 2010/07/10

Functions:
     InvSqMat() AbsDetSqMat()  SubMat() SubMatIndex() MatScalarMul()
     MatMulAB() MatMulAtB() MatMulABt() MatMiuAB() MatAddAB()

	Mat_weighted_mean() Mat_I_nxn() Mat_O_mxn() Mat_y_Ax_b() Mat_A_addcol_B( ) 
	Mat_A_addrow_B( ) Mat_A_B_At( ) showMat() showVar()
Classes: None

Description:
     本程式主要用在建立一套矩陣運算函式庫、定義了矩陣運算之函式，
	 包含反矩陣、行列式值、矩陣加減乘法之運算，
	 各函式與變數之說明請詳見下方宣告與定義處之說明

Note: None
***************************************************************************************************/
#include "stdafx.h"
#include "MatOptInterface.h"



void InvSqMat(double * Data, int NCol)
{
	/******************************************************************
	input: 矩陣本身開頭指標，矩陣column數
	output: 直接儲存在矩陣本身位置

	Note: 輸入一定要是方陣!!
	******************************************************************/

	// 準備丟給CLAPACK的變數 還有記憶體分配
	integer Dim = NCol;
	doublereal *work = new doublereal[Dim*Dim];
    integer INFO; 

    integer * ipiv = new integer[Dim*Dim];

	// CLAPACK 的 LU decomposition
    dgetrf_(&Dim,&Dim,Data,&Dim,ipiv,&INFO);

	if (INFO != 0)
	{
		cout << "Fail to Calculate LU decomposition" << endl;
		system("pause");
	}
	// CLAPACK 的 矩陣 inverse
    dgetri_(&Dim,Data,&Dim,ipiv,work,&Dim,&INFO);

	if (INFO != 0)
	{
		cout << "Fail to Calculate Matrix Inversion" << endl;
		system("pause");
	}

	// 清除動態記憶體
	delete[] work;
	delete[] ipiv;

}


double AbsDetSqMat(double * Data, int NCol)
{
	/******************************************************************
	input: 矩陣本身開頭指標，矩陣column數
	output: 回傳行列式值的絕對值

	Note:
	// 注意!! 此操作只回傳正值，假若Deterimant是負的，會回傳正的
	// 取deterimant 不該破壞矩陣值，要不然之後的Jacobian會用不了
	******************************************************************/

	double detValue = 1;

	// copy the elements of data
	double* LU = new double[NCol*NCol];
	for (int i =0 ; i < NCol*NCol ; i++)
		LU[i] = Data[i];

	// 準備要做CLAPACK的 LU decomposition
	// 分配記憶體與創建變數以供使用
	integer Dim = NCol;
	doublereal *work = new doublereal[Dim*Dim];
    integer lwork = Dim;
    integer INFO; 
    integer * ipiv = new integer[Dim*Dim];

	// LU decomposition
    dgetrf_(&Dim,&Dim,LU,&Dim,ipiv,&INFO);

	// detValue = U矩陣的對角線連乘
	for (int j = 0; j < Dim; j++)
		detValue *= LU[j*Dim+j];

	if (detValue < 0)
		detValue = -detValue;

	// 釋放動態記憶體
	delete[] work;
	delete[] ipiv;
	delete[] LU;

	return detValue;

}


void SubMat(double *MMother, int ColMother, int RowMother, double *MSon, int *CutRegion)
{
	/******************************************************************
	input: *MMother矩陣指標，ColMother母矩陣column數，RowMother母矩陣row數，*MSon子矩陣指標，*CutRegion剪裁區間
	output: void

	Note:
	// 取得要剪下矩陣的微度， CutRegion[0]代表Row的下限  
	   CutRegion[1]代表Row的上限 CutRegion[2]代表Col的下限  CutRegion[3]代表Col的上限
	// 不能跳著剪，要跳著剪，請使用SubMatIndex
	******************************************************************/
	int RowSon = CutRegion[1]-CutRegion[0]+1; 
	int ColSon = CutRegion[3]-CutRegion[2]+1;

	// 矩陣剪取
	for (int i = 0 ; i < RowSon ; i++)
		for (int j = 0 ; j < ColSon ; j++)
			MSon[i*ColSon+j] = MMother[(i+CutRegion[0]-1)*ColMother+CutRegion[2]+j-1];


}


void SubMatIndex(double *MMother, int ColMother, int RowMother, double *MSon, int *ColRegion, int *RowRegion)
{
	/******************************************************************
	input: *MMother矩陣指標，ColMother母矩陣column數，RowMother母矩陣row數，*MSon子矩陣指標，*ColRegion *RowRegion 為剪裁區間
	output: void

	Note:
	      基本上等同於SubMat()，多了可以跳著剪的功能
	******************************************************************/	


	int RowSon = RowRegion[0]; // 丟進去的第一個值，代表剪出來的矩陣長度(有幾Row) 後面的部份代表要選取的Row數，可以跳著填
	int ColSon = ColRegion[0]; // 丟進去的第一個值，代表剪出來的矩陣長度(有幾Col) 後面的部份代表要選取的Col數，可以跳著填

	// 矩陣剪取
	for (int i = 1 ; i <= RowSon ; i++)
		for (int j = 1 ; j <= ColSon ; j++)
			MSon[(i-1)*ColSon+j-1] = MMother[(RowRegion[i]-1)*ColMother+ColRegion[j]-1];

}

void MatMulAB(double* A, int MA, int NA, double* B, int MB, int NB, double* result)
{

	/******************************************************************
	input:  A 矩陣A開頭指標 , MA 矩陣A row數 , NA 矩陣A col數, B 矩陣B開頭指標,  MB 矩陣B row數,  NB 矩陣B col數, result 儲存結果開頭指標
	output: void

	Note:
	// 矩陣乘法 兩個矩陣的維度要對到 不可以是不可乘的矩陣
	******************************************************************/

	if (NA != MB)
	{
		printf("Wrong Matrix Dimension for Mat_Mul");
		system("pause");
	}

	double tempAcc = 0;
	int p = 0;
	int r = 0;

	// 執行矩陣乘法
	for (int i = 0; i < MA*NB ; i++)
	{
		for (int j = 0 ; j < NA; j ++)
		{
			tempAcc += A[p*NA+j]*B[r+j*NB];
		}
		result[i] = tempAcc;
		tempAcc = 0;

		if (r == NB-1)
		{
			r = 0;
			p += 1;
		}
		else
		{
			r += 1;
		}

	}
}

//
//void MatMulAB(double* A, int MA, int NA, double* B, int MB, intNB, double* result)
//{
//
//	if (NA != MB)
//	{
//		printf("Wrong Matrix Dimension for Mat_Mul");
//		system("pause");
//	}
//
//	// C++版本
//	//double tempAcc = 0;
//	//int p, r;
//
//	//for (int i = 0; i < MA*NB ; i++)
//	//{
//	//	p = i/NB;
//	//	r  =i%NB;
//	//	for (int j = 0 ; j < NA; j ++)
//	//	{
//	//		//for (int k = 0 ; k<MB ; k++)
//	//			tempAcc += A[p*NA+j]*B[r+j*NB];
//	//	}
//	//	result[i] = tempAcc;
//	//	tempAcc = 0;
//	//}
//
//	float dummy;
//	int Cnt = MA*NB-1;
//
//	__asm{
//
//
//		
//			mov ch, 8 // ch = 8 because double = 8 bytes
//			mov eax, 8 
//
//			mul NB    // jump a row = 8 * NB
//			mov cl, al // move jump a row to cl
//
//			mov esi, A	 // load matrix
//			mov edi, B	 // load matrix
//
//
//		mul_loop: // 計算矩陣乘法
//
//				mov esi, A	 // load matrix
//				mov edi, B	 // load matrix
//
//			
//				//mov eax, 0
//				mov edx,0 // initialize before div
//				mov eax, Cnt // move the index Cnt to eax, 指向輸入輸出矩陣的位置
//				div NB   // 求得第幾行，第幾列，會被存在 bh, bl
//
//				mov bh, al // al, eax = p
//				mov bl, dl // dl, edx = r
//
//				mul NA
//				add eax, edx
//
//				mov edx, 0
//				mov dl, ch
//				mul edx
//
//				add esi, eax
//
//
//				mov eax, 0 // clear 高位元
//				mov al, dl
//				mul NB
//				add al, bl
//				mul ch
//				add edi, eax
//
//				mov dl, 0  // clear dl, 給 for 迴圈用
//				mov dh, byte ptr NA
//
//				fldz		 // st = 0
//
//			add_loop:
//
//				fld qword ptr [esi]
//				fmul qword ptr [edi]
//				fadd st(1),st(0)
//				fstp dummy
//
//				inc dl  // j ++
//
//				add esi,8
//				mov eax, 0
//				mov al, cl
//				add edi, eax
//
//				cmp dl, dh
//				jnz	add_loop
//
//				mov eax, Cnt
//
//				mov ebx, 8
//				mul ebx
//				add eax, result
//			
//				fstp qword ptr [eax]
//
//				dec Cnt
//				cmp Cnt, -1
//				jz end_this_loop
//
//				jmp mul_loop
//
//			end_this_loop:
//
//	}
//
//}


//void MatMulAB(double* A, int MA, int NA, double* B, int MB, int NB, double* result)
//{
//
//	if (NA != MB)
//	{
//		printf("Wrong Matrix Dimension for Mat_Mul");
//		system("pause");
//	}
//
//	// C++版本
//	//double tempAcc = 0;
//	//int p, r;
//
//	//for (int i = 0; i < MA*NB ; i++)
//	//{
//	//	p = i/NB;
//	//	r  =i%NB;
//	//	for (int j = 0 ; j < NA; j ++)
//	//	{
//	//		//for (int k = 0 ; k<MB ; k++)
//	//			tempAcc += A[p*NA+j]*B[r+j*NB];
//	//	}
//	//	result[i] = tempAcc;
//	//	tempAcc = 0;
//	//}
//
//	float dummy;
//	int Cnt = MA*NB-1;
//
//	__asm{
//
//
//		
//			mov ch, 8 // ch = 8 because double = 8 bytes
//			mov eax, 8 
//
//			mul NB    // jump a row = 8 * NB
//			mov cl, al // move jump a row to cl
//
//			mov esi, A	 // load matrix
//			mov edi, B	 // load matrix
//
//
//		mul_loop: // 計算矩陣乘法
//
//				mov esi, A	 // load matrix
//				mov edi, B	 // load matrix
//
//			
//				//mov eax, 0
//				mov edx,0 // initialize before div
//				mov eax, Cnt // move the index Cnt to eax, 指向輸入輸出矩陣的位置
//				div NB   // 求得第幾行，第幾列，會被存在 bh, bl
//
//				mov bh, al // al, eax = p
//				mov bl, dl // dl, edx = r
//
//				mul NA
//				add eax, edx
//
//				mov edx, 0
//				mov dl, ch
//				mul edx
//
//				add esi, eax
//
//
//				mov eax, 0 // clear 高位元
//				mov al, dl
//				mul NB
//				add al, bl
//				mul ch
//				add edi, eax
//
//				mov dl, 0  // clear dl, 給 for 迴圈用
//				mov dh, byte ptr NA
//
//				fldz		 // st = 0
//
//			add_loop:
//
//				fld qword ptr [esi]
//				fmul qword ptr [edi]
//				fadd st(1),st(0)
//				fstp dummy
//
//				inc dl  // j ++
//
//				add esi,8
//				mov eax, 0
//				mov al, cl
//				add edi, eax
//
//				cmp dl, dh
//				jnz	add_loop
//
//				mov eax, Cnt
//
//				mov ebx, 8
//				mul ebx
//				add eax, result
//			
//				fstp qword ptr [eax]
//
//				dec Cnt
//				cmp Cnt, -1
//				jz end_this_loop
//
//				jmp mul_loop
//
//			end_this_loop:
//
//	}
//
//}


//
//// 第一版，穩定板，可以處理較大矩陣
//void MatMulAB(double* A, int MA, int NA, double* B, int MB, int NB, double* result)
//{
//
//	if (NA != MB)
//	{
//		printf("Wrong Matrix Dimension for Mat_Mul");
//		system("pause");
//	}
//
//	double dummy;
//
//	//int S_D = 8; // size of double
//	int p, r;
//	int warp = NB*8;
//
//	__asm{
//
//			
//
//			mov eax, MA
//			//mov esi, NB  // esi = row 數量
//			mul NB      // eax = eax * NB
//			dec eax
//			mov ecx, eax 
//
//
//		mul_loop: // 計算矩陣乘法
//
//			//fldz		 // st = 0
//			//fld1		 // st = 1 0
//
//			mov esi, A	 // load matrix
//			mov edi, B	 // load matrix
//
//			mov eax, ecx
//
//			mov edx, 0
//			mov ebx, 0
//			div NB		 // eax = 商 edx = 餘數	
//			mov r, edx
//			mov p, eax
//
//				// get p*NA+j
//				//mov eax, p
//				mul NA		 // p*NA
//				add eax, ebx		 // p*NA+j
//				rol eax, 3	// to byte address
//				add esi, eax
//
//				mov eax,ebx // j          
//				mul NB		// j*NB
//				add eax,r // r+ j*NB
//				rol eax, 3		// to byte address
//				add edi, eax 
//				
//				fldz		 // st = 0
//
//			add_loop:
//
//
//				fld qword ptr [esi]
//				fmul qword ptr [edi]
//				fadd st(1),st(0)
//				fstp dummy
//	
//
//				add esi,8
//				add edi,warp
//
//				inc ebx
//				mov eax, NA
//				cmp ebx, eax
//
//				jnz	add_loop
//
//				mov eax, ecx
//				rol eax, 3
//		
//				add eax, result
//				fstp qword ptr [eax]
//
//				dec ecx
//				cmp ecx, -1
//				jnz mul_loop
//
//	}
//
//}

void MatMulAtB(double* A, int MA, int NA, double* B, int MB, int NB, double* result)
{
	/******************************************************************
	input:  A 矩陣A開頭指標 , MA 矩陣A row數 , NA 矩陣A col數, B 矩陣B開頭指標,  MB 矩陣B row數,  NB 矩陣B col數, result 儲存結果開頭指標
	output: void

	Note:
	// 矩陣乘法 兩個矩陣的維度要對到 不可以是不可乘的矩陣
	// 也就是說 transpose(A)*B 要可乘 
	// 並且輸入的MA NA 是A原本的維度!!
	******************************************************************/
	if (MA != MB)
	{
		printf("Wrong Matrix Dimension for Mat_Mul");
		system("pause");
	}

	int p = 0;
	int r = 0;

	double tempAcc = 0;

	// 執行矩陣乘法
	for (int i = 0; i < NA*NB ; i++)
	{
		for (int j = 0 ; j < MA; j ++)
		{
				tempAcc += A[j*NA+p]*B[r+j*NB];
		}
		result[i] = tempAcc;
		tempAcc = 0;

		if (r == NB-1)
		{
			r = 0;
			p += 1;
		}
		else
		{
			r += 1;
		}
	}
}


void MatMulABt(double* A, int MA, int NA, double* B, int MB, int NB, double* result)
{
	/******************************************************************
	input:  A 矩陣A開頭指標 , MA 矩陣A row數 , NA 矩陣A col數, B 矩陣B開頭指標,  MB 矩陣B row數,  NB 矩陣B col數, result 儲存結果開頭指標
	output: void

	Note:
	// 矩陣乘法 兩個矩陣的維度要對到 不可以是不可乘的矩陣
	// 也就是說 A*transpose(B) 要可乘 
	// 並且輸入的MB NB 是B原本的維度!!
	******************************************************************/

	if (NA != NB)
	{
		printf("Wrong Matrix Dimension for Mat_Mul");
		system("pause");
	}

	int p = 0;
	int r = 0;
	double tempAcc = 0;

	// 執行矩陣乘法
	for (int i = 0; i < MA*MB ; i++)
	{
		for (int j = 0 ; j < NA; j ++)
		{
				tempAcc += A[p*NA+j]*B[j+r*NB];
		}
		result[i] = tempAcc;
		tempAcc = 0;

		if (r == MB-1)
		{
			r = 0;
			p += 1;
		}
		else
		{
			r += 1;
		}


	}
}


void MatScalarMul(double* A, int LenA, double* Scalar, double* result)
{
	/******************************************************************
	input:  A 矩陣A開頭指標 , LenA 矩陣A 元素數 , *Scalar 放大倍率, result 儲存結果開頭指標
	output: void

	Note: None
	******************************************************************/

	for (int i=0; i<LenA; i++)
	{
		result[i] = A[i]*Scalar[0];
	}
}

void MatScalarMul(double* A, int LenA, double Scalar, double* result)
{
	/******************************************************************
	input:  A 矩陣A開頭指標 , LenA 矩陣A 元素數 , *Scalar 放大倍率, result 儲存結果開頭指標
	output: void

	Note: None
	******************************************************************/

	for (int i=0; i<LenA; i++)
	{
		result[i] = A[i]*Scalar;
	}
}

void MatMiuAB(double* A, double* B, double* result, int Len) // 矩陣相減
{
	/******************************************************************
	input:  A 矩陣A開頭指標, B 矩陣B開頭指標,  result 儲存結果開頭指標, Len 為兩矩陣之元素數
	output: void

	Note:
	// 矩陣加減法 兩個矩陣的維度要一樣 

	// 原始碼
	for (int i = 0; i < Len ;  i++)
		result[i] = A[i]-B[i];

	******************************************************************/

	// 組合語言
	__asm{

		mov eax, A
		mov ebx, B
		mov edx, result

		mov ecx, Len // 是多少 就會進去幾次回圈

	// for loop
	mat_miu: // ecx 會自動減一直到0 跳出回圈

		fld QWORD PTR [eax] // load A[i] to FPU stack 0 ( st(0) )
		fsub QWORD PTR [ebx] // add B[i] to st(0)
		
		fstp QWORD PTR [edx] // pop the value 

		add eax,8
		add ebx,8
		add edx,8

		dec ecx
		cmp ecx, 0
		jnz mat_miu
		//loop mat_add
	 //for loop

	}



}

void MatAddAB(double* A, double* B, double* result, int Len) // 矩陣相加
{
	/******************************************************************
	input:  A 矩陣A開頭指標, B 矩陣B開頭指標,  result 儲存結果開頭指標, Len 為兩矩陣之元素數
	output: void

	Note:
	// 矩陣加減法 兩個矩陣的維度要一樣 

	// 原始碼
	for (int i = 0; i < Len ;  i++)
		result[i] = A[i]+B[i];

	******************************************************************/

	__asm{

		mov eax, A
		mov ebx, B
		mov edx, result

		mov ecx, Len // 是多少 就會進去幾次回圈

	// for loop
	mat_add: // ecx 會自動減一直到0 跳出回圈

		fld QWORD PTR [eax] // load A[i] to FPU stack 0 ( st(0) )
		fadd QWORD PTR [ebx] // add B[i] to st(0)
		
		fstp QWORD PTR [edx] // pop the value 

		add eax,8
		add ebx,8
		add edx,8

		//loop mat_add

		dec ecx
		cmp ecx, 0
		jnz mat_add

	// for loop

	}

}

void dare_sqrmat_matlab(double * Data, double * Data_U, int NCol)
{
	/******************************************************************
	input: 矩陣本身開頭指標，矩陣column數，
	output: 置放分解後U T矩陣的指標

	SUBROUTINE DGEES( JOBVS, SORT, SELECT, N, A, LDA, SDIM, WR,
	WI, VS, LDVS, WORK, LWORK, BWORK, INFO )
	A = U*T*U'

	Note: 輸入一定要是方陣!!
	******************************************************************/

	// 準備丟給CLAPACK的變數 還有記憶體分配
	integer Dim = NCol;
	doublereal *work = new doublereal[Dim*Dim];
    integer INFO; 

    //integer * ipiv = new integer[Dim*Dim];
	//new add
	integer SDIM;//If SORT = 'N', SDIM = 0.  If SORT = 'S', SDIM = number of eigenvalues (after sorting) for which SELECT is true. (Complex conjugate pairs for which SELECT is true for either eigenvalue count as 2.)
	double * WR = new double[Dim];
	double * WI = new double[Dim];
	L_fp SELECT=0;
	
	char JOBVS;
	char SORT;
	logical * BWORK = new logical[Dim];
	//doublereal *VS = new doublereal[Dim*Dim];//schur vectors U
	integer LWORK = 3*Dim;

	// CLAPACK 的 Schur decomposition
    //dgetrf_(&Dim,&Dim,Data,&Dim,ipiv,&INFO);
	JOBVS = 'V';
	SORT = 'S';//N 不重排特徵值 S要排 排了才會和DARE相同，讓T11在單位圓外
	dgees_(&JOBVS,&SORT,(L_fp) SELECT_eigen,&Dim,Data,&Dim,&SDIM,WR,WI,Data_U,&Dim,work,&LWORK,BWORK,&INFO);

	if (INFO != 0)
	{
		cout << "Fail to Calculate schur decomposition" << endl;
		system("pause");
	}
	// CLAPACK 的 矩陣 inverse
    /*dgetri_(&Dim,Data,&Dim,ipiv,work,&Dim,&INFO);

	if (INFO != 0)
	{
		cout << "Fail to Calculate Matrix Inversion" << endl;
		system("pause");
	}*/

	// 清除動態記憶體
	delete[] work;
	//delete[] ipiv;
	//delete[] VS;
	delete[] WR;
	delete[] WI;
	delete[] BWORK;

}

logical SELECT_eigen(double *ar,double *ai)//test eigenvalue if > 1.0 then move
{
	logical flag;
	//printf("Inside SELECT function!\n");
	//printf("%g, %g\n", *ar, *ai);
	if (abs(*ar) > 1.0)
	{
		flag=1;
	}
	else
	{
		flag=0;
	}
	return flag;//*ar < 0.0;//true;
}

void showMat(double* A,int Nrow, int Ncol, char * MatName )
{
	cout << "-------------------Data:"<< MatName <<"-------------------------- " << endl;

	for (int i = 0 ; i < Nrow ; i++)
	{
		for (int j = 0 ; j < Ncol ; j++)
		{

			if (A[j+i*Ncol] == 0)
				cout<< "\t" << "0" ;
			else if (fabs(A[j+i*Ncol]) < 0.0001)
				cout << "\t"<< "0.0";
			else
				cout << "\t"<< A[j+i*Ncol] ;
		}
		cout << endl;
	}

	cout << endl;
}
void Mat_weighted_mean( double* A, int MA, int NA,double WA, double* B, int MB, int NB,double WB, double* result)
{
	double w_A[3];
	MatScalarMul(A,MA*NA,WA,w_A);

	double w_B[3];
	MatScalarMul(B,MB*NB,WB,w_B);

	MatAddAB(w_A,w_B,result,3);
	MatScalarMul(result,3,1/(WA+WB),result);
}

void Mat_y_Ax_b( double*A, int MA, int NA, double *x, int Mx ,int Nx, double *b, double *result )
{
	MatMulAB(A,MA,NA,x,Mx,Nx,result);
	MatAddAB(result,b,result,MA*Nx);
}

void Mat_I_nxn( double* result, int n )
{
	for ( int i = 0 ; i < n*n ; i++ )
	{
		if ( (i%(n+1)) == 0)
			result[i] = 1; 
		else
			result[i] = 0;
	}
}

void Mat_O_mxn(double* result,int m, int n)
{
	int size = m*n;
	for ( int i =0 ;i<size;i++ )
	{
		result[i] = 0;
	}
}
void Mat_A_addcol_B( double* A, int MA, int NA, double* B, int MB, int NB, double* result )
{
	if (MA!=MB) 
	{
		printf("Wrong Matrix Dimension for Mat_A_addcol_B");
		system("pause");
	}
	int M = MA;
	int N = NA+NB;

	int i_A(0),i_B(0),k(0);

	for(int i =0 ;i<M*N;i++)
	{
		k = i%N;
		if (k < NA )
		{
			result[i] =  A[i_A];
			i_A++;
		} 
		else
		{ 
			result[i] =  B[i_B];
			i_B++;
		}
	}
	
}
void Mat_A_addrow_B(double* A, int MA, int NA, double* B, int MB, int NB, double* result)
{
	if (NA!=NB) 
	{
		printf("Wrong Matrix Dimension for Mat_A_addrow_B");
		system("pause");
	}
	int M = MA+MB;
	int N = NA;
	int i_A(0),i_B(0);

	for (int i = 0; i<M*N ; i++)
	{
		if (i<MA*NA)
		{
			result[i] = A[i_A];
			i_A++;
		} 
		else
		{
			result[i] = B[i_B];
			i_B++;
		}
	}

}

void showVar(double x,char * VarName)
{
	cout << "-------------------Data:"<<VarName<<"---------------------------" << endl;
	cout<<"\t"<<x<<endl;
	cout << endl;
}

void R_rpy( double r,double p, double y,double *result )
{
	double Rzy[9] = { cos(y),-sin(y),0,sin(y),cos(y),0,0,0,1 };
	double Ryp[9] = { cos(p),0 , sin(p),0,1,0, -sin(p),0,cos(p) };
	double Rxr[9] = { 1 ,0 ,0, 0 ,cos(r),-sin(r),0,sin(r),cos(r)};

	double tmp[9];
	MatMulAB(Rzy,3,3,Ryp,3,3,tmp);
	MatMulAB(tmp,3,3,Rxr,3,3,result);
}