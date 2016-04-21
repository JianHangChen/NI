
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include "plane.h"
#include <vector>

#ifndef estimate_H
#define estimate_H
using namespace std;
using namespace Kalman;

class estimate {

public:
    
	estimate(); // 建構子
	~estimate (); // 解構子
	int compute (double cog ,double  *zmp , double cogaccel,int estimatecount,int direction) ;
	
	vector <double>  Cogstatelateral;
	vector <double>  Dcogstatelateral;
	vector <double>  zmpstatelateral;
		
	vector <double>  Cogstatesaggital;
	vector <double>  Dcogstatesaggital;
	vector <double>  zmpstatesaggital;
	
	//double COGsaggital;
	//double DCOGsaggital;
	//double ZMPsaggital;
	//double COGsaggital_before;
	//double DCOGsaggital_before;
	//double ZMPsaggital_before;





};


#endif
