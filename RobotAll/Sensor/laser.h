#include "UrgLaser.h"
#include <iostream>
#include "string.h"
#include <math.h>

#include <cmath>

#include <fstream>

using namespace std;

#pragma once
class laser
{
public:
	int *ya; 
	float *x_data;
	float *y_data;

	int n;
	int j;
	float sample;
	float pi ;
	UrgLaser *urg;
	laser(void);
	~laser(void);

    void transfer(void);
};

