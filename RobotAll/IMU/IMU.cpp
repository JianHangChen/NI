#include "stdafx.h"
#include "IMU.h"
#include "math.h"



IMU ::IMU (void){
 
	absangle  = 0 ;
	
	calicount  = 1000 ;   //cali的筆數
	
	//adamssimcount = 0;

	Ts =0.005;             //IMU sampling freqency 目前為 200 HZ

	microstraincalix = 0;//0.00847 ; //cali值  與 microstrain MIP 比較    直接平移一個量

	microstraincaliy = 0;//0.0171  ; 
	
	



	anglexbias = 0 ;
	angleybias = 0 ; 
	absxbias  = 0 ;
	absybias  = 0 ;


	abscalix = 0 ;
	abscaliy = 0 ;
	anglecalix = 0 ;  
	anglecaliy = 0 ;  
	anglecaliz = 0 ;  

	
	for(int i = 0; i<10000 ; i++)
	{
	
	anglex[i]=0;    // vel 積分
	angley[i]=0;
	
	anglefilterx[i]=0;
	anglefiltery[i]=0;

	anglez[i]=0;
	
	filterx [i]=0;
	filtery [i]=0;
	
	filterabsx [i]=0;
	filterabsy [i]=0;
		
	absanglex [i]=0 ;
	absangley [i]=0 ;
	
	
	finalanglex[i]=0;  //complementary 算出  
	finalangley[i]=0;   
	finalanglez[i]=0;
	
	finalanglex_test[i]=0;

	velx [i]=0;	
	vely [i]=0;	
	velz [i]=0;	
	
	accelx[i]=0;
	accely[i]=0;
	accelz[i]=0;
	
	filteraccelx[i]=0;
	filteraccely[i]=0;
	filteraccelz[i]=0;
	
	
	
	
	
	
	
	
	
	
	
	}






















	
	IMU_Lock=false;
}
	
void  IMU :: integation(float *velx ,float *anglex , float *vely ,float *angley ,float *velz ,float *anglez ){
                   
				   
				   velbiasx1  = (*velx)   -velbiasx;  
			       velbiasx2  = *(velx-1) -velbiasx;
			   
				   velbiasy1  = (*vely) -velbiasy;  
			       velbiasy2  = *(vely-1) -velbiasy;

				   velbiasz1  = (*velz) -velbiasz;  
			       velbiasz2  = *(velz-1) -velbiasz;   // vel  bias 

				    
		      *anglex=  *(anglex-1)  +  0.5*Ts* (  velbiasx1) +   0.5*Ts* (  velbiasx2 )    ; 
              *angley=  *(angley-1)  +  0.5*Ts* (  velbiasy1) +   0.5*Ts* (  velbiasy2 )    ; 
			  *anglez=  *(anglez-1)  +  0.5*Ts* (  velbiasz1) +   0.5*Ts* (  velbiasz2 )  ; 
  }


void SetDataRate(int portNum)
{
	int status = 0, Value = 0, Ival =0;
	int chOption = 0;
	int rDataRate = 0;
	int mDataRate = 0;
	int sDataRate = 0;
	int valCheck = 0;
	BYTE Record[20];
	char AddA = 0xFC;
	char AddB = 0xA2;
    unsigned short nvalS  = 0;
	char Jjunk[20];
	purge_port(portNum);
	Value == 0 ;
	
	
	WriteDataRate(portNum, 200, &Record[0]);

};

void IMU :: complementaryfilter(float   * ang ,  float * abs , float * filterang , float *filterabs , float *filter , float bias  ){

	

	* ang =  *ang - bias ;  //扣掉相對角度的bias



	* filterang=  0.81425 * * ang -3.257 *  * (ang-1) + 4.8855 * * (ang-2) -3.257 *  * (ang-3) + 0.81425 ** (ang-4)
		- ( -3.5897* * (filterang-1) + 4.8513 * * (filterang-2) -2.9241** (filterang-3) +0.66301 ** (filterang-4)) ;



	* filterabs =  3.1239e-005 * * abs + 0.00012496 *  * (abs-1) + 0.00018743 * * (abs-2) +0.00012496 *  * (abs-3) + 3.1239e-005 ** (abs-4)
		- ( -3.5897 * * (filterabs-1) + 4.8513 * * (filterabs-2) -2.9241 ** (filterabs-3) +0.66301 ** (filterabs-4)) ; 



	*filter   =  * filterang + * filterabs ; 


}



void IMU :: ReadContinuousData(int portNum)
{
	
	BOOL bStopContinuous = FALSE;
	DWORD dwcharsRead = 0;
	I3dmgx3Set Record;
	BYTE cmd_return;
	int Curs_posY = 0;
	int Curs_posX = 0;
	int status    = 0;
	int error_record = 0;
	char consoleBuff[60] = {0};
	long valid_rec = 0;
	unsigned char error_cmd;
	
	i3dmgx3_openPort(portNum, 115200, 8, 0, 1, 1024, 1024);//開port
	SetDataRate(portNum) ;//直接設定sampling rate 目前為200 Hz 
	
	
	//put the node in continuous mode
	status =SetContinuousMode(portNum, 0xC2);
	
	

	//printf("_________________________________  angle_________________________________\n");
	//printf("          ROLL                   PITCH                         COUNT         \n");
	//printf("\n");
	//
	//getConXY(&Curs_posX, &Curs_posY); 
	//printf("\n\n\n\n\n");
	//
	

	//continue until the user hits the s key
	while(!bStopContinuous)
	{
		if(ReadNextRecord(portNum, &Record, &cmd_return) != SUCCESS)
			error_record++;

			if (cmd_return == 0xC2){
			   //move to the acceleration position and print the data
		         IMU_Lock=true;
				 		 
				 
				 
				  velx[count] =  Record.setB[0] ;//     單位rad/s

				  vely[count] =  Record.setB[1] ;

		          velz[count] =  Record.setB[2] ;	
				               
				  
				  accelx[count] =  Record.setA[1] ;//   單位m/s^2

				  accely[count] =  Record.setA[0] ;

		          accelz[count] =  Record.setA[2] ;	
		             
				



		          //絕對角度  由accel分量算出
				  absanglex[count]   = - asin(   accelx[count] /9.8) -  abscalix;
				  absangley[count]   =   asin(   accely[count] /9.8) -  abscaliy;  
				  
					
				  //相對角度的bias累積 在cali前
				  anglebiasx    =  anglebiasx  +   anglex[count-1] ;
				  anglebiasy    =  anglebiasy  +   angley[count-1] ;
				  anglebiasz    =  anglebiasz  +   angley[count-1] ;
              
				  
				  //絕對角度的bias累積  在cali前
				  absxbias =  absxbias  +  absanglex[count] ; 
			      absybias =  absybias  +  absangley [count]; 





   //當count = calicount時  將累積bias填入buffer內 						 
   if (count ==0.5*calicount){    
	                          
	     
		  if (absangle == 1){ //使用絕對角度 ==> 不cali

	    abscalix = 0; 
	    abscaliy = 0;
		                    }
	    
		  else {  
		abscalix =  absxbias /(0.5*calicount);   //使用相對角度 ==> cali
	    abscaliy =  absybias /(0.5*calicount);

		       }





							  
		for (int i = 0 ; i<(0.5*calicount);i++){  //將vel累積誤差bias消去
		velbiasx  = velbiasx + velx[i] ;
		velbiasy  = velbiasy + vely[i] ;
		velbiasz  = velbiasz + velz[i] ;
																				
		}
							
	    velbiasx = velbiasx /(calicount*0.5);
	    velbiasy = velbiasy /(calicount*0.5);
	    velbiasz = velbiasz /(calicount*0.5);
        }	



 

   if (count == calicount )//先讓vel bias 消除   往後0.5cali count 再把角度bias 減掉

   
   
   {
   anglecalix  =  anglebiasx /(0.5*calicount);
   anglecaliy   = anglebiasy /(0.5*calicount);
   anglecaliz   = anglebiasz /(0.5*calicount); 
   }




   if (count > 0.5*calicount){

	   integation (velx +count,anglex +count,vely +count,angley +count,velz +count,anglez+count);

	   
	   }	
   
   
   ;					 



   //finalanglex (相對角度)   absanglex (絕對角度) 


   complementaryfilter(  anglex +count  ,  absanglex+count  ,filterx +count, filterabsx +count ,  anglefilterx+count , anglecalix);
   complementaryfilter(  angley +count  ,  absangley+count  ,filtery +count, filterabsy +count ,  anglefiltery+count , anglecaliy);        


   finalanglex[count]  = anglefilterx[count] + microstraincalix ;  //加入 cali值  與 microstrain MIP 比較

   finalangley[count]  = anglefiltery[count] + microstraincaliy ;



   //print 角度 (度度制) cali完後就不顯示 比較不占資源
   


	   //sprintf(consoleBuff, "\t%f\t\t%f\t\t\t%d\n", finalanglex[count]*180/3.14,finalangley[count]*180/3.14 , count );
	   //setConXY(Curs_posX, Curs_posY ,  &consoleBuff[0]);
   




   valid_rec++;
   count++; 
   
   		         IMU_Lock=false;

			}

	}

}


void IMU ::adamssim(int gFlagSimulation,int adamssimcount){
	if (gFlagSimulation ==1 ) //adams sim
	{     
		///motion control thread begins in 1000counts
		  
		   


		          //adams讀入  accel (mm/s)     angvel(rad/s)
            /*     accelx[count] = accelx[count]*0.001;
				 accely[count] = accely[count]*0.001;*/
						



		          //絕對角度  由accel分量算出
				  absanglex[adamssimcount]   =   asin(   0.001*filteraccelx[adamssimcount] /9.9) -  abscalix;
				  absangley[adamssimcount]   =   asin(   0.001*filteraccely[adamssimcount] /9.9) -  abscaliy;  
				  
					
				  //相對角度的bias累積 在cali前
				  anglebiasx    =  anglebiasx  +   anglex[adamssimcount-1] ;
				  anglebiasy    =  anglebiasy  +   angley[adamssimcount-1] ;
				  anglebiasz    =  anglebiasz  +   angley[adamssimcount-1] ;
              
				  
				  //絕對角度的bias累積  在cali前
				  absxbias =  absxbias  +  absanglex[adamssimcount-1] ; 
			      absybias =  absybias  +  absangley [adamssimcount-1]; 





   //當count = calicount時  將累積bias填入buffer內 						 
   if (adamssimcount ==0.5*calicount){    
	                          
	     
		  if (absangle == 1){ //使用絕對角度 ==> 不cali

	    abscalix = 0; 
	    abscaliy = 0;
		                    }
	    
		  else {  
		abscalix =  absxbias /(0.5*calicount);   //使用相對角度 ==> cali
	    abscaliy =  absybias /(0.5*calicount);

		       }





							  
		for (int i = 0 ; i<(0.5*calicount);i++){  //將vel累積誤差bias消去
		velbiasx  = velbiasx + velx[i] ;
		velbiasy  = velbiasy + vely[i] ;
		velbiasz  = velbiasz + velz[i] ;
																				
		}
							
	    velbiasx = velbiasx /(calicount*0.5);
	    velbiasy = velbiasy /(calicount*0.5);
	    velbiasz = velbiasz /(calicount*0.5);
        }	



 

   if (adamssimcount == calicount )//先讓vel bias 消除   往後0.5cali count 再把角度bias 減掉

   
   
   {
   anglecalix  =  anglebiasx /(0.5*calicount);
   anglecaliy   = anglebiasy /(0.5*calicount);
   anglecaliz   = anglebiasz /(0.5*calicount); 
   }




   if (adamssimcount > 0.5*calicount){

	   integation (velx +adamssimcount,anglex +adamssimcount,vely +adamssimcount,angley +adamssimcount,velz +adamssimcount,anglez+adamssimcount);

	   
	   } ;					 



   //finalanglex (相對角度)   absanglex (絕對角度) 


   complementaryfilter(  anglex +adamssimcount  ,  absanglex+adamssimcount  ,filterx +adamssimcount, filterabsx +adamssimcount ,  anglefilterx+adamssimcount , anglecalix);
   complementaryfilter(  angley +adamssimcount  ,  absangley+adamssimcount  ,filtery +adamssimcount, filterabsy +adamssimcount ,  anglefiltery+adamssimcount , anglecaliy);        




   finalanglex[adamssimcount]  = anglefilterx[adamssimcount] + microstraincalix ;  //加入 cali值  與 microstrain MIP 比較

   finalangley[adamssimcount]  = anglefiltery[adamssimcount] + microstraincaliy ;      

   
   
   //0418
   // IMU x sagittal direction(roll)    y lateral direction(pitch) 
   // 需乘上 rotation matrix 將 IMU accel 在 local frame 轉到 world frame 
   
   //ax = ax*cos(pitch) + ay*sin(pitch) sin(roll)  + sin(pitch) sin(roll)
   //ay = cos(roll)*ay - sin(roll)
   //az = 1

   //IMU_roll  =   finalanglex[count] ;
   //IMU_pitch =   finalangley[count] ;



   //accelx[count] =   accelx[count] * cos( IMU_pitch) + accely[count]*sin( IMU_pitch)*sin(IMU_roll)+sin(IMU_pitch)*sin(IMU_roll)*1;
   //accely[count] =   accely[count] * cos(IMU_roll) -sin(IMU_roll)*1;
    
   //cout<<IMU_roll<<" " << IMU_pitch << endl;

  
   
   //0418



    adamssimcount++; 
   
 }

}








IMU ::~IMU (void){} ;

