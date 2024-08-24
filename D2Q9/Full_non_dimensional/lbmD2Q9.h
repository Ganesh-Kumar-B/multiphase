#ifndef _lbmD2Q9_H_
#define _lbmD2Q9_H_

enum velocityDir{ dV_ZERO_ZERO, dV_P1_ZERO, dV_ZERO_P1, dV_M1_ZERO, dV_ZERO_M1, dV_P1_P1, dV_M1_P1, dV_M1_M1, dV_P1_M1};

int oppdV[9] =  {dV_ZERO_ZERO,  dV_M1_ZERO, dV_ZERO_M1, dV_P1_ZERO, dV_ZERO_P1, dV_M1_M1, dV_P1_M1, dV_P1_P1, dV_M1_P1};





template<typename T >
struct lbmD2Q9
{
    int dvN;

     T W[9]; 
     T Cx[9]; 
     T Cy[9]; 
     T Cz[9]; 

    T theta0, thetaInverse ;
    lbmD2Q9(T c); //constructor
    
};





//constructor
template<typename T>
lbmD2Q9<T>::lbmD2Q9(T c1)
{
    dvN = 9;
    theta0 = (c1*c1)/3.0 ;
    thetaInverse = 1.0/theta0;


        W[dV_ZERO_ZERO]=16.0/36.0;	         Cx[dV_ZERO_ZERO]=0.0;	     Cy[dV_ZERO_ZERO]=0.0;	      
        W[dV_P1_ZERO]  =4.0/36.0;	         Cx[dV_P1_ZERO]  = c1;	     Cy[dV_P1_ZERO]  =0.0;	      
        W[dV_ZERO_P1]  =4.0/36.0;	         Cx[dV_ZERO_P1]  =0.0;	     Cy[dV_ZERO_P1]  = c1;	      
        W[dV_M1_ZERO]  =4.0/36.0;	         Cx[dV_M1_ZERO]  =-c1;	     Cy[dV_M1_ZERO]  =0.0;	      
        W[dV_ZERO_M1]  =4.0/36.0;	         Cx[dV_ZERO_M1]  =0.0;	     Cy[dV_ZERO_M1]  =-c1;	      
        W[dV_P1_P1]    =1.0/36.0;		     Cx[dV_P1_P1]    = c1;	     Cy[dV_P1_P1]    = c1;	      
        W[dV_M1_P1]    =1.0/36.0;		     Cx[dV_M1_P1]    =-c1;	     Cy[dV_M1_P1]    = c1;	      
        W[dV_M1_M1]    =1.0/36.0;		     Cx[dV_M1_M1]    =-c1;	     Cy[dV_M1_M1]    =-c1;	      
        W[dV_P1_M1]    =1.0/36.0;		     Cx[dV_P1_M1]    = c1;	     Cy[dV_P1_M1]    =-c1;	      


}



#endif

