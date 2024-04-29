#ifndef _lbmD3Q15_H_
#define _lbmD3Q15_H_

//this here is used for calculating the gradients
enum velocityDir_15{
                dV_15_ZERO_ZERO_ZERO, 
                dV_15_P1_ZERO_ZERO, dV_15_M1_ZERO_ZERO, dV_15_ZERO_P1_ZERO, dV_15_ZERO_M1_ZERO, dV_15_ZERO_ZERO_P1,dV_15_ZERO_ZERO_M1,
                dV_15_PH1_PH1_PH1, dV_15_PH1_MH1_PH1, dV_15_MH1_MH1_PH1, dV_15_MH1_PH1_PH1, dV_15_PH1_PH1_MH1, dV_15_PH1_MH1_MH1, dV_15_MH1_MH1_MH1, dV_15_MH1_PH1_MH1,
                };


int oppdV[15] =  {  
                    dV_15_ZERO_ZERO_ZERO,
                    dV_15_M1_ZERO_ZERO, dV_15_P1_ZERO_ZERO, dV_15_ZERO_M1_ZERO, dV_15_ZERO_P1_ZERO, dV_15_ZERO_ZERO_M1 ,dV_15_ZERO_ZERO_P1,
                    dV_15_MH1_MH1_MH1,dV_15_MH1_PH1_MH1,dV_15_PH1_PH1_MH1,dV_15_PH1_MH1_MH1,dV_15_MH1_MH1_PH1,dV_15_MH1_PH1_PH1,dV_15_PH1_PH1_PH1,dV_15_PH1_MH1_PH1,

                    };







template<typename T >
struct lbmD3Q15
{
    int dvN;

     T W[15]; 

     T Cx[15]; T CxC[15]; T CxF[15]; 
     T Cy[15]; T CyC[15]; T CyF[15]; 
     T Cz[15]; T CzC[15]; T CzF[15]; 

    T theta0, thetaInverse ;
    lbmD3Q15(T c); //constructor
    
};





//constructor
template<typename T>
lbmD3Q15<T>::lbmD3Q15(T c1)
    {
        dvN = 15;
        theta0 =(c1*c1) *1.0/6.0;
        thetaInverse = 1.0/theta0;


W[dV_15_ZERO_ZERO_ZERO ]    =  14.0/36.0;

W[dV_15_P1_ZERO_ZERO   ]    =  1.0/36.0;	            	
W[dV_15_M1_ZERO_ZERO   ]    =  1.0/36.0;	            	
W[dV_15_ZERO_P1_ZERO   ]    =  1.0/36.0;	            	
W[dV_15_ZERO_M1_ZERO   ]    =  1.0/36.0;	            	
W[dV_15_ZERO_ZERO_P1   ]    =  1.0/36.0;		        	
W[dV_15_ZERO_ZERO_M1   ]    =  1.0/36.0;

W[dV_15_PH1_PH1_PH1    ]    = 2.0/36.0 ;
W[dV_15_PH1_MH1_PH1    ]    = 2.0/36.0 ;
W[dV_15_MH1_MH1_PH1    ]    = 2.0/36.0 ;
W[dV_15_MH1_PH1_PH1    ]    = 2.0/36.0 ;
W[dV_15_PH1_PH1_MH1    ]    = 2.0/36.0 ;
W[dV_15_PH1_MH1_MH1    ]    = 2.0/36.0 ;
W[dV_15_MH1_MH1_MH1    ]    = 2.0/36.0 ;
W[dV_15_MH1_PH1_MH1    ]    = 2.0/36.0 ;

    
    
Cx[dV_15_ZERO_ZERO_ZERO ]    =0.0;         Cy[dV_15_ZERO_ZERO_ZERO ]    =0.0;             Cz[dV_15_ZERO_ZERO_ZERO ]    =0.0;

Cx[dV_15_P1_ZERO_ZERO   ]    = c1;         Cy[dV_15_P1_ZERO_ZERO   ]    =0.0;             Cz[dV_15_P1_ZERO_ZERO   ]    =0.0;
Cx[dV_15_M1_ZERO_ZERO   ]    =-c1;         Cy[dV_15_M1_ZERO_ZERO   ]    =0.0;             Cz[dV_15_M1_ZERO_ZERO   ]    =0.0;
Cx[dV_15_ZERO_P1_ZERO   ]    =0.0;         Cy[dV_15_ZERO_P1_ZERO   ]    = c1;             Cz[dV_15_ZERO_P1_ZERO   ]    =0.0;
Cx[dV_15_ZERO_M1_ZERO   ]    =0.0;         Cy[dV_15_ZERO_M1_ZERO   ]    =-c1;             Cz[dV_15_ZERO_M1_ZERO   ]    =0.0;
Cx[dV_15_ZERO_ZERO_P1   ]    =0.0;         Cy[dV_15_ZERO_ZERO_P1   ]    =0.0;             Cz[dV_15_ZERO_ZERO_P1   ]    = c1;
Cx[dV_15_ZERO_ZERO_M1   ]    =0.0;         Cy[dV_15_ZERO_ZERO_M1   ]    =0.0;             Cz[dV_15_ZERO_ZERO_M1   ]    =-c1;


Cx[dV_15_PH1_PH1_PH1   ]    = 0.5*c1;              Cy[dV_15_PH1_PH1_PH1    ]    = 0.5*c1;             Cz[dV_15_PH1_PH1_PH1   ]    = 0.5*c1;
Cx[dV_15_PH1_MH1_PH1   ]    = 0.5*c1;              Cy[dV_15_PH1_MH1_PH1    ]    =-0.5*c1;             Cz[dV_15_PH1_MH1_PH1   ]    = 0.5*c1;
Cx[dV_15_MH1_MH1_PH1   ]    =-0.5*c1;              Cy[dV_15_MH1_MH1_PH1    ]    =-0.5*c1;             Cz[dV_15_MH1_MH1_PH1   ]    = 0.5*c1;
Cx[dV_15_MH1_PH1_PH1   ]    =-0.5*c1;              Cy[dV_15_MH1_PH1_PH1    ]    = 0.5*c1;             Cz[dV_15_MH1_PH1_PH1   ]    = 0.5*c1;
Cx[dV_15_PH1_PH1_MH1   ]    = 0.5*c1;              Cy[dV_15_PH1_PH1_MH1    ]    = 0.5*c1;             Cz[dV_15_PH1_PH1_MH1   ]    =-0.5*c1;
Cx[dV_15_PH1_MH1_MH1   ]    = 0.5*c1;              Cy[dV_15_PH1_MH1_MH1    ]    =-0.5*c1;             Cz[dV_15_PH1_MH1_MH1   ]    =-0.5*c1;
Cx[dV_15_MH1_MH1_MH1   ]    =-0.5*c1;              Cy[dV_15_MH1_MH1_MH1    ]    =-0.5*c1;             Cz[dV_15_MH1_MH1_MH1   ]    =-0.5*c1;
Cx[dV_15_MH1_PH1_MH1   ]    =-0.5*c1;              Cy[dV_15_MH1_PH1_MH1    ]    = 0.5*c1;             Cz[dV_15_MH1_PH1_MH1   ]    =-0.5*c1;
    
    
    
    
CxC[dV_15_PH1_PH1_PH1   ]     =ceil( 0.5*c1);             CyC[dV_15_PH1_PH1_PH1    ]    =ceil( 0.5*c1);             CzC[dV_15_PH1_PH1_PH1    ]    =ceil( 0.5*c1);
CxC[dV_15_PH1_MH1_PH1   ]     =ceil( 0.5*c1);             CyC[dV_15_PH1_MH1_PH1    ]    =ceil(-0.5*c1);             CzC[dV_15_PH1_MH1_PH1    ]    =ceil( 0.5*c1);
CxC[dV_15_MH1_MH1_PH1   ]     =ceil(-0.5*c1);             CyC[dV_15_MH1_MH1_PH1    ]    =ceil(-0.5*c1);             CzC[dV_15_MH1_MH1_PH1    ]    =ceil( 0.5*c1);
CxC[dV_15_MH1_PH1_PH1   ]     =ceil(-0.5*c1);             CyC[dV_15_MH1_PH1_PH1    ]    =ceil( 0.5*c1);             CzC[dV_15_MH1_PH1_PH1    ]    =ceil( 0.5*c1);
CxC[dV_15_PH1_PH1_MH1   ]     =ceil( 0.5*c1);             CyC[dV_15_PH1_PH1_MH1    ]    =ceil( 0.5*c1);             CzC[dV_15_PH1_PH1_MH1    ]    =ceil(-0.5*c1);
CxC[dV_15_PH1_MH1_MH1   ]     =ceil( 0.5*c1);             CyC[dV_15_PH1_MH1_MH1    ]    =ceil(-0.5*c1);             CzC[dV_15_PH1_MH1_MH1    ]    =ceil(-0.5*c1);
CxC[dV_15_MH1_MH1_MH1   ]     =ceil(-0.5*c1);             CyC[dV_15_MH1_MH1_MH1    ]    =ceil(-0.5*c1);             CzC[dV_15_MH1_MH1_MH1    ]    =ceil(-0.5*c1);
CxC[dV_15_MH1_PH1_MH1   ]     =ceil(-0.5*c1);             CyC[dV_15_MH1_PH1_MH1    ]    =ceil( 0.5*c1);             CzC[dV_15_MH1_PH1_MH1    ]    =ceil(-0.5*c1);
    
    
    
CxF[dV_15_PH1_PH1_PH1   ]     =floor( 0.5*c1);             CyF[dV_15_PH1_PH1_PH1    ]    =floor( 0.5*c1);             CzF[dV_15_PH1_PH1_PH1    ]    =floor( 0.5*c1);
CxF[dV_15_PH1_MH1_PH1   ]     =floor( 0.5*c1);             CyF[dV_15_PH1_MH1_PH1    ]    =floor(-0.5*c1);             CzF[dV_15_PH1_MH1_PH1    ]    =floor( 0.5*c1);
CxF[dV_15_MH1_MH1_PH1   ]     =floor(-0.5*c1);             CyF[dV_15_MH1_MH1_PH1    ]    =floor(-0.5*c1);             CzF[dV_15_MH1_MH1_PH1    ]    =floor( 0.5*c1);
CxF[dV_15_MH1_PH1_PH1   ]     =floor(-0.5*c1);             CyF[dV_15_MH1_PH1_PH1    ]    =floor( 0.5*c1);             CzF[dV_15_MH1_PH1_PH1    ]    =floor( 0.5*c1);
CxF[dV_15_PH1_PH1_MH1   ]     =floor( 0.5*c1);             CyF[dV_15_PH1_PH1_MH1    ]    =floor( 0.5*c1);             CzF[dV_15_PH1_PH1_MH1    ]    =floor(-0.5*c1);
CxF[dV_15_PH1_MH1_MH1   ]     =floor( 0.5*c1);             CyF[dV_15_PH1_MH1_MH1    ]    =floor(-0.5*c1);             CzF[dV_15_PH1_MH1_MH1    ]    =floor(-0.5*c1);
CxF[dV_15_MH1_MH1_MH1   ]     =floor(-0.5*c1);             CyF[dV_15_MH1_MH1_MH1    ]    =floor(-0.5*c1);             CzF[dV_15_MH1_MH1_MH1    ]    =floor(-0.5*c1);
CxF[dV_15_MH1_PH1_MH1   ]     =floor(-0.5*c1);             CyF[dV_15_MH1_PH1_MH1    ]    =floor( 0.5*c1);             CzF[dV_15_MH1_PH1_MH1    ]    =floor(-0.5*c1);
    
    
    
    
    
}



#endif

