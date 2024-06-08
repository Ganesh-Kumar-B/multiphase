#pragma once
#include<math.h>
enum velocityDir{dV_ZERO_ZERO_ZERO, 

                dV_P2_ZERO_ZERO, dV_M2_ZERO_ZERO, dV_ZERO_P2_ZERO, dV_ZERO_M2_ZERO, dV_ZERO_ZERO_P2,dV_ZERO_ZERO_M2,
                dV_P1_ZERO_ZERO, dV_M1_ZERO_ZERO, dV_ZERO_P1_ZERO, dV_ZERO_M1_ZERO, dV_ZERO_ZERO_P1,dV_ZERO_ZERO_M1,

                dV_P1_P1_ZERO, dV_M1_P1_ZERO ,dV_M1_M1_ZERO,dV_P1_M1_ZERO , 
                dV_P1_ZERO_P1,dV_M1_ZERO_P1,dV_M1_ZERO_M1,dV_P1_ZERO_M1,
                dV_ZERO_P1_P1,dV_ZERO_M1_P1,dV_ZERO_M1_M1,dV_ZERO_P1_M1,
                


                
                dV_P1_P1_P1, dV_P1_M1_P1, dV_M1_M1_P1, dV_M1_P1_P1, dV_P1_P1_M1, dV_P1_M1_M1, dV_M1_M1_M1, dV_M1_P1_M1,
                dV_PH1_PH1_PH1, dV_PH1_MH1_PH1, dV_MH1_MH1_PH1, dV_MH1_PH1_PH1, dV_PH1_PH1_MH1, dV_PH1_MH1_MH1, dV_MH1_MH1_MH1, dV_MH1_PH1_MH1,


                };

int oppdV[41] =  {  dV_ZERO_ZERO_ZERO,
                    dV_M2_ZERO_ZERO, dV_P2_ZERO_ZERO, dV_ZERO_M2_ZERO, dV_ZERO_P2_ZERO, dV_ZERO_ZERO_M2 ,dV_ZERO_ZERO_P2,
                    dV_M1_ZERO_ZERO, dV_P1_ZERO_ZERO, dV_ZERO_M1_ZERO, dV_ZERO_P1_ZERO, dV_ZERO_ZERO_M1 ,dV_ZERO_ZERO_P1,

                    dV_M1_M1_ZERO,dV_P1_M1_ZERO,dV_P1_P1_ZERO,dV_M1_P1_ZERO,
                    dV_M1_ZERO_M1,dV_P1_ZERO_M1,dV_P1_ZERO_P1,dV_M1_ZERO_P1,
                    dV_ZERO_M1_M1,dV_ZERO_P1_M1,dV_ZERO_P1_P1,dV_ZERO_M1_P1,

                    
                    dV_M1_M1_M1,dV_M1_P1_M1,dV_P1_P1_M1,dV_P1_M1_M1,dV_M1_M1_P1,dV_M1_P1_P1,dV_P1_P1_P1,dV_P1_M1_P1,
                    dV_MH1_MH1_MH1,dV_MH1_PH1_MH1,dV_PH1_PH1_MH1,dV_PH1_MH1_MH1,dV_MH1_MH1_PH1,dV_MH1_PH1_PH1,dV_PH1_PH1_PH1,dV_PH1_MH1_PH1,
                    };







template<typename T >
struct lbmD3Q41
{
    int dvN;

     T W[41]; 
     T Cx[41];  T CxC[41]; T CxF[41]; 
     T Cy[41];  T CyC[41]; T CyF[41]; 
     T Cz[41];  T CzC[41]; T CzF[41]; 

    T theta0, thetaInverse ;
    lbmD3Q41(T c); //constructor
    
};





//constructor
template<typename T>
lbmD3Q41<T>::lbmD3Q41(T c1)
    {   
        dvN = 41;
        theta0 = 0.2948964908710636;
        thetaInverse = 1.0/theta0;


W[dV_ZERO_ZERO_ZERO ]    = 0.1975697820320          ;

W[dV_P1_ZERO_ZERO   ]    = 0.04743040745116578      ;      	
W[dV_M1_ZERO_ZERO   ]    = 0.04743040745116578      ;      	
W[dV_ZERO_P1_ZERO   ]    = 0.04743040745116578      ;      	
W[dV_ZERO_M1_ZERO   ]    = 0.04743040745116578      ;      	
W[dV_ZERO_ZERO_P1   ]    = 0.04743040745116578      ;    	
W[dV_ZERO_ZERO_M1   ]    = 0.04743040745116578      ;

W[dV_P2_ZERO_ZERO   ]    = 0.0016568766450157603    ;      	
W[dV_M2_ZERO_ZERO   ]    = 0.0016568766450157603    ;      	
W[dV_ZERO_P2_ZERO   ]    = 0.0016568766450157603    ;      	
W[dV_ZERO_M2_ZERO   ]    = 0.0016568766450157603    ;      	
W[dV_ZERO_ZERO_P2   ]    = 0.0016568766450157603    ;    	
W[dV_ZERO_ZERO_M2   ]    = 0.0016568766450157603    ;


W[dV_P1_P1_ZERO     ]    =	0.006511753278324646    ;     	
W[dV_M1_M1_ZERO     ]    =	0.006511753278324646    ;     	
W[dV_ZERO_P1_P1     ]    =  0.006511753278324646    ;       
W[dV_ZERO_M1_M1     ]    =  0.006511753278324646    ;       
W[dV_P1_ZERO_P1     ]    =  0.006511753278324646    ;       
W[dV_M1_ZERO_M1     ]    =  0.006511753278324646    ;       
W[dV_P1_M1_ZERO     ]    =  0.006511753278324646    ;       
W[dV_M1_P1_ZERO     ]    =  0.006511753278324646    ;       
W[dV_ZERO_P1_M1     ]    =  0.006511753278324646    ;       
W[dV_ZERO_M1_P1     ]    =  0.006511753278324646    ;   	
W[dV_P1_ZERO_M1     ]    =  0.006511753278324646    ;   	
W[dV_M1_ZERO_P1     ]    =  0.006511753278324646    ;     	


W[dV_P1_P1_P1    ]    =  0.004540878011544409       ;
W[dV_P1_M1_P1    ]    =  0.004540878011544409       ;
W[dV_M1_M1_P1    ]    =  0.004540878011544409       ;
W[dV_M1_P1_P1    ]    =  0.004540878011544409       ;
W[dV_P1_P1_M1    ]    =  0.004540878011544409       ;
W[dV_P1_M1_M1    ]    =  0.004540878011544409       ;
W[dV_M1_M1_M1    ]    =  0.004540878011544409       ;
W[dV_M1_P1_M1    ]    =  0.004540878011544409       ;


W[dV_PH1_PH1_PH1    ]    =    0.04917980624482672  ;
W[dV_PH1_MH1_PH1    ]    =    0.04917980624482672  ;
W[dV_MH1_MH1_PH1    ]    =    0.04917980624482672  ;
W[dV_MH1_PH1_PH1    ]    =    0.04917980624482672  ;
W[dV_PH1_PH1_MH1    ]    =    0.04917980624482672  ;
W[dV_PH1_MH1_MH1    ]    =    0.04917980624482672  ;
W[dV_MH1_MH1_MH1    ]    =    0.04917980624482672  ;
W[dV_MH1_PH1_MH1    ]    =    0.04917980624482672  ;



    
Cx[dV_ZERO_ZERO_ZERO ]    =    0.0;             Cy[dV_ZERO_ZERO_ZERO ]    =    0.0;             Cz[dV_ZERO_ZERO_ZERO ]    =  0.0;
Cx[dV_P1_ZERO_ZERO   ]    = 1.0*c1;             Cy[dV_P1_ZERO_ZERO   ]    =    0.0;             Cz[dV_P1_ZERO_ZERO   ]    =  0.0;
Cx[dV_M1_ZERO_ZERO   ]    =-1.0*c1;             Cy[dV_M1_ZERO_ZERO   ]    =    0.0;             Cz[dV_M1_ZERO_ZERO   ]    =  0.0;
Cx[dV_ZERO_P1_ZERO   ]    =    0.0;             Cy[dV_ZERO_P1_ZERO   ]    = 1.0*c1;             Cz[dV_ZERO_P1_ZERO   ]    =  0.0;
Cx[dV_ZERO_M1_ZERO   ]    =    0.0;             Cy[dV_ZERO_M1_ZERO   ]    =-1.0*c1;             Cz[dV_ZERO_M1_ZERO   ]    =  0.0;
Cx[dV_ZERO_ZERO_P1   ]    =    0.0;             Cy[dV_ZERO_ZERO_P1   ]    =    0.0;             Cz[dV_ZERO_ZERO_P1   ]    = 1.0*c1;
Cx[dV_ZERO_ZERO_M1   ]    =    0.0;             Cy[dV_ZERO_ZERO_M1   ]    =    0.0;             Cz[dV_ZERO_ZERO_M1   ]    =-1.0*c1;

Cx[dV_P2_ZERO_ZERO   ]    = 2.0*c1;             Cy[dV_P2_ZERO_ZERO   ]    =    0.0;             Cz[dV_P2_ZERO_ZERO   ]    =  0.0;
Cx[dV_M2_ZERO_ZERO   ]    =-2.0*c1;             Cy[dV_M2_ZERO_ZERO   ]    =    0.0;             Cz[dV_M2_ZERO_ZERO   ]    =  0.0;
Cx[dV_ZERO_P2_ZERO   ]    =    0.0;             Cy[dV_ZERO_P2_ZERO   ]    = 2.0*c1;             Cz[dV_ZERO_P2_ZERO   ]    =  0.0;
Cx[dV_ZERO_M2_ZERO   ]    =    0.0;             Cy[dV_ZERO_M2_ZERO   ]    =-2.0*c1;             Cz[dV_ZERO_M2_ZERO   ]    =  0.0;
Cx[dV_ZERO_ZERO_P2   ]    =    0.0;             Cy[dV_ZERO_ZERO_P2   ]    =    0.0;             Cz[dV_ZERO_ZERO_P2   ]    = 2.0*c1;
Cx[dV_ZERO_ZERO_M2   ]    =    0.0;             Cy[dV_ZERO_ZERO_M2   ]    =    0.0;             Cz[dV_ZERO_ZERO_M2   ]    =-2.0*c1;

Cx[dV_P1_P1_ZERO     ]    = 1.0*c1;             Cy[dV_P1_P1_ZERO     ]    = 1.0*c1;             Cz[dV_P1_P1_ZERO     ]    = 0.0*c1;
Cx[dV_M1_M1_ZERO     ]    =-1.0*c1;             Cy[dV_M1_M1_ZERO     ]    =-1.0*c1;             Cz[dV_M1_M1_ZERO     ]    = 0.0*c1;
Cx[dV_ZERO_P1_P1     ]    = 0.0*c1;             Cy[dV_ZERO_P1_P1     ]    = 1.0*c1;             Cz[dV_ZERO_P1_P1     ]    = 1.0*c1; 
Cx[dV_ZERO_M1_M1     ]    = 0.0*c1;             Cy[dV_ZERO_M1_M1     ]    =-1.0*c1;             Cz[dV_ZERO_M1_M1     ]    =-1.0*c1;
Cx[dV_P1_ZERO_P1     ]    = 1.0*c1;             Cy[dV_P1_ZERO_P1     ]    = 0.0*c1;             Cz[dV_P1_ZERO_P1     ]    = 1.0*c1;
Cx[dV_M1_ZERO_M1     ]    =-1.0*c1;             Cy[dV_M1_ZERO_M1     ]    = 0.0*c1;             Cz[dV_M1_ZERO_M1     ]    =-1.0*c1;
Cx[dV_P1_M1_ZERO     ]    = 1.0*c1;             Cy[dV_P1_M1_ZERO     ]    =-1.0*c1;             Cz[dV_P1_M1_ZERO     ]    = 0.0*c1;
Cx[dV_M1_P1_ZERO     ]    =-1.0*c1;             Cy[dV_M1_P1_ZERO     ]    = 1.0*c1;             Cz[dV_M1_P1_ZERO     ]    = 0.0*c1;
Cx[dV_ZERO_P1_M1     ]    = 0.0*c1;             Cy[dV_ZERO_P1_M1     ]    = 1.0*c1;             Cz[dV_ZERO_P1_M1     ]    =-1.0*c1;
Cx[dV_ZERO_M1_P1     ]    = 0.0*c1;             Cy[dV_ZERO_M1_P1     ]    =-1.0*c1;             Cz[dV_ZERO_M1_P1     ]    = 1.0*c1;
Cx[dV_P1_ZERO_M1     ]    = 1.0*c1;             Cy[dV_P1_ZERO_M1     ]    = 0.0*c1;             Cz[dV_P1_ZERO_M1     ]    =-1.0*c1;
Cx[dV_M1_ZERO_P1     ]    =-1.0*c1;             Cy[dV_M1_ZERO_P1     ]    = 0.0*c1;             Cz[dV_M1_ZERO_P1     ]    = 1.0*c1;

Cx[dV_PH1_PH1_PH1    ]    = 0.5*c1;             Cy[dV_PH1_PH1_PH1    ]    = 0.5*c1;             Cz[dV_PH1_PH1_PH1    ]    = 0.5*c1;
Cx[dV_PH1_MH1_PH1    ]    = 0.5*c1;             Cy[dV_PH1_MH1_PH1    ]    =-0.5*c1;             Cz[dV_PH1_MH1_PH1    ]    = 0.5*c1;
Cx[dV_MH1_MH1_PH1    ]    =-0.5*c1;             Cy[dV_MH1_MH1_PH1    ]    =-0.5*c1;             Cz[dV_MH1_MH1_PH1    ]    = 0.5*c1;
Cx[dV_MH1_PH1_PH1    ]    =-0.5*c1;             Cy[dV_MH1_PH1_PH1    ]    = 0.5*c1;             Cz[dV_MH1_PH1_PH1    ]    = 0.5*c1;
Cx[dV_PH1_PH1_MH1    ]    = 0.5*c1;             Cy[dV_PH1_PH1_MH1    ]    = 0.5*c1;             Cz[dV_PH1_PH1_MH1    ]    =-0.5*c1;
Cx[dV_PH1_MH1_MH1    ]    = 0.5*c1;             Cy[dV_PH1_MH1_MH1    ]    =-0.5*c1;             Cz[dV_PH1_MH1_MH1    ]    =-0.5*c1;
Cx[dV_MH1_MH1_MH1    ]    =-0.5*c1;             Cy[dV_MH1_MH1_MH1    ]    =-0.5*c1;             Cz[dV_MH1_MH1_MH1    ]    =-0.5*c1;
Cx[dV_MH1_PH1_MH1    ]    =-0.5*c1;             Cy[dV_MH1_PH1_MH1    ]    = 0.5*c1;             Cz[dV_MH1_PH1_MH1    ]    =-0.5*c1;  
 
Cx[dV_P1_P1_P1       ]    = 1.0*c1;             Cy[dV_P1_P1_P1       ]    = 1.0*c1;             Cz[dV_P1_P1_P1       ]    = 1.0*c1;
Cx[dV_P1_M1_P1       ]    = 1.0*c1;             Cy[dV_P1_M1_P1       ]    =-1.0*c1;             Cz[dV_P1_M1_P1       ]    = 1.0*c1;
Cx[dV_M1_M1_P1       ]    =-1.0*c1;             Cy[dV_M1_M1_P1       ]    =-1.0*c1;             Cz[dV_M1_M1_P1       ]    = 1.0*c1;
Cx[dV_M1_P1_P1       ]    =-1.0*c1;             Cy[dV_M1_P1_P1       ]    = 1.0*c1;             Cz[dV_M1_P1_P1       ]    = 1.0*c1;
Cx[dV_P1_P1_M1       ]    = 1.0*c1;             Cy[dV_P1_P1_M1       ]    = 1.0*c1;             Cz[dV_P1_P1_M1       ]    =-1.0*c1;
Cx[dV_P1_M1_M1       ]    = 1.0*c1;             Cy[dV_P1_M1_M1       ]    =-1.0*c1;             Cz[dV_P1_M1_M1       ]    =-1.0*c1;
Cx[dV_M1_M1_M1       ]    =-1.0*c1;             Cy[dV_M1_M1_M1       ]    =-1.0*c1;             Cz[dV_M1_M1_M1       ]    =-1.0*c1;
Cx[dV_M1_P1_M1       ]    =-1.0*c1;             Cy[dV_M1_P1_M1       ]    = 1.0*c1;             Cz[dV_M1_P1_M1       ]    =-1.0*c1;    

    
    
//>ceil_functions
CxC[dV_ZERO_ZERO_ZERO ]    =ceil(    0.0);             CyC[dV_ZERO_ZERO_ZERO ]    =ceil(    0.0);             CzC[dV_ZERO_ZERO_ZERO ]    =ceil(    0.0);

CxC[dV_P1_ZERO_ZERO   ]    =ceil( 1.0*c1);             CyC[dV_P1_ZERO_ZERO   ]    =ceil(    0.0);             CzC[dV_P1_ZERO_ZERO   ]    =ceil(    0.0);
CxC[dV_M1_ZERO_ZERO   ]    =ceil(-1.0*c1);             CyC[dV_M1_ZERO_ZERO   ]    =ceil(    0.0);             CzC[dV_M1_ZERO_ZERO   ]    =ceil(    0.0);
CxC[dV_ZERO_P1_ZERO   ]    =ceil(    0.0);             CyC[dV_ZERO_P1_ZERO   ]    =ceil( 1.0*c1);             CzC[dV_ZERO_P1_ZERO   ]    =ceil(    0.0);
CxC[dV_ZERO_M1_ZERO   ]    =ceil(    0.0);             CyC[dV_ZERO_M1_ZERO   ]    =ceil(-1.0*c1);             CzC[dV_ZERO_M1_ZERO   ]    =ceil(    0.0);
CxC[dV_ZERO_ZERO_P1   ]    =ceil(    0.0);             CyC[dV_ZERO_ZERO_P1   ]    =ceil(    0.0);             CzC[dV_ZERO_ZERO_P1   ]    =ceil( 1.0*c1);
CxC[dV_ZERO_ZERO_M1   ]    =ceil(    0.0);             CyC[dV_ZERO_ZERO_M1   ]    =ceil(    0.0);             CzC[dV_ZERO_ZERO_M1   ]    =ceil(-1.0*c1);

CxC[dV_P2_ZERO_ZERO   ]    =ceil( 2.0*c1);             CyC[dV_P2_ZERO_ZERO   ]    =ceil(    0.0);             CzC[dV_P2_ZERO_ZERO   ]    =ceil(    0.0);
CxC[dV_M2_ZERO_ZERO   ]    =ceil(-2.0*c1);             CyC[dV_M2_ZERO_ZERO   ]    =ceil(    0.0);             CzC[dV_M2_ZERO_ZERO   ]    =ceil(    0.0);
CxC[dV_ZERO_P2_ZERO   ]    =ceil(    0.0);             CyC[dV_ZERO_P2_ZERO   ]    =ceil( 2.0*c1);             CzC[dV_ZERO_P2_ZERO   ]    =ceil(    0.0);
CxC[dV_ZERO_M2_ZERO   ]    =ceil(    0.0);             CyC[dV_ZERO_M2_ZERO   ]    =ceil(-2.0*c1);             CzC[dV_ZERO_M2_ZERO   ]    =ceil(    0.0);
CxC[dV_ZERO_ZERO_P2   ]    =ceil(    0.0);             CyC[dV_ZERO_ZERO_P2   ]    =ceil(    0.0);             CzC[dV_ZERO_ZERO_P2   ]    =ceil( 2.0*c1);
CxC[dV_ZERO_ZERO_M2   ]    =ceil(    0.0);             CyC[dV_ZERO_ZERO_M2   ]    =ceil(    0.0);             CzC[dV_ZERO_ZERO_M2   ]    =ceil(-2.0*c1);

CxC[dV_P1_P1_ZERO     ]    =ceil( 1.0*c1);             CyC[dV_P1_P1_ZERO     ]    =ceil( 1.0*c1);             CzC[dV_P1_P1_ZERO     ]    =ceil( 0.0*c1);
CxC[dV_M1_M1_ZERO     ]    =ceil(-1.0*c1);             CyC[dV_M1_M1_ZERO     ]    =ceil(-1.0*c1);             CzC[dV_M1_M1_ZERO     ]    =ceil( 0.0*c1);
CxC[dV_ZERO_P1_P1     ]    =ceil( 0.0*c1);             CyC[dV_ZERO_P1_P1     ]    =ceil( 1.0*c1);             CzC[dV_ZERO_P1_P1     ]    =ceil( 1.0*c1); 
CxC[dV_ZERO_M1_M1     ]    =ceil( 0.0*c1);             CyC[dV_ZERO_M1_M1     ]    =ceil(-1.0*c1);             CzC[dV_ZERO_M1_M1     ]    =ceil(-1.0*c1);
CxC[dV_P1_ZERO_P1     ]    =ceil( 1.0*c1);             CyC[dV_P1_ZERO_P1     ]    =ceil( 0.0*c1);             CzC[dV_P1_ZERO_P1     ]    =ceil( 1.0*c1);
CxC[dV_M1_ZERO_M1     ]    =ceil(-1.0*c1);             CyC[dV_M1_ZERO_M1     ]    =ceil( 0.0*c1);             CzC[dV_M1_ZERO_M1     ]    =ceil(-1.0*c1);
CxC[dV_P1_M1_ZERO     ]    =ceil( 1.0*c1);             CyC[dV_P1_M1_ZERO     ]    =ceil(-1.0*c1);             CzC[dV_P1_M1_ZERO     ]    =ceil( 0.0*c1);
CxC[dV_M1_P1_ZERO     ]    =ceil(-1.0*c1);             CyC[dV_M1_P1_ZERO     ]    =ceil( 1.0*c1);             CzC[dV_M1_P1_ZERO     ]    =ceil( 0.0*c1);
CxC[dV_ZERO_P1_M1     ]    =ceil( 0.0*c1);             CyC[dV_ZERO_P1_M1     ]    =ceil( 1.0*c1);             CzC[dV_ZERO_P1_M1     ]    =ceil(-1.0*c1);
CxC[dV_ZERO_M1_P1     ]    =ceil( 0.0*c1);             CyC[dV_ZERO_M1_P1     ]    =ceil(-1.0*c1);             CzC[dV_ZERO_M1_P1     ]    =ceil( 1.0*c1);
CxC[dV_P1_ZERO_M1     ]    =ceil( 1.0*c1);             CyC[dV_P1_ZERO_M1     ]    =ceil( 0.0*c1);             CzC[dV_P1_ZERO_M1     ]    =ceil(-1.0*c1);
CxC[dV_M1_ZERO_P1     ]    =ceil(-1.0*c1);             CyC[dV_M1_ZERO_P1     ]    =ceil( 0.0*c1);             CzC[dV_M1_ZERO_P1     ]    =ceil( 1.0*c1);

CxC[dV_PH1_PH1_PH1    ]    =ceil( 0.5*c1);             CyC[dV_PH1_PH1_PH1    ]    =ceil( 0.5*c1);             CzC[dV_PH1_PH1_PH1    ]    =ceil( 0.5*c1);
CxC[dV_PH1_MH1_PH1    ]    =ceil( 0.5*c1);             CyC[dV_PH1_MH1_PH1    ]    =ceil(-0.5*c1);             CzC[dV_PH1_MH1_PH1    ]    =ceil( 0.5*c1);
CxC[dV_MH1_MH1_PH1    ]    =ceil(-0.5*c1);             CyC[dV_MH1_MH1_PH1    ]    =ceil(-0.5*c1);             CzC[dV_MH1_MH1_PH1    ]    =ceil( 0.5*c1);
CxC[dV_MH1_PH1_PH1    ]    =ceil(-0.5*c1);             CyC[dV_MH1_PH1_PH1    ]    =ceil( 0.5*c1);             CzC[dV_MH1_PH1_PH1    ]    =ceil( 0.5*c1);
CxC[dV_PH1_PH1_MH1    ]    =ceil( 0.5*c1);             CyC[dV_PH1_PH1_MH1    ]    =ceil( 0.5*c1);             CzC[dV_PH1_PH1_MH1    ]    =ceil(-0.5*c1);
CxC[dV_PH1_MH1_MH1    ]    =ceil( 0.5*c1);             CyC[dV_PH1_MH1_MH1    ]    =ceil(-0.5*c1);             CzC[dV_PH1_MH1_MH1    ]    =ceil(-0.5*c1);
CxC[dV_MH1_MH1_MH1    ]    =ceil(-0.5*c1);             CyC[dV_MH1_MH1_MH1    ]    =ceil(-0.5*c1);             CzC[dV_MH1_MH1_MH1    ]    =ceil(-0.5*c1);
CxC[dV_MH1_PH1_MH1    ]    =ceil(-0.5*c1);             CyC[dV_MH1_PH1_MH1    ]    =ceil( 0.5*c1);             CzC[dV_MH1_PH1_MH1    ]    =ceil(-0.5*c1); 

CxC[dV_P1_P1_P1       ]    =ceil( 1.0*c1);             CyC[dV_P1_P1_P1       ]    =ceil( 1.0*c1);             CzC[dV_P1_P1_P1       ]    =ceil( 1.0*c1);
CxC[dV_P1_M1_P1       ]    =ceil( 1.0*c1);             CyC[dV_P1_M1_P1       ]    =ceil(-1.0*c1);             CzC[dV_P1_M1_P1       ]    =ceil( 1.0*c1);
CxC[dV_M1_M1_P1       ]    =ceil(-1.0*c1);             CyC[dV_M1_M1_P1       ]    =ceil(-1.0*c1);             CzC[dV_M1_M1_P1       ]    =ceil( 1.0*c1);
CxC[dV_M1_P1_P1       ]    =ceil(-1.0*c1);             CyC[dV_M1_P1_P1       ]    =ceil( 1.0*c1);             CzC[dV_M1_P1_P1       ]    =ceil( 1.0*c1);
CxC[dV_P1_P1_M1       ]    =ceil( 1.0*c1);             CyC[dV_P1_P1_M1       ]    =ceil( 1.0*c1);             CzC[dV_P1_P1_M1       ]    =ceil(-1.0*c1);
CxC[dV_P1_M1_M1       ]    =ceil( 1.0*c1);             CyC[dV_P1_M1_M1       ]    =ceil(-1.0*c1);             CzC[dV_P1_M1_M1       ]    =ceil(-1.0*c1);
CxC[dV_M1_M1_M1       ]    =ceil(-1.0*c1);             CyC[dV_M1_M1_M1       ]    =ceil(-1.0*c1);             CzC[dV_M1_M1_M1       ]    =ceil(-1.0*c1);
CxC[dV_M1_P1_M1       ]    =ceil(-1.0*c1);             CyC[dV_M1_P1_M1       ]    =ceil( 1.0*c1);             CzC[dV_M1_P1_M1       ]    =ceil(-1.0*c1);    
   
    


    
//>floor_functions
CxF[dV_ZERO_ZERO_ZERO ]    =floor(    0.0);             CyF[dV_ZERO_ZERO_ZERO ]    =floor(    0.0);             CzF[dV_ZERO_ZERO_ZERO ]    =floor(    0.0);
CxF[dV_P1_ZERO_ZERO   ]    =floor( 1.0*c1);             CyF[dV_P1_ZERO_ZERO   ]    =floor(    0.0);             CzF[dV_P1_ZERO_ZERO   ]    =floor(    0.0);
CxF[dV_M1_ZERO_ZERO   ]    =floor(-1.0*c1);             CyF[dV_M1_ZERO_ZERO   ]    =floor(    0.0);             CzF[dV_M1_ZERO_ZERO   ]    =floor(    0.0);
CxF[dV_ZERO_P1_ZERO   ]    =floor(    0.0);             CyF[dV_ZERO_P1_ZERO   ]    =floor( 1.0*c1);             CzF[dV_ZERO_P1_ZERO   ]    =floor(    0.0);
CxF[dV_ZERO_M1_ZERO   ]    =floor(    0.0);             CyF[dV_ZERO_M1_ZERO   ]    =floor(-1.0*c1);             CzF[dV_ZERO_M1_ZERO   ]    =floor(    0.0);
CxF[dV_ZERO_ZERO_P1   ]    =floor(    0.0);             CyF[dV_ZERO_ZERO_P1   ]    =floor(    0.0);             CzF[dV_ZERO_ZERO_P1   ]    =floor( 1.0*c1);
CxF[dV_ZERO_ZERO_M1   ]    =floor(    0.0);             CyF[dV_ZERO_ZERO_M1   ]    =floor(    0.0);             CzF[dV_ZERO_ZERO_M1   ]    =floor(-1.0*c1);

CxF[dV_P2_ZERO_ZERO   ]    =floor( 2.0*c1);             CyF[dV_P2_ZERO_ZERO   ]    =floor(    0.0);             CzF[dV_P2_ZERO_ZERO   ]    =floor(    0.0);
CxF[dV_M2_ZERO_ZERO   ]    =floor(-2.0*c1);             CyF[dV_M2_ZERO_ZERO   ]    =floor(    0.0);             CzF[dV_M2_ZERO_ZERO   ]    =floor(    0.0);
CxF[dV_ZERO_P2_ZERO   ]    =floor(    0.0);             CyF[dV_ZERO_P2_ZERO   ]    =floor( 2.0*c1);             CzF[dV_ZERO_P2_ZERO   ]    =floor(    0.0);
CxF[dV_ZERO_M2_ZERO   ]    =floor(    0.0);             CyF[dV_ZERO_M2_ZERO   ]    =floor(-2.0*c1);             CzF[dV_ZERO_M2_ZERO   ]    =floor(    0.0);
CxF[dV_ZERO_ZERO_P2   ]    =floor(    0.0);             CyF[dV_ZERO_ZERO_P2   ]    =floor(    0.0);             CzF[dV_ZERO_ZERO_P2   ]    =floor( 2.0*c1);
CxF[dV_ZERO_ZERO_M2   ]    =floor(    0.0);             CyF[dV_ZERO_ZERO_M2   ]    =floor(    0.0);             CzF[dV_ZERO_ZERO_M2   ]    =floor(-2.0*c1);

CxF[dV_P1_P1_ZERO     ]    =floor( 1.0*c1);             CyF[dV_P1_P1_ZERO     ]    =floor( 1.0*c1);             CzF[dV_P1_P1_ZERO     ]    =floor( 0.0*c1);
CxF[dV_M1_M1_ZERO     ]    =floor(-1.0*c1);             CyF[dV_M1_M1_ZERO     ]    =floor(-1.0*c1);             CzF[dV_M1_M1_ZERO     ]    =floor( 0.0*c1);
CxF[dV_ZERO_P1_P1     ]    =floor( 0.0*c1);             CyF[dV_ZERO_P1_P1     ]    =floor( 1.0*c1);             CzF[dV_ZERO_P1_P1     ]    =floor( 1.0*c1); 
CxF[dV_ZERO_M1_M1     ]    =floor( 0.0*c1);             CyF[dV_ZERO_M1_M1     ]    =floor(-1.0*c1);             CzF[dV_ZERO_M1_M1     ]    =floor(-1.0*c1);
CxF[dV_P1_ZERO_P1     ]    =floor( 1.0*c1);             CyF[dV_P1_ZERO_P1     ]    =floor( 0.0*c1);             CzF[dV_P1_ZERO_P1     ]    =floor( 1.0*c1);
CxF[dV_M1_ZERO_M1     ]    =floor(-1.0*c1);             CyF[dV_M1_ZERO_M1     ]    =floor( 0.0*c1);             CzF[dV_M1_ZERO_M1     ]    =floor(-1.0*c1);
CxF[dV_P1_M1_ZERO     ]    =floor( 1.0*c1);             CyF[dV_P1_M1_ZERO     ]    =floor(-1.0*c1);             CzF[dV_P1_M1_ZERO     ]    =floor( 0.0*c1);
CxF[dV_M1_P1_ZERO     ]    =floor(-1.0*c1);             CyF[dV_M1_P1_ZERO     ]    =floor( 1.0*c1);             CzF[dV_M1_P1_ZERO     ]    =floor( 0.0*c1);
CxF[dV_ZERO_P1_M1     ]    =floor( 0.0*c1);             CyF[dV_ZERO_P1_M1     ]    =floor( 1.0*c1);             CzF[dV_ZERO_P1_M1     ]    =floor(-1.0*c1);
CxF[dV_ZERO_M1_P1     ]    =floor( 0.0*c1);             CyF[dV_ZERO_M1_P1     ]    =floor(-1.0*c1);             CzF[dV_ZERO_M1_P1     ]    =floor( 1.0*c1);
CxF[dV_P1_ZERO_M1     ]    =floor( 1.0*c1);             CyF[dV_P1_ZERO_M1     ]    =floor( 0.0*c1);             CzF[dV_P1_ZERO_M1     ]    =floor(-1.0*c1);
CxF[dV_M1_ZERO_P1     ]    =floor(-1.0*c1);             CyF[dV_M1_ZERO_P1     ]    =floor( 0.0*c1);             CzF[dV_M1_ZERO_P1     ]    =floor( 1.0*c1);

CxF[dV_PH1_PH1_PH1    ]    =floor( 0.5*c1);             CyF[dV_PH1_PH1_PH1    ]    =floor( 0.5*c1);             CzF[dV_PH1_PH1_PH1    ]    =floor( 0.5*c1);
CxF[dV_PH1_MH1_PH1    ]    =floor( 0.5*c1);             CyF[dV_PH1_MH1_PH1    ]    =floor(-0.5*c1);             CzF[dV_PH1_MH1_PH1    ]    =floor( 0.5*c1);
CxF[dV_MH1_MH1_PH1    ]    =floor(-0.5*c1);             CyF[dV_MH1_MH1_PH1    ]    =floor(-0.5*c1);             CzF[dV_MH1_MH1_PH1    ]    =floor( 0.5*c1);
CxF[dV_MH1_PH1_PH1    ]    =floor(-0.5*c1);             CyF[dV_MH1_PH1_PH1    ]    =floor( 0.5*c1);             CzF[dV_MH1_PH1_PH1    ]    =floor( 0.5*c1);
CxF[dV_PH1_PH1_MH1    ]    =floor( 0.5*c1);             CyF[dV_PH1_PH1_MH1    ]    =floor( 0.5*c1);             CzF[dV_PH1_PH1_MH1    ]    =floor(-0.5*c1);
CxF[dV_PH1_MH1_MH1    ]    =floor( 0.5*c1);             CyF[dV_PH1_MH1_MH1    ]    =floor(-0.5*c1);             CzF[dV_PH1_MH1_MH1    ]    =floor(-0.5*c1);
CxF[dV_MH1_MH1_MH1    ]    =floor(-0.5*c1);             CyF[dV_MH1_MH1_MH1    ]    =floor(-0.5*c1);             CzF[dV_MH1_MH1_MH1    ]    =floor(-0.5*c1);
CxF[dV_MH1_PH1_MH1    ]    =floor(-0.5*c1);             CyF[dV_MH1_PH1_MH1    ]    =floor( 0.5*c1);             CzF[dV_MH1_PH1_MH1    ]    =floor(-0.5*c1);   

CxF[dV_P1_P1_P1       ]    =floor( 1.0*c1);             CyF[dV_P1_P1_P1       ]    =floor( 1.0*c1);             CzF[dV_P1_P1_P1       ]    =floor( 1.0*c1);
CxF[dV_P1_M1_P1       ]    =floor( 1.0*c1);             CyF[dV_P1_M1_P1       ]    =floor(-1.0*c1);             CzF[dV_P1_M1_P1       ]    =floor( 1.0*c1);
CxF[dV_M1_M1_P1       ]    =floor(-1.0*c1);             CyF[dV_M1_M1_P1       ]    =floor(-1.0*c1);             CzF[dV_M1_M1_P1       ]    =floor( 1.0*c1);
CxF[dV_M1_P1_P1       ]    =floor(-1.0*c1);             CyF[dV_M1_P1_P1       ]    =floor( 1.0*c1);             CzF[dV_M1_P1_P1       ]    =floor( 1.0*c1);
CxF[dV_P1_P1_M1       ]    =floor( 1.0*c1);             CyF[dV_P1_P1_M1       ]    =floor( 1.0*c1);             CzF[dV_P1_P1_M1       ]    =floor(-1.0*c1);
CxF[dV_P1_M1_M1       ]    =floor( 1.0*c1);             CyF[dV_P1_M1_M1       ]    =floor(-1.0*c1);             CzF[dV_P1_M1_M1       ]    =floor(-1.0*c1);
CxF[dV_M1_M1_M1       ]    =floor(-1.0*c1);             CyF[dV_M1_M1_M1       ]    =floor(-1.0*c1);             CzF[dV_M1_M1_M1       ]    =floor(-1.0*c1);
CxF[dV_M1_P1_M1       ]    =floor(-1.0*c1);             CyF[dV_M1_P1_M1       ]    =floor( 1.0*c1);             CzF[dV_M1_P1_M1       ]    =floor(-1.0*c1);    
    
    
    
}





