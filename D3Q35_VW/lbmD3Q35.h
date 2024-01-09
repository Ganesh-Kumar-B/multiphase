#ifndef _lbmD3Q35_H_
#define _lbmD3Q35_H_

enum velocityDir{dV_ZERO_ZERO_ZERO, 
                dV_P2_ZERO_ZERO, dV_M2_ZERO_ZERO, dV_ZERO_P2_ZERO, dV_ZERO_M2_ZERO, dV_ZERO_ZERO_P2,dV_ZERO_ZERO_M2,

                dV_P1_P1_ZERO, dV_M1_P1_ZERO ,dV_M1_M1_ZERO,dV_P1_M1_ZERO , 
                dV_P1_ZERO_P1,dV_M1_ZERO_P1,dV_M1_ZERO_M1,dV_P1_ZERO_M1,
                dV_ZERO_P1_P1,dV_ZERO_M1_P1,dV_ZERO_M1_M1,dV_ZERO_P1_M1,
                
    

                dV_PH1_PH1_PH1, dV_PH1_MH1_PH1, dV_MH1_MH1_PH1, dV_MH1_PH1_PH1, dV_PH1_PH1_MH1, dV_PH1_MH1_MH1, dV_MH1_MH1_MH1, dV_MH1_PH1_MH1,
                dV_P2_P2_P2, dV_P2_M2_P2, dV_M2_M2_P2, dV_M2_P2_P2, dV_P2_P2_M2, dV_P2_M2_M2, dV_M2_M2_M2, dV_M2_P2_M2,
                };

int oppdV[35] =  {  dV_ZERO_ZERO_ZERO,
                    dV_M2_ZERO_ZERO, dV_P2_ZERO_ZERO, dV_ZERO_M2_ZERO, dV_ZERO_P2_ZERO, dV_ZERO_ZERO_M2 ,dV_ZERO_ZERO_P2,
                    
                    dV_M1_M1_ZERO,  dV_P1_P1_ZERO,  dV_ZERO_M1_M1,dV_ZERO_P1_P1 , dV_M1_ZERO_M1,dV_P1_ZERO_P1,
                    dV_M1_P1_ZERO, dV_P1_M1_ZERO ,dV_ZERO_M1_P1,dV_ZERO_P1_M1 , dV_M1_ZERO_P1,dV_P1_ZERO_M1,

                    };




// const auto e = MyEnum::All;
// for(int dv = 2; dv<6; dv++){
//     int goinout;
//     goinout = (int)*(e+dv);
//     std::cout<<goinout<<std::endl;
// }



template<typename T >
struct lbmD3Q35
{
    int dvN;

     T W[35]; 
     T Cx[35]; 
     T Cy[35]; 
     T Cz[35]; 

    T theta0, thetaInverse ;
    lbmD3Q35(T c,T cs2); //constructor
    
};





//constructor
template<typename T>
lbmD3Q35<T>::lbmD3Q35(T c1,T cs2)
    {
        dvN = 35;
        theta0 = 1.0/3.0;
        thetaInverse = 1.0/theta0;


W[dV_ZERO_ZERO_ZERO ]    =1.0/3.0;

W[dV_P1_ZERO_ZERO   ]    =1.0/18.0;	            	
W[dV_M1_ZERO_ZERO   ]    =1.0/18.0;	            	
W[dV_ZERO_P1_ZERO   ]    =1.0/18.0;	            	
W[dV_ZERO_M1_ZERO   ]    =1.0/18.0;	            	
W[dV_ZERO_ZERO_P1   ]    =1.0/18.0;		        	
W[dV_ZERO_ZERO_M1   ]    =1.0/18.0;

W[dV_P1_P1_ZERO     ]    =1.0/36.0;		        	
W[dV_M1_M1_ZERO     ]    =1.0/36.0;		        	
W[dV_ZERO_P1_P1     ]    =1.0/36.0;	            
W[dV_ZERO_M1_M1     ]    =1.0/36.0;	            
W[dV_P1_ZERO_P1     ]    =1.0/36.0;	            
W[dV_M1_ZERO_M1     ]    =1.0/36.0;	            
W[dV_P1_M1_ZERO     ]    =1.0/36.0;	            
W[dV_M1_P1_ZERO     ]    =1.0/36.0;	            
W[dV_ZERO_P1_M1     ]    =1.0/36.0;	            
W[dV_ZERO_M1_P1     ]    =1.0/36.0;	        	
W[dV_P1_ZERO_M1     ]    =1.0/36.0;	        	
W[dV_M1_ZERO_P1     ]    =1.0/36.0;	        	
    
    
Cx[dV_ZERO_ZERO_ZERO ]    =0.0;         Cy[dV_ZERO_ZERO_ZERO ]    =0.0;             Cz[dV_ZERO_ZERO_ZERO ]    =0.0;

Cx[dV_P1_ZERO_ZERO   ]    = c1;         Cy[dV_P1_ZERO_ZERO   ]    =0.0;             Cz[dV_P1_ZERO_ZERO   ]    =0.0;
Cx[dV_M1_ZERO_ZERO   ]    =-c1;         Cy[dV_M1_ZERO_ZERO   ]    =0.0;             Cz[dV_M1_ZERO_ZERO   ]    =0.0;
Cx[dV_ZERO_P1_ZERO   ]    =0.0;         Cy[dV_ZERO_P1_ZERO   ]    = c1;             Cz[dV_ZERO_P1_ZERO   ]    =0.0;
Cx[dV_ZERO_M1_ZERO   ]    =0.0;         Cy[dV_ZERO_M1_ZERO   ]    =-c1;             Cz[dV_ZERO_M1_ZERO   ]    =0.0;
Cx[dV_ZERO_ZERO_P1   ]    =0.0;         Cy[dV_ZERO_ZERO_P1   ]    =0.0;             Cz[dV_ZERO_ZERO_P1   ]    = c1;
Cx[dV_ZERO_ZERO_M1   ]    =0.0;         Cy[dV_ZERO_ZERO_M1   ]    =0.0;             Cz[dV_ZERO_ZERO_M1   ]    =-c1;

Cx[dV_P1_P1_ZERO     ]    = c1;         Cy[dV_P1_P1_ZERO     ]    = c1;             Cz[dV_P1_P1_ZERO     ]    =0.0;
Cx[dV_M1_M1_ZERO     ]    =-c1;         Cy[dV_M1_M1_ZERO     ]    =-c1;             Cz[dV_M1_M1_ZERO     ]    =0.0;
Cx[dV_ZERO_P1_P1     ]    =0.0;         Cy[dV_ZERO_P1_P1     ]    = c1;             Cz[dV_ZERO_P1_P1     ]    = c1; 
Cx[dV_ZERO_M1_M1     ]    =0.0;         Cy[dV_ZERO_M1_M1     ]    =-c1;             Cz[dV_ZERO_M1_M1     ]    =-c1;
Cx[dV_P1_ZERO_P1     ]    = c1;         Cy[dV_P1_ZERO_P1     ]    =0.0;             Cz[dV_P1_ZERO_P1     ]    = c1;
Cx[dV_M1_ZERO_M1     ]    =-c1;         Cy[dV_M1_ZERO_M1     ]    =0.0;             Cz[dV_M1_ZERO_M1     ]    =-c1;
Cx[dV_P1_M1_ZERO     ]    = c1;         Cy[dV_P1_M1_ZERO     ]    =-c1;             Cz[dV_P1_M1_ZERO     ]    =0.0;
Cx[dV_M1_P1_ZERO     ]    =-c1;         Cy[dV_M1_P1_ZERO     ]    = c1;             Cz[dV_M1_P1_ZERO     ]    =0.0;
Cx[dV_ZERO_P1_M1     ]    =0.0;         Cy[dV_ZERO_P1_M1     ]    = c1;             Cz[dV_ZERO_P1_M1     ]    =-c1;
Cx[dV_ZERO_M1_P1     ]    =0.0;         Cy[dV_ZERO_M1_P1     ]    =-c1;             Cz[dV_ZERO_M1_P1     ]    = c1;
Cx[dV_P1_ZERO_M1     ]    = c1;         Cy[dV_P1_ZERO_M1     ]    =0.0;             Cz[dV_P1_ZERO_M1     ]    =-c1;
Cx[dV_M1_ZERO_P1     ]    =-c1;         Cy[dV_M1_ZERO_P1     ]    =0.0;             Cz[dV_M1_ZERO_P1     ]    = c1;
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    }



#endif

