#ifndef _lbmD2Q9_H_
#define _lbmD2Q9_H_

enum velocityDir{ dV_ZERO_ZERO, dV_P1_ZERO, dV_ZERO_P1, dV_M1_ZERO, dV_ZERO_M1, dV_P1_P1, dV_M1_P1, dV_M1_M1, dV_P1_M1};

int oppdV[9] =  {dV_ZERO_ZERO,  dV_M1_ZERO, dV_ZERO_M1, dV_P1_ZERO, dV_ZERO_P1, dV_M1_M1, dV_P1_M1, dV_P1_P1, dV_M1_P1};


namespace MyEnum //////wrtoes 
{
  enum Type
  {
    dv_L0 = dV_ZERO_ZERO,   dv_L1 = dV_P1_ZERO, dv_L2 = dV_P1_P1,dv_L3 = dV_ZERO_P1,dv_L4 = dV_M1_P1,
                            dv_L5 = dV_M1_ZERO, dv_L6 = dV_M1_M1,dv_L7 = dV_ZERO_M1,dv_L8 = dV_ZERO_M1,
    
  };

  static const Type All[] = { dv_L0, dv_L1, dv_L2,dv_L3,dv_L4,dv_L5,dv_L6,dv_L7,dv_L8 };
}


// const auto e = MyEnum::All;
// for(int dv = 2; dv<6; dv++){
//     int goinout;
//     goinout = (int)*(e+dv);
//     std::cout<<goinout<<std::endl;
// }



template<typename T >
struct lbmD2Q9
{
    int dvN;

     T W[9]; 
     T Cx[9]; 
     T Cy[9]; 
     T Cz[9]; 

    T theta0, thetaInverse ;
    lbmD2Q9(T c,T cs2); //constructor
    
};





//constructor
template<typename T>
lbmD2Q9<T>::lbmD2Q9(T c1,T cs2)
    {
        dvN = 9;
        theta0 = 1.0/3.0;
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

