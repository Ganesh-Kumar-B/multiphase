
/*********************************************************************************************
 * Copyright (c) <2024>, <Santosh Ansumali@JNCASR>                                                                                                             *
 *  All rights reserved.                                                                                                                                                                          *
 *   Redistribution and use in source and binary forms, with or without modification, are                                                                *
 *   permitted provided that the following conditions are met:                                                                                                           *
 *                                                                                            *
 *    1. Redistributions of source code must retain the above copyright notice, this list of                                                               *
 *       conditions and the following disclaimer.                                                                                                                              *
 *    2. Redistributions in binary form must reproduce the above copyright notice, this list                                                              *
 *       of conditions and the following disclaimer in the documentation and/or other                                                                   *
 *       materials provided with the distribution.                                                                                                                               *
 *    3. Neither the name of the <JNCASR> nor the names of its contributors may be used to                                                         *
 *       endorse or promote products derived from this software without specific prior                                                                  *
 *       written permission.                                                                                                                                                                *
 *                                                                                             *
 *       THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND                           *
 *       ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED                                      *
 *       WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.                     *
 *       IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,                           *
 *       INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,                                  *
 *       BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,                                *
 *       DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF                           *
 *       LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE                                      *
 *       OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED                              *
 *       OF THE POSSIBILITY OF SUCH DAMAGE.                                                                                                                          *
 *                                                                                             *
 *       Suggestions:          ansumali@jncasr.ac.in                                                                                                                  *
 *       Bugs:                 ansumali@jncasr.ac.in                                                                                                                       *
 *       Contributors:                                                                                                                                                                           *
 *       Last Change Date:                                                                                                                                                                  *
 *                                                                                             *
 ************************************************************************************************************************************************/
 


// Base Version for understanding
// no of operation is 9* (3*2+1)  +19 = 82  ops
template<typename dataType1>
void getHydroMoment(lbmD2Q9<dataType1> &lbModel, dataType1 *rho, dataType1 *uX, dataType1 *uY, dataType1 *theta)
{

    dataType1 sum(0.0);

    rho = theta = 0.0;
    uX = uY  = 0.0;
    for(int dv=0; dv< lbModel.dvN; dv++)  {
        rho   += lbModel.fTemp[dv];
        uX    += lbModel.fTemp[dv] * lbModel.cX[ dv];
        uY    += lbModel.fTemp[dv] * lbModel.cY[ dv];
        theta += lbModel.fTemp[dv] * lbModel.cSq[dv];
    }
     
    // 9 ops + 1 div = 19 ops
    // const ops does not contribute as it is evaluated once only in pgm lifetime

    const dataType1  oneByTwo = 1.0/2.0;
    dataType1  oneByRho = 1.0/rho;
    uX    *= oneByRho;
    uY    *= oneByRho;
    theta  = oneByTwo * (   theta - rho * (uX * uX + uY * uY )   );
    theta *= oneByRho;
    
   
}

// Better Version for understanding
// no of operation is = 20 + 19 =39  ops
template<typename dataType1>
void getHydroMomentV1(lbmD2Q9<dataType1> &lbModel, dataType1 &rho, dataType1 &uX, dataType1 &uY,   dataType1 &theta)
{
    dataType1 sum(0.0);

    rho = theta = 0.0;
    
    rho = lbModel.fTemp[lbModel.DV_ZERO_ZERO];
    
    dataType1 t1(0.0), t2(0.0), t3(0.0);
    
    // 8  ops in SC
    t1 = lbModel.fTemp[lbModel.DV_P_ZERO] + lbModel.fTemp[lbModel.DV_M_ZERO];
    uX = lbModel.fTemp[lbModel.DV_P_ZERO] - lbModel.fTemp[lbModel.DV_M_ZERO];
     
     
    t1 += lbModel.fTemp[lbModel.DV_ZERO_P] + lbModel.fTemp[lbModel.DV_ZERO_M];
    uY  = lbModel.fTemp[lbModel.DV_ZERO_P] - lbModel.fTemp[lbModel.DV_ZERO_M];
     
   
    theta += 2.0 * t1;
    rho   +=       t1;
    
    // 12 ops in FCC
    
    t1  = lbModel.fTemp[lbModel.DV_P_P] + lbModel.fTemp[lbModel.DV_M_M] ;
    t2  = lbModel.fTemp[lbModel.DV_P_P] - lbModel.fTemp[lbModel.DV_M_M] ;
    t1 += lbModel.fTemp[lbModel.DV_P_M] + lbModel.fTemp[lbModel.DV_M_P] ;
    t3  = lbModel.fTemp[lbModel.DV_P_M] - lbModel.fTemp[lbModel.DV_M_P] ;
    
    theta +=       2.0 * t1 ;
    rho   +=              t1;
    uX    +=       t2 +  t3 ;
    uY    +=       t2 -  t3 ;
    

    // 9 ops + 1 div = 19 ops
    // const ops does not contribute as it is evaluated once only in pgm lifetime
    const dataType1  oneByTwo = 1.0/2.0;
    dataType1  oneByRho = 1.0/rho;
    uX    *= oneByRho;
    uY    *= oneByRho;
    theta  = oneByTwo * (   theta - rho * (uX * uX + uY * uY )   );
    theta  *= oneByRho;
}
   

 


template<typename dataType1>
void getFEqIsoSecond(lbmD2Q9<dataType1> &lbModel,  dataType1 &rho, dataType1 &uX, dataType1 &uY)
{
    
    dataType1 t2 = 1.0 - 1.5*(uX*uX + uY*uY) ;
    dataType1 dot1 =3.0 * uX;
    dataType1 dot2 =3.0 * uY;
    
    for(int dv=0;dv<lbModel.dvN;dv++)
      lbModel.fEq[dv]   = rho * lbModel.wt[dv];
    
    lbModel.fEq[lbModel.DV_ZERO_ZERO] *= t2;
    
    dataType1 t1 = 0.5*dot1*dot1;
    lbModel.fEq[lbModel.DV_P_ZERO] *= (t2 + dot1 + t1);
    lbModel.fEq[lbModel.DV_M_ZERO] *= (t2 - dot1 + t1);
    t1 = 0.5*dot2*dot2;
    lbModel.fEq[lbModel.DV_ZERO_P] *= (t2 + dot2 + t1);
    lbModel.fEq[lbModel.DV_ZERO_M] *= (t2 - dot2 + t1);
    
    t1 = 0.5 * (dot1 + dot2) * (dot1 + dot2);
    lbModel.fEq[lbModel.DV_P_P] *=  (t2 + (dot1 + dot2)+ t1);
    lbModel.fEq[lbModel.DV_M_M] *=  (t2 - (dot1 + dot2)+ t1);

    t1 = 0.5 * (dot1 - dot2) * (dot1 - dot2);
    lbModel.fEq[lbModel.DV_M_P] *=  (t2 - (dot1 - dot2) + t1);
    lbModel.fEq[lbModel.DV_P_M] *=  (t2 + (dot1 - dot2) + t1);
    
     
}



template<typename dataType>
void advectionD2Q9(field2D<dataType,9> &lbmGrid,  lbmD2Q9<dataType>  &d2q9Model ){

    for(int i2 = lbmGrid.n2Begin ; i2 <= lbmGrid.n2End; i2++){
        for(int i1 = lbmGrid.n1Begin ; i1 <= lbmGrid.n1End; i1++){
            lbmGrid(i1, i2, d2q9Model.DV_ZERO_M ) = lbmGrid(i1  , i2+1, d2q9Model.DV_ZERO_M );
            lbmGrid(i1, i2, d2q9Model.DV_M_ZERO ) = lbmGrid(i1  , i2+1, d2q9Model.DV_M_ZERO );
            lbmGrid(i1, i2, d2q9Model.DV_M_M    ) = lbmGrid(i1+1, i2+1, d2q9Model.DV_M_M    );
            lbmGrid(i1, i2, d2q9Model.DV_P_M    ) = lbmGrid(i1-1, i2+1, d2q9Model.DV_P_M    );
            
        }
    }
    
    for(int i2 = lbmGrid.n2End ; i2 >= lbmGrid.n2Begin; i2--){
        for(int i1 = lbmGrid.n1End ; i1 >= lbmGrid.n1Begin; i1--){
            lbmGrid(i1, i2, d2q9Model.DV_ZERO_P ) = lbmGrid(i1  , i2-1, d2q9Model.DV_ZERO_P );
            lbmGrid(i1, i2, d2q9Model.DV_P_ZERO ) = lbmGrid(i1  , i2-1, d2q9Model.DV_P_ZERO );
            lbmGrid(i1, i2, d2q9Model.DV_M_P    ) = lbmGrid(i1+1, i2-1, d2q9Model.DV_M_P );
            lbmGrid(i1, i2, d2q9Model.DV_P_P    ) = lbmGrid(i1-1, i2-1, d2q9Model.DV_P_P );
        }
    }
       
}

template<typename dataType>
void getDenField(field2D<dataType,9> &lbmGrid,  lbmD2Q9<dataType>  &d2q9Model,dataType dt, field2D<dataType,2> denField){
    dataType rho;
    for(int i2 = lbmGrid.n2Begin ; i2 <= lbmGrid.n2End; i2++){
        for(int i1 = lbmGrid.n1Begin ; i1 <= lbmGrid.n1End; i1++){
            rho=0.0;
            for(int dv =0; dv <d2q9Model.dvN; dv++ )
                rho += lbmGrid(i1, i2,dv);
            denField(i1, i2,0) = rho;
        }
    }

    denField.makePeriodicX();
    denField.makePeriodicY();
}



template<typename dataType>
void collideD2Q9(field2D<dataType,9> &lbmGrid,  lbmD2Q9<dataType>  &d2q9Model, dataType beta,dataType dt, field2D<dataType,2> forceField){
   
    dataType rho, uX, uY,theta;
    for(int i2 = lbmGrid.n2Begin ; i2 <= lbmGrid.n2End; i2++){
        for(int i1 = lbmGrid.n1Begin ; i1 <= lbmGrid.n1End; i1++){
            for(int dv =0; dv <d2q9Model.dvN; dv++ )
                d2q9Model.fTemp[dv] = lbmGrid(i1, i2,dv);

            getHydroMomentV1(d2q9Model, rho, uX,  uY,   theta);
            uX += forceField(i1,i2,0) * dt * 0.5;
            uY += forceField(i1,i2,1) * dt * 0.5;
            getFEqIsoSecond(d2q9Model, rho, uX,  uY);
        
            for (int dv = 0; dv< d2q9Model.dvN; dv++){
                dataType dot = forceField(i1,i2,0) * (d2q9Model.cX[dv] - uX) + forceField(i1,i2,1) * (d2q9Model.cY[dv] - uY);
                lbmGrid(i1,i2,dv) +=    2.0* beta*(d2q9Model.fEq[dv] - lbmGrid(i1,i2,dv))+ (1.0 - beta)*dt*d2q9Model.oneByTheta0 *  d2q9Model.fEq[dv] * (dot );
            }
        }
    }
}

template<typename dataType>
void initializeTaylorGreen(field2D<dataType,9> &lbmGrid,  lbmD2Q9<dataType>  &d2q9Model, dataType u_ref){

    dataType rho = 1.0;

    for(int i2 = lbmGrid.n2Begin ; i2 <= lbmGrid.n2End; i2++){
        for(int i1 = lbmGrid.n1Begin ; i1 <= lbmGrid.n1End; i1++){
            
            dataType x = (i1/lbmGrid.n1)*M_PI,   y = (i2/lbmGrid.n2)*M_PI;
            dataType uX = u_ref*sin(x)*cos(y) , uY = u_ref* -cos(x)*sin(y);


            getFEqIsoSecond(d2q9Model, rho,uX,uY);

            for(int dv =0; dv <d2q9Model.dvN; dv++ )
                lbmGrid(i1,i2,dv) = d2q9Model.fEq[dv];

        }
    }

}

template<typename dataType>
void initializeFEq(field2D<dataType,9> &lbmGrid,  lbmD2Q9<dataType>  &d2q9Model){

    dataType rho = 1.0 , uX = 0, uY = 0;

    for(int i2 = lbmGrid.n2Begin ; i2 <= lbmGrid.n2End; i2++){
        for(int i1 = lbmGrid.n1Begin ; i1 <= lbmGrid.n1End; i1++){
          
            getFEqIsoSecond(d2q9Model, rho,uX,uY);

            std::cout<<"rho:    ="<<d2q9Model.fEq[5]<<std::endl;
            for(int dv =0; dv <d2q9Model.dvN; dv++ )
                lbmGrid(i1,i2,dv) = d2q9Model.fEq[dv];

        }
    }

}


template<typename dataType>
void calculateForce(field2D<dataType,9> &lbmGrid,  lbmD2Q9<dataType>  &d2q9Model, field2D<dataType,2> &forceField, dataType g){

    for(int i2 = lbmGrid.n2Begin ; i2 <= lbmGrid.n2End; i2++){
        for(int i1 = lbmGrid.n1Begin ; i1 <= lbmGrid.n1End; i1++){

            dataType x = i1/lbmGrid.n1, y = i2/lbmGrid.n2;

            if(y > 0.5){
                forceField(i1,i2,0) = g;
            }else{
                forceField(i1,i2,0) = -1.0* g ;
            }

            

        }
    }

}




template<typename dataType>
void printVtk(field2D<dataType,9> &lbmGrid,  lbmD2Q9<dataType>  &d2q9Model, field2D<dataType,2> denField, int timeStep){


    std::ofstream file;
    char fileName[250];

    dataType ux, uy ;
    sprintf(fileName,"./Result/velocity_%d.vtk", timeStep) ;
    file.open(fileName);

    file<<"# vtk DataFile Version 3.0\nVelocity\nASCII\nDATASET STRUCTURED_POINTS"<<std::endl;
    
    file<<"DIMENSIONS "<<lbmGrid.n1<<" "<<lbmGrid.n2<<" "<<1<<std::endl;
    
    file<<"ORIGIN "<<0<<" "<<0<<" "<<0<<std::endl;
    file<<"SPACING "<<1<<" "<<1<<" "<<1<<std::endl;


    file<<"POINT_DATA "<<1*lbmGrid.n1*1*lbmGrid.n2<<std::endl;
    file<<"SCALARS density double 1\nLOOKUP_TABLE default"<<std::endl;

    dataType rho, uX, uY,theta;

    for(int i2 = lbmGrid.n2Begin ; i2 <= lbmGrid.n2End; i2++){
        for(int i1 = lbmGrid.n1Begin ; i1 <= lbmGrid.n1End; i1++){
        
            getHydroMomentV1(d2q9Model, rho, uX,  uY,   theta);
        
            file<<rho<<std::endl;

        }
    }

    file<<"VECTORS velocity double"<<std::endl;

    for(int i2 = lbmGrid.n2Begin ; i2 <= lbmGrid.n2End; i2++){
        for(int i1 = lbmGrid.n1Begin ; i1 <= lbmGrid.n1End; i1++){
        
            getHydroMomentV1(d2q9Model, rho, uX,  uY,   theta);
        
            file<<uX<<" "<<uY<<" "<<0.0<<std::endl;;

        }
    }


}











