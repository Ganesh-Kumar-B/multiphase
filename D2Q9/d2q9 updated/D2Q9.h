/*********************************************************************************************
 * Copyright (c) 2024, Santosh Ansumali@JNCASR
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without modification, 
 * are permitted provided that the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice, 
 *    this list of conditions and the following disclaimer in the documentation 
 *    and/or other materials provided with the distribution.
 * 3. Neither the name of JNCASR nor the names of its contributors may be used to 
 *    endorse or promote products derived from this software without specific 
 *    prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE 
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
 * POSSIBILITY OF SUCH DAMAGE.
 * 
 * Suggestions: ansumali@jncasr.ac.in
 * Bugs: ansumali@jncasr.ac.in
 * Contributors:
 * Last Change Date:
 *********************************************************************************************/

#ifndef _MESOSCALE_D2Q9_H_
#define _MESOSCALE_D2Q9_H_
#include<math.h>
#include<string.h>
#include<fstream>
 
template<typename dataType>
struct lbmD2Q9
{
    lbmD2Q9(dataType _speed)
    {
        latticeSpeed  = _speed;
    }
     void setModelParameters(){
         
        theta0          =  1.0/3.0;
        oneByTheta0     =  3.0;
        oneByTheta0Sqr  =  9.0;
        cs = sqrt((5.0/3.0)*theta0);
         
          /******************************************************************************
          *    Constraints  on weights                                                                                                                               *
          *   w0+  4 * wSC + 4 * wBCC = 1                                                                                                                    *
          *       2 * wSC + 4 * wBCC = theta0                                                                                                                *
          *       2 * wSC + 4 * wBCC = 3* theta0^2                                                                                                        *
          *              + 4 * wBCC = theta0^2                                                                                                                     *
          **********************************************************************************************************************/
         
        w0              =  16.0 / 36.0;
        wSC             =  4.0  / 36.0;
        wBCC            =  1.0  / 36.0;
          
        wt[DV_ZERO_ZERO ] = w0;    cX[DV_ZERO_ZERO ]  =   0.0;     cY[DV_ZERO_ZERO ]   =   0.0;
          
        wt[DV_P_ZERO ]    = wSC;   cX[DV_P_ZERO    ]  =   1.0;     cY[DV_P_ZERO    ]   =   0.0;
        wt[DV_M_ZERO ]    = wSC;   cX[DV_M_ZERO    ]  =  -1.0;     cY[DV_M_ZERO    ]   =   0.0;
        wt[DV_ZERO_P ]    = wSC;   cX[DV_ZERO_P    ]  =   0.0;     cY[DV_ZERO_P    ]   =   1.0;
        wt[DV_ZERO_M ]    = wSC;   cX[DV_ZERO_M    ]  =   0.0;     cY[DV_ZERO_M    ]   =  -1.0;
          
        wt[DV_P_P    ]    = wBCC;  cX[DV_P_P       ]  =   1.0;     cY[DV_P_P       ]   =   1.0;
        wt[DV_M_P    ]    = wBCC;  cX[DV_M_P       ]  =  -1.0;     cY[DV_M_P       ]   =   1.0;
        wt[DV_M_M    ]    = wBCC;  cX[DV_M_M       ]  =  -1.0;     cY[DV_M_M       ]   =  -1.0;
        wt[DV_P_M    ]    = wBCC;  cX[DV_P_M       ]  =   1.0;     cY[DV_P_M       ]   =  -1.0;

        
          
           for(int dv=0; dv<  dvN; dv++)  {
                     cX2[dv] =  cX[dv] *  cX[dv]  ;
                     cY2[dv] =  cY[dv] *  cY[dv]  ;
                     cSq[dv] =  cX2[dv] + cY2[dv] ;
      }
             
    }
          
    const static int dvN              =   9;
          
    const static int  DV_ZERO_ZERO    =   0;
    const static int  DV_P_ZERO       =   1;
    const static int  DV_M_ZERO       =   2;
    const static int  DV_ZERO_P       =   3;
    const static int  DV_ZERO_M       =   4;
    const static int  DV_P_P          =   5;
    const static int  DV_M_P          =   6;
    const static int  DV_M_M          =   7;
    const static int  DV_P_M          =   8;

    //field indices
    const static int RHO_INDEX = 0;
    const static int UX_INDEX = 1;
    const static int UY_INDEX = 2;
    const static int THETA_INDEX = 3;
    
    
    // velocity models
    dataType wt[dvN];
    dataType cX[dvN];
    dataType cY[dvN];
    dataType cX2[dvN];
    dataType cY2[dvN];
    dataType cSq[dvN];
    
    //Lattice Parameters
    dataType latticeSpeed;
    dataType theta0; //Reference temperature
    dataType oneByTheta0   ;
    dataType oneByTheta0Sqr;
    dataType cs;
    dataType w0;
    dataType wSC;
    dataType wBCC;
    dataType fTemp[dvN];
    dataType fEq[dvN];
     
   
};

#include"D2Q9.C"
#endif
