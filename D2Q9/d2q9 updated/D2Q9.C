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

#pragma once
#include<iostream>
#include<math.h>
#include "D2Q9.h"


// Base Version for understanding
// no of operation is 9* (3*2+1)  +19 = 82  ops
template <typename dataType1>
void getHydroMoment(lbmD2Q9<dataType1> &lbModel,
                    dataType1 *rho,
                    dataType1 *uX,
                    dataType1 *uY,
                    dataType1 *theta)
{

    dataType1 sum(0.0);

    rho = theta = 0.0;
    uX = uY = 0.0;
    for (int dv = 0; dv < lbModel.dvN; dv++)
    {
        rho += lbModel.fTemp[dv];
        uX += lbModel.fTemp[dv] * lbModel.cX[dv];
        uY += lbModel.fTemp[dv] * lbModel.cY[dv];
        theta += lbModel.fTemp[dv] * lbModel.cSq[dv];
    }

    // 9 ops + 1 div = 19 ops
    // const ops does not contribute as it is evaluated once only in pgm lifetime

    const dataType1 oneByTwo = 1.0 / 2.0;
    dataType1 oneByRho = 1.0 / rho;
    uX *= oneByRho;
    uY *= oneByRho;
    theta = oneByTwo * (theta - rho * (uX * uX + uY * uY));
    theta *= oneByRho;
}

// Better Version for understanding
// no of operation is = 20 + 19 =39  ops
template <typename dataType1>
void getHydroMomentV1(lbmD2Q9<dataType1> &lbModel, dataType1 &rho,
                      dataType1 &uX, dataType1 &uY, dataType1 &theta)
{
    dataType1 sum(0.0);

    rho = theta = 0.0;

    rho = lbModel.fTemp[lbModel.DV_ZERO_ZERO];

    dataType1 t1(0.0), t2(0.0), t3(0.0);

    // 8  ops in SC
    t1 = lbModel.fTemp[lbModel.DV_P_ZERO] + lbModel.fTemp[lbModel.DV_M_ZERO];
    uX = lbModel.fTemp[lbModel.DV_P_ZERO] - lbModel.fTemp[lbModel.DV_M_ZERO];

    t1 += lbModel.fTemp[lbModel.DV_ZERO_P] + lbModel.fTemp[lbModel.DV_ZERO_M];
    uY = lbModel.fTemp[lbModel.DV_ZERO_P] - lbModel.fTemp[lbModel.DV_ZERO_M];

    theta += t1; // changed here
    rho += t1;

    // 12 ops in FCC

    t1 = lbModel.fTemp[lbModel.DV_P_P] + lbModel.fTemp[lbModel.DV_M_M];
    t2 = lbModel.fTemp[lbModel.DV_P_P] - lbModel.fTemp[lbModel.DV_M_M];
    t1 += lbModel.fTemp[lbModel.DV_P_M] + lbModel.fTemp[lbModel.DV_M_P];
    t3 = lbModel.fTemp[lbModel.DV_P_M] - lbModel.fTemp[lbModel.DV_M_P];

    theta += 2.0 * t1;
    rho += t1;
    uX += t2 + t3;
    uY += t2 - t3;

    // 9 ops + 1 div = 19 ops
    // const ops does not contribute as it is evaluated once only in pgm lifetime
    const dataType1 oneByTwo = 1.0 / 2.0;
    dataType1 oneByRho = 1.0 / rho;
    uX *= oneByRho;
    uY *= oneByRho;
    theta = oneByTwo * (theta - rho * (uX * uX + uY * uY));
    theta *= oneByRho;
}

template <typename dataType1>
void getHydroMomentGrid(lbmD2Q9<dataType1> &lbModel,
                        field2D<dataType1, 9> &lbmGrid,
                        field2D<dataType1, 4> &fieldGrid)
{

    dataType1 t1, t2;

    // Initialize the field arrays
    for (int i2 = lbmGrid.n2Begin; i2 <= lbmGrid.n2End; i2++)
    {
        for (int i1 = lbmGrid.n1Begin; i1 <= lbmGrid.n1End; i1++)
        {
            fieldGrid(i1, i2, lbModel.RHO_INDEX) = 0.0;
            fieldGrid(i1, i2, lbModel.UX_INDEX) = 0.0;
            fieldGrid(i1, i2, lbModel.UY_INDEX) = 0.0;
            fieldGrid(i1, i2, lbModel.THETA_INDEX) = 0.0;
        }
    }

    //  Loop over field grid
    for (int i2 = lbmGrid.n2Begin; i2 <= lbmGrid.n2End; i2++)
    {
        for (int i1 = lbmGrid.n1Begin; i1 <= lbmGrid.n1End; i1++)
        {
            // 0 shell
            fieldGrid(i1, i2, lbModel.RHO_INDEX) =
                lbmGrid(i1, i2, lbModel.DV_ZERO_ZERO);

            // SC shell
            t1 = lbmGrid(i1, i2, lbModel.DV_P_ZERO) +
                 lbmGrid(i1, i2, lbModel.DV_M_ZERO);
            t2 = lbmGrid(i1, i2, lbModel.DV_P_ZERO) -
                 lbmGrid(i1, i2, lbModel.DV_M_ZERO);

            fieldGrid(i1, i2, lbModel.UX_INDEX) += t2;

            t1 += lbmGrid(i1, i2, lbModel.DV_ZERO_P) +
                  lbmGrid(i1, i2, lbModel.DV_ZERO_M);
            t2 = lbmGrid(i1, i2, lbModel.DV_ZERO_P) -
                 lbmGrid(i1, i2, lbModel.DV_ZERO_M);

            fieldGrid(i1, i2, lbModel.UY_INDEX) += t2;
            fieldGrid(i1, i2, lbModel.RHO_INDEX) += t1;
            fieldGrid(i1, i2, lbModel.THETA_INDEX) += t1;

            // FCC shell
            t1 = lbmGrid(i1, i2, lbModel.DV_P_P) + lbmGrid(i1, i2, lbModel.DV_M_M);
            t2 = lbmGrid(i1, i2, lbModel.DV_P_P) - lbmGrid(i1, i2, lbModel.DV_M_M);

            fieldGrid(i1, i2, lbModel.UX_INDEX) += t2;
            fieldGrid(i1, i2, lbModel.UY_INDEX) += t2;

            t1 += lbmGrid(i1, i2, lbModel.DV_P_M) + lbmGrid(i1, i2, lbModel.DV_M_P);
            t2 = lbmGrid(i1, i2, lbModel.DV_P_M) - lbmGrid(i1, i2, lbModel.DV_M_P);

            fieldGrid(i1, i2, lbModel.UX_INDEX) += t2;
            fieldGrid(i1, i2, lbModel.UY_INDEX) -= t2;
            fieldGrid(i1, i2, lbModel.RHO_INDEX) += t1;
            fieldGrid(i1, i2, lbModel.THETA_INDEX) += 2.0 * t1;
        }
    }

    for (int i2 = lbmGrid.n2Begin; i2 <= lbmGrid.n2End; i2++)
    {
        for (int i1 = lbmGrid.n1Begin; i1 <= lbmGrid.n1End; i1++)
        {
            dataType1 oneByRho = 1.0 / fieldGrid(i1, i2, lbModel.RHO_INDEX);
            fieldGrid(i1, i2, lbModel.UX_INDEX) *= oneByRho;
            fieldGrid(i1, i2, lbModel.UY_INDEX) *= oneByRho;
            fieldGrid(i1, i2, lbModel.THETA_INDEX) =
                0.5 * (fieldGrid(i1, i2, lbModel.THETA_INDEX) -
                       fieldGrid(i1, i2, lbModel.RHO_INDEX) *
                           (fieldGrid(i1, i2, lbModel.UX_INDEX) *
                                fieldGrid(i1, i2, lbModel.UX_INDEX) +
                            fieldGrid(i1, i2, lbModel.UY_INDEX) *
                                fieldGrid(i1, i2, lbModel.UY_INDEX)));
            fieldGrid(i1, i2, lbModel.THETA_INDEX) *= oneByRho;
        }
    }
}

template <typename dataType1>
void getFEq(lbmD2Q9<dataType1> &lbModel,
            dataType1 &rho,
            dataType1 &uX,
            dataType1 &uY)
{
    for (int dv = 0; dv < 9; dv++)
    {
        dataType1 dot1 = 0.0, dot2 = 0.0;
        dot1 = (lbModel.cX[dv] * uX + lbModel.cY[dv] * uY) * lbModel.oneByTheta0;
        dot2 = ((lbModel.cX[dv] * uX + lbModel.cY[dv] * uY) *
                    (lbModel.cX[dv] * uX + lbModel.cY[dv] * uY) -
                (uX * uX + uY * uY) * lbModel.theta0) *
               0.5 * lbModel.oneByTheta0Sqr;
        lbModel.fEq[dv] = lbModel.wt[dv] * rho * (1.0 + dot1 + dot2);
    }
}

template <typename dataType1>
void getFEqIsoSecond(lbmD2Q9<dataType1> &lbModel,
                     dataType1 &rho,
                     dataType1 &uX,
                     dataType1 &uY)
{

    dataType1 t2 = (uX * uX + uY * uY) * 0.5 * lbModel.oneByTheta0;
    dataType1 dot1 = lbModel.oneByTheta0 * uX;
    dataType1 dot2 = lbModel.oneByTheta0 * uY;

    for (int dv = 0; dv < lbModel.dvN; dv++)
        lbModel.fEq[dv] = rho * lbModel.wt[dv];

    lbModel.fEq[lbModel.DV_ZERO_ZERO] *= (1.0 - t2);

    dataType1 t1 = 0.5 * dot1 * dot1;
    lbModel.fEq[lbModel.DV_P_ZERO] *= (1.0 + dot1 + t1 - t2);
    lbModel.fEq[lbModel.DV_M_ZERO] *= (1.0 - dot1 + t1 - t2);

    t1 = 0.5 * dot2 * dot2;
    lbModel.fEq[lbModel.DV_ZERO_P] *= (1.0 + dot2 + t1 - t2);
    lbModel.fEq[lbModel.DV_ZERO_M] *= (1.0 - dot2 + t1 - t2);

    dataType1 dot1Plusdot2 = dot1 + dot2;
    dataType1 dot1Minusdot2 = dot1 - dot2;

    t1 = 0.5 * dot1Plusdot2 * dot1Plusdot2;
    lbModel.fEq[lbModel.DV_P_P] *= (1.0 + dot1Plusdot2 + t1 - t2);
    lbModel.fEq[lbModel.DV_M_M] *= (1.0 - dot1Plusdot2 + t1 - t2);

    t1 = 0.5 * dot1Minusdot2 * dot1Minusdot2;
    lbModel.fEq[lbModel.DV_M_P] *= (1.0 - dot1Minusdot2 + t1 - t2);
    lbModel.fEq[lbModel.DV_P_M] *= (1.0 + dot1Minusdot2 + t1 - t2);
}

template <typename dataType>
void advectionD2Q9(field2D<dataType, 9> &lbmGrid,
                   lbmD2Q9<dataType> &d2q9Model)
{

    for (int i2 = lbmGrid.n2Begin; i2 <= lbmGrid.n2End; i2++)
    {
        for (int i1 = lbmGrid.n1Begin; i1 <= lbmGrid.n1End; i1++)
        {
            lbmGrid(i1, i2, d2q9Model.DV_ZERO_M) = lbmGrid(i1, i2 + 1, d2q9Model.DV_ZERO_M);
            lbmGrid(i1, i2, d2q9Model.DV_M_ZERO) = lbmGrid(i1 + 1, i2, d2q9Model.DV_M_ZERO);
            lbmGrid(i1, i2, d2q9Model.DV_M_M)    = lbmGrid(i1 + 1, i2 + 1, d2q9Model.DV_M_M);
            lbmGrid(i1, i2, d2q9Model.DV_P_M)    = lbmGrid(i1 - 1, i2 + 1, d2q9Model.DV_P_M);
        }
    }

    for (int i2 = lbmGrid.n2End; i2 >= lbmGrid.n2Begin; i2--)
    {
        for (int i1 = lbmGrid.n1End; i1 >= lbmGrid.n1Begin; i1--)
        {
            lbmGrid(i1, i2, d2q9Model.DV_ZERO_P) = lbmGrid(i1, i2 - 1, d2q9Model.DV_ZERO_P);
            lbmGrid(i1, i2, d2q9Model.DV_P_ZERO) = lbmGrid(i1 - 1, i2, d2q9Model.DV_P_ZERO);
            lbmGrid(i1, i2, d2q9Model.DV_M_P)    = lbmGrid(i1 + 1, i2 - 1, d2q9Model.DV_M_P);
            lbmGrid(i1, i2, d2q9Model.DV_P_P)    = lbmGrid(i1 - 1, i2 - 1, d2q9Model.DV_P_P);
        }
    }
}

template <typename dataType>
void getDenField(field2D<dataType, 9> &lbmGrid,
                 lbmD2Q9<dataType> &d2q9Model,
                 dataType dt,
                 field2D<dataType, 2> denField)
{
    dataType rho;
    for (int i2 = lbmGrid.n2Begin; i2 <= lbmGrid.n2End; i2++)
    {
        for (int i1 = lbmGrid.n1Begin; i1 <= lbmGrid.n1End; i1++)
        {
            rho = 0.0;
            for (int dv = 0; dv < d2q9Model.dvN; dv++)
                rho += lbmGrid(i1, i2, dv);
            denField(i1, i2, 0) = rho;
        }
    }

    denField.makePeriodicX();
    denField.makePeriodicY();
}

template <typename dataType>
void getLaplacianDenField(field2D<dataType, 9> &lbmGrid,
                 lbmD2Q9<dataType> &d2q9Model,
                 dataType dt,
                 field2D<dataType, 2> denField){

    dataType coeff = 2.0/(dt*dt*d2q9Model.theta0);


    for (int i2 = lbmGrid.n2Begin; i2 <= lbmGrid.n2End; i2++)
    {
        for (int i1 = lbmGrid.n1Begin; i1 <= lbmGrid.n1End; i1++)
        {
            denField(i1, i2, 1) = 0.0;

            for (int dv = 0; dv < d2q9Model.dvN; dv++)
                denField(i1, i2, 1) += d2q9Model.wt[dv]*denField(i1 + (int)d2q9Model.cX[dv], i2 + (int)d2q9Model.cY[dv],0 ); 
            
            denField(i1, i2, 1) = coeff*(denField(i1, i2, 1) - denField(i1 ,i2 , 0));

            // std::cout<<denField(i1,i2,1) <<std::endl;

        }
    }

    denField.makePeriodicX();
    denField.makePeriodicY();
}





template <typename dataType>
void getMuNid(field2D<dataType, 9> &lbmGrid,
                 lbmD2Q9<dataType> &d2q9Model,
                 field2D<dataType, 2> denField,
                 field2D<dataType, 1> muNid, dataType a, dataType b, dataType kappa ,dataType dt){

    

    for (int i2 = lbmGrid.n2Begin; i2 <= lbmGrid.n2End; i2++)
    {
        for (int i1 = lbmGrid.n1Begin; i1 <= lbmGrid.n1End; i1++)
        {
            muNid(i1 ,i2 , 0)  = -d2q9Model.theta0*log(1.0 - denField(i1,i2,0) *b );
            muNid(i1 ,i2 , 0) += denField(i1,i2,0)*b*d2q9Model.theta0/(1.0 - denField(i1,i2,0) *b ) ;
            muNid(i1 ,i2 , 0) -= 2.0*denField(i1,i2,0) *a;
            muNid(i1 ,i2 , 0) -= kappa*denField(i1,i2,1);
        }
    }

    muNid.makePeriodicX();
    muNid.makePeriodicY();
   
}



template <typename dataType>
void getFNid(field2D<dataType, 9> &lbmGrid,
                 lbmD2Q9<dataType> &d2q9Model,
                 field2D<dataType, 2> denField,
                 field2D<dataType, 1> FNid, dataType a, dataType b, dataType kappa ,dataType dt){

    

    for (int i2 = lbmGrid.n2Begin; i2 <= lbmGrid.n2End; i2++)
    {
        for (int i1 = lbmGrid.n1Begin; i1 <= lbmGrid.n1End; i1++)
        {
            FNid(i1 ,i2 , 0)  = -a*denField(i1,i2,0)*denField(i1,i2,0) ;
            FNid(i1 ,i2 , 0) -= denField(i1,i2,0)*d2q9Model.theta0*log(1.0 - denField(i1,i2,0) *b );
        }
    }

    FNid.makePeriodicX();
    FNid.makePeriodicY();
   
}


template <typename dataType>
void collideD2Q9(field2D<dataType, 9> &lbmGrid,
                 lbmD2Q9<dataType> &d2q9Model,
                 dataType beta,
                 dataType dt,
                 field2D<dataType, 2> forceField)
{

    dataType rho, uX, uY, theta;
    for (int i2 = lbmGrid.n2Begin; i2 <= lbmGrid.n2End; i2++)
    {
        for (int i1 = lbmGrid.n1Begin; i1 <= lbmGrid.n1End; i1++)
        {
            for (int dv = 0; dv < d2q9Model.dvN; dv++)
            {
                d2q9Model.fTemp[dv] = lbmGrid(i1, i2, dv);
            }

            getHydroMomentV1(d2q9Model, rho, uX, uY, theta);

            uX += forceField(i1, i2, 0) * dt * 0.5;
            uY += forceField(i1, i2, 1) * dt * 0.5;

            getFEqIsoSecond(d2q9Model, rho, uX, uY);

            for (int dv = 0; dv < d2q9Model.dvN; dv++)
            {
                dataType Fi = d2q9Model.wt[dv] * d2q9Model.oneByTheta0 *
                                (forceField(i1, i2, 0) * d2q9Model.cX[dv] + forceField(i1, i2, 1) * d2q9Model.cY[dv]);
                lbmGrid(i1, i2, dv) += 2.0 * beta * (d2q9Model.fEq[dv] - lbmGrid(i1, i2, dv)) + (1.0 - beta) * dt * Fi;
            }
        }
    }
}


template <typename dataType>
void initializeTaylorGreen(field2D<dataType, 9> &lbmGrid,
                           lbmD2Q9<dataType> &d2q9Model,
                           dataType u_ref)
{

    dataType rho = 1.0;

    for (int i2 = lbmGrid.n2Begin; i2 <= lbmGrid.n2End; i2++)
    {
        for (int i1 = lbmGrid.n1Begin; i1 <= lbmGrid.n1End; i1++)
        {

            dataType x =
                (static_cast<dataType>(i1) / static_cast<dataType>(lbmGrid.n1)) * 2.0 *
                M_PI;
            dataType y =
                (static_cast<dataType>(i2) / static_cast<dataType>(lbmGrid.n2)) * 2.0*
                M_PI;

            dataType uX = u_ref * sin(x) * cos(y), uY = u_ref * -cos(x) * sin(y);

            getFEqIsoSecond(d2q9Model, rho, uX, uY);
            // getFEq(d2q9Model, rho, uX, uY);

            for (int dv = 0; dv < d2q9Model.dvN; dv++)
            {
                lbmGrid(i1, i2, dv) = d2q9Model.fEq[dv];
            }
        }
    }
}

template <typename dataType>
void initializeFEq(field2D<dataType, 9> &lbmGrid,
                   lbmD2Q9<dataType> &d2q9Model)
{

    dataType rho = 1.0, uX = 0.0, uY = 0.0;

    for (int i2 = lbmGrid.n2Begin; i2 <= lbmGrid.n2End; i2++)
    {
        for (int i1 = lbmGrid.n1Begin; i1 <= lbmGrid.n1End; i1++)
        {

            // getFEq(d2q9Model, rho,uX,uY);
            getFEqIsoSecond(d2q9Model, rho, uX, uY);

            for (int dv = 0; dv < d2q9Model.dvN; dv++)
                lbmGrid(i1, i2, dv) = d2q9Model.fEq[dv];
        }
    }
}



template <typename dataType>
void initialize1DInterface(field2D<dataType, 9> &lbmGrid,
                   lbmD2Q9<dataType> &d2q9Model, dataType rhoLiq, dataType rhoGas)
{

    dataType rho = 1.0, uX = 0.0, uY = 0.0;

    for (int i2 = lbmGrid.n2Begin; i2 <= lbmGrid.n2End; i2++)
    {
        for (int i1 = lbmGrid.n1Begin; i1 <= lbmGrid.n1End; i1++)
        {

            dataType x = (static_cast<dataType>(i1) / static_cast<dataType>(lbmGrid.n1));
            dataType y = (static_cast<dataType>(i2) / static_cast<dataType>(lbmGrid.n1));
            
            dataType radius = sqrt(x*x + y*y);
            rho = (rhoLiq + rhoGas)*0.5 + (rhoLiq - rhoGas)*0.5*tanh(x); 


            getFEqIsoSecond(d2q9Model, rho, uX, uY);

            for (int dv = 0; dv < d2q9Model.dvN; dv++)
                lbmGrid(i1, i2, dv) = d2q9Model.fEq[dv];
        }
    }
}

template <typename dataType>
void calculateForce(field2D<dataType, 9> &lbmGrid,
                    lbmD2Q9<dataType> &d2q9Model,
                    field2D<dataType, 2> &forceField,
                    field2D<dataType, 2> &denField,
                    field2D<dataType, 1> &muNid,
                    dataType a, dataType b, dataType kappa,dataType dt  )
{

    getDenField(lbmGrid,d2q9Model,dt,denField);
    getLaplacianDenField(lbmGrid,d2q9Model,dt,denField);
    getMuNid(lbmGrid,d2q9Model,denField,muNid,a, b, kappa,dt );

    dataType coeff = 1.0/(dt*d2q9Model.theta0);


    for (int i2 = lbmGrid.n2Begin; i2 <= lbmGrid.n2End; i2++)
    {
        for (int i1 = lbmGrid.n1Begin; i1 <= lbmGrid.n1End; i1++)
        {

            forceField(i1,i2,0) = 0.0;
            forceField(i1,i2,1) = 0.0;
            for (int dv = 0; dv < d2q9Model.dvN; dv++){
                forceField(i1,i2,0) += d2q9Model.wt[dv]*d2q9Model.cX[dv]*muNid(i1 + (int)d2q9Model.cX[dv], i2 + (int)d2q9Model.cY[dv],0 );
                forceField(i1,i2,1) += d2q9Model.wt[dv]*d2q9Model.cY[dv]*muNid(i1 + (int)d2q9Model.cX[dv], i2 + (int)d2q9Model.cY[dv],0 ); 
            }
            
            forceField(i1,i2,0) = -coeff*(forceField(i1,i2,0));
            forceField(i1,i2,1) = -coeff*(forceField(i1,i2,1));

        }
    }

}


template <typename dataType>
void calculateForce1(field2D<dataType, 9> &lbmGrid,
                    lbmD2Q9<dataType> &d2q9Model,
                    field2D<dataType, 2> &forceField,
                    field2D<dataType, 2> &denField,
                    field2D<dataType, 1> &muNid,
                    field2D<dataType, 1> &FNid,

                    dataType a, dataType b, dataType kappa,dataType dt  )
{

    getDenField(lbmGrid,d2q9Model,dt,denField);
    getLaplacianDenField(lbmGrid,d2q9Model,dt,denField);
    getMuNid(lbmGrid,d2q9Model,denField,muNid,a, b, kappa,dt );
    getFNid(lbmGrid,d2q9Model,denField,FNid,a, b, kappa,dt );

    dataType coeff = 1.0/(dt*d2q9Model.theta0);



    for (int i2 = lbmGrid.n2Begin; i2 <= lbmGrid.n2End; i2++)
    {
        for (int i1 = lbmGrid.n1Begin; i1 <= lbmGrid.n1End; i1++)
        {
            
            dataType oneByRho = 1.0/denField(i1, i2, 0);

            forceField(i1,i2,0) = 0.0;
            forceField(i1,i2,1) = 0.0;

            for (int dv = 0; dv < d2q9Model.dvN; dv++){
                forceField(i1,i2,0) += d2q9Model.wt[dv]*d2q9Model.cX[dv]*muNid(i1 + (int)d2q9Model.cX[dv], i2 + (int)d2q9Model.cY[dv],0 );
                forceField(i1,i2,1) += d2q9Model.wt[dv]*d2q9Model.cY[dv]*muNid(i1 + (int)d2q9Model.cX[dv], i2 + (int)d2q9Model.cY[dv],0 ); 
            }

            for (int dv = 0; dv < d2q9Model.dvN; dv++){
                forceField(i1,i2,0) += oneByRho* muNid(i1, i2, 0)*d2q9Model.wt[dv]*d2q9Model.cX[dv]*denField(i1 + (int)d2q9Model.cX[dv], i2 + (int)d2q9Model.cY[dv],0 );
                forceField(i1,i2,1) += oneByRho* muNid(i1, i2, 0)*d2q9Model.wt[dv]*d2q9Model.cY[dv]*denField(i1 + (int)d2q9Model.cX[dv], i2 + (int)d2q9Model.cY[dv],0 ); 
            }

            for (int dv = 0; dv < d2q9Model.dvN; dv++){
                forceField(i1,i2,0) -= oneByRho*d2q9Model.wt[dv]*d2q9Model.cX[dv]*FNid(i1 + (int)d2q9Model.cX[dv], i2 + (int)d2q9Model.cY[dv],0 );
                forceField(i1,i2,1) -= oneByRho*d2q9Model.wt[dv]*d2q9Model.cY[dv]*FNid(i1 + (int)d2q9Model.cX[dv], i2 + (int)d2q9Model.cY[dv],0 ); 
            }
            
            forceField(i1,i2,0) = -coeff*(forceField(i1,i2,0));
            forceField(i1,i2,1) = -coeff*(forceField(i1,i2,1));

        }
    }

}




template <typename dataType>
void calculateMass(field2D<dataType, 9> &lbmGrid)
{
    dataType mass = 0.0;
    for (int i2 = lbmGrid.n2Begin; i2 <= lbmGrid.n2End; i2++)
    {
        for (int i1 = lbmGrid.n1Begin; i1 <= lbmGrid.n1End; i1++)
        {
            for (int dv = 0; dv < 9; dv++)
            {
                mass += lbmGrid(i1, i2, dv);
            }
        }
    }
    std::cout << "total mass is: " << mass << std::endl;
}

template <typename dataType>
void printVtk(field2D<dataType, 9> &lbmGrid,
              lbmD2Q9<dataType> &d2q9Model,
              dataType dt,
              field2D<dataType, 2> &forceField,
              int timeStep)
{
    std::ofstream file;
    char fileName[250];
    dataType ux, uy;
    sprintf(fileName, "./Result/velocity_%d.vtk", timeStep);
    file.open(fileName);

    file << "# vtk DataFile Version 3.0\nVelocity\nASCII\nDATASET "
            "STRUCTURED_POINTS"
         << std::endl;

    file << "DIMENSIONS " << lbmGrid.n1 << " " << lbmGrid.n2 << " " << 1
         << std::endl;
    file << "ORIGIN " << 0 << " " << 0 << " " << 0 << std::endl;
    file << "SPACING " << 1 << " " << 1 << " " << 1 << std::endl;
    file << "POINT_DATA " << 1 * lbmGrid.n1 * 1 * lbmGrid.n2 << std::endl;
    file << "SCALARS density double 1\nLOOKUP_TABLE default" << std::endl;

    dataType rho, uX, uY, theta;

    for (int i2 = lbmGrid.n2Begin; i2 <= lbmGrid.n2End; i2++)
    {
        for (int i1 = lbmGrid.n1Begin; i1 <= lbmGrid.n1End; i1++)
        {

            for (int dv = 0; dv < d2q9Model.dvN; dv++)
            {
                d2q9Model.fTemp[dv] = lbmGrid(i1, i2, dv);

            }

            getHydroMomentV1(d2q9Model, rho, uX, uY, theta);
            file<<rho<<std::endl;
        }
    }

    file << "VECTORS velocity double" << std::endl;
    for (int i2 = lbmGrid.n2Begin; i2 <= lbmGrid.n2End; i2++)
    {
        for (int i1 = lbmGrid.n1Begin; i1 <= lbmGrid.n1End; i1++)
        {

            for (int dv = 0; dv < d2q9Model.dvN; dv++)
            {
                d2q9Model.fTemp[dv] = lbmGrid(i1, i2, dv);
            }
            getHydroMomentV1(d2q9Model, rho, uX, uY, theta);
            uX += 0.5 * dt * forceField(i1, i2, 0);
            uY += 0.5 * dt * forceField(i1, i2, 1);

            file << uX << " " << uY << " " << 0.0 << std::endl;
        }
    }
}

template <typename dataType>
void printVelocityToTxt(field2D<dataType, 9> &lbmGrid,
                        lbmD2Q9<dataType> &d2q9Model, int timeStep)
{
    // Open a file to write the velocity components
    std::ofstream file;
    char fileName[250];
    sprintf(fileName, "./ResultTxt/velocity_%d_a.txt",
            timeStep); // Create a file name based on the timeStep
    file.open(fileName);

    if (!file.is_open())
    {
        std::cerr << "Error: Unable to open file " << fileName << std::endl;
        return;
    }

    // Write the header
    file << "# i1 i2 uX uY" << std::endl;

    // Variables to hold density and velocity components
    dataType rho, uX, uY, theta;

    // Loop over the grid and compute the velocity components at each point
    for (int i2 = lbmGrid.n2Begin; i2 <= lbmGrid.n2End; i2++)
    {
        for (int i1 = lbmGrid.n1Begin; i1 <= lbmGrid.n1End; i1++)
        {
            // Get the hydrodynamic moments (rho, uX, uY, theta) for the current grid
            // point
            getHydroMomentV1(d2q9Model, rho, uX, uY, theta);

            // Write the grid point indices (i1, i2) and the velocity components (uX,
            // uY)
            file << i1 << " " << i2 << " " << uX << " " << uY << std::endl;
        }
    }

    file.close(); // Close the file
}
