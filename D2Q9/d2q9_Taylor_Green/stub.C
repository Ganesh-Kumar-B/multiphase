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

#include "field2D.h"
#include "D2Q9.h"
#include "boundary.cpp"
#include "enstrophy.cpp"
#include "field2D.h"
#include <math.h>

typedef double myReal;
int main() {

    // Parameter related to grid
    int nX(256);
    int nY(256);
    myReal latticeSpeed = 1.0;
    myReal boxLength = 2.0*M_PI;
    myReal dx = boxLength / nY;
    myReal dt = dx / latticeSpeed;

    // defining discrete velocity model
    lbmD2Q9<myReal> d2q9Model(latticeSpeed); 

    d2q9Model.setModelParameters();

    // Parameters related to flow
    myReal knudsenNum = 0.001;
    myReal theta0 = d2q9Model.theta0;
    // myReal tau = knudsenNum * boxLength / sqrt(theta0);
    myReal ReynoldsNumber = 10;

    myReal Ma = 0.05;
    myReal cs = std::sqrt(theta0);
    myReal target_U0 = 0.05;
    

    myReal kinematicVisc = target_U0 * boxLength/ReynoldsNumber;
    myReal tau = kinematicVisc / theta0;


    myReal uRef = (ReynoldsNumber * kinematicVisc) / boxLength;
    std::cout << "u ref: " << uRef << std::endl;

    myReal tauNdim = tau / dt;

    myReal beta = 1.0 / (2.0 * tauNdim + 1.0);
    // myReal beta = 0.9;
    std::cout<<"beta    "<<beta<<std::endl;

    std::cout << "kinematic Viscosity = " << kinematicVisc << std::endl;

    // defining grid
    // defines grid to store populations with 1 ghost index padding
    field2D<myReal, 9> lbmGrid(nX, nY, 1);

    // defines grid to store forces pointwise
    field2D<myReal, 2> forceField(nX, nY, 1); 

    // defines grid to store hydromoments
    field2D<myReal, 4> fieldGrid(nX, nY, 1); 

    // defines grid to store rho and laplacianRho
    field2D<myReal,2> denField  (nX,nY,1);
    field2D<myReal,1> muNid     (nX,nY,1);
    field2D<myReal,1> FNid      (nX,nY,1);


    //defining the VdW paramters and temperature
   

    initializeTaylorGreen(lbmGrid,d2q9Model,uRef);
    // initializeFEq(lbmGrid, d2q9Model);
    myReal radiusC = 0.20;
    std::cout<<"Radius  :"<<radiusC<<std::endl;

    // initializeTanhxRadial(lbmGrid, d2q9Model, rhoLiq, rhoGas,radiusC);

    // initialize1DInterface(lbmGrid, d2q9Model, rhoLiq, rhoGas);

    getHydroMomentGrid(d2q9Model, lbmGrid, fieldGrid);
    calculateMass(lbmGrid);

    std::string name="Results_";
    printVtk(lbmGrid,d2q9Model,dt, forceField, uRef,name,0);



    double convectionTime = (double)nX / (uRef); // based on ref length?
    int iterations = 5000; (int)(15.0 * convectionTime / dt);

    myReal time = 0.0;

    for (int timeStep = 1; timeStep <= iterations; timeStep++) {

        collideD2Q9(lbmGrid, d2q9Model, beta, dt, forceField);

        lbmGrid.makePeriodicX();
        lbmGrid.makePeriodicY();

        // periodicX(lbmGrid, d2q9Model);
        // periodicY(lbmGrid, d2q9Model);

        advectionD2Q9(lbmGrid, d2q9Model);


        time += dt;
        if (timeStep % 20 == 0) {
            printVtk(lbmGrid, d2q9Model, dt, forceField,uRef,name, timeStep);
            std::cout << " time " << time << " ";
            calculateMass(lbmGrid);
            std::cout << "\n";
            // printVelocityToTxt(lbmGrid, d2q9Model, timeStep);
        }
        // if(timeStep % 2 == 0)
        //     printEnstrophy(d2q9Model, lbmGrid, beta, dt, timeStep,
        //     convectionTime, u_ref);
    }

    // printVelocityToTxt(lbmGrid, d2q9Model, iterations);
    return 0;
}
