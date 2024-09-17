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
    int nX(128);
    int nY(128);
    myReal latticeSpeed = 1.0;
    myReal boxLength = nY;
    myReal dx = boxLength / nY;
    myReal dt = dx / latticeSpeed;

    // defining discrete velocity model
    lbmD2Q9<myReal> d2q9Model(latticeSpeed); 

    d2q9Model.setModelParameters();

    // Parameters related to flow
    myReal knudsenNum = 0.001;
    myReal theta0 = d2q9Model.theta0;
    // myReal tau = knudsenNum * boxLength / sqrt(theta0);
    myReal ReynoldsNumber = 200;

    myReal Ma = 0.05;
    myReal cs = std::sqrt(theta0);
    myReal target_U0 = Ma * cs;

    myReal external_force = 8.0*target_U0*target_U0/(3.0*ReynoldsNumber*boxLength/2.0);
    myReal kinematicVisc = target_U0 * boxLength/ReynoldsNumber;
    myReal tau = kinematicVisc / theta0;


    myReal g = external_force;
    std::cout<<"g   "<<g<<std::endl;

    myReal u_ref = (ReynoldsNumber * kinematicVisc) / boxLength;
    std::cout << "u ref: " << u_ref << std::endl;

    myReal tauNdim = tau / dt;

    myReal beta = 1.0 / (2.0 * tauNdim + 1.0);
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
    myReal Tr = 0.95;
    myReal rhoCritical = 1.0, TCritical = d2q9Model.theta0/Tr;
    myReal b = 1.0/(3.0*rhoCritical), a = b*TCritical*27.0/8.0;
    myReal kappa = 0.001;
    myReal rhoLiq = 1.6165, rhoGas = 0.4997;

    // initializeTaylorGreen(lbmGrid,d2q9Model,u_ref);
    // initializeFEq(lbmGrid, d2q9Model);
    initialize1DInterface(lbmGrid, d2q9Model, rhoLiq, rhoGas);

    getHydroMomentGrid(d2q9Model, lbmGrid, fieldGrid);

    calculateMass(lbmGrid);
    printVtk(lbmGrid,d2q9Model,dt, forceField,0);



    double convectionTime = (double)nX / (u_ref); // based on ref length?
    int iterations = 200000; (int)(15.0 * convectionTime / dt);

    myReal time = 0.0;

    for (int timeStep = 1; timeStep <= iterations; timeStep++) {

        // calculateForce  (lbmGrid, d2q9Model, forceField, denField,muNid, a, b, kappa , dt);   // g_alpha = \rho * grad Munid 
        calculateForce1 (lbmGrid, d2q9Model, forceField, denField,muNid,FNid, a, b, kappa , dt); // g_alpha = \rho * grad Munid

        collideD2Q9(lbmGrid, d2q9Model, beta, dt, forceField);

        lbmGrid.makePeriodicX();
        lbmGrid.makePeriodicY();

        // periodicX(lbmGrid, d2q9Model);
        // periodicY(lbmGrid, d2q9Model);

        advectionD2Q9(lbmGrid, d2q9Model);


        time += dt;
            if (timeStep % 5000 == 0) {
            printVtk(lbmGrid, d2q9Model, dt, forceField, timeStep);
            std::cout << " time " << timeStep << " ";
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
