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


#include <sstream>
// Enstrophy
template <typename dataType>
void getEnstrophyMoments(lbmD2Q9<dataType> &lbModel,
                         field2D<dataType, 9> &lbmGrid,
                         int i1, int i2,
                         double &rho, double &ux, double &uy, double &temp,
                         double &pxx, double &pyy, double &pxy, double &sXX,
                         double &sYY, double &sXY, double beta, double dt) {
  // dataType lbModel.fTemp[35];

  for (int dv = 0; dv < lbModel.dvN; dv++) {
    lbModel.fTemp[dv] = lbmGrid(i1, i2, dv);
  }

  // getHydroMomentV1()

  double rhoZero, rhoSC, rhoFCC, uxSC, uxFCC, uySC, uyFCC, uzSC, uzFCC;

  rhoZero = lbModel.fTemp[lbModel.DV_ZERO_ZERO];
  rhoSC = lbModel.fTemp[lbModel.DV_P_ZERO] + lbModel.fTemp[lbModel.DV_M_ZERO] +
          lbModel.fTemp[lbModel.DV_ZERO_P] + lbModel.fTemp[lbModel.DV_ZERO_M];
  rhoFCC = lbModel.fTemp[lbModel.DV_P_P] + lbModel.fTemp[lbModel.DV_M_M] +
           lbModel.fTemp[lbModel.DV_P_M] + lbModel.fTemp[lbModel.DV_M_P];

  uxSC = (lbModel.fTemp[lbModel.DV_P_ZERO] - lbModel.fTemp[lbModel.DV_M_ZERO]);
  uySC = (lbModel.fTemp[lbModel.DV_ZERO_P] - lbModel.fTemp[lbModel.DV_ZERO_M]);

  uxFCC = lbModel.fTemp[lbModel.DV_P_P] - lbModel.fTemp[lbModel.DV_M_M] +
          lbModel.fTemp[lbModel.DV_P_M] - lbModel.fTemp[lbModel.DV_M_P];
  uyFCC = lbModel.fTemp[lbModel.DV_P_P] - lbModel.fTemp[lbModel.DV_M_M] -
          lbModel.fTemp[lbModel.DV_P_M] + lbModel.fTemp[lbModel.DV_M_P];

  pxy = lbModel.fTemp[lbModel.DV_P_P] + lbModel.fTemp[lbModel.DV_M_M] -
        lbModel.fTemp[lbModel.DV_P_M] - lbModel.fTemp[lbModel.DV_M_P];

  pxx = (lbModel.fTemp[lbModel.DV_P_ZERO] + lbModel.fTemp[lbModel.DV_M_ZERO]) +
        lbModel.fTemp[lbModel.DV_P_P] + lbModel.fTemp[lbModel.DV_M_M] +
        lbModel.fTemp[lbModel.DV_P_M] + lbModel.fTemp[lbModel.DV_M_P];
  pyy = (lbModel.fTemp[lbModel.DV_ZERO_P] + lbModel.fTemp[lbModel.DV_ZERO_M]) +
        lbModel.fTemp[lbModel.DV_P_P] + lbModel.fTemp[lbModel.DV_M_M] +
        lbModel.fTemp[lbModel.DV_P_M] + lbModel.fTemp[lbModel.DV_M_P];

  rho = rhoZero + rhoSC + rhoFCC;
  ux = uxSC + uxFCC;
  uy = uySC + uyFCC;

  double rhoInv = 1.0 / rho;
  ux *= rhoInv;
  uy *= rhoInv;

  temp = rhoSC + rhoFCC * 2.0;
  temp = (temp - rho * (ux * ux + uy * uy)) * rhoInv * 0.5;

  double factor = beta / (dt * rho * temp);

  sXX = (rho * ux * ux + rho * temp - pxx) * factor;
  sYY = (rho * uy * uy + rho * temp - pyy) * factor;
  sXY = (rho * ux * uy - pxy) * factor;
}

template <typename dataType>
void printEnstrophy(lbmD2Q9<dataType> &lbModel,
                    field2D<dataType, 9> &lbmGrid,
                    double beta, double dt, int t,
                    double convectionTime,
                    double U0) {

  double rho, ux, uy, uz, theta, pxx, pyy, pzz, pxy, pxz, pyz, sXX, sYY, sZZ,
      sXY, sYZ, sXZ;

  double energy = 0.0, rhosum = 0.0, internalEnergy(0.0), kineticEnergy(0.0),
         thetasum(0.0), enstrophySum(0.0), enstrophy(0.0), enstrophyMax(0.0);

  for (int i2 = lbmGrid.n2Begin; i2 <= lbmGrid.n2End; i2++) {
    for (int i1 = lbmGrid.n1Begin; i1 <= lbmGrid.n1End; i1++) {
      getEnstrophyMoments(lbModel, lbmGrid, i1, i2, rho, ux, uy, theta, pxx,
                          pyy, pxy, sXX, sYY, sXY, beta, dt);
      rhosum += rho;
      thetasum += theta;
      internalEnergy += 1.5 * rho * theta;
      kineticEnergy += rho * (ux * ux + uy * uy) * 0.5;
      enstrophy = sXX * sXX + sYY * sYY + 2.0 * sXY * sXY;
      enstrophySum += enstrophy;
    }
  }

  int totalLatticePoints = lbmGrid.numPoints;

  kineticEnergy /= U0 * U0 * totalLatticePoints;
  internalEnergy /= U0 * U0 * totalLatticePoints;
  enstrophySum /= U0 * U0 * totalLatticePoints;
  rhosum /= totalLatticePoints;
  thetasum /= totalLatticePoints;

  // std::cout<<"time: "<<t<<"  mass: "<<rhosum<<" "<<rhosum * 2.0 *
  // totalLatticePoints<<std::endl;

  // Calculate t/convectionTime
  double time = static_cast<dataType>(t) * dt / (1.0 / U0);

  // Create a filename based on Re, Mach, and globalSize
  std::ostringstream filename;
  filename << "enstrophy_data.txt";

  // Convert the filename to a std::string
  std::string filenameStr = filename.str();

  // Open the file in append mode
  std::ofstream outputFile(filenameStr.c_str(), std::ios_base::app);
  if (!outputFile.is_open()) {
    std::cerr << "Error: Unable to open file: " << filename.str() << std::endl;
    return;
  }

  // Write to file
  outputFile << time << " "
             << rhosum << " "
             << thetasum / lbModel.theta0 << " "
             << internalEnergy << "   "
             << kineticEnergy << "   "
             << internalEnergy + kineticEnergy
             << "  " << enstrophySum
             << std::endl;

  outputFile.close(); // Close the file
}