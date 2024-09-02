#include"field2D.h"
#include"D2Q9.h"
#include<math.h>
typedef  float myReal;
int main(){
    
    
    // Parameter related to grid
    int nX(128);
    int nY(128);
    myReal latticeSpeed = 1.0;
    myReal boxLength = 128.0;
    myReal dx = boxLength/nX;
    myReal dt = dx/latticeSpeed;

    
    myReal knudsenNum = 0.001;
    
    myReal theta0 = latticeSpeed * latticeSpeed / 3.0;
    myReal tau    =  knudsenNum * boxLength/sqrt(theta0);
    myReal kinematicVisc = tau*theta0;
    myReal ReynoldsNumber = 200;


    myReal g = 0.0001;
    
    myReal u_ref = (ReynoldsNumber*kinematicVisc)/ boxLength;
    std::cout<<"kinematic Viscosity = "<<kinematicVisc<<std::endl;
    myReal tauNdim = tau/dt;
    
    myReal beta  =  1.0/(2.0*tauNdim+1.0);
         
        
    
    field2D<myReal, 9> lbmGrid(nX,nY,1);
    
    // To store forces pointwise
    field2D<myReal,2> forceField(nX,nY,1);

    // to store rho and laplacianRho
    field2D<myReal,2> denField(nX,nY,1);
    
    lbmD2Q9<myReal>  d2q9Model(latticeSpeed);

    // initializeTaylorGreen(lbmGrid,d2q9Model,u_ref);
    initializeFEq(lbmGrid,d2q9Model);

    printVtk(lbmGrid,d2q9Model,denField,0.0);

    // calculateForce(lbmGrid,   d2q9Model, forceField, g);
     
    myReal time =0.0;
    for(int timeStep = 1; timeStep <= 5;timeStep++){
        

        collideD2Q9(lbmGrid,   d2q9Model,   beta,  dt,forceField );
        lbmGrid.makePeriodicX();
        lbmGrid.makePeriodicY();

        advectionD2Q9(lbmGrid, d2q9Model );
      
        time += dt;
        if(timeStep % 1 == 0)
            printVtk(lbmGrid,d2q9Model,denField,timeStep);
    }
    return 0;
}
