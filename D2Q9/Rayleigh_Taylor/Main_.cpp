#include<iostream>
#include<fstream>
#include<cmath>
#include"collison.h"
#include"advection.h"
#include"print.h"
#include<iomanip>
#include<sstream>





int main()
{

    int Nx = 64 ;int Ny = 128;


    Grid_N_C_2D<real> grid                  (Nx,Ny,1,9);
    Grid_N_C_2D<real> Force                 (Nx,Ny,1,2);




    lbmD2Q9<real> d2q9(1.0,(1.0/3.0));
    
    real cs = sqrt(d2q9.theta0);
    std::cout<<"theta= "<<d2q9.theta0<<std::endl;




    real Re = 2048;
    std::cout<<"Re      = "<<Re<<std::endl;
    real L  = Nx;


    real u0 = 0.2;
    std::cout<<"u0      = "<<u0<<std::endl;

    real g = (u0*u0)/L;
    std::cout<<"g       = "<<g<<std::endl;


    real Kin_Vis = u0*(L)/Re;
    std::cout<<"Kin_Vis = "<<Kin_Vis<<std::endl;



    real tau = Kin_Vis/(cs*cs);
    std::cout<<"tau     = "<<tau<<std::endl;

    real beta = 1.0/(2.0*tau + 1.0);
    std::cout<<"beta    = "<<beta<<std::endl;
   




    real Rho_mean = 1.0;
    real rho_liq = 1.63;
    real rho_gas = 0.58;


    real TbyTc = 0.96  ;       ;
    std::cout<<"T/T0    = "<<TbyTc<<std::endl;
    real kappa = 0.0;



    //fixed ------------------------------Main code--------------------------//

    initialization_equilibrium_profile_y(grid,d2q9,rho_liq, rho_gas);


    print_vtk(d2q9,grid,0.0,u0,TbyTc,kappa,Force);

    // exit(0);

    int sim_time = 50*20*Nx/u0;

    std::cout<<"Simulation time "<< sim_time<<std::endl;


    for(int t = 1; t <=10;t++){
        std::cout<<"done"<<std::endl;

        collide (grid,d2q9,beta,tau,TbyTc,kappa, t,Force,g);

        Periodic_left_Right(grid);
        // Periodic_top_bottom(grid);
      

        // BB_left     (grid,d2q9,u0);
        // BB_right    (grid,d2q9,u0);

        BB_top      (grid,d2q9,u0);
        BB_bottom   (grid,d2q9,u0);







        advection_D2Q9(grid);


        if(t%1== 0){
            std::cout<<t<<" ";
            printMass(grid);
            print_vtk(d2q9,grid,t,u0,TbyTc,kappa,Force);
        }
    }


}
   














;