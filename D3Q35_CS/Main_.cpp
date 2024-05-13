#include<iostream>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<sstream>

#include"collison.h"
#include"advection.h"
#include"print.h"
#include"boundary.h"




int main()
{


    int Nx =50;int Ny = 50; int Nz = 5;
    std::cout<<" domain size Nx =  "<<Nx<<" Ny = "<<Ny<<" Nz = "<< Nz<< std::endl;

    Grid_N_C_3D<real> grid            (Nx,Ny,Nz,2,35);
    Grid_N_C_3D<real> grid2           (Nx,Ny,Nz,2,35);
    Grid_N_C_3D<real> Force           (Nx,Ny,Nz,2,3);






    real b = 4.0;
    real a = 1.0;

	real TbyTc = 0.90;
    real T_critical = (0.377332*a)/b;
    real T_actual = TbyTc* T_critical;
    std::cout<<"T/T0        = "<<TbyTc<<std::endl;


    real rho_critical = 0.521772/b;
    real rho_by_rhoc  = 1.0;
    real Rho_mean = rho_by_rhoc*rho_critical;
    std::cout<<"Rho_mean    = "<<Rho_mean<<std::endl;




    real c = sqrt( T_actual *  (1.0/((31.0 + sqrt(7009))/252.0)));
    std::cout<<T_actual<<" "    <<(1.0/((31.0 + sqrt(7009))/252.0))<<std::endl;
    std::cout<<"c           = " <<c<<std::endl;


    lbmD3Q35<real> d3q35(c);

    //check this can c be sccaled with the same theta0 for the 15 and 35
    lbmD3Q15<real> d3q15(c);
    

    real cs = sqrt(d3q35.theta0);
    std::cout<<"theta       = "<<d3q35.theta0<<std::endl;


    real Re = 10;
    real L  = Ny;
    real Kn = 0.002;
    real Ma = Kn * Re;
    real u0 = Ma * cs;
    std ::cout<<"u0          = "<<u0<<std::endl;

    real dx = pow(3.0*b,1.0/3.0);

    real dt = dx/c;
    std::cout<<"dt          = "<<dt<<std::endl;


    real Kin_Vis = u0*(L)/Re;
    real tau = Kin_Vis/(cs*cs);
    std::cout<<"tau         = "<<tau<<std::endl;


    real tauNdim = tau /dt;
    //  real beta = 1.0/(2.0*tauNdim + 1);
    real beta = 0.5;
    std::cout<<"beta        = "<<beta<<std::endl;




    
    real kappabar   = 0.0625;
    real kappa      = kappabar*a * dx *dx;
    std::cout<<"kappa       = "<<kappa<<std::endl;





    //:fixed ------------------------------Main code--------------------------//
    
    //      initialization(grid,d3q35,Rho_mean,0.0,0.0);
    initialization_equilibrium_profile(grid,d3q35,Rho_mean);
    // initialization_2D_droplet(grid,d3q35,Rho_mean);





    print_vtk(d3q35,grid,0,u0,TbyTc,Force,dt);
    printMass(grid);
    int sim_time = 20*Nx/u0;

    std::cout<<"simulation started and Simulation time "<< sim_time<<std::endl;
    
    for(int t = 1; t <=10000;t++){

        // Periodic(grid);
        collide (grid,d3q35,d3q15,beta,tau,TbyTc,kappa, t,Force, dt ,dx,a , b);

        Periodic(grid);

        //  Diffuse_35(grid,d3q35,u0,0.0);
        //  BB_wall(grid,d3q35,u0,0.0);

        advection(grid);
        // stationary_correction(grid);

        if(t%5== 0){
            std::cout<<t<<" ";
            printMass(grid);
            // print_vtk(d3q35,grid,t,u0,TbyTc,Force,dt);
        }
    }


    // for(int t = 5001; t <=40000;t++){

    //     // Periodic(grid);
    //     collide(grid,d3q35,beta,tau,TbyTc,kappa, t);

    //     Periodic(grid);

    //     Diffuse_35(grid,d3q35,u0,0.0);
    //     //   BB_wall(grid,d3q35,u0,0.0);

    //     advection(grid);
    //     // stationary_correction(grid);


    //     if(t%500== 0){
    //         std::cout<<t<<" ";
    //         printMass(grid);
    //         print_vtk(d3q35,grid,t,u0,TbyTc);
    //     }
    // }



}
;