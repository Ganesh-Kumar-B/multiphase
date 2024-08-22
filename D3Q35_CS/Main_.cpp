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

    int Nx = 256 ;int Ny = 256; int Nz = 2;
    std::cout<<" domain size Nx =  "<<Nx<<" Ny = "<<Ny<<" Nz = "<< Nz<< std::endl;

    Grid_N_C_3D<real> grid            (Nx,Ny,Nz,2,35);
    Grid_N_C_3D<real> grid2           (Nx,Ny,Nz,2,35);
    Grid_N_C_3D<real> Force           (Nx,Ny,Nz,2,3);

    real c = 1.0;

    
    lbmD3Q35<real> d3q35(c);
    lbmD3Q15<real> d3q15(c);


    real cs = sqrt(d3q35.theta0);
    std::cout<<"theta   = "<<d3q35.theta0<<std::endl;


    
    real dx = 64.0/Nx;
    std::cout<<"dx      = "<<dx<<std::endl;
    real dt = dx/c;


    real Re = 100;
    std::cout<<"Re      = "<<Re<<std::endl;
    real L  = 128;


    real u0 = 0.04;
    std::cout<<"u0      = "<<u0<<std::endl;

    
    real Ma = u0/cs;
    std::cout<<"Ma      = "<<Ma<<std::endl;


    real g = (u0*u0)/L;
    g = 0;
    std::cout<<"g       = "<<g<<std::endl;


    real Kin_Vis = u0*(L)/Re;
    std::cout<<"Kin_Vis = "<<Kin_Vis<<std::endl;

    real tau = Kin_Vis/(cs*cs);
    real tauNdim = tau/dt;


    // real tau = Kin_Vis/(cs*cs) +0.5;
    std::cout<<"tau     = "<<tau<<std::endl;


    real beta = 1.0/(2.0*tauNdim + 1.0);
    std::cout<<"beta    = "<<beta<<std::endl;





    real Rho_mean = 1.0;
    real rho_liq = 1.60163;
    real rho_gas = 0.534782;


    real TbyTc = 0.954  ;       ;
    std::cout<<"T/T0    = "<<TbyTc<<std::endl;
    real kappa = 0.01*dx*dx;
    std::cout<<"kappa   = "<<kappa<<std::endl;
    real sigma = 0.0;

    
    real rho_critical = 1.0, T_critical = d3q35.theta0/TbyTc ; 
    real b = 0.521772/(rho_critical), a = b*T_critical/0.377332;



    //:fixed ------------------------------Main code--------------------------//
    
    // initialization(grid,d3q35,Rho_mean,0.0,0.0);

    // initialization_equilibrium_profile_x(grid,d3q35,Rho_mean);
    // initialization_equilibrium_profile_y(grid,d3q35,rho_liq, rho_gas);


    real R = 0.20;
    std::cout<<"circle Radius: "<<R*Nx<<std::endl;
    initialization_2D_droplet(grid,d3q35,Rho_mean,rho_liq,rho_gas, R);





    std::string name="Result_256_0.20_0.01";

    print_vtk(d3q35,grid,0,u0,TbyTc,kappa, a, b,Force,name, dx, dt);
    printMass(grid);
    int sim_time = 20*Nx/u0;

    std::cout<<"simulation started and Simulation time "<< sim_time<<std::endl;
    
    for(int t = 1; t <=25000;t++){

        // Periodic(grid);
        collide (grid,d3q35,d3q15,beta,tau,TbyTc,kappa, a, b,sigma, t,Force,g, dx, dt);

        Periodic_x(grid);
        Periodic_y(grid);
        Periodic_z(grid); //#fixeed for all

        // Diffuse_35(grid,d3q35,u0,0.0);
        // BB_wall_top     (grid,d3q35,u0,0.0);
        // BB_wall_bottom  (grid,d3q35,u0,0.0);    
        // BB_wall_left    (grid,d3q35,u0,0.0);
        // BB_wall_right   (grid,d3q35,u0,0.0);

        advection(grid);
        // stationary_correction(grid);

        if(t%1000== 0){
            std::cout<<sigma<<std::endl;
            std::cout<<t<<" ";
            printMass(grid);
            print_vtk(d3q35,grid,t,u0,TbyTc,kappa, a, b,Force,name, dx, dt);
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