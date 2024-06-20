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

    int Nx =50 ;    int Ny = 100; int Nz = 2;
    std::cout<<" domain size Nx =  "<<Nx<<" Ny = "<<Ny<<" Nz = "<< Nz<< std::endl;

    Grid_N_C_3D<real> grid            (Nx,Ny,Nz,2,35);
    Grid_N_C_3D<real> grid2           (Nx,Ny,Nz,2,35);
    Grid_N_C_3D<real> Force           (Nx,Ny,Nz,2,3);

    
    
    lbmD3Q35<real> d3q35(1.0);
    lbmD3Q15<real> d3q15(1.0);


    real cs = sqrt(d3q35.theta0);
    std::cout<<"theta= "<<d3q35.theta0<<std::endl;


    real Re = 400;
    std::cout<<"Re = "<<Re<<std::endl;
    real L  = Ny;
    real Kn =0.00005;
    // real Ma = Kn * Re;
    real Ma = 0.05;
    std::cout<<"Ma = "<<Ma<<std::endl;

    real u0 = Ma * cs;
    std::cout<<"u0 = "<<u0<<std::endl;


    real Kin_Vis = u0*(L)/Re;
    // real tau = Kin_Vis/(cs*cs);
    // std::cout<<"tau "<<tau<<std::endl;

    // real beta = 1.0/(2.0*tau + 1.0);
    real beta = 0.5;
    real tau = (1.0-beta)/beta *0.5;
    std::cout<<"tau "<<tau<<std::endl;
    std ::cout<<"beta"<<beta<<std::endl;


    real Rho_mean = 1.0;
    real rho_liq = 1.53867;
    real rho_gas = 0.552003;


    real TbyTc = 0.96  ;       ;
    std::cout<<"T/T0 = "<<TbyTc<<std::endl;
    real kappa = 0.0;



    //:fixed ------------------------------Main code--------------------------//
    
    // initialization(grid,d3q35,Rho_mean,0.0,0.0);

    // initialization_equilibrium_profile_x(grid,d3q35,Rho_mean);
    initialization_equilibrium_profile_y(grid,d3q35,rho_liq, rho_gas);

    // initialization_2D_droplet(grid,d3q35,Rho_mean);







    print_vtk(d3q35,grid,0,u0,TbyTc,kappa,Force);
    printMass(grid);
    int sim_time = 20*Nx/u0;

    std::cout<<"simulation started and Simulation time "<< sim_time<<std::endl;
    
    for(int t = 1; t <=50000;t++){

        // Periodic(grid);
        collide (grid,d3q35,d3q15,beta,tau,TbyTc,kappa, t,Force);

        // Periodic_x(grid);
        // Periodic_y(grid);
        Periodic_z(grid); //#fixeed for all

        // Diffuse_35(grid,d3q35,u0,0.0);
        BB_wall_top     (grid,d3q35,u0,0.0);
        BB_wall_bottom  (grid,d3q35,u0,0.0);    
        BB_wall_left    (grid,d3q35,u0,0.0);
        BB_wall_right   (grid,d3q35,u0,0.0);

        advection(grid);
        // stationary_correction(grid);

        if(t%500== 0){
            std::cout<<t<<" ";
            printMass(grid);
            print_vtk(d3q35,grid,t,u0,TbyTc,kappa,Force);
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