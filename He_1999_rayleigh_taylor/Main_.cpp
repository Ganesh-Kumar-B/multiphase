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

    int Nx =100;int Ny = 100; int Nz = 5;
    std::cout<<" domain size Nx =  "<<Nx<<" Ny = "<<Ny<<" Nz = "<< Nz<< std::endl;

    Grid_N_C_3D<real> gridf                             (Nx,Ny,Nz,2,35);
    Grid_N_C_3D<real> gridg                             (Nx,Ny,Nz,2,35);
    Grid_N_C_3D<real> Force                             (Nx,Ny,Nz,2,3);

    Grid_N_C_3D<real>  rho                              (Nx,Ny,Nz,2,1);
    Grid_N_C_3D<real>  phi                              (Nx,Ny,Nz,2,1);   
    Grid_N_C_3D<real>  mu                               (Nx,Ny,Nz,2,1);   
    Grid_N_C_3D<real>  laplacian_phi                    (Nx,Ny,Nz,2,1);   
    
    lbmD3Q35<real> d3q35(1.0,0.33333333333333);

    real cs = sqrt(d3q35.theta0);
    std::cout<<"theta= "<<d3q35.theta0<<std::endl;


    real Re = 10;
    real L  = Ny;
    real Kn =0.002;
    real Ma = Kn * Re;
    real u0 = Ma * cs;
    std::cout<<"u0 = "<<u0<<std::endl;



    real Kin_Vis = u0*(L)/Re;


    real tau = Kin_Vis/(cs*cs);

    real tauphi = 1.0; 
    std::cout<<"tau "<<tau<<std::endl;


    real beta = 1.0/(2.0*tau + 1);
    std::cout<<"beta"<<beta<<std::endl;


    real Rho_mean = 1.0;



    real TbyTc = 0.82      ;
    std::cout<<"T/T0 = "<<TbyTc<<std::endl;

    real kappa = 0.00318246;
    real gamma_s = 0.1 ;
    real A = 0.003535;


    


    //:fixed ------------------------------Main code--------------------------//
    
    //      initialization(gridf,d3q35,Rho_mean,0.0,0.0);

    //  initialization_equilibrium_profile(gridf,d3q35,Rho_mean);

    initialization_2D_droplet(gridf,gridg,rho,phi,mu,laplacian_phi,d3q35,Rho_mean,kappa, gamma_s, A);

    print_vtk(d3q35,gridf,gridg,0,u0,TbyTc,Force);
    printMass(gridf);
    int sim_time = 20*Nx/u0;

    std::cout<<"simulation started and Simulation time "<< sim_time<<std::endl;
    
    for(int t = 1; t <=10000;t++){

        // Periodic(gridf);
        collide (gridf,gridg,rho,phi,mu,laplacian_phi,d3q35,beta,tau,tauphi,TbyTc, t,Force, kappa ,gamma_s,A);

        Periodic(gridf);
        Periodic(gridg);

        //  Diffuse_35(gridf,d3q35,u0,0.0);
        //  BB_wall(gridf,d3q35,u0,0.0);

        advection(gridf);
        advection(gridg);
        // stationary_correction(gridf);

        if(t%1== 0){
            std::cout<<t<<" ";
            printMass(gridf);
            print_vtk(d3q35,gridf,gridg,t,u0,TbyTc,Force);
        }
    }


    // for(int t = 5001; t <=40000;t++){

    //     Periodic(gridf);
    //     collide(gridf,d3q35,beta,tau,TbyTc,kappa, t);

    //     Periodic(gridf);

    //     Diffuse_35(gridf,d3q35,u0,0.0);
    //     //   BB_wall(gridf,d3q35,u0,0.0);

    //     advection(gridf);
    //     // stationary_correction(gridf);


    //     if(t%500== 0){
    //         std::cout<<t<<" ";
    //         printMass(gridf);
    //         print_vtk(d3q35,gridf,t,u0,TbyTc);
    //     }
    // }



}
;