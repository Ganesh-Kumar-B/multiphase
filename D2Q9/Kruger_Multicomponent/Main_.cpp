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

    int Nx = 128 ;int Ny = 512;

    Grid_N_C_2D<real> gridf                             (Nx,Ny,1,9);
    Grid_N_C_2D<real> gridg                             (Nx,Ny,1,9);

    Grid_N_C_2D<real> Force                             (Nx,Ny,1,2);
    Grid_N_C_2D<real>  rho                              (Nx,Ny,1,1);
    Grid_N_C_2D<real>  phi                              (Nx,Ny,1,1);   
    Grid_N_C_2D<real>  mu                               (Nx,Ny,1,1);   
    Grid_N_C_2D<real>  laplacian_phi                    (Nx,Ny,1,1);   
    
    lbmD2Q9<real> d2q9(1.0,(1.0/3.0));

    real cs = sqrt(d2q9.theta0);
    std::cout<<"theta= "<<d2q9.theta0<<std::endl;


    real Re = 1000;
    real L  = Nx;
    

    real u0 =0.04;
    std::cout<<"u0 = "<<u0<<std::endl;



    real Kin_Vis = u0*(L)/Re;


    real tau = Kin_Vis/(cs*cs) +0.5;

    real tauphi = 1.0; 
    std::cout<<"tau "<<tau<<std::endl;


    real beta = 1.0/(2.0*tau + 1);
    std::cout<<"beta"<<beta<<std::endl;


    real Rho_mean = 1.0;




    real kappa = 0.00318246;
    real gamma_s = 1 ;
    real A = 0.003535;

    real g = u0*u0/L;

    //fixed ------------------------------Main code--------------------------//
    // initialization(grid,d2q9,Rho_mean);
    // initialization_2D_droplet(gridf,gridg,rho,phi,mu,laplacian_phi,d2q9,Rho_mean,kappa, gamma_s, A);

    initialization_y_RT(gridf,gridg,rho,phi,mu,laplacian_phi,d2q9,Rho_mean,kappa, gamma_s, A);



    std::string name="Result";
    print_vtk(d2q9,gridf,gridg,0,u0,Force);


    printMass(gridf);

    // exit(0);


    int sim_time = 50*20*Nx/u0;

    std::cout<<"Simulation time "<< sim_time<<std::endl;


    for(int t = 1; t <=50000;t++){

        collide (gridf,gridg,rho,phi,mu,laplacian_phi,d2q9,beta,tau,tauphi, t,Force, kappa ,gamma_s,A,g);

        Periodic_left_Right(gridf);
        // Periodic_top_bottom(gridf);
      
        Periodic_left_Right(gridg);
        // Periodic_top_bottom(gridg);


        // BB_left     (gridf,d2q9,u0);
        // BB_right    (gridf,d2q9,u0);
        BB_top      (gridf,d2q9,u0);
        BB_bottom   (gridf,d2q9,u0);

        // BB_left     (gridg,d2q9,u0);
        // BB_right    (gridg,d2q9,u0);
        BB_top      (gridg,d2q9,u0);
        BB_bottom   (gridg,d2q9,u0);


        advection_D2Q9(gridf);
        advection_D2Q9(gridg);


        if(t%500== 0){
            std::cout<<t<<" ";
            printMass(gridf);
            print_vtk(d2q9,gridf,gridg,t,u0,Force);
        }
    }


}
   






// working for this parameters

// real kappa = 0.00318246;
//     real gamma_s = 0.1 ;
//     real A = 0.003535;

//     real g = 0.000025;





;





















    // real Re = 10;
    // real L  = Ny;
    // real Kn =0.002;
    // real Ma = Kn * Re;
    // real u0 = Ma * cs;
    // std::cout<<"u0 = "<<u0<<std::endl;



    // real Kin_Vis = u0*(L)/Re;


    // real tau = Kin_Vis/(cs*cs);

    // real tauphi = 1.0; 
    // std::cout<<"tau "<<tau<<std::endl;


    // real beta = 1.0/(2.0*tau + 1);
    // std::cout<<"beta"<<beta<<std::endl;


    // real Rho_mean = 1.0;



    // real TbyTc = 0.82      ;
    // std::cout<<"T/T0 = "<<TbyTc<<std::endl;

    // real kappa = 0.00318246;
    // real gamma_s = 1 ;
    // real A = 0.003535;

    // real g = 0.000025;