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

    int Nx = 256 ;int Ny = 256;


    Grid_N_C_2D<real> grid                  (Nx,Ny,1,9);
    Grid_N_C_2D<real> Force                 (Nx,Ny,1,2);
    Grid_N_C_2D<real> P_tensor              (Nx,Ny,1,4);   //4 components - 0-xx, 1xy, 2yx, 3yy


    real c =1.0;
    lbmD2Q9<real> d2q9(c);
    
    real cs = sqrt(d2q9.theta0);
    std::cout<<"theta   = "<<d2q9.theta0<<std::endl;


    real dx = 1.0/Nx;
    real dt = dx/c;



    real Re = 1000;
    std::cout<<"Re      = "<<Re<<std::endl;
    real L  = Nx;


    real u0 = 0.001;
    std::cout<<"u0      = "<<u0<<std::endl;


    real Ma = u0/cs;
    std::cout<<"Ma      = "<<Ma<<std::endl;


    // real g =1.0*(u0*u0)/1.0;
    real g = 0.0;
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
    real rho_liq =  1.5819;
    real rho_gas =  0.46344;





    real TbyTc = 0.95  ;       ;
    std::cout<<"T/T0    = "<<TbyTc<<std::endl;
    real kappa = 0.006;
    std::cout<<"kappa   = "<<kappa<<std::endl;



    //fixed ------------------------------Main code--------------------------//
    // initialization(grid,d2q9,Rho_mean, rho_liq, rho_gas);
    // initialization_equilibrium_profile_y(grid,d2q9,rho_liq, rho_gas);

    // initialization_ellipse(grid,d2q9,rho_liq, rho_gas); // with bubble on top domain

    initialization_circle(grid,d2q9,rho_liq, rho_gas);


    std::string name="Result_k_06";
    print_vtk(d2q9,grid,0.0,u0,TbyTc,kappa,Force,P_tensor,name, dx ,dt);



    // exit(0);


    int sim_time = 2*Nx/u0;

    std::cout<<"Simulation time "<< sim_time<<std::endl;


    for(int t = 1; t <=200000;t++){

        collide (grid,d2q9,beta,tau,TbyTc,kappa, t,Force,P_tensor,g, dx, dt);

        Periodic_left_Right(grid);
        Periodic_top_bottom(grid);


        // BB_left     (grid,d2q9,u0);
        // BB_right    (grid,d2q9,u0);

        // BB_top      (grid,d2q9,u0);
        // BB_bottom   (grid,d2q9,u0);



        advection_D2Q9(grid);


        if(t%2000== 0){
            std::cout<<t<<" ";
            printMass(grid);
            print_vtk(d2q9,grid,t,u0,TbyTc,kappa,Force,P_tensor,name,dx,dt);
        }
    }


}









    //for the drop acoustic thing
    // real Re = 2048;
    // std::cout<<"Re      = "<<Re<<std::endl;
    // real L  = Nx;


    // real u0 = 0.01;
    // std::cout<<"u0      = "<<u0<<std::endl;

    // real g = 10.0*(u0*u0)/L;
    // // real g = 0.0;
    // std::cout<<"g       = "<<g<<std::endl;



    // real Kin_Vis = u0*(L)/Re;
    // std::cout<<"Kin_Vis = "<<Kin_Vis<<std::endl;


    // real tau = Kin_Vis/(cs*cs);
    // // real tau = Kin_Vis/(cs*cs) +0.5;
    // std::cout<<"tau     = "<<tau<<std::endl;

    // real beta = 1.0/(2.0*tau + 1.0);
    // std::cout<<"beta    = "<<beta<<std::endl;





    // real Rho_mean = 1.0;
    // real rho_liq =  1.61;
    // real rho_gas =0.50;


    // real TbyTc = 0.95  ;       ;
    // std::cout<<"T/T0    = "<<TbyTc<<std::endl;
    // real kappa = 0.00;







;