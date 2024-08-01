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

    std::cout<<"Nx = "<<Nx<<"Ny = "<<Ny<<std::endl;
    Grid_N_C_2D<real> grid                  (Nx,Ny,1,9);
    Grid_N_C_2D<real> Force                 (Nx,Ny,1,2);
    Grid_N_C_2D<real> P_tensor              (Nx,Ny,1,4);   //4 components - 0-xx, 1xy, 2yx, 3yy



    real c = 1.0;
    lbmD2Q9<real> d2q9(c);
    
    real cs = sqrt(d2q9.theta0);
    std::cout<<"theta   = "<<d2q9.theta0<<std::endl;


    real dx = 0.5;
    real dt = dx/c;
    std::cout<<"dx      = "<<dx<<std::endl;


    real Re = 200;
    std::cout<<"Re      = "<<Re<<std::endl;
    real L  = 128;


    real u0 = 0.05;
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


    // real Rho_mean = 1.0;
    // real rho_liq =  1.5819;
    // real rho_gas =  0.46344;



    real rho_liq =  1.36861;
    real rho_gas =  0.688708;



    real TbyTc = 0.98  ;       ;
    std::cout<<"T/T0    = "<<TbyTc<<std::endl;
    real kappa = 0.001;
    std::cout<<"kappabar   = "<<kappa<<std::endl;
    real sigma  = 0;



    real rho_critical = 1.0, T_critical = d2q9.theta0/TbyTc ; 

    real b = 0.521772/(rho_critical), a = b*T_critical/0.377332;    //CS

    // double b = 1.0/(3.0*rho_critical), a = b*T_critical*27.0/8.0;   //VW

    std::cout<<"kappa   ="<<kappa<<std::endl;





    //fixed ------------------------------Main code--------------------------//
    //  initialization(grid,d2q9,Rho_mean, rho_liq, rho_gas);
    //  initialization_equilibrium_profile_y(grid,d2q9,rho_liq, rho_gas);

    //  initialization_ellipse(grid,d2q9,rho_liq, rho_gas); // with bubble on top domain

    real R = 0.20;
    std::cout<<"circle Radius: "<<R*Nx<<std::endl;
    initialization_circle(grid,d2q9,rho_liq, rho_gas,R);


    std::string name="Result_256_0.20_yes";
    print_vtk(d2q9,grid,0.0,u0,TbyTc,kappa,Force,P_tensor,name, dx ,dt);



    // exit(0);


    int sim_time = 2*Nx/u0;

    std::cout<<"Simulation time "<< sim_time<<std::endl;


    for(int t = 1; t <=50000;t++){


        collide (grid,d2q9,beta,tau,TbyTc, a, b,kappa,sigma, t,Force,P_tensor,g, dx, dt);



        Periodic_left_Right(grid);
        Periodic_top_bottom(grid);


        // BB_left     (grid,d2q9,u0);
        // BB_right    (grid,d2q9,u0);

        // BB_top      (grid,d2q9,u0);
        // BB_bottom   (grid,d2q9,u0);



        advection_D2Q9(grid);


        if(t%500== 0 ){
            std::cout<<t<<"     =";
            std::cout<<"sigma "<<sigma<<std::endl;
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