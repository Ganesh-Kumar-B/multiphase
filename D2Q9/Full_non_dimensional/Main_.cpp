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

    int Nx = 128 ;int Ny = 128;

    std::cout<<"Nx = "<<Nx<<"Ny = "<<Ny<<std::endl;
    Grid_N_C_2D<real> grid                  (Nx,Ny,1,9);
    Grid_N_C_2D<real> Force                 (Nx,Ny,1,2);
    Grid_N_C_2D<real> P_tensor              (Nx,Ny,1,4);   //4 components - 0-xx, 1xy, 2yx, 3yy



    real c = 1.0;
    lbmD2Q9<real> d2q9(c);
    
    real cs = sqrt(d2q9.theta0);
    std::cout<<"theta   = "<<d2q9.theta0<<std::endl;



    real L  = 128;
    real dx = L/Nx ;
    real dt = dx/c;
    std::cout<<"dx      = "<<dx<<std::endl;


    real Re = 100;
    std::cout<<"Re      = "<<Re<<std::endl;


    real u0 = 0.05;
    std::cout<<"u0      = "<<u0<<std::endl;


    real Ma = u0/cs;
    std::cout<<"Ma      = "<<Ma<<std::endl;


    real g = 0.0;
    std::cout<<"g       = "<<g<<std::endl;



    real Kin_Vis = u0*(L)/Re;
    std::cout<<"Kin_Vis = "<<Kin_Vis<<std::endl;



    real tau = Kin_Vis/(cs*cs);
    real tauNdim = tau/dt;

    std::cout<<"tau     = "<<tau<<std::endl;



    real beta = 1.0/(2.0*tauNdim + 1.0);
    std::cout<<"beta    = "<<beta<<std::endl;





    real rho_liq =  1.6165;
    real rho_gas =  0.49947;
    std::cout<<"rhol    = "<<rho_liq<<std::endl;
    std::cout<<"rhog    = "<<rho_gas<<std::endl;



    real theta_r = 0.95  ;       ;
    std::cout<<"T/T0    = "<<theta_r<<std::endl;



    real rho_critical = 1.0, T_critical = d2q9.theta0/theta_r ; 
    std::cout<<"rho_c    = "<<rho_critical<<std::endl;

    // real b = 0.521772/(rho_critical), a = b*T_critical/0.377332;    //CS

    // real b = 1.0/(3.0*rho_critical), a = b*T_critical*27.0/8.0;   //VW

    real kappa = 1.0;
    std::cout<<"kappa   = "<<kappa<<std::endl;
    real sigma  = 0;






    // ------------------------------Main code--------------------------//
    

    real R = 0.2;  // --> for the finite  
    std::cout<<"circle Radius: "<<R*L<<std::endl;
    initialization_circle(grid,d2q9,rho_liq, rho_gas,rho_critical,R);


    std::string name="Results_0.2";
    // print_vtk(d2q9,grid,0.0,u0,theta_r,kappa, a, b,Force,P_tensor,name, dx ,dt);



    int sim_time = 2*Nx/u0;

    std::cout<<"Simulation time "<< sim_time<<std::endl;


    for(int t = 1; t <=50000;t++){


        collide (grid,d2q9,beta,tau,theta_r,kappa,sigma, t,Force,P_tensor,g, dx, dt);



        Periodic_left_Right(grid);
        Periodic_top_bottom(grid);


        advection_D2Q9(grid);


        if(t%1== 0 ){
            std::cout<<t<<"     =";
            printMass(grid);
            // print_vtk(d2q9,grid,t,u0,TbyTc,kappa, a,b,Force,P_tensor,name,dx,dt);
        }
    }


}




;