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


    Grid_N_C_2D<double> gridf               (Nx,Ny,1,9);
    Grid_N_C_2D<double> gridg               (Nx,Ny,1,9);
    Grid_N_C_2D<double> grad_psi_rho        (Nx,Ny,1,2); // $psi 
    Grid_N_C_2D<double> Force               (Nx,Ny,1,2); // $psi 



    lbmD2Q9<double> d2q9(1.0,(1.0/3.0));
    
    double cs = sqrt(d2q9.theta0);
    std::cout<<"theta= "<<d2q9.theta0<<std::endl;




    double Re = 2048;
    double L = Nx;

    double u0 = 0.04;

    double g =(u0*u0)/L;
    std::cout<<"g       = "<<g<<std::endl;


    double Kin_Vis = u0*(L)/Re;

    std::cout<<"Kin_Vis = "<<Kin_Vis<<std::endl;


    double tau = Kin_Vis/(cs*cs) + 0.5;
    std::cout<<"tau     = "<<tau<<std::endl;



    double beta = 1.0/(2.0*tau + 1.0);
    std::cout<<"beta    ="<<beta<<std::endl;
    

    double rho_l = 0.50, rho_h = 1.61;
    double phi_l = -1.0, phi_h = 1.0;



    double TbyTc = 0.95;
    double rho_critical = 1.0, T_critical = d2q9.theta0/TbyTc ; 
    double b = 0.521772/(rho_critical), a = b*T_critical/0.377332;

    double kappa = 0;


    //fixed ------------------------------Main code--------------------------//
    initialization(gridf,gridg,d2q9,phi_l,phi_h,rho_l, rho_h, a,b   );
    print_vtk(d2q9,gridf,gridg,grad_psi_rho,0.0,u0, kappa,phi_l,phi_h,rho_l, rho_h,Force);

    // exit(0);
    printMass(gridf,gridg);

    int sim_time = 50*20*Nx/u0;

    std::cout<<"Simulation time "<< sim_time<<std::endl;


    for(int t = 1; t <=50000;t++){


        collide(gridf,gridg,grad_psi_rho,d2q9,beta,tau,kappa,g, phi_l, phi_h,rho_l, rho_h,  a, b, Force);

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


        if(t%25== 0){
            std::cout<<t<<" ";
            printMass(gridf,gridg);
            print_vtk(d2q9,gridf,gridg,grad_psi_rho,t,u0, kappa,phi_l,phi_h,rho_l, rho_h, Force );
        }
    }


}
   














;