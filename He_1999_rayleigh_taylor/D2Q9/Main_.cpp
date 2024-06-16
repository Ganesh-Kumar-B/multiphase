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

    int Nx = 125 ;int Ny = 125;
    Grid_N_C_2D<double> gridf(Nx,Ny,1,9);
    Grid_N_C_2D<double> gridg(Nx,Ny,1,9);


    double kappa = 0.5;
    double theta = 0.69;

    lbmD2Q9<double> d2q9(1.0,(1.0/3.0));
    // double beta = 0.53;
    // double tau =0.5*(1/beta-1);
    // double Re=100.0;
    double cs = sqrt(d2q9.theta0);
    std::cout<<"theta= "<<d2q9.theta0<<std::endl;
    // std::cout<<cs<<std::endl;
    // double u0 = 0.05   *cs;
    // double gx = 8.0*u0*u0/(40*Re);
    // double gy = 0.0*u0*u0/(40*Re);



    double Re = 400;
    double L = Ny;
    double Kn =0.0005;
    double Ma = 0.05;
    double u0 = Ma * cs;
    std::cout<<"u0 = "<<u0<<std::endl;



    double Kin_Vis = u0*(L)/Re;
    double tau = Kin_Vis/(cs*cs);
    std::cout<<"tau "<<tau<<std::endl;

    double beta = 1.0/(2.0*tau + 1.0);
    std::cout<<"beta"<<beta<<std::endl;
    double gx = 0.0*u0*u0/(L*Re);
    std::cout<<"gx = "<<gx<<std::endl;
    double gy = 0.0*u0*u0/(L*Re);

    //fixed ------------------------------Main code--------------------------//
    initialization(gridf,gridg,d2q9,0.0,u0);
    print_vtk(d2q9,gridf,gridg,0.0,u0, kappa , theta);

    printMass(gridf,gridg);

int sim_time = 50*20*Nx/u0;

std::cout<<"Simulation time "<< sim_time<<std::endl;


    for(int t = 1; t <=500000;t++){

        collide(gridf,gridg,d2q9,beta,tau,gx,gy);
        
        // Periodic_left_Right(gridf,gridg);
        // Periodic_top_bottom(gridf,gridg);


        BB_top(gridf,gridg,d2q9,u0);
        BB_bottom(gridf,gridg,d2q9,u0);
        BB_left(gridf,gridg,d2q9,u0);
        BB_right(gridf,gridg,d2q9,u0);



        advection_D2Q9(gridf,gridg);

        

        if(t%5000== 0){
            std::cout<<t<<" ";
            printMass(gridf,gridg);
            print_vtk(d2q9,gridf,gridg,0.0,u0, kappa , theta);
        }
    }


}
   














;