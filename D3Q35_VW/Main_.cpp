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

   int Nx =200;int Ny = 200; int Nz = 5;
   std::cout<<" domain size Nx =  "<<Nx<<" Ny = "<<Ny<<" Nz = "<< Nz<< std::endl;

   Grid_N_C_3D<real> grid            (Nx,Ny,Nz,2,35);
  
  
   
   lbmD3Q35<real> d3q35(1.0,0.33333333333333);

    
   real cs = sqrt(d3q35.theta0);
   std::cout<<"theta= "<<d3q35.theta0<<std::endl;

   real Re = 1;
   real L = Ny;
   real Kn =0.002;
   real Ma = Kn * Re;
   real u0 = Ma * cs;
   std::cout<<"u0 = "<<u0<<std::endl;



   real Kin_Vis = u0*(L)/Re;
   real tau = Kin_Vis/(cs*cs);
   std::cout<<"tau "<<tau<<std::endl;


   real beta = 1.0/(2.0*tau + 1);
   std::cout<<"beta"<<beta<<std::endl;


   real Rho_mean = 1.0;



   real TbyTc = 0.92;
   std::cout<<"T/T0 = "<<TbyTc<<std::endl;
   real kappa = 0.000125;





   //:fixed ------------------------------Main code--------------------------//
   

   // initialization_equilibrium_profile(grid,d3q35,Rho_mean);
   initialization_2D_droplet(grid,d3q35,Rho_mean);
   // exit(0);

   print_vtk(d3q35,grid,0,u0,TbyTc);
   printMass(grid);
   int sim_time = 20*Nx/u0;

   std::cout<<"simulation started and Simulation time "<< sim_time<<std::endl;
   
   for(int t = 1; t <=10000;t++){

      // Periodic(grid);
      collide(grid,d3q35,beta,tau,TbyTc,kappa, t);

      Periodic(grid);

      advection(grid);

      if(t%500== 0){
         std::cout<<t<<" ";
         printMass(grid);
         print_vtk(d3q35,grid,t,u0,TbyTc);
      }
   }
}
;