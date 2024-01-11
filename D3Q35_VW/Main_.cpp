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

   int Nx =30;int Ny = 30; int Nz = 30;
   std::cout<<"domain size Nx =  "<<Nx<<" Ny = "<<Ny<<" Nz = "<< Nz<< std::endl;

   Grid_N_C_3D<double> grid            (Nx,Ny,Nz,2,35);
   Grid_N_C_3D<double> pnid            (Nx,Ny,Nz,2,1);
   Grid_N_C_3D<double> munid           (Nx,Ny,Nz,2,1);
   Grid_N_C_3D<double> fnid            (Nx,Ny,Nz,2,1);
   Grid_N_C_3D<double> rho             (Nx,Ny,Nz,2,1);
   Grid_N_C_3D<double> laplacian_rho   (Nx,Ny,Nz,2,1);

   
   lbmD3Q35<double> d3q35(1.0,0.33333333333333);


   double cs = sqrt(d3q35.theta0);
   std::cout<<"theta= "<<d3q35.theta0<<std::endl;


   double Re = 1;
   double L = Ny;
   double Kn =0.002;
   double Ma = Kn * Re;
   double u0 = Ma * cs;
   std::cout<<"u0 = "<<u0<<std::endl;



   double Kin_Vis = u0*(L)/Re;
   double tau = Kin_Vis/(cs*cs);
   std::cout<<"tau "<<tau<<std::endl;


   double beta = 1.0/(2.0*tau + 1);
   std::cout<<"beta"<<beta<<std::endl;


   double Rho_mean = 1.0;

   double TbyTc = 0.90;
   std::cout<<"T/T0 = "<<TbyTc<<std::endl;
   double kappa = 0.0625;



   //:fixed ------------------------------Main code--------------------------//
   double check = 0;
   for(int dv = 0; dv<grid.d_v; dv++){
      check += d3q35.W[dv]; 
   }
   std::cout<<check<<std::endl;
   initialization(grid,d3q35,Rho_mean,0.001,2);
   printdata(d3q35,grid,0,u0);
   printMass(grid);

   int sim_time = 20*Nx/u0;

   std::cout<<"Simulation time "<< sim_time<<std::endl;
 

   for(int t = 1; t <=25000;t++){

      // Periodic(grid);
      collide(grid,rho,pnid,fnid,munid,laplacian_rho,d3q35,beta,tau,TbyTc,kappa, t);

      Periodic(grid);

      advection(grid);

      if(t%1000== 0){
         std::cout<<t<<" ";
         printMass(grid);
         printdata(d3q35,grid,t,u0);
      }
   }

    
}
   


;