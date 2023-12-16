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

   int Nx =100;int Ny = 100;

   Grid_N_C_2D<double> grid            (Nx,Ny,1,9);
   Grid_N_C_2D<double> pnid            (Nx,Ny,1,1);
   Grid_N_C_2D<double> munid           (Nx,Ny,1,1);
   Grid_N_C_2D<double> fnid            (Nx,Ny,1,1);
   Grid_N_C_2D<double> rho             (Nx,Ny,1,1);
   Grid_N_C_2D<double> laplacian_rho   (Nx,Ny,1,1);

   
   lbmD2Q9<double> d2q9(1.0,0.33333333333333);


   double cs = sqrt(d2q9.theta0);
   std::cout<<"theta= "<<d2q9.theta0<<std::endl;


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

   double TbyTc = 0.88;
   std::cout<<"T/T0 = "<<TbyTc<<std::endl;
   double kappa = 0.0625;


   //fixed ------------------------------Main code--------------------------//
   initialization(grid,d2q9,Rho_mean,0.001,2);
   printdata(d2q9,grid,0,u0);
   printMass(grid);

   int sim_time = 20*Nx/u0;

   std::cout<<"Simulation time "<< sim_time<<std::endl;
 

   for(int t = 1; t <=50000;t++){

      Periodic(grid);
      
      collide(grid,rho,pnid,fnid,munid,laplacian_rho,d2q9,beta,tau,TbyTc,kappa, t);

      Periodic(grid);

      advection_D2Q9(grid);
      if(t%1000== 0){
         std::cout<<t<<" ";
         printMass(grid);
         printdata(d2q9,grid,t,u0);
      }
   }

    
}
   


;