#pragma once

#include<iostream>
#include<cmath>
#include<vector>
#include<math.h>
#include<fstream>
#include<algorithm>
#include <sstream>
#include<string>
#include "lbmD3Q35.h"
#include "GRID_3D.h"
#include "multiphase.h"
#define PI 3.14159265



template<typename T, typename T1>
void collide(Grid_N_C_3D<T> &gridf,Grid_N_C_3D<T> &gridg,Grid_N_C_3D<T> &rho, Grid_N_C_3D<T> &phi,Grid_N_C_3D<T> &mu,Grid_N_C_3D<T> &laplacian_phi,
            lbmD3Q35<T1> &lb,real beta,real tau,real tauphi, real TbyTc, int t,Grid_N_C_3D<T> &Force, real kappa,real gamma_s, real A ){

    

    real feq_Node[35] = {0}, feq_Cell[35] = {0}, ux = 0, uy = 0, uz = 0;

    real eta =0;   //   0 ----> fourth order   1-----> second order 
    

    Multiphase_terms(gridf, gridg,lb,rho,phi , mu,laplacian_phi,kappa, A );

    //< Collision
    for(int i = 0 + gridf.noghost; i < gridf.n_x_node - (gridf.noghost) ; i++){
        for(int j = 0 + gridf.noghost;j < gridf.n_y_node - (gridf.noghost) ; j++){
            for(int k = 0 + gridf.noghost;k < gridf.n_z_node - (gridf.noghost) ; k++){
                
                real Rho = 0.0;

                Multiphase_Force_Node(gridf,phi , mu,lb,Force,i,j,k );             										
                
                get_moments_Node_f(gridf, lb,  ux, uy, uz,Rho, i, j, k,Force);            //for the node
                get_equi_f(feq_Node,lb, ux, uy,uz, Rho, phi.Node(i,j,k),mu.Node(i,j,k));

                // // //> normal
                for (int dv = 0; dv< gridf.d_v; dv++){
                    gridf.Node(i,j,k,dv) =  gridf.Node(i,j,k,dv) + (1.0/tau )* beta*(feq_Node[dv] - gridf.Node(i,j,k,dv))
                                        +(1.0 - (1.0 / (2.0* tau)  ) )*lb.thetaInverse * rho.Node(i,j,k)* lb.W[dv] * (Force.Node(i,j,k,0) * lb.Cx[dv] + Force.Node(i,j,k,1) * lb.Cy[dv] + Force.Node(i,j,k,1) * lb.Cz[dv])
                                        ;
                }

                //< CELLS    
                Rho = 0.0;

                Multiphase_Force_Cell(gridf,rho, mu,lb,Force,i,j,k );                                               

                get_moments_Cell_f(gridf, lb,  ux, uy, uz,Rho, i, j, k,Force);            //for the node
                get_equi_f(feq_Cell,lb, ux, uy,uz, Rho,phi.Cell(i,j,k),mu.Node(i,j,k));

                //> normal
                for (int dv = 0; dv< gridf.d_v; dv++){
                    gridf.Cell(i,j,k,dv) =  gridf.Cell(i,j,k,dv) +(1.0/tau )*(feq_Cell[dv] - gridf.Cell(i,j,k,dv))
                                        + (1.0 - (1.0 / (2.0* tau)  ) )*lb.thetaInverse * rho.Cell(i,j,k)* lb.W[dv] * (Force.Cell(i,j,k,0) * lb.Cx[dv] + Force.Cell(i,j,k,1) * lb.Cy[dv] + Force.Cell(i,j,k,1) * lb.Cz[dv])
                                        ;
                }       
            }
        }
    }


    for(int i = 0 + gridg.noghost; i < gridg.n_x_node - (gridg.noghost) ; i++){
        for(int j = 0 + gridg.noghost;j < gridg.n_y_node - (gridg.noghost) ; j++){
            for(int k = 0 + gridg.noghost;k < gridg.n_z_node - (gridg.noghost) ; k++){
                
                real Rho = 0.0;
                
                get_moments_Node_f(gridf, lb,  ux, uy, uz,Rho, i, j, k,Force);  

                get_equi_g(feq_Node,lb, ux, uy,uz, phi.Node(i,j,k),gamma_s,mu.Node(i,j,k));

                // // //> normal
                for (int dv = 0; dv< gridg.d_v; dv++){
                    gridg.Node(i,j,k,dv) =  gridg.Node(i,j,k,dv) + (1.0/tauphi)*(feq_Node[dv] - gridg.Node(i,j, k,dv))

                                        ;
                }

                //< CELLS    

                get_moments_Cell_f(gridf, lb,  ux, uy, uz,Rho, i, j, k,Force); 

                get_equi_g(feq_Cell,lb, ux, uy,uz, phi.Cell(i,j,k),gamma_s,mu.Cell(i,j,k));

                //> normal
                for (int dv = 0; dv< gridg.d_v; dv++){
                    gridg.Cell(i,j,k,dv) =  gridg.Cell(i,j,k,dv) + (1.0/tauphi)*(feq_Cell[dv] - gridg.Cell(i,j,k,dv))

                                        ;
                }       


            }
        }
    }

    






}
    


template<typename T>
void calculateAlpha(lbmD3Q35<T> &lbModel,T* x_i,T* f_i,T beta,T& alpha)
{
  real a(0.0), b(0.0), c(0.0),oneBySix(1.0/6.0);
  real ximin(0.0), ximax(0.0);

  alignas(32) T xSq [35];
  alignas(32) T fxSq[35];

  for(int dv = 0;dv<35;dv++)
  { 
    ximin = std::min(ximin, x_i[dv])  ;
    ximax = std::max(ximax, x_i[dv])  ;
  }

  for(int dv = 0;dv<35;dv++)
  {
    xSq[dv]  = x_i[dv]*x_i[dv] ;
    fxSq[dv] = f_i[dv]*x_i[dv]*x_i[dv] ;

    if(x_i[dv]<0.0)
      a += fxSq[dv]*x_i[dv]*0.5 ;

    b += fxSq[dv]*0.5 ;
    c += fxSq[dv]/(1.0 + 0.5*x_i[dv]) ;
  }

  T alphaMax = -1.0/(beta*ximin);

  T k;
  if(a<0 && b>0 && c>0)
    k = (b-sqrt(b*b - 4.0*a*c))/(2.0*a);
  else
    k = 1.5;

  a = 0.0;b=0.0;c=0.0;
//   T kBeta   = k*beta;
  T beta2   = beta*beta;
  T fourByK = 4.0/k;
  T hBeta = 0.0;

  for(int dv = 0; dv < 35; dv++)
  {
    if(x_i[dv]<0.0)
    {
      a += fxSq[dv]*x_i[dv]*beta2*( 1.0/6.0 - hBeta*x_i[dv]/12.0 + hBeta*hBeta*x_i[dv]*x_i[dv]/20 - hBeta*hBeta*hBeta*x_i[dv]*x_i[dv]*x_i[dv]/5.0 );
      b += fxSq[dv]*0.5;
    }

    if(x_i[dv]>0.0)
      b += f_i[dv]*( (x_i[dv]*x_i[dv]*0.5) - beta2*(x_i[dv]*x_i[dv]*x_i[dv]/15.0)* ( (4.0/(fourByK+x_i[dv])) + (2.0/(fourByK+2.0*x_i[dv])) + (4.0/(fourByK+3.0*x_i[dv])) ));

    c += f_i[dv]*(60.0*x_i[dv]*x_i[dv] + 60.0*x_i[dv]*x_i[dv]*x_i[dv] + 11.0*x_i[dv]*x_i[dv]*x_i[dv]*x_i[dv])/( 60.0 + 90.0*x_i[dv] + 36.0*x_i[dv]*x_i[dv] + 3.0*x_i[dv]*x_i[dv]*x_i[dv]);
  }

  T  h;
  if(a<0 && b>0 && c>0)
    h = (b-std::sqrt(b*b - 4.0*a*c))/(2.0*a);
  else
    h = 2.1;

  a = 0.0;
  hBeta = h*beta;

  for(int dv = 0; dv < 35; dv++)
  {
    if(x_i[dv]<0.0)
    {
      a += f_i[dv]*x_i[dv]*x_i[dv]*x_i[dv]*beta2*( 1.0/6.0 - hBeta*x_i[dv]/12.0 + hBeta*hBeta*x_i[dv]*x_i[dv]/20 - hBeta*hBeta*hBeta*x_i[dv]*x_i[dv]*x_i[dv]/5.0 );
    }
  }

  if(a<0 && b>0 && c>0)
    alpha = 2.0*c/(b+sqrt(b*b - 4.0*a*c));

  if(alpha > alphaMax)
  {
    if ( alphaMax > 1.0)
      alpha = 0.5*(1.0+alphaMax) ;
    else
      alpha = 0.95*alphaMax;
  }
}



template<typename T>
void get_equi_f(real *feq , lbmD3Q35<T> &lb, real ux, real uy, real uz, real &rho, real &phi, real &mu){

    real u2 = ux*ux + uy*uy + uz*uz;
    real a1=0;
    real first,second, third, chem_pot ,feq0=0;
    real sum = 0;
    for (int dv = 1; dv< 35; dv++){

        feq0 = rho*lb.W[dv];

        chem_pot = (phi* mu) / (rho* lb.theta0) ;
        chem_pot = 0;

        first  = (ux*lb.Cx[dv] + uy*lb.Cy[dv] + uz*lb.Cz[dv])*lb.thetaInverse;
        second = 0.5*(first * first);
        third = -0.5*u2*lb.thetaInverse;
        feq[dv] = feq0*(1+ chem_pot +first + second + third);    

        sum += feq[dv];
    }

    feq[0] = rho - sum;
         
}

template<typename T>
void get_equi_g(real *feq , lbmD3Q35<T> &lb, real ux, real uy, real uz, real &phi,  real &gamma, real &mu){


    real u2 = ux*ux + uy*uy + uz*uz;
    real first,second, third, chem_pot ,feq0=0;
    real sum = 0;

    for (int dv = 1; dv< 35; dv++){

        feq0 = lb.W[dv];

        chem_pot = (gamma* mu) / (lb.theta0) ;

        first  = (ux*lb.Cx[dv] + uy*lb.Cy[dv] + uz*lb.Cz[dv])*lb.thetaInverse;
        second = 0.5*(first * first) ;
        third = - 0.5*u2*lb.thetaInverse;
        feq[dv] = feq0*( chem_pot + phi*  first +phi * second + phi * third);    

        sum += feq[dv];
    }
    
    
    feq[0] = phi - sum;

}





template<typename T,typename T1>
void get_moments_Node_f(Grid_N_C_3D<T> &grid, lbmD3Q35<T1> &lb,real &Ux, real &Uy, real &Uz,real &Rho,  int X, int Y, int Z, Grid_N_C_3D<T> &Force){ ///node or cell 0-Node 1- cell
    Ux  = 0.0;
    Uy  = 0.0;
    Uz  = 0.0;
    Rho = 0.0;

    for(int dv = 0; dv <grid.d_v; dv++){
        Ux  += grid.Node(X,Y,Z,dv)*lb.Cx[dv];
        Uy  += grid.Node(X,Y,Z,dv)*lb.Cy[dv];
        Uz  += grid.Node(X,Y,Z,dv)*lb.Cz[dv];
        Rho += grid.Node(X,Y,Z,dv);
    }  

    Ux = Ux/Rho + 0.5*Force.Node(X,Y,Z,0);
    Uy = Uy/Rho + 0.5*Force.Node(X,Y,Z,1);
    Uz = Uz/Rho + 0.5*Force.Node(X,Y,Z,2);  

}

template<typename T,typename T1>
void get_moments_Cell_f(Grid_N_C_3D<T> &grid, lbmD3Q35<T1> &lb,real &Ux, real &Uy, real &Uz,real &Rho,  int X, int Y, int Z,Grid_N_C_3D<T> &Force){ ///node or cell 0-Cell 1- cell
    Ux  = 0.0;
    Uy  = 0.0;
    Uz  = 0.0;
    Rho = 0.0;



    for(int dv = 0; dv <grid.d_v; dv++){
        Ux  += grid.Cell(X,Y,Z,dv)*lb.Cx[dv];
        Uy  += grid.Cell(X,Y,Z,dv)*lb.Cy[dv];
        Uz  += grid.Cell(X,Y,Z,dv)*lb.Cz[dv];
        Rho += grid.Cell(X,Y,Z,dv);
    }  

    Ux = Ux/Rho + 0.5*Force.Cell(X,Y,Z,0);
    Uy = Uy/Rho + 0.5*Force.Cell(X,Y,Z,1);
    Uz = Uz/Rho + 0.5*Force.Cell(X,Y,Z,2);  

}




template<typename T,typename T1>
void get_moments_Node_g(Grid_N_C_3D<T> &grid, lbmD3Q35<T1> &lb,real &phi,  int X, int Y, int Z){ ///node or cell 0-Node 1- cell
    
    phi = 0.0;


    for(int dv = 0; dv <grid.d_v; dv++){
        
        phi += grid.Node(X,Y,Z,dv);
    }  

}

template<typename T,typename T1>
void get_moments_Cell_g(Grid_N_C_3D<T> &grid, lbmD3Q35<T1> &lb,real &phi,  int X, int Y, int Z){ ///node or cell 0-Cell 1- cell

    phi = 0.0;

    for(int dv = 0; dv <grid.d_v; dv++){
        
        phi += grid.Cell(X,Y,Z,dv);
    }  


}


//period no of waves
template<typename T, typename T1>
void initialization(Grid_N_C_3D<T> &grid,lbmD3Q35<T1> &lb,real Rho_mean ,real amplitude, real period){


    real Feq_node[35] = {0},Feq_cell[35] = {0},Rho = 0.0;
    real x,y,z
           ;    ///distance between nodes 
    
    real  x_0 = 0.5;
    real  y_0 = 0.5;
    real  z_0 = 0.5;

    real ux_node = 0., uy_node = 0, uz_node = 0;
    
    real k_w = 2*M_PI* period;
    
    for(int i = 0 + grid.noghost; i < grid.n_x_node - (grid.noghost); i++){
        for(int j = 0 + grid.noghost; j < grid.n_y_node - (grid.noghost); j++){
            for(int k = 0 + grid.noghost; k < grid.n_z_node - (grid.noghost); k++){

                x = ((real)i)/ grid.n_x;
                y = ((real)j)/ grid.n_y;
                z = ((real)k)/ grid.n_z;

				Rho = 1.0;
				get_equi(Feq_node,lb,ux_node,uy_node,uz_node,Rho_mean);


				for (int dv = 0; dv<grid.d_v; dv++)
					grid.Node(i,j,k,dv) = Feq_node[dv];


                x = (((real)i) + 0.5)/ grid.n_x;
                y = (((real)j) + 0.5)/ grid.n_y;
                z = (((real)k) + 0.5)/ grid.n_z;

            
				Rho = 1.0;


				get_equi(Feq_node,lb,ux_node,uy_node,uz_node,Rho_mean);

				for (int dv = 0; dv<grid.d_v; dv++)
					grid.Cell(i,j,k,dv) = Feq_node[dv];
				


                
            }
        }
    }
}



template<typename T, typename T1>
void initialization_equilibrium_profile(Grid_N_C_3D<T> &grid,lbmD3Q35<T1> &lb,real Rho_mean ){

	real Feq_node[35] = {0},Feq_cell[35] = {0},Rho = 0.0;
    real x,y,z
           ;    ///distance between nodes 
    
    real  x_0 = 0.5;
    real  y_0 = 0.5;
    real  z_0 = 0.5;

    real ux_node = 0., uy_node = 0, uz_node = 0;
    
    for(int i = 0 + grid.noghost; i < grid.n_x_node - (grid.noghost); i++){
        for(int j = 0 + grid.noghost; j < grid.n_y_node - (grid.noghost); j++){
            for(int k = 0 + grid.noghost; k < grid.n_z_node - (grid.noghost); k++){

                x = ((real)i)/ grid.n_x - x_0;
                y = ((real)j)/ grid.n_y - y_0;
                z = ((real)k)/ grid.n_z - z_0;

				double rho_gas = 0.4227;
				double rho_liq = 1.6572;

				Rho = (rho_liq + rho_gas)* 0.5 + (rho_liq - rho_gas) *0.5* tanh(x);

				get_equi(Feq_node,lb,ux_node,uy_node,uz_node,Rho);

				for (int dv = 0; dv<grid.d_v; dv++)
					grid.Node(i,j,k,dv) = Feq_node[dv];

				x = ((real)i+0.5)/ grid.n_x - x_0;
                y = ((real)j+0.5)/ grid.n_y - y_0;
                z = ((real)k+0.5)/ grid.n_z - z_0;

				Rho = (rho_liq + rho_gas)* 0.5 + (rho_liq - rho_gas) *0.5* tanh(x);

				get_equi(Feq_node,lb,ux_node,uy_node,uz_node,Rho);

				for (int dv = 0; dv<grid.d_v; dv++)
					grid.Cell(i,j,k,dv) = Feq_node[dv];


            }
        }
    }
}
        


template<typename T, typename T1>
void initialization_2D_droplet(Grid_N_C_3D<T> &gridf,Grid_N_C_3D<T> &gridg,Grid_N_C_3D<T> &rho,
Grid_N_C_3D<T> &phi,Grid_N_C_3D<T> &mu,Grid_N_C_3D<T> &laplacian_phi,lbmD3Q35<T1> &lb,real Rho_mean, real kappa, real gamma_s, real A ){

	real Feq_node[35] = {0},Feq_cell[35] = {0};
    real x,y,z;     ///distance between nodes 
    

    real  x_0 = 0.5;
    real  y_0 = 0.5;
    real  z_0 = 0.5;

    real ux_node = 0.0, uy_node = 0.0, uz_node = 0.0;
    

    for(int i = 0 + gridf.noghost; i < gridf.n_x_node - (gridf.noghost); i++){
        for(int j = 0 + gridf.noghost; j < gridf.n_y_node - (gridf.noghost); j++){
            for(int k = 0 + gridf.noghost; k < gridf.n_z_node - (gridf.noghost); k++){

                x = ((real)i)/ gridf.n_x - x_0;
                y = ((real)j)/ gridf.n_y - y_0;
                z = ((real)k)/ gridf.n_z - z_0;
                
                phi.Node(i,j,k) = -1.0;

                if( x * x  + y * y < 0.25*0.25 ){

                    phi.Node(i,j,k) = 1.0;

                }

                x = ((real)i + 0.5)/ gridf.n_x - x_0;
                y = ((real)j + 0.5)/ gridf.n_y - y_0;
                z = ((real)k + 0.5)/ gridf.n_z - z_0;

                phi.Cell(i,j,k) = -1.0;

                if( x * x  + y * y < 0.25*0.25 ){

                    phi.Cell(i,j,k) = 1.0;

                }

            }
        }
    }

    Periodic(phi);


    for(int i = 0 + gridf.noghost; i < gridf.n_x_node - (gridf.noghost); i++){
        for(int j = 0 + gridf.noghost; j < gridf.n_y_node - (gridf.noghost); j++){
            for(int k = 0 + gridf.noghost; k < gridf.n_z_node - (gridf.noghost); k++){


                //> laplacian of phi  
                
                real del_t = 1.0;
                real Coeff = (2.0/(del_t*del_t*lb.theta0));
                
                laplacian_phi.Node(i,j,k)  = 0.0;

                for(int dv = 0; dv< 27; dv++)
                    laplacian_phi.Node(i,j,k) += lb.W[dv]*phi.Node( i+ (int)lb.Cx[dv]  , j + (int)lb.Cy[dv] , k + (int)lb.Cz[dv] ) ;

                for(int dv = 27; dv< 35;dv++)
                    laplacian_phi.Node(i,j,k) += lb.W[dv]*phi.Cell( i+ (int)lb.CxF[dv] , j + (int)lb.CyF[dv], k + (int)lb.CzF[dv]) ;

                laplacian_phi.Node(i,j,k) = Coeff * ( laplacian_phi.Node(i,j,k) - phi.Node(i,j,k));

                laplacian_phi.Cell(i,j,k)  = 0.0;
                
                for(int dv = 0; dv< 27; dv++)
                    laplacian_phi.Cell(i,j,k) += lb.W[dv]*phi.Cell( i+ (int)lb.Cx[dv]  , j + (int)lb.Cy[dv] , k + (int)lb.Cz[dv] ) ;

                for(int dv = 27; dv< 35;dv++)
                    laplacian_phi.Cell(i,j,k) += lb.W[dv]*phi.Node( i+ (int)lb.CxC[dv] , j + (int)lb.CyC[dv], k + (int)lb.CzC[dv]) ;

                laplacian_phi.Cell(i,j,k) = Coeff * ( laplacian_phi.Cell(i,j,k) - phi.Cell(i,j,k));
                
                //>-------------------------------------------<\\


                //$------------------Node

                real    mu = - A * phi.Node(i,j,k) +  A * phi.Node(i,j,k) *  phi.Node(i,j,k) * phi.Node(i,j,k)   ;
                        mu -= kappa*laplacian_phi.Node(i,j,k);


				get_equi_f(Feq_node,lb,ux_node,uy_node,uz_node  ,Rho_mean, phi.Node(i,j,k),mu);

				for (int dv = 0; dv<gridf.d_v; dv++)
					gridf.Node(i,j,k,dv) = Feq_node[dv];

				get_equi_g(Feq_node,lb,ux_node,uy_node,uz_node  ,phi.Node(i,j,k), gamma_s, mu);

				for (int dv = 0; dv<gridf.d_v; dv++)
					gridg.Node(i,j,k,dv) = Feq_node[dv];
                

                //$-------------------Cell
                mu = - A * phi.Cell(i,j,k) +  A * phi.Cell(i,j,k) *  phi.Cell(i,j,k) * phi.Cell(i,j,k)   ;
                mu -= kappa*laplacian_phi.Cell(i,j,k);

				get_equi_f(Feq_node,lb,ux_node,uy_node,uz_node  ,Rho_mean, phi.Cell(i,j,k),mu);

				for (int dv = 0; dv<gridf.d_v; dv++)
					gridf.Cell(i,j,k,dv) = Feq_node[dv];

				get_equi_g(Feq_node,lb,ux_node,uy_node,uz_node  ,phi.Cell(i,j,k), gamma_s, mu);

				for (int dv = 0; dv<gridf.d_v; dv++)
					gridg.Cell(i,j,k,dv) = Feq_node[dv];
                
                
            
            }
        }
    }




}
        



template<typename T>
void printMass(Grid_N_C_3D<T> &grid){    
    real a = 0;
    for(int i = 0 + grid.noghost; i < grid.n_x_node - (grid.noghost); i++){
        for(int j = 0 + grid.noghost; j < grid.n_y_node - (grid.noghost); j++){
            for(int k = 0 + grid.noghost; k < grid.n_z_node - (grid.noghost); k++){

                for (int dv = 0; dv< grid.d_v; dv++){
                    a += grid.Node(i,j,k,dv) ;
                    a += grid.Cell(i,j,k,dv) ;

                }
            }
        }
    }
std::cout<<"Total mass =  "<<a<<std::endl;

}


;














                // //> with entropic 
                // //<this works properly
                // //$for node
                // real x_i[35];
                // for(int dv = 0; dv< gridf.d_v; dv++)
                //     x_i[dv] = feq_Node[dv]/gridf.Node(i,j,dv) - 1.0;
                
                // real alpha = 0;
                // for(int dv = 0; dv< gridf.d_v; dv++){
                //     alpha = 2.0;
                //     if( std::fabs(x_i[dv]) > 0.0001){
                //         calculateAlpha(lb,x_i,&gridf.Node(i,j,0),beta,alpha);
                //         break;
                //     }
                // }

                // for (int dv = 0; dv< 35; dv++){
                //     gridf.Node(i,j,k,dv) =  gridf.Node(i,j,k,dv) + alpha* beta*(feq_Node[dv] - gridf.Node(i,j,k,dv))
                //                         + (1 - 0.5*alpha*beta)*lb.thetaInverse * rho.Node(i,j,k)* lb.W[dv] * (Force.Node(i,j,k,0) * lb.Cx[dv] + Force.Node(i,j,k,1) * lb.Cy[dv] + Force.Node(i,j,k,2) * lb.Cz[dv] );
                // }    
    




