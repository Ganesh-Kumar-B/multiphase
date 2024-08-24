#pragma once

#include<iostream>
#include<cmath>
#include<vector>
#include<math.h>
#include<fstream>
#include<algorithm>
#include <sstream>
#include<string>
#include "lbmD2Q9.h"
#include "GRID_2D.h"
#include "multiphase.h"
#define PI 3.14159265



template<typename T, typename T1>
void collide(Grid_N_C_2D<T> &grid,
            lbmD2Q9<T1> &lb9, real beta,real tau, real theta_r,real kappa, real &sigma, int t,Grid_N_C_2D<T> &Force, Grid_N_C_2D<T> &P_tensor ,real g, real dx , real dt ){

    Grid_N_C_2D<T>  rho                             (grid.n_x,grid.n_y,1,1);   
    Grid_N_C_2D<T>  munid                           (grid.n_x,grid.n_y,1,1);   
    Grid_N_C_2D<T>  laplacian_rho                   (grid.n_x,grid.n_y,1,1);   



    real feq_Node[9] = {0},
    ux = 0, uy = 0;


    sigma = 0;

    for(int i = 0 + grid.noghost; i < grid.n_x_node - (grid.noghost) ; i++){
        for(int j = 0 + grid.noghost;j < grid.n_y_node - (grid.noghost) ; j++){
            
            get_moments_Node(grid, lb9,  ux, uy,  rho.Node(i,j), i, j , Force,dx,dt); 

        }
    }

    Periodic_left_Right(rho);
    Periodic_top_bottom(rho);


    for(int i = 0 + grid.noghost; i < grid.n_x_node - (grid.noghost) ; i++){
        for(int j = 0 + grid.noghost;j < grid.n_y_node - (grid.noghost) ; j++){

            

            //> laplacian of Rho  
            real Coeff = (2.0/(dt*dt*lb9.theta0));
            
            laplacian_rho.Node(i,j)  = 0.0;

            for(int dv = 0; dv<  grid.d_v; dv++)
                laplacian_rho.Node(i,j) += lb9.W[dv]*rho.Node( i+ (int)lb9.Cx[dv]  , j + (int)lb9.Cy[dv]  ) ;


            laplacian_rho.Node(i,j) = Coeff * ( laplacian_rho.Node(i,j) - rho.Node(i,j));

        



            //  //>------
            real kappa_node = kappa;

            //  //> for CS
            // kappa_node = kappa - (1.0/2.0)*dt*dt* lb9.theta0* (
            //                                         -32.0 *b* lb9.theta0*(-16.0 + b* rho.Node(i,j)) / (pow(-4.0 + b*rho.Node(i,j) , 4.0))
            //                                             - 2.0*a
            //                                         ) 
            //                                     ;
            
            // //  //> for VdW
            // kappa_node = kappa - (1.0/2.0)*dt*dt*lb9.theta0* (
            //                                                     rho.Node(i,j,k) * b * b * lb.theta0/(pow(1.0 - rho.Node(i,j,k)*b ,2)) 
            //                                                     +(b*lb.theta0)/( 1.0 - rho.Node(i,j,k) * b  )   
            //                                                     - 2.0*a
            //                                                 )
            //                                     ;


            //> for CS
            // munid.Node(i,j) = -2.0*rho.Node(i,j)*a;
            // munid.Node(i,j) += lb9.theta0*(3.0*eta*eta*eta - 9.0*eta*eta + 8.0*eta) / pow(1.0 - eta, 3.0);
            // munid.Node(i,j) -= kappa_node*laplacian_rho.Node(i,j);

            //> for VdW
            munid.Node(i,j)  = -(9.0/4.0)*rho.Node(i,j);
            munid.Node(i,j) += theta_r*(-log(3.0 - rho.Node(i,j)) + log(3.0) + (rho.Node(i,j))/(3.0 -rho.Node(i,j)) );
            munid.Node(i,j) -= kappa_node*(laplacian_rho.Node(i,j));
            


            
        }    
    }


    Periodic_left_Right(laplacian_rho);
    Periodic_top_bottom(laplacian_rho);


    Periodic_left_Right(munid);
    Periodic_top_bottom(munid);

    //---------------------------------------------------------------------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------MAIN LOOP-------------------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------------------------------------------------------------------//

    for(int i = 0 + grid.noghost; i < grid.n_x_node - (grid.noghost) ; i++){
        for(int j = 0 + grid.noghost;j < grid.n_y_node - (grid.noghost) ; j++){
            
        
            real Rho = 0.0;

            Force.Node(i,j,0) = 0.0;
            Force.Node(i,j,1) = 0.0;


            real Coeff_grad = (1.0/(dt*lb9.theta0));

            for(int dv = 0; dv< 9; dv++){
                Force.Node(i,j,0) += lb9.W[dv]*lb9.Cx[dv]*munid.Node( i+ (int)lb9.Cx[dv] , j + (int)lb9.Cy[dv]) ;
                Force.Node(i,j,1) += lb9.W[dv]*lb9.Cy[dv]*munid.Node( i+ (int)lb9.Cx[dv] , j + (int)lb9.Cy[dv]) ;
            }

            Force.Node(i,j,0) = - Coeff_grad*(Force.Node(i,j,0)) ;
            Force.Node(i,j,1) = - Coeff_grad*(Force.Node(i,j,1)) ;
                
            
            get_moments_Node(grid, lb9,  ux, uy,Rho, i, j, Force,dx,dt );            //for the node
            get_equi(feq_Node ,lb9, ux, uy, Rho);


            for (int dv = 0; dv< grid.d_v; dv++){
                grid.Node(i,j,dv) =  grid.Node(i,j,dv) + 2.0* beta*(feq_Node[dv] - grid.Node(i,j,dv))
                                    + (1.0 - beta)*dt*lb9.thetaInverse *  feq_Node[dv] * (Force.Node(i,j,0) * (lb9.Cx[dv] - ux) + Force.Node(i,j,1) * (lb9.Cy[dv] - uy) );
                                    ;
            }


        
        }
    }


}
    


template<typename T, typename T1>
void initialization_circle(Grid_N_C_2D<T> &grid,lbmD2Q9<T1> &lb,real rho_liq,real rho_gas,real rho_c, real R ){

	real Feq_node[9] = {0},Rho = 0.0;
    real x,y
           ;    ///distance between nodes 
    
    real  x_0 = 0.0;
    real  y_0 = 0.0;

    real phi_l = -1.0;
    real phi_h = +1.0;

    rho_liq = rho_liq*rho_c;
    rho_gas = rho_gas*rho_c;

    real phi;
    real ux_node = 0.0, uy_node = 0.0;


    x_0 = 0.5, y_0 = 0.5;
    for(int i = 0 + grid.noghost; i < grid.n_x_node - (grid.noghost); i++){
        for(int j = 0 + grid.noghost; j < grid.n_y_node - (grid.noghost); j++){
            
            
            
            x = ((real)i)/ grid.n_x - x_0;
            y = ((real)j)/ grid.n_x - y_0;
            

            Rho = rho_gas;
            if(x*x + y*y < R*R){
                Rho = rho_liq;
            }
            

            get_equi(Feq_node,lb,ux_node,uy_node,Rho);

            for (int dv = 0; dv<grid.d_v; dv++)
                grid.Node(i,j,dv) = Feq_node[dv];


        }
    }




}



template<typename T>
void get_equi(real *feq , lbmD2Q9<T> &lb, real ux, real uy, real rho){


    real u2 = ux*ux + uy*uy ;
    real a1 = 0;
    real first,second, third,feq0=0;
    for (int dv = 0; dv< 9; dv++){

        feq0 = rho*lb.W[dv];

        first  = (ux*lb.Cx[dv] + uy*lb.Cy[dv] )*lb.thetaInverse;
        second = 0.5*(first * first);
        third = -0.5*u2*lb.thetaInverse;
        feq[dv] = feq0*(1+ first + second + third);    

        
    }
         
}





template<typename T,typename T1>
void get_moments_Node(Grid_N_C_2D<T> &grid, lbmD2Q9<T1> &lb,real &Ux, real &Uy, real &Rho,  int X, int Y, Grid_N_C_2D<T> &Force,real dx, real dt){ ///node or cell 0-Node 1- cell
    Ux  = 0.0;
    Uy  = 0.0;
    Rho = 0.0;


    for(int dv = 0; dv <grid.d_v; dv++){
        Ux  += grid.Node(X,Y,dv)*lb.Cx[dv];
        Uy  += grid.Node(X,Y,dv)*lb.Cy[dv];
        Rho += grid.Node(X,Y,dv);
    }  

    Ux = Ux/Rho + 0.5*Force.Node(X,Y,0)*dt;
    Uy = Uy/Rho + 0.5*Force.Node(X,Y,1)*dt;


}











template<typename T>
void printMass(Grid_N_C_2D<T> &gridf){    
    real a = 0;
    for(int i = 0 + gridf.noghost; i < gridf.n_x_node - (gridf.noghost) ; i++){
        for(int j = 0 + gridf.noghost;j < gridf.n_y_node - (gridf.noghost) ; j++){
            for (int dv = 0; dv< 9; dv++){
                a += gridf.Node(i,j,dv) ;
            }
        }
    }
    std::cout<<"   "<<a<<std::endl;

}
