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
#include "collison.h"

template<typename T, typename T1>
void Multiphase_terms(  Grid_N_C_2D<T> &grid, Grid_N_C_2D<T> &Force, Grid_N_C_2D<T> &rho, Grid_N_C_2D<T> &pnid, Grid_N_C_2D<T> &fnid, Grid_N_C_2D<T> &munid,
                        Grid_N_C_2D<T> &laplacian_rho,Grid_N_C_2D<T> &laplacian_munid,  Grid_N_C_2D<T> &laplacian_fnid, Grid_N_C_2D<T> &gradient_rho, Grid_N_C_2D<T> &P_tensor,
            lbmD2Q9<T1> &lb9, real TbyTc, real kappa,real &sigma, real a , real b, real dx , real dt ){

    real ux = 0, uy = 0;


    for(int i = 0 + grid.noghost; i < grid.n_x_node - (grid.noghost) ; i++){
        for(int j = 0 + grid.noghost;j < grid.n_y_node - (grid.noghost) ; j++){
            
            get_moments_Node(grid, lb9,  ux, uy,  rho.Node(i,j), i, j , Force,dx,dt); 

            real eta = rho.Node(i,j)*b/4.0;           

            
        }
    }


    Periodic_left_Right(rho);
    Periodic_top_bottom(rho);

    // Grad_zero_left_Right(rho);
    // Grad_zero_top_bottom(rho);



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

            double eta = rho.Node(i,j)*b/4.0;

            //> for CS
            // munid.Node(i,j) = -2.0*rho.Node(i,j)*a;
            // munid.Node(i,j) += lb9.theta0*(3.0*eta*eta*eta - 9.0*eta*eta + 8.0*eta) / pow(1.0 - eta, 3.0);
            // munid.Node(i,j) -= kappa_node*laplacian_rho.Node(i,j);

            //> for VdW
            munid.Node(i,j)  = -lb9.theta0*log(1.0 - rho.Node(i,j)*b) ;               
            munid.Node(i,j) += rho.Node(i,j)*b*lb9.theta0/(1.0 - rho.Node(i,j)*b);   
            munid.Node(i,j) -= 2.0*rho.Node(i,j)*a;                                   
            munid.Node(i,j) -= kappa_node*laplacian_rho.Node(i,j);                   


            
        }    
    }


    Periodic_left_Right(laplacian_rho);
    Periodic_top_bottom(laplacian_rho);

    // Grad_zero_left_Right(laplacian_rho);
    // Grad_zero_top_bottom(laplacian_rho);


    Periodic_left_Right(munid);
    Periodic_top_bottom(munid);

    // Grad_zero_left_Right(munid);
    // Grad_zero_top_bottom(munid);





}



template<typename T, typename T1>
void Multiphase_Force_Node(Grid_N_C_2D<T> &grid, Grid_N_C_2D<T> &rho, Grid_N_C_2D<T> &pnid, Grid_N_C_2D<T> &fnid, Grid_N_C_2D<T> &munid, Grid_N_C_2D<T> &laplacian_rho,
            lbmD2Q9<T1> &lb9, Grid_N_C_2D<T> &Force, int i, int j, real kappa, real a , real b, real g, real dx , real dt){


    //> CHEMICAL POTENTIAL FORMULATION 

    Force.Node(i,j,0) = 0.0;
    Force.Node(i,j,1) = 0.0;


    real Coeff_grad = (1.0/(dt*lb9.theta0));

    for(int dv = 0; dv< 9; dv++){
        Force.Node(i,j,0) += lb9.W[dv]*lb9.Cx[dv]*munid.Node( i+ (int)lb9.Cx[dv] , j + (int)lb9.Cy[dv]) ;
        Force.Node(i,j,1) += lb9.W[dv]*lb9.Cy[dv]*munid.Node( i+ (int)lb9.Cx[dv] , j + (int)lb9.Cy[dv]) ;
    }

    Force.Node(i,j,0) = - Coeff_grad*(Force.Node(i,j,0))       ;
    Force.Node(i,j,1) = - Coeff_grad*(Force.Node(i,j,1)) - g   ;

}























































