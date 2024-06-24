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
void Multiphase_terms(Grid_N_C_2D<T> &grid, Grid_N_C_2D<T> &Force, Grid_N_C_2D<T> &rho, Grid_N_C_2D<T> &pnid, Grid_N_C_2D<T> &fnid, Grid_N_C_2D<T> &munid, Grid_N_C_2D<T> &laplacian_rho,  Grid_N_C_2D<T> &laplacian_fnid, Grid_N_C_2D<T> &gradient_rho, 
            lbmD2Q9<T1> &lb9, real TbyTc, real kappa, real a , real b ){

    real ux = 0, uy = 0;



    for(int i = 0 + grid.noghost; i < grid.n_x_node - (grid.noghost) ; i++){
        for(int j = 0 + grid.noghost;j < grid.n_y_node - (grid.noghost) ; j++){
            
            get_moments_Node(grid, lb9,  ux, uy,  rho.Node(i,j), i, j , Force); 


            real eta = rho.Node(i,j)*b/4.0;

            //> pnid Node
            pnid.Node(i,j) = (rho.Node(i,j)*lb9.theta0*(1.0 + eta + eta* eta - eta*eta*eta) )/pow(1.0 - eta,  3.0)  - 
                                a * rho.Node(i,j)*rho.Node(i,j) ;

            //> Fnid_node
            fnid.Node(i,j) = -1.0*(rho.Node(i,j)*lb9.theta0*(3.0*eta*eta - 4.0 * eta) )/pow(1.0 - eta,  2.0)  - 
                                a * rho.Node(i,j)*rho.Node(i,j) ;

        }
    }



    Periodic_left_Right(rho);
    Periodic_top_bottom(rho);

    Grad_zero_left_Right(rho);
    Grad_zero_top_bottom(rho);

    //>  -----------------------
    Periodic_left_Right(pnid);
    Periodic_top_bottom(pnid);

    Grad_zero_left_Right(pnid);
    Grad_zero_top_bottom(pnid);

    //>--------------------------
    Periodic_left_Right(fnid);
    Periodic_top_bottom(fnid);
    
    Grad_zero_left_Right(fnid);
    Grad_zero_top_bottom(fnid);





    for(int i = 0 + grid.noghost; i < grid.n_x_node - (grid.noghost) ; i++){
        for(int j = 0 + grid.noghost;j < grid.n_y_node - (grid.noghost) ; j++){



            //> laplacian of Rho  
            real del_t = 1.0;
            real Coeff = (2.0/(del_t*del_t*lb9.theta0));
            
            laplacian_rho.Node(i,j)  = 0.0;

            for(int dv = 0; dv<  grid.d_v; dv++)
                laplacian_rho.Node(i,j) += lb9.W[dv]*rho.Node( i+ (int)lb9.Cx[dv]  , j + (int)lb9.Cy[dv] , k + (int)lb9.Cz[dv] ) ;


            laplacian_rho.Node(i,j) = Coeff * ( laplacian_rho.Node(i,j) - rho.Node(i,j));






            // //>------
            real kappa_node = kappa;

            // kappa_node = kappa - (1.0/2.0)* lb9.theta0* (
            //                                         -32.0 *b* lb9.theta0*(-16.0 + b* rho.Node(i,j)) / (pow(-4.0 + b*rho.Node(i,j) , 4.0))
            //                                             - 2.0*a
            //                                         ) 
            //                                     ;

            double eta = rho.Node(i,j)*b/4.0;

            //> munid
            munid.Node(i,j) = -2.0*rho.Node(i,j)*a;
            munid.Node(i,j) += lb9.theta0*(3.0*eta*eta*eta - 9.0*eta*eta + 8.0*eta) / pow(1.0 - eta, 3.0);
            munid.Node(i,j) -= kappa_node*laplacian_rho.Node(i,j);





            
        }    
    }



    Periodic_left_Right(laplacian_rho);
    Periodic_top_bottom(laplacian_rho);

    Grad_zero_left_Right(laplacian_rho);
    Grad_zero_top_bottom(laplacian_rho);


    Periodic_left_Right(laplacian_fnid);
    Periodic_top_bottom(laplacian_fnid);

    Grad_zero_left_Right(laplacian_fnid);
    Grad_zero_top_bottom(laplacian_fnid);

    Periodic_left_Right(munid);
    Periodic_top_bottom(munid);

    Grad_zero_left_Right(munid);
    Grad_zero_top_bottom(munid);


}




template<typename T, typename T1>
void Multiphase_Force_Node(Grid_N_C_2D<T> &grid, Grid_N_C_2D<T> &rho, Grid_N_C_2D<T> &pnid, Grid_N_C_2D<T> &fnid, Grid_N_C_2D<T> &munid, Grid_N_C_2D<T> &laplacian_rho,
            lbmD2Q9<T1> &lb9, Grid_N_C_2D<T> &Force, int i, int j, real kappa, real a , real b, real g){


    //> CHEMICAL POTENTIAL FORMULATION 
    real grad_mux = 0.0, grad_muy = 0.0;

    Force.Node(i,j,k,0) = 0.0;
    Force.Node(i,j,k,1) = 0.0;


    real del_t = 1.0;
    real Coeff_grad = (1.0/(del_t*lb9.theta0));

    for(int dv = 0; dv< 9; dv++){
        grad_mux += lb9.W[dv]*lb9.Cx[dv]*munid.Node( i+ (int)lb9.Cx[dv] , j + (int)lb9.Cy[dv]) ;
        grad_muy += lb9.W[dv]*lb9.Cy[dv]*munid.Node( i+ (int)lb9.Cx[dv] , j + (int)lb9.Cy[dv]) ;
    }


    Force.Node(i,j,k,0) = - Coeff_grad*(grad_mux)        ;
    Force.Node(i,j,k,1) = - Coeff_grad*(grad_muy)  -g    ;


}


;