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
void Multiphase_terms(  Grid_N_C_2D<T> &gridf,Grid_N_C_2D<T> &gridg, lbmD2Q9<T1> &lb9,Grid_N_C_2D<T> &rho,Grid_N_C_2D<T> &phi, Grid_N_C_2D<T> &laplacian_rho,
                        Grid_N_C_2D<T> &psi_phi, Grid_N_C_2D<T> &psi_rho,Grid_N_C_2D<T> &grad_psi_phi, Grid_N_C_2D<T> &grad_psi_rho,
                        double rho_l, double rho_h, double phi_l, double phi_h, double a, double b ){



   // need to define a and b and check hte definitions of psi_phi and psi_phi;

    for(int i = 0 + gridf.nbx ; i <= gridf.nex; i++){
        for(int j = 0 + gridf.nby ; j <= gridf.ney; j++){

            get_phi(gridf,lb9,phi.Node(i,j),i,j);

            rho.Node(i,j) = rho_l + ((phi.Node(i,j) - phi_l)/(phi_h - phi_l)) *(rho_h - rho_l);

            double eta = rho.Node(i,j)*b/4.0;

            psi_phi.Node(i,j) = rho.Node(i,j) *lb9.theta0*(1 + eta + eta*eta - eta*eta*eta)/pow(1 - eta, 3)  - a *rho.Node(i,j)*rho.Node(i,j) - phi.Node(i,j)*lb9.theta0;

            psi_rho.Node(i,j) = rho.Node(i,j) *lb9.theta0*(1 + eta + eta*eta - eta*eta*eta)/pow(1 - eta, 3)  - a *rho.Node(i,j)*rho.Node(i,j) - rho.Node(i,j)*lb9.theta0;


        }
    }

    Periodic_left_Right(rho);
    // Periodic_top_bottom(rho);

    Periodic_left_Right(psi_phi);
    // Periodic_top_bottom(psi_phi);

    Periodic_left_Right(psi_rho);
    // Periodic_top_bottom(psi_rho);



    // Grad_zero_left_Right(rho);
    Grad_zero_top_bottom(rho);

    // Grad_zero_left_Right(psi_phi);
    Grad_zero_top_bottom(psi_phi);

    // Grad_zero_left_Right(psi_rho);
    Grad_zero_top_bottom(psi_rho);


    //  $ LAPLACIAN OF RHO
    for(int i = 0 + gridf.nbx ; i <= gridf.nex; i++){
        for(int j = 0 + gridf.nby ; j <= gridf.ney; j++){
        
            //> laplacian of Rho  
            double del_t = 1.0;
            double Coeff = (2.0/(del_t*del_t*lb9.theta0));
            
            laplacian_rho.Node(i,j)  = 0.0;

            for(int dv = 0; dv< gridf.d_v; dv++)
                laplacian_rho.Node(i,j) += lb9.W[dv]*rho.Node( i+ (int)lb9.Cx[dv]  , j + (int)lb9.Cy[dv] ) ;

            laplacian_rho.Node(i,j) = Coeff * ( laplacian_rho.Node(i,j) - rho.Node(i,j));

     
        }
    }


    Periodic_left_Right(laplacian_rho);
    // Periodic_top_bottom(laplacian_rho);

    // Grad_zero_left_Right(laplacian_rho);
    Grad_zero_top_bottom(laplacian_rho);



    double del_t = 1.0;
    double Coeff_grad = (1.0/(del_t*lb9.theta0));

    for(int i = 0 + gridf.nbx ; i <= gridf.nex; i++){
        for(int j = 0 + gridf.nby ; j <= gridf.ney; j++){


            //$ grad_psi_phi
            grad_psi_phi.Node(i,j,0) = 0;
            grad_psi_phi.Node(i,j,1) = 0;

            for(int dv = 0; dv< gridf.d_v; dv++){

                grad_psi_phi.Node(i,j,0) += lb9.W[dv]*lb9.Cx[dv]*psi_phi.Node( i+ (int)lb9.Cx[dv] , j + (int)lb9.Cy[dv] );
                grad_psi_phi.Node(i,j,1) += lb9.W[dv]*lb9.Cy[dv]*psi_phi.Node( i+ (int)lb9.Cx[dv] , j + (int)lb9.Cy[dv] );

            }
            grad_psi_phi.Node(i,j,0) =  Coeff_grad*(grad_psi_phi.Node(i,j,0));
            grad_psi_phi.Node(i,j,1) =  Coeff_grad*(grad_psi_phi.Node(i,j,1));


            //$   grad_psi_rho
            grad_psi_rho.Node(i,j,0) = 0;
            grad_psi_rho.Node(i,j,1) = 0;

            for(int dv = 0; dv< gridf.d_v; dv++){

                grad_psi_rho.Node(i,j,0) += lb9.W[dv]*lb9.Cx[dv]*psi_rho.Node( i+ (int)lb9.Cx[dv] , j + (int)lb9.Cy[dv] );
                grad_psi_rho.Node(i,j,1) += lb9.W[dv]*lb9.Cy[dv]*psi_rho.Node( i+ (int)lb9.Cx[dv] , j + (int)lb9.Cy[dv] );

            }
            grad_psi_rho.Node(i,j,0) =  Coeff_grad*(grad_psi_rho.Node(i,j,0));
            grad_psi_rho.Node(i,j,1) =  Coeff_grad*(grad_psi_rho.Node(i,j,1));



        }
    }


    Periodic_left_Right(grad_psi_phi);
    // Periodic_top_bottom(grad_psi_phi);

    Periodic_left_Right(grad_psi_rho);
    // Periodic_top_bottom(grad_psi_rho);



    // Grad_zero_left_Right(grad_psi_phi);
    Grad_zero_top_bottom(grad_psi_phi);

    // Grad_zero_left_Right(grad_psi_rho);
    Grad_zero_top_bottom(grad_psi_rho);

}



template<typename T,typename T1>
void get_Force_and_gravity( Grid_N_C_2D<T> &gridf, lbmD2Q9<T1> &lb9,Grid_N_C_2D<T> &laplacian_rho,
                            double kappa, double g,  Grid_N_C_2D<T> &Force,  int i, int j){ ///node or cell 0-Node 1- cell

    double del_t = 1.0;
    double Coeff_grad = (1.0/(del_t*lb9.theta0));

    Force.Node(i,j,0) = 0;
    Force.Node(i,j,1) = 0;

    for(int dv = 0; dv< gridf.d_v; dv++){

        Force.Node(i,j,0)  += lb9.W[dv]*lb9.Cx[dv]*laplacian_rho.Node( i+ (int)lb9.Cx[dv] , j + (int)lb9.Cy[dv] );
        Force.Node(i,j,1)  += lb9.W[dv]*lb9.Cy[dv]*laplacian_rho.Node( i+ (int)lb9.Cx[dv] , j + (int)lb9.Cy[dv] );

    }

    Force.Node(i,j,0) =  kappa*Force.Node(i,j,0)    ;
    Force.Node(i,j,1) =  kappa*Force.Node(i,j,1) - g;



}