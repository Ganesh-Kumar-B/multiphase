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

enum coodinates{X,Y,Z};

template<typename T, typename T1>
void Multiphase_terms(Grid_N_C_3D<T> &gridf,Grid_N_C_3D<T> &gridg, lbmD3Q35<T1> &lb,  Grid_N_C_3D<T> &rho, Grid_N_C_3D<T> &phi, Grid_N_C_3D<T> &munid, 
            Grid_N_C_3D<T> &laplacian_phi, real kappa, real A
            ){


    for(int i = 0 + grid.noghost; i < grid.n_x_node - (grid.noghost) ; i++){
        for(int j = 0 + grid.noghost;j < grid.n_y_node - (grid.noghost) ; j++){
            for(int k = 0 + grid.noghost;k < grid.n_z_node - (grid.noghost) ; k++){
            
                get_moments_Node_g(gridg, lb, phi.Node(i,j,k), i, j ,k); 

                get_moments_Cell_g(gridg, lb, phi.Cell(i,j,k), i, j ,k); 

            }
        }
    }


    Periodic(phi);



    for(int i = 0 + grid.noghost; i < grid.n_x_node - (grid.noghost) ; i++){
        for(int j = 0 + grid.noghost;j < grid.n_y_node - (grid.noghost) ; j++){
            for(int k = 0 + grid.noghost;k < grid.n_z_node - (grid.noghost) ; k++){


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
                //>---------------------






                //> munid
                munid.Node(i,j,k) = - A * phi.Node(i,j,k) +  A * phi.Node(i,j,k) *  phi.Node(i,j,k) * phi.Node(i,j,k)   ;
                
                munid.Node(i,j,k) -= kappa*laplacian_phi.Node(i,j,k);

    

                munid.Cell(i,j,k) = - A * phi.Cell(i,j,k) +  A * phi.Cell(i,j,k) *  phi.Cell(i,j,k) * phi.Cell(i,j,k)   ;

                munid.Cell(i,j,k) -= kappa*laplacian_phi.Cell(i,j,k);

            }
        }    
    }

    Periodic(munid); 

}






template<typename T, typename T1>
void Multiphase_Force_Node(Grid_N_C_3D<T> &grid, Grid_N_C_3D<T> &phi,   Grid_N_C_3D<T> &mu, 
            lbmD3Q35<T1> &lb, Grid_N_C_3D<T> &Force, int i, int j,int k){

    //> CHEMICAL POTENTIAL FORMULATION 
    real grad_mux = 0.0, grad_muy = 0.0, grad_muz = 0.0;

    Force.Node(i,j,k,0) = 0.0;
    Force.Node(i,j,k,1) = 0.0;
    Force.Node(i,j,k,2) = 0.0;


    real del_t = 1.0;
    real Coeff_grad = (1.0/(del_t*lb.theta0));

    for(int dv = 0; dv< 27; dv++){
        grad_mux += lb.W[dv]*lb.Cx[dv]*phi.Node( i+ (int)lb.Cx[dv]  , j + (int)lb.Cy[dv]  , k + (int)lb.Cz[dv]) ;
        grad_muy += lb.W[dv]*lb.Cy[dv]*phi.Node( i+ (int)lb.Cx[dv]  , j + (int)lb.Cy[dv]  , k + (int)lb.Cz[dv]) ;
        grad_muz += lb.W[dv]*lb.Cz[dv]*phi.Node( i+ (int)lb.Cx[dv]  , j + (int)lb.Cy[dv]  , k + (int)lb.Cz[dv]) ;
    }

    for(int dv = 27; dv<grid.d_v; dv++){
        grad_mux += lb.W[dv]*lb.Cx[dv]*phi.Cell( i+ (int)lb.CxF[dv] , j + (int)lb.CyF[dv] , k + (int)lb.CzF[dv]) ;
        grad_muy += lb.W[dv]*lb.Cy[dv]*phi.Cell( i+ (int)lb.CxF[dv] , j + (int)lb.CyF[dv] , k + (int)lb.CzF[dv]) ;
        grad_muz += lb.W[dv]*lb.Cz[dv]*phi.Cell( i+ (int)lb.CxF[dv] , j + (int)lb.CyF[dv] , k + (int)lb.CzF[dv]) ;
    }

    Force.Node(i,j,k,0) = -mu.Node(i,j,k)*Coeff_grad*(grad_mux);
    Force.Node(i,j,k,1) = -mu.Node(i,j,k)*Coeff_grad*(grad_muy);  
    Force.Node(i,j,k,2) = -mu.Node(i,j,k)*Coeff_grad*(grad_muz);  

}

//this is from the 41 paper formulation
template<typename T, typename T1>
void Multiphase_Force_eta_Node(Grid_N_C_3D<T> &grid, Grid_N_C_3D<T> &rho, Grid_N_C_3D<T> &pnid, Grid_N_C_3D<T> &fnid, Grid_N_C_3D<T> &munid, Grid_N_C_3D<T> &laplacian_rho,  Grid_N_C_3D<T> &laplacian_fnid, Grid_N_C_3D<T> &gradient_rho, 
            lbmD3Q35<T1> &lb, real &Fx, real &Fy, real &Fz, int i, int j,int k, real kappa, real eta){

    //> CHEMICAL POTENTIAL FORMULATION 
    real grad_mux = 0.0, grad_muy = 0.0, grad_muz = 0.0;
    Fx = 0, Fy = 0, Fz = 0.0;
    real del_t = 1.0;
    real Coeff_grad = (1.0/(del_t*lb.theta0));

    //#change the munid to bulk remove the laplacian rho when using this
    //
    for(int dv = 0; dv< 27; dv++){
        grad_mux += rho.Node(i,j,k)*Coeff_grad*lb.W[dv]*lb.Cx[dv]*munid.Node( i+ (int)lb.Cx[dv] , j + (int)lb.Cy[dv], k + (int)lb.Cz[dv]) ;
        grad_muy += rho.Node(i,j,k)*Coeff_grad*lb.W[dv]*lb.Cy[dv]*munid.Node( i+ (int)lb.Cx[dv] , j + (int)lb.Cy[dv], k + (int)lb.Cz[dv]) ;
        grad_muz += rho.Node(i,j,k)*Coeff_grad*lb.W[dv]*lb.Cz[dv]*munid.Node( i+ (int)lb.Cx[dv] , j + (int)lb.Cy[dv], k + (int)lb.Cz[dv]) ;
    }
    for(int dv = 27; dv<grid.d_v; dv++){
        grad_mux += rho.Node(i,j,k)*Coeff_grad*lb.W[dv]*lb.Cx[dv]*munid.Cell( i+ (int)lb.CxF[dv] , j + (int)lb.CyF[dv], k + (int)lb.CzF[dv]) ;
        grad_muy += rho.Node(i,j,k)*Coeff_grad*lb.W[dv]*lb.Cy[dv]*munid.Cell( i+ (int)lb.CxF[dv] , j + (int)lb.CyF[dv], k + (int)lb.CzF[dv]) ;
        grad_muz += rho.Node(i,j,k)*Coeff_grad*lb.W[dv]*lb.Cz[dv]*munid.Cell( i+ (int)lb.CxF[dv] , j + (int)lb.CyF[dv], k + (int)lb.CzF[dv]) ;
    }

    //munid*grad(rho)
        grad_mux += munid.Node(i,j,k)*gradient_rho.Node(i,j,k,0);
        grad_muy += munid.Node(i,j,k)*gradient_rho.Node(i,j,k,1);
        grad_muz += munid.Node(i,j,k)*gradient_rho.Node(i,j,k,2);
        
    //fnid
    for(int dv = 0; dv< 27; dv++){
        grad_mux = grad_mux - eta* ( Coeff_grad*lb.W[dv]*lb.Cx[dv]*fnid.Node( i+ (int)lb.Cx[dv] , j + (int)lb.Cy[dv], k + (int)lb.Cz[dv])) ;
        grad_muy = grad_muy - eta* ( Coeff_grad*lb.W[dv]*lb.Cy[dv]*fnid.Node( i+ (int)lb.Cx[dv] , j + (int)lb.Cy[dv], k + (int)lb.Cz[dv])) ;
        grad_muz = grad_muz - eta* ( Coeff_grad*lb.W[dv]*lb.Cz[dv]*fnid.Node( i+ (int)lb.Cx[dv] , j + (int)lb.Cy[dv], k + (int)lb.Cz[dv])) ;
    }

    for(int dv = 27; dv<grid.d_v; dv++){
        grad_mux = grad_mux - eta* (Coeff_grad*lb.W[dv]*lb.Cx[dv]*fnid.Cell( i+ (int)lb.CxF[dv] , j + (int)lb.CyF[dv], k + (int)lb.CzF[dv])) ;
        grad_muy = grad_muy - eta* (Coeff_grad*lb.W[dv]*lb.Cy[dv]*fnid.Cell( i+ (int)lb.CxF[dv] , j + (int)lb.CyF[dv], k + (int)lb.CzF[dv])) ;
        grad_muz = grad_muz - eta* (Coeff_grad*lb.W[dv]*lb.Cz[dv]*fnid.Cell( i+ (int)lb.CxF[dv] , j + (int)lb.CyF[dv], k + (int)lb.CzF[dv])) ;
    }

    //fnid
   // fourth oder corrections
    real coeff_grad_4th = 0.5*lb.theta0*del_t*del_t   ;
    for(int dv = 0; dv< 27; dv++){
        grad_mux =  grad_mux - (1 - eta) * ((Coeff_grad*lb.W[dv]*lb.Cx[dv]*fnid.Node( i+ (int)lb.Cx[dv] , j + (int)lb.Cy[dv], k + (int)lb.Cz[dv]) - coeff_grad_4th*Coeff_grad*lb.W[dv]*lb.Cx[dv]*laplacian_fnid.Node( i+ (int)lb.Cx[dv] , j + (int)lb.Cy[dv], k + (int)lb.Cz[dv]) )   ) ;                    
        grad_muy =  grad_muy - (1 - eta) * ((Coeff_grad*lb.W[dv]*lb.Cy[dv]*fnid.Node( i+ (int)lb.Cx[dv] , j + (int)lb.Cy[dv], k + (int)lb.Cz[dv]) - coeff_grad_4th*Coeff_grad*lb.W[dv]*lb.Cy[dv]*laplacian_fnid.Node( i+ (int)lb.Cx[dv] , j + (int)lb.Cy[dv], k + (int)lb.Cz[dv]) )   ) ;                    
        grad_muz =  grad_muz - (1 - eta) * ((Coeff_grad*lb.W[dv]*lb.Cz[dv]*fnid.Node( i+ (int)lb.Cx[dv] , j + (int)lb.Cy[dv], k + (int)lb.Cz[dv]) - coeff_grad_4th*Coeff_grad*lb.W[dv]*lb.Cz[dv]*laplacian_fnid.Node( i+ (int)lb.Cx[dv] , j + (int)lb.Cy[dv], k + (int)lb.Cz[dv]) )   ) ;                    
    }

    for(int dv = 27; dv<grid.d_v; dv++){
        grad_mux =  grad_mux -  (1 - eta) *((Coeff_grad*lb.W[dv]*lb.Cx[dv]*fnid.Cell( i+ (int)lb.CxF[dv] , j + (int)lb.CyF[dv], k + (int)lb.CzF[dv]) - coeff_grad_4th*Coeff_grad*lb.W[dv]*lb.Cx[dv]*laplacian_fnid.Cell( i+ (int)lb.CxF[dv] , j + (int)lb.CyF[dv], k + (int)lb.CzF[dv])) );
        grad_muy =  grad_muy -  (1 - eta) *((Coeff_grad*lb.W[dv]*lb.Cy[dv]*fnid.Cell( i+ (int)lb.CxF[dv] , j + (int)lb.CyF[dv], k + (int)lb.CzF[dv]) - coeff_grad_4th*Coeff_grad*lb.W[dv]*lb.Cy[dv]*laplacian_fnid.Cell( i+ (int)lb.CxF[dv] , j + (int)lb.CyF[dv], k + (int)lb.CzF[dv])) );
        grad_muz =  grad_muz -  (1 - eta) *((Coeff_grad*lb.W[dv]*lb.Cz[dv]*fnid.Cell( i+ (int)lb.CxF[dv] , j + (int)lb.CyF[dv], k + (int)lb.CzF[dv]) - coeff_grad_4th*Coeff_grad*lb.W[dv]*lb.Cz[dv]*laplacian_fnid.Cell( i+ (int)lb.CxF[dv] , j + (int)lb.CyF[dv], k + (int)lb.CzF[dv])) );
    }




    //
    for(int dv = 0; dv< 27; dv++){

        int ic = i+ (int)lb.Cx[dv]; int jc = j + (int)lb.Cy[dv];int kc = k + (int)lb.Cz[dv];

        grad_mux += kappa*Coeff_grad*lb.W[dv]*lb.Cx[dv]*( 0.5*pow(gradient_rho.Node( ic ,jc, kc,X), 2 ) - 0.5*pow(gradient_rho.Node( ic ,jc, kc,Y), 2 )  - 0.5*pow(gradient_rho.Node( ic ,jc, kc,Z), 2 ) - rho.Node(ic ,jc, kc) *laplacian_rho.Node(ic ,jc, kc)) 
                 +  kappa*Coeff_grad*lb.W[dv]*lb.Cy[dv]*( gradient_rho.Node( ic ,jc, kc,X) * gradient_rho.Node( ic ,jc, kc,Y)) 
                 +  kappa*Coeff_grad*lb.W[dv]*lb.Cz[dv]*( gradient_rho.Node( ic ,jc, kc,X) * gradient_rho.Node( ic ,jc, kc,Z)) ;



        grad_muy += kappa*Coeff_grad*lb.W[dv]*lb.Cx[dv]*( gradient_rho.Node( ic ,jc, kc,X) * gradient_rho.Node( ic ,jc, kc,Y)) 
                 +  kappa*Coeff_grad*lb.W[dv]*lb.Cy[dv]*( 0.5*pow(gradient_rho.Node( ic ,jc, kc,Y), 2 ) - 0.5*pow(gradient_rho.Node( ic ,jc, kc,X), 2 )  - 0.5*pow(gradient_rho.Node( ic ,jc, kc,Z), 2 ) - rho.Node(ic ,jc, kc) *laplacian_rho.Node(ic ,jc, kc)) 
                 +  kappa*Coeff_grad*lb.W[dv]*lb.Cz[dv]*( gradient_rho.Node( ic ,jc, kc,Y) * gradient_rho.Node( ic ,jc, kc,Z)) ;


        grad_muz += kappa*Coeff_grad*lb.W[dv]*lb.Cx[dv]*( gradient_rho.Node( ic ,jc, kc,Z) * gradient_rho.Node( ic ,jc, kc,X)) 
                 +  kappa*Coeff_grad*lb.W[dv]*lb.Cy[dv]*( gradient_rho.Node( ic ,jc, kc,Z) * gradient_rho.Node( ic ,jc, kc,Y)) 
                 +  kappa*Coeff_grad*lb.W[dv]*lb.Cz[dv]*( 0.5*pow(gradient_rho.Node( ic ,jc, kc,Z), 2 ) - 0.5*pow(gradient_rho.Node( ic ,jc, kc,X), 2 )  - 0.5*pow(gradient_rho.Node( ic ,jc, kc,Y), 2 ) - rho.Node(ic ,jc, kc) *laplacian_rho.Node(ic ,jc, kc)) ;

    }




    for(int dv = 27; dv<grid.d_v; dv++){

        int ic = i+ (int)lb.CxF[dv]; int jc = j + (int)lb.CyF[dv];int kc = k + (int)lb.CzF[dv];


        grad_mux += kappa*Coeff_grad*lb.W[dv]*lb.Cx[dv]*( 0.5*pow(gradient_rho.Cell( ic ,jc, kc,X), 2 ) - 0.5*pow(gradient_rho.Cell( ic ,jc, kc,Y), 2 )  - 0.5*pow(gradient_rho.Cell( ic ,jc, kc,Z), 2 ) - rho.Cell(ic ,jc, kc) *laplacian_rho.Cell(ic ,jc, kc)) 
                 +  kappa*Coeff_grad*lb.W[dv]*lb.Cy[dv]*( gradient_rho.Cell( ic ,jc, kc,X) * gradient_rho.Cell( ic ,jc, kc,Y)) 
                 +  kappa*Coeff_grad*lb.W[dv]*lb.Cz[dv]*( gradient_rho.Cell( ic ,jc, kc,X) * gradient_rho.Cell( ic ,jc, kc,Z)) ;



        grad_muy += kappa*Coeff_grad*lb.W[dv]*lb.Cx[dv]*( gradient_rho.Cell( ic ,jc, kc,X) * gradient_rho.Cell( ic ,jc, kc,Y)) 
                 +  kappa*Coeff_grad*lb.W[dv]*lb.Cy[dv]*( 0.5*pow(gradient_rho.Cell( ic ,jc, kc,Y), 2 ) - 0.5*pow(gradient_rho.Cell( ic ,jc, kc,X), 2 )  - 0.5*pow(gradient_rho.Cell( ic ,jc, kc,Z), 2 ) - rho.Cell(ic ,jc, kc) *laplacian_rho.Cell(ic ,jc, kc)) 
                 +  kappa*Coeff_grad*lb.W[dv]*lb.Cz[dv]*( gradient_rho.Cell( ic ,jc, kc,Y) * gradient_rho.Cell( ic ,jc, kc,Z)) ;


        grad_muz += kappa*Coeff_grad*lb.W[dv]*lb.Cx[dv]*( gradient_rho.Cell( ic ,jc, kc,Z) * gradient_rho.Cell( ic ,jc, kc,X)) 
                 +  kappa*Coeff_grad*lb.W[dv]*lb.Cy[dv]*( gradient_rho.Cell( ic ,jc, kc,Z) * gradient_rho.Cell( ic ,jc, kc,Y)) 
                 +  kappa*Coeff_grad*lb.W[dv]*lb.Cz[dv]*( 0.5*pow(gradient_rho.Cell( ic ,jc, kc,Z), 2 ) - 0.5*pow(gradient_rho.Cell( ic ,jc, kc,X), 2 )  - 0.5*pow(gradient_rho.Cell( ic ,jc, kc,Y), 2 ) - rho.Cell(ic ,jc, kc) *laplacian_rho.Cell(ic ,jc, kc)) ;


    }

    Fx = - grad_mux/rho.Node(i,j,k);
    Fy = - grad_muy/rho.Node(i,j,k);  
    Fz = - grad_muz/rho.Node(i,j,k);  

}




















template<typename T, typename T1>
void Multiphase_Force_Cell(Grid_N_C_3D<T> &grid, Grid_N_C_3D<T> &phi,   Grid_N_C_3D<T> &mu, 
            lbmD3Q35<T1> &lb, Grid_N_C_3D<T> &Force, int i, int j,int k){

    //> CHEMICAL POTENTIAL FORMULATION 
    real grad_mux = 0.0, grad_muy = 0.0, grad_muz = 0.0;
    
    Force.Cell(i,j,k,0) = 0.0;
    Force.Cell(i,j,k,1) = 0.0;
    Force.Cell(i,j,k,2) = 0.0;

    real del_t = 1.0;
    real Coeff_grad = (1.0/(del_t*lb.theta0));
    
    
    //> CHEMICAL POTENTIAL FORMULATION 
    grad_mux = 0.0, grad_muy = 0.0, grad_muz = 0.0;




    del_t = 1.0;
    Coeff_grad = (1.0/(del_t*lb.theta0));

    for(int dv = 0; dv< 27; dv++){
        grad_mux += lb.W[dv]*lb.Cx[dv]*phi.Cell( i+ (int)lb.Cx[dv] , j + (int)lb.Cy[dv], k + (int)lb.Cz[dv]) ;
        grad_muy += lb.W[dv]*lb.Cy[dv]*phi.Cell( i+ (int)lb.Cx[dv] , j + (int)lb.Cy[dv], k + (int)lb.Cz[dv]) ;
        grad_muz += lb.W[dv]*lb.Cz[dv]*phi.Cell( i+ (int)lb.Cx[dv] , j + (int)lb.Cy[dv], k + (int)lb.Cz[dv]) ;
    }

    for(int dv = 27; dv< 35; dv++){
        grad_mux += lb.W[dv]*lb.Cx[dv]*phi.Node( i+ (int)lb.CxC[dv] , j + (int)lb.CyC[dv], k + (int)lb.CzC[dv]) ;
        grad_muy += lb.W[dv]*lb.Cy[dv]*phi.Node( i+ (int)lb.CxC[dv] , j + (int)lb.CyC[dv], k + (int)lb.CzC[dv]) ;
        grad_muz += lb.W[dv]*lb.Cz[dv]*phi.Node( i+ (int)lb.CxC[dv] , j + (int)lb.CyC[dv], k + (int)lb.CzC[dv]) ;
    }


    Force.Cell(i,j,k,0) = -mu.Cell(i,j,k)* Coeff_grad*(grad_mux);
    Force.Cell(i,j,k,1) = -mu.Cell(i,j,k)* Coeff_grad*(grad_muy);  
    Force.Cell(i,j,k,2) = -mu.Cell(i,j,k)* Coeff_grad*(grad_muz);  

    
    }



    //> gradient of Rho 
    //! was used in the free energy butnot needed
    // real grad_rhox = 0.0,grad_rhoy = 0.0,grad_rhoz = 0.0;
    // del_t = 1.0;
    // Coeff = (1.0/(del_t*lb.theta0));

    // for(int dv = 0; dv< grid.d_v; dv++){
    //     grad_rhox += lb.W[dv]*lb.Cx[dv]*rho.Node( i+ lb.Cx[dv] , j + lb.Cy[dv], k + lb.Cz[dv]) ;
    //     grad_rhoy += lb.W[dv]*lb.Cy[dv]*rho.Node( i+ lb.Cx[dv] , j + lb.Cy[dv], k + lb.Cz[dv]) ;
    //     grad_rhoz += lb.W[dv]*lb.Cz[dv]*rho.Node( i+ lb.Cx[dv] , j + lb.Cy[dv], k + lb.Cz[dv]) ;
    // }
    // real grad_rho = Coeff* (grad_rhox + grad_rhoy + grad_rhoz);





    //>laplacian of rho
    // for(int i = 0 + grid.noghost; i < grid.n_x_node - (grid.noghost) ; i++){
    //     for(int j = 0 + grid.noghost;j < grid.n_y_node - (grid.noghost) ; j++){
    //         for(int k = 0 + grid.noghost;k < grid.n_z_node - (grid.noghost) ; k++){

    //             //>node
    //             laplacian_pnidplusfnidbyrho.Node(i,j,k)  = 0.0;
    //             real del_t = 1.0;
    //             real Coeff = (2.0/(del_t*del_t*lb.theta0));

    //             for(int dv = 0; dv< grid.d_v; dv++)
    //                 laplacian_pnidplusfnidbyrho.Node(i,j,k) += lb.W[dv]*((pnid.Node( i+ lb.Cx[dv] , j + lb.Cy[dv], k + lb.Cz[dv]) + fnid.Node(i+ lb.Cx[dv] , j + lb.Cy[dv], k + lb.Cz[dv])) / rho.Node(i+ lb.Cx[dv] , j + lb.Cy[dv], k + lb.Cz[dv])) ;

    //             laplacian_pnidplusfnidbyrho.Node(i,j,k) = Coeff * ( laplacian_pnidplusfnidbyrho.Node(i,j,k) - ((pnid.Node(i,j,k) + fnid.Node(i,j,k)) / rho.Node(i,j,k)));

    //             //>cell
    //             laplacian_pnidplusfnidbyrho.Cell(i,j,k)  = 0.0;
                
    //             for(int dv = 0; dv< grid.d_v; dv++)
    //                 laplacian_pnidplusfnidbyrho.Cell(i,j,k) += lb.W[dv]*((pnid.Cell( i+ lb.Cx[dv] , j + lb.Cy[dv], k + lb.Cz[dv]) + fnid.Cell(i+ lb.Cx[dv] , j + lb.Cy[dv], k + lb.Cz[dv])) / rho.Cell(i+ lb.Cx[dv] , j + lb.Cy[dv], k + lb.Cz[dv])) ;

    //             laplacian_pnidplusfnidbyrho.Cell(i,j,k) = Coeff * ( laplacian_pnidplusfnidbyrho.Cell(i,j,k) - ((pnid.Cell(i,j,k) + fnid.Cell(i,j,k)) / rho.Cell(i,j,k)));



    //         }
    //     }
    // }

    // Periodic(laplacian_pnidplusfnidbyrho);





    template<typename T, typename T1>
void Multiphase_Force_eta_Cell(Grid_N_C_3D<T> &grid, Grid_N_C_3D<T> &rho, Grid_N_C_3D<T> &pnid, Grid_N_C_3D<T> &fnid, Grid_N_C_3D<T> &munid, Grid_N_C_3D<T> &laplacian_rho, Grid_N_C_3D<T> &laplacian_fnid, Grid_N_C_3D<T> &gradient_rho, 
            lbmD3Q35<T1> &lb, real &Fx, real &Fy, real &Fz, int i, int j,int k, real kappa, real eta){

    //> CHEMICAL POTENTIAL FORMULATION 
    real grad_mux = 0.0, grad_muy = 0.0, grad_muz = 0.0;
    Fx = 0, Fy = 0, Fz = 0.0;
    real del_t = 1.0;
    real Coeff_grad = (1.0/(del_t*lb.theta0));

    //#change the munid to bulk remove the laplacian rho when using this
    //
    for(int dv = 0; dv< 27; dv++){
        grad_mux += rho.Cell(i,j,k)*Coeff_grad*lb.W[dv]*lb.Cx[dv]*munid.Cell( i+ (int)lb.Cx[dv] , j + (int)lb.Cy[dv], k + (int)lb.Cz[dv]) ;
        grad_muy += rho.Cell(i,j,k)*Coeff_grad*lb.W[dv]*lb.Cy[dv]*munid.Cell( i+ (int)lb.Cx[dv] , j + (int)lb.Cy[dv], k + (int)lb.Cz[dv]) ;
        grad_muz += rho.Cell(i,j,k)*Coeff_grad*lb.W[dv]*lb.Cz[dv]*munid.Cell( i+ (int)lb.Cx[dv] , j + (int)lb.Cy[dv], k + (int)lb.Cz[dv]) ;
    }

    for(int dv = 27; dv<grid.d_v; dv++){
        grad_mux += rho.Cell(i,j,k)*Coeff_grad*lb.W[dv]*lb.Cx[dv]*munid.Node( i+ (int)lb.CxC[dv] , j + (int)lb.CyC[dv], k + (int)lb.CzC[dv]) ;
        grad_muy += rho.Cell(i,j,k)*Coeff_grad*lb.W[dv]*lb.Cy[dv]*munid.Node( i+ (int)lb.CxC[dv] , j + (int)lb.CyC[dv], k + (int)lb.CzC[dv]) ;
        grad_muz += rho.Cell(i,j,k)*Coeff_grad*lb.W[dv]*lb.Cz[dv]*munid.Node( i+ (int)lb.CxC[dv] , j + (int)lb.CyC[dv], k + (int)lb.CzC[dv]) ;
    }

    //munid*grad(rho)
        grad_mux += munid.Cell(i,j,k)*gradient_rho.Cell(i,j,k,0);
        grad_muy += munid.Cell(i,j,k)*gradient_rho.Cell(i,j,k,1);
        grad_muz += munid.Cell(i,j,k)*gradient_rho.Cell(i,j,k,2);

    //fnid
    for(int dv = 0; dv< 27; dv++){
        grad_mux = grad_mux - eta* ( Coeff_grad*lb.W[dv]*lb.Cx[dv]*fnid.Cell( i+ (int)lb.Cx[dv] , j + (int)lb.Cy[dv], k + (int)lb.Cz[dv])) ;
        grad_muy = grad_muy - eta* ( Coeff_grad*lb.W[dv]*lb.Cy[dv]*fnid.Cell( i+ (int)lb.Cx[dv] , j + (int)lb.Cy[dv], k + (int)lb.Cz[dv])) ;
        grad_muz = grad_muz - eta* ( Coeff_grad*lb.W[dv]*lb.Cz[dv]*fnid.Cell( i+ (int)lb.Cx[dv] , j + (int)lb.Cy[dv], k + (int)lb.Cz[dv])) ;
    }

    for(int dv = 27; dv<grid.d_v; dv++){
        grad_mux =  grad_mux -  eta *(Coeff_grad*lb.W[dv]*lb.Cx[dv]*fnid.Node( i+ (int)lb.CxC[dv] , j + (int)lb.CyC[dv], k + (int)lb.CzC[dv])) ;
        grad_muy =  grad_muy -  eta *(Coeff_grad*lb.W[dv]*lb.Cy[dv]*fnid.Node( i+ (int)lb.CxC[dv] , j + (int)lb.CyC[dv], k + (int)lb.CzC[dv])) ;
        grad_muz =  grad_muz -  eta *(Coeff_grad*lb.W[dv]*lb.Cz[dv]*fnid.Node( i+ (int)lb.CxC[dv] , j + (int)lb.CyC[dv], k + (int)lb.CzC[dv])) ;
    }


    //fnid
   // fourth oder corrections

    real coeff_grad_4th = 0.5*lb.theta0*del_t*del_t   ;
    for(int dv = 0; dv< 27; dv++){
        grad_mux = grad_mux - (1 - eta) * ((Coeff_grad*lb.W[dv]*lb.Cx[dv]*fnid.Cell( i+ (int)lb.Cx[dv] , j + (int)lb.Cy[dv], k + (int)lb.Cz[dv]) - coeff_grad_4th*Coeff_grad*lb.W[dv]*lb.Cx[dv]*laplacian_fnid.Cell( i+ (int)lb.Cx[dv] , j + (int)lb.Cy[dv], k + (int)lb.Cz[dv])));                    
        grad_muy = grad_muy - (1 - eta) * ((Coeff_grad*lb.W[dv]*lb.Cy[dv]*fnid.Cell( i+ (int)lb.Cx[dv] , j + (int)lb.Cy[dv], k + (int)lb.Cz[dv]) - coeff_grad_4th*Coeff_grad*lb.W[dv]*lb.Cy[dv]*laplacian_fnid.Cell( i+ (int)lb.Cx[dv] , j + (int)lb.Cy[dv], k + (int)lb.Cz[dv])));                    
        grad_muz = grad_muz - (1 - eta) * ((Coeff_grad*lb.W[dv]*lb.Cz[dv]*fnid.Cell( i+ (int)lb.Cx[dv] , j + (int)lb.Cy[dv], k + (int)lb.Cz[dv]) - coeff_grad_4th*Coeff_grad*lb.W[dv]*lb.Cz[dv]*laplacian_fnid.Cell( i+ (int)lb.Cx[dv] , j + (int)lb.Cy[dv], k + (int)lb.Cz[dv])));                    
    }

    for(int dv = 27; dv<grid.d_v; dv++){
        grad_mux = grad_mux - (1 - eta) * ((Coeff_grad*lb.W[dv]*lb.Cx[dv]*fnid.Node( i+ (int)lb.CxC[dv] , j + (int)lb.CyC[dv], k + (int)lb.CzC[dv]) - coeff_grad_4th*Coeff_grad*lb.W[dv]*lb.Cx[dv]*laplacian_fnid.Node( i+ (int)lb.CxC[dv] , j + (int)lb.CyC[dv], k + (int)lb.CzC[dv]))) ;
        grad_muy = grad_muy - (1 - eta) * ((Coeff_grad*lb.W[dv]*lb.Cy[dv]*fnid.Node( i+ (int)lb.CxC[dv] , j + (int)lb.CyC[dv], k + (int)lb.CzC[dv]) - coeff_grad_4th*Coeff_grad*lb.W[dv]*lb.Cy[dv]*laplacian_fnid.Node( i+ (int)lb.CxC[dv] , j + (int)lb.CyC[dv], k + (int)lb.CzC[dv]))) ;
        grad_muz = grad_muz - (1 - eta) * ((Coeff_grad*lb.W[dv]*lb.Cz[dv]*fnid.Node( i+ (int)lb.CxC[dv] , j + (int)lb.CyC[dv], k + (int)lb.CzC[dv]) - coeff_grad_4th*Coeff_grad*lb.W[dv]*lb.Cz[dv]*laplacian_fnid.Node( i+ (int)lb.CxC[dv] , j + (int)lb.CyC[dv], k + (int)lb.CzC[dv]))) ;
    }



    //
    for(int dv = 0; dv< 27; dv++){

        int ic = i+ (int)lb.Cx[dv]; int jc = j + (int)lb.Cy[dv];int kc = k + (int)lb.Cz[dv];

        grad_mux += kappa*Coeff_grad*lb.W[dv]*lb.Cx[dv]*( 0.5*pow(gradient_rho.Cell( ic ,jc, kc,X), 2 ) - 0.5*pow(gradient_rho.Cell( ic ,jc, kc,Y), 2 )  - 0.5*pow(gradient_rho.Cell( ic ,jc, kc,Z), 2 ) - rho.Cell(ic ,jc, kc) *laplacian_rho.Cell(ic ,jc, kc)) 
                 +  kappa*Coeff_grad*lb.W[dv]*lb.Cy[dv]*( gradient_rho.Cell( ic ,jc, kc,X) * gradient_rho.Cell( ic ,jc, kc,Y)) 
                 +  kappa*Coeff_grad*lb.W[dv]*lb.Cz[dv]*( gradient_rho.Cell( ic ,jc, kc,X) * gradient_rho.Cell( ic ,jc, kc,Z)) ;



        grad_muy += kappa*Coeff_grad*lb.W[dv]*lb.Cx[dv]*( gradient_rho.Cell( ic ,jc, kc,X) * gradient_rho.Cell( ic ,jc, kc,Y)) 
                 +  kappa*Coeff_grad*lb.W[dv]*lb.Cy[dv]*( 0.5*pow(gradient_rho.Cell( ic ,jc, kc,Y), 2 ) - 0.5*pow(gradient_rho.Cell( ic ,jc, kc,X), 2 )  - 0.5*pow(gradient_rho.Cell( ic ,jc, kc,Z), 2 ) - rho.Cell(ic ,jc, kc) *laplacian_rho.Cell(ic ,jc, kc)) 
                 +  kappa*Coeff_grad*lb.W[dv]*lb.Cz[dv]*( gradient_rho.Cell( ic ,jc, kc,Y) * gradient_rho.Cell( ic ,jc, kc,Z)) ;


        grad_muz += kappa*Coeff_grad*lb.W[dv]*lb.Cx[dv]*( gradient_rho.Cell( ic ,jc, kc,Z) * gradient_rho.Cell( ic ,jc, kc,X)) 
                 +  kappa*Coeff_grad*lb.W[dv]*lb.Cy[dv]*( gradient_rho.Cell( ic ,jc, kc,Z) * gradient_rho.Cell( ic ,jc, kc,Y)) 
                 +  kappa*Coeff_grad*lb.W[dv]*lb.Cz[dv]*( 0.5*pow(gradient_rho.Cell( ic ,jc, kc,Z), 2 ) - 0.5*pow(gradient_rho.Cell( ic ,jc, kc,X), 2 )  - 0.5*pow(gradient_rho.Cell( ic ,jc, kc,Y), 2 ) - rho.Cell(ic ,jc, kc) *laplacian_rho.Cell(ic ,jc, kc)) ;




    }

    for(int dv = 27; dv<grid.d_v; dv++){
 
        int ic = i+ (int)lb.CxC[dv]; int jc = j + (int)lb.CyC[dv];int kc = k + (int)lb.CzC[dv];


        grad_mux += kappa*Coeff_grad*lb.W[dv]*lb.Cx[dv]*( 0.5*pow(gradient_rho.Node( ic ,jc, kc,X), 2 ) - 0.5*pow(gradient_rho.Node( ic ,jc, kc,Y), 2 )  - 0.5*pow(gradient_rho.Node( ic ,jc, kc,Z), 2 ) - rho.Node(ic ,jc, kc) *laplacian_rho.Node(ic ,jc, kc)) 
                 +  kappa*Coeff_grad*lb.W[dv]*lb.Cy[dv]*( gradient_rho.Node( ic ,jc, kc,X) * gradient_rho.Node( ic ,jc, kc,Y)) 
                 +  kappa*Coeff_grad*lb.W[dv]*lb.Cz[dv]*( gradient_rho.Node( ic ,jc, kc,X) * gradient_rho.Node( ic ,jc, kc,Z)) ;



        grad_muy += kappa*Coeff_grad*lb.W[dv]*lb.Cx[dv]*( gradient_rho.Node( ic ,jc, kc,X) * gradient_rho.Node( ic ,jc, kc,Y)) 
                 +  kappa*Coeff_grad*lb.W[dv]*lb.Cy[dv]*( 0.5*pow(gradient_rho.Node( ic ,jc, kc,Y), 2 ) - 0.5*pow(gradient_rho.Node( ic ,jc, kc,X), 2 )  - 0.5*pow(gradient_rho.Node( ic ,jc, kc,Z), 2 ) - rho.Node(ic ,jc, kc) *laplacian_rho.Node(ic ,jc, kc)) 
                 +  kappa*Coeff_grad*lb.W[dv]*lb.Cz[dv]*( gradient_rho.Node( ic ,jc, kc,Y) * gradient_rho.Node( ic ,jc, kc,Z)) ;


        grad_muz += kappa*Coeff_grad*lb.W[dv]*lb.Cx[dv]*( gradient_rho.Node( ic ,jc, kc,Z) * gradient_rho.Node( ic ,jc, kc,X)) 
                 +  kappa*Coeff_grad*lb.W[dv]*lb.Cy[dv]*( gradient_rho.Node( ic ,jc, kc,Z) * gradient_rho.Node( ic ,jc, kc,Y)) 
                 +  kappa*Coeff_grad*lb.W[dv]*lb.Cz[dv]*( 0.5*pow(gradient_rho.Node( ic ,jc, kc,Z), 2 ) - 0.5*pow(gradient_rho.Node( ic ,jc, kc,X), 2 )  - 0.5*pow(gradient_rho.Node( ic ,jc, kc,Y), 2 ) - rho.Node(ic ,jc, kc) *laplacian_rho.Node(ic ,jc, kc)) ;


    }





    Fx = - grad_mux/rho.Cell(i,j,k);
    Fy = - grad_muy/rho.Cell(i,j,k);  
    Fz = - grad_muz/rho.Cell(i,j,k);  


    


}
