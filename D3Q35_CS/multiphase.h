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
void Multiphase_terms(Grid_N_C_3D<T> &grid, Grid_N_C_3D<T> &Force, Grid_N_C_3D<T> &rho, Grid_N_C_3D<T> &pnid, Grid_N_C_3D<T> &fnid, Grid_N_C_3D<T> &munid, Grid_N_C_3D<T> &laplacian_rho,  Grid_N_C_3D<T> &laplacian_fnid, Grid_N_C_3D<T> &gradient_rho, 
            lbmD3Q35<T1> &lb35,lbmD3Q15<T1> &lb15, real TbyTc, real kappa, real a , real b, real dx, real dt ){

    real ux = 0, uy = 0, uz = 0;



    for(int i = 0 + grid.noghost; i < grid.n_x_node - (grid.noghost) ; i++){
        for(int j = 0 + grid.noghost;j < grid.n_y_node - (grid.noghost) ; j++){
            for(int k = 0 + grid.noghost;k < grid.n_z_node - (grid.noghost) ; k++){
            
                get_moments_Node(grid, lb35,  ux, uy, uz, rho.Node(i,j,k), i, j ,k, Force, dx, dt); 

                get_moments_Cell(grid, lb35,  ux, uy, uz, rho.Cell(i,j,k), i, j ,k, Force, dx , dt); 

                real eta = rho.Node(i,j,k)*b/4.0;

                //> pnid Node
                pnid.Node(i,j,k) = (rho.Node(i,j,k)*lb35.theta0*(1.0 + eta + eta* eta - eta*eta*eta) )/pow(1.0 - eta,  3.0)  - 
                                    a * rho.Node(i,j,k)*rho.Node(i,j,k) ;

                //> Fnid_node
                fnid.Node(i,j,k) = -1.0*(rho.Node(i,j,k)*lb35.theta0*(3.0*eta*eta - 4.0 * eta) )/pow(1.0 - eta,  2.0)  - 
                                    a * rho.Node(i,j,k)*rho.Node(i,j,k) ;

                eta  = rho.Cell(i,j,k)*b/4.0;
                //> pnid Cell
                pnid.Cell(i,j,k) = (rho.Cell(i,j,k)*lb35.theta0*(1.0 + eta + eta* eta - eta*eta*eta) )/pow(1.0 - eta,  3.0)  - 
                                    a * rho.Cell(i,j,k)*rho.Cell(i,j,k) ;

                //>Fnid_Cell
                fnid.Cell(i,j,k) =  -1.0*(rho.Cell(i,j,k)*lb35.theta0*(3.0*eta*eta - 4.0 * eta) )/pow(1.0 - eta,  2.0)  - 
                                    a * rho.Cell(i,j,k)*rho.Cell(i,j,k) ;

            }
        }
    }



     Periodic_x(rho);
     Periodic_y(rho);
     Periodic_z(rho);


    // grad_zero_x(rho);
    // grad_zero_y(rho);
    // grad_zero_z(rho);



     Periodic_x(pnid);
     Periodic_y(pnid);
     Periodic_z(pnid);

    // grad_zero_x(pnid);
    // grad_zero_y(pnid);
    // grad_zero_z(pnid);




     Periodic_x(fnid);
     Periodic_y(fnid);
     Periodic_z(fnid);
    
    // grad_zero_x(fnid);
    // grad_zero_y(fnid);
    // grad_zero_z(fnid);





    for(int i = 0 + grid.noghost; i < grid.n_x_node - (grid.noghost) ; i++){
        for(int j = 0 + grid.noghost;j < grid.n_y_node - (grid.noghost) ; j++){
            for(int k = 0 + grid.noghost;k < grid.n_z_node - (grid.noghost) ; k++){



                //> laplacian of Rho  
                real Coeff = (2.0/(dt*dt*lb35.theta0));
                
                laplacian_rho.Node(i,j,k)  = 0.0;

                for(int dv = 0; dv< 27; dv++)
                    laplacian_rho.Node(i,j,k) += lb35.W[dv]*rho.Node( i+ (int)lb35.Cx[dv]  , j + (int)lb35.Cy[dv] , k + (int)lb35.Cz[dv] ) ;

                for(int dv = 27; dv< 35;dv++)
                    laplacian_rho.Node(i,j,k) += lb35.W[dv]*rho.Cell( i+ (int)lb35.CxF[dv] , j + (int)lb35.CyF[dv], k + (int)lb35.CzF[dv]) ;

                laplacian_rho.Node(i,j,k) = Coeff * ( laplacian_rho.Node(i,j,k) - rho.Node(i,j,k));


                laplacian_rho.Cell(i,j,k)  = 0.0;
                
                for(int dv = 0; dv< 27; dv++)
                    laplacian_rho.Cell(i,j,k) += lb35.W[dv]*rho.Cell( i+ (int)lb35.Cx[dv]  , j + (int)lb35.Cy[dv] , k + (int)lb35.Cz[dv] ) ;

                for(int dv = 27; dv< 35;dv++)
                    laplacian_rho.Cell(i,j,k) += lb35.W[dv]*rho.Node( i+ (int)lb35.CxC[dv] , j + (int)lb35.CyC[dv], k + (int)lb35.CzC[dv]) ;

                laplacian_rho.Cell(i,j,k) = Coeff * ( laplacian_rho.Cell(i,j,k) - rho.Cell(i,j,k));
                //>---------------------




                // //>------
                real kappa_node = kappa;

                // kappa_node = kappa - (1.0/2.0)*dx*dx* lb35.theta0* (
                //                                         -32.0 *b* lb35.theta0*(-16.0 + b* rho.Node(i,j,k)) / (pow(-4.0 + b*rho.Node(i,j,k) , 4.0))
                //                                             - 2.0*a
                //                                         ) 
                //                                     ;

                double eta = rho.Node(i,j,k)*b/4.0;

                //> munid
                munid.Node(i,j,k) = -2.0*rho.Node(i,j,k)*a;
                munid.Node(i,j,k) += lb35.theta0*(3.0*eta*eta*eta - 9.0*eta*eta + 8.0*eta) / pow(1.0 - eta, 3.0);
                munid.Node(i,j,k) -= kappa_node*laplacian_rho.Node(i,j,k);





                real kappa_cell = kappa;
                // kappa_cell = kappa - (1.0/2.0)*dx*dx* lb35.theta0* (
                //                                        -32.0 *b* lb35.theta0*(-16.0 + b* rho.Cell(i,j,k)) / (pow(-4.0 + b*rho.Cell(i,j,k) , 4.0))
                //                                             - 2.0*a
                //                                         ) 
                //                                     ;


                eta = rho.Cell(i,j,k)*b/4.0;

                munid.Cell(i,j,k) = -2.0*rho.Cell(i,j,k)*a;
                munid.Cell(i,j,k) += lb35.theta0*(3.0*eta*eta*eta - 9.0*eta*eta + 8.0*eta) / pow(1.0 - eta, 3.0);
                munid.Cell(i,j,k) -= kappa_cell*laplacian_rho.Cell(i,j,k);

            }
        }    
    }



    Periodic_x(laplacian_rho);
    Periodic_y(laplacian_rho);
    Periodic_z(laplacian_rho);

    // grad_zero_x(laplacian_rho);
    // grad_zero_y(laplacian_rho);
    // grad_zero_z(laplacian_rho);


    Periodic_x(laplacian_fnid);
    Periodic_y(laplacian_fnid);
    Periodic_z(laplacian_fnid);

    // grad_zero_x(laplacian_fnid);
    // grad_zero_y(laplacian_fnid);
    // grad_zero_z(laplacian_fnid);

    Periodic_x(munid);
    Periodic_y(munid);
    Periodic_z(munid);

    // grad_zero_x(munid);
    // grad_zero_y(munid);
    // grad_zero_z(munid);


}





template<typename T, typename T1>
void Multiphase_Force_Node(Grid_N_C_3D<T> &grid, Grid_N_C_3D<T> &rho, Grid_N_C_3D<T> &pnid, Grid_N_C_3D<T> &fnid, Grid_N_C_3D<T> &munid, Grid_N_C_3D<T> &laplacian_rho,
            lbmD3Q35<T1> &lb35,lbmD3Q15<T1> &lb15, Grid_N_C_3D<T> &Force, int i, int j,int k, real kappa, real a , real b, real g, real dx, real dt){

    

    //> CHEMICAL POTENTIAL FORMULATION 

    Force.Node(i,j,k,0) = 0.0;
    Force.Node(i,j,k,1) = 0.0;
    Force.Node(i,j,k,2) = 0.0;


    real Coeff_grad = (1.0/(dt*lb35.theta0));

    for(int dv = 0; dv< 27; dv++){
        Force.Node(i,j,k,0)+= lb35.W[dv]*lb35.Cx[dv]*munid.Node( i+ (int)lb35.Cx[dv] , j + (int)lb35.Cy[dv], k + (int)lb35.Cz[dv]) ;
        Force.Node(i,j,k,1)+= lb35.W[dv]*lb35.Cy[dv]*munid.Node( i+ (int)lb35.Cx[dv] , j + (int)lb35.Cy[dv], k + (int)lb35.Cz[dv]) ;
        Force.Node(i,j,k,2)+= lb35.W[dv]*lb35.Cz[dv]*munid.Node( i+ (int)lb35.Cx[dv] , j + (int)lb35.Cy[dv], k + (int)lb35.Cz[dv]) ;
    }

    for(int dv = 27; dv<grid.d_v; dv++){
        Force.Node(i,j,k,0) += lb35.W[dv]*lb35.Cx[dv]*munid.Cell( i+ (int)lb35.CxF[dv] , j + (int)lb35.CyF[dv], k + (int)lb35.CzF[dv]) ;
        Force.Node(i,j,k,1) += lb35.W[dv]*lb35.Cy[dv]*munid.Cell( i+ (int)lb35.CxF[dv] , j + (int)lb35.CyF[dv], k + (int)lb35.CzF[dv]) ;
        Force.Node(i,j,k,2) += lb35.W[dv]*lb35.Cz[dv]*munid.Cell( i+ (int)lb35.CxF[dv] , j + (int)lb35.CyF[dv], k + (int)lb35.CzF[dv]) ;
    }

    Force.Node(i,j,k,0) = - Coeff_grad*(Force.Node(i,j,k,0))        ;
    Force.Node(i,j,k,1) = - Coeff_grad*(Force.Node(i,j,k,1))  -g    ;
    Force.Node(i,j,k,2) = - Coeff_grad*(Force.Node(i,j,k,2))        ;

    // // // > direct              

    // Force.Node(i,j,k,0) = 0.0;
    // Force.Node(i,j,k,1) = 0.0;
    // Force.Node(i,j,k,2) = 0.0;


    // real del_t = 1.0;
    // real Coeff_grad = (1.0/(del_t*lb35.theta0));


    // real kappa_node = kappa;

    //             // kappa_node = kappa - (1.0/6.0)*  (
    //             //                                        -32.0 *b* lb35.theta0*(-16.0 + b* rho.Node(i,j,k)) / (pow(-4.0 + b*rho.Node(i,j,k) , 4.0))
    //             //                                             - 2.0*a
    //             //                                         ) 
    //             //                                     ;





    // for(int dv = 0; dv < 27; dv++){

    //     Force.Node(i,j,k,0) += Coeff_grad*(               
    //                         lb35.W[dv]*lb35.Cx[dv]*( (pnid.Node( i+ (int)lb35.Cx[dv] , j + (int)lb35.Cy[dv], k + (int)lb35.Cz[dv]) + fnid.Node( i+ (int)lb35.Cx[dv] , j + (int)lb35.Cy[dv], k + (int)lb35.Cz[dv]))/ rho.Node( i+ (int)lb35.Cx[dv] ,j + (int)lb35.Cy[dv], k + (int)lb35.Cz[dv])  )
    //                         )
    //                         - kappa_node*Coeff_grad*( lb35.W[dv]*lb35.Cx[dv]*( laplacian_rho.Node( i+ (int)lb35.Cx[dv] , j + (int)lb35.Cy[dv], k + (int)lb35.Cz[dv])))
    //                         ;

    //     Force.Node(i,j,k,1) += Coeff_grad*(               
    //                         lb35.W[dv]*lb35.Cy[dv]*( (pnid.Node( i+ (int)lb35.Cx[dv] , j + (int)lb35.Cy[dv], k + (int)lb35.Cz[dv]) + fnid.Node( i+ (int)lb35.Cx[dv] , j + (int)lb35.Cy[dv], k + (int)lb35.Cz[dv]))/ rho.Node( i+ (int)lb35.Cx[dv] ,j + (int)lb35.Cy[dv], k + (int)lb35.Cz[dv])  )
    //                         )
    //                         - kappa_node*Coeff_grad*( lb35.W[dv]*lb35.Cy[dv]*( laplacian_rho.Node( i+ (int)lb35.Cx[dv] , j + (int)lb35.Cy[dv], k + (int)lb35.Cz[dv])))
    //                         ;

    //     Force.Node(i,j,k,2) += Coeff_grad*(               
    //                         lb35.W[dv]*lb35.Cz[dv]*( (pnid.Node( i+ (int)lb35.Cx[dv] , j + (int)lb35.Cy[dv], k + (int)lb35.Cz[dv]) + fnid.Node( i+ (int)lb35.Cx[dv] , j + (int)lb35.Cy[dv], k + (int)lb35.Cz[dv]))/ rho.Node( i+ (int)lb35.Cx[dv] ,j + (int)lb35.Cy[dv], k + (int)lb35.Cz[dv])  )
    //                         )
    //                         - kappa_node*Coeff_grad*( lb35.W[dv]*lb35.Cz[dv]*( laplacian_rho.Node( i+ (int)lb35.Cx[dv] , j + (int)lb35.Cy[dv], k + (int)lb35.Cz[dv])))
    //                         ;

    // }

    // for(int dv = 27; dv < 35; dv++){

    //     Force.Node(i,j,k,0) += Coeff_grad*(               
    //                         lb35.W[dv]*lb35.Cx[dv]*( (pnid.Cell( i+ (int)lb35.CxF[dv] , j + (int)lb35.CyF[dv], k + (int)lb35.CzF[dv]) + fnid.Cell( i+ (int)lb35.CxF[dv] , j + (int)lb35.CyF[dv], k + (int)lb35.CzF[dv]))/ rho.Cell( i+ (int)lb35.CxF[dv] ,j + (int)lb35.CyF[dv], k + (int)lb35.CzF[dv])  )
    //                         )
    //                         - kappa_node*Coeff_grad*( lb35.W[dv]*lb35.Cx[dv]*( laplacian_rho.Cell( i+ (int)lb35.CxF[dv] , j + (int)lb35.CyF[dv], k + (int)lb35.CzF[dv])))
    //                         ;

    //     Force.Node(i,j,k,1) += Coeff_grad*(               
    //                         lb35.W[dv]*lb35.Cy[dv]*( (pnid.Cell( i+ (int)lb35.CxF[dv] , j + (int)lb35.CyF[dv], k + (int)lb35.CzF[dv]) + fnid.Cell( i+ (int)lb35.CxF[dv] , j + (int)lb35.CyF[dv], k + (int)lb35.CzF[dv]))/ rho.Cell( i+ (int)lb35.CxF[dv] ,j + (int)lb35.CyF[dv], k + (int)lb35.CzF[dv])  )
    //                         )
    //                         - kappa_node*Coeff_grad*( lb35.W[dv]*lb35.Cy[dv]*( laplacian_rho.Cell( i+ (int)lb35.CxF[dv] , j + (int)lb35.CyF[dv], k + (int)lb35.CzF[dv])))
    //                         ;

    //     Force.Node(i,j,k,2) += Coeff_grad*(               
    //                         lb35.W[dv]*lb35.Cz[dv]*( (pnid.Cell( i+ (int)lb35.CxF[dv] , j + (int)lb35.CyF[dv], k + (int)lb35.CzF[dv]) + fnid.Cell( i+ (int)lb35.CxF[dv] , j + (int)lb35.CyF[dv], k + (int)lb35.CzF[dv]))/ rho.Cell( i+ (int)lb35.CxF[dv] ,j + (int)lb35.CyF[dv], k + (int)lb35.CzF[dv])  )
    //                         )
    //                         - kappa_node*Coeff_grad*( lb35.W[dv]*lb35.Cz[dv]*( laplacian_rho.Cell( i+ (int)lb35.CxF[dv] , j + (int)lb35.CyF[dv], k + (int)lb35.CzF[dv])))
    //                         ;

    // }


    
    
    // // // 4th order corrections
    // // real coeff_4th_grad = (lb35.theta0*del_t*del_t)/2.0;

    // // for(int dv = 0; dv < grid.d_v; dv++){
    // //     Force.Node(i,j,k,0) -= coeff_4th_grad* Coeff_grad*(lb35.W[dv]*lb35.Cx[dv]*laplacian_pnidplusfnidbyrho.Node( i+ lb35.Cx[dv] , j + lb35.Cy[dv],k + lb35.Cz[dv])) ;
    // //     Force.Node(i,j,k,1) -= coeff_4th_grad* Coeff_grad*(lb35.W[dv]*lb35.Cy[dv]*laplacian_pnidplusfnidbyrho.Node( i+ lb35.Cx[dv] , j + lb35.Cy[dv],k + lb35.Cz[dv])) ;
    // //     Force.Node(i,j,k,2) -= coeff_4th_grad* Coeff_grad*(lb35.W[dv]*lb35.Cz[dv]*laplacian_pnidplusfnidbyrho.Node( i+ lb35.Cx[dv] , j + lb35.Cy[dv],k + lb35.Cz[dv])) ;
    // // }

    // Force.Node(i,j,k,0) = -1.0*Force.Node(i,j,k,0);
    // Force.Node(i,j,k,1) = -1.0*Force.Node(i,j,k,1);
    // Force.Node(i,j,k,2) = -1.0*Force.Node(i,j,k,2);
    // //>--------------------------


}


template<typename T, typename T1>
void Multiphase_Force_Cell(Grid_N_C_3D<T> &grid, Grid_N_C_3D<T> &rho, Grid_N_C_3D<T> &pnid, Grid_N_C_3D<T> &fnid, Grid_N_C_3D<T> &munid, Grid_N_C_3D<T> &laplacian_rho,
            lbmD3Q35<T1> &lb35, lbmD3Q15<T1> &lb15, Grid_N_C_3D<T> &Force, int i, int j,int k, real kappa, real a, real b, real g, real dx ,real dt){

    //> CHEMICAL POTENTIAL FORMULATION 
    
    Force.Cell(i,j,k,0) = 0.0;
    Force.Cell(i,j,k,1) = 0.0;
    Force.Cell(i,j,k,2) = 0.0;

    
    //> CHEMICAL POTENTIAL FORMULATION 
    real Coeff_grad = (1.0/(dt*lb35.theta0));

    for(int dv = 0; dv< 27; dv++){
        Force.Cell(i,j,k,0)  += lb35.W[dv]*lb35.Cx[dv]*munid.Cell( i+ (int)lb35.Cx[dv] , j + (int)lb35.Cy[dv], k + (int)lb35.Cz[dv]) ;
        Force.Cell(i,j,k,1)  += lb35.W[dv]*lb35.Cy[dv]*munid.Cell( i+ (int)lb35.Cx[dv] , j + (int)lb35.Cy[dv], k + (int)lb35.Cz[dv]) ;
        Force.Cell(i,j,k,2)  += lb35.W[dv]*lb35.Cz[dv]*munid.Cell( i+ (int)lb35.Cx[dv] , j + (int)lb35.Cy[dv], k + (int)lb35.Cz[dv]) ;
    }

    for(int dv = 27; dv< 35; dv++){
        Force.Cell(i,j,k,0)  += lb35.W[dv]*lb35.Cx[dv]*munid.Node( i+ (int)lb35.CxC[dv] , j + (int)lb35.CyC[dv], k + (int)lb35.CzC[dv]) ;
        Force.Cell(i,j,k,1)  += lb35.W[dv]*lb35.Cy[dv]*munid.Node( i+ (int)lb35.CxC[dv] , j + (int)lb35.CyC[dv], k + (int)lb35.CzC[dv]) ;
        Force.Cell(i,j,k,2)  += lb35.W[dv]*lb35.Cz[dv]*munid.Node( i+ (int)lb35.CxC[dv] , j + (int)lb35.CyC[dv], k + (int)lb35.CzC[dv]) ;
    }


    Force.Cell(i,j,k,0) = - Coeff_grad*(Force.Cell(i,j,k,0) )       ;
    Force.Cell(i,j,k,1) = - Coeff_grad*(Force.Cell(i,j,k,1) ) - g   ;  
    Force.Cell(i,j,k,2) = - Coeff_grad*(Force.Cell(i,j,k,2) )       ;  



//     //> MECHANICAL FORMULATION
//     // // > direct              
//     Force.Cell(i,j,k,0) = 0.0;
//     Force.Cell(i,j,k,1) = 0.0;
//     Force.Cell(i,j,k,2) = 0.0;


//     del_t = 1.0;
//     Coeff_grad = (1.0/(del_t*lb35.theta0));


//     real kappa_cell = kappa;

//                 // kappa_cell = kappa - (1.0/6.0)*  (
//                 //                                        -32.0 *b* lb35.theta0*(-16.0 + b* rho.Cell(i,j,k)) / (pow(-4.0 + b*rho.Cell(i,j,k) , 4.0))
//                 //                                             - 2.0*a
//                 //                                         ) 
//                 //                                     ;





//     for(int dv = 0; dv < 27; dv++){

//             Force.Cell(i,j,k,0) += Coeff_grad*(               
//                             lb35.W[dv]*lb35.Cx[dv]*( (pnid.Cell( i+ (int)lb35.Cx[dv] , j + (int)lb35.Cy[dv], k + (int)lb35.Cz[dv]) + fnid.Cell( i+ (int)lb35.Cx[dv] , j + (int)lb35.Cy[dv], k + (int)lb35.Cz[dv]))/ rho.Cell( i+ (int)lb35.Cx[dv] ,j + (int)lb35.Cy[dv], k + (int)lb35.Cz[dv])  )
//                             )
//                             - kappa_cell*Coeff_grad*( lb35.W[dv]*lb35.Cx[dv]*( laplacian_rho.Cell( i+ (int)lb35.Cx[dv] , j + (int)lb35.Cy[dv], k + (int)lb35.Cz[dv])))
//                             ;

//             Force.Cell(i,j,k,1) += Coeff_grad*(               
//                             lb35.W[dv]*lb35.Cy[dv]*( (pnid.Cell( i+ (int)lb35.Cx[dv] , j + (int)lb35.Cy[dv], k + (int)lb35.Cz[dv]) + fnid.Cell( i+ (int)lb35.Cx[dv] , j + (int)lb35.Cy[dv], k + (int)lb35.Cz[dv]))/ rho.Cell( i+ (int)lb35.Cx[dv] ,j + (int)lb35.Cy[dv], k + (int)lb35.Cz[dv])  )
//                             )
//                             - kappa_cell*Coeff_grad*( lb35.W[dv]*lb35.Cy[dv]*( laplacian_rho.Cell( i+ (int)lb35.Cx[dv] , j + (int)lb35.Cy[dv], k + (int)lb35.Cz[dv])))
//                             ;

//             Force.Cell(i,j,k,2) += Coeff_grad*(               
//                             lb35.W[dv]*lb35.Cz[dv]*( (pnid.Cell( i+ (int)lb35.Cx[dv] , j + (int)lb35.Cy[dv], k + (int)lb35.Cz[dv]) + fnid.Cell( i+ (int)lb35.Cx[dv] , j + (int)lb35.Cy[dv], k + (int)lb35.Cz[dv]))/ rho.Cell( i+ (int)lb35.Cx[dv] ,j + (int)lb35.Cy[dv], k + (int)lb35.Cz[dv])  )
//                             )
//                             - kappa_cell*Coeff_grad*( lb35.W[dv]*lb35.Cz[dv]*( laplacian_rho.Cell( i+ (int)lb35.Cx[dv] , j + (int)lb35.Cy[dv], k + (int)lb35.Cz[dv])))
//                             ;

//     }
    




//     for(int dv = 27; dv < 35; dv++){

//             Force.Cell(i,j,k,0) += Coeff_grad*(               
//                             lb35.W[dv]*lb35.Cx[dv]*( (pnid.Node( i+ (int)lb35.CxC[dv] , j + (int)lb35.CyC[dv], k + (int)lb35.CzC[dv]) + fnid.Node( i+ (int)lb35.CxC[dv] , j + (int)lb35.CyC[dv], k + (int)lb35.CzC[dv]))/ rho.Node( i+ (int)lb35.CxC[dv] ,j + (int)lb35.CyC[dv], k + (int)lb35.CzC[dv])  )
//                             )
//                             - kappa_cell*Coeff_grad*( lb35.W[dv]*lb35.Cx[dv]*( laplacian_rho.Node( i+ (int)lb35.CxC[dv] , j + (int)lb35.CyC[dv], k + (int)lb35.CzC[dv])))
//                             ;

//             Force.Cell(i,j,k,1) += Coeff_grad*(               
//                             lb35.W[dv]*lb35.Cy[dv]*( (pnid.Node( i+ (int)lb35.CxC[dv] , j + (int)lb35.CyC[dv], k + (int)lb35.CzC[dv]) + fnid.Node( i+ (int)lb35.CxC[dv] , j + (int)lb35.CyC[dv], k + (int)lb35.CzC[dv]))/ rho.Node( i+ (int)lb35.CxC[dv] ,j + (int)lb35.CyC[dv], k + (int)lb35.CzC[dv])  )
//                             )
//                             - kappa_cell*Coeff_grad*( lb35.W[dv]*lb35.Cy[dv]*( laplacian_rho.Node( i+ (int)lb35.CxC[dv] , j + (int)lb35.CyC[dv], k + (int)lb35.CzC[dv])))
//                             ;

//             Force.Cell(i,j,k,2) += Coeff_grad*(               
//                             lb35.W[dv]*lb35.Cz[dv]*( (pnid.Node( i+ (int)lb35.CxC[dv] , j + (int)lb35.CyC[dv], k + (int)lb35.CzC[dv]) + fnid.Node( i+ (int)lb35.CxC[dv] , j + (int)lb35.CyC[dv], k + (int)lb35.CzC[dv]))/ rho.Node( i+ (int)lb35.CxC[dv] ,j + (int)lb35.CyC[dv], k + (int)lb35.CzC[dv])  )
//                             )
//                             - kappa_cell*Coeff_grad*( lb35.W[dv]*lb35.Cz[dv]*( laplacian_rho.Node( i+ (int)lb35.CxC[dv] , j + (int)lb35.CyC[dv], k + (int)lb35.CzC[dv])))
//                             ;

//     }
    


//     // // //> 4th order corrections
//     // //  coeff_4th_grad = (lb35.theta0*del_t*del_t)/2.0;

//     // // for(int dv = 0; dv < grid.d_v; dv++){
//     // //     Force.Cell(i,j,k,0) -= coeff_4th_grad* Coeff_grad*(lb35.W[dv]*lb35.Cx[dv]*laplacian_pnidplusfnidbyrho.Cell( i+ lb35.Cx[dv] , j + lb35.Cy[dv],k + lb35.Cz[dv])) ;
//     // //     Force.Cell(i,j,k,1) -= coeff_4th_grad* Coeff_grad*(lb35.W[dv]*lb35.Cy[dv]*laplacian_pnidplusfnidbyrho.Cell( i+ lb35.Cx[dv] , j + lb35.Cy[dv],k + lb35.Cz[dv])) ;
//     // //     Force.Cell(i,j,k,2) -= coeff_4th_grad* Coeff_grad*(lb35.W[dv]*lb35.Cz[dv]*laplacian_pnidplusfnidbyrho.Cell( i+ lb35.Cx[dv] , j + lb35.Cy[dv],k + lb35.Cz[dv])) ;
//     // // }

//    Force.Cell(i,j,k,0) = -1.0*Force.Cell(i,j,k,0);
//    Force.Cell(i,j,k,1) = -1.0*Force.Cell(i,j,k,1);
//    Force.Cell(i,j,k,2) = -1.0*Force.Cell(i,j,k,2);
//     //> ------------------------------

    
    
    }














//%                 this is from the 41 paper formulation
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
