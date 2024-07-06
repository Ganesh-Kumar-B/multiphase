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
void collide(Grid_N_C_2D<T> &gridf,Grid_N_C_2D<T> &gridg, Grid_N_C_2D<T> &grad_psi_rho,
            lbmD2Q9<T1> &lb9,double beta,double tau, double kappa, double g ,
            double phi_l, double phi_h,double rho_l, double rho_h,  double a , double b, Grid_N_C_2D<T> &Force){


    double  feq_Node[9] = {0},geq_Node[9]={0},gamma_Node[9]={0}, 
            ux = 0, uy = 0, p = 0;


    Grid_N_C_2D<T> phi(gridf.n_x,gridf.n_y,1,1);
    Grid_N_C_2D<T> rho(gridf.n_x,gridf.n_y,1,1);

    Grid_N_C_2D<T> laplacian_rho(gridf.n_x,gridf.n_y,1,1);

    Grid_N_C_2D<T> grad_psi_phi(gridf.n_x,gridf.n_y,1,2); // $phi 

    Grid_N_C_2D<T> psi_rho(gridf.n_x,gridf.n_y,1,1); // $psi 
    Grid_N_C_2D<T> psi_phi(gridf.n_x,gridf.n_y,1,1); // $phi 


    Multiphase_terms(gridf,gridg,lb9 ,rho, phi, laplacian_rho,psi_phi,psi_rho, grad_psi_phi, grad_psi_rho, rho_l, rho_h, phi_l, phi_h, a, b);

    //  // first the population of nodes are resetted and second  the population of the cells are resetted
    for(int i = 0 + gridf.noghost; i < gridf.n_x_node - (gridf.noghost) ; i++){
        for(int j = 0 + gridf.noghost;j < gridf.n_y_node - (gridf.noghost) ; j++){
            



            
            get_Force_and_gravity(gridf, lb9, laplacian_rho, kappa, g, Force,i,j);

            
            get_P_and_u(gridg, grad_psi_rho,lb9,  ux, uy,p,rho.Node(i,j), i, j, Force);            //for the node
            
            
            get_equi_f      (feq_Node   ,lb9, ux, uy, phi.Node(i,j)  );   //need phi, u 
            get_equi_g      (geq_Node   ,lb9, ux, uy, rho.Node(i,j),p);                               // need p and rho and u
            get_equi_gamma  (gamma_Node ,lb9, ux, uy, phi.Node(i,j)  );                           //need u onlu


            for (int dv = 0; dv< gridf.d_v; dv++){      
                gridf.Node(i,j,dv) =  gridf.Node(i,j,dv) + (1.0/tau)*(feq_Node[dv] - gridf.Node(i,j,dv)) 
                                    - (1 - 1.0/(2.0*tau))*lb9.thetaInverse *(
                                        gamma_Node[dv]*((lb9.Cx[dv] - ux)*grad_psi_phi.Node(i,j,0) + (lb9.Cy[dv] - uy)*grad_psi_phi.Node(i,j,1)  )
                                        )    ;
            }



            for (int dv = 0; dv< gridg.d_v; dv++){      
                gridg.Node(i,j,dv) =  gridg.Node(i,j,dv) +(1.0/tau)*(geq_Node[dv] - gridg.Node(i,j,dv)) 
                                            + (1 - 1.0/(2.0*tau))*lb9.thetaInverse *(
                                            gamma_Node[dv]*((lb9.Cx[dv] - ux)*Force.Node(i,j,0) + (lb9.Cy[dv] - uy)*Force.Node(i,j,1) ) 
                                            -   (gamma_Node[dv] - lb9.W[dv])*((lb9.Cx[dv] - ux)*grad_psi_rho.Node(i,j,0) + (lb9.Cy[dv] - uy)*grad_psi_rho.Node(i,j,1) ) 
                                    )   ;
            }
        }
    }


}
    




template<typename T>
void get_equi_f(double feq[9], lbmD2Q9<T> &lb9, double ux, double uy, double phi){


    double u2 = ux*ux + uy*uy;
    double a1=0;
    double first,second, third,feq0=0;
    for (int dv = 0; dv< 9; dv++){

        feq0 = phi*lb9.W[dv];

        first  = (ux*lb9.Cx[dv] + uy*lb9.Cy[dv])*lb9.thetaInverse;
        second = 0.5*(first * first);
        third = -0.5*u2*lb9.thetaInverse;
        feq[dv] = feq0*(1+ first + second + third);    


    }
}




template<typename T>
void get_equi_g(double geq[9], lbmD2Q9<T> &lb9, double ux, double uy, double rho, double p){


    double u2 = ux*ux + uy*uy;
    double a1=0;
    double first,second, third,geq0=0;

    for (int dv = 0; dv< 9; dv++){

        geq0 = lb9.W[dv];

        first  = (ux*lb9.Cx[dv] + uy*lb9.Cy[dv])*lb9.thetaInverse;
        second = 0.5*(first * first);
        third = -0.5*u2*lb9.thetaInverse;

        geq[dv] = geq0*(p + rho*lb9.theta0*(first + second + third));    
   
    }
}



template<typename T>
void get_equi_gamma(double gamma[9], lbmD2Q9<T> &lb9, double ux, double uy, double phi){


    double u2 = ux*ux + uy*uy;
    double a1=0;
    double first,second, third,gamma0=0;
    for (int dv = 0; dv< 9; dv++){

        gamma0 =  lb9.W[dv];

        first  = (ux*lb9.Cx[dv] + uy*lb9.Cy[dv])*lb9.thetaInverse;
        second = 0.5*(first * first);
        third = -0.5*u2*lb9.thetaInverse;
        gamma[dv] = gamma0*(1+ first + second + third);    
   
    }
}





template<typename T,typename T1>
void get_P_and_u(Grid_N_C_2D<T> &gridg,Grid_N_C_2D<T> &grad_psi_rho, lbmD2Q9<T1> &lb9,double &Ux, double &Uy,double &p,double rho , int i,int j, Grid_N_C_2D<T> &Force){ ///node or cell 0-Node 1- cell
    Ux  = 0.0;
    Uy  = 0.0;
    p   = 0;

    double sum_g = 0;
    for(int dv = 0; dv <9; dv++){

        Ux  += gridg.Node(i,j,dv)*lb9.Cx[dv];
        Uy  += gridg.Node(i,j,dv)*lb9.Cy[dv];
        p   += gridg.Node(i,j,dv);

    }  

    Ux = Ux/(rho*lb9.theta0) + 0.5*(Force.Node(i,j,0)); //dont include the density in F
    Uy = Uy/(rho*lb9.theta0) + 0.5*(Force.Node(i,j,1)); //dont include the density in F

    p = p - 0.5* (Ux* grad_psi_rho.Node(i,j,0) + Uy* grad_psi_rho.Node(i,j,1));

}


template<typename T,typename T1>
void get_phi(Grid_N_C_2D<T> &lbgrid, lbmD2Q9<T1> &lb9,double &phi,  int X, int Y){ ///node or cell 0-Node 1- cell
   
    phi = 0.0;
    for(int dv = 0; dv <lbgrid.d_v; dv++){

        phi += lbgrid.Node(X,Y,dv);

    }  
}



template<typename T, typename T1>
void initialization(Grid_N_C_2D<T> &gridf,Grid_N_C_2D<T> &gridg,lbmD2Q9<T1> &lb9, double phi_l ,double phi_h, double rho_l , double rho_h, double a, double b){

    double Feq_node[9] = {0},Geq_node[9] = {0},rho = 1.0;
    double x,y;

    double phi = 0;
    
    double ux = 0,  uy =  0;
    
    double u1 = 0.0, u2 = 0.0;

    double p = 0;
    
    for(int i = gridf.nbx; i <= gridf.nex ; i++){
        for(int j = gridf.nby;j <= gridf.ney ; j++){
            
            x = ((double)i)/ gridf.n_x ;
            y = ((double)j)/ gridf.n_x ;



            phi = (tanh((y - 2.0 - 0.1*cos(2.0*M_PI*x))/(sqrt(4.0) * (1.0/gridf.n_x) )));
            
            get_equi_f(Feq_node   ,lb9, ux, uy, phi);   //need phi, u 

            for (int dv = 0; dv< gridf.d_v; dv++){
                gridf.Node(i,j,dv) = Feq_node[dv];
            }

            rho = rho_l + ((phi - phi_l)/(phi_h - phi_l)) *(rho_h - rho_l);

            double eta = rho*b/4.0;

            p = rho *lb9.theta0*(1 + eta + eta*eta - eta*eta*eta)/pow(1 - eta, 3)  - a *rho*rho ;

            get_equi_g (Geq_node   ,lb9, ux, uy, rho,p);                               // need p and rho and u

            for(int dv = 0; dv< gridg.d_v; dv++){
                gridg.Node(i,j,dv) = Geq_node[dv];
            }

        }
    }
}







template<typename T, typename T1>
void initialization_ellipse(Grid_N_C_2D<T> &gridf,Grid_N_C_2D<T> &gridg,lbmD2Q9<T1> &lb9, double phi_l ,double phi_h, double rho_l , double rho_h, double a, double b){

    double Feq_node[9] = {0},Geq_node[9] = {0},rho = 1.0;
    double x,y;

    double phi = 0;
    
    double ux = 0,  uy =  0;
    
    double u1 = 0.0, u2 = 0.0;

    double p = 0;
    
    for(int i = gridf.nbx; i <= gridf.nex ; i++){
        for(int j = gridf.nby;j <= gridf.ney ; j++){
            
            x = ((double)i)/ gridf.n_x ;
            y = ((double)j)/ gridf.n_x ;



            phi = tanh( (0.2 - sqrt( (x -0.5)*(x - 0.5 ) + 0.5*(y - 0.5)*(y - 0.5) )  )/
                        (sqrt(2.0) * (1.5/ gridf.n_x) )  
                        );



            get_equi_f(Feq_node   ,lb9, ux, uy, phi);   //need phi, u 

            for (int dv = 0; dv< gridf.d_v; dv++){
                gridf.Node(i,j,dv) = Feq_node[dv];
            }

            rho = rho_l + ((phi - phi_l)/(phi_h - phi_l)) *(rho_h - rho_l);

            double eta = rho*b/4.0;

            p = rho *lb9.theta0*(1 + eta + eta*eta - eta*eta*eta)/pow(1 - eta, 3)  - a *rho*rho ;

            get_equi_g (Geq_node   ,lb9, ux, uy, rho,p);                               // need p and rho and u

            for(int dv = 0; dv< gridg.d_v; dv++){
                gridg.Node(i,j,dv) = Geq_node[dv];
            }

        }
    }
}




















template<typename T>
void printMass(Grid_N_C_2D<T> &gridf,Grid_N_C_2D<T> &gridg){    
    double a = 0;
    for(int i = 0 + gridf.noghost; i < gridf.n_x_node - (gridf.noghost) ; i++){
        for(int j = 0 + gridf.noghost;j < gridf.n_y_node - (gridf.noghost) ; j++){
            for (int dv = 0; dv< 9; dv++){
                a += gridf.Node(i,j,dv) ;
            }
        }
    }
    std::cout<<"   "<<a<<std::endl;

}
