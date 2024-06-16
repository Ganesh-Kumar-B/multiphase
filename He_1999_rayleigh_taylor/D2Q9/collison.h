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
void collide(Grid_N_C_2D<T> &gridf,Grid_N_C_2D<T> &gridg,lbmD2Q9<T1> &lbD2Q9,double beta,double tau, double gx, double gy ){
    double feq_Node[9] = {0},feq_Cell[9]={0}, ux = 0, uy = 0, rho = 0, G[9]={0};
    int a = 0,b = 1; //just for representation of node or cell 

    Grid_N_C_2D<T> phi(gridf.n_x_node,gridf.n_y_node,1,9);
    Grid_N_C_2D<T> rho(gridf.n_x_node,gridf.n_y_node,1,9);



    Multiphase_terms(gridf,gridg, phi, rho);


    //  // first the population of nodes are resetted and second  the population of the cells are resetted
    for(int i = 0 + gridf.noghost; i < gridf.n_x_node - (gridf.noghost) ; i++){
        for(int j = 0 + gridf.noghost;j < gridf.n_y_node - (gridf.noghost) ; j++){
            
            
            get_moments(gridf, lbD2Q9,  ux, uy,rho, i, j);            //for the node
            get_equi(feq_Node,lbD2Q9, ux, uy, rho);


            for (int dv = 0; dv< 9; dv++){      
                gridf.Node(i,j,dv) =  gridf.Node(i,j,dv) + 2.0* beta*(feq_Node[dv] - gridf.Node(i,j,dv)) 
                                    - (1 - 1.0/(2.0*tau))    ;
            }




        }
    }


}
    





template<typename T>
void get_equi(double feq[9], lbmD2Q9<T> &lbD2Q9, double ux, double uy, double rho){


    double u2 = ux*ux + uy*uy;
     double a1=0;
    double first,second, third,feq0=0;
    for (int dv = 0; dv< 9; dv++){

        feq0 = rho*lbD2Q9.W[dv];

        first  = (ux*lbD2Q9.Cx[dv] + uy*lbD2Q9.Cy[dv])*lbD2Q9.thetaInverse;
        second = 0.5*(first * first);
        third = -0.5*u2*lbD2Q9.thetaInverse;
        feq[dv] = feq0*(1+ first + second + third);    
   
    }
}





template<typename T,typename T1>
void get_moments(Grid_N_C_2D<T> &lbgrid, lbmD2Q9<T1> &lbD2Q9,double &Ux, double &Uy,double &Rho,  int X, int Y){ ///node or cell 0-Node 1- cell
   Ux  = 0.0;
   Uy  = 0.0;
   Rho = 0.0;


    for(int dv = 0; dv <9; dv++){

        Ux  += lbgrid.Node(X,Y,dv)*lbD2Q9.Cx[dv];
        Uy  += lbgrid.Node(X,Y,dv)*lbD2Q9.Cy[dv];
        Rho += lbgrid.Node(X,Y,dv);

    }  

    Ux = Ux/Rho;
    Uy = Uy/Rho;  

}


template<typename T, typename T1>
void initialization(Grid_N_C_2D<T> &gridf,Grid_N_C_2D<T> &gridg,lbmD2Q9<T1> &lbD2Q9, double U0,double u1){

    double  Feq_node[9] = {0},Rho = 1.0;
    double  x,y,
            n_to_n_dist = (2*PI)/gridf.n_x;    ///distance between nodes 
    
    double ux_node = 0,uy_node =0;

    get_equi(Feq_node,lbD2Q9,ux_node,uy_node,Rho);
    
    for(int i = gridf.nbx; i <= gridf.nex ; i++){
        for(int j = gridf.nby;j <= gridf.ney ; j++){
            
            for(int dv = 0; dv< 9; dv++){
                gridf.Node(i,j,dv) = Feq_node[dv];
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
