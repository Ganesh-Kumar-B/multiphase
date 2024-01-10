#pragma once

#include "lbmD3Q35.h"
#include "GRID_3D.h"
#include <fstream>
#include <iostream>








template<typename T, typename T1>
void printdata(lbmD3Q35<T1> &lbModel,  Grid_N_C_3D<T> &gridLB,  int step, double u0)
{

    double u_inv,nx_inv, ny_inv;
    u_inv = 1/ u0;
    
    nx_inv = 1/(double)gridLB.n_x;
    ny_inv = 1/(double)gridLB.n_y;


    std::vector<int> line;bool isPresent;
    T u1,u2,u3,u4,um, rho1,rho2,del=0.05;
    std::ofstream file;
    char fileName[250];
    sprintf(fileName,"./Result/velocity_%d.txt",step) ;
    file.open(fileName) ;
    file<<"x,y,z,ux,uy,uz,rho"<<std::endl;

    for(int i = 0 + gridLB.noghost; i < gridLB.n_x_node - (gridLB.noghost); i++){ 
	    for (int j = 0 + gridLB.noghost; j < gridLB.n_y_node - (gridLB.noghost); j++){
            for(int k = 0 + gridLB.noghost; k < gridLB.n_z_node - (gridLB.noghost); k++){ 


                get_moments(gridLB,lbModel,u1, u2,u3, rho1, i,j,k);
                


                // file<<(double)(i-gridLB.noghost)*ny_inv  <<"," <<(double)(j-gridLB.noghost) *ny_inv   <<","<<u1*u_inv<<","<<u2*u_inv<<","<<rho1<<std::endl;
                file<<i<<"," <<j  <<","<<k  <<","<<u1<<","<<u2<<","<<u3<<","<<rho1<<std::endl;
    

            }
        }
    }
}


