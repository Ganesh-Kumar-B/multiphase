#ifndef PRINT
#define PRINT

#include "lbmD2Q9.h"
#include "GRID_2D.h"
#include <fstream>
#include <iostream>

#include <sys/stat.h>
#include <sys/types.h>





template<typename T, typename T1>
void print_vtk(lbmD2Q9<T1> &lb,  Grid_N_C_2D<T> &grid,  int step, real u0, real theta,real kappa, Grid_N_C_2D<T> &Force, const std::string &name)
{
    T u1,u2,u3,u4,um, rho1,rho2,del=0.05;

    std::ofstream file;
    char fileName[250];
    char foldername[250];
    sprintf(foldername,"%s_%.2f",name.c_str(),theta);
    mkdir(foldername,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    sprintf(fileName,"./%s_%.2f/velocity_%d.vtk",name.c_str(), theta,step) ;
    file.open(fileName);
    // file<<"# vtk DataFile Version 3.0\nVelocity\nASCII\nDATASET STRUCTURED_GRID"<<std::endl;
    file<<"# vtk DataFile Version 3.0\nVelocity\nASCII\nDATASET STRUCTURED_POINTS"<<std::endl;
    
    file<<"DIMENSIONS "<<1*grid.n_x<<" "<<1*grid.n_y<<" "<<1<<std::endl;
    
    file<<"ORIGIN "<<0<<" "<<0<<" "<<0<<std::endl;
    file<<"SPACING "<<1<<" "<<1<<" "<<1<<std::endl;


    // file<<"POINTS "<<1*grid.n_x*1*grid.n_y*1*grid.n_z<<" double"<<std::endl;

    
    // for(int k = 0 + grid.noghost; k < grid.n_z_node - (grid.noghost); k++){
	//     for (int j = 0 + grid.noghost; j < grid.n_y_node - (grid.noghost); j++){
    //         for(int i = 0 + grid.noghost; i < grid.n_x_node - (grid.noghost); i++){ 

    //             file<<i<<" "<<j     <<" "<<k     <<std::endl;
    //             file<<i +0.5<<" "<<j +0.5<<" "<<k +0.5<<std::endl;


    //         }
    //     }
    // }
        

    file<<"POINT_DATA "<<1*grid.n_x*1*grid.n_y<<std::endl;
    file<<"SCALARS density double 1\nLOOKUP_TABLE default"<<std::endl;

    for (int j = 0 + grid.noghost; j < grid.n_y_node - (grid.noghost); j++){
        for(int i = 0 + grid.noghost; i < grid.n_x_node - (grid.noghost); i++){ 

            get_moments_Node(grid,lb,u1, u2, rho1, i,j,Force);
            file<<rho1<<std::endl;

        }
    }
    

    file<<"VECTORS velocity double"<<std::endl;

    for (int j = 0 + grid.noghost; j < grid.n_y_node - (grid.noghost); j++){
        for(int i = 0 + grid.noghost; i < grid.n_x_node - (grid.noghost); i++){ 

            get_moments_Node(grid,lb,u1, u2, rho1, i,j,Force);

            file<<u1<<" "<<u2<<" "<<0.0<<std::endl;

        }
    }
    



}





#endif