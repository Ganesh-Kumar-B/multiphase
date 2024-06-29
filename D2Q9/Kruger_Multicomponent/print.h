#ifndef PRINT
#define PRINT

#include "lbmD2Q9.h"
#include "GRID_2D.h"
#include <fstream>
#include <iostream>

#include <sys/stat.h>
#include <sys/types.h>




template<typename T, typename T1>
void print_vtk(lbmD2Q9<T1> &lb,  Grid_N_C_2D<T> &gridf ,Grid_N_C_2D<T> &gridg,  int step, real u0,  Grid_N_C_2D<T> &Force)
{
   
    T u1,u2,um, rho1,rho2,del=0.05;

    std::ofstream file;
    char fileName[250];
    char foldername[250];
    sprintf(foldername,"Result_%.2f",0.1);
    mkdir(foldername,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    sprintf(fileName,"./Result_%.2f/velocity_%d.vtk", 0.1,step) ;
    file.open(fileName);

    file<<"# vtk DataFile Version 3.0\nVelocity\nASCII\nDATASET STRUCTURED_POINTS"<<std::endl;
    
    file<<"DIMENSIONS "<<1*gridf.n_x<<" "<<1*gridf.n_y<<" "<<1<<std::endl;
    
    file<<"ORIGIN "<<0<<" "<<0<<" "<<0<<std::endl;
    file<<"SPACING "<<1<<" "<<1<<" "<<1<<std::endl;

    file<<"POINT_DATA "<<1*gridf.n_x*1*gridf.n_y<<std::endl;
    file<<"SCALARS density double 1\nLOOKUP_TABLE default"<<std::endl;

    for (int j = 0 + gridf.noghost; j < gridf.n_y_node - (gridf.noghost); j++){
        for(int i = 0 + gridf.noghost; i < gridf.n_x_node - (gridf.noghost); i++){ 

            get_moments_Node_f(gridf,lb,u1, u2, rho1, i,j,Force);
            file<<rho1<<std::endl;

        }
    }
    

    

    file<<"SCALARS phi double 1\nLOOKUP_TABLE default"<<std::endl;


    for (int j = 0 + gridf.noghost; j < gridf.n_y_node - (gridf.noghost); j++){
        for(int i = 0 + gridf.noghost; i < gridf.n_x_node - (gridf.noghost); i++){ 
            
            real phi = 0;

            get_moments_Node_g(gridg,lb,phi,i,j);

            file<< phi <<std::endl;

        }
    }





    file<<"VECTORS velocity double"<<std::endl;

    for (int j = 0 + gridf.noghost; j < gridf.n_y_node - (gridf.noghost); j++){
        for(int i = 0 + gridf.noghost; i < gridf.n_x_node - (gridf.noghost); i++){ 

            get_moments_Node_f(gridf,lb,u1, u2, rho1, i,j,Force);

            file<<u1<<" "<<u2<<" "<<0.0<<std::endl;
            
        }
    }



}



#endif