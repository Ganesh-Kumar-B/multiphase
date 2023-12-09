#ifndef PRINT
#define PRINT

#include "lbmD2Q9.h"
#include "GRID_2D.h"
#include <fstream>
#include <iostream>








template<typename T, typename T1>
void printdata(lbmD2Q9<T1> &lbModel,  Grid_N_C_2D<T> &gridLB,  int step, double u0)
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
  file<<"x,y,ux,uy,rho"<<std::endl;
  for(int i = 0 + gridLB.noghost; i < gridLB.n_x_node - (gridLB.noghost); i++) 
	for (int j = 0 + gridLB.noghost; j < gridLB.n_y_node - (gridLB.noghost); j++)
		{{

    get_moments(gridLB,lbModel,u1, u2, rho1, i,j);
    // if(fabs(u1) < 0.0000001){u1 = 0; }
    // if(fabs(u2) < 0.0000001){u2 = 0; }

    //  get_moments(gridLB,lbModel,u3,u4, rho2, i,j,1);
    // if(fabs(u3) < 0.00000000000001){u3 = 0;}
    // if(fabs(u4) < 0.000000000001){u4 = 0;}
  file<<(double)(i-gridLB.noghost)*ny_inv  <<"," <<(double)(j-gridLB.noghost) *ny_inv   <<","<<u1*u_inv<<","<<u2*u_inv<<","<<rho1<<std::endl;
    

}
}
}


#endif