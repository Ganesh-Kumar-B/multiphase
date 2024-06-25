#ifndef Grid_2D_D2Q9
#define Grid_2D_D2Q9



#include<iostream>
typedef double real;


template<typename T>
class Grid_N_C_2D{
    public:
int nNodes ,dv,total_mem_node, noghost;

int n_x_node,n_y_node,
    n_x, n_y,
    nbx, nby,
    nex, ney,
    d_v;

T *data_node;

Grid_N_C_2D(int nx, int ny,int no_ghost, int dv){
    
    n_x = nx;
    n_y = ny;

    noghost = no_ghost;
  


    n_x_node = nx + 2*no_ghost;    //2 ghost nodes in each direction
    n_y_node = ny + 2*no_ghost;
   
    d_v = dv;                      //No of population
    


    nbx = no_ghost;                     //physical domain starting
    nby = no_ghost;

    nex = n_x_node - no_ghost -1;       //physical domain ending
    ney = n_y_node - no_ghost -1;


    nNodes = (n_x_node)*(n_y_node); 
   

    total_mem_node = nNodes*dv;
   
    data_node = new T[total_mem_node]{0};
}


~Grid_N_C_2D(){
    delete [] data_node;
}
T Node(const int i,const int j,const int dv = 0) const{return data_node[(i*n_y_node +j)*d_v + dv ];} //start 0 from ghost node
T& Node(const int i, const int j,const  int dv = 0) {return data_node[(i*n_y_node +j)*d_v + dv ];}   // start 0 from ghost node



};





#endif