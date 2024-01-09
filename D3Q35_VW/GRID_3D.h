#ifndef Grid_3D_D2Q19
#define Grid_3D_D2Q19



#include<iostream>


template<typename T>
class Grid_N_C_3D{
    public:
int nNodes ,dv,total_mem_node, noghost;

int n_x_node,n_y_node,n_z_node,
    n_x     ,n_y     ,n_z     ,
    d_v ;

T *data_node, *data_cell;

Grid_N_C_3D(int nx, int ny, int nz,int no_ghost, int dv){
    
    n_x = nx;
    n_y = ny;
    n_z = nz;

    noghost = no_ghost;
  


    n_x_node = nx + 2*no_ghost;    //2 ghost nodes in each direction
    n_y_node = ny + 2*no_ghost;
    n_z_node = nz + 2*no_ghost;

    d_v = dv;                      //No of population
    



    nNodes = n_x_node*n_y_node*n_z_node; 

    total_mem_node = nNodes*d_v;
   
    data_node = new T[total_mem_node]{0};
    data_cell = new T[total_mem_node]{0};
    
}


~Grid_N_C_3D(){
    delete [] data_node;
    delete [] data_cell;
}
T Node(const int i,const int j,const int k,const int dv = 0) const  {return data_node[((i* n_y_node + j)*n_z_node + k )*d_v + dv ];} //start 0 from ghost node
T& Node(const int i, const int j,const int k,const  int dv = 0)     {return data_node[((i* n_y_node + j)*n_z_node + k )*d_v + dv ];}   // start 0 from ghost node

T Cell(const int i,const int j,const int k,const int dv = 0) const  {return data_cell[((i* n_y_node + j)*n_z_node + k )*d_v + dv ];} //start 0 from ghost node
T& Cell(const int i, const int j,const int k,const  int dv = 0)     {return data_cell[((i* n_y_node + j)*n_z_node + k )*d_v + dv ];}   // start 0 from ghost node

};





#endif