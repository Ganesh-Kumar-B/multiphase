#ifndef ADVECTION
#define ADVECTION
//till second order
// // enum velocityDir{ 
// // dV_ZERO_ZERO,
// // dV_P1_ZERO, dV_ZERO_P1, dV_M1_ZERO, dV_ZERO_M1,
// // dV_P3_ZERO, dV_ZERO_P3, dV_M3_ZERO, dV_ZERO_M3,
// // dV_PH1_PH1, dV_MH1_PH1, dV_MH1_MH1, dV_PH1_MH1
// // }; 
#include"lbmD3Q19.h"
#include"GRID_3D.h"
#include<iomanip>




template<typename T>
void advection(Grid_N_C_3D<T> &grid){


    for(int i = 0 + grid.noghost ; i < grid.n_x_node - (grid.noghost);i++){
        for(int j = 0 + grid.noghost ; j < grid.n_y_node - (grid.noghost);j++){ 
            for(int k = 0 + grid.noghost ; k < grid.n_z_node - (grid.noghost);k++){ 

            
                grid.Node(i,j,k,  dV_M1_ZERO_ZERO   )     = grid.Node(i+1,j  ,k  , dV_M1_ZERO_ZERO ); 
                grid.Node(i,j,k,  dV_ZERO_M1_ZERO   )     = grid.Node(i  ,j+1,k  , dV_ZERO_M1_ZERO ); 
                grid.Node(i,j,k,  dV_ZERO_ZERO_M1   )     = grid.Node(i  ,j  ,k+1, dV_ZERO_ZERO_M1 ); 
                grid.Node(i,j,k,  dV_M1_M1_ZERO  )        = grid.Node(i+1,j+1,k  , dV_M1_M1_ZERO   ); 
                grid.Node(i,j,k,  dV_M1_P1_ZERO  )        = grid.Node(i+1,j-1,k  , dV_M1_P1_ZERO   ); 
                grid.Node(i,j,k,  dV_ZERO_M1_M1  )        = grid.Node(i  ,j+1,k+1, dV_ZERO_M1_M1   ); 
                grid.Node(i,j,k,  dV_ZERO_M1_P1  )        = grid.Node(i  ,j+1,k-1, dV_ZERO_M1_P1   ); 
                grid.Node(i,j,k,  dV_M1_ZERO_M1  )        = grid.Node(i+1,j  ,k+1, dV_M1_ZERO_M1   ); 
                grid.Node(i,j,k,  dV_M1_ZERO_P1  )        = grid.Node(i+1,j  ,k-1, dV_M1_ZERO_P1   ); 

            
                
            }
        }
    }



    for(int i = grid.n_x_node - (grid.noghost +1); i> grid.noghost-1; i--){                                                                                        //  #pragma omp parallel for shared(b)
        for(int j = grid.n_y_node - (grid.noghost +1); j> grid.noghost-1; j--){
            for(int k = grid.n_z_node - (grid.noghost +1); k> grid.noghost-1; k--){

                grid.Node(i,j,k,dV_P1_ZERO_ZERO   ) = grid.Node(i-1,j  ,k  , dV_P1_ZERO_ZERO );      
                grid.Node(i,j,k,dV_ZERO_P1_ZERO   ) = grid.Node(i  ,j-1,k  , dV_ZERO_P1_ZERO );      
                grid.Node(i,j,k,dV_ZERO_ZERO_P1   ) = grid.Node(i  ,j  ,k-1, dV_ZERO_ZERO_P1 );      
                grid.Node(i,j,k,dV_P1_M1_ZERO   )   = grid.Node(i-1,j+1,k  , dV_P1_M1_ZERO   );    
                grid.Node(i,j,k,dV_ZERO_P1_P1   )   = grid.Node(i  ,j-1,k-1, dV_ZERO_P1_P1   );    
                grid.Node(i,j,k,dV_ZERO_P1_M1   )   = grid.Node(i  ,j-1,k+1, dV_ZERO_P1_M1   );    
                grid.Node(i,j,k,dV_P1_P1_ZERO   )   = grid.Node(i-1,j-1,k  , dV_P1_P1_ZERO   );    
                grid.Node(i,j,k,dV_P1_ZERO_P1   )   = grid.Node(i-1,j  ,k-1, dV_P1_ZERO_P1   );    
                grid.Node(i,j,k,dV_P1_ZERO_M1   )   = grid.Node(i-1,j  ,k+1, dV_P1_ZERO_M1   );    
            }       
        }
    }

}





template<typename T>
void Periodic(Grid_N_C_3D<T> &grid){
    
    for(int i = 0 ; i < grid.n_x_node; i++)
        for(int j = 0 ; j < grid.n_y_node;j++ ){
            for(int dv = 0; dv<grid.d_v; dv++){

                grid.Node(i ,j  ,0                            , dv)   = grid.Node(i,j,grid.n_z_node-grid.noghost -1 ,dv);
                grid.Node(i ,j  ,grid.n_z_node - grid.noghost , dv)   = grid.Node(i,j,1                             ,dv);

            }
        }

    for(int i = 0 ; i < grid.n_x_node;i++ ){
        for(int k = 0 ; k < grid.n_z_node;k++ ){
            for(int dv = 0; dv<grid.d_v; dv++){

                grid.Node(i ,0                            ,k, dv)   = grid.Node(i,grid.n_y_node-grid.noghost -1 ,k,dv);
                grid.Node(i ,grid.n_y_node - grid.noghost ,k, dv)   = grid.Node(i,1                             ,k,dv);

            }
        }
    }

    for(int j = 0 ; j < grid.n_y_node;j++ ){
        for(int k = 0 ; k < grid.n_z_node;k++ ){
            for(int dv = 0; dv<grid.d_v; dv++){

                grid.Node(0                            ,j,k, dv)   = grid.Node(grid.n_x_node-grid.noghost -1 ,j,k,dv);
                grid.Node(grid.n_x_node - grid.noghost ,j,k, dv)   = grid.Node(1                             ,j,k,dv);

            }
        }
    }

}


















;



#endif