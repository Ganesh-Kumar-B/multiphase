#pragma once

#include"lbmD3Q35.h"
#include"GRID_3D.h"
#include"collison.h"
#include<iomanip>



template<typename T, typename T1>
void Diffuse_35(Grid_N_C_3D<T> &grid,  lbmD3Q35<T1> &lb, real utop, real ubottom){

        real feq[35] ={0};
        real rho_wall_top = 1.0, rho_wall_bottom = 1.0;

        // #  ------------------------------------------------------------------------------=--top wall

        int topC_last[]         =   { 
                                    dV_ZERO_P2_ZERO ,
                                    dV_P1_P1_ZERO   ,dV_M1_P1_ZERO  ,dV_ZERO_P1_P1  ,dV_ZERO_P1_M1  ,
                                    dV_P2_P2_P2     ,dV_M2_P2_P2    ,dV_P2_P2_M2    ,dV_M2_P2_M2    ,
                                    dV_PH1_PH1_PH1  ,dV_MH1_PH1_PH1 ,dV_PH1_PH1_MH1 , dV_MH1_PH1_MH1
                                    };

        int topC_2nd_last[]     =   {
                                    dV_ZERO_P2_ZERO ,
                                    dV_P2_P2_P2     ,dV_M2_P2_P2    ,dV_P2_P2_M2    ,dV_M2_P2_M2,
                                    };
                                    
        int topN_last[]         =   {
                                    dV_ZERO_P2_ZERO ,
                                    dV_P1_P1_ZERO   ,dV_M1_P1_ZERO  ,dV_ZERO_P1_P1  ,dV_ZERO_P1_M1  ,
                                    dV_P2_P2_P2     ,dV_M2_P2_P2    ,dV_P2_P2_M2    ,dV_M2_P2_M2    ,
                                    };

        int topN_2nd_last[]   =     {
                                    dV_ZERO_P2_ZERO,
                                    dV_P2_P2_P2     ,dV_M2_P2_P2    ,dV_P2_P2_M2    ,dV_M2_P2_M2,  
                                    }
                                    ;

    get_equi(feq,lb,utop,0.0,0.0,rho_wall_top);


    for(int i = 0 + grid.nbx ; i <= grid.nex;i++){
        for(int k = 0 + grid.nbz;k <= grid.nez; k++){ 

            //top wall
            real G_num_cell = 0, G_den_cell = 0, G_num_node = 0, G_den_node = 0,G_node = 0,G_cell = 0;

            for(int l = 0; l<sizeof(topC_last)/sizeof(topC_last[0]); l++){
                G_num_cell += grid.Cell(i   ,grid.ney   ,k  , topC_last[l]  ) * abs(lb.Cy[topC_last[l]]);
                G_den_cell += feq[topC_last[l]] *abs(lb.Cy[topC_last[l]]);       
            }
                
            for(int l = 0; l<sizeof(topC_2nd_last)/sizeof(topC_2nd_last[0]); l++){
                G_num_cell += grid.Cell(i   ,grid.ney -1    ,k  , topC_2nd_last[l]) * abs(lb.Cy[topC_2nd_last[l]]);
                G_den_cell += feq[topC_2nd_last[l]] *abs(lb.Cy[topC_2nd_last[l]]) ;       
            }

            G_cell = G_num_cell/G_den_cell;

            for(int l = 0; l<sizeof(topN_last)/sizeof(topN_last[0]); l++){
                G_num_node += grid.Node(i   ,grid.ney   ,k  , topN_last[l]) * abs(lb.Cy[topN_last[l]]);
                G_den_node += feq[topN_last[l]] *abs(lb.Cy[topN_last[l]]) ;       
            }

            for(int l = 0; l<sizeof(topN_2nd_last)/sizeof(topN_2nd_last[0]); l++){
                G_num_node += grid.Node(i   ,grid.ney -1    ,k  , topN_2nd_last[l]) * abs(lb.Cy[topN_2nd_last[l]]);    
                G_den_node += feq[topN_2nd_last[l]] *abs(lb.Cy[topN_2nd_last[l]]) ;       
            }

            G_node = G_num_node / G_den_node;

            // pushing to the last node and cell
            int j = grid.ney;

            for(int l = 0; l < 9; l++) 
                grid.Cell(i + (int)lb.Cx[topC_last[l]], j +  (int)lb.Cy[topC_last[l]],  k +  (int)lb.Cz[topC_last[l]] , oppdV[topC_last[l]]) = G_cell * feq[oppdV[topC_last[l]]];

            for(int l = 9; l < 13;l++)
                grid.Node(i + (int)lb.CxC[topC_last[l]], j + (int)lb.CyC[topC_last[l]],  k +  (int)lb.CzC[topC_last[l]] ,oppdV[topC_last[l]]) = G_cell * feq[oppdV[topC_last[l]]];

            for(int l = 0; l < 9;l++)
                grid.Node(i + (int)lb.Cx[topN_last[l]], j + (int)lb.Cy[topN_last[l]] ,  k +  (int)lb.Cz[topN_last[l]] ,oppdV[topN_last[l]]) = G_node * feq[oppdV[topN_last[l]]];


            //pushing to the last cell
            j = grid.ney - 1;

            for(int l = 0; l<5; l++)
                grid.Cell(i + (int)lb.Cx[topC_2nd_last[l]], j +  (int)lb.Cy[topC_2nd_last[l]],  k +  (int)lb.Cz[topC_2nd_last[l]] ,oppdV[topC_2nd_last[l]]) = G_cell * feq[oppdV[topC_2nd_last[l]]];

            for(int l = 0; l<5; l++)
                grid.Node(i + (int)lb.Cx[topN_2nd_last[l]], j +  (int)lb.Cy[topN_2nd_last[l]],  k +  (int)lb.Cz[topN_2nd_last[l]] ,oppdV[topN_2nd_last[l]]) = G_node * feq[oppdV[topN_2nd_last[l]]];



        }
    }    


    //     //#--bottom wall

    //     int botomC_last[]         =   { 

    //                                 dV_ZERO_M2_ZERO ,
    //                                 dV_P1_M1_ZERO   ,dV_M1_M1_ZERO  ,dV_ZERO_M1_P1  ,dV_ZERO_M1_M1  ,
    //                                 dV_P2_M2_P2     ,dV_M2_M2_P2    ,dV_P2_M2_M2    ,dV_M2_M2_M2    ,

    //                                 };

    //     int botomC_2nd_last[]     =   {
    //                                 dV_ZERO_M2_ZERO ,
    //                                 dV_P2_M2_P2     ,dV_M2_M2_P2    ,dV_P2_M2_M2    ,dV_M2_M2_M2,
    //                                 };
                                    


    //     int botomN_last[]         =   {
    //                                 dV_ZERO_M2_ZERO ,
    //                                 dV_P1_M1_ZERO   ,dV_M1_M1_ZERO  ,dV_ZERO_M1_P1  ,dV_ZERO_M1_M1  ,
    //                                 dV_P2_M2_P2     ,dV_M2_M2_P2    ,dV_P2_M2_M2    ,dV_M2_M2_M2    ,
    //                                 dV_PH1_MH1_PH1  ,dV_MH1_MH1_PH1 ,dV_PH1_MH1_MH1 , dV_MH1_MH1_MH1
    //                                 };

    //     int botomN_2nd_last[]   =   {
    //                                 dV_ZERO_M2_ZERO ,
    //                                 dV_P2_M2_P2     ,dV_M2_M2_P2    ,dV_P2_M2_M2    ,dV_M2_M2_M2,  
    //                                 }
    //                                 ;


    // get_equi(feq,lb,ubottom,0.0,0.0,rho_wall_bottom);


    // for(int i = 0 + grid.nbx ; i <= grid.nex;i++){
    //     for(int k = 0 + grid.nbz;k <= grid.nez; k++){ 

    //         //botom wall
    //         real G_num_cell = 0, G_den_cell = 0, G_num_node = 0, G_den_node = 0,G_node = 0,G_cell = 0;

    //         for(int l = 0; l<sizeof(botomC_last)/sizeof(botomC_last[0]); l++){
    //             G_num_cell += grid.Cell(i   ,grid.nby   ,k  , botomC_last[l]  ) * abs(lb.Cy[botomC_last[l]]);
    //             G_den_cell += feq[botomC_last[l]] *abs(lb.Cy[botomC_last[l]]);       
    //         }
                
    //         for(int l = 0; l<sizeof(botomC_2nd_last)/sizeof(botomC_2nd_last[0]); l++){
    //             G_num_cell += grid.Cell(i   ,grid.nby +1     ,k  , botomC_2nd_last[l]) * abs(lb.Cy[botomC_2nd_last[l]]);
    //             G_den_cell += feq[botomC_2nd_last[l]] *abs(lb.Cy[botomC_2nd_last[l]]) ;       
    //         }

    //         G_cell = G_num_cell / G_den_cell;

    //         for(int l = 0; l<sizeof(botomN_last)/sizeof(botomN_last[0]); l++){
    //             G_num_node += grid.Node(i   ,grid.nby   ,k  , botomN_last[l]) * abs(lb.Cy[botomN_last[l]]);
    //             G_den_node += feq[botomN_last[l]] *abs(lb.Cy[botomN_last[l]]) ;       
    //         }

    //         for(int l = 0; l<sizeof(botomN_2nd_last)/sizeof(botomN_2nd_last[0]); l++){
    //             G_num_node += grid.Node(i   ,grid.nby +1     ,k  , botomN_2nd_last[l]) * abs(lb.Cy[botomN_2nd_last[l]]);    
    //             G_den_node += feq[botomN_2nd_last[l]] *abs(lb.Cy[botomN_2nd_last[l]]) ;       
    //         }

    //         G_node = G_num_node / G_den_node;

    //         // pushing to the last node and cell
    //         int j = grid.nby;

    //         for(int l = 0; l < 9; l++) 
    //             grid.Cell(i + (int)lb.Cx[botomC_last[l]], j +  (int)lb.Cy[botomC_last[l]],  k +  (int)lb.Cz[botomC_last[l]] , oppdV[botomC_last[l]]) = G_cell * feq[oppdV[botomC_last[l]]];

    //         for(int l = 0; l < 9;l++)
    //             grid.Node(i + (int)lb.Cx[botomN_last[l]], j + (int)lb.Cy[botomN_last[l]] ,  k +  (int)lb.Cz[botomN_last[l]] ,oppdV[botomN_last[l]]) = G_node * feq[oppdV[botomN_last[l]]];

    //          for(int l = 9; l < 13;l++)
    //             grid.Cell(i + (int)lb.CxF[botomN_last[l]], j + (int)lb.CyF[botomN_last[l]] ,  k +  (int)lb.CzF[botomN_last[l]] ,oppdV[botomN_last[l]]) = G_node * feq[oppdV[botomN_last[l]]];


    //         //pushing to the last cell
    //         j = grid.nby + 1;

    //         for(int l = 0; l<5; l++)
    //             grid.Cell(i + (int)lb.Cx[botomC_2nd_last[l]], j +  (int)lb.Cy[botomC_2nd_last[l]],  k +  (int)lb.Cz[botomC_2nd_last[l]] ,oppdV[botomC_2nd_last[l]]) = G_cell * feq[oppdV[botomC_2nd_last[l]]];

    //         for(int l = 0; l<5; l++)
    //             grid.Node(i + (int)lb.Cx[botomN_2nd_last[l]], j +  (int)lb.Cy[botomN_2nd_last[l]],  k +  (int)lb.Cz[botomN_2nd_last[l]] ,oppdV[botomN_2nd_last[l]]) = G_node * feq[oppdV[botomN_2nd_last[l]]];



    //     }
    // }    
}   



template<typename T, typename T1>
void BB_wall_bottom(Grid_N_C_3D<T> &grid,  lbmD3Q35<T1> &lb, real utop, real ubottom){

    //#--bottom wall

    int botomC_last[]       =   { 
                                dV_ZERO_M2_ZERO ,

                                dV_P1_M1_ZERO   ,dV_M1_M1_ZERO  ,dV_ZERO_M1_P1  ,dV_ZERO_M1_M1  ,
                                dV_P2_M2_P2     ,dV_M2_M2_P2    ,dV_P2_M2_M2    ,dV_M2_M2_M2    ,
                                };



    int botomC_2nd_last[]   =   {
                                dV_ZERO_M2_ZERO ,
                                dV_P2_M2_P2     ,dV_M2_M2_P2    ,dV_P2_M2_M2    ,dV_M2_M2_M2,
                                };
                                
    

    int botomN_last[]       =   {
                                dV_ZERO_M2_ZERO ,
                                dV_P1_M1_ZERO   ,dV_M1_M1_ZERO  ,dV_ZERO_M1_P1  ,dV_ZERO_M1_M1  ,
                                dV_P2_M2_P2     ,dV_M2_M2_P2    ,dV_P2_M2_M2    ,dV_M2_M2_M2    ,
                                dV_PH1_MH1_PH1  ,dV_MH1_MH1_PH1 ,dV_PH1_MH1_MH1 , dV_MH1_MH1_MH1
                                };

    int botomN_2nd_last[]   =   {
                                dV_ZERO_M2_ZERO ,
                                dV_P2_M2_P2     ,dV_M2_M2_P2    ,dV_P2_M2_M2    ,dV_M2_M2_M2,  
                                }
                                ;




    for(int i = 0 + grid.nbx ; i <= grid.nex;i++){
        for(int k = 0 + grid.nbz;k <= grid.nez; k++){ 
            
            // last
            int j = grid.nby;

            for(int l = 0; l < 9; l++) 
                grid.Cell(i + (int)lb.Cx[botomC_last[l]], j +  (int)lb.Cy[botomC_last[l]],  k +  (int)lb.Cz[botomC_last[l]] , oppdV[botomC_last[l]]) = grid.Cell(i,j,k, botomC_last[l]);

            for(int l = 0; l < 9;l++)
                grid.Node(i + (int)lb.Cx[botomN_last[l]], j + (int)lb.Cy[botomN_last[l]] ,  k +  (int)lb.Cz[botomN_last[l]] ,oppdV[botomN_last[l]]) = grid.Node(i,j,k, botomN_last[l]);

             for(int l = 9; l < 13;l++)
                grid.Cell(i + (int)lb.CxF[botomN_last[l]], j + (int)lb.CyF[botomN_last[l]] ,  k +  (int)lb.CzF[botomN_last[l]] ,oppdV[botomN_last[l]]) = grid.Node(i,j,k, botomN_last[l]);

            //second last
            j = grid.nby + 1;

            for(int l = 0; l<5; l++)
                grid.Cell(i + (int)lb.Cx[botomC_2nd_last[l]], j +  (int)lb.Cy[botomC_2nd_last[l]],  k +  (int)lb.Cz[botomC_2nd_last[l]] ,oppdV[botomC_2nd_last[l]]) = grid.Cell(i,j,k,botomC_2nd_last[l]);

            for(int l = 0; l<5; l++)
                grid.Node(i + (int)lb.Cx[botomN_2nd_last[l]], j +  (int)lb.Cy[botomN_2nd_last[l]],  k +  (int)lb.Cz[botomN_2nd_last[l]] ,oppdV[botomN_2nd_last[l]]) = grid.Node(i,j,k,botomN_2nd_last[l]);

        }
    }

}




template<typename T, typename T1>
void BB_wall_left(Grid_N_C_3D<T> &grid,  lbmD3Q35<T1> &lb, real utop, real ubottom){

    //#--bottom wall

    int leftC_last[]       =   { 
                                dV_M2_ZERO_ZERO ,
                                
                                dV_M1_P1_ZERO   ,dV_M1_M1_ZERO  ,dV_M1_ZERO_P1  ,dV_M1_ZERO_M1  ,

                                dV_M2_P2_P2     ,dV_M2_P2_M2    ,dV_M2_M2_P2    ,dV_M2_M2_M2    ,
                                };



    int leftC_2nd_last[]   =   {
                                dV_M2_ZERO_ZERO ,

                                dV_M2_P2_P2     ,dV_M2_P2_M2    ,dV_M2_M2_P2    ,dV_M2_M2_M2    ,
                                };
                                


    int leftN_last[]       =   {
                                dV_M2_ZERO_ZERO ,
                                
                                dV_M1_P1_ZERO   ,dV_M1_M1_ZERO  ,dV_M1_ZERO_P1  ,dV_M1_ZERO_M1  ,

                                dV_M2_P2_P2     ,dV_M2_P2_M2    ,dV_M2_M2_P2    ,dV_M2_M2_M2    ,


                                dV_MH1_PH1_PH1     ,dV_MH1_PH1_MH1    ,dV_MH1_MH1_PH1    ,dV_MH1_MH1_MH1    
                                };

    int leftN_2nd_last[]   =   {
                                dV_M2_ZERO_ZERO ,

                                dV_M2_P2_P2     ,dV_M2_P2_M2    ,dV_M2_M2_P2    ,dV_M2_M2_M2    ,  
                                }
                                ;




    for(int j = 0 + grid.nby ; j <= grid.ney;j++){
        for(int k = 0 + grid.nbz;k <= grid.nez; k++){ 
            
            // last
            int i = grid.nbx;

            for(int l = 0; l < 9; l++) 
                grid.Cell(i + (int)lb.Cx[leftC_last[l]], j +  (int)lb.Cy[leftC_last[l]],  k +  (int)lb.Cz[leftC_last[l]] , oppdV[leftC_last[l]]) = grid.Cell(i,j,k, leftC_last[l]);

            for(int l = 0; l < 9;l++)
                grid.Node(i + (int)lb.Cx[leftN_last[l]], j + (int)lb.Cy[leftN_last[l]] ,  k +  (int)lb.Cz[leftN_last[l]] ,oppdV[leftN_last[l]]) = grid.Node(i,j,k, leftN_last[l]);

             for(int l = 9; l < 13;l++)
                grid.Cell(i + (int)lb.CxF[leftN_last[l]], j + (int)lb.CyF[leftN_last[l]] ,  k +  (int)lb.CzF[leftN_last[l]] ,oppdV[leftN_last[l]]) = grid.Node(i,j,k, leftN_last[l]);

            //second last
            i = grid.nbx + 1;

            for(int l = 0; l<5; l++)
                grid.Cell(i + (int)lb.Cx[leftC_2nd_last[l]], j +  (int)lb.Cy[leftC_2nd_last[l]],  k +  (int)lb.Cz[leftC_2nd_last[l]] ,oppdV[leftC_2nd_last[l]]) = grid.Cell(i,j,k,leftC_2nd_last[l]);

            for(int l = 0; l<5; l++)
                grid.Node(i + (int)lb.Cx[leftN_2nd_last[l]], j +  (int)lb.Cy[leftN_2nd_last[l]],  k +  (int)lb.Cz[leftN_2nd_last[l]] ,oppdV[leftN_2nd_last[l]]) = grid.Node(i,j,k,leftN_2nd_last[l]);

        }
    }

}




template<typename T, typename T1>
void BB_wall_top(Grid_N_C_3D<T> &grid,  lbmD3Q35<T1> &lb, real utop, real ubottom){

    //#--top wall

    int topC_last[]       =   { 
                                dV_ZERO_P2_ZERO ,
                                dV_P1_P1_ZERO   ,dV_M1_P1_ZERO  ,dV_ZERO_P1_P1  ,dV_ZERO_P1_M1  ,
                                dV_P2_P2_P2     ,dV_M2_P2_P2    ,dV_P2_P2_M2    ,dV_M2_P2_M2    ,
                                dV_PH1_PH1_PH1  ,dV_MH1_PH1_PH1 ,dV_PH1_PH1_MH1 , dV_MH1_PH1_MH1

                                };

                                
    int topC_2nd_last[]   =   {
                                dV_ZERO_P2_ZERO ,
                                dV_P2_P2_P2     ,dV_M2_P2_P2    ,dV_P2_P2_M2    ,dV_M2_P2_M2,  
                                };
                                


    int topN_last[]       =   {  dV_ZERO_P2_ZERO ,
                                dV_P1_P1_ZERO   ,dV_M1_P1_ZERO  ,dV_ZERO_P1_P1  ,dV_ZERO_P1_M1  ,
                                dV_P2_P2_P2     ,dV_M2_P2_P2    ,dV_P2_P2_M2    ,dV_M2_P2_M2    ,
                                
                                };

    int topN_2nd_last[]   =   {
                                dV_ZERO_P2_ZERO ,
                                dV_P2_P2_P2     ,dV_M2_P2_P2    ,dV_P2_P2_M2    ,dV_M2_P2_M2,
                                }
                                ;



    for(int i = 0 + grid.nbx ; i <= grid.nex;i++){
        for(int k = 0 + grid.nbz;k <= grid.nez; k++){ 
            
            // last
            int j = grid.ney;

            for(int l = 0; l < 9; l++) 
                grid.Cell(i + (int)lb.Cx[topC_last[l]], j +  (int)lb.Cy[topC_last[l]],  k +  (int)lb.Cz[topC_last[l]] , oppdV[topC_last[l]]) = grid.Cell(i,j,k, topC_last[l]);

            for(int l = 0; l < 9;l++)
                grid.Node(i + (int)lb.Cx[topN_last[l]], j + (int)lb.Cy[topN_last[l]] ,  k +  (int)lb.Cz[topN_last[l]] ,oppdV[topN_last[l]]) = grid.Node(i,j,k, topN_last[l]);

            for(int l = 9; l < 13;l++)
                grid.Node(i + (int)lb.CxC[topC_last[l]], j + (int)lb.CyC[topC_last[l]] ,  k +  (int)lb.CzC[topC_last[l]] ,oppdV[topC_last[l]]) = grid.Cell(i,j,k, topC_last[l]);

            //second last
            j = grid.ney - 1;

            for(int l = 0; l<5; l++)
                grid.Cell(i + (int)lb.Cx[topC_2nd_last[l]], j +  (int)lb.Cy[topC_2nd_last[l]],  k +  (int)lb.Cz[topC_2nd_last[l]] ,oppdV[topC_2nd_last[l]]) = grid.Cell(i,j,k,topC_2nd_last[l]);

            for(int l = 0; l<5; l++)
                grid.Node(i + (int)lb.Cx[topN_2nd_last[l]], j +  (int)lb.Cy[topN_2nd_last[l]],  k +  (int)lb.Cz[topN_2nd_last[l]] ,oppdV[topN_2nd_last[l]]) = grid.Node(i,j,k,topN_2nd_last[l]);

        }
    }

}



template<typename T, typename T1>
void BB_wall_right(Grid_N_C_3D<T> &grid,  lbmD3Q35<T1> &lb, real utop, real ubottom){

    //#--bottom wall

     int rightC_last[]       =   {  dV_P2_ZERO_ZERO ,
                                
                                    dV_P1_P1_ZERO   ,dV_P1_M1_ZERO  ,dV_P1_ZERO_P1  ,dV_P1_ZERO_M1  ,
                                    dV_P2_P2_P2     ,dV_P2_P2_M2    ,dV_P2_M2_P2    ,dV_P2_M2_M2    ,
                                    dV_MH1_PH1_PH1     ,dV_MH1_PH1_MH1    ,dV_MH1_MH1_PH1    ,dV_MH1_MH1_MH1    

                               
                                };



    int rightC_2nd_last[]   =   {
                                dV_P2_ZERO_ZERO ,

                                dV_P2_P2_P2     ,dV_P2_P2_M2    ,dV_P2_M2_P2    ,dV_P2_M2_M2    ,
                                };
                                


    int rightN_last[]       =   { dV_P2_ZERO_ZERO ,
                                
                                dV_P1_P1_ZERO   ,dV_P1_M1_ZERO  ,dV_P1_ZERO_P1  ,dV_P1_ZERO_M1  ,

                                dV_P2_P2_P2     ,dV_P2_P2_M2    ,dV_P2_M2_P2    ,dV_P2_M2_M2    ,
                                
                                };

    int rightN_2nd_last[]   =   {
                                dV_P2_ZERO_ZERO ,

                                dV_P2_P2_P2     ,dV_P2_P2_M2    ,dV_P2_M2_P2    ,dV_P2_M2_M2    ,  
                                }
                                ;






    for(int j = 0 + grid.nby ; j <= grid.ney;j++){
        for(int k = 0 + grid.nbz;k <= grid.nez; k++){ 
            
           
            
            // last
            int i = grid.nex;

            for(int l = 0; l < 9; l++) 
                grid.Cell(i + (int)lb.Cx[rightC_last[l]], j +  (int)lb.Cy[rightC_last[l]],  k +  (int)lb.Cz[rightC_last[l]] , oppdV[rightC_last[l]]) = grid.Cell(i,j,k, rightC_last[l]);

            for(int l = 0; l < 9;l++)
                grid.Node(i + (int)lb.Cx[rightN_last[l]], j + (int)lb.Cy[rightN_last[l]] ,  k +  (int)lb.Cz[rightN_last[l]] ,oppdV[rightN_last[l]]) = grid.Node(i,j,k, rightN_last[l]);

            for(int l = 9; l < 13;l++)
                grid.Node(i + (int)lb.CxC[rightC_last[l]], j + (int)lb.CyC[rightC_last[l]] ,  k +  (int)lb.CzC[rightC_last[l]] ,oppdV[rightC_last[l]]) = grid.Cell(i,j,k, rightC_last[l]);

            //second last
            i = grid.nex - 1;

            for(int l = 0; l<5; l++)
                grid.Cell(i + (int)lb.Cx[rightC_2nd_last[l]], j +  (int)lb.Cy[rightC_2nd_last[l]],  k +  (int)lb.Cz[rightC_2nd_last[l]] ,oppdV[rightC_2nd_last[l]]) = grid.Cell(i,j,k,rightC_2nd_last[l]);

            for(int l = 0; l<5; l++)
                grid.Node(i + (int)lb.Cx[rightN_2nd_last[l]], j +  (int)lb.Cy[rightN_2nd_last[l]],  k +  (int)lb.Cz[rightN_2nd_last[l]] ,oppdV[rightN_2nd_last[l]]) = grid.Node(i,j,k,rightN_2nd_last[l]);

        }
    }


}






template<typename T>
void stationary_correction(Grid_N_C_3D<T> &grid){
    
    //refactor ---------TOP -----------------//

    for(int i = 0 + grid.nbx ; i <= grid.nex;i++){
        for(int k = 0 + grid.nbz;k <= grid.nez; k++){ 
            for(int j = grid.ney - 1; j <= grid.ney; j++){

                double sum = 0;
                for(int l = 1; l<grid.d_v; l++)
                    sum += grid.Node(i,j,k,l);

                grid.Node(i,j,k,dV_ZERO_ZERO_ZERO) = 1 - sum;
                
                sum = 0;
                for(int l = 1; l<grid.d_v; l++)
                    sum += grid.Cell(i,j,k,l);

                grid.Cell(i,j,k,dV_ZERO_ZERO_ZERO) = 1 - sum;

            }
        }
    }
   
   
    //refactor  ------BOTTOM-------------//
    
    for(int i = 0 + grid.nbx ; i <= grid.nex;i++){
        for(int k = 0 + grid.nbz;k <= grid.nez; k++){ 
            for(int j = grid.nby - 1; j <= grid.nby +1 ; j++){

                double sum = 0;
                for(int l = 1; l<grid.d_v; l++)
                    sum += grid.Node(i,j,k,l);

                grid.Node(i,j,k,dV_ZERO_ZERO_ZERO) = 1 - sum;
                
                sum = 0;
                for(int l = 1; l<grid.d_v; l++)
                    sum += grid.Cell(i,j,k,l);

                grid.Cell(i,j,k,dV_ZERO_ZERO_ZERO) = 1 - sum;

            }
        }
    }


}

