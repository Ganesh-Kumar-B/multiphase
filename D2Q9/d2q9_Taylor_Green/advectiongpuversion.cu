#include"assert.h"
#define CUDA_CALL(expr) do {                  \
  cudaError_t err;                            \
  err = expr;                                 \
  assert(err == cudaSuccess);                 \
} while(0)

#define CHECK_ERRORS {\
   cudaError_t code =   cudaGetLastError();\
   if ( cudaSuccess != code )\
    printf( "Error at line %s, %d  %s\n", __FILE__, __LINE__,  cudaGetErrorString(code));\
   }


__global__ void VoxelAdv(double *a,int Nx,int Ny,int Nz,int NX,int NY,int NZ,int NE1,int NE2,int NE3,int NB1,int NB2,int NB3,int LT,int NT,int VecSize,int numFields)
//__host__ void VoxelAdv(double *a,int Nx,int Ny,int Nz,int NE1,int NE2,int NE3,int NB1,int NB2,int NB3,int LT,int NT,int VecSize,int numFields)
{

   int voxelXIndex=threadIdx.x+blockIdx.x*blockDim.x; 
   int voxelYIndex=threadIdx.y+blockIdx.y*blockDim.y; 
   int voxelZIndex=threadIdx.z+blockIdx.z*blockDim.z;
        if(voxelXIndex<NX && voxelYIndex<NY && voxelZIndex<NZ)
        {
           long int index=(((voxelZIndex*Ny+voxelYIndex)*Nx)+voxelXIndex)*NT;
            for(int cellType = 0; cellType < LT; cellType++){
                for(int k = NE3; k >= NB3; k--){
    			for(int j = NE2; j >= NB2; j--){
    				for(int i = NE1; i >= NB1; i--){
    				a[index+cellType*numFields*VecSize+DV_P2_ZERO_ZERO*VecSize+((k*Ny+j)*Nx)+i]=a[index+cellType*numFields*VecSize+DV_P2_ZERO_ZERO*VecSize+(((k)*Ny+(j))*Nx)+(i-2)] ;
                                a[index+cellType*numFields*VecSize+DV_ZERO_P2_ZERO*VecSize+((k*Ny+j)*Nx)+i]=a[index+cellType*numFields*VecSize+DV_ZERO_P2_ZERO*VecSize+(((k)*Ny+(j-2))*Nx)+(i)] ;
                               
                                a[index+cellType*numFields*VecSize+DV_ZERO_ZERO_P2*VecSize+((k*Ny+j)*Nx)+i]=a[index+cellType*numFields*VecSize+DV_ZERO_ZERO_P2*VecSize+(((k-2)*Ny+(j))*Nx)+(i)];                                
                             
                                a[index+cellType*numFields*VecSize+DV_P1_P1_ZERO*VecSize+((k*Ny+j)*Nx)+i]=a[index+cellType*numFields*VecSize+DV_P1_P1_ZERO*VecSize+(((k)*Ny+(j-1))*Nx)+(i-1)] ;
    				a[index+cellType*numFields*VecSize+DV_M1_P1_ZERO*VecSize+((k*Ny+j)*Nx)+i]=a[index+cellType*numFields*VecSize+DV_M1_P1_ZERO*VecSize+(((k)*Ny+(j-1))*Nx)+(i+1)] ;
    				a[index+cellType*numFields*VecSize+DV_P1_ZERO_P1*VecSize+((k*Ny+j)*Nx)+i]=a[index+cellType*numFields*VecSize+DV_P1_ZERO_P1*VecSize+(((k-1)*Ny+(j))*Nx)+(i-1)] ;                               
    				a[index+cellType*numFields*VecSize+DV_ZERO_P1_P1*VecSize+((k*Ny+j)*Nx)+i]=a[index+cellType*numFields*VecSize+DV_ZERO_P1_P1*VecSize+(((k-1)*Ny+(j-1))*Nx)+(i)] ;
    				a[index+cellType*numFields*VecSize+DV_ZERO_M1_P1*VecSize+((k*Ny+j)*Nx)+i]=a[index+cellType*numFields*VecSize+DV_ZERO_M1_P1*VecSize+(((k-1)*Ny+(j+1))*Nx)+(i)] ;
                                a[index+cellType*numFields*VecSize+DV_M1_ZERO_P1*VecSize+((k*Ny+j)*Nx)+i]=a[index+cellType*numFields*VecSize+DV_M1_ZERO_P1*VecSize+(((k-1)*Ny+(j))*Nx)+(i+1)] ;
    			//---------------------------------------------------------------------------------------------------------------------------------//	

    				a[index+cellType*numFields*VecSize+DV_P2_P2_P2*VecSize+((k*Ny+j)*Nx)+i]=a[index+cellType*numFields*VecSize+DV_P2_P2_P2*VecSize+(((k-2)*Ny+(j-2))*Nx)+(i-2)] ;                               

    				a[index+cellType*numFields*VecSize+DV_P2_M2_P2*VecSize+((k*Ny+j)*Nx)+i]=a[index+cellType*numFields*VecSize+DV_P2_M2_P2*VecSize+(((k-2)*Ny+(j+2))*Nx)+(i-2)] ;                               

    				a[index+cellType*numFields*VecSize+DV_M2_P2_P2*VecSize+((k*Ny+j)*Nx)+i]=a[index+cellType*numFields*VecSize+DV_M2_P2_P2*VecSize+(((k-2)*Ny+(j-2))*Nx)+(i+2)] ;                               
    				
                                a[index+cellType*numFields*VecSize+DV_M2_M2_P2*VecSize+((k*Ny+j)*Nx)+i]=a[index+cellType*numFields*VecSize+DV_M2_M2_P2*VecSize+(((k-2)*Ny+(j+2))*Nx)+(i+2)] ;                               
                            }
    			}
    		}


    		for(int k = NB3; k <= NB3; k++){
    			for(int j = NB2; j <= NB2; j++){
    				for(int i = NB1; i <= NE1; i++){
   				a[index+cellType*numFields*VecSize+DV_ZERO_ZERO_M2*VecSize+((k*Ny+j)*Nx)+i]=a[index+cellType*numFields*VecSize+DV_ZERO_ZERO_M2*VecSize+(((k+2)*Ny+(j))*Nx)+(i)] ;
    				a[index+cellType*numFields*VecSize+DV_ZERO_P1_M1*VecSize+((k*Ny+j)*Nx)+i]=a[index+cellType*numFields*VecSize+DV_ZERO_P1_M1*VecSize+(((k+1)*Ny+(j-1))*Nx)+(i)] ;
    				a[index+cellType*numFields*VecSize+DV_ZERO_M1_M1*VecSize+((k*Ny+j)*Nx)+i]=a[index+cellType*numFields*VecSize+DV_ZERO_M1_M1*VecSize+(((k+1)*Ny+(j+1))*Nx)+(i)] ;
                    	        a[index+cellType*numFields*VecSize+DV_P1_ZERO_M1*VecSize+((k*Ny+j)*Nx)+i]=a[index+cellType*numFields*VecSize+DV_ZERO_P1_M1*VecSize+(((k+1)*Ny+(j))*Nx)+(i-1)] ;
                                a[index+cellType*numFields*VecSize+DV_M1_ZERO_M1*VecSize+((k*Ny+j)*Nx)+i]=a[index+cellType*numFields*VecSize+DV_M1_ZERO_M1*VecSize+(((k+1)*Ny+(j))*Nx)+(i+1)] ;
                              //----------------------------------------------------------------------------------------------------------------------	
    				a[index+cellType*numFields*VecSize+DV_M2_M2_M2*VecSize+((k*Ny+j)*Nx)+i]=a[index+cellType*numFields*VecSize+DV_M2_M2_M2*VecSize+(((k+2)*Ny+(j+2))*Nx)+(i+2)] ;
                                a[index+cellType*numFields*VecSize+DV_M2_P2_M2*VecSize+((k*Ny+j)*Nx)+i]=a[index+cellType*numFields*VecSize+DV_M2_P2_M2*VecSize+(((k+2)*Ny+(j-2))*Nx)+(i+2)] ;
                                a[index+cellType*numFields*VecSize+DV_P2_M2_M2*VecSize+((k*Ny+j)*Nx)+i]=a[index+cellType*numFields*VecSize+DV_P2_M2_M2*VecSize+(((k+2)*Ny+(j+2))*Nx)+(i-2)] ;
    				a[index+cellType*numFields*VecSize+DV_P2_P2_M2*VecSize+((k*Ny+j)*Nx)+i]=a[index+cellType*numFields*VecSize+DV_P2_P2_M2*VecSize+(((k+2)*Ny+(j-2))*Nx)+(i-2)] ;
                                a[index+cellType*numFields*VecSize+DV_M2_ZERO_ZERO*VecSize+((k*Ny+j)*Nx)+i]=a[index+cellType*numFields*VecSize+DV_M2_ZERO_ZERO*VecSize+(((k)*Ny+(j))*Nx)+(i+2)] ;
 			        a[index+cellType*numFields*VecSize+DV_ZERO_M2_ZERO*VecSize+((k*Ny+j)*Nx)+i]=a[index+cellType*numFields*VecSize+DV_ZERO_M2_ZERO*VecSize+(((k)*Ny+(j+2))*Nx)+(i)];                                
    	                        a[index+cellType*numFields*VecSize+DV_M1_M1_ZERO*VecSize+((k*Ny+j)*Nx)+i]=a[index+cellType*numFields*VecSize+DV_M1_M1_ZERO*VecSize+(((k)*Ny+(j+1))*Nx)+(i+1)] ;
 			        a[index+cellType*numFields*VecSize+DV_P1_M1_ZERO*VecSize+((k*Ny+j)*Nx)+i]=a[index+cellType*numFields*VecSize+DV_P1_M1_ZERO*VecSize+(((k)*Ny+(j+1))*Nx)+(i-1)] ;


                               }
    			}
    		}

    }
 
int node = 0;
int cell = 1;

        // remove ifs
    for(int k = NE3; k>= NB3; k--){
                for(int j = NE2; j>= NB2; j--){
                        for(int i = NE1; i >= NB1; i--){
                                for(int cellType = 1; cellType >= 0; cellType--){
                                        if(cellType == cell){
                                       //currentVoxel(i, j, k, cellType, lbModel.DV_PH_PH_PH) = currentVoxel(i, j, k, node, lbModel.DV_PH_PH_PH);       
                                                      a[index+cellType*numFields*VecSize+DV_PH_PH_PH*VecSize+((k*Ny+j)*Nx)+i]=a[index+node*numFields*VecSize+DV_PH_PH_PH*VecSize+(((k)*Ny+(j))*Nx)+(i)] ;
                      }
                                        if(cellType == node){
                                           //     currentVoxel(i, j, k, cellType, lbModel.DV_PH_PH_PH) = currentVoxel(i-1, j-1, k-1, cell, lbModel.DV_PH_PH_PH);
                                                a[index+cellType*numFields*VecSize+DV_PH_PH_PH*VecSize+((k*Ny+j)*Nx)+i]=a[index+cell*numFields*VecSize+DV_PH_PH_PH*VecSize+(((k-1)*Ny+(j-1))*Nx)+(i-1)] ;
                                        }
                                }
                        }
                }
        }

   for(int k =NB3; k <= NE3; k++){
        for(int j = NB2; j <= NE2; j++){
                for(int i = NB1; i <=NE1; i++){
                        for(int cellType = 0; cellType <= 1; cellType++){
                                if(cellType == node){
                                      //  currentVoxel(i, j, k, cellType, lbModel.DV_MH_MH_MH) = currentVoxel(i, j, k, cell, lbModel.DV_MH_MH_MH);
                                                a[index+cellType*numFields*VecSize+DV_MH_MH_MH*VecSize+((k*Ny+j)*Nx)+i]=a[index+cell*numFields*VecSize+DV_MH_MH_MH*VecSize+(((k)*Ny+(j))*Nx)+(i)] ;
                                }
                                if(cellType == cell){
                                      //  currentVoxel(i, j, k, cellType, lbModel.DV_MH_MH_MH) = currentVoxel(i+1, j+1, k+1, node, lbModel.DV_MH_MH_MH);
                                      a[index+cellType*numFields*VecSize+DV_MH_MH_MH*VecSize+((k*Ny+j)*Nx)+i]=a[index+node*numFields*VecSize+DV_MH_MH_MH*VecSize+(((k+1)*Ny+(j+1))*Nx)+(i+1)] ;
                                     }
                        }
                }
        }
    }
  for(int k = NE3; k>= NB3; k--){
                for(int j = NB2; j<= NE2; j++){
                        for(int i = NE1; i >= NB1; i--){
                                for(int cellType = 1; cellType >= 0; cellType--){
                                        if(cellType == cell){
                                        // currentVoxel(i, j, k, cellType, lbModel.DV_PH_MH_PH) = currentVoxel(i, j+1, k, node, lbModel.DV_PH_MH_PH);
                                      a[index+cellType*numFields*VecSize+DV_PH_MH_PH*VecSize+((k*Ny+j)*Nx)+i]=a[index+node*numFields*VecSize+DV_PH_MH_PH*VecSize+(((k)*Ny+(j+1))*Nx)+(i)] ;
                                        }
                                        if(cellType == node){
                               // currentVoxel(i, j, k, cellType, lbModel.DV_PH_MH_PH) = currentVoxel(i-1, j, k-1, cell, lbModel.DV_PH_MH_PH);
                                                a[index+cellType*numFields*VecSize+DV_PH_MH_PH*VecSize+((k*Ny+j)*Nx)+i]=a[index+cell*numFields*VecSize+DV_PH_MH_PH*VecSize+(((k-1)*Ny+(j))*Nx)+(i-1)] ;
                                        }
                                }
                        }
                }
        }


        for(int k =NB3; k <=NE3; k++){
        for(int j = NE2; j >= NB2; j--){
                for(int i = NB1; i <= NE1; i++){
                        for(int cellType = 0; cellType <= 1; cellType++){
                                if(cellType == node){
                                // currentVoxel(i, j, k, cellType, lbModel.DV_MH_PH_MH) = currentVoxel(i, j-1, k, cell, lbModel.DV_MH_PH_MH);
                                 a[index+cellType*numFields*VecSize+DV_MH_PH_MH*VecSize+((k*Ny+j)*Nx)+i]=a[index+cell*numFields*VecSize+DV_MH_PH_MH*VecSize+(((k)*Ny+(j-1))*Nx)+(i)] ;
                                    }
                                if(cellType == cell){
                              //  currentVoxel(i, j, k, cellType, lbModel.DV_MH_PH_MH) = currentVoxel(i+1, j, k+1, node, lbModel.DV_MH_PH_MH);
                                 a[index+cellType*numFields*VecSize+DV_MH_PH_MH*VecSize+((k*Ny+j)*Nx)+i]=a[index+node*numFields*VecSize+DV_MH_PH_MH*VecSize+(((k+1)*Ny+(j))*Nx)+(i+1)] ;
                                }
                        }
                }
        }
    }

 for(int k = NE3; k>= NB3; k--){
                for(int j = NB2; j<= NE2; j++){
                        for(int i = NB1; i <= NE1; i++){
                                for(int cellType = 1; cellType >= 0; cellType--){
                                        if(cellType == cell){
//                                currentVoxel(i, j, k, cellType, lbModel.DV_MH_MH_PH) = currentVoxel(i+1, j+1, k, node, lbModel.DV_MH_MH_PH);
                                 a[index+cellType*numFields*VecSize+DV_MH_MH_PH*VecSize+((k*Ny+j)*Nx)+i]=a[index+node*numFields*VecSize+DV_MH_MH_PH*VecSize+(((k)*Ny+(j+1))*Nx)+(i+1)] ;
                                        }
                                        if(cellType == node){
                            //    currentVoxel(i, j, k, cellType, lbModel.DV_MH_MH_PH) = currentVoxel(i, j, k-1, cell, lbModel.DV_MH_MH_PH);                              
                              a[index+cellType*numFields*VecSize+DV_MH_MH_PH*VecSize+((k*Ny+j)*Nx)+i]=a[index+cell*numFields*VecSize+DV_MH_MH_PH*VecSize+(((k-1)*Ny+(j))*Nx)+(i)] ;
                                        }
                                }
                        }
                }
        }


  for(int k = NB3; k <= NE3; k++){
        for(int j = NE2; j >= NB2; j--){
                for(int i =NE1; i >= NB1; i--){
                        for(int cellType = 0; cellType <= 1; cellType++){
                                if(cellType == node){
                               // currentVoxel(i, j, k, cellType, lbModel.DV_PH_PH_MH) = currentVoxel(i-1, j-1, k, cell, lbModel.DV_PH_PH_MH);
                                 a[index+cellType*numFields*VecSize+DV_PH_PH_MH*VecSize+((k*Ny+j)*Nx)+i]=a[index+cell*numFields*VecSize+DV_PH_PH_MH*VecSize+(((k)*Ny+(j-1))*Nx)+(i-1)] ;
                               }
                                if(cellType == cell){
//                                currentVoxel(i, j, k, cellType, lbModel.DV_PH_PH_MH) = currentVoxel(i, j, k+1, node, lbModel.DV_PH_PH_MH);
                                 a[index+cellType*numFields*VecSize+DV_PH_PH_MH*VecSize+((k*Ny+j)*Nx)+i]=a[index+node*numFields*VecSize+DV_PH_PH_MH*VecSize+(((k+1)*Ny+(j))*Nx)+(i)] ;
                                }
                        }
                }
        }
    }


for(int k = NE3; k>= NB3; k--){
                for(int j = NE2; j>= NB2; j--){
                        for(int i = NB1; i <= NE1; i++){
                                for(int cellType = 1; cellType >= 0; cellType--){
                                        if(cellType == cell){
//                                currentVoxel(i, j, k, cellType, lbModel.DV_MH_PH_PH) = currentVoxel(i+1, j, k, node, lbModel.DV_MH_PH_PH);
                                 a[index+cellType*numFields*VecSize+DV_MH_PH_PH*VecSize+((k*Ny+j)*Nx)+i]=a[index+node*numFields*VecSize+DV_MH_PH_PH*VecSize+(((k)*Ny+(j))*Nx)+(i+1)] ;
                                        }
                                        if(cellType == node){
                               // currentVoxel(i, j, k, cellType, lbModel.DV_MH_PH_PH) = currentVoxel(i, j-1, k-1, cell, lbModel.DV_MH_PH_PH);
                                 a[index+cellType*numFields*VecSize+DV_MH_PH_PH*VecSize+((k*Ny+j)*Nx)+i]=a[index+cell*numFields*VecSize+DV_MH_PH_PH*VecSize+(((k-1)*Ny+(j-1))*Nx)+(i)] ;
                                        }
                                }
                        }
                }
        }


 for(int k = NB3; k <= NE3; k++){
        for(int j = NB2; j <= NE2; j++){
                for(int i = NE1; i >=NB1; i--){
                        for(int cellType = 0; cellType <= 1; cellType++){
                                if(cellType == node){
                              //  currentVoxel(i, j, k, cellType, lbModel.DV_PH_MH_MH) = currentVoxel(i-1, j, k, cell, lbModel.DV_PH_MH_MH);
                           a[index+cellType*numFields*VecSize+DV_PH_MH_MH*VecSize+((k*Ny+j)*Nx)+i]=a[index+cell*numFields*VecSize+DV_PH_MH_MH*VecSize+(((k)*Ny+(j))*Nx)+(i-1)] ;
                                }
                                if(cellType == cell){
                              //  currentVoxel(i, j, k, cellType, lbModel.DV_PH_MH_MH) = currentVoxel(i, j+1, k+1, node, lbModel.DV_PH_MH_MH);
                           a[index+cellType*numFields*VecSize+DV_PH_MH_MH*VecSize+((k*Ny+j)*Nx)+i]=a[index+node*numFields*VecSize+DV_PH_MH_MH*VecSize+(((k+1)*Ny+(j+1))*Nx)+(i)] ;
                          }
                        }
                }
        }
    }
}
}



template<typename dataType, int numFields1>
void advectionGrid(lbmD3Q35<dataType> &lbModel, voxelGrid<dataType, numFields1> &mainGrid){ //FLAG
                       int  nx= mainGrid.tempVoxel.n1;
                       int  NX= mainGrid.n1;
                       int  ny= mainGrid.tempVoxel.n2;
                       int  NY= mainGrid.n2;
                       int  nz= mainGrid.tempVoxel.n3;
                       int  NZ= mainGrid.n3;
                       int  NB1= mainGrid.tempVoxel.nB1;
                       int  NB2= mainGrid.tempVoxel.nB2;
                       int  NB3= mainGrid.tempVoxel.nB3;
                       int  NE1= mainGrid.tempVoxel.nE1;
                       int  NE2= mainGrid.tempVoxel.nE2;
                       int  NE3= mainGrid.tempVoxel.nE3;
                       int LT= mainGrid.tempVoxel.LatticeType;     
                       //int Size_of_Grid= NX*NY*NZ;  
                       int Size_of_Voxel= nx*ny*nz*LT*numFields1; 
                     //  printf("%d\n",Size_of_Voxel); 
                       int Size_of_Grid= NX*NY*NZ*Size_of_Voxel;  
                       int VecSize= nx*ny*nz;  
                         
                       int vecsize= NX*NY*NZ;  
                       int i=0;
             //--------------------File For Printing and Testing the Data----------------------------------------------------//
                       FILE *fp,*fb;
                       fp=fopen("VoxelA.txt","w+"); 
                       fb=fopen("VoxelB.txt","w+"); 
            //-------------------------------Declaration of Variables for CUDA------------------------------------------------------------------------------------//
                                dataType *dev_Voxel;   
	        
           //-----------------------------------Memory Allocation for CUDA Variables----------------------------------------------------------------//
                            CUDA_CALL(cudaMalloc((void **) &dev_Voxel,Size_of_Grid*sizeof(dataType))) ;
                       	        
           //----------------Declaration and Memory Allocation of Host array to Store the Entire Grid----------------------------------------------------------//
                                 dataType *host_Voxel;
	                        host_Voxel=(dataType*)calloc(Size_of_Grid,sizeof(dataType));
                           //     printf("Grid=%dVoxel=%d",sizeof(dataType)*Size_of_Grid,sizeof(dataType)*Size_of_Voxel);
            //----------------------------------------Passing Element by Element from Main Grid to Host Array-------------------------------------------------------------------------------//

              for(int voxelZIndex = mainGrid.nB3; voxelZIndex <= mainGrid.nE3; voxelZIndex++){
		for(int voxelYIndex = mainGrid.nB2; voxelYIndex <= mainGrid.nE2; voxelYIndex++){
			for(int voxelXIndex = mainGrid.nB1; voxelXIndex <= mainGrid.nE1; voxelXIndex++){
      	               for(int cellType =0;  cellType<2; cellType++){
	               for(int dv = 0;  dv<35; dv++){
	               for(int k =0 ;  k< mainGrid.tempVoxel.n3; k++){
		        for(int j =0; j < mainGrid.tempVoxel.n2; j++){
			 for(int i =0 ; i < mainGrid.tempVoxel.n1; i++){
                                host_Voxel[(((voxelZIndex*ny+voxelYIndex)*nx)+voxelXIndex)*Size_of_Voxel+cellType*numFields1*VecSize+dv*VecSize+((k*ny+j)*nx)+i]=mainGrid(voxelXIndex,voxelYIndex,voxelZIndex).fData[cellType*numFields1*VecSize+dv*VecSize+((k*ny+j)*nx)+i];
			}
		}
	}
     }
   }
   }
  } 
 }             
 
//-----------------------------------------------Copying Data from the Host to Device Varaible-------------------------------------------------------------------------------//

                     cudaMemcpy(dev_Voxel,host_Voxel,Size_of_Grid*sizeof(dataType),cudaMemcpyHostToDevice);
                       //   cudaDeviceSynchronize();
                                                                                                            
//----------------------------------------------------------------------------------------------------------------------------------------------//
                                    dim3 grid(1,1,NZ);

                                    dim3 block(NX,NY,1);	    
                  VoxelAdv<<<grid,block>>>(dev_Voxel,nx,ny,nz,NX,NY,NZ,NE1,NE2,NE3,NB1,NB2,NB3,LT,Size_of_Voxel,VecSize,numFields1) ;
			CHECK_ERRORS;  
                     cudaMemcpy(host_Voxel,dev_Voxel,Size_of_Grid*sizeof(dataType),cudaMemcpyDeviceToHost);
              for(int voxelZIndex = mainGrid.nB3; voxelZIndex <= mainGrid.nE3; voxelZIndex++){
		for(int voxelYIndex = mainGrid.nB2; voxelYIndex <= mainGrid.nE2; voxelYIndex++){
			for(int voxelXIndex = mainGrid.nB1; voxelXIndex <= mainGrid.nE1; voxelXIndex++){
      	               for(int cellType =0;  cellType<2; cellType++){
	               for(int dv = 0;  dv<35; dv++){
	               for(int k =0 ;  k< mainGrid.tempVoxel.n3; k++){
		        for(int j =0; j < mainGrid.tempVoxel.n2; j++){
			 for(int i =0 ; i < mainGrid.tempVoxel.n1; i++){
                              mainGrid(voxelXIndex,voxelYIndex,voxelZIndex).fData[cellType*numFields1*VecSize+dv*VecSize+((k*ny+j)*nx)+i]=host_Voxel[(((voxelZIndex*ny+voxelYIndex)*nx)+voxelXIndex)*Size_of_Voxel+cellType*numFields1*VecSize+dv*VecSize+((k*ny+j)*nx)+i];
                                            if(i==3 && j==3 && k==3 && cellType==0 && dv==1) 
                     fprintf(fp,"VoxelXIndex=%d\tVoxelYIndex=%d\tVoxelZIndex=%d\tIndex=%d\tMG=%lf\tMG2=%lf\tDV=%lf\n",voxelXIndex,voxelYIndex,voxelZIndex,(((voxelZIndex*ny+voxelYIndex)*nx)+voxelXIndex)*Size_of_Voxel+cellType*numFields1*VecSize+dv*VecSize+((k*ny+j)*nx)+i,mainGrid(voxelXIndex,voxelYIndex,voxelZIndex).fData[cellType*numFields1*VecSize+dv*VecSize+((k*ny+j)*nx)+i],mainGrid(voxelXIndex,voxelYIndex,voxelZIndex)(i,j,k,cellType,dv),host_Voxel[(((voxelZIndex*ny+voxelYIndex)*nx)+voxelXIndex)*Size_of_Voxel+cellType*numFields1*VecSize+dv*VecSize+((k*ny+j)*nx)+i]);     
			}
		}
	}
     }
   }
   }
  } 
 }             
 

 }






                                 // if(voxelXIndex==0 && voxelYIndex==0 && voxelZIndex==1)
                                // { 
                   //  fprintf(fb,"VoxelXIndex=%d\tVoxelYIndex=%d\tVoxelZIndex=%d\tIndex=%d\tHT2=%lf\tMG=%lf\n",voxelXIndex,voxelYIndex,voxelZIndex,(((voxelZIndex*ny+voxelYIndex)*nx)+voxelXIndex)*70000+cellType*numFields1*VecSize+dv*VecSize+((k*ny+j)*nx)+i,host_Voxel[(((voxelZIndex*ny+voxelYIndex)*nx)+voxelXIndex)*70000+cellType*numFields1*VecSize+dv*VecSize+((k*ny+j)*nx)+i],mainGrid(voxelXIndex,voxelYIndex,voxelZIndex).fData[cellType*numFields1*VecSize+dv*VecSize+((k*ny+j)*nx)+i]);    
                     
                                           // if(i==3 && j==3 && k==3 && cellType==0 && dv==1) 
                               // fprintf(fb,"VoxelXIndex=%d\tVoxelYIndex=%d\tVoxelZIndex=%d\tIndex=%d\tMG=%lf\tMG2=%lf\tDV=%lf\n",voxelXIndex,voxelYIndex,voxelZIndex,(((voxelZIndex*ny+voxelYIndex)*nx)+voxelXIndex)*70000+cellType*numFields1*VecSize+dv*VecSize+((k*ny+j)*nx)+i,mainGrid(voxelXIndex,voxelYIndex,voxelZIndex).fData[cellType*numFields1*VecSize+dv*VecSize+((k*ny+j)*nx)+i],mainGrid(voxelXIndex,voxelYIndex,voxelZIndex)(i,j,k,cellType,dv),host_Voxel[(((voxelZIndex*ny+voxelYIndex)*nx)+voxelXIndex)*70000+cellType*numFields1*VecSize+dv*VecSize+((k*ny+j)*nx)+i]);     
                               // } 
             /* for(int voxelZIndex = mainGrid.nB3+1; voxelZIndex <= mainGrid.nB3+1; voxelZIndex++){
		for(int voxelYIndex = mainGrid.nB2; voxelYIndex <= mainGrid.nB2; voxelYIndex++){
			for(int voxelXIndex = mainGrid.nB1; voxelXIndex <= mainGrid.nB1; voxelXIndex++){
      	               for(int cellType =0;  cellType<2; cellType++){
	               for(int dv = 0;  dv<=26; dv++){
	               for(int k =0 ;  k< mainGrid.tempVoxel.n3; k++){
		        for(int j =0; j < mainGrid.tempVoxel.n2; j++){
			 for(int i =0 ; i < mainGrid.tempVoxel.n1; i++){
                     //fprintf(fp,"VoxelXIndex=%d\tVoxelYIndex=%d\tVoxelZIndex=%d\tIndex=%d\tDV=%lf\tHV=%lf\tMG=%lf\n",voxelXIndex,voxelYIndex,voxelZIndex,(((voxelZIndex*ny+voxelYIndex)*nx)+voxelXIndex)*70000+cellType*numFields1*VecSize+dv*VecSize+((k*ny+j)*nx)+i,dev_Voxel[(((voxelZIndex*ny+voxelYIndex)*nx)+voxelXIndex)*70000+cellType*numFields1*VecSize+dv*VecSize+((k*ny+j)*nx)+i],host_Voxel[(((voxelZIndex*ny+voxelYIndex)*nx)+voxelXIndex)*70000+cellType*numFields1*VecSize+dv*VecSize+((k*ny+j)*nx)+i],mainGrid(voxelXIndex,voxelYIndex,voxelZIndex).fData[cellType*numFields1*VecSize+dv*VecSize+((k*ny+j)*nx)+i]);     
                     //fprintf(fp,"VoxelXIndex=%d\tVoxelYIndex=%d\tVoxelZIndex=%d\tIndex=%d\tMG=%lf\n",voxelXIndex,voxelYIndex,voxelZIndex,(((voxelZIndex*ny+voxelYIndex)*nx)+voxelXIndex)*70000+cellType*numFields1*VecSize+dv*VecSize+((k*ny+j)*nx)+i,mainGrid(voxelXIndex,voxelYIndex,voxelZIndex).fData[cellType*numFields1*VecSize+dv*VecSize+((k*ny+j)*nx)+i]);     
                                            if(i==3 && j==3 && k==3 && cellType==0 && dv==1) 
                     fprintf(fp,"VoxelXIndex=%d\tVoxelYIndex=%d\tVoxelZIndex=%d\tIndex=%d\tMG=%lf\tMG2=%lf\tDV=%lf\n",voxelXIndex,voxelYIndex,voxelZIndex,(((voxelZIndex*ny+voxelYIndex)*nx)+voxelXIndex)*70000+cellType*numFields1*VecSize+dv*VecSize+((k*ny+j)*nx)+i,mainGrid(voxelXIndex,voxelYIndex,voxelZIndex).fData[cellType*numFields1*VecSize+dv*VecSize+((k*ny+j)*nx)+i],mainGrid(voxelXIndex,voxelYIndex,voxelZIndex)(i,j,k,cellType,dv),host_Voxel[(((voxelZIndex*ny+voxelYIndex)*nx)+voxelXIndex)*70000+cellType*numFields1*VecSize+dv*VecSize+((k*ny+j)*nx)+i]);     
			}
		}
	}
     }
   }
   }
  } 
 }*/             
                      //    cudaDeviceSynchronize();
//                   VoxelAdv(dev_Voxel,nx,ny,nz,NX,NY,NZ,NE1,NE2,NE3,NB1,NB2,NB3,LT,Size_of_Voxel,VecSize,numFields1);

                      //    printf("Bf=%lf\n",mainGrid(0,0,1)(3,3,3,0,1)); 
			//	advectionVoxel(lbModel,mainGrid(0,0,1),fp);
                        //  printf("Af=%lf\n",mainGrid(0,0,1)(3,3,3,0,1)); 
