/*
  A Simple 1D D1Q3 code for MultiPhase flow
 */
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <math.h>
#include <stddef.h>
/*****Max Size of Array *******/
#define MAX_SIZE  5000 
#define N_DV 3
 
enum direction
{
    ZERO, /* Alias for 0 */
    DX, /* Alias for 1 */
    DMX, /* Alias for -1 */
    D3X, /* Alias for 3 */
    DM3X, /* Alias for -3 */
};

typedef  struct nodeD1Q3
{
    double rho;
    double rhoOld;
    double pNid;
    double FNid;
    double chi;
    double chiPredictor;
    double chiCorrector;
    double muR;
    double muA;
    double dmuA;
    double dmuR;
    double tauM;
    double beta2;
    double betaM;
    double oneM2beta;
    double vel;
    //new foHelmholtz
    double fHelmholtz;
    //define a forcing term
    double Force;
    //define pressure
    double surface;
    double f[N_DV]; 
    double fEq[N_DV]; 
    double f0[N_DV]; 
    double Fh[N_DV]; 
    double df[N_DV];
} nodeD1Q3;

typedef  nodeD1Q3 latticeArr[MAX_SIZE]  ;
 
typedef  struct latticeD1Q3
{
    double weight[N_DV];
    double dvD1Q3[N_DV];
    int sgn[N_DV];
    double T0;
    double T0Inv;
    double c;
    double c2;
} latticeD1Q3;


typedef  struct nonIdealParam  
{
    double a;
    double b;
    double alpha1;
    double kappa;
    double rho0;
} nonIdealParam ;

void getLatticeD1Q3(double c, latticeD1Q3 *myD1Q3);
void printRho(latticeArr lattice, latticeD1Q3 myD1Q3, nonIdealParam myVDW,int nX, double beta,int iX_begin, int iX_end, int time,double c);
void getFeqP2(double fEq[N_DV], latticeD1Q3 myD1Q3,  double rho,  double vel);
void getFeqP(double fEq[N_DV], latticeD1Q3 myD1Q3,  double rho,  double vel);
void getFeqPQuad(double fEq[N_DV], latticeD1Q3 myD1Q3,  double rho,  double vel);

void getFeqThermal(double fEq[N_DV], latticeD1Q3 myD1Q3,  double rho,  double vel, double T);
void getFeqThermal2(double fEq[N_DV], latticeD1Q3 myD1Q3,  double rho,  double vel, double T);
void initializePerturbPeriodic(latticeArr lattice, latticeD1Q3 myD1Q3,  int nX, double vel, double rhoMean, double ampDisturb, double periodDisturb);
void createBoundaryPeriodic(latticeArr lattice, int nX );
void collide(latticeArr myLattice, latticeD1Q3 myD1Q3, nonIdealParam myVDW, int nX, double beta,double tau, int time);
void collideCS(latticeArr myLattice, latticeD1Q3 myD1Q3, nonIdealParam myVDW, int nX, double beta,double tau, int time);
void collideWorkingOld(latticeArr myLattice, latticeD1Q3 myD1Q3, nonIdealParam myVDW, int nX, double beta,double tau, int time);
void collideWorking (latticeArr myLattice, latticeD1Q3 myD1Q3, nonIdealParam myVDW, int nX, double beta,double tau, int time,double c);
void advect(latticeArr lattice, int nX);
void calculateAlpha( latticeArr myLattice, double& alpha, int i, int time);

void initializeGrandPotential(latticeArr myLattice, latticeD1Q3 myD1Q3,  nonIdealParam myVDW, int nX, double vel, double rhoMean, double ampDisturb, double periodDisturb)
{
    int iX;
    double coord, rho, kx;
    for(iX =3; iX<=nX+2;iX++)
    {
        rho = myLattice[iX].f[ZERO] + myLattice[iX].f[DX] + myLattice[iX].f[DMX];
        myLattice[iX].chi = rho*myVDW.b/(1.0-myLattice[iX].rho*myVDW.b) - myVDW.a*rho/myD1Q3.T0 ;
    }
	myLattice[0     ].chi  = myLattice[nX     ].chi ;
    myLattice[1     ].chi  = myLattice[nX + 1 ].chi ;
    myLattice[2     ].chi  = myLattice[nX + 2 ].chi ;
    myLattice[nX + 3].chi  = myLattice[3      ].chi ;
    myLattice[nX + 4].chi  = myLattice[4      ].chi ;
    myLattice[nX + 5].chi  = myLattice[5      ].chi ;
    createBoundaryPeriodic(  myLattice,   nX );
}
 
void relaxGrandPnid(latticeArr myLattice, latticeD1Q3 myD1Q3, nonIdealParam myVDW,  int nX, double dt, double tau, double dx)
{
//     dt *= 0.1;
    double dtBy2tau = dt/tau*0.5;
    double oneMinusdtBy2tau = 1-dtBy2tau;
    double twoBeta = 2*dt/(2*tau+dt);
 
    for( int iX = nX+2  ; iX >=3 ; iX--)
    {
       double chieqnid = myLattice[iX].rho*myVDW.b/(1-myLattice[iX].rho*myVDW.b) - myVDW.a*myLattice[iX].rho;///myD1Q3.T0;
//        myLattice[iX].chi = (1-2*dtBy2tau)*myLattice[iX].chi + 2*dtBy2tau*chieqnid ;
       myLattice[iX].chi = (myLattice[iX].chi + dtBy2tau*chieqnid)/(1+dtBy2tau) ;
//        myLattice[iX].chi = (myLattice[iX].chi*(1-0.5*dtBy2tau) + dtBy2tau*chieqnid)/(1+0.5*dtBy2tau) ;
    }
	myLattice[0     ].chi  = myLattice[nX     ].chi ;
    myLattice[1     ].chi  = myLattice[nX + 1 ].chi ;
    myLattice[2     ].chi  = myLattice[nX + 2 ].chi ;
    myLattice[nX + 3].chi  = myLattice[3      ].chi ;
    myLattice[nX + 4].chi  = myLattice[4      ].chi ;
    myLattice[nX + 5].chi  = myLattice[5      ].chi ;
}
    

void advectGrandPnid (latticeArr myLattice, latticeD1Q3 myD1Q3, int nX, double dt, double dx )
{
//     dt *= 0.1;
    
    // // // Predictor // // //
    for( int iX = nX+2  ; iX >=3 ; iX--)
       myLattice[iX].chiPredictor = myLattice[iX].chi - dt/dx*myLattice[iX].vel*( myLattice[iX+1].chi - myLattice[iX].chi ) ;
	myLattice[0     ].chiPredictor  = myLattice[nX     ].chiPredictor ;
    myLattice[1     ].chiPredictor  = myLattice[nX + 1 ].chiPredictor ;
    myLattice[2     ].chiPredictor  = myLattice[nX + 2 ].chiPredictor ;
    myLattice[nX + 3].chiPredictor  = myLattice[3      ].chiPredictor ;
    myLattice[nX + 4].chiPredictor  = myLattice[4      ].chiPredictor ;
    myLattice[nX + 5].chiPredictor  = myLattice[5      ].chiPredictor ;

    // // // Corrector // // //
    for( int iX = nX+2  ; iX >=3 ; iX--)
       myLattice[iX].chiCorrector = 0.5*( myLattice[iX].chiPredictor + myLattice[iX].chi ) - dt/dx*myLattice[iX].vel*( myLattice[iX].chiPredictor - myLattice[iX-1].chiPredictor ) ;

          
    for( int iX = nX+2  ; iX >=3 ; iX--) {
//        myLattice[iX].pNid = 0.5*(myLattice[iX].chiCorrector + myLattice[iX].chiPredictor);
       myLattice[iX].chi  = myLattice[iX].chiCorrector ;
    }
    
}

main()
{
  
    latticeArr myLattice;
	latticeD1Q3 myD1Q3;
	double kappabar, TbyTc,rho0byrhoc;
	double c, rho,u, beta,  T_critical, tau, dt, dX, rho_critical;
	nonIdealParam myVDW;
	int iX,nX,time, finalTime;
	c = sqrt(3.0);
	getLatticeD1Q3(c, &myD1Q3);
 
    nX = 1000;
	beta = 0.6;
	TbyTc = 0.93;
	rho0byrhoc = 0.9;
	kappabar = 0.0625;//0.0625;
	dX = 1.0/(nX-1.0);
	dt =  dX/c;
	
	tau = (1.0-beta)/beta *dt*0.5;   
	printf("\n beta=%lf kn =%lf \n",beta, tau);
	
    T_critical = myD1Q3.T0/TbyTc ;
    rho_critical = 1.0;

    myVDW.b = 0.521772/(rho_critical);
    myVDW.a = myVDW.b*T_critical/0.377332  ;

	myVDW.kappa = kappabar*myVDW.a*dX*dX;
	myVDW.rho0 = rho0byrhoc*rho_critical;	

    finalTime = 990000    ;

	initializePerturbPeriodic(  myLattice,   myD1Q3,    nX, 0.0, myVDW.rho0, 0.0005,2 );//0.001,2
	initializeGrandPotential (  myLattice,   myD1Q3,  myVDW,  nX, 0.0, myVDW.rho0, 0.001, 2 );//0.001,2
	createBoundaryPeriodic(  myLattice,  nX );
    
	for(time =0; time<= finalTime ; time++)
    {      
       for( int iX = nX+2  ; iX >=3 ; iX--) 
       {
          myLattice[iX].rho = myLattice[iX].f[ZERO] + myLattice[iX].f[DX] + myLattice[iX].f[DMX];
//           myLattice[iX].vel = (myD1Q3.dvD1Q3[DX]*(myLattice[iX].f[DX]-myLattice[iX].f[DMX]))/myLattice[iX].rho ; 
        //   double Force_i    = -(myLattice[iX+1].pNid - myLattice[iX-1].pNid)/(dX*2.0);
        //   myLattice[iX].vel = (myLattice[iX].vel + 0.5*dt*Force_i)/myLattice[iX].rho ; //the velocity is actually not needed this thing is for computing the kappa term gradient
       }
       /* Set Periodicity for computing gradients */
	   myLattice[0     ].vel  = myLattice[nX     ].vel ;
       myLattice[1     ].vel  = myLattice[nX + 1 ].vel ;
       myLattice[2     ].vel  = myLattice[nX + 2 ].vel ;
       myLattice[nX + 3].vel  = myLattice[3      ].vel ;
       myLattice[nX + 4].vel  = myLattice[4      ].vel ;
       myLattice[nX + 5].vel  = myLattice[5      ].vel ;        

       myLattice[0].rho = myLattice[nX].rho ;
       myLattice[1].rho = myLattice[nX + 1].rho;
       myLattice[2].rho = myLattice[nX + 2].rho;
       myLattice[nX + 3].rho = myLattice[3].rho;
       myLattice[nX + 4].rho = myLattice[4].rho;
       myLattice[nX + 5].rho = myLattice[5].rho;
        
    //    relaxGrandPnid  ( myLattice, myD1Q3, myVDW, nX, dt, tau, dX);   
    //    advectGrandPnid ( myLattice, myD1Q3, nX, dt, dX );                    
    //    relaxGrandPnid  ( myLattice, myD1Q3, myVDW, nX, dt, tau, dX);           
        
	   collideWorking (  myLattice,   myD1Q3, myVDW,   nX,  beta,tau,time,c);  
       
	   createBoundaryPeriodic ( myLattice,  nX );
	   advect(   myLattice,       nX);

       if(time % 9900 == 0) {
         printRho(myLattice,myD1Q3,myVDW,nX,beta,3,nX+2,time,c);
       }
/*  ______________________________________________________________   */	  
	}
    printRho(myLattice,myD1Q3,myVDW,nX,beta,3,nX+2,finalTime,c);
}


void   getLatticeD1Q3(double c, latticeD1Q3 *myD1Q3)
{
        double tmp;
		tmp =sqrt(10.0);
		myD1Q3->T0 = c*c/3.0;
		printf("\nTemp = %lf\n", myD1Q3->T0);
		myD1Q3->T0Inv = 1.0/myD1Q3->T0;
		myD1Q3->c = c; 
	   	myD1Q3->c2 = c*c;

		myD1Q3->weight[ZERO] =  2.0/3.0;
		myD1Q3->weight[DX]   =myD1Q3->weight[DMX]  = 1.0/6.0;
		myD1Q3->weight[ZERO] = 1.0-2.0*(myD1Q3->weight[DX]);		


		myD1Q3->dvD1Q3[ZERO] = 0.0;
		myD1Q3->dvD1Q3[DX]   = c;
		myD1Q3->dvD1Q3[DMX]  = -1.0 * c;
		myD1Q3->dvD1Q3[D3X]   = 3.0*c;
		myD1Q3->dvD1Q3[DM3X]  = -3.0 * c;
		myD1Q3->sgn[ZERO] = 0;
		myD1Q3->sgn[DX]   = 1;
		myD1Q3->sgn[DMX]  = -1 ;
		myD1Q3->sgn[D3X]   = 1;
		myD1Q3->sgn[DM3X]  = -1;
		printf("\n w0=%lf w1= %lf  \n", myD1Q3->weight[ZERO], myD1Q3->weight[DX]);
		printf("\n rho=%.16lf  T0=%.16lf\n", myD1Q3->weight[ZERO]+2.0*(myD1Q3->weight[DX]),myD1Q3->T0 );
		return;
}
		
		
		                    /*                 Periodic boundary conditions                      */
/*
 *  Periodicity:
 *
 *         NX-2.5  NX-1.5  NX-0.5   0.5     1.5    2.5    NX-2.5    NX-1.5  NX-0.5   0.5      1.5       2.5
 *          |-------|-------|--------|------|------|--------|--------|-------|--------|--------|---------|  
 * array:   0       1       2        3      4      5       NX       NX+1    NX+2     NX+3     NX+4      NX+5
 *                                           
 *
 *   *********************************************************************************/       
	
void createBoundaryPeriodic(latticeArr myLattice, int nX )
{
    int iX,  dv; 
    /* Periodic BC */
    for (dv = 0; dv < N_DV; dv++) 
    {
        myLattice[0].f[dv] = myLattice[nX     ].f[dv];
        myLattice[1].f[dv] = myLattice[nX + 1 ].f[dv];
        myLattice[2].f[dv] = myLattice[nX + 2 ].f[dv];

        myLattice[nX + 3].f[dv] = myLattice[3].f[dv];
        myLattice[nX + 4].f[dv] = myLattice[4].f[dv];
        myLattice[nX + 5].f[dv] = myLattice[5].f[dv];   
    }
    return;
}
	
void advect(latticeArr myLattice, int nX)
{
    int iX, end;
    end =nX+2;
    for(iX =3; iX<=end;iX++)
    {
        myLattice[iX] .f[DMX] =myLattice[iX+1] .f[DMX]; 
    }
    for(iX =end; iX>=3;iX--)
    {
        myLattice[iX] .f[DX] =myLattice[iX-1] .f[DX]; 
    }
    return;
} 

void initializePerturbPeriodic(latticeArr myLattice, latticeD1Q3 myD1Q3,  int nX, double vel, double rhoMean, double ampDisturb, double periodDisturb)
{
    int iX;
    double coord, rhoIn, kx;
    FILE *fpt;
    fpt = fopen("iC.dat","w");
    /*Initialise the lattice */
    kx = 2*M_PI* periodDisturb;
    for(iX =3; iX<=nX+2;iX++)
    {
        coord= (iX-2.5)/nX;
	    rhoIn =rhoMean + ampDisturb*sin(kx*coord);
        myLattice[iX].rho =rhoIn;
        /*Inititalise f using equilibrium*/
        getFeqPQuad(myLattice[iX].f,myD1Q3, rhoIn, vel);
        fprintf(fpt,"%.10lf  %.10lf \n", coord, rhoIn);
    }
    createBoundaryPeriodic(  myLattice,   nX );
    fclose(fpt);
 }

void printRho(latticeArr myLattice, latticeD1Q3 myD1Q3, nonIdealParam myVDW, int nX, double beta,int iX_begin, int iX_end, int time, double c)
{
    FILE *fpt;
    int iX;
    char *fileName;
    double rho, vel;
    fileName = ( char  * ) malloc( (size_t) (45 * sizeof( char ) ) );
  
    if (!fileName)
    { 
        fprintf(stderr, "Memory allocation failure");
        exit(0);
    }
    sprintf(fileName,"./outputs/rho%d.dat", time);
    
    fpt=fopen(fileName,"w");

    fprintf(fpt,  " # iX \t rho \t vel \t P \t mu \n");
    
    for(iX =iX_begin; iX<iX_end;iX++)
    {
        double P = 0.0, mu = 0.0, dx=0.0, temp;
        dx = 1.0/(nX-1.0);
    
        P = myLattice[iX].rho*myD1Q3.T0/(1.0-myLattice[iX].rho*myVDW.b);
        P += -myVDW.a*myLattice[iX].rho*myLattice[iX].rho;
        P += -myVDW.kappa*myLattice[iX].rho*(myLattice[iX+1].rho+myLattice[iX-1].rho-2.0*myLattice[iX].rho)/(dx*dx);
        P += -0.5*myVDW.kappa*(myLattice[iX+1].rho - myLattice[iX-1].rho)*(myLattice[iX+1].rho - myLattice[iX-1].rho)/(4.0*dx*dx);
        //correction
        //P += -(0.5/beta - 0.25)*temp*temp/myLattice[iX].rho; //term hurts
        //P += 0.25*myVDW.kappa*(myLattice[iX+1].P+myLattice[iX-1].P-2.0*myLattice[iX].P)/(dx*dx); //term may be fine
        //P += -1.0/12.0*myVDW.kappa*(myLattice[iX+1].rho+myLattice[iX-1].rho-2.0*myLattice[iX].rho)/(dx*dx); // this term is good

        mu = myD1Q3.T0*log(myLattice[iX].rho) ;
        mu += -1.5*myD1Q3.T0*log(2.0*M_PI*myD1Q3.T0);
        mu += -2.0*myVDW.a*myLattice[iX].rho -myD1Q3.T0*log(1.0-myLattice[iX].rho*myVDW.b);
        mu += myLattice[iX].rho*myD1Q3.T0*myVDW.b/(1.0-myLattice[iX].rho*myVDW.b);
        mu += - myVDW.kappa*(myLattice[iX-1].rho + myLattice[iX+1].rho - 2.0*myLattice[iX].rho)/(dx*dx);

        fprintf(fpt,  "%d %.7lf  %.7lf %.7lf  %.7lf  %.7lf \n",iX, myLattice[iX].rho, myLattice[iX].vel, P, mu,myLattice[iX].chi);
    }
    fclose(fpt);
}

void collideWorking(latticeArr myLattice, latticeD1Q3 myD1Q3, nonIdealParam myVDW,  int nX, double beta,double tau, int time, double c)
{
    int iX,dv;
    double  fact,quad, quad1, tmp;
    double fEq[N_DV],fc[N_DV], dRho,dPnid;
    double  vel,alpha, lapRho, dmuA,dmuR, mass,sum, df, dt,dx;
    double rhoReduced ,tmp2,fact2,tauM,betaM,rhoRedby4,g;
    mass = 0.0;
    dx = 1.0/(nX-1.0);
    dt = dx / c;

//    	//define pressure
//     for( iX = nX+2 ; iX >=3 ; iX--)
//     {
// //         myLattice[iX].pNid = myVDW.b*myLattice[iX].rho*myLattice[iX].rho*myD1Q3.T0/(1-myVDW.b*myLattice[iX].rho) - myLattice[iX].rho*myLattice[iX].rho*myVDW.a ;
//        myLattice[iX].pNid = myLattice[iX].chi*myLattice[iX].rho*myD1Q3.T0 ;
//         //surface energy
//     } 

    //$pressure from Pnid carnahan starling EOS
        for( iX = nX+2 ; iX >=3 ; iX--)
    {
        double eta_EOS = myLattice[iX].rho * myVDW.b / 4.0;
        myLattice[iX].pNid = (myLattice[iX].rho*myD1Q3.T0  ) *(1+eta_EOS + eta_EOS*eta_EOS - pow(eta_EOS, 3) )/pow((1.0 - eta_EOS), 3)  - myVDW.a * myLattice[iX].rho*myLattice[iX].rho ;
                  // pnid =  (rho**2 b theta) / (1 - rho * b) - a rho**2  
    }

    //$ F_nid  free energy
        for( iX = nX+2 ; iX >=3 ; iX--)
    {   
        double eta_EOS = myLattice[iX].rho * myVDW.b / 4.0;
        myLattice[iX].FNid =  - myVDW.a * myLattice[iX].rho*myLattice[iX].rho   -   myLattice[iX].rho*myD1Q3.T0 *(3.0 * eta_EOS*eta_EOS - 4.0*eta_EOS) / pow(1 - eta_EOS,2)
                            - myVDW.kappa * 0.5*pow((myLattice[iX+1].rho -myLattice[iX-1].rho)/(dx*2.0),2)
        ;
    }

    myLattice[0     ].FNid = myLattice[nX       ].FNid ;
    myLattice[1     ].FNid = myLattice[nX + 1   ].FNid ;
    myLattice[2     ].FNid = myLattice[nX + 2   ].FNid ;
    myLattice[nX + 3].FNid = myLattice[3        ].FNid ;
    myLattice[nX + 4].FNid = myLattice[4        ].FNid ;
    myLattice[nX + 5].FNid = myLattice[5        ].FNid ;    
    
    myLattice[0].pNid = myLattice[nX].pNid ;
    myLattice[1].pNid = myLattice[nX + 1 ].pNid ;
    myLattice[2].pNid = myLattice[nX + 2 ].pNid ;
    myLattice[nX + 3].pNid = myLattice[3].pNid ;
    myLattice[nX + 4].pNid = myLattice[4].pNid ;
    myLattice[nX + 5].pNid = myLattice[5].pNid ;    

   	//calculate chemical potential
   	for( iX = nX+2; iX >=3; iX--)
    {
        //from continuous derivative
        double eta_EOS = myLattice[iX].rho * myVDW.b / 4.0;

        myLattice[iX].muA = myD1Q3.T0*(3.0*pow(eta_EOS,3) - 9.0 *eta_EOS*eta_EOS + 8.0 *eta_EOS )/(pow(1 - eta_EOS,3));
        myLattice[iX].muA -= 2.0*myLattice[iX].rho*myVDW.a ;

        myLattice[iX].muA -= myVDW.kappa*(myLattice[iX-1].rho + myLattice[iX+1].rho - 2.0*myLattice[iX].rho)/(dx*dx) ;

    }
	myLattice[0].muA = myLattice[nX].muA    ;
    myLattice[1].muA = myLattice[nX + 1 ].muA ;
    myLattice[2].muA = myLattice[nX + 2 ].muA ;
    myLattice[nX + 3].muA = myLattice[3].muA ;
    myLattice[nX + 4].muA = myLattice[4].muA ;
    myLattice[nX + 5].muA = myLattice[5].muA ;

    //calcualte the surface terms contribution from kappa
    for( iX = nX+2; iX >=3; iX--)
    {
        myLattice[iX].surface = (myLattice[iX-1].rho + myLattice[iX+1].rho - 2.0*myLattice[iX].rho)/(dx*dx) ;
    }
    myLattice[0].surface = myLattice[nX].surface ;
    myLattice[1].surface = myLattice[nX + 1 ].surface ;
    myLattice[2].surface = myLattice[nX + 2 ].surface ;
    myLattice[nX + 3].surface = myLattice[3].surface ;
    myLattice[nX + 4].surface = myLattice[4].surface ;
    myLattice[nX + 5].surface = myLattice[5].surface ;
    
    
    double eta = 0;

    //calculate force
    for( iX = nX+2; iX >=3 ; iX--)
    {
        mass += myLattice[iX].rho;
		
        myLattice[iX].Force = myLattice[iX].rho*(myLattice[iX+1].muA -myLattice[iX-1].muA)/(dx*2.0);
        // myLattice[iX].Force = (myLattice[iX+1].pNid - myLattice[iX-1].pNid)/(dx*2.0);
        
        // // //# second order
        // double del_muA_second_order =  (1.0/ (2.0 *dx))* ( 
        //                 (eta)*
        //                 (
        //                 (myLattice[iX+1].pNid  + myLattice[iX+1].FNid)* (1.0 /myLattice[iX].rho ) + (1.0/myLattice[iX +1].rho )* (myLattice[iX].pNid  + myLattice[iX].FNid)
        //                 -(myLattice[iX].pNid  + myLattice[iX].FNid) * (1.0 /myLattice[iX - 1].rho ) - (1.0 /myLattice[iX].rho )* (myLattice[iX-1].pNid  + myLattice[iX-1].FNid)
        //                 )
        
        
        //                 -(myVDW.kappa* (myLattice[iX + 1 ].surface - myLattice[iX -1 ].surface) )
        
        //                 );

        // // // //#fourth order
        // double del_muA_fourth_order = (1.0/ (12.0 *dx))*( 
                        
        //                 (1.0 - eta)*  
        //                 (
        //                 8*(myLattice[iX+1].pNid  + myLattice[iX+1].FNid) * (1.0 /myLattice[iX   ].rho ) +
        //                 8*(myLattice[iX  ].pNid  + myLattice[iX  ].FNid) * (1.0 /myLattice[iX +1].rho ) +
        //                   (myLattice[iX-2].pNid  + myLattice[iX-2].FNid) * (1.0 /myLattice[iX   ].rho ) +
        //                   (myLattice[iX  ].pNid  + myLattice[iX  ].FNid) * (1.0 /myLattice[iX -2].rho ) - 
                          
        //                   (myLattice[iX  ].pNid  + myLattice[iX  ].FNid) * (1.0 /myLattice[iX +2].rho ) - 
        //                   (myLattice[iX+2].pNid  + myLattice[iX+2].FNid) * (1.0 /myLattice[iX   ].rho ) - 
        //                 8*(myLattice[iX  ].pNid  + myLattice[iX  ].FNid) * (1.0 /myLattice[iX -1].rho ) - 
        //                 8*(myLattice[iX-1].pNid  + myLattice[iX-1].FNid) * (1.0 /myLattice[iX   ].rho )
        //                 )
        //                 )

        //                 // - (1.0/ (2.0 *dx))* myVDW.kappa* (myLattice[iX + 2].surface + 8 * myLattice[iX +1].surface - 8 *myLattice[iX - 1 ].surface + myLattice[iX -2 ].surface)

        //                 ;
        // myLattice[iX].Force = myLattice[iX].rho*(del_muA_second_order + del_muA_fourth_order); //total force

//   //#fourth order to second order smoothly
        double a = 4.0/3.0; double b = 1.0 - a;
        double del_muA_fourth_order =   (myLattice[iX  ].pNid  + myLattice[iX  ].FNid)*
                                        (   (a/(2.0*dx)) *( (1.0 /myLattice[iX +1].rho ) - (1.0 /myLattice[iX -1].rho )) +   
                                            (b/(4.0*dx)) *( (1.0 /myLattice[iX +2].rho ) - (1.0 /myLattice[iX -2].rho ))  ) +
                                       
                                        (1.0 /myLattice[iX   ].rho )*
                                        (   (a/(2.0*dx)) *( ( myLattice[iX+1].pNid  + myLattice[iX+1].FNid ) - (myLattice[iX-1].pNid  + myLattice[iX-1].FNid)) +   
                                            (b/(4.0*dx)) *( ( myLattice[iX+2].pNid  + myLattice[iX+2].FNid ) - (myLattice[iX-2].pNid  + myLattice[iX-2].FNid))  ) 





                        - (1.0/ (2.0 *dx))*(myVDW.kappa* (myLattice[iX + 1 ].surface - myLattice[iX -1 ].surface) )
                        ;

        myLattice[iX].Force = myLattice[iX].rho*(del_muA_fourth_order); //total force


        // # from direct continuous derivative
        // myLattice[iX].Force = myLattice[iX].rho*(myLattice[iX+1].muA -myLattice[iX-1].muA)/(dx*2.0);
        // std::cout<<myLattice[iX].Force<<std::endl;

        myLattice[iX].Force = -myLattice[iX].Force /myLattice[iX].rho; //force density

    	myLattice[iX].vel = (myD1Q3.dvD1Q3[DX]*(myLattice[iX].f[DX]-myLattice[iX].f[DMX])) ; 
        myLattice[iX].vel = myLattice[iX].vel/myLattice[iX].rho +  0.5*dt*(myLattice[iX].Force);
        
        getFeqPQuad( myLattice[iX].fEq ,   myD1Q3,    myLattice[iX].rho,   myLattice[iX].vel); 
    } 

    for( iX = nX+2  ; iX >=3 ; iX--)
    {
	    alpha = 2.0;
//       	    calculateAlpha( myLattice, alpha, iX, time);
	
	    double chi = 0.0;

        rhoReduced = myVDW.b*myLattice[iX].rho;
        chi = rhoReduced/(1.0-rhoReduced) - myVDW.a*myLattice[iX].rho/myD1Q3.T0;
        tauM = tau;//*myVDW.rho0/myLattice[iX].rho *(1 + rhoReduced*(5.0/8.0 + rhoReduced*(0.2869+ rhoReduced*(0.1103+ 0.0386*rhoReduced))));
//             tauM = tauM/(1.0+chi);
	    betaM = dt/(2.0*tauM + dt);
      
        for(dv = 0; dv < N_DV; dv++)
        {
            myLattice[iX].f[dv]  = myLattice[iX] .f[dv] + 2.0*betaM*(myLattice[iX].fEq[dv] - myLattice[iX] .f[dv]); 
            myLattice[iX].f[dv] += alpha*betaM*tauM*myD1Q3.T0Inv*myLattice[iX].rho*myD1Q3.weight[dv]*myD1Q3.dvD1Q3[dv]*myLattice[iX].Force;
        }      
    }
      
    if(time % 9900 == 0)
    printf("\nAt time =%d mass = %.16lf ", time, mass);
}

void getFeqPQuad(double fEq[N_DV], latticeD1Q3 myD1Q3,  double rho,  double vel)
{
   
    double  fact,quad, quad1, tmp, cubic;
    double  ratio, cubic1, vel1;
    ratio = myD1Q3.T0Inv*myD1Q3.c2;
    quad = myD1Q3.T0Inv *0.5*  vel*vel ;
  
    fEq[ZERO]= myD1Q3.weight[ZERO] *rho*  (1.0  -  quad) ;

    fact =  myD1Q3.T0Inv * vel*myD1Q3.dvD1Q3[DX];
    quad1  = quad *(ratio  -1.0); 
    tmp = rho* myD1Q3.weight[DX];
  
    fEq[DX]= tmp*(1.0 + fact + quad1   ); 
    fEq[DMX]= tmp*(1.0 - fact + quad1 ); 
  
    return;
}


void calculateAlpha( latticeArr myLattice, double& alpha, int i,int time)
{
      alpha = 2.0; 
      
      double a(0.), b(0.), c(0.), ximax;
      double x_i[N_DV];
      
      for(int dv = 0; dv<N_DV; dv++)
      {
	x_i[dv] = myLattice[i].fEq[dv]/myLattice[i].f[dv]-1.0;
      }
	
      ximax = 0.0;
      for(int dv = 0; dv<N_DV; dv++)
      {
	ximax = std::max(fabs(ximax),x_i[dv]);
      }
      
      if(ximax>0.01)
      {
	  for(int dv = 0; dv<N_DV; dv++)
	  {
	      if(x_i[dv]<0.0)
	      {
		  a += myLattice[i].f[dv]*x_i[dv]*x_i[dv]*x_i[dv]*0.5;
		  c += myLattice[i].f[dv]*x_i[dv]*x_i[dv];
	      }
	      
	      else
	      {
		  c += myLattice[i].f[dv]*x_i[dv]*x_i[dv]/(1.0+x_i[dv]);
	      }
	      
	      b += myLattice[i].f[dv]*x_i[dv]*x_i[dv]*0.5;
	  }
	  alpha = (b-sqrt(b*b - 4.0*a*c))/(2.0*a);       
      }
      
//       if( (time==3260) && (i==46)) 
//       {
// 	printf("\n f_i  %.16lf  %.16lf %.16lf   ", myLattice[i].f[DX],myLattice[i].f[ZERO],myLattice[i].f[DMX]);
// 	printf("\n fEq  %.16lf  %.16lf %.16lf   ", myLattice[i].fEq[DX],myLattice[i].fEq[ZERO],myLattice[i].fEq[DMX]);
// 	printf("\n x_i  %.16lf  %.16lf %.16lf   ", x_i[DX],x_i[ZERO],x_i[DMX]);
// 	printf("\n f*x  %.16lf  %.16lf %.16lf   ", myLattice[i].f[DX]*x_i[DX],myLattice[i].f[ZERO]*x_i[ZERO],myLattice[i].f[DMX]*x_i[DMX]);
// 	printf("\n abc  %.16lf  %.16lf %.16lf  alpha %.16lf ", a, b ,c, alpha);
// 	printf("\n %16lf diff in f and feq",myLattice[i].f[ZERO] + myLattice[i].f[DX] + myLattice[i].f[DMX] - myLattice[i].fEq[ZERO] - myLattice[i].fEq[DX] - myLattice[i].fEq[DMX]);
// // 	std::cout << "\n"<< b*b - 4.0*a*c << "\n";
//       }
	for(int dv = 0; dv<N_DV; dv++)
	{
	  if( (myLattice[i].f[dv]<0.0) && (a>=0.0) )
	  {
	    alpha = 2.0;
	    printf(" forcing alpha=2 here %d \t %d \t %.16lf \t %.16lf \n",i,dv, myLattice[i].f[dv],a);
	    break;
	  }
// 	  if (time==3427) std::cout << myLattice[i] .f[dv]<< "  " << myLattice[i] .fEq[dv]<< "  " << alpha << "  "<< dv << "  "<<i << "\n";
	}
     
}
