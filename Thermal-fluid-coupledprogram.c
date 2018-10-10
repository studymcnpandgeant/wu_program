static char help[] = "Nonlinear driven cavity with multigrid in 2d.\n \
  \n\
The 2D driven cavity problem is solved in a velocity-vorticity formulation.\n\
The flow can be driven with the lid or with bouyancy or both:\n\
  -lidvelocity &ltlid&gt, where &ltlid&gt = dimensionless velocity of lid\n\
  -grashof &ltgr&gt, where &ltgr&gt = dimensionless temperature gradent\n\
  -prandtl &ltpr&gt, where &ltpr&gt = dimensionless thermal/momentum diffusity ratio\n\
 -contours : draw contour plots of solution\n\n";
/* in HTML, '&lt' = '<' and '&gt' = '>' */

/*
      See src/ksp/ksp/examples/tutorials/ex45.c
*/

/*T
   Concepts: SNES^solving a system of nonlinear equations (parallel multicomponent example);
   Concepts: DMDA^using distributed arrays;
   Concepts: multicomponent
   Processors: n
T*/


/*F-----------------------------------------------------------------------

    We thank David E. Keyes for contributing the driven cavity discretization within this example code.

    This problem is modeled by the partial differential equation system

\begin{eqnarray}
        - \triangle U - \nabla_y \Omega & = & 0  \\
        - \triangle V + \nabla_x\Omega & = & 0  \\
        - \triangle \Omega + \nabla \cdot ([U*\Omega,V*\Omega]) - GR* \nabla_x T & = & 0  \\
        - \triangle T + PR* \nabla \cdot ([U*T,V*T]) & = & 0
\end{eqnarray}

    in the unit square, which is uniformly discretized in each of x and y in this simple encoding.

    No-slip, rigid-wall Dirichlet conditions are used for $ [U,V]$.
    Dirichlet conditions are used for Omega, based on the definition of
    vorticity: $ \Omega = - \nabla_y U + \nabla_x V$, where along each
    constant coordinate boundary, the tangential derivative is zero.
    Dirichlet conditions are used for T on the left and right walls,
    and insulation homogeneous Neumann conditions are used for T on
    the top and bottom walls.

    A finite difference approximation with the usual 5-point stencil
    is used to discretize the boundary value problem to obtain a
    nonlinear system of equations.  Upwinding is used for the divergence
    (convective) terms and central for the gradient (source) terms.

    The Jacobian can be either
      * formed via finite differencing using coloring (the default), or
      * applied matrix-free via the option -snes_mf
        (for larger grid problems this variant may not converge
        without a preconditioner due to ill-conditioning).

  ------------------------------------------------------------------------F*/

/*
   Include "petscdmda.h" so that we can use distributed arrays (DMDAs).
   Include "petscsnes.h" so that we can use SNES solvers.  Note that this
   file automatically includes:
     petscsys.h       - base PETSc routines   petscvec.h - vectors
     petscmat.h - matrices
     petscis.h     - index sets            petscksp.h - Krylov subspace methods
     petscviewer.h - viewers               petscpc.h  - preconditioners
     petscksp.h   - linear solvers
*/
#if defined(PETSC_APPLE_FRAMEWORK)
#import <PETSc/petscsnes.h>
#import <PETSc/petscdmda.h>
#else
#include <petscsnes.h>
#include <petscdm.h>
#include <petscdmda.h>
#endif

/*
   User-defined routines and data structures
*/
typedef struct {
  PetscScalar p,u,tempc;/* Wu-Three physical field, Pressure , y-velocity , Fluid Temperature*/
} Field;

PetscErrorCode FormFunctionLocal(DMDALocalInfo*,Field**,Field**,void*);

typedef struct {
  PetscReal   lidvelocity,prandtl,grashof;  /* physical parameters */
  PetscReal   PboundaryN,PboundaryS,Tcboundary;                    /* Wu- Pressure, Fluid Temperature Dirichlet Boundary*/
  PetscReal   alpha, AA, g, w, epsilon;     /*Wu- const parameters*/
  PetscBool   draw_contours;                /* flag - 1 indicates drawing contours */
} AppCtx;

extern PetscErrorCode FormInitialGuess(AppCtx*,DM,Vec);
extern PetscErrorCode NonlinearGS(SNES,Vec,Vec,void*);

int main(int argc,char **argv)
{
  AppCtx         user;                /* user-defined work context */
  PetscInt       mx,my,its;
  PetscErrorCode ierr;
  MPI_Comm       comm;
  SNES           snes;
  DM             da;
  Vec            x;

  ierr = PetscInitialize(&argc,&argv,(char*)0,help);if (ierr) return(1);

  PetscFunctionBeginUser;
  comm = PETSC_COMM_WORLD;
  ierr = SNESCreate(comm,&snes);CHKERRQ(ierr);

  /*
      Create distributed array object to manage parallel grid and vectors
      for principal unknowns (x) and governing residuals (f)
  */
  ierr = DMDACreate2d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_STAR,10,10,PETSC_DECIDE,PETSC_DECIDE,3,1,0,0,&da);CHKERRQ(ierr);/*Wu- three dof*/
  ierr = DMSetFromOptions(da);CHKERRQ(ierr);
  ierr = DMSetUp(da);CHKERRQ(ierr);
  ierr = SNESSetDM(snes,(DM)da);CHKERRQ(ierr);
  ierr = SNESSetNGS(snes, NonlinearGS, (void*)&user);CHKERRQ(ierr);

  ierr = DMDAGetInfo(da,0,&mx,&my,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);CHKERRQ(ierr);
  /*
     Problem parameters (velocity of lid, prandtl, and grashof numbers)
  */
  user.lidvelocity = 1.0/(mx*my);
  user.prandtl     = 1.0;
  user.grashof     = 1.0;
  user.PboundaryN   = 5000.0;/*Wu- Define Pressure Boundary */
  user.PboundaryS   = 6500000.0;/*Wu- Define Pressure Boundary */
  user.Tcboundary   = 290;
  user.w           = 0.0134*5*1000/11.8/2.0; /*Wu- Define ficition coefficient*/
  user.alpha       = 41050;
  user.AA          = 187.99;
  user.g           = 9.8;
  user.epsilon     = 0.4465;

  ierr = PetscOptionsGetReal(NULL,NULL,"-lidvelocity",&user.lidvelocity,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-prandtl",&user.prandtl,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-grashof",&user.grashof,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsHasName(NULL,NULL,"-contours",&user.draw_contours);CHKERRQ(ierr);

  /*Wu- Set DMDA Field Name, Three Dofs*/
  ierr = DMDASetFieldName(da,0,"Pressure");CHKERRQ(ierr);
  ierr = DMDASetFieldName(da,1,"Y-Velocity");CHKERRQ(ierr);
  ierr = DMDASetFieldName(da,2,"Fluid_temperature");CHKERRQ(ierr);
  //ierr = DMDASetFieldName(da,3,"temperature");CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create user context, set problem data, create vector data structures.
     Also, compute the initial guess.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create nonlinear solver context

     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = DMSetApplicationContext(da,&user);CHKERRQ(ierr);
  ierr = DMDASNESSetFunctionLocal(da,INSERT_VALUES,(PetscErrorCode (*)(DMDALocalInfo*,void*,void*,void*))FormFunctionLocal,&user);CHKERRQ(ierr);
  ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);
  ierr = PetscPrintf(comm,"lid velocity = %g, prandtl # = %g, grashof # = %g\n",(double)user.lidvelocity,(double)user.prandtl,(double)user.grashof);CHKERRQ(ierr);


  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Solve the nonlinear system
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = DMCreateGlobalVector(da,&x);CHKERRQ(ierr);
  ierr = FormInitialGuess(&user,da,x);CHKERRQ(ierr);

  ierr = SNESSolve(snes,NULL,x);CHKERRQ(ierr);

  ierr = SNESGetIterationNumber(snes,&its);CHKERRQ(ierr);
  ierr = PetscPrintf(comm,"Number of SNES iterations = %D\n", its);CHKERRQ(ierr);

  /*
     Visualize solution
  */
  if (user.draw_contours) {
    ierr = VecView(x,PETSC_VIEWER_DRAW_WORLD);CHKERRQ(ierr);
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Free work space.  All PETSc objects should be destroyed when they
     are no longer needed.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = VecDestroy(&x);CHKERRQ(ierr);
  ierr = DMDestroy(&da);CHKERRQ(ierr);
  ierr = SNESDestroy(&snes);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return ierr;
}

/* ------------------------------------------------------------------- */

/*
   FormInitialGuess - Forms initial approximation.

   Input Parameters:
   user - user-defined application context
   X - vector

   Output Parameter:
   X - vector
*/
PetscErrorCode FormInitialGuess(AppCtx *user,DM da,Vec X)
{
  PetscInt       i,j,mx,xs,ys,xm,ym;
  PetscErrorCode ierr;
  PetscReal      grashof,dx;
  PetscReal      PboundaryN;/*Wu - Fetch Pboundary*/
  Field          **x;

  PetscFunctionBeginUser;
  grashof = user->grashof;
  PboundaryN = user->PboundaryN; /*Wu -Fetch Pboundary*/

  ierr = DMDAGetInfo(da,0,&mx,0,0,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
  dx   = 1.0/(mx-1);

  /*
     Get local grid boundaries (for 2-dimensional DMDA):
       xs, ys   - starting grid indices (no ghost points)
       xm, ym   - widths of local grid (no ghost points)
  */
  ierr = DMDAGetCorners(da,&xs,&ys,NULL,&xm,&ym,NULL);CHKERRQ(ierr);

  /*
     Get a pointer to vector data.
       - For default PETSc vectors, VecGetArray() returns a pointer to
         the data array.  Otherwise, the routine is implementation dependent.
       - You MUST call VecRestoreArray() when you no longer need access to
         the array.
  */
  ierr = DMDAVecGetArray(da,X,&x);CHKERRQ(ierr);

  /*
     Compute initial guess over the locally owned part of the grid
     Initial condition is motionless fluid and equilibrium temperature
  */
  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      x[j][i].p     = 20000.0;
      x[j][i].u     = 5.0;
      x[j][i].tempc = 290.0;
      //x[j][i].temp  = (grashof>0)*i*dx;
    }
  }

  /*
     Restore vector
  */
  ierr = DMDAVecRestoreArray(da,X,&x);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PETSC_STATIC_INLINE PetscScalar LiquDen(PetscScalar Tc)
{
  return -0.00836 * Tc * Tc +2.932 * Tc + 598.24;
}
PETSC_STATIC_INLINE PetscScalar LiquCond(PetscScalar Tc)
{
  return -8.82e-6 * Tc *Tc + 0.0034 * Tc + 0.3351;
}
PETSC_STATIC_INLINE PetscScalar LiquHetM(PetscScalar Tc)
{
  return 0.25393 * Tc * Tc - 126.86 * Tc + 20685;
}

PetscErrorCode FormFunctionLocal(DMDALocalInfo *info,Field **x,Field **f,void *ptr)
{
  AppCtx         *user = (AppCtx*)ptr;
  PetscErrorCode ierr;
  PetscInt       xints,xinte,yints,yinte,i,j;
  PetscInt       xmax,ymax;/*Wu-defined*/
  PetscReal      hx,hy,dhx,dhy,hxdhy,hydhx;
  PetscReal      grashof,prandtl,lid;
  PetscReal      PboundaryN,PboundaryS,w,Tcboundary;/*Wu-Fetch Pressure Boundary*/
  PetscReal      alpha,AA,g,epsilon;
  PetscScalar    center,east,west,s,n,uxx,uyy,vx,vy,avx,avy,vxp,vxm,vyp,vym;
  PetscScalar    a,b,c,d,e,rightside; /*Coefficient Matrix*/
  PetscScalar    rho_E,rho_W,rho_S,rho_N,rho_C;/*Material coefficient*/
  PetscScalar    Cp_E,Cp_W,Cp_S,Cp_N,Cp_C;
  PetscScalar    lambda_E,lambda_N,lambda_W,lambda_S,lambda_C;

  PetscFunctionBeginUser;
  grashof = user->grashof;
  prandtl = user->prandtl;
  lid     = user->lidvelocity;
  PboundaryN = user->PboundaryN;/*Wu-Fetch Pressure Boundary*/
  PboundaryS = user->PboundaryS;/*Wu-Fetch Pressure Boundary*/
  Tcboundary = user->Tcboundary;
  alpha   = user->alpha;
  w       = user->w;/*Wu-Fetch friction coefficient*/
  AA      = user->AA;
  g       = user->g;
  epsilon = user->epsilon;

  /*
     Define mesh intervals ratios for uniform grid.

     Note: FD formulae below are normalized by multiplying through by
     local volume element (i.e. hx*hy) to obtain coefficients O(1) in two dimensions.


  */
  dhx   = (PetscReal)(info->mx-1);  dhy = (PetscReal)(info->my-1);
  /*Interchange the hx and the hy */
  hx    = 3.0/dhy;                   hy = 2.0/dhx;
  hxdhy = hx*dhy;                 hydhx = hy*dhx;

  xints = info->xs; xinte = info->xs+info->xm; yints = info->ys; yinte = info->ys+info->ym;
  xmax = xinte -1; ymax = info->my -1;

  /* Test whether we are on the bottom edge of the global array */
  if (yints == 0) {
    j     = 0;
    yints = yints + 1;
    /* bottom edge */
    for (i=info->xs; i<info->xs+info->xm; i++) {
      /*fluid field*/
      rightside = Tcboundary;
      b         = 0.0;
      c         = 0.0;
      d         = 0.0;
      e         = 0.0;
      a         = 1.0;
      f[j][i].tempc = a * x[j][i].tempc - rightside;

      /*velocity field*/
      center = x[j][i].tempc;
      rho_C = LiquDen( center );
      rightside = epsilon * rho_C * g -epsilon * ( x[j+1][i].p - x[j][i].p ) / hx;
      b         = 0.0;
      c         = 0.0;
      d         = 0.0;
      e         = 0.0;
      a         = w * rho_C;
      f[j][i].u = a * x[j][i].u - rightside;

      /*pressure field*/
      center = x[j][i].tempc;
      n = x[j+1][i].tempc;
      rho_C = LiquDen( center );
      rho_N = LiquDen( n );
      rightside = ((rho_N * g / w)-(rho_C * g / w )) / hx - ( epsilon * rho_C * g - w * rho_C * 5.0) * hx / epsilon /(w * hx * hx);
      b         = 0.0;
      c         = 1 / ( w * hx * hx );
      d         = 0.0;
      e         = 0.0;
      a         = c;
      f[j][i].p = a*x[j][i].p - c*x[j+1][i].p -rightside;
    }
  }

  /* Test whether we are on the top edge of the global array */
  if (yinte == info->my) {
    j     = info->my - 1;
    yinte = yinte - 1;
    /* top edge there is no j+1 */
    for (i=info->xs; i<info->xs+info->xm; i++) {
      /*Up-left*/
      if ( i == 0 ){
           /*fluid field*/
      center = x[j][i].tempc;
      lambda_C = LiquCond( center );
      rho_C = LiquDen( center );
      Cp_C = LiquHetM( center );
      s = x[j-1][i].tempc;
      lambda_S = LiquCond( s );
      rho_S = LiquDen( s );
      Cp_S = LiquHetM( s );
      east = x[j][i+1].tempc;
      lambda_E = LiquCond( east );
      rightside = alpha * AA * 300.0;
      b         = 1.0 / ( (hx*hx) / (2*epsilon*lambda_C) + (hx*hx) / (2*epsilon*lambda_S) ) + ( epsilon * rho_S * Cp_S * x[j-1][i].u ) / hx;
      c         = 0.0;
      d         = 0.0;/*Left corner there is no i-1 */
      e         = 1.0 /( (hy*hy) / (2*epsilon*lambda_C) + (hy*hy) / (2*epsilon*lambda_E) );
      a         = 1.0 /( (hx*hx) / (2*epsilon*lambda_C) + (hx*hx) / (2*epsilon*lambda_S) ) + e + ( epsilon * rho_C * Cp_C * x[j][i].u ) / hx + alpha * AA;
      f[j][i].tempc = a * x[j][i].tempc - b * x[j-1][i].tempc - e * x[j][i+1].tempc - rightside;

      /*velocity field*/
      center = x[j][i].tempc;
      rho_C = LiquDen( center );
      rightside = epsilon * rho_C * g -epsilon * ( x[j][i].p - x[j-1][i].p ) / hx;
      b         = 0.0;
      c         = 0.0;
      d         = 0.0;
      e         = 0.0;
      a         = w * rho_C;
      f[j][i].u = a * x[j][i].u - rightside;

      /*pressure field*/
      center = x[j][i].tempc;
      rho_C = LiquDen( center );
      s = x[j-1][i].tempc;
      rho_S = LiquDen( s );
      rightside = ((rho_C * g / w)-(rho_S * g / w )) / hx - ( PboundaryN ) /(w * hx * hx);
      b         = 1.0 / ( w * hx * hx );
      c         = 0.0;
      d         = 0.0;
      e         = 0.0;
      a         = b + 1.0 / ( w * hx * hx );
      f[j][i].p = a*x[j][i].p - b * x[j-1][i].p -rightside;
      }
      /*Up-right there is no i+1*/
      else if ( i == info->xs+info->xm-1 ){
           /*fluid field*/
      center = x[j][i].tempc;
      lambda_C = LiquCond( center );
      rho_C = LiquDen( center );
      Cp_C = LiquHetM( center );
      s = x[j-1][i].tempc;
      lambda_S = LiquCond( s );
      rho_S = LiquDen( s );
      Cp_S = LiquHetM( s );
      west = x[j][i-1].tempc;
      lambda_W = LiquCond( west );
      rightside = alpha * AA * 300.0;
      b         = 1.0 /( (hx*hx) / (2*epsilon*lambda_C) + (hx*hx) / (2*epsilon*lambda_S) ) + ( epsilon * rho_S * Cp_S * x[j-1][i].u ) / hx;
      c         = 0.0;
      d         = 1.0 /( (hy*hy) / (2*epsilon*lambda_C) + (hy*hy) / (2*epsilon*lambda_W) );
      e         = 0.0;
      a         = 1.0 /( (hx*hx) / (2*epsilon*lambda_C) + (hx*hx) / (2*epsilon*lambda_S) ) + d + ( epsilon * rho_C * Cp_C * x[j][i].u ) / hx + alpha * AA;
      f[j][i].tempc = a * x[j][i].tempc - b * x[j-1][i].tempc - d * x[j][i-1].tempc - rightside;

      /*velocity field*/
      center = x[j][i].tempc;
      rho_C = LiquDen( center );
      rightside = epsilon * rho_C * g -epsilon * ( x[j][i].p - x[j-1][i].p ) / hx;
      b         = 0.0;
      c         = 0.0;
      d         = 0.0;
      e         = 0.0;
      a         = w * rho_C;
      f[j][i].u = a * x[j][i].u - rightside;

      /*pressure field*/
      center = x[j][i].tempc;
      rho_C = LiquDen( center );
      s = x[j-1][i].tempc;
      rho_S = LiquDen( s );
      rightside = ((rho_C * g / w)-(rho_S * g / w )) / hx - ( PboundaryN ) /(w * hx * hx);
      b         = 1.0 / ( w * hx * hx );
      c         = 0.0;
      d         = 0.0;
      e         = 0.0;
      a         = b + 1.0 / ( w * hx * hx );
      f[j][i].p = a*x[j][i].p - b * x[j-1][i].p -rightside;
      }
      /*Top boundary*/
      else{
            /*fluid field*/
      center = x[j][i].tempc;
      lambda_C = LiquCond( center );
      rho_C = LiquDen( center );
      Cp_C = LiquHetM( center );
      s = x[j-1][i].tempc;
      lambda_S = LiquCond( s );
      rho_S = LiquDen( s );
      Cp_S = LiquHetM( s );
      east = x[j][i+1].tempc;
      lambda_E = LiquCond( east );
      west = x[j][i-1].tempc;
      lambda_W = LiquCond( west );
      rightside = alpha * AA * 300.0;
      b         = 1/( (hx*hx) / (2*epsilon*lambda_C) + (hx*hx) / (2*epsilon*lambda_S) ) + ( epsilon * rho_S * Cp_S * x[j-1][i].u ) / hx;
      c         = 0.0;
      d         = 1/( (hy*hy) / (2*epsilon*lambda_C) + (hy*hy) / (2*epsilon*lambda_W) );
      e         = 1/( (hy*hy) / (2*epsilon*lambda_C) + (hy*hy) / (2*epsilon*lambda_E) );
      a         = 1/( (hx*hx) / (2*epsilon*lambda_C) + (hx*hx) / (2*epsilon*lambda_S) ) + d + e + ( epsilon * rho_C * Cp_C * x[j][i].u ) / hx + alpha * AA;
      f[j][i].tempc = a * x[j][i].tempc - b * x[j-1][i].tempc - d * x[j][i-1].tempc - e * x[j][i+1].tempc - rightside;

      /*velocity field*/
      center = x[j][i].tempc;
      rho_C = LiquDen( center );
      rightside = epsilon * rho_C * g -epsilon * ( x[j][i].p - x[j-1][i].p ) / hx;
      b         = 0.0;
      c         = 0.0;
      d         = 0.0;
      e         = 0.0;
      a         = w * rho_C;
      f[j][i].u = a * x[j][i].u - rightside;

      /*pressure field*/
      center = x[j][i].tempc;
      rho_C = LiquDen( center );
      s = x[j-1][i].tempc;
      rho_S = LiquDen( s );
      rightside = ((rho_C * g / w)-(rho_S * g / w )) / hx - ( PboundaryN ) /(w * hx * hx);
      b         = 1.0 / ( w * hx * hx );
      c         = 0.0;
      d         = 0.0;
      e         = 0.0;
      a         = b + 1.0 / ( w * hx * hx );
      f[j][i].p = a*x[j][i].p - b * x[j-1][i].p -rightside;
          }
    }
  }

  /* Test whether we are on the left edge of the global array */
  if (xints == 0) {
    i     = 0;
    xints = xints + 1;
    /* left edge  there is no i-1*/
    for (j=info->ys; j<info->ys+info->ym; j++) {
      /*left-bottom  there is no j-1*/
      if ( j == 0 ){
           /*fluid field*/
      rightside = Tcboundary;
      b         = 0.0;
      c         = 0.0;
      d         = 0.0;
      e         = 0.0;
      a         = 1.0;
      f[j][i].tempc = a * x[j][i].tempc - rightside;

      /*velocity field*/
      center = x[j][i].tempc;
      rho_C = LiquDen( center );
      rightside = epsilon * rho_C * g -epsilon * ( x[j+1][i].p - x[j][i].p ) / hx;
      b         = 0.0;
      c         = 0.0;
      d         = 0.0;
      e         = 0.0;
      a         = w * rho_C;
      f[j][i].u = a * x[j][i].u - rightside;

      /*pressure field*/
      center = x[j][i].tempc;
      n = x[j+1][i].tempc;
      rho_C = LiquDen( center );
      rho_N = LiquDen( n );
      rightside = ((rho_N * g / w)-(rho_C * g / w )) / hx - ( epsilon * rho_C * g - w * rho_C * 5.0) * hx / epsilon /(w * hx * hx);
      b         = 0.0;
      c         = 1.0 / ( w * hx * hx );
      d         = 0.0;
      e         = 0.0;
      a         = c;
      f[j][i].p = a*x[j][i].p - c*x[j+1][i].p -rightside;
      }
      /*Up-left there is no j+1*/
      else if ( j == info->ys+info->ym-1 ){
           /*fluid field*/
      center = x[j][i].tempc;
      lambda_C = LiquCond( center );
      rho_C = LiquDen( center );
      Cp_C = LiquHetM( center );
      s = x[j-1][i].tempc;
      lambda_S = LiquCond( s );
      rho_S = LiquDen( s );
      Cp_S = LiquHetM( s );
      east = x[j][i+1].tempc;
      lambda_E = LiquCond( east );
      rightside = alpha * AA * 300.0;
      b         = 1.0 /( (hx*hx) / (2*epsilon*lambda_C) + (hx*hx) / (2*epsilon*lambda_S) ) + ( epsilon * rho_S * Cp_S * x[j-1][i].u ) / hx;
      c         = 0.0;
      d         = 0.0;/*Left corner there is no i-1 */
      e         = 1.0 /( (hy*hy) / (2*epsilon*lambda_C) + (hy*hy) / (2*epsilon*lambda_E) );
      a         = 1.0 /( (hx*hx) / (2*epsilon*lambda_C) + (hx*hx) / (2*epsilon*lambda_S) ) + e + ( epsilon * rho_C * Cp_C * x[j][i].u ) / hx + alpha * AA;
      f[j][i].tempc = a * x[j][i].tempc - b * x[j-1][i].tempc - e * x[j][i+1].tempc - rightside;

      /*velocity field*/
      center = x[j][i].tempc;
      rho_C = LiquDen( center );
      rightside = epsilon * rho_C * g -epsilon * ( x[j][i].p - x[j-1][i].p ) / hx;
      b         = 0.0;
      c         = 0.0;
      d         = 0.0;
      e         = 0.0;
      a         = w * rho_C;
      f[j][i].u = a * x[j][i].u - rightside;

      /*pressure field*/
      center = x[j][i].tempc;
      rho_C = LiquDen( center );
      s = x[j-1][i].tempc;
      rho_S = LiquDen( s );
      rightside = ((rho_C * g / w)-(rho_S * g / w )) / hx - ( PboundaryN ) /(w * hx * hx);
      b         = 1.0 / ( w * hx * hx );
      c         = 0.0;
      d         = 0.0;
      e         = 0.0;
      a         = b + 1.0 / ( w * hx * hx );
      f[j][i].p = a*x[j][i].p - b * x[j-1][i].p -rightside;
      }
      /*Left boundary*/
      else {
            /*fluid field*/
      center = x[j][i].tempc;
      lambda_C = LiquCond( center );
      rho_C = LiquDen( center );
      Cp_C = LiquHetM( center );
      s = x[j-1][i].tempc;
      lambda_S = LiquCond( s );
      rho_S = LiquDen( s );
      Cp_S = LiquHetM( s );
      east = x[j][i+1].tempc;
      lambda_E = LiquCond( east );
      n= x[j+1][i].tempc;
      lambda_N = LiquCond( n );
      rightside = alpha * AA * 300.0;
      b         = 1.0 /( (hx*hx) / (2*epsilon*lambda_C) + (hx*hx) / (2*epsilon*lambda_S) ) + ( epsilon * rho_S * Cp_S * x[j-1][i].u ) / hx;
      c         = 1.0 /( (hx*hx) / (2*epsilon*lambda_C) + (hx*hx) / (2*epsilon*lambda_N) );
      d         = 0.0;
      e         = 1.0 /( (hy*hy) / (2*epsilon*lambda_C) + (hy*hy) / (2*epsilon*lambda_E) );
      a         = 1.0 /( (hx*hx) / (2*epsilon*lambda_C) + (hx*hx) / (2*epsilon*lambda_S) ) + c + e + ( epsilon * rho_C * Cp_C * x[j][i].u ) / hx + alpha * AA;
      f[j][i].tempc = a * x[j][i].tempc - b * x[j-1][i].tempc - c * x[j+1][i].tempc - e * x[j][i+1].tempc - rightside;

      /*velocity field*/
      center = x[j][i].tempc;
      rho_C = LiquDen( center );
      rightside = epsilon * rho_C * g -epsilon * ( x[j+1][i].p - x[j][i].p ) / hx;
      b         = 0.0;
      c         = 0.0;
      d         = 0.0;
      e         = 0.0;
      a         = w * rho_C;
      f[j][i].u = a * x[j][i].u - rightside;

      /*pressure field*/
      center = x[j][i].tempc;
      rho_C = LiquDen( center );
      s = x[j-1][i].tempc;
      rho_S = LiquDen( s );
      rightside = ((rho_C * g / w)-(rho_S * g / w )) / hx;
      b         = 1.0 / ( w * hx * hx );
      c         = 1.0 / ( w * hx * hx );
      d         = 0.0;
      e         = 0.0;
      a         = b + c;
      f[j][i].p = a*x[j][i].p - b * x[j-1][i].p - c * x[j+1][i].p -rightside;
      }
    }
  }

  /* Test whether we are on the right edge of the global array */
  if (xinte == info->mx) {
    i     = info->mx - 1;
    xinte = xinte - 1;
    /* right edge there is no i+1*/
    for (j=info->ys; j<info->ys+info->ym; j++) {
      /*right bottom there is no j-1*/
      if ( j == 0 ){
           /*fluid field*/
      rightside = Tcboundary;
      b         = 0.0;
      c         = 0.0;
      d         = 0.0;
      e         = 0.0;
      a         = 1.0;
      f[j][i].tempc = a * x[j][i].tempc - rightside;

      /*velocity field*/
      center = x[j][i].tempc;
      rho_C = LiquDen( center );
      rightside = epsilon * rho_C * g -epsilon * ( x[j+1][i].p - x[j][i].p ) / hx;
      b         = 0.0;
      c         = 0.0;
      d         = 0.0;
      e         = 0.0;
      a         = w * rho_C;
      f[j][i].u = a * x[j][i].u - rightside;

      /*pressure field*/
      center = x[j][i].tempc;
      n = x[j+1][i].tempc;
      rho_C = LiquDen( center );
      rho_N = LiquDen( n );
      rightside = ((rho_N * g / w)-(rho_C * g / w )) / hx - ( epsilon * rho_C * g - w * rho_C * 5.0) * hx / epsilon /(w * hx * hx);
      b         = 0.0;
      c         = 1.0 / ( w * hx * hx );
      d         = 0.0;
      e         = 0.0;
      a         = c;
      f[j][i].p = a*x[j][i].p - c*x[j+1][i].p -rightside;
      }
      /*Up-right there is no j+1*/
      else if ( j == info->ys+info->ym-1 ){
             /*fluid field*/
      center = x[j][i].tempc;
      lambda_C = LiquCond( center );
      rho_C = LiquDen( center );
      Cp_C = LiquHetM( center );
      s = x[j-1][i].tempc;
      lambda_S = LiquCond( s );
      rho_S = LiquDen( s );
      Cp_S = LiquHetM( s );
      west = x[j][i-1].tempc;
      lambda_W = LiquCond( west );
      rightside = alpha * AA * 300.0;
      b         = 1.0 /( (hx*hx) / (2*epsilon*lambda_C) + (hx*hx) / (2*epsilon*lambda_S) ) + ( epsilon * rho_S * Cp_S * x[j-1][i].u ) / hx;
      c         = 0.0;
      d         = 1.0 /( (hy*hy) / (2*epsilon*lambda_C) + (hy*hy) / (2*epsilon*lambda_W) );
      e         = 0.0;
      a         = 1.0 /( (hx*hx) / (2*epsilon*lambda_C) + (hx*hx) / (2*epsilon*lambda_S) ) + d + ( epsilon * rho_C * Cp_C * x[j][i].u ) / hx + alpha * AA;
      f[j][i].tempc = a * x[j][i].tempc - b * x[j-1][i].tempc - d * x[j][i-1].tempc - rightside;

      /*velocity field*/
      center = x[j][i].tempc;
      rho_C = LiquDen( center );
      rightside = epsilon * rho_C * g -epsilon * ( x[j][i].p - x[j-1][i].p ) / hx;
      b         = 0.0;
      c         = 0.0;
      d         = 0.0;
      e         = 0.0;
      a         = w * rho_C;
      f[j][i].u = a * x[j][i].u - rightside;

      /*pressure field*/
      center = x[j][i].tempc;
      rho_C = LiquDen( center );
      s = x[j-1][i].tempc;
      rho_S = LiquDen( s );
      rightside = ((rho_C * g / w)-(rho_S * g / w )) / hx - ( PboundaryN ) /(w * hx * hx);
      b         = 1.0 / ( w * hx * hx );
      c         = 0.0;
      d         = 0.0;
      e         = 0.0;
      a         = b + 1.0 / ( w * hx * hx );
      f[j][i].p = a*x[j][i].p - b * x[j-1][i].p -rightside;
      }
      /*right*/
      else {
           /*fluid field*/
      center = x[j][i].tempc;
      lambda_C = LiquCond( center );
      rho_C = LiquDen( center );
      Cp_C = LiquHetM( center );
      s = x[j-1][i].tempc;
      lambda_S = LiquCond( s );
      rho_S = LiquDen( s );
      Cp_S = LiquHetM( s );
      west = x[j][i-1].tempc;
      lambda_W = LiquCond( west );
      n = x[j+1][i].tempc;
      lambda_N = LiquCond( n );
      rightside = alpha * AA * 300.0;
      b         = 1.0 /( (hx*hx) / (2*epsilon*lambda_C) + (hx*hx) / (2*epsilon*lambda_S) ) + ( epsilon * rho_S * Cp_S * x[j-1][i].u ) / hx;
      c         = 1.0 /( (hx*hx) / (2*epsilon*lambda_C) + (hx*hx) / (2*epsilon*lambda_N) );
      d         = 1.0 /( (hy*hy) / (2*epsilon*lambda_C) + (hy*hy) / (2*epsilon*lambda_W) );
      e         = 0.0;
      a         = 1.0 /( (hx*hx) / (2*epsilon*lambda_C) + (hx*hx) / (2*epsilon*lambda_S) ) + c + d + ( epsilon * rho_C * Cp_C * x[j][i].u ) / hx + alpha * AA;
      f[j][i].tempc = a * x[j][i].tempc - b * x[j-1][i].tempc - c * x[j+1][i].tempc- d * x[j][i-1].tempc - rightside;

      /*velocity field*/
      center = x[j][i].tempc;
      rho_C = LiquDen( center );
      rightside = epsilon * rho_C * g -epsilon * ( x[j+1][i].p - x[j][i].p ) / hx;
      b         = 0.0;
      c         = 0.0;
      d         = 0.0;
      e         = 0.0;
      a         = w * rho_C;
      f[j][i].u = a * x[j][i].u - rightside;

      /*pressure field*/
      center = x[j][i].tempc;
      rho_C = LiquDen( center );
      s = x[j-1][i].tempc;
      rho_S = LiquDen( s );
      rightside = ((rho_C * g / w)-(rho_S * g / w )) / hx;
      b         = 1.0 / ( w * hx * hx );
      c         = 1.0 / ( w * hx * hx );
      d         = 0.0;
      e         = 0.0;
      a         = b + c;
      f[j][i].p = a*x[j][i].p - b * x[j-1][i].p - c * x[j+1][i].p -rightside;
      }
    }
  }

  /* Compute over the interior points */
  for (j=yints; j<yinte; j++) {
    for (i=xints; i<xinte; i++) {
        /*fluid temperature field*/
      center = x[j][i].tempc;
      lambda_C = LiquCond( center );
      rho_C = LiquDen( center );
      Cp_C = LiquHetM( center );
      s = x[j-1][i].tempc;
      lambda_S = LiquCond( s );
      rho_S = LiquDen( s );
      Cp_S = LiquHetM( s );
      west = x[j][i-1].tempc;
      lambda_W = LiquCond( west );
      n = x[j+1][i].tempc;
      lambda_N = LiquCond( n );
      east = x[j][i+1].tempc;
      lambda_E = LiquCond( east );
      rightside = alpha * AA * 350.0;
      b         = 1.0 /( (hx*hx) / (2*epsilon*lambda_C) + (hx*hx) / (2*epsilon*lambda_S) ) + ( epsilon * rho_S * Cp_S * x[j-1][i].u ) / hx;
      c         = 1.0 /( (hx*hx) / (2*epsilon*lambda_C) + (hx*hx) / (2*epsilon*lambda_N) );
      d         = 1.0 /( (hy*hy) / (2*epsilon*lambda_C) + (hy*hy) / (2*epsilon*lambda_W) );
      e         = 1.0 /( (hy*hy) / (2*epsilon*lambda_C) + (hy*hy) / (2*epsilon*lambda_E) );
      a         = 1.0 /( (hx*hx) / (2*epsilon*lambda_C) + (hx*hx) / (2*epsilon*lambda_S) ) + c + d + e + ( epsilon * rho_C * Cp_C * x[j][i].u ) / hx + alpha * AA;
      f[j][i].tempc = a * x[j][i].tempc - b * x[j-1][i].tempc - c * x[j+1][i].tempc- d * x[j][i-1].tempc -e * x[j][i+1].tempc - rightside;
      /*velocity field*/
      center = x[j][i].tempc;
      rho_C = LiquDen( center );
      rightside = epsilon * rho_C * g -epsilon * ( x[j+1][i].p - x[j][i].p ) / hx;
      b         = 0.0;
      c         = 0.0;
      d         = 0.0;
      e         = 0.0;
      a         = w * rho_C;
      f[j][i].u = a * x[j][i].u - rightside;
      /*pressure field*/
      center = x[j][i].tempc;
      rho_C = LiquDen( center );
      s = x[j-1][i].tempc;
      rho_S = LiquDen( s );
      rightside = ((rho_C * g / w)-(rho_S * g / w )) / hx;
      b         = 1.0 / ( w * hx * hx );
      c         = 1.0 / ( w * hx * hx );
      d         = 0.0;
      e         = 0.0;
      a         = b + c;
      f[j][i].p = a*x[j][i].p - b * x[j-1][i].p - c * x[j+1][i].p -rightside;
    }
  }

  /*
     Flop count (multiply-adds are counted as 2 operations)
  */
  ierr = PetscLogFlops(84.0*info->ym*info->xm);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode NonlinearGS(SNES snes, Vec X, Vec B, void *ctx)
{
  DMDALocalInfo  info;
  Field          **x,**b;
  PetscErrorCode ierr;
  Vec            localX, localB;
  DM             da;
  PetscInt       xints,xinte,yints,yinte,i,j,k,l;
  PetscInt       max_its,tot_its;
  PetscInt       sweeps;
  PetscReal      rtol,atol,stol;
  PetscReal      hx,hy,dhx,dhy,hxdhy,hydhx;
  PetscReal      grashof,prandtl,lid;
  PetscScalar    u,uxx,uyy,vx,vy,avx,avy,vxp,vxm,vyp,vym;
  PetscScalar    fu, fv, fomega, ftemp;
  PetscScalar    dfudu;
  PetscScalar    dfvdv;
  PetscScalar    dfodu, dfodv, dfodo;
  PetscScalar    dftdu, dftdv, dftdt;
  PetscScalar    yu=0, yv=0, yo=0, yt=0;
  PetscScalar    bjiu, bjiv, bjiomega, bjitemp;
  PetscBool      ptconverged;
  PetscReal      pfnorm,pfnorm0,pynorm,pxnorm;
  AppCtx         *user = (AppCtx*)ctx;

  PetscFunctionBeginUser;

  PetscFunctionReturn(0);
}

/*TEST

   test:
      nsize: 1
      args: $PETSC_DIR/lib/petsc/bin/petscmpiexec -n 1 ./ex213 -snes_fd -pc_type lu -snes_view -snes_converged_reason -snes_monitor -ksp_converged_reason -ksp_monitor_true_residual



TEST*/