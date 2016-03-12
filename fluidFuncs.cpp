
void initializeGrid()
// initialize the grid cell coordinates
{
    for ( int i=1 ; i<=N ; i++ ) {
        Nx[i] = i;
        Ny[i] = i;
    }
}


void initializeFluid(float *u,float *v,float *dens)
{
    // initialize the fluid state (velocities and density) at t=0
#pragma omp parallel for
    for ( int i=1 ; i<=N ; i++ ) {
        for ( int j=1 ; j<=N ; j++ ) {
            u[IX(i,j)] = 0.01; // u velocity at t=0
            v[IX(i,j)] = 0.01; // v velocity at t=0
            dens[IX(i,j)] = 1.; // density at t=0
        }
    }
}

void diffuse( int N, int b, float * x, float * x0, float diff, float dt )
{
    // diffusion step is obtained by Gauss-Seidel relaxation equation system solver
    // used for density, u-component and v-component of velocity field separately

    float a=dt*diff*N*N;

#pragma omp parallel for
    for (int k=0 ; k<20 ; k++ ) {
        for (int i=1 ; i<=N ; i++ ) {
            for (int j=1 ; j<=N ; j++ ) {
                x[IX(i,j)] = (x0[IX(i,j)] + a*(x[IX(i-1,j)]+x[IX(i+1,j)]+x[IX(i,j-1)]+x[IX(i,j+1)]))/(1+4*a);
            }
        }
        set_bnd( N, b, x );
    }
}


void advect( int N, int b, float * d, float * d0, float * u, float * v, float dt )
{
    // calculate the advection of density in velocity field and velocity field along itself

    int i, j, i0, j0, i1, j1;
    float x, y, s0, t0, s1, t1, dt0;
    bool occ;
    dt0 = dt*N;

    for ( i=1 ; i<=N ; i++ ) {
        for ( j=1 ; j<=N ; j++ ) {

            x = i-dt0*u[IX(i,j)];
            y = j-dt0*v[IX(i,j)];

            if (x<0.5) x=0.5;
            if (x>N+0.5) x=N+ 0.5; i0=(int)x; i1=i0+ 1;
            if (y<0.5) y=0.5;
            if (y>N+0.5) y=N+ 0.5; j0=(int)y; j1=j0+1;
            s1 = x-i0;
            s0 = 1-s1;
            t1 = y-j0;
            t0 = 1-t1;

            occ = occupiedGrid[IX(i,j)];

            if(occ == 0) d[IX(i,j)] = s0*(t0*d0[IX(i0,j0)]+t1*d0[IX(i0,j1)])+s1*(t0*d0[IX(i1,j0)]+t1*d0[IX(i1,j1)]);
            else
                d[IX(i,j)] = 0;
        }
    }

    set_bnd( N, b, d );
}



void add_source ( int N, float * x, float * s, float dt )
{
    // add sources for velocity field or density field
    int i, size=(N+2)*(N+2);
    for ( i=0 ; i<size ; i++ ){
        x[i] += dt*s[i];
    }
}

void dens_step (int N, float * x, float * x0, float * u, float * v, float diff,float dt)
{
    // executes all routines for motion of density field in one time step
    add_source( N, x, x0, dt );
    SWAP ( x0,x);
    diffuse( N, 0, x, x0, diff, dt );
    SWAP ( x0,x);
    advect( N, 0, x, x0, u, v, dt );
}


void vel_step ( int N, float * u, float * v, float *  u0, float * v0,float visc, float dt )
{
    // executes all routines for motion of velocity field in one time step
    add_source ( N, u, u0, dt );
    SWAP ( u0, u );
    diffuse(N, 1, u, u0, visc, dt);

    add_source ( N, v, v0, dt );
    SWAP ( v0, v );
    diffuse(N, 2, v, v0, visc, dt);

    project ( N, u, v, u0, v0 );

    SWAP ( u0, u );
    SWAP ( v0, v );

    advect(N, 1, u, u0, u0, v0, dt );
    advect(N, 2, v, v0, u0, v0, dt );

    project ( N, u, v, u0, v0 );
}

void project ( int N, float * u, float * v, float * p, float * div )
{
    // force routing to be mass conserving (use "hodge decomposition" for obtained velocity field and
    // eliminate gradient field)
    // this will make the velocity field to have fluid-like swirls as desired

    float h;
    h = 1.0/N;

#pragma omp parallel for
    for (int i=1 ; i<=N ; i++ ) {
        for (int j=1 ; j<=N ; j++ ) {
            div[IX(i,j)] = -0.5*h*(u[IX(i+1,j)]-u[IX(i-1,j)]+v[IX(i,j+1)]-v[IX(i,j-1)]);
            p[IX(i,j)] = 0;
        }
    }

    set_bnd( N, 0, div );
    set_bnd( N, 0, p );

#pragma omp parallel for
    for (int k=0 ; k<20 ; k++ ) {
        for (int i=1 ; i<=N ; i++ ) {
            for (int j=1 ; j<=N ; j++ ) {
                p[IX(i,j)] = (div[IX( i,j)]+p[IX(i-1,j)]+p[IX(i+1,j)]+p[IX(i,j-1)]+p[IX(i,j+1)])/4;
            }
        }
        set_bnd( N, 0, p );
    }

#pragma omp parallel for
    for (int i=1 ; i<=N ; i++ ) {
        for (int j=1 ; j<=N ; j++ ) {
            u[IX(i,j)] -= 0.5*(p[IX(i+1,j)]-p[IX(i-1,j)])/h;
            v[IX(i,j)] -= 0.5*(p[IX(i,j+1)]-p[IX(i,j-1)])/h;
        }
    }
    set_bnd( N, 1, u );
    set_bnd( N, 2, v );
}








void set_bnd(int N, int b, float * x)
{
    // define boundary values for velocity and density

    for (int i=1 ; i<=N; i++ ) {

        // left border

        x[IX(0 ,i)] = b==1 ? -x[IX(1,i)] : x[IX(1,i)]; // bounded box boundary conditions



        // right border
        x[IX(N+1,i)] = b==1 ? -x[IX(N,i)] : x[IX(N,i)]; // bounded box boundary conditions


        // bottom border
        //        x[IX(i,0 )] = b==2 ? -x[IX(i,1)] : x[IX(i,1)]; // bounded box boundary conditions

        // make jets at the bottom
        x[IX(i,0 )] = 0;

        if ((i > 4.5*N/10 && i < 5.5*N/10))
            x[IX(i,0 )] += 0.4;

        else if ((i > 0.5*N/10 && i < 1.8*N/10))
            x[IX(i,0 )] += 0.4;
        else if ((i > 7.8*N/10 && i < 8.9*N/10))
            x[IX(i,0 )] += 0.4;


        // upper border
        x[IX(i,N+1)] = b==2 ? -x[IX(i,N)] : x[IX(i,N)]; // bounded box boundary conditions
    }



    // implementing internal flow obstacles
    if(b != 0) {  // only changed boundaries for flow -> b = 1,2

#pragma omp parallel for
        for ( int i=1 ; i<=N ; i++ ) {
            for ( int j=1 ; j<=N ; j++ ) {
                bool occ = occupiedGrid[IX(i,j)];
                if(occ == 1){
                    x[IX(i-1,j)] = b==1 ? -x[IX(i,j)] : x[IX(i,j)];
                    x[IX(i+1,j)] = b==1 ? -x[IX(i,j)] : x[IX(i,j)];
                    x[IX(i,j-1 )] = b==2 ? -x[IX(i,j)] : x[IX(i,j)];
                    x[IX(i,j+1)] = b==2 ? -x[IX(i,j)] : x[IX(i,j)];
                }
            }
        }
    }

    // define edge cells as median of neighborhood
    x[IX(0 ,0 )] = 0.5*(x[IX(1,0 )]+x[IX(0 ,1)]);
    x[IX(0 ,N+1)] = 0.5*(x[IX(1,N+1)]+x[IX(0 ,N )]);
    x[IX(N+1,0 )] = 0.5*(x[IX(N,0 )]+x[IX(N+1,1)]);
    x[IX(N+1,N+1)] = 0.5*(x[IX(N,N+1)]+x[IX(N+1,N )]);
}



void get_from_UI(float *dens, float *u, float *v, float *dens_prev, float *u_prev, float *v_prev, float t, Gnuplot &densPlot)
{
    // adds density and/or velocity in this function

#pragma omp parallel for
    for ( int i=1 ; i<=N ; i++ ) {
        for ( int j=1 ; j<=N ; j++ ) {


            // Density input for nice jets
            if ((i > 4.6*N/10 && i < 5.4*N/10) && (j<3))
                dens[IX(i,j)] = 10.0;
            else if ((i > 0.6*N/10 && i < 1.7*N/10)&& (j<3))
                dens[IX(i,j)] = 10.0;
            else if ((i > 7.9*N/10 && i < 8.8*N/10)&& (j<3))
                dens[IX(i,j)] = 10.0;

        }
    }
}





