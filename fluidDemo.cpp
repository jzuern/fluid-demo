#include "defines.h"
#include "fluidFuncs.cpp"
#include "Obstacle.cpp"

using namespace std;

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%     EXPLANATION     %%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// s : Source for density
// i -> Index right/left
// j -> Index up/down
// u -> velocity right/left
// v -> velocity up/down

int main() // MAIN FUNCTION
{

    //densPlot << "set term wxt\n";
    //clock setup
    double duration;
    double start = omp_get_wtime();
    double tick = omp_get_wtick();

    // Initialization
    initializeGrid();
    Obstacle particle = defineObstacles();
    initializeFluid(u,v,dens);

densPlot << "set term x11 0\n";
densPlot << "plot sin(x)\n";
cout << "hit enter to proceed!\n";
std::cin.ignore(); // wait for user to hit enter
    densPlot << "set term xterm\n";
    //densPlot << "set term wxt\n";

    // Console output of parameters
    cout<< "time step dt:"<< dt << endl;
    cout<< "grid size:"<< N << " x " << N << " cells." << endl;
    cout<< "::::::::::::::::::::::::::::::::::::::" << "\n" << "\n" << endl;

    while (t<tmax ) // simulation loop
    {
        cout<< "Iteration " << counter << ", Zeit t = "<< t << endl;

        get_from_UI( dens,u,v,dens_prev, u_prev, v_prev ,t,densPlot); // external influence on fluid

        // NAVIER-STOKES SOLUTION: VELOCITY FIELD AND DENSITY FIELD SEPARATELY SOLVED
        vel_step( N, u, v, u_prev, v_prev, visc, dt);
        dens_step( N, dens, dens_prev, u, v, diff, dt);


        // MOVE OBSTACLE IN FLUID AROUND
        particle = obstacle_step(particle);
        updateGrid(particle);

        // VISUALIZATION
//                drawVel();
//                drawAbsVel();
        drawDens();



        // Measurement of performance
        //        if (counter % 10 == 0) {
        //            double end = omp_get_wtime();
        //            duration = (end-start);
        //            cout << "Mean CPU time per Iteration: " << duration/10 << " sec/iteration\n";
        //            cout << tick << endl;
        //            start = omp_get_wtime();
        //        }



        t += dt;
        ++counter;

    } // end simulation loop

    return 0;
}



Obstacle defineObstacles()
{


    double mx=0.5, my=0.5; // initialize mouse button variables
    int mb=1;

    densPlot << "set title 'Please select obstacle positions\n";
    densPlot << "set cbrange [0:1]\n";
    densPlot << "set yrange [0:" << N << "]\n";
    densPlot << "set xrange [0:" << N << "]\n";
    densPlot << "set nokey\n";
    densPlot << "set pal gray\n";
    densPlot << "plot -1 with lines\n";
    densPlot << "set size square\n";

    densPlot.getMouse(mx, my, mb, "Left click to select area for obstacle, right click to exit.");

    Obstacle particle(mx,my); // initialize Obstacle with coordinates of mouse pointer click


    // update visualization
    vector<vector<double> > grid;
    for(int j=0; j<N; j++) {
        vector<double> row;
        for(int i=0; i<N; i++) {
            bool occ = (i-mx)*(i-mx) + (j-my)*(j-my) < particle.getRad()*particle.getRad(); // determine occupation
            occupiedGrid[IX(i,j)] = occupiedGrid[IX(i,j)] || occ; // update occupied grid
            row.push_back((double)occupiedGrid[IX(i,j)]);
        }
        grid.push_back(row); // build image of grid
    }
    densPlot << "plot '-' binary" << densPlot.binFmt2d(grid, "array") << "with image\n";
    densPlot.sendBinary2d(grid);


    if(mb < 0) {
        printf("The gnuplot window was closed.\n");
    }


#pragma omp parallel for
        for ( int i=1 ; i<=N ; i++ ) {
            for ( int j=1 ; j<=N ; j++ ) {
                bool occ = (i-mx)*(i-mx) + (j-my)*(j-my) < particle.getRad()*particle.getRad();
                occupiedGrid[IX(i,j)] = occupiedGrid[IX(i,j)] || occ;
            }
        }

    return particle;
}

void drawVel()
{
    velPlot << "set xrange [0:" << N << "]\n";
    velPlot << "set yrange [0:" << N << "]\n";

    vector<boost::tuple<double, double, double, double> > vecs;
    for(int i=0; i<maxDraw; i++) {
        for(int j=0; j<maxDraw; j++){
            float norm = sqrt(u[IX(drawRatio*i,drawRatio*j)]*u[IX(drawRatio*i,drawRatio*j)]+v[IX(drawRatio*i,drawRatio*j)]*v[IX(drawRatio*i,drawRatio*j)]);
            int ir = drawRatio*i;
            int jr = drawRatio*j;
            vecs.push_back(boost::make_tuple(
                               (double)ir,      (double)jr,
                               (double)u[IX(ir,jr)]/norm*drawRatio, (double)v[IX(ir,jr)]/norm*drawRatio
                    ));
        }
    }

    velPlot << "plot '-' with vectors\n";
    velPlot.send1d(vecs);
}


void drawAbsVel()
{

    absVelPlot << "set cbrange [0:500]\n";
    absVelPlot << "set xrange [0:" << N << "]\n";
    absVelPlot << "set yrange [0:" << N << "]\n";
    absVelPlot << "set title 'Fluid Speed'\n";
    absVelPlot << "set pal defined\n";
    absVelPlot << "set size square\n";

    vector<vector<double> > image;
    for(int j=0; j<N; j++) {
        vector<double> row;
        for(int i=0; i<N; i++) {
            float abs = sqrt(u[IX(i,j)]*u[IX(i,j)]+v[IX(i,j)]*v[IX(i,j)]);
            double z = (double)abs;
            row.push_back(z);
        }
        image.push_back(row);
    }

    absVelPlot << "plot '-' binary" << absVelPlot.binFmt2d(image, "array") << "with image\n";
    absVelPlot.sendBinary2d(image);
}


void drawDens()
{
	//densPlot << "set term wxt\n";
    densPlot << "set cbrange [0:15]\n";
    densPlot << "set xrange [0:" << N << "]\n";
    densPlot << "set yrange [0:" << N << "]\n";
    densPlot << "set title 'Density Field'\n";
    densPlot << "set pal defined\n";
    densPlot << "set size square\n";


    // convert density array to double vector in order to use Gnuplot.binFmt2d() function for fast image plotting
    vector<vector<double> > image;
    for(int j=0; j<N; j++) {
        vector<double> row;
        for(int i=0; i<N; i++) {
            double z = (double)dens[IX(i,j)];
            row.push_back(z);
        }
        image.push_back(row);
    }

    densPlot << "plot '-' binary" << densPlot.binFmt2d(image, "array") << "with image\n";
    densPlot.sendBinary2d(image);
}



Obstacle obstacle_step(Obstacle &particle)
{
    particle.updatePos(u,v,dens,dt);
    return particle;
}


void updateGrid(Obstacle &particle)
{

    float posX = particle.getPosX();
    float posY= particle.getPosY();
    float rad = particle.getRad();

    for(int j=0; j<N; j++) {
        for(int i=0; i<N; i++) {
            occupiedGrid[IX(i,j)] = ((i-posX)*(i-posX) + (j-posY)*(j-posY) < rad*rad) ? 1:0;
        }
    }

}
