#include "Obstacle.h"
#include "defines.h"


// Obstacle constructor
Obstacle::Obstacle(double mx, double my)
{
    m_posx = mx;
    m_posy = my;
    m_velx = 0;
    m_vely = 0;
    m_mass = 3; // initialize with any mass
    m_radius = 15;// initialize with any radius [px]
}


// Obstacle position update according to momentum conservation
void Obstacle::updatePos(const float* u, const float*v, const float*dens, const float dt)
{

    int intX = m_posx; // make integers of position (for better handling)
    int intY = m_posy;


    // MOMENTUM BALANCE - CALCULATE RESULTING FORCE FROM OBSTACLE "SURFACE-INTEGRAL"

    double Fx = 0; // initialize force variables
    double Fy = 0;

    int tempX,tempY;
    int nRound = 10; // use 10 points to extrapolate a surface integral of resulting force on object

    for (int i  = 0; i < nRound; i++)
    {
        tempX = intX + (int)(m_radius+2)*sin(2*PI*i/nRound);
        tempY = intY + (int)(m_radius+2)*cos(2*PI*i/nRound);

        if (tempX > 0 && tempX < N && tempY > 0 && tempY < N){ // make sure we use no index out of simulation grid
            Fx += dens[IX(tempX,tempY)]*(u[IX(tempX,tempY)]);
            Fy += dens[IX(tempX,tempY)]*(v[IX(tempX,tempY)]);
        }
    }


    Fy -= 9.81*m_mass*50; // influence of gravity, I use some scaling to get a better behavior


    // update position and velocities according to external forces

    m_velx += Fx/m_mass *10000*dt; // time scaling for realistic behavior of the object
    m_vely += Fy/m_mass *10000*dt;

    if(m_posx < m_radius || m_posx > (N-m_radius)) m_velx = -m_velx; // repel from solid walls, no damping
    if(m_posy < m_radius || m_posy > (N-m_radius)) m_vely = -m_vely;

    m_posx += m_velx*dt; // use explicit Euler method for time integration
    m_posy += m_vely*dt;

}




void Obstacle::setX(float x)
{
    m_posx = x;
}

void Obstacle::setY(float y)
{
    m_posy = y;
}

void Obstacle::testMember()
{
    std::cout << "testMember: m_posx:" << m_posx<< ", m_posy:" << m_posy << " velx: " << m_velx << " vely: " << m_vely << std::endl;
}