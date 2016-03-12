#ifndef OBSTACLE_H
#define OBSTACLE_H


class Obstacle
{
private:
    double m_posx; // x Position of Obstacle
    double m_posy; // y Position of Obstacle
    double m_vely; // x Velocity of Obstacle
    double m_velx; // y Velocity of Obstacle

    double m_mass; // mass of Obstacle
    double m_radius;  // Obstacle radius

    Obstacle() { } // private default constructor

public:
    Obstacle(double mx, double my); // constructor

    void updatePos(const float *u, const float *v, const float *dens, const float dt); // Obstacle position update

    float getPosX() { return m_posx; } // getter and setter functions for private members
    float getPosY()  { return m_posy; }
    float getRad()  { return m_radius; }
    void setX(float x);
    void setY(float y);

    void testMember(); // testfunction for debugging
};

#endif
