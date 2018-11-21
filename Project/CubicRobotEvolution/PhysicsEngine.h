#ifndef _PHYSICSENGINE_H_
#define _PHYSICSENGINE_H_

#include <vector>
#include <chrono>
using std::vector;


void PhysicsEngine(void);
void EnergyRecord(void);
void TrajectoryRecord(void);

extern int nt;
extern double t;
extern double dt;

#endif
