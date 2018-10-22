#include <iostream>
#include <GL/glut.h>
#include <GL/freeglut.h>
#include <sys/stat.h>
#include "opengl.h"
#include "cube_model.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;

// Time
int frame = 0 ;
double t = 0.0;
int nt = 0;

double dt = 0.001;
double dt_ms = 1;

//----------------------------------------------------
// Function prototype
//----------------------------------------------------

void timerDisplay(int value);
void timerPhysicsEngine(int value);

//----------------------------------------------------
// Main
//----------------------------------------------------
int main(int argc, char *argv[]){

  //Create a folder for output files.
  mkdir("png", 0775);
  mkdir("log", 0775);

  // GLUT initial setting
  glutInit(&argc, argv);                                    //Initialize environment.
  glutInitWindowPosition(WindowPositionX, WindowPositionY); //Define window position.
  glutInitWindowSize(WindowWidth, WindowHeight);            //Define window size.
  glutCreateWindow(WindowTitle);                            //Create window.
  glutInitDisplayMode(GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE);//Define display mode.

  // Assign functions.
  glutDisplayFunc(Display);     // Display function.
  glutKeyboardFunc(Keyboard);   // Keyboard function.
  glutSpecialFunc(SpecialKey);  // Special Keyboard Function.
  glutIdleFunc(Idle);           // Idle function.
  glutMouseFunc(mouse_on);      // Mouse on function
  glutMotionFunc(mouse_motion); // Mouse drag function
  glutTimerFunc(dt_ms, timerPhysicsEngine, 0); // Timer

  // Initialization
  InitializeCube();  //initialize cube robot setting.
  Initialize();      //initialize display setting.

  // Main loop start.
  glutMainLoop();

  return 0;
}

// ---------------------------------------------------
// Timer
// ---------------------------------------------------
void timerPhysicsEngine(int value){

  PhysicsEngine();

  // Call timer function after 100ms.
  glutTimerFunc(dt_ms, timerPhysicsEngine, 0);

  if (nt%10 == 0){
    glutPostRedisplay();
  }

  t = t + dt;
  nt += 1;

  if (t>100.0) exit(0);
}
