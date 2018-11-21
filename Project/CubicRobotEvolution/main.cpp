#include <iostream>
#include <fstream>
#include <GL/glut.h>
#include <GL/freeglut.h>
#include <sys/stat.h>
#include <string>
#include <mpi.h>
#include "test.h"
#include "OpenGL.h"
#include "PhysicsEngine.h"
#include "CubeGenerator.h"
#include "utility.h"
#include "GeneticOperators.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;

// Time
int frame = 0;
double dt_ms = 1;

//----------------------------------------------------
// Function prototype
//----------------------------------------------------

void timerDisplay(int value);
void timerPhysicsEngine(int value);
void DisplayCube(int, char**);

//----------------------------------------------------
// Main
//----------------------------------------------------
int main(int argc, char **argv){
  int rank, size;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size( MPI_COMM_WORLD, &size );

  //Create a folder for output files.
  mkdir("png", 0775);
  mkdir("log", 0775);

  GenerateCube();  //initialize cube robot setting.


  /* Evolve from scratch & display */
  if (*argv[1]=='e'){
    EvolveCube();
    MPI_Finalize();
    if (rank==0){
      cout << "Press enter to see the best motion." << endl;
      getchar();
      DisplayCube(argc, argv);
    }
  }

  /* Restart & display */
  else if (*argv[1]=='r'){
    Restart();
    MPI_Finalize();
    if (rank==0){
      cout << "Press enter to see the best motion." << endl;
      getchar();
      DisplayCube(argc, argv);
    }
  }

  /* Display */
  else if (*argv[1]=='d'){
    MPI_Finalize();
    if (rank==0){
      cout << "Display mode" << endl;
      DisplayCube(argc, argv);
    }
  }

  /* For test */
  else if (*argv[1]=='t'){
    test();
    MPI_Finalize();
    exit(0);
  }

  else{
    MPI_Finalize();
    cout << "Invalid arguments" << endl;
    exit(0);
  }

  return 0;
}

void DisplayCube(int argc, char **argv){

  InitializeCube();
  setBestIndividual();

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
  Initialize();      //initialize display setting.

  // Main loop start.
  glutMainLoop();
}

// ---------------------------------------------------
// Timer
// ---------------------------------------------------
void timerPhysicsEngine(int value){

  PhysicsEngine();

  // Call timer function after 100ms.
  glutTimerFunc(dt_ms, timerPhysicsEngine, 0);

  if (nt%60 == 0){
    glutPostRedisplay();
    TrajectoryRecord();
    EnergyRecord();
  }

  t = t + dt;
  nt += 1;

  //if (t>0.5) _Damping = true;
  //if (t>10.0) exit(0);
}
