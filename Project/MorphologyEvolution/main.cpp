#include <iostream>
#include <fstream>
#include <GL/glut.h>
#include <GL/freeglut.h>
#include <sys/stat.h>
#include <string>
#include <mpi.h>
#include <algorithm>

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
vector<Mass> robots[10];
vector<vector<double>> df_robots[10];
int N_ROBOTS;
double END_TIME = 100.0;
bool save=false;
int slot;

//----------------------------------------------------
// Function prototype
//----------------------------------------------------

void timerPhysicsEngine(int value);
void timerMultiRobots( int value);
void DisplayCube(int, char**);
void readParameterFile(string, int);
void readMotionFile(char*);

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
  mkdir("log/parameter", 0775);
  mkdir("log/trajectory", 0775);
  mkdir("log/fitness", 0775);
  mkdir("log/robots", 0775);

  //GenerateTetra();  //initialize cube robot setting.
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
      cout << "Save? (T/F)" << endl;
      char char_save;
      std::cin >> char_save;
      if (char_save=='T'){
        save = true;
        cout << "Slot? (0-9)" << endl;
        std::cin >> slot;
      }
      DisplayCube(argc, argv);
    }
  }

  /* Multi-robot display */
  else if (*argv[1]=='m'){
    MPI_Finalize();
    if (rank==0){
      cout << "Multi-robots display mode" << endl;
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

void DisplayCube(int argc, char** argv){

  InitializeCube();

  if (*argv[1]=='e' || *argv[1]=='r'){
    readParameterFile("./log/parameter/para_" + std::to_string(gen-1) + ".csv", 0);
    robots[0] = mass;
  }else if (*argv[1]=='d'){
    if(argc==3){
      int id;
      cout << "Robot ID?"<< endl;
      std::cin >> id;
      readParameterFile(argv[2], id);
    }
    robots[0] = mass;
  }else if(*argv[1]=='m'){
    readMotionFile(argv[2]);
    for(int i=0; i<N_ROBOTS; i++) robots[i] = mass;
  }

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

  if(*argv[1]!='m'){
    glutTimerFunc(dt_ms, timerPhysicsEngine, 0); // Timer
  }else{
    glutTimerFunc(dt_ms, timerMultiRobots, 1); // Timer
  }

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
  robots[0] = mass;

  // Call timer function after 100ms.
  glutTimerFunc(dt_ms, timerPhysicsEngine, 0);

  if (nt%60 == 0){
    glutPostRedisplay();
    if (save){
      TrajectoryRecord(slot);
      EnergyRecord();
    }
  }

  if (!_Stop){
    t = t + dt;
    nt += 1;
  }

  if (t>END_TIME) exit(0);
}

// ---------------------------------------------------
// Timer
// ---------------------------------------------------
void timerMultiRobots(int value){

  int INI = nt*N_MASS;
  int FIN = INI + N_MASS;

  for (int i=0; i<N_ROBOTS; i++){
    vector<double> offset = {0.0, (double)i*0.5, 0.0};
    for (int j=INI; j<FIN; j++){
      robots[i][j-INI].p = {df_robots[i][j][2] + offset[0],
                            df_robots[i][j][3] + offset[1],
                            df_robots[i][j][4] + offset[2]};
    }
  }

  t = df_robots[0][INI][0];


  // Call timer function after 100ms.
  glutTimerFunc(30, timerMultiRobots, 0);
  glutPostRedisplay();

  if (!_Stop) nt += 1;
  if (t >= END_TIME) exit(0);
}

void readParameterFile(string filename, int id){

  Individual ind;
  //int id = 0;

  cout << "Inputfile: " << filename << endl;
  vector<vector<string>> strvec = read_csv_string(filename, 0, 0);
  const char* rep = strvec[0][0].c_str();
  setRepresentation(rep[0]);

  int INI = 2 + SIZE_OF_CHROMOSOME * id;
  int FIN = INI + SIZE_OF_CHROMOSOME - 1;

  vector<vector<double>> parameter_table = read_csv(filename, INI, FIN);
  vector<int> num_of_parts = Decorder(parameter_table, strvec[0][0]);

  cout << "Num. of Masses : " << num_of_parts[0] << endl;
  cout << "Num. of Springs: " << num_of_parts[1] << endl;
}

void readMotionFile(char* folder){
  vector<double> time_list;
  vector<string> filenames;
  vector<string> strvec = get_filename(folder);

  for (string str: strvec) {
    string keyword = "robot";
    if (std::equal(keyword.begin(), keyword.end(), str.begin())){
      filenames.push_back(string(folder) + str);
    }
  }

  N_ROBOTS = filenames.size();
  cout << "Num of robots: " << N_ROBOTS << endl;
  for (int i=0; i<N_ROBOTS; i++){
    df_robots[i] = read_csv(filenames[i], 1, -1);
    double end_time = df_robots[i].back()[0];
    time_list.push_back(end_time);
    cout << filenames[i] << " | END TIME: "<< end_time << endl;
  }

  END_TIME = *std::min_element(time_list.begin(), time_list.end());

}
