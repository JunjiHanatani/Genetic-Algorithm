#include <math.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string.h>
#include <vector>
#include <GL/glut.h>
#include <sys/stat.h>
#include <png.h>
using std::cout;
using std::endl;
using std::string;
using std::vector;

#define SCALE 500.0
#define N_MASS 8
#define N_SPRING 28
double PI = acos(-1.0);

// Window Setting
int WindowPositionX = 200;  //Window x position
int WindowPositionY = 200;  //Window y position
int WindowWidth = 512;    //Window width
int WindowHeight = 512;    //Window height
char WindowTitle[] = "The beginning of the world";  //Window title

// Output video
bool _Png = false;

// Time
int tn = 0 ;
double t = 0.0;
double dt = 0.001;
double dt_ms = 1;

// Masses
struct {
  double m = 0.1;
  vector<double> p = {0.0, 0.0, 0.0};
  vector<double> v = {0.0, 0.0, 0.0};
  vector<double> a = {0.0, 0.0, -9.81};
}mass[N_MASS];
double MassR = 0.01;

// Springs
struct {
  double k = 10000.0;
  double l0 = 0.1;
  vector<int> masses={0, 0};
}spring[N_SPRING];

// Ball (test)
struct {
  double x, y, z;
  double vx, vy, vz;
}p[100];
int pn = 0;
double ax = 0.0 , ay = 0.0 , az = -9.81;
double vx = 0.05 , vy = 0.1 , vz = 1.5;
double xini = 0.0 , yini = 0.0 , zini = 0.1;
double BallR = 0.004;
double hanpatu = 0.95;

//----------------------------------------------------
// Initial Viewpoint
//----------------------------------------------------
double InitialViewPointX = 0.0;
double InitialViewPointY = -500.0;
double InitialViewPointZ = 100.0;
double ViewPointR = sqrt(InitialViewPointX * InitialViewPointX + InitialViewPointY * InitialViewPointY);
double ViewPointTheta = atan2(InitialViewPointY, InitialViewPointX);
double omega = 2.0 * PI * 0.1;

//----------------------------------------------------
// Rotation by mouse action
//----------------------------------------------------

int cx, cy;                // Drag start position
double sx = 1.0 / (double)512; // Transform coefficient from absolute position of mouse to relative position in the window
double sy = 1.0 / (double)512;
double cq[4] = { 1.0, 0.0, 0.0, 0.0 };  // Initial orientation (Quaternion)
double tq[4];              // Orientation during dragging (Quaternion)
double rt[16];              // Rotation transform matrix

unsigned int listNumber;
float camera_z_pos =50.0;

//----------------------------------------------------
// Define a box.
//----------------------------------------------------
GLdouble vertex[][3] = {//Define vertices.
  { 0.0, 0.0, 0.0 },
  { 2.0, 0.0, 0.0 },
  { 2.0, 2.0, 0.0 },
  { 0.0, 2.0, 0.0 },
  { 0.0, 0.0, 30.0 },
  { 2.0, 0.0, 30.0 },
  { 2.0, 2.0, 30.0 },
  { 0.0, 2.0, 30.0 }
};
int face[][4] = {//Define surfaces.
  { 0, 1, 2, 3 },
  { 1, 5, 6, 2 },
  { 5, 4, 7, 6 },
  { 4, 0, 3, 7 },
  { 4, 5, 1, 0 },
  { 3, 2, 6, 7 }
};
GLdouble normal[][3] = {//Normal vector of the surface.
  { 0.0, 0.0,-1.0 },
  { 1.0, 0.0, 0.0 },
  { 0.0, 0.0, 1.0 },
  {-1.0, 0.0, 0.0 },
  { 0.0,-1.0, 0.0 },
  { 0.0, 1.0, 0.0 }
};

//----------------------------------------------------
// Texture
//----------------------------------------------------
struct MaterialStruct {
  GLfloat ambient[4];
  GLfloat diffuse[4];
  GLfloat specular[4];
  GLfloat shininess;
};
//jade
MaterialStruct ms_jade = {
  {0.135,     0.2225,   0.1575,   1.0},
  {0.54,      0.89,     0.63,     1.0},
  {0.316228,  0.316228, 0.316228, 1.0},
  12.8};
//ruby
MaterialStruct ms_ruby  = {
  {0.1745,   0.01175,  0.01175,   1.0},
  {0.61424,  0.04136,  0.04136,   1.0},
  {0.727811, 0.626959, 0.626959,  1.0},
  76.8};

//----------------------------------------------------
// Color
//----------------------------------------------------
GLfloat red[] = { 0.8, 0.2, 0.2, 1.0 };
GLfloat green[] = { 0.2, 0.8, 0.2, 1.0 };
GLfloat blue[] = { 0.2, 0.2, 0.8, 1.0 };
GLfloat yellow[] = { 0.8, 0.8, 0.2, 1.0 };
GLfloat white[] = { 1.0, 1.0, 1.0, 1.0 };
GLfloat shininess = 30.0;

//----------------------------------------------------
// Function prototype
//----------------------------------------------------
void Initialize(void);
void Idle(void);
void Display(void);
void Ground(void);
void StaticObjects(void);
void PopUpBall(void);
void Keyboard(unsigned char key, int x, int y);
void DRAW_STRING(int x, int y, string str, void *font = GLUT_BITMAP_TIMES_ROMAN_24);
void DISPLAY_TEXT(int x, int y, string str);
void capture(const string fname);
void qmul(double r[], const double p[], const double q[]);
void qrot(double r[], double q[]);
void mouse_motion(int x, int y);
void mouse_on(int button, int state, int x, int y);
void mouse_wheel(float z);
void timer(int value);
double calc_distance(vector<double>p0, vector<double>p1);
void CubicRobot();
void InitializeCube();

//----------------------------------------------------
// Main
//----------------------------------------------------
int main(int argc, char *argv[]){
  if(_Png) mkdir("png", 0775); //Create a folder for png files.

  glutInit(&argc, argv);//Initialize environment.
  glutInitWindowPosition(WindowPositionX, WindowPositionY);//Define window position.
  glutInitWindowSize(WindowWidth, WindowHeight); //Define window size.
  glutCreateWindow(WindowTitle);  //Create window.
  glutInitDisplayMode(GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE);//Define display mode.

  glutDisplayFunc(Display);     // Assign a display function.(Function: Display)
  glutKeyboardFunc(Keyboard);   // Assign a keyboard function. (Function: Keyboard)
  glutIdleFunc(Idle);           // Assign an idle function. (Function: Idle)
  glutMouseFunc(mouse_on);      // Assign a mouse on function
  glutMotionFunc(mouse_motion); // Assign a mouse drag function
  glutTimerFunc(dt_ms, timer, 0);
  Initialize(); //Call initial setting.
  glutMainLoop();
  return 0;
}

//----------------------------------------------------
// Initialize function
//----------------------------------------------------
void Initialize(void){

  glClearColor(1.0, 1.0, 1.0, 1.0); //Assign background color. (R, G, B, alpha)
  glEnable(GL_DEPTH_TEST);//Use depth buffer: (Assign GLUT_DEPTH using glutInitDisplayMode())

  //Set a light source --------------------------------------
  GLfloat light_position0[] = { -50.0, -50.0, 20.0, 1.0 }; //Coordinate of a light source 0
  glLightfv(GL_LIGHT0, GL_POSITION, light_position0); //

  // Create display list
  listNumber = glGenLists(1);
  glNewList( listNumber, GL_COMPILE );
  glEndList();

  // Initialize rotate matrix
  qrot(rt, cq);

  // Initialize CubeRobot
  InitializeCube();

}

//----------------------------------------------------
// Display function
//----------------------------------------------------
void Display(void) {
  //Clear buffer.
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  // Set a perspective projection matrix ------------------------------
  glMatrixMode(GL_PROJECTION);//Matrix mode (GL PROJECTION: Perspective projection matrix / GL_MODELVIEW: model view transformation matrix)
  glLoadIdentity();//Initialize a matrix
  gluPerspective(30.0, (double)WindowWidth/(double)WindowHeight, 0.1, 1000.0); //Apparent volume gluPerspactive(th, w/h, near, far);

  // Set a model view matrix --------------------------
  glMatrixMode(GL_MODELVIEW);//
  glLoadIdentity(); //Initialize a matrix
  glViewport(0, 0, WindowWidth, WindowHeight);

  // Set a viewpoint------------------------------------------------
  double ViewPointX = ViewPointR * cos( omega * t + ViewPointTheta);
  double ViewPointY = ViewPointR * sin( omega * t + ViewPointTheta);
  double ViewPointZ = InitialViewPointZ;
  gluLookAt(
      ViewPointX, ViewPointY, ViewPointZ, // Viewpoint: x,y,z;
      0.0,        0.0,        0.0,        // Reference point: x,y,z
      0.0,        0.0,        1.0 );      //Vector: x,y,z

  // Rotation -------------------------
  glMultMatrixd(rt);

  //Shadow ON-----------------------------
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);//Use a light source 0

  // Objects -------------------------
  //PopUpBall();
  //StaticObjects();
  CubicRobot();

  // Shadow OFF -----------------------------
  glDisable(GL_LIGHTING);
  // -----------------------------------

  // Draw ground. --------------------
  Ground();
  glPopMatrix();
  // ---------------------------------

  // Draw characters. -------------------------------------------------------
  string t_str = "t = " + std::to_string(t);
  DISPLAY_TEXT(5, 95, t_str);
  string f_str = "frame = " + std::to_string(tn);
  DISPLAY_TEXT(5, 90, f_str);
  string p_str = "p = " + std::to_string(pn);
  DISPLAY_TEXT(5, 85, p_str);
  // ------------------------------------------------------------------------

  // Output png. -------------------------------------------------------
  if(_Png){
    int tt = tn +10000;
    string fname = "./png/" + std::to_string(tt) + ".png";//Output file name.
    capture(fname);
  }
  // ------------------------------------------------------------------------

  //glutInitDisplayMode(GLUT_DOUBLE) enables "double buffering"
  glutSwapBuffers();

  tn++ ;
}

//----------------------------------------------------
// Idling function
//----------------------------------------------------
void Idle(){
  glutPostRedisplay(); //Execute glutDisplayFunc() once.
}

//----------------------------------------------------
// Keyboard function.
//----------------------------------------------------
void Keyboard(unsigned char key, int x, int y){
  switch ( key )
  {
  case 'a':
    pn++;
    p[pn].x = xini;
    p[pn].y = yini;
    p[pn].z = zini;
    p[pn].vx = vx * ( (double)rand()/(double)RAND_MAX - (double)rand()/(double)RAND_MAX );
    p[pn].vy = vy * ( (double)rand()/(double)RAND_MAX );
    p[pn].vz = vz * ( (double)rand()/(double)RAND_MAX );
    break;
  case 'q':
    exit(0);
    break;

  default:
    break;
  }
}

void InitializeCube(void){

  mass[0].p = {0.0, 0.0, 0.0};
  mass[1].p = {0.1, 0.0, 0.0};
  mass[2].p = {0.1, 0.1, 0.0};
  mass[3].p = {0.0, 0.1, 0.0};
  mass[4].p = {0.0, 0.0, 0.1};
  mass[5].p = {0.1, 0.0, 0.1};
  mass[6].p = {0.1, 0.1, 0.1};
  mass[7].p = {0.0, 0.1, 0.1};

  for (int i=0; i<N_SPRING; i++){
    if(i<12){
        spring[i].l0 = 0.1;
    }else if(i<24){
        spring[i].l0 = 0.1 * sqrt(2);
    }else{
        spring[i].l0 = 0.1 * sqrt(3);
    }
  }

  int k=0;
  for (int i=0; i<N_MASS; i++){
    for (int j=i+1; j<N_MASS; j++){
      spring[k].masses = {i, j};
      spring[k].l0 = calc_distance(mass[i].p, mass[j].p);
      k += 1;
    }
  }

}

//----------------------------------------------------
// Ground function
//----------------------------------------------------

void Ground(void) {
    double ground_max_x = 300.0;
    double ground_max_y = 300.0;
    glColor3d(0.8, 0.8, 0.8);  // Color of the ground
    glLineWidth(1.0d);

    // Draw lines.
    glBegin(GL_LINES);
    // Horizontal lines.
    for(double ly = -ground_max_y ;ly <= ground_max_y; ly+=10.0){
      glVertex3d(-ground_max_x, ly,0);
      glVertex3d(ground_max_x, ly,0);
    }
    // Vertical lines.
    for(double lx = -ground_max_x ;lx <= ground_max_x; lx+=10.0){
      glVertex3d(lx, ground_max_y,0);
      glVertex3d(lx, -ground_max_y,0);
    }
    glEnd();
}

// ---------------------------------------------------
// Timer
// ---------------------------------------------------

void timer(int value){

  for(int i=1; i<=pn; i++){
    p[i].vx += ax * dt;
    p[i].vy += ay * dt;
    p[i].vz += az * dt;
    p[i].x += p[i].vx * dt;
    p[i].y += p[i].vy * dt;
    p[i].z += p[i].vz * dt;
    if(p[i].z < BallR){
      p[i].z = BallR;
      p[i].vz = -hanpatu * p[i].vz;
    }
  }

  // Update the display
  glutPostRedisplay();

  // Call timer function after 100ms.
  glutTimerFunc(dt_ms, timer, 0);

  t = t + dt;
}


//----------------------------------------------------
// Static Objects function
//----------------------------------------------------

void CubicRobot(void){
  // Masses
  for(int i=0; i<N_MASS; i++){
    glPushMatrix();
    glMaterialfv(GL_FRONT, GL_AMBIENT, ms_ruby.ambient);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, ms_ruby.diffuse);
    glMaterialfv(GL_FRONT, GL_SPECULAR, ms_ruby.specular);
    glMaterialfv(GL_FRONT, GL_SHININESS, &ms_ruby.shininess);
    glTranslated(mass[i].p[0]*SCALE , mass[i].p[1]*SCALE , mass[i].p[2]*SCALE );
    glutSolidSphere(MassR*SCALE, 20, 20);
    glPopMatrix();
  }

  // Springs
  vector<double> begin_pt;
  vector<double> end_pt;
  glMaterialfv(GL_FRONT, GL_DIFFUSE, green);
  glColor3d(0.0, 1.0, 0.0);//Color
  glLineWidth(3.0d);
  glBegin(GL_LINES);
  for(int i=0; i<N_SPRING; i++){
    begin_pt = mass[spring[i].masses[0]].p;
    end_pt = mass[spring[i].masses[1]].p;
    glVertex3d( begin_pt[0]*SCALE, begin_pt[1]*SCALE, begin_pt[2]*SCALE );
    glVertex3d( end_pt[0]*SCALE, end_pt[1]*SCALE, end_pt[2]*SCALE );
  }
  glEnd();
}


//----------------------------------------------------
// Static Objects function
//----------------------------------------------------

void StaticObjects(void){
  // Sphere
  glPushMatrix();
  glColor3d(1.0, 0.0, 0.0); //Color
  glMaterialfv(GL_FRONT, GL_AMBIENT, ms_ruby.ambient);
  glMaterialfv(GL_FRONT, GL_DIFFUSE, ms_ruby.diffuse);
  glMaterialfv(GL_FRONT, GL_SPECULAR, ms_ruby.specular);
  glMaterialfv(GL_FRONT, GL_SHININESS, &ms_ruby.shininess);
  glTranslated(0.0, 10.0, 20.0);//Translational position
  glutSolidSphere(4.0, 20, 20);//Radius, division around z-axis, division along z-axis)
  glPopMatrix();

  //Cube
  glPushMatrix();
  glMaterialfv(GL_FRONT, GL_DIFFUSE, green);
  glColor3d(0.0, 1.0, 0.0);//Color
  glTranslated(-20.0, 0.0, 20.0);//Translational position
  glutSolidCube(10.0);//Length of an edge
  glPopMatrix();

  //Cone
  glPushMatrix();
  glColor3d(0.0, 0.0, 1.0);//Color
  glMaterialfv(GL_FRONT, GL_DIFFUSE, blue);
  glTranslated(20.0, 100.0, 0.0);//Translational position
  glutSolidCone(5.0,10.0,20,20);//Radius, Height, Division around z-axis, Division along z-axis
  glPopMatrix();

  //Vertex
  glPushMatrix();
  glMaterialfv(GL_FRONT, GL_AMBIENT, ms_jade.ambient);
  glMaterialfv(GL_FRONT, GL_DIFFUSE, ms_jade.diffuse);
  glMaterialfv(GL_FRONT, GL_SPECULAR, ms_jade.specular);
  glMaterialfv(GL_FRONT, GL_SHININESS, &ms_jade.shininess);
  glColor3d(0.0, 1.0, 1.0);//Color
  glTranslated(30.0, 50.0, 0.0);//Translational position
  glBegin(GL_QUADS);
  for (int j = 0; j < 6; ++j) {
    glNormal3dv(normal[j]); //Normal vector
    for (int i = 0; i < 4; ++i) {
      glVertex3dv(vertex[face[j][i]]);
    }
  }
  glEnd();
  glPopMatrix();
}

// ---------------------------------------------------
// Pop Up Ball
// ---------------------------------------------------

void PopUpBall(void){
  for(int i=1; i<=pn; i++){
    glPushMatrix();
    glMaterialfv(GL_FRONT, GL_AMBIENT, ms_ruby.ambient);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, ms_ruby.diffuse);
    glMaterialfv(GL_FRONT, GL_SPECULAR, ms_ruby.specular);
    glMaterialfv(GL_FRONT, GL_SHININESS, &ms_ruby.shininess);
    glTranslated(p[i].x*SCALE , p[i].y*SCALE , p[i].z*SCALE );
    glutSolidSphere(BallR*SCALE, 20, 20);
    glPopMatrix();
  }
}

//----------------------------------------------------
// Draw Characters.
//----------------------------------------------------

void DISPLAY_TEXT(int x, int y, string str){
  static int list=0;

  glDisable(GL_LIGHTING);
  glDisable(GL_LIGHT0);

  glPushAttrib(GL_ENABLE_BIT);
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  gluOrtho2D(0, 100, 0, 100);
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  glColor3f(0.0, 0.0, 0.0);
  glCallList(list);
  glPopMatrix();
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glPopAttrib();
  glMatrixMode(GL_MODELVIEW);
  list=glGenLists(1);
  glNewList(list,GL_COMPILE);

  DRAW_STRING(x, y, str , GLUT_BITMAP_TIMES_ROMAN_24);
  glEndList();

glEnable(GL_LIGHTING);
glEnable(GL_LIGHT0);
}

void DRAW_STRING(int x, int y, string str, void *font){
  int len, i;
  glRasterPos2f(x, y);
  len = str.size();
  char t_char[len];
  str.copy(t_char, len);
  for (i = 0; i < len; i++){
    glutBitmapCharacter(font, t_char[i]);
  }
}

// -----------------------------------------------------------------------
// Frame capture -> png file
// -----------------------------------------------------------------------

void capture(const string fname)
{
    //string fname1 = "./png/.png";
    int len = fname.size();
    //cout << len << endl;
    char filepath[len];
    fname.copy(filepath, len);

    //for (int i; i<len; i++) cout << filepath[i];
    //cout << endl;
    //getchar();

    //const char filepath[] = "./output.png";

    png_bytep raw1D;
    png_bytepp raw2D;
    int i;
    int width = glutGet(GLUT_WINDOW_WIDTH);
    int height = glutGet(GLUT_WINDOW_HEIGHT);

    // Create structure.
    FILE *fp = fopen(filepath, "wb");
    png_structp pp = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    png_infop ip = png_create_info_struct(pp);

    // Pre-process
    png_init_io(pp, fp);
    png_set_IHDR(pp, ip, width, height,
        8, // 8bit
        PNG_COLOR_TYPE_RGBA, // RGBA
        PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);

    // Pixel area
    raw1D = (png_bytep)malloc(height * png_get_rowbytes(pp, ip));
    raw2D = (png_bytepp)malloc(height * sizeof(png_bytep));
    for (i = 0; i < height; i++)
        raw2D[i] = &raw1D[i*png_get_rowbytes(pp, ip)];

    // Capture.
    glReadBuffer(GL_FRONT);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1); // Initial value = 4.
    glReadPixels(0, 0, width, height,
            GL_RGBA, // RGBA
            GL_UNSIGNED_BYTE, // 8bit
            (void*)raw1D);

    // Flip.
    for (i = 0; i < height/ 2; i++){
        png_bytep swp = raw2D[i];
        raw2D[i] = raw2D[height - i - 1];
        raw2D[height - i - 1] = swp;
    }

    // Write.
    png_write_info(pp, ip);
    png_write_image(pp, raw2D);
    png_write_end(pp, ip);

    // Close.
    png_destroy_write_struct(&pp, &ip);
    fclose(fp);
    free(raw1D);
    free(raw2D);

    printf("write out screen capture to '%s'\n", filepath);
}

//------------------------------------------------------------------------
// Viewpoint rotation by mouse drag
//------------------------------------------------------------------------

// --- Drag mouse action
void mouse_motion(int x, int y){
  double dx, dy, distance;

  // Displacement of the mouse pointer from the drag start point.
  dx = (x - cx) * sx;
  dy = (y - cy) * sy;
  distance = sqrt(dx * dx + dy * dy);

  if( distance != 0.0 )
  {
    // Rotation by mouse drag (Quaternion dq).
    double ar = distance * 2.0 * PI * 0.5;
    double as = sin(ar) / distance;
    double dq[4] = { cos(ar), dy * as, dx * as, 0.0 };

    // Initial orientation cq * quarternion dq
    qmul(tq, dq, cq);

    // Quaternion -> Rotation matrix
    qrot(rt, tq);
  }
}

// --- On mouse action
void mouse_on(int button, int state, int x, int y)
{
  switch (button) {
  case 0:
    switch (state) {
    case 0:
      // Drag starting point
      cx = x;
      cy = y;
      break;
    case 1:
      // Save initial orientation.
      cq[0] = tq[0];
      cq[1] = tq[1];
      cq[2] = tq[2];
      cq[3] = tq[3];
      break;
    default:
      break;
    }
    break;
  default:
    break;
  }
  cout << x << " " << y<<endl;
}

// --- multiplication of quaternions r <- p x q
void qmul(double r[], const double p[], const double q[])
{
  r[0] = p[0] * q[0] - p[1] * q[1] - p[2] * q[2] - p[3] * q[3];
  r[1] = p[0] * q[1] + p[1] * q[0] + p[2] * q[3] - p[3] * q[2];
  r[2] = p[0] * q[2] - p[1] * q[3] + p[2] * q[0] + p[3] * q[1];
  r[3] = p[0] * q[3] + p[1] * q[2] - p[2] * q[1] + p[3] * q[0];
}

// --- transform matrix r <- quaternion q
void qrot(double r[], double q[]){
  double x2 = q[1] * q[1] * 2.0;
  double y2 = q[2] * q[2] * 2.0;
  double z2 = q[3] * q[3] * 2.0;
  double xy = q[1] * q[2] * 2.0;
  double yz = q[2] * q[3] * 2.0;
  double zx = q[3] * q[1] * 2.0;
  double xw = q[1] * q[0] * 2.0;
  double yw = q[2] * q[0] * 2.0;
  double zw = q[3] * q[0] * 2.0;

  r[ 0] = 1.0 - y2 - z2;
  r[ 1] = xy + zw;
  r[ 2] = zx - yw;
  r[ 4] = xy - zw;
  r[ 5] = 1.0 - z2 - x2;
  r[ 6] = yz + xw;
  r[ 8] = zx + yw;
  r[ 9] = yz - xw;
  r[10] = 1.0 - x2 - y2;
  r[ 3] = r[ 7] = r[11] = r[12] = r[13] = r[14] = 0.0;
  r[15] = 1.0;
}

void mouse_wheel(float z){
  camera_z_pos += z;

  if( camera_z_pos < 0.0f )
  {
    camera_z_pos = 0.0f;
  }
}

double calc_distance(vector<double>p0, vector<double>p1){
    double length;
    double sum = 0.0;
    for (int i=0; i<3; i++){
        sum += (p0[i]- p1[i]) * (p0[i]- p1[i]);
    }
    length = sqrt(sum);
    return length;
}
