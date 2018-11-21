#ifndef _OPENGL_H_
#define _OPENGL_H_

#include <string>
#include <vector>
using std::string;
using std::vector;

void Initialize(void);
void Display(void);
void Idle(void);
void Ground(void);
void StaticObjects(void);
void PopUpBall(void);
void Keyboard(unsigned char, int, int);
void SpecialKey(int, int, int);
void DRAW_STRING(int, int, string, void*);
void DISPLAY_TEXT(int, int, string);
void capture(const string);
void qmul(double[], const double[], const double[]);
void qrot(double[], double[]);
void mouse_motion(int, int);
void mouse_on(int, int, int, int);
void CubicRobot(void);
void CubicRobotSolid(void);
void fixed_view(int);
void glDrawAxisd(double);
void getMatrix(void);

// Window Setting
extern const int WindowPositionX;  //Window x position
extern const int WindowPositionY;  //Window y position
extern const int WindowWidth;      //Window width
extern const int WindowHeight;     //Window height
extern const char WindowTitle[];   //Window title

// Control
extern bool _Friction;
extern bool _Damping;
extern bool _Breathe;

extern int frame;
extern double dt_ms;

extern vector<double> nudgeForce;
extern int nudgeID;
#endif
