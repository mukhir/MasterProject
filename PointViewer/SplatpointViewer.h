/*
  
*/
#ifndef COLORCUBE_H
#define COLORCUBE_H

// Defining the ESCAPE Key Code
#define ESCAPE 27
// Defining the DELETE Key Code
#define DELETE 127


enum MouseState { NONE, LEFT, MIDDLE, RIGHT };

// Translation Parameters
GLfloat xpos=0.0,ypos=0.0,zpos=0.0;
// Rotation Parameters
GLfloat xrot=0.0,yrot=0.0,zrot=0.0;

//Running variable to toggle culling on/off
bool enable_culling=true;
//Running variable to toggle wireframe/solid modelling
bool solid=true;



struct splatpoint
{
    float x, y, z;
    float rad;
    float nx, ny, nz;
    float r, g, b;
    
    void operator = (const splatpoint& sp)
    {
        x = sp.x;
        y = sp.y;
        z = sp.z;
        rad = sp.rad;
        nx = sp.nx;
        ny = sp.ny;
        nz = sp.nz;
        r = sp.r;
        g = sp.g;
        b = sp.b;
    }  
};

GLint startX = 0, startY = 0;
MouseState _mstate;
bool move = false;


GLuint displayList = 0;
size_t pointCount = 0;
float* dataArray;

GLuint vs_glsl;
GLuint fs_glsl;
GLuint sp_glsl;
GLint CameraVector_;
GLint KFOV_;

const static size_t maxVBOSize = 150000;
const static size_t maxVBOs = 100; 

float camera[] = {0.0, 0.0, 10.0};
float cameraN[] = {1.0, 1.0, 1.0};
float lookAt[] = {0.0, 0.0, 0.0};

GLuint vboId[maxVBOs] = {0}, nboId[maxVBOs] = {0}, cboId[maxVBOs] = {0};
size_t vboSizes[maxVBOs] = {0};
size_t vboStart[maxVBOs] = {0};

size_t vboCount;

int widthW = 1024, heightW = 1024;

float minX, minY, minZ;
float maxX, maxY, maxZ;

float quanta = 10.0;

bool hasColors = false;
//-------------------------------------------------------------------------

void InitShadersGLSL();
void ReadFile(const char* modelFile = NULL);
void Destructor();
void CreateVBO();
void DeleteVBO();
void InitLightsGL();

GLvoid InitGL(const char* modelFile = NULL);
GLvoid ReshapeGL (GLsizei Width, GLsizei Height);
GLvoid DisplayGL(GLvoid);
GLvoid KeyPressedGL(unsigned char key, GLint x, GLint y);
GLvoid SpecialKeyPressedGL(GLint key, GLint x, GLint y);
GLvoid RenderGL(int argc, char** argv);

//-------------------------------------------------------------------------

#endif
