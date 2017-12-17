/*
  A program which opens a window and draws the "color cube."

  Use the arrow keys and +/-, PgUp,PgDn, Home,End, Ins,Del 
  keys to make the cube move.

  Use w/W to toggle between wireframe and solid models
  Use c/C to toggle backface culling on/off

  Written by - 
               Parag Chaudhuri
	       Research Scholar,
	       DCSE, IITD.
*/

#include <cstdlib>
//Including only glut.h is sufficient as it includes GL.h and GLX.h for us
#include <GLUT/glut.h>
#include "SplatpointViewer.h"
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

inline float reverseFloat (char *c) 
{
    float i;
    char *p = (char *)&i;
	
    p[0] = c[3];
    p[1] = c[2];
    p[2] = c[1];
    p[3] = c[0];
	
    return i;
}


void ReadFile(const char* modelFile)
{
    string filename = "armadillo.points";
    
    if(modelFile)
    {
        filename = modelFile;
    }
/*        rfile.open(modelFile);
    else
        rfile.open(filename.c_str());
  */      
    ifstream rfile(filename.c_str());
    if(!rfile)
    {
        cout << "File : " << filename << " could not be opened " << endl;
    }
    
    // Calculate file size and splatpoint count
    rfile.seekg(0, ios::end);
    size_t size = rfile.tellg();
    size /= sizeof(splatpoint);
    rfile.seekg(0, ios::beg);
	//size /= 10.0;
	cout << size << " " << sizeof(splatpoint) << endl;
	
    // Malloc array and read file
    // 10 attributes in each splatpoint
    dataArray = new float[10*size];
    
    minX = dataArray[0];
    minY = dataArray[1];
    minZ = dataArray[2];
    
    maxX = dataArray[0];
    maxY = dataArray[1];
    maxZ = dataArray[2];  
    
    for(size_t i=0; i<size; i++)
    {
        splatpoint sp;
        rfile.read((char*)&sp, sizeof(splatpoint));
        
		sp.x = reverseFloat((char*)&sp.x);
		sp.y = reverseFloat((char*)&sp.y);
		sp.z = reverseFloat((char*)&sp.z);
		sp.rad = reverseFloat((char*)&sp.rad);
		
		sp.nx = reverseFloat((char*)&sp.nx);
		sp.ny = reverseFloat((char*)&sp.ny);
		sp.nz = reverseFloat((char*)&sp.nz);
		
        dataArray[4*i] = sp.x;
        dataArray[4*i+1] = sp.y;
        dataArray[4*i+2] = sp.z;
        dataArray[4*i+3] = sp.rad;
        
        dataArray[4*size + 3*i] = sp.nx;
        dataArray[4*size + 3*i+1] = sp.ny;
        dataArray[4*size + 3*i+2] = sp.nz;
        
        minX = fmin(minX, dataArray[4*i]);
        minY = fmin(minY, dataArray[4*i+1]);
        minZ = fmin(minZ, dataArray[4*i+2]);
        
        maxX = fmax(maxX, dataArray[4*i]);
        maxY = fmax(maxY, dataArray[4*i+1]);
        maxZ = fmax(maxZ, dataArray[4*i+2]);
    }
    rfile.close();
	
	cout<<minX<<" "<<minY<<" "<<minZ<<" "<<maxX<<" "<<maxY<<" "<<maxZ<<endl;
    
	float midX, midY, midZ;
	
	midX = (minX + maxX)/2.0;
	midY = (minY + maxY)/2.0;
	midZ = (minZ + maxZ)/2.0;
	
   for(size_t i=0; i<size; i++)
    {
		dataArray[4*i] -= midX;
        dataArray[4*i+1] -= midY;
        dataArray[4*i+2] -= midZ;
	}

	lookAt[0] = midX;
	lookAt[1] = midY;
	lookAt[2] = midZ;
	
	camera[0] = midX; 
	camera[1] = midY; 
    camera[2] = maxZ + 500; 
	
    cout<<minX<<" "<<minY<<" "<<minZ<<endl;;
    pointCount = size;
}


uint32_t ReadFileIntoString(const string& filename, unsigned char*& str)
{
	uint32_t size;
	ifstream rfile(filename.c_str());
	rfile.seekg(0,ios::end);
	size = rfile.tellg();
	rfile.seekg(0,ios::beg);
	rfile.read((char*)str, size);
	rfile.close();
	return size;
}

void CreateVBO()
{
    float ratio = ceil((float)pointCount/(float)maxVBOSize);
    vboCount = ratio;
    
    vboCount = min((int)ratio, (int)maxVBOs);
    
    size_t start = 0;
    size_t localCount = maxVBOSize;
    
    for(size_t i=0; i<vboCount; i++)
    {
        localCount = ((start + maxVBOSize) >= pointCount ) ? (pointCount-start) :
                                                        maxVBOSize;
        
        vboSizes[i] = localCount;
        vboStart[i] = start;                                                
         
        glGenBuffers(1, &vboId[i]);
        glBindBufferARB( GL_ARRAY_BUFFER_ARB, vboId[i] );
        glBufferDataARB( GL_ARRAY_BUFFER_ARB, 4*localCount*sizeof(GL_FLOAT), 
                                    &dataArray[4*start], GL_STATIC_DRAW_ARB );
    
        glGenBuffersARB( 1, &nboId[i] );
        glBindBufferARB( GL_ARRAY_BUFFER_ARB, nboId[i] );
        glBufferDataARB( GL_ARRAY_BUFFER_ARB, 3*localCount*sizeof(GL_FLOAT), 
                                    &dataArray[4*pointCount+3*start], GL_STATIC_DRAW_ARB );
        
        if(hasColors)
        {
            /*glGenBuffersARB( 1, &nboId[i] );
            glBindBufferARB( GL_ARRAY_BUFFER_ARB, nboId[i] );
            glBufferDataARB( GL_ARRAY_BUFFER_ARB, 3*localCount*sizeof(GL_FLOAT), 
                                         &dataArray[4*pointCount+3*start], GL_STATIC_DRAW_ARB );*/
        }
        
        start += localCount;
    }
    
    if(dataArray)
        delete[] dataArray;
    cout << vboCount << "VBOs created for point count " << pointCount << endl;
}

void DeleteVBO()
{
    for(size_t i=0; i<vboCount; i++)
    {
        if(vboId[i] > 0)
            glDeleteBuffers(1, &vboId[i]);
            
        if(nboId[i] > 0)
            glDeleteBuffers(1, &nboId[i]);
            
        if(cboId[i] > 0)
            glDeleteBuffers(1, &cboId[i]);
    }   
}

void InitShadersGLSL()
{
	unsigned char* str1 = new unsigned char[65536];
	unsigned char* str2 = new unsigned char[65536];
	char message[1024];
	uint32_t size=0;
	GLint vertCompiled, fragCompiled;
	int logsize, flag;
	
	vs_glsl = glCreateShader(GL_VERTEX_SHADER);
	size = ReadFileIntoString("vertex.glsl", str1);
	const char* vv = (char*)str1;
	glShaderSource(vs_glsl, 1, &vv,NULL);
	
	glCompileShader(vs_glsl);
	glGetShaderiv(vs_glsl, GL_COMPILE_STATUS, &vertCompiled);
	if(!vertCompiled)
	{
		cout << "Vertex shader failed to compiled " << endl;
		glGetShaderiv(vs_glsl, GL_INFO_LOG_LENGTH, &logsize);
		glGetShaderInfoLog(vs_glsl, logsize, &flag, message);
		cout << message << endl;
		exit(1);
	}
	
	fs_glsl = glCreateShader(GL_FRAGMENT_SHADER);
 	size = ReadFileIntoString("fragment.glsl", str2);
	const char* ff = (char*)str2;
	glShaderSource(fs_glsl, 1, &ff,NULL);
	
	glCompileShader(fs_glsl);
	glGetShaderiv(fs_glsl, GL_COMPILE_STATUS, &fragCompiled);
	if(!fragCompiled)
	{
		cout << "Fragment shader failed to compiled " << endl;
		glGetShaderiv(fs_glsl, GL_INFO_LOG_LENGTH, &logsize);
		glGetShaderInfoLog(fs_glsl, logsize, &flag, message);
		cout << message << endl;
		exit(1);
	}
	
	sp_glsl = glCreateProgram();
	
	glAttachShader(sp_glsl, vs_glsl);
	glAttachShader(sp_glsl, fs_glsl);
	
	GLint status=1;
	glLinkProgram(sp_glsl);
	glGetProgramiv( sp_glsl, GL_LINK_STATUS, &status );
	if(!status)
	{
		cout << "Failed to link shader program" << endl;
		glGetShaderiv(sp_glsl, GL_INFO_LOG_LENGTH, &logsize);
		glGetShaderInfoLog(sp_glsl, logsize, &flag, message);
		cout << message << " " << flag << endl;
		exit(1);
	}
	glUseProgram(sp_glsl);
	
	CameraVector_ = glGetUniformLocation( sp_glsl, "cameravector" );
	KFOV_ = glGetUniformLocation( sp_glsl, "KFOV" );

	
	delete[] str1;
	delete[] str2;
}

void InitLightsGL()
{
   GLfloat light_ambient[4] = {0.3, 0.3, 0.3, 1.0};
   GLfloat light_diffuse[4] = {0.8, 0.8, 0.8, 1.0};
   GLfloat light_specular[4] = {0.7, 0.7, 0.7, 1.0};
   GLfloat light_position[4] = {0.0, 0.0, 1.0, 0.0};
   GLfloat mat_amb[4] = {0.65, 0.65, 0.65, 1.0};
   GLfloat mat_diffuse[4] = {0.7, 0.7, 0.7, 1.0};
   GLfloat mat_specular[4] = { 0.5, 0.5, 0.5, 1.0 };
   GLfloat mat_shininess[1] = {60.0};

   glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, mat_amb);
   glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_diffuse);
   glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
   glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);
   glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
   glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
   glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
   glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	
   glEnable(GL_LIGHTING);
   glEnable(GL_LIGHT0);
   //glEnable(GL_SMOOTH);

   glEnable(GL_DEPTH_TEST);
   //glEnable(GL_COLOR_MATERIAL);
   glDisable(GL_LIGHTING);
   glDisable(GL_COLOR_MATERIAL);
}



GLvoid UpdateCamera()
{

	float modelview[16] = {0};
	float matrix[16] = {0};
	glGetFloatv( GL_MODELVIEW_MATRIX, modelview );
		
	// calculations of the original matrix.
	matrix[0]  = modelview[0]; matrix[1] = modelview[4]; matrix[2]  = modelview[8];
	matrix[4]  = modelview[1]; matrix[5] = modelview[5]; matrix[6]  = modelview[9];
	matrix[8]  = modelview[2]; matrix[9] = modelview[6]; matrix[10] = modelview[10];
	matrix[3]  = 0.0f; matrix[7] = 0.0f; matrix[11] = 0.0f;
	matrix[15] = 1.0f;

	matrix[12] = -(modelview[12] * modelview[0]) - (modelview[13] * modelview[1]) - (modelview[14] * modelview[2]);
	matrix[13] = -(modelview[12] * modelview[4]) - (modelview[13] * modelview[5]) - (modelview[14] * modelview[6]);
	matrix[14] = -(modelview[12] * modelview[8]) - (modelview[13] * modelview[9]) - (modelview[14] * modelview[10]);
	
	float x = matrix[0]*camera[0] + matrix[4]*camera[1] + matrix[8]*camera[2] + matrix[12];
	float y = matrix[1]*camera[0] + matrix[5]*camera[1] + matrix[9]*camera[2] + matrix[13];
	float z = matrix[2]*camera[0] + matrix[6]*camera[1] + matrix[10]*camera[2] + matrix[14];
	float w = matrix[3]*camera[0] + matrix[7]*camera[1] + matrix[11]*camera[2] + matrix[15];	
	if( fabs(w) < 0.0001 ) 
		w = 1.0;
		
	cameraN[0] = x/w;
	cameraN[1] = y/w;
	cameraN[2] = z/w;
}


//Function to initialize some basic parameters
GLvoid InitGL(const char* modelFile)
{
  //Setting the color used to clear the framebuffer
  glClearColor(0.0, 0.0, 0.0, 1.0f);
  //Setting the depth used to clear the depthbuffer
  glClearDepth(1.0);
  //Enabling Z-buffering
  glEnable(GL_DEPTH_TEST);
  //Enabling Smooth Shading 
  glShadeModel(GL_SMOOTH);
  //Setup backface culling
  glCullFace(GL_BACK);
  
  ReadFile(modelFile);
  InitShadersGLSL();
  CreateVBO();
  
  InitLightsGL();
}

//-------------------------------------------------------------------------

//GL reshape callback
GLvoid ReshapeGL (GLsizei Width, GLsizei Height)
{
  //Prevent a divison by zero
  if (Height == 0) Height=1;
  //Set the viewport to the whole of the current window
  glViewport (0, 0, Width, Height);
  //Change the matrix mode to projection
  glMatrixMode (GL_PROJECTION);
  //Load the identity matrix
  glLoadIdentity ();
  //Setup a perspective projection
  //Parameters are in the order -
  //vertical field of view angle - 60 degrees
  //Aspect Ratio i.e. width/height of window - 1.0
  //Near clipping plane distance - 0.1
  //Far clipping plane distance - 1000.0
  gluPerspective(45.0f,(GLfloat)Width/(GLfloat)Height, 1.0, 10000.0);
  //Change the matrix mode to modelview
  glMatrixMode (GL_MODELVIEW);
  
  widthW = Width;
  heightW = Height;
}

//-------------------------------------------------------------------------

//GL display callback - does all the drawing
GLvoid DisplayGL(GLvoid)
{
  //Clear the frame and depth buffers
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  //Setup the camera
  //Camera center at 0.0,5.0,0.0
  //LookAt 0.0,0.0,0.0
  //Up Vector 0.0,1.0,0.0
  gluLookAt(camera[0], camera[1], camera[2],
            lookAt[0], lookAt[1], lookAt[2], 0.0,1.0,0.0);
  //Translate as per the current translation parameters 
  glTranslatef(xpos,ypos,zpos);
  //Rotate as per the current rotation parameters
  glRotatef(xrot,1.0,0.0,0.0);
  glRotatef(yrot,0.0,1.0,0.0);
  glRotatef(zrot,0.0,0.0,1.0);

  //Toggle Backface culling
  /*if (enable_culling) 
    glEnable(GL_CULL_FACE);
  else
    glDisable(GL_CULL_FACE);

  //Toggle solid/wireframe drawing
  if(solid)
    glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
  else
    glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
  
  //Draw the color cube
  glDrawElements(GL_QUADS, 24, GL_UNSIGNED_BYTE, allIndices);*/
  
  glEnableClientState( GL_VERTEX_ARRAY );
  glEnableClientState( GL_NORMAL_ARRAY );
  
  glEnable( GL_VERTEX_PROGRAM_POINT_SIZE );
	
  float f = widthW / tan( 0.5 * 0.0174603174603 * 45.0f );
  glUseProgram(sp_glsl);
  glUniform3fARB(CameraVector_, cameraN[0], cameraN[1], cameraN[2]);
  glUniform4fARB(KFOV_, f, f, f, f);
 
  
  for(size_t i=0; i<vboCount; i++)
  {
      assert(vboId[i] > 0);
      glBindBufferARB( GL_ARRAY_BUFFER_ARB, vboId[i] );
      glVertexPointer( 4, GL_FLOAT, 0, NULL );

      assert(nboId[i] > 0);
      glBindBufferARB( GL_ARRAY_BUFFER_ARB, nboId[i] );
      glNormalPointer( GL_FLOAT, 0, NULL );
      
      glDrawArrays( GL_POINTS, 0, vboSizes[i]);
  }
 
  
  glDisableClientState( GL_VERTEX_ARRAY );
  glDisableClientState( GL_NORMAL_ARRAY );
  
  UpdateCamera();
  
  //Swap the double buffers
  glutSwapBuffers();
}

//-------------------------------------------------------------------------

//GL keyboard callback
GLvoid KeyPressedGL(unsigned char key, GLint x, GLint y) 
{
    switch (key)
      {    
	//quit
        case ESCAPE: 
                exit(1);                   	
	             break;
	//move along Z
        case '+' :zpos+=quanta; 
	          glutPostRedisplay();
		  break;
	//move along Z
        case '-' :zpos-=quanta; 
	          glutPostRedisplay();
		  break;
	//rotate left about Z
        case DELETE:yrot-=1.0;
	            glutPostRedisplay();
		    break;
	//toggle wireframe/solid
        case 'w':
        case 'W':solid = !solid;
	         glutPostRedisplay();
		 break;
	//Toggle culling
        case 'c':
        case 'C':enable_culling = !enable_culling;
	         glutPostRedisplay();
		 break;
         default:
		  break;
      }	

}

//-------------------------------------------------------------------------

//GL Keyboard callback for special keys
GLvoid SpecialKeyPressedGL(GLint key, GLint x, GLint y) 
{
 switch (key)
      {   
       //Move along Y
       case GLUT_KEY_UP:ypos+=quanta; 
	                glutPostRedisplay();
		        break;
       //Move along Y
       case GLUT_KEY_DOWN:ypos-=quanta; 
	                  glutPostRedisplay();
		          break;
       //Move along X
       case GLUT_KEY_LEFT:xpos+=quanta; 
	                glutPostRedisplay();
		        break;
       //Move along X
       case GLUT_KEY_RIGHT:xpos-=quanta; 
	                  glutPostRedisplay();
		          break;
       //Rotate right about Z
       case GLUT_KEY_PAGE_DOWN:yrot+=1.0;
	                    glutPostRedisplay();
			    break;
       //Rotate about X
       case GLUT_KEY_HOME:xrot+=1.0;
	                  glutPostRedisplay();
			  break;
       //Rotate about X
       case GLUT_KEY_END:xrot-=1.0;
	                    glutPostRedisplay();
			    break;
       //Rotate about Y
       case GLUT_KEY_INSERT:zrot+=1.0;
	                    glutPostRedisplay();
			    break;
       //Rotate about Y
       case GLUT_KEY_PAGE_UP:zrot-=1.0;
	                    glutPostRedisplay();
			    break;
       default:
	        break;
      }	
}

GLvoid MouseCallBackGL( int button, int state, int x, int y )
{
	switch( button )
	{
		case GLUT_LEFT_BUTTON:
		if( state == GLUT_DOWN )
		{
			_mstate = LEFT;
			startX = x;
			startY = y;
			move = true;
		}
		else
		{
			_mstate = NONE;
			move = false;
		}
		glutPostRedisplay();
		break;
		
		case GLUT_MIDDLE_BUTTON:
		if( state == GLUT_DOWN )
		{
			_mstate = MIDDLE;
			startX = x;
			startY = y;
			move = true;
		}
		else
		{
			_mstate = NONE;
			move = false;
		}
		break;
		
		case GLUT_RIGHT_BUTTON:
		if( state == GLUT_DOWN )
		{
			_mstate = RIGHT;
			move = true;
			
			startX = x;
			startY = y;
		}
		else
		{
			_mstate = NONE;
			move = false;
		}	
			
		break;
		
		default:
		_mstate = NONE;
		break;
	}
}

GLvoid MotionCallBackGL( int x, int y )
{

	if( move )
	{
		if( _mstate == LEFT )
		{
			zrot += (float)( x- startX ) / 10.0;
			xrot += (float)( y- startY ) / 10.0;

		}
		else if( _mstate == MIDDLE )
		{
		    //yrot += (float)( x - startX ) / 5.0;
			xpos += (float)( x - startX ) / 2.0;
			ypos -= (float)( y - startY ) / 2.0;
		    
		    startX = x;
    		startY = y;
		}
		else if( _mstate == RIGHT )
		{
			zpos += 1.0*(y-startY);

		}
			
        startX = x;
		startY = y;

	}
	glutPostRedisplay();
}

//-------------------------------------------------------------------------

// The main rendering function
GLvoid RenderGL(int argc, char** argv)
{       
  //Initialize GLUT
  glutInit(&argc, argv);
  //Initialize GLUT display
  glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH  | GLUT_STENCIL);
  //Window top-left corner position
  glutInitWindowPosition(0,0 );
  //Window size
  glutInitWindowSize( widthW, heightW );
  //Window title
  glutCreateWindow("The Color Cube");

  //Our Init function
  if(argc > 1)
      InitGL(argv[1]);
  else
      InitGL(NULL);
      
  if(argc > 2)
      hasColors = atoi(argv[2]);

  //Register the callbacks  
  glutDisplayFunc(&DisplayGL);
  glutReshapeFunc(&ReshapeGL);  
  glutKeyboardFunc(&KeyPressedGL);
  glutSpecialFunc(&SpecialKeyPressedGL);
  glutMouseFunc( &MouseCallBackGL );	
  glutMotionFunc( &MotionCallBackGL );
  
  //Start the GLUT event handling loop
  glutMainLoop();
}

//-------------------------------------------------------------------------

int main(int argc, char** argv)
{   
  RenderGL(argc,argv);
}

//-------------------------------------------------------------------------

