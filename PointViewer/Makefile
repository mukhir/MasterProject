CC=g++

GL_INC_DIR=/opt/local/include
GL_LIB_DIR=/opt/local/lib
X_LIB_DIR=/usr/X11R6/lib

GL_LIBS=-L$(GL_LIB_DIR) -framework GLUT -framework OpenGL
X_LIBS=-L$(X_LIB_DIR) -lXmu -lX11

CFLAGS = -Wall 

.PHONY:clean

OBJ = SplatpointViewer.o
LIBS = -L$(GL_LIBS)

all: $(OBJ) SplatpointViewer

$(OBJ): %.o : %.cpp
	$(CC) -c $(CFLAGS) $< -o $@

SplatpointViewer: $(OBJ) 
	$(CC) $(OBJ) $(LIBS) -o $@

clean:
	rm -f  ./*~ ./core $(OBJ) 

#DO NOT DELETE
SplatpointViewer.o : SplatpointViewer.h
