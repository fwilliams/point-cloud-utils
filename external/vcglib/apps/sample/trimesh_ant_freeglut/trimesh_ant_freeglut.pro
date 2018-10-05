# Base options
TEMPLATE = app
LANGUAGE  = C++

# Executable name
TARGET = trimesh_ant_freeglut

# Directories
DESTDIR = .
OBJECTS_DIR = build/obj

# Lib headers
INCLUDEPATH += .
INCLUDEPATH += ../../..

# Lib sources
SOURCES += ../../../wrap/ply/plylib.cpp
SOURCES += ../../../wrap/gui/trackball.cpp
SOURCES += ../../../wrap/gui/trackmode.cpp


# Compile glew
DEFINES += GLEW_STATIC
SOURCES += ../../../../code/lib/glew/src/glew.c

# Awful problem with windows..
win32{
  DEFINES += NOMINMAX
  INCLUDEPATH += ../../../../code/lib/glew/include
  INCLUDEPATH += ../../../../code/lib/AntTweakBar/include
  INCLUDEPATH += ../../../../code/lib/freeglut/include
  LIBS +=../../../../code/lib/AntTweakBar/lib/AntTweakBar.lib
  LIBS +=../../../../code/lib/freeglut/lib/freeglut.lib
}
unix {
  LIBS += -lfreeglut
  }

# Input
SOURCES += main.cpp
