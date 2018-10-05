VCGLIBDIR = ../../../../vcglib
GLEWDIR   = ../../../../code/lib/glew/
ANTDIR    = ../../../../code/lib/AntTweakBar1.14/

HEADERS       = glwidget.h
SOURCES       = glwidget.cpp \
                main.cpp
QT           += opengl

# Compile glew
DEFINES += GLEW_STATIC
INCLUDEPATH += $$GLEWDIR/include
SOURCES += $$GLEWDIR/src/glew.c

INCLUDEPATH += $$VCGLIBDIR
INCLUDEPATH += $$GLEWDIR/include
INCLUDEPATH += $$ANTDIR/include

SOURCES += $$VCGLIBDIR/wrap/ply/plylib.cpp
SOURCES += $$VCGLIBDIR/wrap/gui/trackball.cpp
SOURCES += $$VCGLIBDIR/wrap/gui/trackmode.cpp
SOURCES += $$VCGLIBDIR/wrap/qt/anttweakbarMapper.cpp


# Awful problem with windows..
win32{
  DEFINES += NOMINMAX
  LIBS +=$$ANTDIR/lib/AntTweakBar.lib
}

mac{
# Mac specific Config required to avoid to make application bundles
  CONFIG -= app_bundle
  LIBS +=$$ANTDIR/lib/libAntTweakBar.dylib
  QMAKE_POST_LINK ="cp -P ../../../../code/lib/AntTweakBar1.14/lib/libAntTweakBar.dylib . ; install_name_tool -change ../lib/libAntTweakBar.dylib ./libAntTweakBar.dylib $$TARGET"
}

unix:!macx{
  LIBS +=$$ANTDIR/lib/libAntTweakBar.so -lGLU
}
