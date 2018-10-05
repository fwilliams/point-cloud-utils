INCLUDEPATH += . ../../.. ../../../../code/lib ../../../../code/lib/glew/include
HEADERS       = glwidget.h \
                window.h \
		mesh_type.h	
SOURCES       = glwidget.cpp \
                main.cpp \
                window.cpp\
		 ../../../../code/lib/glew/src/glew.c \
		../../../wrap/ply/plylib.cpp\
		../../../wrap/gui/trackmode.cpp\
		../../../wrap/gui/trackball.cpp
QT           += opengl

# install
target.path = $$./debug
sources.files = $$SOURCES $$HEADERS $$RESOURCES $$FORMS trimesh_pos_demo.pro
sources.path =  ./
INSTALLS += target sources

