include(../common.pri)
QT += opengl svg
TARGET = trimesh_texture
SOURCES += trimesh_texture.cpp  \
../../../../vcglib/wrap/qt/Outline2ToQImage.cpp \
../../../../vcglib/wrap/qt/outline2_rasterizer.cpp \
../../../wrap/ply/plylib.cpp
