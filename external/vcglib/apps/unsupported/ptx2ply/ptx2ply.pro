
TARGET = ptx2ply
LIBPATH += 
DEPENDPATH += . 
INCLUDEPATH += . ../..
CONFIG += console stl
TEMPLATE = app
#HEADERS += defs.h sampling.h mesh_type.h
SOURCES += ptx2ply.cpp ../../wrap/ply/plylib.cpp
# Mac specific Config required to avoid to make application bundles
CONFIG -= app_bundle