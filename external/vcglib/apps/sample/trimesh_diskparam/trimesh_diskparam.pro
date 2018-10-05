#DEFINES += VCG_USE_EIGEN

TARGET = trimesh_diskparam
DEPENDPATH += . ../../..

INCLUDEPATH += . ../../..

CONFIG += console stl
TEMPLATE = app
SOURCES += trimesh_diskparam.cpp ../../../wrap/ply/plylib.cpp

# Mac specific Config required to avoid to make application bundles
CONFIG -= app_bundle
