#DEFINES += VCG_USE_EIGEN
TARGET = texdeci
#LIBPATH +=
DEPENDPATH += . ../..
INCLUDEPATH += . ../..
CONFIG += console stl debug_and_release
TEMPLATE = app
HEADERS += 
SOURCES += texdeci.cpp ../../wrap/ply/plylib.cpp

QT-=gui


# Mac specific Config required to avoid to make application bundles
CONFIG -= app_bundle

