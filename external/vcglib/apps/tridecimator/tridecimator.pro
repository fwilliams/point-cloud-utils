
TARGET = tridecimator
DEPENDPATH += ../..
INCLUDEPATH += . ../..
CONFIG += console stl  c++11 debug_and_release
TEMPLATE = app
HEADERS += 
SOURCES += tridecimator.cpp ../../wrap/ply/plylib.cpp


# Mac specific Config required to avoid to make application bundles
CONFIG -= app_bundle 
