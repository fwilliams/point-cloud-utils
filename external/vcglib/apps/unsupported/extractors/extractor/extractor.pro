INCLUDEPATH += . ../../..
CONFIG += console stl
TEMPLATE = app

TARGET = 
DEPENDPATH += .
INCLUDEPATH += .

# Input
HEADERS += Definitions.h \
           Implicit.h \
           ImplicitSphere.h \
           SphereDifference.h \
           SphereUnion.h \
           Volume.h \
           Walker.h
SOURCES += main.cpp
