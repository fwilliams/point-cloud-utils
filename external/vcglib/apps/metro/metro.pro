
TARGET = metro
INCLUDEPATH += . ../..
CONFIG += console stl
TEMPLATE = app
HEADERS += sampling.h
SOURCES += metro.cpp ../../wrap/ply/plylib.cpp

# Mac specific Config required to avoid to make application bundles
CONFIG -= app_bundle
