INCLUDEPATH += ../../../
CONFIG += console stl c++11
CONFIG -= qt
TEMPLATE = app
SOURCES += main.cpp

# Mac specific Config required to avoid to make application bundles
CONFIG -= app_bundle
