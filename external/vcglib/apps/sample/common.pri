DEPENDPATH += . ../../..
INCLUDEPATH += . ../../..
CONFIG += console c++11
TEMPLATE = app
# Mac specific Config required to avoid to make application bundles
CONFIG -= app_bundle

QMAKE_CXXFLAGS += -std=c++11
