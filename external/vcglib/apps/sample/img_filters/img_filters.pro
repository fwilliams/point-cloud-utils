# debugging
CONFIG += debug

# Base options
TEMPLATE = app
LANGUAGE = C++

# Executable name
TARGET = img_filters

# STL support is enabled
CONFIG += stl

# enable console
CONFIG += console

# Awful..
win32{
  DEFINES += NOMINMAX
}

# The following define is needed in gcc to remove the asserts
win32-g++:DEFINES += NDEBUG
CONFIG(debug, debug|release) {
	win32-g++:release:DEFINES -= NDEBUG
}

#include current path
INCLUDEPATH += . 

#include lib path
INCLUDEPATH += ../../..

SOURCES += img_filters.cpp
