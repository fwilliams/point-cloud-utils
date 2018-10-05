include(../common.pri)
TARGET		 = kdTree_test

HEADERS		 =  nanoflann.hpp

SOURCES		 =  trimesh_indexing.cpp \
				../../../wrap/ply/plylib.cpp

win32-msvc2010:QMAKE_CXXFLAGS   += /openmp
win32-msvc2012:QMAKE_CXXFLAGS   += /openmp

win32-g++:QMAKE_CXXFLAGS    += -fopenmp
win32-g++:QMAKE_LIB   		+= -lgomp

mac-g++:QMAKE_CXXFLAGS    	+= -fopenmp
mac-g++:QMAKE_LIB   		+= -lgomp


				
win32{
  DEFINES += NOMINMAX
}