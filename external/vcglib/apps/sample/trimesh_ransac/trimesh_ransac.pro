include(../common.pri)
TARGET = trimesh_ransac 
SOURCES +=  ../../../wrap/ply/plylib.cpp \
  ../../../vcg/complex/algorithms/ransac_matching.h \
  trimesh_ransac.cpp
