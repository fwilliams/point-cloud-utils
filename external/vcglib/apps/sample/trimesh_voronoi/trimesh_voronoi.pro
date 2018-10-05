include(../common.pri)
TARGET = trimesh_voronoi
SOURCES += trimesh_voronoi.cpp ../../../wrap/ply/plylib.cpp

CONFIG(release, debug | release){
DEFINES += NDEBUG
}

