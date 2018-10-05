include(../common.pri)
TARGET = trimesh_harmonic
SOURCES += trimesh_harmonic.cpp ../../../wrap/ply/plylib.cpp

CONFIG(release, debug | release){
DEFINES += NDEBUG
}

