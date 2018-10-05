# IMPORTANT
# We assume that you have unzipped and compiled somewhere the COMISO solver
# put here that path and and the name of the build dir (from which the lib are copied).
COMISODIR = ../../../code/lib/CoMISo
COMISOBUILDDIR = $$COMISODIR/buildMACOSX

TEMPLATE = app
TARGET = quadrangulator

# This is needed by CoMISo
DEFINES += INCLUDE_TEMPLATES
INCLUDEPATH += $$COMISODIR/include
INCLUDEPATH += $$COMISODIR/Solver
INCLUDEPATH += $$COMISODIR/..

# just a shortcut
VCGLIBDIR = ../../../vcglib

# include of vcg library
INCLUDEPATH += $$VCGLIBDIR

# CORE
HEADERS += $$VCGLIBDIR/wrap/miq/core/vertex_indexing.h
HEADERS += $$VCGLIBDIR/wrap/miq/core/poisson_solver.h
HEADERS += $$VCGLIBDIR/wrap/miq/quadrangulator.h
HEADERS += $$VCGLIBDIR/wrap/miq/core/param_stats.h
HEADERS += $$VCGLIBDIR/wrap/miq/MIQ.h
HEADERS += $$VCGLIBDIR/wrap/miq/core/seams_initializer.h
HEADERS += $$VCGLIBDIR/wrap/miq/core/stiffening.h

# VCG
HEADERS += $$VCGLIBDIR/vcg/complex/algorithms/parametrization/tangent_field_operators.h
HEADERS += $$VCGLIBDIR/vcg/complex/algorithms/parametrization/distortion.h
HEADERS += $$VCGLIBDIR/wrap/io_trimesh/import_field.h
HEADERS += $$VCGLIBDIR/wrap/io_trimesh/export_field.h
SOURCES += $$VCGLIBDIR/wrap/ply/plylib.cpp
SOURCES += quadrangulator.cpp

win32{
# Awful problem with windows..
  DEFINES += NOMINMAX
}

mac{
  CONFIG   += console
# Mac specific Config required to avoid to make application bundles
  CONFIG   -= app_bundle
  LIBS += -L/opt/local/lib -lamd -lcamd -lccolamd -lcholmod -lcolamd -lcxsparse -lblas -framework accelerate
}

#Comiso
mac{
LIBS += -L $$COMISOBUILDDIR/Build/lib/CoMISo/ -lCoMISo
  QMAKE_POST_LINK +="install_name_tool -change @executable_path/../lib/CoMISo/libCoMISo.dylib @executable_path/libCoMISo.dylib $$TARGET ; "
  QMAKE_POST_LINK +="cp -P $$COMISOBUILDDIR/Build/lib/CoMISo/libCoMISo.dylib . ; "
}

