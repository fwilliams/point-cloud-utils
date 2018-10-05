include(../common.pri)
TARGET = space_packer
SOURCES +=  space_packer.cpp \
            ../../../../vcglib/wrap/qt/Outline2ToQImage.cpp \
            ../../../../vcglib/wrap/qt/outline2_rasterizer.cpp
QT += opengl svg
