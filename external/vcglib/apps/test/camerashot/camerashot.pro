TARGET = camerashot_test
INCLUDEPATH += . ../../..
CONFIG += console stl
TEMPLATE = app
SOURCES += camerashot_test.cpp

win32: CONFIG += NOMINMAX

# Mac specific Config required to avoid to make application bundles
CONFIG -= app_bundle
