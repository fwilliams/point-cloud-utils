
#include <QApplication>
#include "window.h"

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);
    Window wdw;
    wdw.show();
    return app.exec();
}
