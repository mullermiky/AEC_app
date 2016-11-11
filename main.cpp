#include "aec_app.h"
#include <QtWidgets/QApplication>

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);
	AEC_app w;

	w.setFixedHeight(776);
	w.setFixedWidth(771);
	w.show();
	qRegisterMetaType<vector<float> >("vector<float>");
	return a.exec();
}
