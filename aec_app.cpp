#include "aec_app.h"
#include "portaudio.h"
#include <qdebug.h>
#include "audioProcesser.h"

#define FRAMES_PER_BUFFER (1024)

AEC_app::AEC_app(QWidget *parent)
	: QMainWindow(parent)
{
	ui.setupUi(this);
	PaError err = Pa_Initialize();
	if (err != paNoError)
	{
		qDebug() << "Could not initialize PortAudio.";
		exit(2);
	}
	int deviceNum = Pa_GetDeviceCount();

	const PaDeviceInfo *deviceInfo, **AsioDeviceInfo;
	AsioDeviceInfo = new const PaDeviceInfo *[deviceNum];
	int * AsioDeviceIndeces = new int[deviceNum];
	int AsioDeviceNum = 0;
	int k = 0;

	//searching for asio devices and adding them to combo box
	for (int i = 0; i < deviceNum; i++)
	{
		deviceInfo = Pa_GetDeviceInfo(i);
		if (deviceInfo->hostApi == Pa_HostApiTypeIdToHostApiIndex(paASIO))
		{
			if (deviceInfo->maxInputChannels > 0)
			{
				ui.comboBox_input->addItem(QString("[ASIO] %1").arg(deviceInfo->name));
				ui.comboBox_input->setItemData(k, QVariant(i));
				//qDebug() << i;
			}
			if (deviceInfo->maxOutputChannels > 0)
			{
				ui.comboBox_output->addItem(QString("[ASIO] %1").arg(deviceInfo->name));
				ui.comboBox_output->setItemData(k, QVariant(i));
				//qDebug() << i;
			}
			k++;
			connect(ui.comboBox_input, SIGNAL(activated(int)), this, SLOT(showChannels()));
			connect(ui.comboBox_output, SIGNAL(activated(int)), this, SLOT(showChannels()));

			AsioDeviceIndeces[AsioDeviceNum] = i;
			AsioDeviceInfo[AsioDeviceNum++] = deviceInfo;
		}
	}

	/*ui.comboBox_method->addItem(QString("Record/Play"));
	ui.comboBox_method->setItemData(0, QVariant(0));
	ui.comboBox_method->addItem(QString("LMS Frequncy"));
	ui.comboBox_method->setItemData(1, QVariant(1));*/

	delete[] AsioDeviceIndeces;
	delete[] AsioDeviceInfo;

	//signals for buttons
	QPushButton* buttonStart = ui.pushButton_start;
	QPushButton* buttonStop = ui.pushButton_stop;
	connect(buttonStart, SIGNAL(clicked()), this, SLOT(startWorking()));
	connect(buttonStop, SIGNAL(clicked()), this, SLOT(stopWorking()));
	ui.statusBar->showMessage(tr("Ready"));
	buttonStop->setDisabled(true);

	ui.lineEdit_numIn->setDisabled(true);
	ui.lineEdit_numOut->setDisabled(true);

	ui.graphPlot->addGraph();
	ui.graphPlot->addGraph();
	ui.graphPlot->addGraph();
	ui.graphPlot->xAxis->setRange(0, FRAMES_PER_BUFFER);
	ui.graphPlot->yAxis->setRange(-0.001, 0.001);
	ui.graphPlot->replot();

	ui.graphPlot_2->addGraph();
	ui.graphPlot_2->xAxis->setRange(0, FRAMES_PER_BUFFER);
	ui.graphPlot_2->yAxis->setRange(-0.005, 0.005);
	ui.graphPlot_2->replot();

	/*ui.sliderMu->setDisabled(true);
	ui.sliderResponse->setDisabled(true);

	ui.doubleSpinBox_alpha->setDisabled(true);
	ui.doubleSpinBox_mu->setDisabled(true);*/
	showChannels();
}

AEC_app::~AEC_app()
{

}

//slot for showing basic info about selected devices
void AEC_app::showChannels() 
{
	const PaDeviceInfo *deviceInfo;
	int deviceIn = ui.comboBox_input->itemData(ui.comboBox_input->currentIndex()).toInt();
	int deviceOut = ui.comboBox_output->itemData(ui.comboBox_output->currentIndex()).toInt();

	QString s;
	qDebug() << "in: " << deviceIn << "out" << deviceOut;

	deviceInfo = Pa_GetDeviceInfo(deviceIn);
	s = QString("%1").arg(deviceInfo->maxInputChannels);
	ui.lineEdit_numIn->setDisabled(true);
	ui.lineEdit_numIn->setText(s);
	deviceInfo = Pa_GetDeviceInfo(deviceOut);
	s = QString("%1").arg(deviceInfo->maxOutputChannels);
	ui.lineEdit_numOut->setDisabled(true);
	ui.lineEdit_numOut->setText(s);

	ui.statusBar->showMessage(QString("Ready <in: %1 out: %2>").arg(deviceIn).arg(deviceOut));
}

//this runs audio processing
void AEC_app::startWorking(){ 
	int inDevice = 0;
	int outDevice = 0;
	int method = 0;

	inDevice = ui.comboBox_input->itemData(ui.comboBox_input->currentIndex()).toInt();
	outDevice = ui.comboBox_output->itemData(ui.comboBox_output->currentIndex()).toInt();
	ui.pushButton_start->setDisabled(true);
	ui.pushButton_stop->setEnabled(true);

	_settings = new Settings();
	_settings->inDevice = inDevice;
	_settings->outDevice = outDevice;
	_settings->usedInChannels = 2;
	_settings->usedOutChannels = 2;
	_settings->method = method;

	qDebug() << "used method" << _settings->method;
	qDebug() << "From main thread: " << QThread::currentThreadId();
	ui.statusBar->showMessage(QString("Working..."));
	audio = new audioProcesser();
	connect(audio, SIGNAL(showChannelData(vector<float>, vector<float>, vector<float>)), this, SLOT(dataToshow(vector<float>, vector<float>, vector<float>)));

	audio->open(_settings);
	audio->start();
}

//this stops audio processing
void AEC_app::stopWorking() 
{
	ui.statusBar->showMessage(QString("Ready"));
	ui.pushButton_stop->setDisabled(true);
	ui.pushButton_start->setEnabled(true);

	audio->stop();
	audio->wait();
	delete audio;
}