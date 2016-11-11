#ifndef AEC_APP_H
#define AEC_APP_H

#include <QtWidgets/QMainWindow>
#include "ui_aec_app.h"
#include "qdebug.h"
#include <stdio.h>

#include "audioProcesser.h"

#define FRAMES_PER_BUFFER (1024)

class AEC_app : public QMainWindow
{
	Q_OBJECT

public:
	AEC_app(QWidget *parent = 0);
	~AEC_app();

private:
	Ui::AEC_appClass ui;
	audioProcesser* audio;
	Settings* _settings;

public slots:
	void dataToshow(vector<float> channels, vector<float> response, vector<float> computedResult)
	{

		//qDebug() << "emited";
		QVector<double> channelOneVec(FRAMES_PER_BUFFER), channelTwoVec(FRAMES_PER_BUFFER), y(FRAMES_PER_BUFFER),
			impR(FRAMES_PER_BUFFER), computedData(FRAMES_PER_BUFFER);


		/*if (!impulseResponse.empty() && !result.empty()){*/
			for (int i = 0; i < FRAMES_PER_BUFFER; i++)
			{
				y[i] = i;
				channelOneVec[i] = channels[2 * i];    //jeste mozna upravit, kazdy kanal extra vector
				channelTwoVec[i] = channels[2 * i + 1];
				impR[i] = response[i];
				computedData[i] = computedResult[2 * i];
			}
		//}


		/*ui.graphPlot->graph(0)->setPen(QPen(Qt::blue));
		ui.graphPlot->graph(0)->setData(y, channelOneVec);*/
		ui.graphPlot->graph(0)->setPen(QPen(Qt::red));
		ui.graphPlot->graph(0)->setData(y, channelTwoVec);
		ui.graphPlot->graph(1)->setPen(QPen(Qt::green));
		ui.graphPlot->graph(1)->setData(y, computedData);

		ui.graphPlot->replot();

		ui.graphPlot_2->graph(0)->setData(y, impR);
		ui.graphPlot_2->replot();
	};

private slots:
	void showChannels();
	void startWorking();
	void stopWorking();
	//void showSliderValues();

};

#endif // AEC_APP_H
