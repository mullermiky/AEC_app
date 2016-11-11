#pragma once
#include <portaudio.h>
#include "pa_asio.h"
#include <qthread.h>
#include "channelProcesser.h"
#include <stdio.h>

struct Settings
{
	int inDevice;
	int outDevice;
	int usedInChannels;
	int usedOutChannels;
	int method;
};

class audioProcesser : public QThread{
	Q_OBJECT
public:
	audioProcesser();
	~audioProcesser();

	bool open(Settings* settings);

	public slots:
	void stop();

signals:
	void showChannelData(vector<float> channelsBuffer, vector<float> response, vector<float> computedResult);

private:
	static int mujCallback(const void*, void*, unsigned long, const PaStreamCallbackTimeInfo*, PaStreamCallbackFlags, void*);
	void run();

	int first;
	Settings* mySettings;
	QMutex a_mutex;
	bool a_stop;
	channelProcesser *channelProc;
	float* _buf;
	vector<float> response;
	vector<float> computedResult;
	vector<float> channelsBuffer;

	float* impulseResp;
	float* computed;
	float parameterMU = 0.00059;
	float parameterALPHA = 0.004;

	PaStream * stream;
	PaStreamParameters inputParameters;
	PaStreamParameters outputParameters;
};