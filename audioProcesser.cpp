#include "audioProcesser.h"
#include "pa_asio.h"
#include "portaudio.h"
#include "channelProcesser.h"
#include <qdebug.h>
#include "qdebug.h"
#include <time.h>
//#include <vld.h>

#define SAMPLE_RATE (32000)
#define FRAMES_PER_BUFFER (1024)
#define PA_SAMPLE_TYPE  paFloat32
#define NUM_CHANNELS (2)

static void PrintSupportedStandardSampleRates(
	const PaStreamParameters *inputParameters,
	const PaStreamParameters *outputParameters)
{
	static double standardSampleRates[] = {
		8000.0, 9600.0, 11025.0, 12000.0, 16000.0, 22050.0, 24000.0, 32000.0,
		44100.0, 48000.0, 88200.0, 96000.0, 192000.0, -1 /* negative terminated  list */
	};
	int     i, printCount;
	PaError err;

	printCount = 0;
	for (i = 0; standardSampleRates[i] > 0; i++)
	{
		err = Pa_IsFormatSupported(inputParameters, outputParameters, standardSampleRates[i]);
		if (err == paFormatIsSupported)
		{
			qDebug() << standardSampleRates[i];
		}
	}
}

audioProcesser::audioProcesser()
{
	PaError err;
	err = Pa_Initialize();
	if (err != paNoError)
	{
		qDebug() << "Could not initialize PortAudio.";
		exit(2);
	}
	a_stop = false;
	first = 0;
}

audioProcesser::~audioProcesser()
{
	//delete channelProc;
}

void audioProcesser::run()
{

	if (stream)
	{
		qDebug() << "Pa_StartStream( stream ) ...";

		response.resize(FRAMES_PER_BUFFER);
		computedResult.resize(FRAMES_PER_BUFFER*2);
		channelsBuffer.resize(FRAMES_PER_BUFFER*2);
		PaError err = Pa_StartStream(stream);
		if (err != paNoError) {
			const PaHostErrorInfo* info;
			info = Pa_GetLastHostErrorInfo();
			Pa_Terminate();
			qDebug() << "chyba zabijim PA v run";
			qDebug() << "chyba: " << Pa_GetErrorText(err);
			qDebug() << "detaily: " << info->errorCode << "->" << info->errorText << "(" << info->hostApiType << ")";
		}
	}
	clock_t lastTime;
	lastTime = clock();

	while (1)
	{
		if (double(clock() - lastTime) / CLOCKS_PER_SEC > 0.05 && !response.empty() && !computedResult.empty())
		{
			//qDebug() << "time now";
			lastTime = clock();			
			emit showChannelData(channelsBuffer, response, computedResult);
			if (a_stop) break;
		}

		msleep(1);
		if (a_stop) break;
	}

	if (stream)
	{
		PaError err = Pa_StopStream(stream);
		if (err != paNoError) {
			qDebug() << "PortAudio_processer chyba: " << Pa_GetErrorText(err);
		}
		Pa_CloseStream(stream);
		Pa_Terminate();

	}

	delete[] mySettings;
}

void audioProcesser::stop()
{
	qDebug() << "Thread::stop called from main thread: " << currentThreadId();
	//QMutexLocker locker(&a_mutex);
	a_stop = true;
	delete channelProc;

}

//this will open audio stream
bool audioProcesser::open(Settings* _settings)
{
	PaError err = Pa_Initialize();
	if (err != paNoError)
	{
		qDebug() << "Could not initialize PortAudio.";
		exit(2);
	}

	mySettings = _settings;

	int inputDevice = mySettings->inDevice;
	int outputDevice = mySettings->outDevice;
	int inputChannels = mySettings->usedInChannels;
	int outputChannels = mySettings->usedOutChannels;

	const PaDeviceInfo *inDevInfo = Pa_GetDeviceInfo(inputDevice);
	const PaDeviceInfo *outDevInfo = Pa_GetDeviceInfo(outputDevice);

	inputParameters.device = inputDevice;
	inputParameters.channelCount = inputChannels;
	inputParameters.sampleFormat = PA_SAMPLE_TYPE;
	inputParameters.suggestedLatency = inDevInfo->defaultLowInputLatency;
	inputParameters.hostApiSpecificStreamInfo = NULL;

	outputParameters.device = outputDevice;
	outputParameters.channelCount = outputChannels;
	outputParameters.sampleFormat = PA_SAMPLE_TYPE;
	outputParameters.suggestedLatency = outDevInfo->defaultLowOutputLatency;
	outputParameters.hostApiSpecificStreamInfo = NULL;

	qDebug() << "suggestedLatency: " << inputParameters.suggestedLatency;
	qDebug() << "inputchannels: " << inputParameters.channelCount;
	qDebug() << "outputchannels: " << outputParameters.channelCount;
	qDebug() << "inputdevice: " << inputParameters.device;
	qDebug() << "outputdevice: " << outputParameters.device;
	qDebug() << "sample rate: " << SAMPLE_RATE;
	qDebug() << "frames per buffer: " << FRAMES_PER_BUFFER;
	qDebug() << "supported rate: ";
	//PrintSupportedStandardSampleRates(&inputParameters, &outputParameters);

	PaAsioStreamInfo asioOutputInfo, asioInputInfo;

	asioOutputInfo.size = sizeof(PaAsioStreamInfo);
	asioOutputInfo.hostApiType = paASIO;
	asioOutputInfo.version = 1;
	asioOutputInfo.flags = paAsioUseChannelSelectors;
	int outputChannelSelectors[2];
	outputChannelSelectors[0] = 0;
	outputChannelSelectors[1] = 1;
	asioOutputInfo.channelSelectors = outputChannelSelectors;

	asioInputInfo = asioOutputInfo;
	int inputChannelSelectors[2];
	inputChannelSelectors[0] = 0;
	inputChannelSelectors[1] = 1;
	asioInputInfo.channelSelectors = inputChannelSelectors;

	inputParameters.hostApiSpecificStreamInfo = &asioOutputInfo;
	outputParameters.hostApiSpecificStreamInfo = &asioInputInfo;

	channelProc = new channelProcesser();

	err = Pa_OpenStream(&stream, &inputParameters, &outputParameters, SAMPLE_RATE,
		FRAMES_PER_BUFFER, paClipOff, mujCallback, this);

	const PaHostErrorInfo* info;
	info = Pa_GetLastHostErrorInfo();


	if (err != paNoError) {
		Pa_Terminate();
		qDebug() << "chyba zabijim PA v open";
		qDebug() << "chyba: " << Pa_GetErrorText(err) << "(" << err << ")";
		qDebug() << "detaily: " << info->errorCode << "->" << info->errorText << "(" << info->hostApiType << ")";
		return false;
	}
	return true;
}

int audioProcesser::mujCallback(const void *inputBuffer, void *outputBuffer,
	unsigned long framesPerBuffer,
	const PaStreamCallbackTimeInfo* timeInfo,
	PaStreamCallbackFlags statusFlags,
	void *userData)
{
	float* out = (float *)outputBuffer;
	float* in = (float *)inputBuffer;
	float* buf = (float *)inputBuffer;

	audioProcesser *audioProc = (audioProcesser*)(userData);
	Settings *set = audioProc->mySettings;

	unsigned long i;
	int finished = 0;
	std::vector<float> channelOne(framesPerBuffer);
	std::vector<float> channelTwo(framesPerBuffer);
	std::vector<std::vector<float>> result(2, std::vector<float>(2*framesPerBuffer));

	//move data to vectors
	for (i = 0; i < framesPerBuffer; i++)
	{
		channelOne[i] = in[i*NUM_CHANNELS];
		channelTwo[i] = in[i*NUM_CHANNELS + 1];
		audioProc->channelsBuffer[i*NUM_CHANNELS] = in[i*NUM_CHANNELS];
		audioProc->channelsBuffer[i*NUM_CHANNELS + 1] = in[i*NUM_CHANNELS + 1];
	}

	/*qDebug() << "left:" << channelOne[512];
	qDebug() << "right:" << channelTwo[512];*/

	
	channelProcesser* _channelProc = audioProc->channelProc;

	_channelProc->computeData(channelOne, channelTwo, result);
	for (i = 0; i < framesPerBuffer; i++)
	{
		buf[i*NUM_CHANNELS] = result[0][i*NUM_CHANNELS];
		buf[i*NUM_CHANNELS + 1] = result[0][i*NUM_CHANNELS + 1];
	}
	audioProc->computedResult = result[0];
	audioProc->response = result[1];
	for (i = 0; i < 2 * framesPerBuffer; i++)
	{
		*out++ = *buf++; //prehrava zpracovany
		//*out++ = *in++; //prehrava origo
	}
	
	audioProc->first = 1;

	return finished;
}