#pragma once
#include "qdebug.h"
#include <stdio.h>
#include <stdlib.h>
#include <array>
#include "complexNum.h"
using namespace std;

#define FRAMES_PER_BUFFER (1024)
#define NUM_CHANNELS (2)
#define NUM_DELAYED_SAMPLES (0)

class channelProcesser{
public:
	channelProcesser();
	~channelProcesser();

	//float** computeData(float* buffer);
	void computeData(vector<float> channelOne, vector<float> channelTwo, std::vector<std::vector<float>>& result);

protected:


private:
	int first;
	int N;
	float lambda = 0.8;

	vector<float> istftSaved;
	vector<float> weigths;
	vector<float> delayedSamples;
	vector<float> savedChannelOneBuffer;
	vector<float> savedChannelTwoBuffer;
	vector<float> PSDreference;
	vector<float> PSDerror;
	vector<float> U;
	vector<complexNum> CPSDrefErr;


	vector<complexNum> impulseResponseFD;
	
	void hamming(vector<float>& window);
	void fourier(vector<float>& data, int N, int isign);
	void stft(vector<float> actualData, vector<float> storedData, int N, vector<vector <complexNum>>& output);
	void istft(vector<vector <complexNum>> data, int N, vector<float>& outputData);

	void makeCpx(float* data, int N, vector<complexNum>& output);
	void makeReal(vector<complexNum> data, int N, vector<float>& output);
	void normalize(vector<complexNum>& data);
	void channelProcesser::makeStereo(vector<float> data, vector<float>& stereo);

	vector<complexNum> LMSFD(vector<complexNum> mixture, vector<complexNum> reference, int N);
};