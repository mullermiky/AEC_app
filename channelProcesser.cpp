#include "channelProcesser.h"
#include <qdebug.h>
#include <qmath.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include "complexNum.h"
#include <qmath.h>
#include <complex>

using namespace std;

#define _USE_MATH_DEFINES
#define PI M_PI

#define TWOPI (2.0*PI)
#define NUM_DELAYED_SAMPLES (0)

channelProcesser::channelProcesser()
{
	//qDebug() << "channel";
	impulseResponseFD.resize(FRAMES_PER_BUFFER);
	delayedSamples.resize(NUM_DELAYED_SAMPLES);
	savedChannelOneBuffer.resize(FRAMES_PER_BUFFER, 0);
	savedChannelTwoBuffer.resize(FRAMES_PER_BUFFER, 0);
	weigths.resize(FRAMES_PER_BUFFER);
	istftSaved.resize(FRAMES_PER_BUFFER);

	PSDreference.resize(FRAMES_PER_BUFFER);
	PSDerror.resize(FRAMES_PER_BUFFER);
	CPSDrefErr.resize(FRAMES_PER_BUFFER);
	U.resize(FRAMES_PER_BUFFER, 0.00199);

	first = 0;
}

channelProcesser::~channelProcesser()
{
	
	qDebug() << "konec";
}

void channelProcesser::computeData(vector<float> _channelOne, vector<float> _channelTwo, std::vector<std::vector<float>>& out)
{
	vector<float> channelOne(_channelOne);
	vector<float> channelTwo(_channelTwo);
	vector<float> savedChannelOne(savedChannelOneBuffer);
	vector<float> savedChannelTwo(savedChannelTwoBuffer);
	vector<float> delayedChannel(FRAMES_PER_BUFFER);
	vector<float> result(FRAMES_PER_BUFFER);
	vector<float> stereoResult(2*FRAMES_PER_BUFFER);
	vector<float> impulseResponseTEMP(2 * FRAMES_PER_BUFFER + 1);
	vector<float> impulseResponseTD(FRAMES_PER_BUFFER);

	
	vector<complexNum> impulseResponseTEMPcpx(FRAMES_PER_BUFFER);
	vector<vector<complexNum>> cpxChannelOne(2, std::vector<complexNum>(FRAMES_PER_BUFFER));
	vector<vector<complexNum>> cpxChannelTwo(2, std::vector<complexNum>(FRAMES_PER_BUFFER));
	vector<vector<complexNum>> cpxResult(2, std::vector<complexNum>(FRAMES_PER_BUFFER));

	for (int i = 0; i < FRAMES_PER_BUFFER; ++i)
	{
		if (i <= NUM_DELAYED_SAMPLES - 1)
		{
			delayedChannel[i] = delayedSamples[i];
		}
		else
		{
			delayedChannel[i] = channelOne[i - NUM_DELAYED_SAMPLES];
		}
	}

	savedChannelOneBuffer = delayedChannel;
	savedChannelTwoBuffer = channelTwo;

	
	//---for testing fourier function---
	/*vector<float> vec(17, 1);
	for (int i = 0; i < vec.size(); i++)
		vec[i] = i;
	vector<complexNum> vecc(8);
	fourier(vec, 8, 1);
	fourier(vec, 8, -1);
	for (int i = 0; i < vecc.size(); i++)
		vecc[i].setData(vec[2 * i + 1], vec[2 * i + 2]);
	normalize(vecc);*/
	//----------------------------------

	stft(delayedChannel, savedChannelOne, FRAMES_PER_BUFFER, cpxChannelOne);
	stft(channelTwo, savedChannelTwo, FRAMES_PER_BUFFER, cpxChannelTwo);

	//istft(cpxChannelOne, FRAMES_PER_BUFFER, result);	//for testing stft/istft only

	if (first == 2)
	{
		cpxResult[0] = LMSFD(cpxChannelTwo[0], cpxChannelOne[0], FRAMES_PER_BUFFER);
		cpxResult[1] = LMSFD(cpxChannelTwo[1], cpxChannelOne[1], FRAMES_PER_BUFFER);
		istft(cpxResult, FRAMES_PER_BUFFER, result);
	}
	else if (first == 1)
	{
		cpxResult[0] = LMSFD(cpxChannelTwo[0], cpxChannelOne[0], FRAMES_PER_BUFFER);
		istft(cpxResult, FRAMES_PER_BUFFER, result);
	}
	else
	{
		istft(cpxChannelTwo, FRAMES_PER_BUFFER, result);
	}
	
	
	makeStereo(result, stereoResult);
	out[0] = stereoResult;

	for (int i = 0; i < impulseResponseFD.size(); i++)
	{
		impulseResponseTEMP[2 * i + 1] = impulseResponseFD[i].getReal();
		impulseResponseTEMP[2 * i + 2] = impulseResponseFD[i].getImag();
	}

	fourier(impulseResponseTEMP, FRAMES_PER_BUFFER, -1);
	for (int i = 0; i < impulseResponseTEMPcpx.size(); i++)
		impulseResponseTEMPcpx[i].setData(impulseResponseTEMP[2 * i + 1], impulseResponseTEMP[2 * i + 2]);

	normalize(impulseResponseTEMPcpx);
	makeReal(impulseResponseTEMPcpx, FRAMES_PER_BUFFER, impulseResponseTD);
	out[1] = impulseResponseTD;
	//first = 1;
	(first == 1) ? first = 2 : 0;
	(first == 0) ? first = 1 : 0;
	
}

//---------NEW METHODs----------

void channelProcesser::fourier(vector<float>& data, int N, int isign)
{
	int n, mmax, m, j, istep, i;
	float wtemp, wr, wpr, wpi, wi, theta;
	float tempr, tempi;

	n = N << 1;
	j = 1;
	for (i = 1; i < n; i += 2) {
		if (j > i) {
			tempr = data[j];     data[j] = data[i];     data[i] = tempr;
			tempr = data[j + 1]; data[j + 1] = data[i + 1]; data[i + 1] = tempr;
		}
		m = n >> 1;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	mmax = 2;
	while (n > mmax) {
		istep = 2 * mmax;
		theta = TWOPI / (isign*mmax);
		wtemp = sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi = sin(theta);
		wr = 1.0;
		wi = 0.0;
		for (m = 1; m < mmax; m += 2) {
			for (i = m; i <= n; i += istep) {
				j = i + mmax;
				tempr = wr*data[j] - wi*data[j + 1];
				tempi = wr*data[j + 1] + wi*data[j];
				data[j] = data[i] - tempr;
				data[j + 1] = data[i + 1] - tempi;
				data[i] += tempr;
				data[i + 1] += tempi;
			}
			wr = (wtemp = wr)*wpr - wi*wpi + wr;
			wi = wi*wpr + wtemp*wpi + wi;
		}
		mmax = istep;
	}
}

void channelProcesser::stft(vector<float> actualData, vector<float> storedData, int N, vector<vector <complexNum>>& CpxOut)
{

	int windowLength = N;

	vector<float> window(windowLength);
	vector<vector<float>> data(2, vector<float>(2 * N + 1));

	hamming(window);

	if (this->first == 0)
	{
		for (int i = 0; i < windowLength; i++)
		{
			CpxOut[0][i].setReal(window[i] * actualData[i]);
		}
	}
	else
	{
		for (int i = 0; i < windowLength; i++)
		{
			if (i < windowLength / 2)
			{
				CpxOut[0][i].setReal(window[i] * storedData[i + N / 2]);
			}
			else
			{
				CpxOut[0][i].setReal(window[i] * actualData[i - N / 2]);
			}
			CpxOut[1][i].setReal(window[i] * actualData[i]);
		}
	}
	for (int i = 0; i < CpxOut[0].size(); i++)
	{
		data[0][2 * i + 1] = CpxOut[0][i].getReal();
		data[0][2 * i + 2] = CpxOut[0][i].getImag();
		data[1][2 * i + 1] = CpxOut[1][i].getReal();
		data[1][2 * i + 2] = CpxOut[1][i].getImag();
	}

	fourier(data[0], N, 1);
	fourier(data[1], N, 1);

	for (int i = 0; i < CpxOut[0].size(); i++)
	{
		CpxOut[0][i].setReal(data[0][2 * i + 1]);
		CpxOut[0][i].setImag(data[0][2 * i + 2]);
		CpxOut[1][i].setReal(data[1][2 * i + 1]);
		CpxOut[1][i].setImag(data[1][2 * i + 2]);
	}
}

void channelProcesser::istft(vector<vector <complexNum>> CpxData, int N, vector<float>& flData)

{
	
	vector<float> realData_1(N);
	vector<float> realData_2(N);
	vector<vector<float>> data(2, vector<float>(2 * N + 1));
	for (int i = 0; i < CpxData[0].size(); i++)
	{
		data[0][2 * i + 1] = CpxData[0][i].getReal();
		data[0][2 * i + 2] = CpxData[0][i].getImag();
		data[1][2 * i + 1] = CpxData[1][i].getReal();
		data[1][2 * i + 2] = CpxData[1][i].getImag();
	}

	if (first == 0)
	{
		fourier(data[0], N, -1); //prvni cely blok
		for (int i = 0; i < CpxData[0].size(); i++)
		{
			CpxData[0][i].setReal(data[0][2 * i + 1]);
			CpxData[0][i].setImag(data[0][2 * i + 2]);
		}
		makeReal(CpxData[0], N, istftSaved);
		makeReal(CpxData[0], N, realData_2);
	}
	else
	{
		fourier(data[0], N, -1); //hybrid 1/2 z prvniho a 1/2 z druhyho
		fourier(data[1], N, -1); //dalsi blok
		for (int i = 0; i < CpxData[0].size(); i++)
		{
			CpxData[0][i].setReal(data[0][2 * i + 1]);
			CpxData[0][i].setImag(data[0][2 * i + 2]);
			CpxData[1][i].setReal(data[1][2 * i + 1]);
			CpxData[1][i].setImag(data[1][2 * i + 2]);
		}
		normalize(CpxData[0]);
		normalize(CpxData[1]);
		makeReal(CpxData[0], N, realData_1);
		makeReal(CpxData[1], N, realData_2);
		for (int i = 0; i < N / 2; i++)
		{
			istftSaved[i + FRAMES_PER_BUFFER / 2] = (istftSaved[i + FRAMES_PER_BUFFER / 2] + realData_1[i]) / 2;
			realData_2[i] = (realData_2[i] + realData_1[i + FRAMES_PER_BUFFER / 2]) / 2;
		}

	}
	flData = istftSaved;
	istftSaved = realData_2;
}

void channelProcesser::makeReal(vector<complexNum> data, int N, vector<float>& output)
{
	for (int i = 0; i < N; i++)
	{
		output[i] = data[i].getReal();
	}
}

void channelProcesser::makeCpx(float* data, int N, vector<complexNum>& output)
{
	vector<complexNum> result(N);

	for (int i = 0; i < N; i++)
	{
		output[i].setData(data[i], 0.0);
	}
}

void channelProcesser::normalize(vector<complexNum>& data)
{
	N = data.size()-1/2;
	for (int i = 0; i < N; i++)
	{
		data[i].setData(data[i].getReal() / N, data[i].getImag() / N);
	}
}

void channelProcesser::hamming(vector<float>& window)
{
	for (int i = 0; i < window.size(); i++)
	{
		window[i] = 0.53836 - 0.46164 * cos(2.0 * M_PI * i / float(window.size() - 1));
	}
}

vector<complexNum> channelProcesser::LMSFD(vector<complexNum> mixture, vector<complexNum> reference, int N)
{
	vector<complexNum> E(N);
	vector<complexNum> H(N);
	vector<complexNum> Y(N);
	vector<float> Htemp(2 * N + 1);
	vector<complexNum> cpxHtemp(N);
	vector<float> realH(N);
	vector<complexNum> X;
	vector<complexNum> R;
	vector<float> delenec(N);
	vector<float> delitel(N);
	vector<float> vydeleno(N);
	vector<float> rozdil1(N);
	vector<complex<double>> odmocnina(N);
	vector<float> rozdil2(N);
	vector<float> odmocnina2(N);
	vector<float> nasobek(N);

	X = mixture;
	R = reference;


	float c;
	H = impulseResponseFD;
	float mu = 0.00199;
	float beta = 0.1;

	for (int k = 0; k < N; k++)
	{
		Y[k] = (R[k] * H[k]);
		E[k] = X[k] - Y[k];
		H[k] = H[k] + ((R[k].getConjugate() * E[k]) * U[k]*beta);

		PSDreference[k] = (lambda * PSDreference[k]) + ((1 - lambda) * (powf((R[k].absoluteValue()), 2)));
		PSDerror[k]		= (lambda * PSDerror[k])	 + ((1 - lambda) * (powf((E[k].absoluteValue()), 2)));

		CPSDrefErr[k]	= ((R[k].getConjugate() * E[k]) * (1 - lambda)) + (lambda * PSDerror[k]);

		//TODO ---predelat zpátky
		delenec[k] = powf(CPSDrefErr[k].absoluteValue(), 2);
		delitel[k] = (PSDreference[k] * PSDerror[k]) + 0.000001;
		vydeleno[k] = delenec[k] / delitel[k];
		rozdil1[k] = 1 - vydeleno[k];
		odmocnina[k] = std::sqrt(std::complex<double>(rozdil1[k]));
		(odmocnina[k].real() == 0) ? odmocnina2[k] = odmocnina[k].imag() : odmocnina2[k] = odmocnina[k].real();
		rozdil2[k] = 1 - odmocnina2[k];
		nasobek[k] = 2 / PSDreference[k] + 0.000001;

		U[k] = nasobek[k] * rozdil2[k];
		//U[k] = (2 / PSDreference[k] + 0.001) * (1 - sqrt(1 - ((powf(CPSDrefErr[k].absoluteValue(), 2)) / ((PSDreference[k] * PSDerror[k]) + 0.001))));
	}
	//-----prevod H do casove oblasti----
	for (int i = 0; i < H.size(); i++)
	{
		Htemp[2 * i + 1] = H[i].getReal();
		Htemp[2 * i + 2] = H[i].getImag();
	}
	fourier(Htemp, N, -1);
	for (int i = 0; i < H.size(); i++)
	{
		cpxHtemp[i].setReal(Htemp[2 * i + 1]);
		cpxHtemp[i].setImag(Htemp[2 * i + 2]);
	}
	normalize(cpxHtemp);
	makeReal(cpxHtemp, N, realH);
	//-----------------------------------
	float max = 0;
	for (int i = 0; i < realH.size(); i++)
		(max < realH[i]) ? max = realH[i] : max = max;

	(max < 0.0009) ? impulseResponseFD = H : H;
	//-----------------------------------------
	
	//impulseResponseFD = H;
	
	return E;
}

void channelProcesser::makeStereo(vector<float> data, vector<float>& stereo)
{
	for (int i = 0; i < data.size(); i++)
	{
		stereo[2 * i] = data[i];
		stereo[2 * i + 1] = data[i];
	}
}