#pragma once

class complexNum{
private:
	float real;
	float imag;

public:
	complexNum(float r = 0.0, float im = 0.0);
	complexNum(const complexNum&);

	complexNum operator +(complexNum);
	complexNum operator +(float);
	complexNum operator -(complexNum);
	complexNum operator *(complexNum);
	complexNum operator *(float);
	complexNum operator /(complexNum);
	complexNum getConjugate();
	void setData(float, float);
	float getReal();
	float getImag();
	void setReal(float);
	void setImag(float);
	bool operator ==(complexNum);
	void operator =(complexNum);
	float absoluteValue();

	
};