#include "complexNum.h"
#include "qmath.h";

complexNum::complexNum(float r, float im)
{
	real = r;
	imag = im;
}

complexNum::complexNum(const complexNum &c)
{
	real = c.real;
	imag = c.imag;
}

void complexNum::operator =(complexNum c)
{
	real = c.real;
	imag = c.imag;
}

complexNum complexNum::operator +(complexNum c)
{
	complexNum result;
	result.real = this->real + c.real;
	result.imag = this->imag + c.imag;
	return result;
}

complexNum complexNum::operator +(float f)
{
	complexNum result;
	result.real = this->real + f;
	result.imag = this->imag;
	return result;
}

complexNum complexNum::operator -(complexNum c)
{
	complexNum result;
	result.real = this->real - c.real;
	result.imag = this->imag - c.imag;
	return result;
}

complexNum complexNum::operator *(complexNum c)
{
	complexNum result;
	result.real = (this->real*c.real) - (this->imag*c.imag);
	result.imag = (this->real*c.imag) + (this->imag*c.real);
	return result;
}

complexNum complexNum::operator *(float f)
{
	complexNum result;
	result.real = this->real*f;
	result.imag = this->imag*f;
	return result;
}

complexNum complexNum::operator /(complexNum c)
{
	float div = (c.real*c.real) + (c.imag*c.imag);
	complexNum result;
	result.real = (this->real*c.real) + (this->imag*c.imag);
	result.real /= div;
	result.imag = (this->imag*c.real) - (this->real*c.imag);
	result.imag /= div;
	return result;
}

complexNum complexNum::getConjugate()
{
	complexNum result;
	result.real = this->real;
	result.imag = this->imag*-1;
	return result;
}

void complexNum::setData(float r, float im)
{
	this->real = r;
	this->imag = im;
}

void complexNum::setReal(float r)
{
	real = r;
}

void complexNum::setImag(float im)
{
	imag = im;
}

float complexNum::getReal()
{
	return real;
}

float complexNum::getImag()
{
	return imag;
}

bool complexNum::operator==(complexNum c)
{
	return (this->real == c.real) && (this->imag == c.imag) ? 1 : 0;
}

float complexNum::absoluteValue()
{
	return sqrtf(powf(this->real, 2) + powf(this->imag, 2));
}