#ifdef _WIN32
	#include "headers\\classes.hpp"
#elif _WIN64
	#include "headers\\classes.hpp"
#else
	#include "headers/classes.hpp"
#endif

Gene::Gene()
{
	generatePosition();
}

void Gene::generatePosition()
{
	random_device rd;
	mt19937 gen{rd()};
	uniform_int_distribution<int> distribution(0,1);
	value = distribution(gen);
}

int Gene::getValue() const
{
	return value;
}

void Gene::mutatePosition()
{
	if(value==0)
		value = 1;
	else
		value = 0;

}