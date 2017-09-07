#include "testUtilities.hpp"

#include "../../host_functionality/logger.hpp"

void printKernelInformation(std::string name)
{
  logger.info() << "Test kernel\t\"" << name << "\"\tagainst reference value";
}

double sumOfIntegers(const int start, const int end, const int increment) noexcept
{
	// One could also implement some variant of Faulhabers Formula here to save the loop
	double sum = 0.;
	for(int iteration = start; iteration <= end; iteration += increment)
	{
		sum += iteration;
	}
	return sum;
}

double sumOfIntegersSquared(const int end) noexcept
{
	return (2*end*end*end + 3*end*end + end) / 6.; // Faulhaber`s formula
}


