#include "generationExecutable.h"

generationExecutable::generationExecutable(int argc, const char* argv[]) : generalExecutable(argc, argv)
{

}

void generationExecutable::setIterationParameters()
{
	thermalizationSteps = parameters.get_thermalizationsteps();
	writeFrequency = parameters.get_writefrequency();
	saveFrequency = parameters.get_savefrequency();
}
