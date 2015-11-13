#ifndef FERMIONTESTER_HPP
#define FERMIONTESTER_HPP

#include "SpinorTester.hpp"
#include "GaugefieldTester.hpp"

//todo: is there a better solution to this?
struct FermionTestParameters : public NonEvenOddSpinorTestParameters, GaugefieldTestParameters
{
	FermionTestParameters(const ReferenceValues referenceValuesIn, const LatticeExtents latticeExtentsIn, const SpinorFillTypes spinorFillTypesIn, const GaugefieldFillType gaugefieldFillTypesIn) :
		NonEvenOddSpinorTestParameters(referenceValuesIn, latticeExtentsIn, spinorFillTypesIn),
		GaugefieldTestParameters(referenceValuesIn, latticeExtentsIn, gaugefieldFillTypesIn) {};
	FermionTestParameters() : NonEvenOddSpinorTestParameters(), GaugefieldTestParameters() {};
};

class FermionTester : public SpinorTester
{
public:
	FermionTester(std::string kernelName, std::string inputfileIn, int numberOfValues = 1):
		SpinorTester(kernelName, getSpecificInputfile(inputfileIn), numberOfValues), gaugefieldBuffer(nullptr)
	{
			code = SpinorTester::device->getFermionCode();
			physics::lattices::GaugefieldParametersImplementation params(&SpinorTester::system->get_inputparameters());
			gaugefield = new physics::lattices::Gaugefield(*SpinorTester::system, &params, *prng);
	}
	FermionTester(std::string kernelName, std::vector<std::string> parameterStrings, int numberOfValues = 1):
	SpinorTester(kernelName, parameterStrings, numberOfValues), gaugefieldBuffer(nullptr)
	{
			code = SpinorTester::device->getFermionCode();
			physics::lattices::GaugefieldParametersImplementation params(&SpinorTester::system->get_inputparameters());
			gaugefield = new physics::lattices::Gaugefield(*SpinorTester::system, &params, *prng);
	}
	FermionTester(std::string kernelName, const ParameterCollection & pC, const FermionTestParameters & testParameters):
		SpinorTester(kernelName, pC, testParameters), gaugefield(nullptr)
	{
		gaugefieldBuffer = new hardware::buffers::SU3( calculateGaugefieldSize(testParameters.SpinorTestParameters::latticeExtents), device);
		const Matrixsu3 * gf_host = createGaugefield(calculateGaugefieldSize(testParameters.SpinorTestParameters::latticeExtents), testParameters.GaugefieldTestParameters::fillType);
		device->getGaugefieldCode()->importGaugefield(gaugefieldBuffer, gf_host);
		delete[] gf_host;

		code = SpinorTester::device->getFermionCode();
	}

	FermionTester(std::string kernelName, const ParameterCollection & pC, const SpinorTestParameters & testParameters):
		SpinorTester(kernelName, pC, testParameters), gaugefield(nullptr), gaugefieldBuffer(nullptr)
	{
		code = SpinorTester::device->getFermionCode();
	}

	~FermionTester()
	{
		if (gaugefield)
		{
			delete gaugefield;
		}
	}
	
protected:
	const hardware::code::Fermions * code;
	const hardware::buffers::SU3 * gaugefieldBuffer;
	physics::lattices::Gaugefield * gaugefield; //todo: remove
	
	//todo: remove both fcts.
	std::string getSpecificInputfile(std::string inputfileIn)
	{
		//todo: this is ugly, find a better solution.
		// The problem is that the parent class calls a similar fct.
		return "../fermions/" + inputfileIn;
	}
	
	const hardware::buffers::SU3* getGaugefieldBuffer() {
		return gaugefield->get_buffers()[0];
	}
};

#endif
