#ifndef CLAS12PHOTONSAMPLITUDE
#define CLAS12PHOTONSAMPLITUDE

#include "IUAmpTools/Amplitude.h"
#include "IUAmpTools/UserAmplitude.h"
#include "IUAmpTools/AmpParameter.h"
#include "GPUManager/GPUCustomTypes.h"
#include "TLorentzVector.h"

#include <utility>
#include <string>
#include <complex>
#include <vector>

#ifdef GPU_ACCELERATION

void GPUConstant_exec( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO,
		GDouble module, GDouble phase);

#endif // GPU_ACCELERATION

using std::complex;
using namespace std;

//Unfortunately, I need to use a structure to move the ElectronScatteringTerm from the calcElectronScattering to the calcAmplitude, because all members that I derive a const!

typedef struct {

	complex<GDouble> JP;
	complex<GDouble> J0;
	complex<GDouble> JM;

} ElectronScatteringTerm;

class Kinematics;

template<class T> class Clas12PhotonsAmplitude: public UserAmplitude<T> { //Inherits from UserAmplitude.

public:

	Clas12PhotonsAmplitude<T>() :
			UserAmplitude<T>() {

		One.real(1.);
		One.imag(0.);
		I.real(0.);
		I.imag(1.);

	}
	Clas12PhotonsAmplitude<T>(const vector<string>& args);

	virtual ~Clas12PhotonsAmplitude<T>() {
	}

	string name() const {
		return "Clas12PhotonsAmplitude";
	}

	complex<GDouble> calcAmplitude(GDouble** pKin) const;

	int calcElectronScattering(GDouble** pKin, ElectronScatteringTerm &ElectronScattering) const;
	virtual complex<GDouble> calcHelicityAmplitude(int helicity, GDouble** pKin) const = 0; //this will be derived by the user in his amplitude!!!

private:
	complex<GDouble> One;
	complex<GDouble> I;

#ifdef GPU_ACCELERATION

	void launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const;

	bool isGPUEnabled() const {return true;}

#endif // GPU_ACCELERATION

private:

protected:
	//it's important that they're protected, NOT private.
	int m_helicity_beam;     //beam helicity
	int m_helicity_electron; //scattered electron helicity

};

#include "Clas12PhotonsAmplitude.tpp"

#endif

