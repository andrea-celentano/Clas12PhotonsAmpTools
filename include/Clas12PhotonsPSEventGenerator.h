#ifndef CLAS12PHOTONSPSEVENTGENERATOR
#define CLAS12PHOTONSPSEVENTGENERATOR

#include <vector>
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"
#include "TRandom3.h"
using namespace std;

class TH1D;
class ReactionInfo;
class TDatabasePDG;

class Clas12PhotonsPSEventGenerator {

public:

	Clas12PhotonsPSEventGenerator();

	int getNp() const{
		return m_Np;
	}

	double getThetaMax() const {
		return m_thetaMax;
	}

	void setThetaMax(double thetaMax) {
		m_thetaMax = thetaMax;
	}

	double getThetaMin() const {
		return m_thetaMin;
	}

	void setThetaMin(double thetaMin) {
		m_thetaMin = thetaMin;
	}

	double getEbeam() const {
		return m_Ebeam;
	}

	void setEbeam(double ebeam) {
		m_Ebeam = ebeam;
		m_beam.SetXYZT(0, 0, m_Ebeam, m_Ebeam);
		m_P0 = m_beam + m_target;

		m_Wmax = sqrt(m_target.M2() + 2 * m_Ebeam * m_target.M());
		if (m_Wdistr != 0) this->computeWdistr(); //W distr. was already calculated, need to do again
	}

	double getSeed() const {
		return m_seed;
	}

	void setSeed(double seed) {
		m_seed = seed;
		gRandom->SetSeed(m_seed);
	}

	const TH1D* getWdistr() const {
		return m_Wdistr;
	}

	double getEprimeMax() const {
		return m_EprimeMax;
	}

	void setEprimeMax(double eprimeMax) {
		m_EprimeMax = eprimeMax;
	}

	double getEprimeMin() const {
		return m_EprimeMin;
	}

	void setEprimeMin(double eprimeMin) {
		m_EprimeMin = eprimeMin;
	}

	void setReaction(ReactionInfo *reaction);

	void Generate();
	TLorentzVector GetDecay(int ip) {
		return m_vP[ip];
	}
	vector<TLorentzVector> GetFinalStateParticles(){
		return m_vP;
	}

	 /*Returns in the order required by AmpTools: beam,e',target,other particles*/
	vector<TLorentzVector> GetAllParticlesAmpToolsOrder(){
		vector<TLorentzVector> v;
		v.push_back(m_beam);
		v.push_back(m_vP[0]);
		v.push_back(m_target);
		for (int ip=1;ip<m_Np;ip++) v.push_back(m_vP[ip]);
		return v;
	}



private:

	void computeWdistr();

	//the reaction
	ReactionInfo *m_reaction;

	//helper DB
	TDatabasePDG *m_dbPDG;

	//seed
	double m_seed;

	//beam energy
	double m_Ebeam;

	//scattered e- angles and energy
	double m_thetaMin; //in radians!
	double m_thetaMax; //in radians!

	double m_EprimeMin;
	double m_EprimeMax;

	//W distribution
	TH1D* m_Wdistr;

	double m_Wmax; //the physical maximum value of W
	double m_Wmin; //the physical minimum value of W

	//4-vectors of particles
	TLorentzVector m_beam;
	TLorentzVector m_target;
	TLorentzVector m_P0;
	TLorentzVector m_recoilElectron;
	TLorentzVector m_virtualPhoton;

	//particles in the final state
	int m_Np;
	vector<string> m_pname;
	vector<int> m_pid;
	vector<double> m_pmass;
	vector<TLorentzVector> m_vP;

	TGenPhaseSpace m_generator;
	double m_generatorMaxWt;
};

#endif
