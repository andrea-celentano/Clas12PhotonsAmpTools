#ifndef CLAS12PHOTONSAMPLITUDEEVENTGENERATOR
#define CLAS12PHOTONSAMPLITUDEEVENTGENERATOR

#include <vector>
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TGenPhaseSpace.h"
#include "TRandom3.h"

#include "IUAmpTools/Kinematics.h"

using namespace std;

class TH1D;
class ReactionInfo;
class TDatabasePDG;
class Clas12PhotonsPSEventGenerator;
class AmpToolsInterface;

class Clas12PhotonsAmplitudeEventGenerator {

public:

	Clas12PhotonsAmplitudeEventGenerator(const string &cfgfile, int Nevents);
	~Clas12PhotonsAmplitudeEventGenerator();

	void setSeed(double seed) {
		m_seed = seed;
		gRandom->SetSeed(m_seed);
	}

	void setEbeam(double ebeam);
	void GenerateEvents();
	void GenerateEvents(int Nevents);


	void setEfficiencySaverdMin(int n) {
		m_savedMin = n;
	}
	void computeEfficiency();
	double GetEfficiency();

	inline void SetNt(int nt){m_Nt=nt;}
	void EnableTweight();
	void DisableTweight();

	void setSafectyFactor(int f) {
		m_safetyFactor = f;
	}

	 TLorentzVector GetDecay(int evt,int ip);
	 double GetWeight(int evt);
	 vector<TLorentzVector> GetAllParticlesAmpToolsOrder(int evt);
	 vector<TLorentzVector> GetFinalStateParticles(int evt);


private:

	//helper DB
	TDatabasePDG *m_dbPDG;

	//reaction
	ReactionInfo* m_reaction;

	//AmpTools interface
	AmpToolsInterface *m_ATI;
	//PS event generator
	Clas12PhotonsPSEventGenerator *m_PSgenerator;

	//How many events to generate
	int m_Nevents;

	//seed
	double m_seed;

	//efficiency of the computation
	double m_efficiency;
	bool m_EfficiencyDone;
	int m_savedMin;

	//beam energy
	double m_Ebeam;

	int m_safetyFactor;
	bool m_GenerationDone;

	bool m_doTweight;
	int m_Nt;
	double m_wtMax;
	TH1D *m_hTweight;

	//particles in the final state
	int m_Np;
	vector<TLorentzVector> m_vP;
	vector<Kinematics> m_kinVPS;
	vector<Kinematics> m_kinVGenerated;
	vector<double> m_intensities;

};

#endif
