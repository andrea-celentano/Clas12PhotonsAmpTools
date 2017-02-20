#include <iostream>

#include "Clas12PhotonsAmplitudeEventGenerator.h"
#include "Clas12PhotonsPSEventGenerator.h"

#include "IUAmpTools/ConfigurationInfo.h"
#include "IUAmpTools/ConfigFileParser.h"
#include "IUAmpTools/AmpToolsInterface.h"

#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "TH1D.h"
#include "TCanvas.h"

Clas12PhotonsAmplitudeEventGenerator::Clas12PhotonsAmplitudeEventGenerator(const string &cfgfile, int Nevents) :
		m_dbPDG(0), m_Ebeam(11.), m_seed(0), m_Np(0), m_reaction(0), m_GenerationDone(false), m_Nevents(Nevents), m_ATI(0), m_hTweight(0) {
	//init the DB
	m_dbPDG = TDatabasePDG::Instance();
	if (m_dbPDG == 0) {
		cerr << "Can't instantiate new TDatabasePDG " << endl;
		cerr << "Exit " << endl;
		exit(1);
	}

	//Setup the PS event generator
	ConfigFileParser parser(cfgfile);
	ConfigurationInfo* cfgInfo = parser.getConfigurationInfo();
	cfgInfo->display();

	m_reaction = cfgInfo->reactionList()[0];
	m_PSgenerator = new Clas12PhotonsPSEventGenerator();
	m_PSgenerator->setReaction(m_reaction);

	//Setup the AmpTools interface
	m_ATI = new AmpToolsInterface(cfgInfo);

	m_doTweight = false;
	m_Nt = 100E3;

	m_EfficiencyDone = false;
	m_efficiency = 0;

	m_savedMin = 100;
	m_safetyFactor = 2;

	m_Np = m_PSgenerator->getNp();

	gRandom->SetSeed(m_seed);
}

Clas12PhotonsAmplitudeEventGenerator::~Clas12PhotonsAmplitudeEventGenerator() {
	if (m_ATI) delete m_ATI;
	if (m_PSgenerator) delete m_PSgenerator;
	if (m_hTweight) delete m_hTweight;
}

void Clas12PhotonsAmplitudeEventGenerator::setEbeam(double ebeam) {
	double s;
	double M;

	M = m_dbPDG->GetParticle("proton")->Mass();
	s = M * M + 2 * M * m_Ebeam;
	m_PSgenerator->setEbeam(ebeam);

	if (m_doTweight) {
		if (m_hTweight != 0) delete m_hTweight;
		m_hTweight = new TH1D("m_hTweight", "m_hTweight", 1000, 0, s);
	}

	m_EfficiencyDone = false;
	m_GenerationDone = false;
}

void Clas12PhotonsAmplitudeEventGenerator::EnableTweight() {
	double s;
	double M;

	M = m_dbPDG->GetParticle("proton")->Mass();
	s = M * M + 2 * M * m_Ebeam;
	if (m_hTweight != 0) delete m_hTweight;
	m_hTweight = new TH1D("m_hTweight", "m_hTweight", 1000, 0, s);

	m_doTweight = true;
}

void Clas12PhotonsAmplitudeEventGenerator::DisableTweight() {
	m_doTweight = false;
}

void Clas12PhotonsAmplitudeEventGenerator::GenerateEvents(int Nevents) {
	this->m_Nevents = Nevents;
	this->GenerateEvents();
}

void Clas12PhotonsAmplitudeEventGenerator::GenerateEvents() {
	int N_PS;
	int it_generation;
	int saved;

	double t, wt;
	double intensity, maxIntensity;

	if (!m_EfficiencyDone) this->computeEfficiency();

	N_PS = int(m_Nevents * m_safetyFactor / m_efficiency);
	it_generation = 1;
	cout << " Will start generating " << N_PS << " PS events" << endl;
	m_kinVPS.clear();
	m_kinVGenerated.clear();




	while (1) {
		Info("GenerateEvents", "Generation iteration %i : generate PS events", it_generation);
		m_ATI->clearEvents();
		//Load previous events - if any
		for (int it = 1; it < it_generation; it++) {
			for (int i = 0; i < N_PS; i++) {
				m_ATI->loadEvent(&m_kinVPS[(it - 1) * N_PS + i], (it - 1) * N_PS + i, N_PS * it_generation);
			}
		}
		for (int i = 0; i < N_PS; i++) {
			m_PSgenerator->Generate();
			if (m_doTweight) {
				t = -(m_PSgenerator->GetAllParticlesAmpToolsOrder()[2] - m_PSgenerator->GetAllParticlesAmpToolsOrder()[3]).M2();
				wt = m_hTweight->GetBinContent(m_hTweight->FindBin(t));
				if (wt < gRandom->Uniform(0, m_wtMax)) {
					i--;
					continue;
				}
			}
			if (i % (N_PS / 10) == 0) cout << " PS event: " << i << endl;
			Kinematics m_kin(m_PSgenerator->GetAllParticlesAmpToolsOrder());
			m_kinVPS.push_back(m_kin); //save the event also in this vector, for later use if this loop is entered again
			m_ATI->loadEvent(&m_kin, i, N_PS * it_generation);
		}
		Info("GenerateEvents", "Generation iteration %i : computing intensity ", it_generation);
		maxIntensity = m_ATI->processEvents(m_reaction->reactionName());
		Info("GenerateEvents", "Done. Max Intensity is: %f", maxIntensity);
		Info("GenerateEvents", "Doing hit-or-miss");

		saved = 0;
		m_kinVGenerated.clear();

		for (int i = 0; i < N_PS; i++) {
			if ((i % 10000) == 0) cout << "Event " << i << " intensity: " << m_ATI->intensity(i) << endl;
			intensity = m_ATI->intensity(i);
			if (intensity > gRandom->Uniform(0, maxIntensity)) {  //if intensity is bigger than random number between 0 and max, keep it.
				saved++;
				if (m_doTweight) {
					t = -(m_ATI->kinematics(i)->particleList()[2] - m_ATI->kinematics(i)->particleList()[3]).M2(); //A.C. definitively need to do this better
					wt = m_hTweight->GetBinContent(m_hTweight->FindBin(t));
				}
				else{
					wt=1;
				}
				Kinematics *m_kin = m_ATI->kinematics(i);
				m_kin->setWeight(wt);
				m_kinVGenerated.push_back(*m_kin);
			}
		}
		if (saved >= m_Nevents) {
			Info("GenerateEvents", "Enough events were generated in iteration: %i", it_generation);
			m_GenerationDone = true;
			break;
		} else {
			Info("GenerateEvents", "NOT enough events were generated in iteration: %i, only %i out of %i. Repeat", it_generation, saved, m_Nevents);
			it_generation++;
		}

	}
}

double Clas12PhotonsAmplitudeEventGenerator::GetEfficiency() {
	if (m_EfficiencyDone == false) {
		Info("GetEfficiency", "computeEfficiency was not called yet. Doing so now!");
		this->computeEfficiency();
	}
	return m_efficiency;
}

void Clas12PhotonsAmplitudeEventGenerator::computeEfficiency() {
	int it_efficiency = 0;
	int saved;
	int N_PS;
	double intensity, maxIntensity;

	double t, wt;

	if (m_doTweight) {
		Info("computeEfficiency", "doing pre-calculation with t-weight from amplitude");
		m_ATI->clearEvents();
		for (int i = 0; i < m_Nt; i++) {
			m_PSgenerator->Generate();
			Kinematics m_kin(m_PSgenerator->GetAllParticlesAmpToolsOrder());
			m_ATI->loadEvent(&m_kin, i, m_Nt);
		}
		maxIntensity = m_ATI->processEvents(m_reaction->reactionName());
		for (int i = 0; i < m_Nt; i++) {
			t = -(m_ATI->kinematics(i)->particleList()[2] - m_ATI->kinematics(i)->particleList()[3]).M2(); //A.C. definitively need to do this better
			intensity = m_ATI->intensity(i);
			m_hTweight->Fill(t, intensity);
		}
		m_hTweight->Scale(1. / m_Nt);
		m_wtMax = m_hTweight->GetMaximum();
		Info("computeEfficiency", "done");
	}

	while (1) {
		Info("computeEfficiency", "Start efficiency computation iteration %i", it_efficiency);
		m_ATI->clearEvents();
		N_PS = m_Nevents * pow(10, it_efficiency);
		Info("computeEfficiency", "Generating %i PS events", N_PS);
		for (int i = 0; i < N_PS; i++) {
			m_PSgenerator->Generate();
			if (m_doTweight) {
				t = -(m_PSgenerator->GetAllParticlesAmpToolsOrder()[2] - m_PSgenerator->GetAllParticlesAmpToolsOrder()[3]).M2(); //A.C. definitively need to do this better
				wt = m_hTweight->GetBinContent(m_hTweight->FindBin(t));
				if (wt < gRandom->Uniform(0, m_wtMax)) {
					i--;
					continue;
				}

			}
			if (i % (N_PS / 10) == 0) cout << " PS event: " << i << endl;
			Kinematics m_kin(m_PSgenerator->GetAllParticlesAmpToolsOrder());
			m_ATI->loadEvent(&m_kin, i, N_PS);
		}


		Info("computeEfficiency", " Efficiency iteration %i, PS events generated. Computing intensity", it_efficiency);
		maxIntensity = m_ATI->processEvents(m_reaction->reactionName());
		Info("computeEfficiency", " Intensity computation done. Max intensity is %f", maxIntensity);
		Info("computeEfficiency", "doing accept/reject...");

		saved = 0;

		for (int i = 0; i < N_PS; i++) {
			if ((i % (N_PS / 10)) == 0) cout << "Event " << i << " intensity: " << m_ATI->intensity(i) << endl;
			intensity = m_ATI->intensity(i);
			if (intensity > gRandom->Uniform(0, maxIntensity)) {  //if intensity if bigger than random number between 0 and max, keep it.
				saved++;
			}
		}

		Info("computeEfficiency", "Done: obtained events are: %i ", saved);
		if (saved > m_savedMin) {
			m_efficiency = 1. * saved / N_PS;
			Info("computeEfficiency", "Efficiency was computed: %f", m_efficiency);
			break;
		} else {
			Info("computeEfficiency", "Too few events. Efficiency was not computed. Repeat!");
			it_efficiency++;
		}
	}
	m_EfficiencyDone = true;
}

double Clas12PhotonsAmplitudeEventGenerator::GetWeight(int evt){
	if (m_GenerationDone == false) {
			Info("GetDecay", "Need to generate events first. Doing so now");
			this->GenerateEvents();
		}
		if (evt >= m_Nevents) {
			Error("GetDecay", "Request for evt %i, there are only %i events", evt, m_Nevents);
		}
		return m_kinVGenerated[evt].weight();

}

TLorentzVector Clas12PhotonsAmplitudeEventGenerator::GetDecay(int evt, int ip) {
	if (m_GenerationDone == false) {
		Info("GetDecay", "Need to generate events first. Doing so now");
		this->GenerateEvents();
	}
	if (evt >= m_Nevents) {
		Error("GetDecay", "Request for evt %i, there are only %i events", evt, m_Nevents);
	}

	return m_kinVGenerated[evt].particleList()[ip + 2];

}

vector<TLorentzVector> Clas12PhotonsAmplitudeEventGenerator::GetFinalStateParticles(int evt) {
	vector<TLorentzVector> v;

	v.push_back(m_kinVGenerated[evt].particleList()[1]); //scattered e'

	for (int ip = 0; ip < (m_Np-1); ip++) {
		v.push_back(m_kinVGenerated[evt].particleList()[ip + 3]); //all the others
	}
	return v;
}

vector<TLorentzVector> Clas12PhotonsAmplitudeEventGenerator::GetAllParticlesAmpToolsOrder(int evt) {
	vector<TLorentzVector> v;

	for (int ip = 0; ip < m_Np + 2; ip++) { //plus2 because of initial state
		v.push_back(m_kinVGenerated[evt].particleList()[ip]);
	}
	return v;

}

