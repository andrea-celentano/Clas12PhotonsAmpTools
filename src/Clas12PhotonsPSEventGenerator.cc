#include <iostream>

#include "IUAmpTools/ConfigurationInfo.h"
#include "Clas12PhotonsPSEventGenerator.h"

#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "TH1D.h"

Clas12PhotonsPSEventGenerator::Clas12PhotonsPSEventGenerator() :
		m_dbPDG(0), m_Ebeam(11.), m_Wdistr(0), m_reaction(0), m_seed(0), m_Np(0), m_generatorMaxWt(0) {
	//init the DB
	m_dbPDG = TDatabasePDG::Instance();
	if (m_dbPDG == 0) {
		cerr << "Can't instantiate new TDatabasePDG " << endl;
		cerr << "Exit " << endl;
		exit(1);
	}

	//set initial state 4-momenta
	m_beam.SetXYZT(0, 0, m_Ebeam, m_Ebeam);
	m_target.SetXYZT(0, 0, 0, m_dbPDG->GetParticle(2212)->Mass()); //proton target
	m_P0 = m_beam + m_target;

	//W2=Mp * Mp + 2E0 * Mp- 2E' * Mp- 4EE' sin^2(theta'/2)
	//max value:Mp*Mp + 2E0 * Mp
	m_Wmax = sqrt(m_target.M2() + 2 * m_Ebeam * m_target.M());
	m_Wmin = m_target.M();

	//now set the FT nominal acceptance (can always be change with the proper method!)
	m_thetaMax = 4.5 * TMath::DegToRad();
	m_thetaMin = 2.5 * TMath::DegToRad();

	m_EprimeMin = 0.5;
	m_EprimeMax = 4.5;

}

void Clas12PhotonsPSEventGenerator::setReaction(ReactionInfo *reaction) {

	TParticlePDG *particle;
	double Qtot;

	m_reaction = reaction;

	m_Np = m_reaction->particleList().size() - 2 ; //first two entries are beam and target
	m_Wmin = 0;
	//add the final state e'
	m_pname.push_back("e-");
	m_pmass.push_back(m_dbPDG->GetParticle("e-")->Mass());
	m_pid.push_back(11);
	Qtot = -1;

	//use the TDatabasePDG to search for these particles and assign masses - pids
	//also do a consistency check for charges

	for (int ip = 3; ip < m_reaction->particleList().size(); ip++) { //starting from the fourth entry, the first 3 are beam, target, recoil
		particle = 0;
		particle = m_dbPDG->GetParticle(m_reaction->particleList()[ip].c_str());
		if (!particle) {
			cout << "Can't find particle: " << m_pname[ip] << endl;
			cout << "Exit " << endl;
			exit(1);
		} else {
			m_pname.push_back(m_reaction->particleList()[ip]);
			m_pmass.push_back(particle->Mass());
			m_pid.push_back(particle->PdgCode());
			m_Wmin += particle->Mass();
			Qtot += particle->Charge() / 3; //return in units |e| / 3
			cout << " Particle: " << particle->GetName() << " found in pdg DB " << endl;
			cout << " Mass: " << particle->Mass() << " Charge: " << particle->Charge() / 3 << endl;
		}
	}

	//The full reaction is: e p --> e' p X.
	//In the MesonEx context, this is g* p -> p X
	//The ReactionInfo part only carry information about "p X" final state, i.e. the hadronic part
	if (Qtot != 0) {
		cerr << "Error, the total charge of the hadronic part of final state is: " << Qtot << endl;
		cerr << "Should be 0. Exit" << endl;
		exit(1);
	}
}

void Clas12PhotonsPSEventGenerator::computeWdistr() {
	double Wt, Wval;
	TLorentzVector Pw;

	m_generator.SetDecay(m_P0, m_Np, &(m_pmass[0])); //c++ standard guarantees this!
	m_generatorMaxWt = m_generator.GetWtMax();
	if (m_Wdistr) delete m_Wdistr;
	m_Wdistr = new TH1D("Wdistr", "Wdistr", 1000, m_Wmin, m_Wmax);
	Info("computeWdistr", "Start computing the events for the W-distr sampling");
	for (int ievt = 0; ievt < 100000; ievt++) {
		Wt = m_generator.Generate();
		if (Wt > m_generatorMaxWt) m_generatorMaxWt = Wt; //should not happen
		if (Wt < gRandom->Uniform(0, m_generatorMaxWt)) {
			ievt--;
			continue;
		}
		Pw.SetXYZT(0, 0, 0, 0);
		for (int ip = 1; ip < m_Np; ip++) {
			Pw += *(m_generator.GetDecay(ip));
		}
		Wval = Pw.M();
		m_Wdistr->Fill(Wval);
	}
	Info("computeWdistr", "Done");
}

void Clas12PhotonsPSEventGenerator::Generate() {

	TLorentzVector Peprime, Pw;
	double u, u_min, u_max;
	double ctheta, ctheta_min, ctheta_max, phi;
	double M, E0;
	double Wval, WminGen, WmaxGen;
	double Eprime;
	double Wt;
	if (m_Wdistr == 0) {
		Info("Generate", "W distribution not yet sampled. Doing so now");
		this->computeWdistr();
	}
	M = m_target.M();
	E0 = m_Ebeam;

	m_vP.clear();

	//First part of the computation: pseudo 2-body reaction e p -> e (W), with W all the other particles in final state
	//See A. Celentano PhD thesis, p.112

	//1-A: extract theta value using inversion!
	//note the order!
	ctheta_min = cos(m_thetaMax);
	ctheta_max = cos(m_thetaMin);

	u_min = M / 2 * (ctheta_min + 1) / (M + E0 * (1 - ctheta_min));
	u_max = M / 2 * (ctheta_max + 1) / (M + E0 * (1 - ctheta_max));

	u = gRandom->Uniform(u_min, u_max);

	ctheta = (2 * u * (E0 + M) - M) / (M + 2 * u * E0);

	//1-B: extract W.
	//Since W*2 = M*M + 2*M*(E0-E')-2EE'(1-ctheta), and now theta is fixed, a lower (upper) limit on E' is an upper (lower) limit on W

	WmaxGen = sqrt(M * M + 2 * M * (E0 - m_EprimeMin) - 2 * E0 * m_EprimeMin * (1 - ctheta));
	WminGen = sqrt(M * M + 2 * M * (E0 - m_EprimeMax) - 2 * E0 * m_EprimeMax * (1 - ctheta));

	if (WmaxGen > m_Wmax) {
		Warning("Generate", "The physical max value of W is: %f, but from e' scattering and limits on the energy the limit is: %f ", m_Wmax, WmaxGen);
		WmaxGen = m_Wmax;
	}
	if (WminGen < m_Wmin) {
		Warning("Generate", "The physical min value of W is: %f, but from e' scattering and limits on the energy the limit is: %f ", m_Wmin, WminGen);
		WminGen = m_Wmin;
	}

	while (1) {
		Wval = m_Wdistr->GetRandom();
		if ((Wval > WminGen) && (Wval < WmaxGen)) break;
	}

	Eprime = (-Wval * Wval + M * M + 2 * M * E0) / (2 * M + 2 * E0 * (1 - ctheta));

	//1-C: fix the kinematics of scattered e' in the LAB frame
	phi = gRandom->Uniform(0, TMath::TwoPi());
	Peprime.SetXYZT(Eprime * sqrt(1 - ctheta * ctheta) * sin(phi), Eprime * sqrt(1 - ctheta * ctheta) * cos(phi), Eprime * ctheta, Eprime);
	m_vP.push_back(Peprime);
	//1-D: fix the kinematics of the pseudo-particle "W" in the LAB frame: this is simply P0-Peprime
	Pw=m_P0-Peprime;


	//2: now handle the "decay" process W->other final state particles
	m_generator.SetDecay(Pw, m_Np - 1, &(m_pmass[1]));
	m_generatorMaxWt = m_generator.GetWtMax()*2; //*2 is temporary

	while (1) {
		Wt = m_generator.Generate();
		if (Wt > gRandom->Uniform(0, m_generatorMaxWt)) {
			for (int ip = 0; ip < m_Np - 1; ip++) {
				m_vP.push_back(*(m_generator.GetDecay(ip)));
			}
			break;
		}
	}

}

