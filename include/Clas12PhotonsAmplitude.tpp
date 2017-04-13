/* This is a derived class from AmpTools UserAmplitude.
 It's written for the quasi-real photoproduction physics of MesonEx in CLAS12.
 Basically what it does is to perform the calculation of the leptonic current, event by event, summing over the helicity of the quasi-real photon.
 The spin of the initial and final state electron are two protected members of this class, that are initialized to their value by the derived amplitude from the user!
 The user has to derive its own amplitude from this class, specifing the photoproduction amplitude as a function of the helicity.
 */
 
 
 /*Important: the formulas that follows are fine when used in a frame where:
 
 1) The virtual photon is moving along +z
 2) The hadronic plane is xz
*/
template<class T> Clas12PhotonsAmplitude<T>::Clas12PhotonsAmplitude(const vector<string>& args) :
		UserAmplitude<T>(args) {
	assert(args.size() >= 2); //helicity beam, helicity scattered electron,- -- then others.
	m_helicity_beam = atoi(args[0].c_str());
	m_helicity_electron = atoi(args[1].c_str());

	One.real(1.);
	One.imag(0.);
	I.real(0.);
	I.imag(1.);

}

//the order of the particles is supposed to be:
//0: beam  
//1: scattered e-
//2: target
//3: recoil baryon
//4.. n: all the other particles

//This ordering is convenient for the helicity!!

template<class T> int Clas12PhotonsAmplitude<T>::calcElectronScattering(GDouble** pKin, ElectronScatteringTerm &ElectronScattering) const {

	const int Ibeam = 0;
	const int Iscattered = 1;
	const int Itarget = 2;
	const int Irecoil = 3;

	TLorentzVector beam;
	TLorentzVector electron;
	TLorentzVector gamma;
	TLorentzVector target;
	TLorentzVector recoil;


	double E1, E2, Eg;
	double theta1, theta2, thetag, phi1, phi2, phig;
	double Pg;
	double Mg;
	double Q2;


	double thetaRot1,thetaRot2,thetaRot3;

	//define here all the relevant variables
	beam.SetPxPyPzE(pKin[Ibeam][1], pKin[Ibeam][2], pKin[Ibeam][3], pKin[Ibeam][0]);
	electron.SetPxPyPzE(pKin[Iscattered][1], pKin[Iscattered][2], pKin[Iscattered][3], pKin[Iscattered][0]);
	target.SetPxPyPzE(pKin[Itarget][1], pKin[Itarget][2], pKin[Itarget][3], pKin[Itarget][0]);
	recoil.SetPxPyPzE(pKin[Irecoil][1], pKin[Irecoil][2], pKin[Irecoil][3], pKin[Irecoil][0]);

	gamma = beam - electron;

	//First, rotate along z to have Pgamma in the xy plane
	 thetaRot1=atan2(-gamma.Y(),gamma.X());
	 beam.RotateZ(thetaRot1);
	 electron.RotateZ(thetaRot1);
	 target.RotateZ(thetaRot1);
	 recoil.RotateZ(thetaRot1);
	 gamma = beam - electron;

	 //then, rotate along y to have Pgamma along z
	 thetaRot2=atan2(-gamma.X(),gamma.Z());
	 beam.RotateY(thetaRot2);
	 electron.RotateY(thetaRot2);
 	 target.RotateY(thetaRot2);
 	 recoil.RotateY(thetaRot2);
	 gamma = beam - electron;
	 

	 //then, rotate along z to have recoil (and hance target) in xz. This last rotation ensures the hadronic plane is in xz, and phi is the angle between hadronic and leptonic plane!
	 //Recoil IS necessary, since the target itself (0,0,0,Mproton) does not define any plane (in the lab frame)
	 thetaRot3=atan2(-recoil.Y(),recoil.X());
	 beam.RotateZ(thetaRot3);
	 electron.RotateZ(thetaRot3);
 	 target.RotateZ(thetaRot3);
 	 recoil.RotateZ(thetaRot3);
	 gamma = beam - electron;
	 


	
	//now I can use Vincent's formulas
	E1 = beam.E();
	E2 = electron.E();
	Eg = E1 - E2;
	Q2 = -(beam - electron).M2();

	theta1 = beam.Theta();
	theta2 = electron.Theta();
	phi1 = beam.Phi();
	phi2 = electron.Phi();

	thetag = gamma.Theta();
	phig = gamma.Phi();

	Pg = gamma.Z(); //quasi real photon momentum (Purely along z)
	Mg = sqrt(Pg * Pg - Eg * Eg);

	if ((m_helicity_beam == 1) && (m_helicity_electron == 1)) {

		ElectronScattering.JP = 2. * sqrt(2. * E1 * E2) * cos(theta1 / 2) * sin(theta2 / 2) * exp(-I * phi2);
		ElectronScattering.JM = -2. * sqrt(2. * E1 * E2) * cos(theta2 / 2) * sin(theta1 / 2) * exp(I * phi1);

		ElectronScattering.J0 = 0.;

	} else if ((m_helicity_beam == 1) && (m_helicity_electron == -1)) {
		ElectronScattering.JP = 0;
		ElectronScattering.JM = 0;
		ElectronScattering.J0 = 0;
	} else if ((m_helicity_beam == -1) && (m_helicity_electron == 1)) {
		ElectronScattering.JP = 0;
		ElectronScattering.JM = 0;
		ElectronScattering.J0 = 0;
	}

	else if ((m_helicity_beam == -1) && (m_helicity_electron == -1)) {

		ElectronScattering.JP = 2. * sqrt(2. * E1 * E2) * cos(theta2 / 2) * sin(theta1 / 2) * exp(-I * phi1);
		ElectronScattering.JM = -2. * sqrt(2. * E1 * E2) * cos(theta1 / 2) * sin(theta2 / 2) * exp(I * phi2);
		ElectronScattering.J0 = 0.;
	}
	ElectronScattering.JP /= Q2;
	ElectronScattering.J0 /= Q2;
	ElectronScattering.JM /= Q2;



return 0;
}

template<class T> complex<GDouble> Clas12PhotonsAmplitude<T>::calcAmplitude(GDouble** pKin) const {

	complex < GDouble > helP;
	complex < GDouble > hel0;
	complex < GDouble > helM;
	complex < GDouble > amp;

	//calculate the photoproduction part of the amplitude trough the code provided by the user
	helP = calcHelicityAmplitude(1, pKin);
	hel0 = calcHelicityAmplitude(0, pKin);
	helM = calcHelicityAmplitude(-1, pKin);

	//trigger the calculation of the electron scattering part.
	//This fills the JP,J0,JM complex numbers
	ElectronScatteringTerm ElectronScattering;
	calcElectronScattering(pKin, ElectronScattering);

	//Note that the 1/Q2 factor is already included in JP,J0,JM
	amp = helP * ElectronScattering.JP + hel0 * ElectronScattering.J0 + helM * ElectronScattering.JM;

/*	cout<<"HEL:: "<<helP<<" "<< ElectronScattering.JP <<" "<<hel0<<" "<< ElectronScattering.J0 << " "<<helM <<" "<<ElectronScattering.JM<<endl;
	cout<<"HEL:: "<<amp<<" "<< m_helicity_beam<<" "<< m_helicity_electron<<" "<<endl;
	cout<<"HEL2:: "<<std::norm(amp)<<endl;
	cin.get();*/
	return amp;

}

