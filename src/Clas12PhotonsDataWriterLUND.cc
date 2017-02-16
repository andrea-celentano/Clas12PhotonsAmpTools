/*
 * JPsiDataWriteLUND.cc
 *
 *  Created on: Nov 22, 2016
 *      Author: celentan
 */

#include "Clas12PhotonsDataWriterLUND.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include <sstream>   


Clas12PhotonsDataWriterLUND::Clas12PhotonsDataWriterLUND(const string& outFile) {
	// TODO Auto-generated constructor stub
	m_outFile.open(outFile.c_str());
	m_eventCounter = 0;   
	m_PDGdb=TDatabasePDG::Instance();

	m_PDGparticle=0;

}

Clas12PhotonsDataWriterLUND::~Clas12PhotonsDataWriterLUND() {
	// TODO Auto-generated destructor stub
	m_outFile.close();
}

void Clas12PhotonsDataWriterLUND::writeEvent(const Kinematics& kin,const vector<TVector3>& vertex,int *pid,int *status){
	m_kin=kin;
	m_nP=m_kin.particleList().size();
	if (vertex.size()!=m_nP){
	  Error("writeEvent","vertex entries are:%i while 4-momenta entries are: %i",vertex.size(),m_nP);
	  return;
	}

	m_vertex=vertex;
	m_pid=pid;
	m_status=status;
	


	/*start writing output
	//header line has 10 entries. Only first one is meaningfull
	NumberParticles
	Number of target nucleons
	Number of target protons
	Target Polarization
	Beam Polarization
	x
	y
	W
	Q2
	nu
	//Then N=NumberParticles lines follows, each with 14 values. Only those marked with "X" are used
	index 
	charge
	type(1 is active)
	particleID PDG format (X)
	parent index
	daugther index
	momentum px(GeV) (X)
	momentum py(GeV) (X)
	momentum pz(GeV) (X)
	energy (GeV)
	mass (GeV)
	vertex x (cm)
	vertex y (cm)
	vertex z (cm)
	*/
	this->header_line();
	m_outFile<<m_header<<endl;
	for (int ip=0;ip<m_nP;ip++){
	  this->particle_line(ip);
	  m_outFile<<m_particle<<endl;
	}
}

void Clas12PhotonsDataWriterLUND::writeEvent(const vector<TLorentzVector>& P,const vector<TVector3>& vertex,int *pid,int *status,double weight){
	Kinematics kin(P,weight);
	this->writeEvent(kin,vertex,pid,status);
}

void Clas12PhotonsDataWriterLUND::header_line(){
  ostringstream stream;
  stream<<m_nP<<" "; //number of particles
  stream<<"0 0 0 0 0 0 0 0 "; //8 fields not used
  stream<<m_kin.weight(); //save the event weight in the last header field

  m_header=stream.str();
}

void Clas12PhotonsDataWriterLUND::particle_line(int np){
  ostringstream stream;
  int status;
  if (m_pid==0){
    Error("particle_line","no pid is provided - pointer of m_pid is 0!");
    return;
  }
  if (m_status==0) status=1;
  else status=m_status[np];


 
 
  m_PDGparticle=m_PDGdb->GetParticle(m_pid[np]);

  //index
  stream<<np<<" "; 
  //charge
  stream<<m_PDGparticle->Charge()/3<<" "; //Charge() returns in units of |e|/3
  //status
  stream<<status<<" ";
  //PID
  stream<<m_pid[np]<<" ";
  //parent and daughter index - not used
  stream<<"0 0 ";
  //momentum px py pz (GeV)
  stream<<m_kin.particle(np).Px()<<" "<<m_kin.particle(np).Py()<<" "<<m_kin.particle(np).Pz()<<" ";
  //energy
  stream<<m_kin.particle(np).E()<<" ";
  //mass
  stream<<m_PDGparticle->Mass()<<" ";
  //vertex
  stream<<m_vertex[np].X()<<" "<<m_vertex[np].Y()<<" "<<m_vertex[np].Z()<<endl;

  
  m_particle=stream.str();
  



}
