/*
 * JPsiDataWriteLUND.hh
 *
 *  Created on: Nov 22, 2016
 *      Author: celentan
 */

#ifndef CLAS12PHOTONS_DATAWRITERLUND_H_
#define CLAS12PHOTONS_DATAWRITERLUND_H_


#include "IUAmpTools/Kinematics.h"

#include "TLorentzVector.h"
#include "TVector3.h"
#include <fstream>


class TDatabasePDG;
class TParticlePDG;

class Clas12PhotonsDataWriterLUND {
public:
	Clas12PhotonsDataWriterLUND(const string& outFile);
	virtual ~Clas12PhotonsDataWriterLUND();

	void writeEvent( const Kinematics& kin,const vector<TVector3>& vertex,int *pid=0,int *status=0);
	void writeEvent( const vector<TLorentzVector>& P,const vector<TVector3>& vertex,int *pid=0,int *status=0,double weight=1);
	int eventCounter() const { return m_eventCounter; }

private:

	  std::ofstream m_outFile;
	  int m_eventCounter;


	  /*Variables of current event*/
	  Kinematics m_kin;                //holds 4-vectors and weight
	  vector<TVector3> m_vertex; //holds vertexes;
	  int m_nP;                        //number of particles in current event
	  int *m_pid;                     //holds particles PID
	  int *m_status;                  //holds status
	  
	  /*strings and methods for these*/
	  string m_header;
	  string m_particle;


	  void header_line();
	  void particle_line(int np);
	  
	  /*helper DB*/
	  TDatabasePDG *m_PDGdb;
	  TParticlePDG *m_PDGparticle;

};

#endif /* JPSIIO_JPSIDATAWRITERLUND_H_ */
