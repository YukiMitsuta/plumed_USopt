/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2021 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "Bias.h"
#include "ActionRegister.h"

namespace PLMD {
namespace bias {

//+PLUMEDOC BIAS RESTRAINTMATRIX
/*
Adds harmonic and/or linear restraints on one or more variables.

Either or both
of SLOPE and KAPPA must be present to specify the linear and harmonic force constants
respectively.  The resulting potential is given by:
\f[
  \transpose{x-a} \frac{K}{2} (x-a)
  = \sum_ij \frac{K_ij}{2} (x_i-a_i) (x_j-a_j)
\f].

The number of components for any vector of force constants must be equal to the number
of arguments to the action.

\par Examples

The following input tells plumed to restrain the two dihedral angles of alanine dipeptide (phi, psi).
\plumedfile
TORSION LABEL=phi ATOMS=5,7,9,15  #CV0
TORSION LABEL=psi ATOMS=7,9,15,17  #CV1
PRINT STRIDE=10 ARG=phi,psi FILE=COLVAR
RESTRAINTMATRIX ...
LABEL=Rest-MAT
ARG=phi,psi
AT=-2.313075,2.334774
KAPPA0=50.107449,-0.077964
KAPPA1=-0.077964,50.182735
... RESTRAINTMATRIX
\endplumedfile

*/
//+ENDPLUMEDOC

class RestraintMatrix : public Bias {
  std::vector<double> at;
  std::vector<std::vector<double> > kappa;
  std::vector<double> slope;
  Value* valueForce2;
public:
  explicit RestraintMatrix(const ActionOptions&);
  void calculate() override;
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(RestraintMatrix,"RESTRAINTMATRIX")

void RestraintMatrix::registerKeywords(Keywords& keys) {
  Bias::registerKeywords(keys);
  keys.use("ARG");
  keys.add("compulsory","SLOPE","0.0","specifies that the restraint is linear and what the values of the force constants on each of the variables are");
  keys.add("numbered","KAPPA","KAPPA\\f$x\\f$ is equal to the value of the force constants at initial.");
  keys.reset_style("KAPPA","compulsory");
  keys.add("compulsory","AT","the position of the restraint");
  keys.addOutputComponent("force2","default","the instantaneous value of the squared force due to this bias potential");
}

RestraintMatrix::RestraintMatrix(const ActionOptions&ao):
  PLUMED_BIAS_INIT(ao),
  at(getNumberOfArguments()),
  //kappa(getNumberOfArguments(),0.0),
  slope(getNumberOfArguments(),0.0)
{
  std::vector<double> kk( getNumberOfArguments() );
  parseVector("SLOPE",slope);
  for(int i=0;; i++) {
    // Try to read kappa
    if( !parseNumberedVector("KAPPA",i,kk) ) break;
    kappa.push_back(kk);
  }
  parseVector("AT",at);
  checkRead();

  log.printf("  at");
  for(unsigned i=0; i<at.size(); i++) log.printf(" %f",at[i]);
  log.printf("\n");
  log.printf("  with harmonic force constant");
  //for(unsigned i=0; i<kappa.size(); i++) log.printf(" %f",kappa[i]);
  for(unsigned i=0; i<kappa[0].size(); i++) {
    for(unsigned j=0; j<kappa[0].size(); j++) {
      log.printf(" %f",kappa[i][j]);
    }
    log.printf("\n");
  }
  log.printf("\n");
  log.printf("  and linear force constant");
  for(unsigned i=0; i<slope.size(); i++) log.printf(" %f",slope[i]);
  log.printf("\n");

  addComponent("force2");
  componentIsNotPeriodic("force2");
  valueForce2=getPntrToComponent("force2");
}


void RestraintMatrix::calculate() {
  double ene=0.0;
  double totf2=0.0;
  double f=0.0;
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    const double cv_i=difference(i,at[i],getArgument(i));
    f=0.0;
    const double m=slope[i];
    f -= m;
    for(unsigned j=0; j<getNumberOfArguments(); ++j) {
      const double cv_j=difference(j,at[j],getArgument(j));
      const double k=kappa[i][j];
      if (i == j) {
        f-=k*cv_i;
      } else {
        f-=0.5*(k*cv_j);
        const double kji=kappa[j][i];
        f-=0.5*kji*cv_j;
      }
      ene+=0.5*k*cv_i*cv_j;
    }
    ene+= m*cv_i;
    setOutputForce(i,f);
    totf2+=f*f;
  }
  setBias(ene);
  valueForce2->set(totf2);
}

}


}
