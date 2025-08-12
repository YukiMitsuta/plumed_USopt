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
  \sum_i \frac{k_i}{2} (x_i-a_i)^2 + m_i*(x_i-a_i)
\f].

The number of components for any vector of force constants must be equal to the number
of arguments to the action.

Additional material and examples can be also found in the tutorial \ref lugano-2

\par Examples

The following input tells plumed to restrain the distance between atoms 3 and 5
and the distance between atoms 2 and 4, at different equilibrium
values, and to print the energy of the restraint
\plumedfile
DISTANCE ATOMS=3,5 LABEL=d1
DISTANCE ATOMS=2,4 LABEL=d2
RESTRAINTMATRIX ARG=d1,d2 AT=1.0,1.5 KAPPA=150.0,150.0 LABEL=restraint
PRINT ARG=restraint.bias
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
/*
       if (i == j) {
         f-=k*cv_j;
       } else {
         f-=0.5*(k*cv_j);
      }
*/
      f-=k*cv_j;
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
