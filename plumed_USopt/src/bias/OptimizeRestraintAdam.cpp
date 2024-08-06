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

//+PLUMEDOC BIAS OPTIMIZERESTRAINTADAM
/*
Add a optimization method of harmonic restraint on one or more variables.

The harmonic restraint on your system is given by:

\f[
V(\vec{s},t) = \frac{1}{2} \kappa(t) ( \vec{s} - \vec{s}_0(t) )^2
\f]

The time dependence of \f$\kappa\f$ and \f$\vec{s}_0\f$ are specified by a list of
STEP, KAPPA and AT keywords.  These keywords tell plumed what values \f$\kappa\f$ and \f$\vec{s}_0\f$
should have at the time specified by the corresponding STEP keyword.  In between these times
the values of \f$\kappa\f$ and \f$\vec{s}_0\f$ are linearly interpolated.

Additional material and examples can be also found in the tutorial \ref belfast-5

\par Examples

The following input is dragging the two dihedral angles of alanine dipeptide (phi, psi).
\plumedfile
TORSION LABEL=phi ATOMS=5,7,9,15  #CV0
TORSION LABEL=psi ATOMS=7,9,15,17  #CV1
PRINT STRIDE=10 ARG=phi,psi FILE=COLVAR_npt
OPTIMIZERESTRAINTADAM ...
OFFDIAGONAL
LABEL=OptRest
ARG=phi,psi
AT=-2.3,2.4
TARGET=-2.3,2.4
KAPPA0=50.0,0.0
KAPPA1=0.0,50.0
OPTSTRIDE=5000
CVSTEPSIZE=0.01
KAPPASTEPSIZE=0.01
OPTMETHOD=ADABELIEF
BETA1=0.9
BETA2=0.999
EPSILON=1.0e-8
... OPTIMIZERESTRAINTADAM
PRINT ...
STRIDE=500
FILE=COLVAR_kappa
ARG=OptRest.phi_phi_kappa,OptRest.phi_psi_kappa,OptRest.psi_phi_kappa,OptRest.psi_psi_kappa
... PRINT
PRINT ...
STRIDE=500
FILE=COLVAR_cntr
ARG=OptRest.phi_cntr,OptRest.psi_cntr
... PRINT
PRINT ...
STRIDE=500
FILE=COLVAR_mean
ARG=OptRest.phi_mean,OptRest.psi_mean
... PRINT
\endplumedfile

\attention Work is not computed properly when KAPPA is time dependent.

*/
//+ENDPLUMEDOC


class OptimizeRestraintAdam : public Bias {
  std::vector<std::vector<double> > kappa;

  std::vector<double> at;
  std::vector<double> target;
  //std::vector<long int> step;
  std::vector<double> oldaa;
  std::vector<double> oldk;
  std::vector<double> olddpotdk;
  std::vector<double> oldf;
  //std::vector<std::string> verse;
  std::vector<double> work;
  std::vector<double> stepsize_cv;
  std::vector<double> stepsize_kappa;
  std::vector<double> mean;
  std::vector<long int> stride;
  std::vector<std::string> optmethod;

  std::vector<double> m_cv;
  std::vector<std::vector<double>> m_kappa;
  double v_cv;
  double v_kappa;
  //std::vector<double> grad_cv;
  //std::vector<std::vector<double>> grad_kappa;
  std::vector<double> beta1;
  std::vector<double> beta2;
  std::vector<double> epsilon;

  unsigned N;
  int t;
  bool offdiagonal;

  double tot_work;
public:
  explicit OptimizeRestraintAdam(const ActionOptions&);
  void calculate() override;
  static void registerKeywords( Keywords& keys );
};

PLUMED_REGISTER_ACTION(OptimizeRestraintAdam,"OPTIMIZERESTRAINTADAM")

void OptimizeRestraintAdam::registerKeywords( Keywords& keys ) {
  Bias::registerKeywords(keys);
  keys.use("ARG");
  /*
  keys.add("compulsory","VERSE","B","Tells plumed whether the restraint is only acting for CV larger (U) or smaller (L) than "
           "the restraint or whether it is acting on both sides (B)");
  keys.add("numbered","STEP","This keyword appears multiple times as STEP\\f$x\\f$ with x=0,1,2,...,n. Each value given represents "
           "the MD step at which the restraint parameters take the values KAPPA\\f$x\\f$ and AT\\f$x\\f$.");
  keys.reset_style("STEP","compulsory");
  */
  keys.add("compulsory","OPTMETHOD","ADAM"," The method of optimization (default = Adam)");
  keys.add("compulsory","BETA1"," The parameter of beta1");
  keys.add("compulsory","BETA2"," The parameter of beta2");
  keys.add("compulsory","EPSILON"," The parameter of epsilon");

  keys.add("compulsory","AT","AT\\f$x\\f$ is equal to the initial position of the restraint.");
  //keys.reset_style("AT","compulsory");
  keys.add("numbered","KAPPA","KAPPA\\f$x\\f$ is equal to the value of the force constants at initial.");
  keys.reset_style("KAPPA","compulsory");
  keys.add("compulsory","TARGET","TARGET\\f$x\\f$ is the target point of optimization.");
  //keys.reset_style("TARGET","compulsory");
  keys.add("compulsory","CVSTEPSIZE","The stepsize of optimization");
  keys.add("compulsory","KAPPASTEPSIZE","The stepsize of optimization");
  keys.add("compulsory","OPTSTRIDE","The stride of optimization");
  keys.addFlag("OFFDIAGONAL",false," evaluate the kappa with offdiagonal");

  keys.addOutputComponent("work","default","the total work performed changing this restraint");
  keys.addOutputComponent("force2","default","the instantaneous value of the squared force due to this bias potential");
  keys.addOutputComponent("_cntr","default","one or multiple instances of this quantity can be referenced elsewhere in the input file. "
                          "these quantities will named with  the arguments of the bias followed by "
                          "the character string _cntr. These quantities give the instantaneous position "
                          "of the center of the harmonic potential.");
  keys.addOutputComponent("_work","default","one or multiple instances of this quantity can be referenced elsewhere in the input file. "
                          "These quantities will named with the arguments of the bias followed by "
                          "the character string _work. These quantities tell the user how much work has "
                          "been done by the potential in dragging the system along the various colvar axis.");
  keys.addOutputComponent("_kappa","default","one or multiple instances of this quantity can be referenced elsewhere in the input file. "
                          "These quantities will named with the arguments of the bias followed by "
                          "the character string _kappa. These quantities tell the user the time dependent value of kappa.");
  keys.addOutputComponent("_mean","default","mean of cvs. ");
}

OptimizeRestraintAdam::OptimizeRestraintAdam(const ActionOptions&ao):
  PLUMED_BIAS_INIT(ao),
  //verse(getNumberOfArguments())
  stepsize_cv(1),
  stepsize_kappa(1),
  optmethod(1),
  beta1(1),
  beta2(1),
  epsilon(1),
  stride(1),
  mean(getNumberOfArguments()),
  m_cv(getNumberOfArguments()),
  at(getNumberOfArguments()),
  target(getNumberOfArguments()),
  offdiagonal(false)
{
  //parseVector("VERSE",verse);

  parseVector("OPTMETHOD",optmethod);
  parseVector("BETA1",beta1);
  parseVector("BETA2",beta2);
  parseVector("EPSILON",epsilon);

  // Try to read step size
  parseVector("CVSTEPSIZE",stepsize_cv);
  parseVector("KAPPASTEPSIZE",stepsize_kappa);
  // Try to read stride
  parseVector("OPTSTRIDE",stride);
  // Try to read offdiagonal
  parseFlag("OFFDIAGONAL",offdiagonal);
  //std::vector<long int> ss(1); ss[0]=-1;
  std::vector<double> kk( getNumberOfArguments() ), aa( getNumberOfArguments() );
  //std::vector<double> targetp( getNumberOfArguments() );
  // Now read AT and TARGET
  parseVector("AT",at);
  parseVector("TARGET",target);

  for(int i=0;; i++) {
    // Try to read kappa
    //if( !parseNumberedVector("KAPPA",i,kk) ) kk=kappa[i-1];
    if( !parseNumberedVector("KAPPA",i,kk) ) break;
    kappa.push_back(kk);
    m_kappa.push_back(kk);
  }
  checkRead();
/*
  for(unsigned i=0; i<step.size(); i++) {
    log.printf("  step%u %ld\n",i,step[i]);
    log.printf("  at");
    for(unsigned j=0; j<at[i].size(); j++) log.printf(" %f",at[i][j]);
    log.printf("\n");
    log.printf("  with force constant");
    for(unsigned j=0; j<kappa[i].size(); j++) log.printf(" %f",kappa[i][j]);
    log.printf("\n");
  };
  this showld be add the output of log file about Kappa and at (future plan)*/
  log.printf("  reference point");
  for(unsigned j=0; j<at.size(); j++) log.printf(" %f",at[j]);
  log.printf("\n");
  log.printf("  target point");
  for(unsigned j=0; j<target.size(); j++) log.printf(" %f",target[j]);
  log.printf("\n");
  log.printf("  with force constant");
  if (offdiagonal) {
    for(unsigned i=0; i<kappa[0].size(); i++) {
      for(unsigned j=0; j<kappa[0].size(); j++) {
        log.printf(" %f",kappa[i][j]);
      }
      log.printf("\n");
    }
  } else {
    for(unsigned j=0; j<kappa[0].size(); j++) log.printf(" %f",kappa[0][j]);
  }
  log.printf("\n");
  log.printf("  with step size in CV");
  log.printf(" %f\n",stepsize_cv[0]);
  log.printf("  with step size in KAPPA");
  log.printf(" %f\n",stepsize_kappa[0]);
  log.printf("  with stride");
  log.printf(" %ld\n",stride[0]);

  addComponent("force2"); componentIsNotPeriodic("force2");

  // add the centers of the restraint as additional components that can be retrieved (useful for debug)

  std::string comp;
  for(unsigned i=0; i< getNumberOfArguments() ; i++) {
    comp=getPntrToArgument(i)->getName()+"_cntr"; // each spring has its own center
    addComponent(comp); componentIsNotPeriodic(comp);
    comp=getPntrToArgument(i)->getName()+"_work"; // each spring has its own work
    addComponent(comp); componentIsNotPeriodic(comp);
    if (offdiagonal) {
      for(unsigned j=0; j< getNumberOfArguments() ; j++) {
        comp=getPntrToArgument(i)->getName()+"_"+getPntrToArgument(j)->getName()+"_kappa"; // each spring has its own kappa
        addComponent(comp); componentIsNotPeriodic(comp);
      }
    } else{
      comp=getPntrToArgument(i)->getName()+"_kappa"; // each spring has its own kappa
      addComponent(comp); componentIsNotPeriodic(comp);
    }
    comp=getPntrToArgument(i)->getName()+"_mean"; // each spring has its own mean
    addComponent(comp); componentIsNotPeriodic(comp);
    work.push_back(0.); // initialize the work value
  }
  addComponent("work"); componentIsNotPeriodic("work");
  tot_work=0.0;
  N=0;
  t=1;
  v_cv = 0.0;
  v_kappa = 0.0;
  unsigned i=0;
  unsigned j=0;
  for(i=0; i<getNumberOfArguments(); ++i) {
    mean[i] = 0.0;
    m_cv[i] = 0.0;
    for(j=0; j<getNumberOfArguments(); ++j) {
      m_kappa[i][j] = 0.0;
    };
  };



  //log<<"  Bibliography ";
  //log<<cite("Grubmuller, Heymann, and Tavan, Science 271, 997 (1996)")<<"\n";


}


void OptimizeRestraintAdam::calculate() {
  double ene=0.0;
  double totf2=0.0;
  unsigned narg=getNumberOfArguments();
  long int now=getStep();
  std::vector<double> kk(narg),aa(narg),f(narg),dpotdk(narg);
  std::vector<double> aa_before(narg),kk_before(narg);
  std::vector<double> grad_cv(narg),grad_kappa_i(narg);
  std::vector<std::vector<double> > kk_offd(narg), kk_before_offd(narg);
  std::vector<std::vector<double> > grad_kappa(narg);

  //std::vector<double> mean(narg);
  //unsigned N;

  unsigned i=0;
  unsigned j=0;
  if(now==0) {
    if (offdiagonal) {
      kk_offd=kappa;
    } else {
      kk=kappa[0];
    }
    aa=at;
    oldaa=at;
    oldk=kappa[0];
    olddpotdk.resize(narg);
    oldf.resize(narg);
    N = 0;
    
    for(i=0; i<narg; ++i) {mean[i] = 0.0;};
  } else if (now % stride[0] == 0){
    // optimization of the reference point and kappa
    aa_before=at;
    for(i=0; i<narg; ++i) {
      mean[i] += aa_before[i];
    }
    if (optmethod[0]=="SGD") {
      if (offdiagonal) {
        kk_before_offd=kappa;
        for(i=0; i<narg; ++i) {
          const double mean_target_i=-difference(i,mean[i],target[i]); // this gives: mean[i] - target[i]
          const double ref_target_i=-difference(i,aa_before[i],target[i]); // this gives: aa[i] - target[i]
          for(j=0; j<narg; ++j) {
            const double mean_target_j=-difference(j,mean[j],target[j]); // this gives: mean[i] - target[i]
            const double ref_mean_j=-difference(j,aa_before[j],mean[j]); // this gives: aa[i] - mean[i]
            at[i] -= stepsize_cv[0]*mean_target_i*kk_before_offd[i][j];
            grad_kappa[i][j]  -= mean_target_j*ref_target_i*0.5; 
            grad_kappa[i][j]  -= mean_target_i*ref_mean_j*0.5; 
          }
        }
      } else {
        kk_before=kappa[0];
        for(i=0; i<narg; ++i) {
          const double mean_target_i=-difference(i,mean[i],target[i]); // this gives: mean[i] - target[i]
          const double ref_target_i=-difference(i,aa_before[i],target[i]); // this gives: aa[i] - target[i]
          const double ref_mean_i=-difference(i,aa_before[i],mean[i]); // this gives: aa[i] - target[i]
          at[i] -= stepsize_cv[0]*mean_target_i*kk_before[i];
          kappa[0][i] -= stepsize_kappa[0]*mean_target_i*(ref_target_i+ref_mean_i)*0.5;
        };
      } // end optimization of Steepest Gradient Descent
    } else {
      if (offdiagonal) {
        kk_before_offd=kappa;
        grad_kappa=kappa;
        double gradgrad_cv = 0.0;
        double gradgrad_kappa = 0.0;
        for(i=0; i<narg; ++i) {
          grad_cv[i] = 0.0;
          const double mean_target_i=-difference(i,mean[i],target[i]); // this gives: mean[i] - target[i]
          const double ref_target_i=-difference(i,aa_before[i],target[i]); // this gives: aa[i] - target[i]
          for(j=0; j<narg; ++j) {
            const double mean_target_j=-difference(j,mean[j],target[j]); // this gives: mean[i] - target[i]
            const double ref_mean_j=-difference(j,aa_before[j],mean[j]); // this gives: aa[i] - mean[i]

            grad_cv[i] += mean_target_j*kk_before_offd[i][j];
            grad_kappa[i][j]  = mean_target_j*ref_target_i*0.5; 
            grad_kappa[i][j]  += mean_target_i*ref_mean_j*0.5; 
            //if (t==1){
              //m_kappa[i][j] = grad_kappa[i][j];
            //} else {
              m_kappa[i][j] = m_kappa[i][j]*beta1[0]+(1.0-beta1[0])*grad_kappa[i][j];
            //}
            if (optmethod[0]=="ADABELIEF") {
              gradgrad_kappa += (grad_kappa[i][j]-m_kappa[i][j])*(grad_kappa[i][j]-m_kappa[i][j]);
            } else {
              gradgrad_kappa += grad_kappa[i][j]*grad_kappa[i][j];
            }
          }
          //if (t==1){
            //m_cv[i] = grad_cv[i];
          //}else{
            m_cv[i] = m_cv[i]*beta1[0]+(1.0-beta1[0])*grad_cv[i];
          //}
          if (optmethod[0]=="ADABELIEF") {
            gradgrad_cv += (grad_cv[i]-m_cv[i])*(grad_cv[i]-m_cv[i]);
          } else {
            gradgrad_cv += grad_cv[i]*grad_cv[i];
          }
        }
        //if (t==1){
          //v_cv = gradgrad_cv;
          //v_kappa = gradgrad_kappa;
        //}else{
          v_cv = v_cv*beta2[0]+(1.0-beta2[0])*gradgrad_cv;
          v_kappa = v_kappa*beta2[0]+(1.0-beta2[0])*gradgrad_kappa;
        //}
        if ((optmethod[0]=="ADAM") || (optmethod[0]=="NADAM") || (optmethod[0]=="ADABELIEF")) {
          double b = sqrt(1.0-pow(beta2[0],t))/(1.0-pow(beta1[0],t));
          for(i=0; i<narg; ++i) {
            if ((optmethod[0]=="ADAM") || (optmethod[0]=="ADABELIEF")) {
              at[i] -= stepsize_cv[0]*m_cv[i]/sqrt(v_cv+epsilon[0])*b;
            } else {
              at[i] -= stepsize_cv[0]*(beta1[0]*m_cv[i]+(1.0-beta1[0])*grad_cv[i])/sqrt(v_cv+epsilon[0])*b;
            }
            for(j=0; j<narg; ++j) {
              //double mcap = m_kappa[i][j]/(1.0-pow(beta1[0],t));
              //double vcap = v_kappa/(1.0-pow(beta2[0],t));
              if ((optmethod[0]=="ADAM") || (optmethod[0]=="ADABELIEF")) {
                kappa[i][j] -= stepsize_kappa[0]*m_kappa[i][j]/sqrt(v_kappa+epsilon[0])*b;
              } else {
                kappa[i][j] -= stepsize_kappa[0]*(beta1[0]*m_kappa[i][j]+(1.0-beta1[0])*grad_kappa[i][j])/sqrt(v_kappa+epsilon[0])*b;
              }
            }
          }
        } else if (optmethod[0]=="MOMENTUM") {
          for(i=0; i<narg; ++i) {
            at[i] -= stepsize_cv[0]*m_cv[i];
            for(j=0; j<narg; ++j) {
              kappa[i][j] -= stepsize_kappa[0]*m_kappa[i][j];
            }
          }
        } else if (optmethod[0]=="RMSPROP") {
          for(i=0; i<narg; ++i) {
            at[i] -= stepsize_cv[0]*grad_cv[i]/sqrt(v_cv+epsilon[0]);
            for(j=0; j<narg; ++j) {
              kappa[i][j] -= stepsize_kappa[0]*grad_kappa[i][j]/sqrt(v_kappa+epsilon[0]);
            }
          }
        } else error("OPTMETHOD cannot understand");
      } else {
        kk_before=kappa[0];
        double gradgrad_cv = 0.0;
        double gradgrad_kappa = 0.0;
        for(i=0; i<narg; ++i) {
          const double mean_target_i=-difference(i,mean[i],target[i]); // this gives: mean[i] - target[i]
          const double ref_target_i=-difference(i,aa_before[i],target[i]); // this gives: aa[i] - target[i]
          const double ref_mean_i=-difference(i,aa_before[i],mean[i]); // this gives: aa[i] - target[i]
          grad_cv[i] = mean_target_i*kk_before[i];
          grad_kappa[0][i] = mean_target_i*(ref_target_i+ref_mean_i)*0.5;
          //if (t==1){
            //m_cv[i] = grad_cv[i];
            //m_kappa[0][i] = grad_kappa[0][i];
          //} else { 
            m_cv[i] = m_cv[i]*beta1[0]+(1.0-beta1[0])*grad_cv[i];
            m_kappa[0][i] = m_kappa[0][i]*beta1[0]+(1.0-beta1[0])*grad_kappa[0][i];
          //}
          if (optmethod[0]=="ADABELIEF") {
            gradgrad_cv += (grad_cv[i]-m_cv[i])*(grad_cv[i]-m_cv[i]);
            gradgrad_kappa += (grad_kappa[0][i]-m_kappa[0][i])*(grad_kappa[0][i]-m_kappa[0][i]);
          } else {
            gradgrad_cv += grad_cv[i]*grad_cv[i];
            gradgrad_kappa += grad_kappa[0][i]*grad_kappa[0][i];
          }
        };
        v_cv = v_cv*beta2[0]+(1.0-beta2[0])*gradgrad_cv;
        v_kappa = v_kappa*beta2[0]+(1.0-beta2[0])*gradgrad_kappa;
        if ((optmethod[0]=="ADAM") || (optmethod[0]=="NADAM") || (optmethod[0]=="ADABELIEF")) {
          double b = sqrt(1.0-pow(beta2[0],t))/(1.0-pow(beta1[0],t));
          for(i=0; i<narg; ++i) {
            if ((optmethod[0]=="ADAM") || (optmethod[0]=="ADABELIEF")) {
              at[i] -= stepsize_cv[0]*m_cv[i]/sqrt(v_cv+epsilon[0])*b;
              kappa[0][i] -= stepsize_kappa[0]*m_kappa[0][i]/sqrt(v_kappa+epsilon[0])*b;
            } else {
              at[i] -= stepsize_cv[0]*m_cv[i]/sqrt(v_cv+epsilon[0])*b;
              kappa[0][i] -= stepsize_kappa[0]*(beta1[0]*m_kappa[0][i]+(1.0-beta1[0]*grad_kappa[0][i]))/sqrt(v_kappa+epsilon[0])*b;
            }
          }
        } else if (optmethod[0]=="MOMENTUM") {
          for(i=0; i<narg; ++i) {
            at[i] -= stepsize_cv[0]*m_cv[i];
            kappa[0][i] -= stepsize_kappa[0]*m_kappa[0][i];
          }
        } else if (optmethod[0]=="RMSPROP") {
          for(i=0; i<narg; ++i) {
            at[i] -= stepsize_cv[0]*grad_cv[i]/sqrt(v_cv+epsilon[0]);
            kappa[0][i] -= stepsize_kappa[0]*grad_kappa[0][i]/sqrt(v_kappa+epsilon[0]);
          }
        } else error("OPTMETHOD cannot understand");
      t += 1;
      } // end optimization of Adam
    }
    for(i=0; i<narg; ++i) {mean[i] = 0.0;};
    N = 0;
  }
  tot_work=0.0;
  N += 1;
  for(unsigned i=0; i<narg; ++i) {
    aa=at;
    const double cv=difference(i,aa[i],getArgument(i)); // this gives: getArgument(i) - at[0][i]
    //mean[i] = mean[i]*(N-1)/N+(cv+aa[i])/N;
    mean[i] = mean[i]*(N-1)/N+cv/N;
    getPntrToComponent(getPntrToArgument(i)->getName()+"_cntr")->set(aa[i]);
    getPntrToComponent(getPntrToArgument(i)->getName()+"_mean")->set(mean[i]+aa[i]);
    if (offdiagonal) {
      f[i]=0.0;
      for(unsigned j=0; j<narg; ++j) {
        const double cv_j=difference(j,aa[j],getArgument(j)); // this gives: getArgument(i) - at[0][i]
        const double k=kappa[i][j];
        f-=k*cv_j;
        ene+=0.5*k*cv*cv_j;
        getPntrToComponent(getPntrToArgument(i)->getName()+"_"+getPntrToArgument(j)->getName()+"_kappa")->set(kappa[i][j]);
      }

    } else {
      const double k=kappa[0][i];
      f[i]=-k*cv;
      dpotdk[i]=0.5*cv*cv;
      if(oldaa.size()==aa.size() && oldf.size()==f.size()) work[i]+=0.5*(oldf[i]+f[i])*(aa[i]-oldaa[i]) + 0.5*( dpotdk[i]+olddpotdk[i] )*(kk[i]-oldk[i]);
      getPntrToComponent(getPntrToArgument(i)->getName()+"_kappa")->set(kappa[0][i]);
      ene+=0.5*k*cv*cv;
    }
    getPntrToComponent(getPntrToArgument(i)->getName()+"_work")->set(work[i]);
    tot_work+=work[i];
    setOutputForce(i,f[i]);
    totf2+=f[i]*f[i];
  };
  getPntrToComponent("work")->set(tot_work);
  oldf=f;
  oldaa=aa;
  oldk=kk;
  olddpotdk=dpotdk;
  setBias(ene);
  getPntrToComponent("force2")->set(totf2);
}

}
}


