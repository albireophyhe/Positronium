#include <iostream>
#include <stdlib.h>

#include <TRandom3.h>
#include <TVector3.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TMath.h>

using namespace std;

const Int_t MAX_NGAMMA = 3;
const Int_t CHANNEL = 6;
const Float_t WholeEne = 1022.;
const Float_t OPsratio = 0.7365;
const Float_t sig = 4.2;
const TVector3 v0;
const Float_t pcomp = 0.8;
const Float_t Z = 64;
Float_t PMTphi[6] = {-TMath::Pi()*5./6.,-TMath::Pi()/2.,-TMath::Pi()/6.,TMath::Pi()/6.,TMath::Pi()/2.,TMath::Pi()*5./6.};
TRandom3 erand(0);
//TRandom3 therand(1);

int DecidePMT(TVector3 v){
  Int_t nPMT = 6;
  for (int nP = 0; nP < 6; nP++) {
    TVector3 vPMT;
    vPMT.SetMagThetaPhi(1.,TMath::Pi()/2.,PMTphi[nP]);
    if(v.Angle(vPMT)<TMath::Pi()/6.){
      nPMT = nP;
    }
  }
  return nPMT;
}


//Systematic error from PMT
float PMTEne(Float_t meanE,Float_t sigma){
  while(true){
    Float_t E = erand.Gaus(meanE,sigma);
    if(E>0){
      return E;
      break;
    }
  }
  return 0;
}

//rest Positronium number by lifetime
float Restnump(Float_t tim){
  Float_t lf=0.125;
  Float_t Rp;
  Float_t insu;
  insu =-1* tim/lf;
  Rp = TMath::Exp(insu);
  //cout << "restp:" << Rp << endl;
  return Rp;
}

float Restnumo(Float_t tim){
  Float_t lf=142.0;
  Float_t Rp;
  Float_t insu;
  insu =-1* tim/lf;
  Rp = TMath::Exp(insu);
  //cout << "resto:" << Rp << endl;
  return Rp;
}

//Count compton scattering
float Compton(Float_t Ei){
  Float_t a;
  a = Ei/511.;
  Float_t Tth;
  Float_t Tcomp;
  Float_t the;
  Float_t ratio;
  Float_t para;
  //  int num=0;
  Tcomp = TMath::Abs(((1+a)/pow(a,2)*(2*(1+a)/(1+2*a)-TMath::Log(1+2*a)/a)+TMath::Log(1+2*a)/(2*a)-(1+3*a)/pow((1+2*a),2)));
  while(true){
    //num +=1;
    the = erand.Rndm(1)*TMath::Pi();
    para = erand.Rndm(1);
    Tth =0.5*pow((1./(1+a*(1-TMath::Cos(the)))),2)*(1+pow(TMath::Cos(the),2)+(pow(a,2)*(pow(1-TMath::Cos(the),2)))/(1+a*(1-TMath::Cos(the))))*TMath::Sin(the);
    ratio = Tth/Tcomp;
    //cout<< ratio<<","<<para<< endl;
    //cout << Tcomp << endl;
    //cout << Tth <<endl;
    if(ratio>1){
      cout << ratio << endl;
    }
    if(para<ratio){
      Float_t Ef;
      Ef =Ei- Ei/(a*(1-TMath::Cos(the))+1);
      return Ef;
      break;
    }
  }
  //cout << "CompEnd" << endl;
  return 0;
}

//This Photon reaction probability may not accurate, so all comment out in main code
float photon(Float_t Ei){
    Double_t Tpe;
    Double_t A;
    Double_t B;
    Double_t C;
    Double_t Tcomp;
    Double_t pcomp;
    Double_t a;
    a=Ei/511.;
    Tcomp = Z*0.49893*TMath::Abs(((1+a)/pow(a,2)*(2*(1+a)/(1+2*a)-TMath::Log(1+2*a)/a)+TMath::Log(1+2*a)/(2*a)-(1+3*a)/pow((1+2*a),2)));
    
    if(Ei>10){
        if(Ei<50){
            A = (1.53*TMath::Exp(-0.0361*Z)+0.510*TMath::Exp(-0.176*Z))*pow(10,-4);
            B = 0.0151*Z-0.18*TMath::Exp(-0.215*Z)-3.32;
            Tpe = A*pow(Z,5)*pow((Ei/30.),B);
            pcomp = Tcomp/(Tpe+Tcomp);
            return pcomp;
        }else if(Ei<250){
            A = (2.73*TMath::Exp(-0.2113*Z)+0.860*TMath::Exp(-0.126*Z))/(1-(2.63/Z+0.473)*(Ei-100)*pow(10,-3))*pow(10,-6);
            B = 0.0080*Z-0.18*TMath::Exp(-0.107*Z)-3.32;
            Tpe = A*pow(Z,5)*pow((Ei/30.),B);
            pcomp = Tcomp/(Tpe+Tcomp);
            return pcomp;
        }else{
            A=-18.849+5.798*TMath::Log(Z)-0.1754*pow(TMath::Log(Z),2);
            B=-2.6204+1.612*pow(10,-2)*Z+2.041*pow(10,-5)*pow(Z,2);
            C=0.4069-4.796*pow(10,-4)*Z-1.837*pow(10,-5)*pow(Z,2);
            Tpe = TMath::Exp(A+B*Ei+C*pow(Ei,2));
            pcomp = Tcomp/(Tpe+Tcomp);
            cout << A << " " << B << " " << C << " " << Tpe << Tcomp << endl;
            cout << pcomp <<endl;
            return pcomp;
        }
    }
    return 0;
}

int main(int argc, char **argv) {
  //Initialize
  Int_t nRun = 20;
  Int_t nEvents = 10000;
  if(argc>1){
    nRun = atoi(argv[1]);
    if(argc>2){
      nEvents = atoi(argv[2]);
    }
  }
  cout << "Run:" << nRun << endl;
  cout << "Events:" << nRun*nEvents << endl;

  //Prepare output ntup
  TFile fout("tgamma.root","recreate");
  TTree tree("pos","Positronium time event tree");
  Int_t nGamma;
  Float_t ch[6];
  Float_t t;
  Float_t pr[3];
  Float_t phi[3];
  Float_t theta[3];

  tree.Branch("nGamma",&nGamma,"nGamma/I");
  tree.Branch("ch",&ch,"ch[6]/F");
  tree.Branch("time",&t,"t/F");
  tree.Branch("pr",&pr,"pr[3]/F");
  tree.Branch("phi",&phi,"phi[3]/F");
  tree.Branch("theta",&theta,"theta[3]/F");

  //main loop
  cout << "sim start" <<endl;
  Int_t posi[nEvents][2];
  //TRandom3 rand(nRun);
  Float_t OorP;
  for(int nposi = 0; nposi<nEvents; nposi++){
    OorP = erand.Rndm(1);
    if(OPsratio<OorP){
      posi[nposi][0] = 0;
    }else{
      posi[nposi][0] = 1;
    }
    posi[nposi][1] = 1;
  }
  cout << "o,p decided" << endl;

  Float_t tin=0.;//time unit nsec
  Int_t nDone=0;
  while(true){
    Float_t EorP;
    Int_t PMTv1;
    Int_t PMTv2;
    Int_t PMTv3;
    Float_t restp=1.;
    Float_t resto=1.;
    Float_t Pdecay;
    Float_t lft;
    Float_t tstep=0.001;
    cout << tin << endl;
    if(tin==0){
      restp = 1;
      resto = 1;
    }else{
      restp = Restnump(tin)/Restnump(tin-tstep);
      resto = Restnumo(tin)/Restnumo(tin-tstep);
    }
    cout << restp << "," << resto << endl;
    for(int i = 0; i<nEvents; i++){
      if(posi[i][1]==1){
        Float_t EPMT=0;
        Pdecay=erand.Rndm(1);
	//cout << Pdecay << endl;
        if(posi[i][0]==0){
          lft = restp;
          if(Pdecay>lft){
            posi[i][1]=0;
	    nDone+=1;
            nGamma = 2;
            t=tin;
            TVector3 v1;
            v1.SetMagThetaPhi(511.,TMath::ACos(2*(erand.Rndm(1)-0.5)),2*(erand.Rndm(1)-0.5)*TMath::Pi());
            TVector3 v2;
            v2 = v0 - v1;
            PMTv1 = DecidePMT(v1);
            PMTv2 = DecidePMT(v2);
            if(PMTv1<6){
              phi[0] = v1.Phi();
              theta[0] = v1.Theta();
              pr[0] = v1.Mag();
              EorP = erand.Rndm(1);
              EPMT=v1.Mag();
              //porc = photon(EPMT);
              if(EorP<pcomp){
                EPMT=Compton(EPMT);
              }
              ch[PMTv1] = PMTEne(EPMT,sig);
            }
            if(PMTv2<6){
              phi[1] = v2.Phi();
              theta[1] = v2.Theta();
              pr[1] = v2.Mag();
              EorP = erand.Rndm(1);
              EPMT= v2.Mag();
              //porc = photon(EPMT);
              if(EorP<pcomp){
                EPMT=Compton(EPMT);
              }
              ch[PMTv2]=PMTEne(EPMT,sig);
            }
            //cout << "2" << endl;
            tree.Fill();

          }
        }else{
          lft = resto;
          if(Pdecay>lft){
            posi[i][1]=0;
	    nDone+=1;
            t=tin;
            nGamma = 3;
            TVector3 v1;
            TVector3 v2;
            TVector3 v3;
            Double_t E1;
            Double_t E2;
            Double_t E3;

              
            while(1){
                while(1){
                    E2=rand.Rndm(1)*511.;
                    E3=rand.Rndm(1)*511.;
                    E1=1022.-E2-E3;
                    if(E1<511.){
                        break;
                    }
                }
                Double_t dif3;
                dif3 = E2*E2+E3*E3-E1*E1;
                if(TMath::Abs(dif3/(2*E2*E3))<1.){
                    th = TMath::Pi()-TMath::ACos(dif3/(2*E2*E3));
                    v2.SetMagThetaPhi(E2,TMath::ACos(2*(erand.Rndm(1)-0.5)),2*(erand.Rndm(1)-0.5)*TMath::Pi());
                    if(v2.Phi()-th>0){
                        v3.SetMagThetaPhi(E3,v2.Theta()-th,v2.Phi());
                    }else{
                        v3.SetMagThetaPhi(E3,v2.Theta()+th,v2.Phi());
                    }
                    v3.Rotate(2*(erand.Rndm(1)-0.5)*TMath::Pi(),v2);
                    v1 = v0-(v2+v3);
                    break;
                }
            }

            PMTv1 = DecidePMT(v1);
            PMTv2 = DecidePMT(v2);
            PMTv3 = DecidePMT(v3);

            if((PMTv1!=PMTv2)&&(PMTv2!=PMTv3)&&(PMTv1!=PMTv3)){
              if(PMTv1<6){
                phi[0] = v1.Phi();
                theta[0] = v1.Theta();
                pr[0] = v1.Mag();
                EorP = erand.Rndm(1);
                EPMT=v1.Mag();
                //porc=photon(EPMT);
                if(EorP<pcomp && EPMT>10){
                  EPMT=Compton(EPMT);
                }
                //cout << "31" << endl;
                ch[PMTv1]=PMTEne(EPMT,sig);
              }
              if(PMTv2<6){
                phi[1] = v2.Phi();
                theta[1] = v2.Theta();
                pr[1] = v2.Mag();
                EorP = erand.Rndm(1);
                EPMT=v2.Mag();
                //porc=photon(EPMT);
                if(EorP<pcomp && EPMT>10){
                  EPMT=Compton(EPMT);
                }
                //cout << "32" << endl;
                ch[PMTv2]=PMTEne(EPMT,sig);
              }
              if(PMTv3<6){
                phi[2] = v3.Phi();
                theta[2] = v3.Theta();
                pr[2] = v3.Mag();
                EorP = erand.Rndm(1);
                EPMT = v3.Mag();
                //porc=photon(EPMT);
                if(EorP<pcomp &&EPMT>10){
                  EPMT=Compton(EPMT);
                }
                //cout << "33" << endl;
                ch[PMTv3]=PMTEne(EPMT,sig);
              }
              //cout << "3" << endl;
              tree.Fill();
            }else if(PMTv1==PMTv2){
              if(PMTv1<6){
                phi[0] = v1.Phi();
                theta[0] = v1.Theta();
                pr[0] = v1.Mag();
                EorP = erand.Rndm(1);
                EPMT=v1.Mag();
                phi[1] = v2.Phi();
                theta[1] = v2.Theta();
                pr[1] = v2.Mag();
                //porc=photon(EPMT);
                if(EorP<pcomp && EPMT>10){
                  EPMT=Compton(EPMT);
                }
                EorP = erand.Rndm(1);
                if(EorP<pcomp && v2.Mag()>10){
                  EPMT+=Compton(v2.Mag());
                }else{
                  EPMT+=v2.Mag();
                }
                ch[PMTv1]=PMTEne(EPMT,sig);
              }
              if(PMTv3<6){
                phi[2] = v3.Phi();
                theta[2] = v3.Theta();
                pr[2] = v3.Mag();
                EorP = erand.Rndm(1);
                EPMT = v3.Mag();
                //porc=photon(EPMT);
                if(EorP<pcomp &&EPMT>10){
                  EPMT=Compton(EPMT);
                }
                //cout << "33" << endl;
                ch[PMTv3]=PMTEne(EPMT,sig);
              }
              tree.Fill();
            }else if(PMTv1==PMTv3){
              if(PMTv1<6){
                phi[0] = v1.Phi();
                theta[0] = v1.Theta();
                pr[0] = v1.Mag();
                EorP = erand.Rndm(1);
               EPMT=v1.Mag();
               phi[2] = v3.Phi();
               theta[2] = v3.Theta();
               pr[2] = v3.Mag();
               //porc=photon(EPMT);
               if(EorP<pcomp && EPMT>10){
                 EPMT=Compton(EPMT);
                }
                EorP = erand.Rndm(1);
                if(EorP<pcomp && v3.Mag()>10){
                  EPMT+=Compton(v3.Mag());
                }else{
                  EPMT+=v3.Mag();
                }
                ch[PMTv1]=PMTEne(EPMT,sig);
              }
              if(PMTv2<6){
                phi[1] = v2.Phi();
                theta[1] = v2.Theta();
                pr[1] = v2.Mag();
                EorP = erand.Rndm(1);
                EPMT = v2.Mag();
                //porc=photon(EPMT);
                if(EorP<pcomp &&EPMT>10){
                  EPMT=Compton(EPMT);
                }
                //cout << "33" << endl;
                ch[PMTv2]=PMTEne(EPMT,sig);
              }
              tree.Fill();
            }else if(PMTv2==PMTv3){
              if(PMTv2<6){
                phi[1] = v2.Phi();
                theta[1] = v2.Theta();
                pr[1] = v2.Mag();
                EorP = erand.Rndm(1);
                EPMT=v2.Mag();
                phi[2] = v3.Phi();
                theta[2] = v3.Theta();
                pr[2] = v3.Mag();
                //porc=photon(EPMT);
                if(EorP<pcomp && EPMT>10){
                  EPMT=Compton(EPMT);
                }
                EorP = erand.Rndm(1);
                if(EorP<pcomp && v3.Mag()>10){
                  EPMT+=Compton(v3.Mag());
                }else{
                  EPMT+=v3.Mag();
                }
                ch[PMTv2]=PMTEne(EPMT,sig);
              }
              if(PMTv1<6){
                phi[0] = v1.Phi();
                theta[0] = v1.Theta();
                pr[0] = v1.Mag();
                EorP = erand.Rndm(1);
                EPMT = v1.Mag();
                //porc=photon(EPMT);
                if(EorP<pcomp &&EPMT>10){
                  EPMT=Compton(EPMT);
                }
                //cout << "33" << endl;
                ch[PMTv1]=PMTEne(EPMT,sig);
              }
              tree.Fill();
            }
          }
        }
      }
    }

    cout << tin << endl;
    tin=tin+tstep;
    if(tin>400){
      break;
    }
    if(nDone==nEvents){
      break;
    }
  }


  cout << "sim over" << endl;
  tree.Write();
  fout.Close();

  return 0;
}
