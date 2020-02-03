/**
 * @file input.cc
 * @author Luca Maccione, Daniele Gaggero
 * @email luca.maccione@desy.de
 * @email daniele.gaggero@sissa.it
 * @brief File where class Input is implemented.
 */

#include "input.h"
#include "constants.h"
#include "tinyxml.h"
#include "errorcode.h"

#include <string>
#include <stdlib.h>
#include <algorithm>
#include <sstream>
#include <vector>

using namespace std;

vector<double> Input::ReadInjectionStructure(string name, TiXmlElement* el1) {
   
  double result = 1.;
  int rho_counter = 0;
  vector<double> hold;
  while (result > 0) {
    stringstream rho_name;
    rho_name << name << rho_counter;
    result = QueryDoubleAttribute(rho_name.str(),el1);
    if (result>0) {
      cout << "Reading " << rho_name.str() << " = " << result << endl;
      hold.push_back(result);
      rho_counter++;
    }
  }
  return hold;
}

int Input::QueryIntAttribute(string obj, TiXmlElement* el) {
  TiXmlElement* el1 = el->FirstChildElement(obj);
  int result = -1;
  if (el1) el1->QueryIntAttribute("value", &result);
  return result;
}

int Input::QueryIntAttributeWithDefault(string obj, TiXmlElement* el, int def) {
  TiXmlElement* el1 = el->FirstChildElement(obj);
  int result = def;
  if (el1) el1->QueryIntAttribute("value", &result);
  return result;
}

double Input::QueryDoubleAttribute(string obj, TiXmlElement* el) {
  TiXmlElement* el1 = el->FirstChildElement(obj);
  double result = -1.0;
  if (el1) el1->QueryDoubleAttribute("value", &result);
  return result;
}

double Input::QueryDoubleAttributeWithDefault(string obj, TiXmlElement* el, double def) {
  TiXmlElement* el1 = el->FirstChildElement(obj);
  double result = def;
  if (el1) el1->QueryDoubleAttribute("value", &result);
  return result;
}

string Input::QueryStringAttribute(string obj, TiXmlElement* el) {
  TiXmlElement* el1 = el->FirstChildElement(obj);
  string result ("error");
  if (el1) el1->QueryStringAttribute("value", &result);
  return result;
}

int Input::LoadFile(const string inputfilename) {
   
  TiXmlDocument* doc = new TiXmlDocument(inputfilename.c_str());
   
  ringmin = 0;
  ringmax = 0;
  TiXmlElement* el = NULL;
  TiXmlElement* el1 = NULL;
  TiXmlElement* el2 = NULL;
   
  if (doc->LoadFile(inputfilename.c_str(), TIXML_DEFAULT_ENCODING)) {
      
    // I/O
    el = doc->FirstChildElement("Output");
      
    int found = inputfilename.rfind('/');
    if (found == string::npos) found = -1;
    filename = string(inputfilename, found+1, inputfilename.size()-4-(found+1));
      
    partialstore = ( el->FirstChildElement("partialstore") );
    fullstore = ( el->FirstChildElement("fullstore") );
    if (!(fullstore || partialstore)) cerr << "Warning: your results will not be saved." << endl;
      
    //MW130801
    timestepstore = ( el->FirstChildElement("timestepstore") );
    feedback = QueryIntAttribute("feedback", el);

    // Grid
    Rmin = 0.0;
    el = doc->FirstChildElement("Grid");
    if (el) {el->QueryStringAttribute("type", &gridtype);
      
      el1 = el->FirstChildElement("Observer");
      if (el1){
	xobs = QueryDoubleAttribute("x", el1);
	yobs = QueryDoubleAttribute("y", el1);
	zobs = QueryDoubleAttribute("z", el1);
	robs = sqrt(xobs*xobs + yobs*yobs);
	if (feedback >0) cout << "Observer is sitting at x = " << xobs << " y = " << yobs << " z = " << zobs << endl;
      }
      else{ //DG28.11.2013
	if (feedback >0) cout<<"Using default values for the position of the observer: x = 8.3 kpc; y = 0; z= 0 " << endl;
	if (feedback >0) cout<<"You can specify a <Observer> block with x y and z coordinates if you want to change these values" << endl;
	xobs = 8.3; yobs = 0.; zobs = 0.; robs=xobs;
	//exit(110);
      }
      
      Rmax = QueryDoubleAttribute("Rmax", el);
      zmax = QueryDoubleAttribute("L", el);
        
      numr = QueryIntAttribute("DimR", el);
      numz = QueryIntAttribute("DimZ", el);
        
      dimz_equidistant = true;
      el1 = el->FirstChildElement("DimZ_division");
      if (el1)
	{
	  dimz_equidistant = false;
	  el1->QueryStringAttribute("points", &divz);
	  numz = count(divz.begin(), divz.end(), ';') + 1; //reset dimension
	}
      
      if (gridtype == "3D") {
	// Matze: nonequidistant z grid
	numx = QueryIntAttribute("DimX", el);
	dimx_equidistant = true;
	el1 = el->FirstChildElement("DimX_division");
	if (el1)
	  {
	    dimx_equidistant = false;
	    el1->QueryStringAttribute("points", &divx);
	    numx = std::count(divx.begin(), divx.end(), ';') + 1;
	  }
            
	numy = QueryIntAttribute("DimY", el);
	dimy_equidistant = true;
	el1 = el->FirstChildElement("DimY_division");
	if (el1)
	  {
	    dimy_equidistant = false;
	    el1->QueryStringAttribute("points", &divy);
	    numy = std::count(divy.begin(), divy.end(), ';') + 1;
	  }

      }
      else
	{
	  dimr_equidistant = true;
	  el1 = el->FirstChildElement("DimR_division");
	  if (el1)
            {
	      dimr_equidistant = false;
	      el1->QueryStringAttribute("points", &divr);
	      numr = std::count(divr.begin(), divr.end(), ';') + 1;
            }
	}
      
      Ekmin = QueryDoubleAttribute("Ekmin", el);
      Ekmax = QueryDoubleAttribute("Ekmax", el);
      Ekfact = QueryDoubleAttribute("Ekfactor", el);
      
      el1 = el->FirstChildElement("NuclearChain");
      if (el1) {
	Zmax = QueryIntAttribute("Zmax", el1);
	Zmin = QueryIntAttribute("Zmin", el1);
        if ( Zmax != 0 ) prop_nuclei = true; else prop_nuclei = false; // CHECK
	prop_secap = ( el1->FirstChildElement("PropSecAntiProton") );
	prop_secdeuteron = ( el1->FirstChildElement("PropSecAntiDeuteron") );
	prop_lep = ( el1->FirstChildElement("PropLepton") );
	prop_TPP = ( el1->FirstChildElement("PropTPP") );
	prop_extracomp = ( el1->FirstChildElement("PropExtraComponent") );
	Kcapture = ( el1->FirstChildElement("Kcapture") );
      }
      else {
	cerr << "You must specify the Nuclear Chain!" << endl;
	exit(NONUCCHAIN);
      }
         
      DoubleRun = false;
      el1 = el->FirstChildElement("DoubleRun");
      if  (el1) {
	DoubleRun = true;
	if (feedback >0) cout << "A second iteration of the nuclear chain will be performed in order to take the beta- decays into account properly. " << endl;
      }
         
      LB_shape = "Cuboid";
      LB_ax = LB_ay = LB_az = -1; //MW: if not set, let grid.cc decide over default values
      el1 = el->FirstChildElement("LocalBubble");
      if (el1) {
	el1->QueryStringAttribute("type", &LB_shape);
	LB_ax = QueryDoubleAttribute("SemiX", el1);
	LB_ay = QueryDoubleAttribute("SemiY", el1);
	LB_az = QueryDoubleAttribute("SemiZ", el1);
            
	LB_smearing = "None";
	el2 = el1->FirstChildElement("Smearing");
	if (el2) el2->QueryStringAttribute("type", &LB_smearing);
      }
         
    }
    else {
      cerr << "You must specify a GRID!" << endl;
      exit(NOGRID);
    }
      
      
    // Algorithm
    TESTMODE = false;
    el = doc->FirstChildElement("Algorithm");
    if (!el) {
      cerr << "You did not set specific algorithm. It will be set to default" << endl;
      Nrept = 20;
      dtfactor = 0.25;
      dtmin = 1.e-6;
      dtmax = 1.e3;
      ADI = false;
      OpSplit = true;
      tol = 1.e-4;
      alpha = 0.09;
    }
    else {
      el1 = el->FirstChildElement("OpSplit");
      OpSplit = (el1);
      if (OpSplit) {
	Nrept = QueryIntAttribute("Nrept", el1);
	dtfactor = QueryDoubleAttribute("Dtfactor", el1);
	dtmin = QueryDoubleAttribute("Dtmin", el1);
	dtmax = QueryDoubleAttribute("Dtmax", el1);
      }
      el1 = el->FirstChildElement("ADI");
      ADI = (el1);
      if (ADI) {
	tol = QueryDoubleAttribute("Tolerance", el1);
	alpha = QueryDoubleAttribute("Alpha", el1);
            
	//MW: not functional at the moment
	cerr << "ADI is not functional at the moment! Because of nonequidistant Grid. Greets, Matthias." << endl;
	ADI=false;
      }
      el1 = el->FirstChildElement("TestMode");
      if (el1) TESTMODE = true;
      else     TESTMODE = false;
      el1 = el->FirstChildElement("MovingSource");
      if (el1) MOVING = true;
      else     MOVING = false;
      if (MOVING) {
	source_x0 = QueryDoubleAttribute("x0", el1);
	source_y0 = QueryDoubleAttribute("y0", el1);
	source_z0 = QueryDoubleAttribute("z0", el1);
	source_vx = QueryDoubleAttribute("source_vx", el1);
	source_vy = QueryDoubleAttribute("source_vy", el1);
	source_vz = QueryDoubleAttribute("source_vz", el1);
      }
    }
      
    // Galaxy
    el = doc->FirstChildElement("Galaxy");
      
    string gmod;
      
    // Gas
    el1 = el->FirstChildElement("Gas");
    if (el1) el1->QueryStringAttribute("type", &gmod);
    if (gmod == "BronfFerr") gas_model = BronfFerr;
    else if (gmod == "NS") gas_model = NS;
    else if (gmod == "Galprop") gas_model = Galprop;
    else if (gmod == "Uniform") gas_model = UniformGas;
    else cerr << "Your gas model is not implemented. You chose: " << gmod << endl;
      
    // ISRF
    ISRF_model = GalpropISRF; //default value
    el1 = el->FirstChildElement("ISRF");
    if (el1) el1->QueryStringAttribute("type", &gmod);
    if      (gmod == "Galprop") ISRF_model = GalpropISRF;
    else if (gmod == "Uniform") ISRF_model = UniformISRF;
    else {
      cout << "implementing the standard Galprop-like ISRF " << endl;
      ISRF_model = GalpropISRF;
    }
    
    // SNR
    el1 = el->FirstChildElement("SNR");
    if (el1) el1->QueryStringAttribute("type", &gmod);
    else {
      cerr << "You have to specify a SNR distribution model" << endl;
      exit(NOSNR);
    }
    if (gmod == "Lorimer") SNR_model = Lorimer;
    else if (gmod == "Galprop") SNR_model = Galprop_;
    else if (gmod == "Ferriere") SNR_model = Ferriere;
    else if (gmod == "CaseBhattacharya") SNR_model = CaseBhattacharya;
    else if (gmod == "GiguereKaspi") SNR_model = GiguereKaspi;
    else if (gmod == "PointSource"){
      SNR_model = PointSource;
      
      pointsrc_x = QueryDoubleAttributeWithDefault("PointX", el1, 0);
      pointsrc_y = QueryDoubleAttributeWithDefault("PointY", el1, 0);
      pointsrc_z = QueryDoubleAttributeWithDefault("PointZ", el1, 0);

      cout <<" PointSource " << pointsrc_x << " " << pointsrc_y << " " << pointsrc_z << endl;


    }
    else if (gmod == "OneRing") {
      SNR_model = Ring;
      ringmin = QueryDoubleAttribute("RingMin", el1);
      ringmax = QueryDoubleAttribute("RingMax", el1);
      if (ringmin >= ringmax) {
	cerr << "Please reconsider your settings for the Ring SNR distribution. Ringmin = " << ringmin << " Ringmax = " << ringmax << endl;
	exit(PROBRING);
      }
    }
    else if (gmod == "Rings") {
      SNR_model = Rings;
      rings_period = QueryDoubleAttribute("RingsPeriod", el1);
      rings_phase = QueryDoubleAttribute("RingsPhase", el1);
    }
    else if (gmod == "Blasi") SNR_model = BlasiSmooth;
    else if (gmod == "OnlyExtra") //MW130828
      {
	if(prop_extracomp) SNR_model = OnlyExtra;
      }
    else {
      cerr << "Your SNR model is not implemented. You chose: " << gmod << endl;
      exit(SNRMODEL);
    }
    
    // SNR - extra component
    el1 = el->FirstChildElement("SNR_Extra");
    if (el1)
      {
	el1->QueryStringAttribute("type", &gmod);
	if (gmod == "Lorimer") SNR_model_Extra = Lorimer;
	else if (gmod == "Ferriere") SNR_model_Extra = Ferriere;
	else if (gmod == "CaseBhattacharya") SNR_model_Extra = CaseBhattacharya;
	else if (gmod == "GiguereKaspi") SNR_model_Extra = GiguereKaspi;
	else if (gmod == "PointSource"){
	  SNR_model_Extra = PointSource;
	  
	  pointsrc_Extra_x = QueryDoubleAttributeWithDefault("PointX", el1, 0);
	  pointsrc_Extra_y = QueryDoubleAttributeWithDefault("PointY", el1, 0);
	  pointsrc_Extra_z = QueryDoubleAttributeWithDefault("PointZ", el1, 0);

          cout <<"SNR_Extra" << " " << SNR_model_Extra << " " << pointsrc_Extra_x << " " << pointsrc_Extra_y << " " << pointsrc_Extra_z << endl;

	}
	else if (gmod == "Ring") {
	  SNR_model_Extra = Ring;
	  ringmin_extra = QueryDoubleAttribute("RingMin", el1);
	  ringmax_extra = QueryDoubleAttribute("RingMax", el1);
	  if (ringmin_extra >= ringmax_extra) {
	    cerr << "Please reconsider your settings for the Ring SNR distribution for the extra component. Ringmin = " << ringmin_extra << " Ringmax = " << ringmax_extra << endl;
	    exit(PROBRING);
	  }
	}
	else if (gmod == "Blasi") SNR_model_Extra = BlasiSmooth;
	else {
	  cerr << "Your SNR model for the extra component is not implemented. You chose: " << gmod << endl;
	}
      }
    else {
      if(SNR_model == OnlyExtra)
	{
	  cerr << "You have to specify a SNR distribution model. You chose OnlyExtra but then gave no extra component." << endl;
	  exit(NOSNR);
	}
    }
    
    // XCO
    el1 = el->FirstChildElement("XCOmode");
    if (el1) el1->QueryStringAttribute("type", &gmod);
    else {
      cerr << "You have to specify a XCO distribution model" << endl;
      exit(NOXCO);
    }
    if (gmod == "SM96") xco_mode = SM96; //Strong And Mattox 1996 -- constant value equal to 1.9 in standard units
    else if (gmod == "galprop_2004") xco_mode = galprop_2004;
    else if (gmod == "galprop_2010") xco_mode = galprop_2010;
    else if (gmod == "constant")  {
      xco_mode = constant;
      xco_constant = QueryDoubleAttribute("XCO_constant", el1);
    }
    else if (gmod == "dragon") {
      xco_mode = dragon;
      xco_inner = QueryDoubleAttribute("XCO_inner", el1);
      xco_outer = QueryDoubleAttribute("XCO_outer", el1);
    }
    else {
      cerr << "Your XCO model is not implemented. You chose: " << gmod << endl;
      exit(XCOMODEL);
    }
      
    // Diffusion
    el1 = el->FirstChildElement("Diffusion");
      
    if (el1)
      {
	DiffT = Isotropic;
	el2 = el1->FirstChildElement("AnisotropicDiffusion");
	if (el2) {
	  DiffT = Anisotropic;
	  Dpar = QueryDoubleAttribute("DPar_1e28", el2)*1.e28/kpc/kpc*Myr;
	  etaTpar = QueryDoubleAttribute("etaTPar", el2);
	  DeltaPar = QueryDoubleAttribute("DeltaPar", el2);
	  Dperp = QueryDoubleAttribute("DPerp_1e28", el2)*1.e28/kpc/kpc*Myr;
	  etaTperp = QueryDoubleAttribute("etaTPerp", el2);
	  DeltaPerp = QueryDoubleAttribute("DeltaPerp", el2);
	}
	
	el1->QueryStringAttribute("type", &gmod);
	if (gmod == "Constant") set_profile = Constant;
	else if (gmod == "Exp") set_profile = Exp;
	else if (gmod == "Blasi") set_profile = Blasi;
	else if (gmod == "Expr") set_profile = Expr;
	else if (gmod == "Qtau") set_profile = Qtau;
	else if (gmod == "Test") set_profile = Test;
	else if (gmod == "ExpRadial") set_profile = ExpRadial;
	else {
	  cerr << "Your Diffusion model is not implemented. You chose: " << gmod << endl;
	  exit(DIFFMODEL);
	}
         
	D0 = QueryDoubleAttribute("D0_1e28", el1)*1.e28/kpc/kpc*Myr;
	D_ref_rig = QueryDoubleAttribute("DiffRefRig", el1);       
	delta = QueryDoubleAttribute("Delta", el1);
	delta_h = QueryDoubleAttribute("Delta_h", el1);
	if(delta_h==-1) delta_h = delta; //default value

	VariableDelta = QueryIntAttributeWithDefault("VariableDelta",el1,0);
	diff_threshold = QueryDoubleAttributeWithDefault("DiffusionThreshold",el1,11.);
	if (VariableDelta == 1) {
	  delta_A = QueryDoubleAttributeWithDefault("deltaA",el1,0.);		
	  delta_B = QueryDoubleAttributeWithDefault("deltaB",el1,delta);		
	}

	rho_b = QueryDoubleAttribute("rho_b", el1);
	if(rho_b==-1) rho_b = 1000.; //default value
	zt = QueryDoubleAttribute("zt", el1);
	etaT = QueryDoubleAttribute("etaT", el1);
	index_radial = QueryDoubleAttribute("index_radial", el1);
	if(index_radial==-1) index_radial = 0; //default value
      }
    else {
      cerr << "No Diffusion " << gmod << endl;
      exit(DIFFMODEL);
    }
      
    // Reacceleration
    el1 = el->FirstChildElement("Reacceleration");
    if (!el1) REACC = false;
    else {
      REACC = true;
      el1->QueryStringAttribute("type", &gmod);
      if (gmod == "Ptuskin94") diff_reacc = Ptuskin94;
      else if (gmod == "Ptuskin03") diff_reacc = Ptuskin03;
      else {
	cerr << "Your Reacceleration model is not implemented. You chose: " << gmod << endl;
	exit(REACCMODEL);
      }
      //            cout << QueryDoubleAttribute("vA_kms", el1) << endl;
      vAlfven = QueryDoubleAttribute("vA_kms", el1)*km*Myr/kpc;
    }
      
    // Convection
    el1 = el->FirstChildElement("Convection");
    if (!el1) CONVECTION = false;
    else {
      //f_b = 0;//default value
      //z_k =0.1;//default value
      conv_index_radial = 0; //default value
      conv_threshold = 0;//default value
      CONVECTION = true;  	    //
      el1->QueryStringAttribute("type", &gmod);
      if (gmod == "Constant")    set_profile_conv = Constant;
      else if (gmod == "Radial") set_profile_conv = Radial;
      else set_profile_conv = Constant;
      //
      v0 = QueryDoubleAttribute("v0_kms", el1)*km*Myr/kpc;
      //f_b = QueryDoubleAttribute("f_b", el1);
      //if(f_b == -1) f_b = 0; //default value
      //vb=v0*f_b;
      //
      //DG28.11.2013 solved little bug pointed out by Simon	
      f_b = QueryDoubleAttributeWithDefault("f_b", el1, 1.);
      vb=v0*f_b;
      //z_k = QueryDoubleAttribute("z_k_kpc", el1);
      //if(z_k == -1) z_k =0.1; //default value
      z_k = QueryDoubleAttributeWithDefault("z_k_kpc", el1, 0.1); 
      //     
      dvdz = QueryDoubleAttribute("dvdz_kmskpc", el1)*km*Myr/kpc;
      conv_index_radial = QueryDoubleAttribute("conv_index_radial",el1);
      if(conv_index_radial == -1) conv_index_radial = 0; //default value
      conv_threshold = QueryDoubleAttribute("conv_threshold",el1);
      if(conv_threshold == -1) conv_threshold = 0; //default value
    }
    if (feedback >0) cout << "Convection settings read " << endl;
      
    // Magnetic field
    string bf;
    el1 = el->FirstChildElement("MagneticField");
    if (!el1) cerr << "You did not specify a magnetic field model. A null magnetic field will be used." << endl;
    else {
      el1->QueryStringAttribute("type", &bf);
      if (bf == "Pshirkov") BM = Pshirkov;
      //else if (bf == "Azimuthal") BM = Azimuthal;
      //else if (bf == "Spiral") BM = Spiral;
      //else if (bf == "Prouda") BM = Prouda;
      //else if (bf == "Wmap") BM = Wmap;
      else if (bf == "Farrar") BM = Farrar;
      else if (bf == "Uniform") BM = Uniform;
      else if (bf == "Simple") BM = Simple;
      //else if (bf == "Tkachev") BM = Tkachev;
      else if (bf == "ToyModel") { BM = ToyModel; //useful to test anisotropic diffusion
	cout << "ToyModel mag field selected in input file" << endl;
      }
      else {
	cerr << "Your specified Magnetic Field model " << bf << " has not been implemented yet." << endl;
	exit(NO_BMODEL);
      }
         
      B0disk = QueryDoubleAttribute("B0disk", el1);
      B0halo = QueryDoubleAttribute("B0halo", el1);
      B0turb = QueryDoubleAttribute("B0turb", el1);
      if (feedback >0) cout << "B0 disk/halo/turb = " << B0disk << " " << B0halo << " " << B0turb << endl;	
      betaFarrar = QueryDoubleAttribute("betaFarrar", el1);
      b0 = QueryDoubleAttribute("b0", el1);
      if(b0 == -1) b0=6.1e-6; //default value
      if (bf == "ToyModel") {
	bx = QueryDoubleAttributeWithDefault("bx", el1, 1.e-6);
	by = QueryDoubleAttributeWithDefault("by", el1, 0.);
	bz = QueryDoubleAttributeWithDefault("bz", el1, 0.);
	bturb = QueryDoubleAttributeWithDefault("bturb", el1, 0.);
      }
    }
    
    el1 = el->FirstChildElement("CrossSection");
    if (el1) {
      
      el1->QueryStringAttribute("type", &gmod);
      if (gmod == "GalpropXSec")  spallationxsec = GalpropXSec;
      else if (gmod == "Webber03") spallationxsec = Webber03;
      else if (gmod == "DRAGON2") spallationxsec = DRAGON2;
      else if (gmod == "Fluka") {
	spallationxsec = Fluka;
        // cout << "Fluka E grid " << Ekmin << " " << Ekmax << endl;
	// Ekfact = pow(10., 1./16.);
	// if(Ekmin <= 0.01) {
	//   Ekmin = 0.01;
	// } else {
	//   int iemin = int(log10(Ekmin/0.01)*16.);
	//   Ekmin = 0.01*pow(10., float(iemin)/16.);
	// }
	// if(Ekmax >= 0.01*pow(Ekfact, 160.)){
	//   Ekmax = 0.01*pow(Ekfact, 160.);
	// } else {
	//   int iemax = int(log10(Ekmax/0.01)*16.)+1;
	//   if(iemax>160) {
	//     iemax = 160;
	//   }
	//   cout << "iemax = " << iemax << endl;
	//   Ekmax = 0.01*pow(10., float(iemax)/16.);
	// }
        // cout << "Fluka E grid " << Ekmin << " " << Ekmax << endl;
      }
      else {
	cerr << "Your cross section model has not been implemented yet. You chose " << gmod << endl;
	//exit(XSECMODEL);
	cerr << "Setting the default XSec model: DRAGON2 algorithm" << endl;
	spallationxsec = DRAGON2;
      }
         
      el1->QueryStringAttribute("leptopt", &gmod);
      if (gmod == "GalpropTable") ly = GalpropTable;
      else if (gmod == "Kamae") ly = Kamae;
      else if (gmod == "Pohl") ly = Pohl;
      else if (gmod == "Fluka") ly = FlukaLep;
      else {
	if (feedback >0) cout << "You did not choose any lepton yield. It will be set to Kamae." << endl;
	ly = Kamae;
      }
         
      el1->QueryIntAttribute("ApCs", &antiproton_cs);
      //antiproton_cs = QueryIntAttribute("ApCs", el);
      if ( antiproton_cs != 2 ){ 
	cout<<"ApCs value not valid!"<<endl;
	exit(111);
      }
      
      el1->QueryStringAttribute("apopt", &gmod);
      if (gmod == "GalpropFunction") apy = GalpropFunction;
      else if (gmod == "QGSJET") apy = QGSJET; //  M.Kachelriess, S.Ostapchenko, submitted to PRD
      else if (gmod == "Fluka") apy = FlukaAp;
      else if (gmod == "Winkler") apy = Winkler;
      else {
	if (feedback >0) cout << "You did not choose any antiproton yield. It will be set to Winkler." << endl;
	apy = Winkler;
      }
    }
    else {
      cerr << "You have to choose a Cross section model." << endl;
      exit(NOXSEC);
    }
      
    // Matthias: factors for Local Bubble and Spiral Arms
    el1 = el->FirstChildElement("LocalBubble");
    LB_diff = LB_delta = LB_convec = LB_vA = LB_source = LB_DMsource = LB_gastotal = LB_gasH2 = LB_gasHI = LB_gasHII = LB_MagField = LB_ISRF = 1;
    if(el1)
      {
	LB_diff     = QueryDoubleAttributeWithDefault("FactorDiffCoeff"    , el1, 1);
	LB_delta    = QueryDoubleAttributeWithDefault("FactorDiffDelta"    , el1, 1); //MW: not implemented yet!
	LB_convec   = QueryDoubleAttributeWithDefault("FactorConvection"   , el1, 1);
	LB_vA       = QueryDoubleAttributeWithDefault("FactorReaccCoeff"   , el1, 1);
	LB_source   = QueryDoubleAttributeWithDefault("FactorSource"       , el1, 1);
	LB_DMsource = QueryDoubleAttributeWithDefault("FactorDMSource"     , el1, 1);
	LB_gastotal = QueryDoubleAttributeWithDefault("FactorGasTotal"     , el1, 1);
	LB_gasH2    = QueryDoubleAttributeWithDefault("FactorGasH2"        , el1, 1);
	LB_gasHI    = QueryDoubleAttributeWithDefault("FactorGasHI"        , el1, 1);
	LB_gasHII   = QueryDoubleAttributeWithDefault("FactorGasHII"       , el1, 1);
	LB_MagField = QueryDoubleAttributeWithDefault("FactorMagneticField", el1, 1);
	LB_ISRF     = QueryDoubleAttributeWithDefault("FactorISRF"         , el1, 1);
      }
      
    SA_type = "None";
    el1 = el->FirstChildElement("SpiralArms");
    SA_diff = SA_convec = SA_vA = SA_source = SA_DMsource = SA_gasH2 = SA_gasHI = SA_gasHII = SA_MagField = SA_ISRFStar = SA_ISRFDust = 0;
    SA_cut_diff = SA_cut_convec = SA_cut_vA = 100;
    SA_cut_source = SA_cut_DMsource = SA_cut_gasH2 = SA_cut_gasHI = SA_cut_gasHII = SA_cut_MagField = SA_cut_ISRF = 1e16;
    //MW130625: have to separate ISRF contribution: Starlight and Dust
      
    if(el1)
      {
	el1->QueryStringAttribute("type", &SA_type); //"CII" or "NII", or "None"
         
	//MW130623: Trying to coin these powers "Kappa"
	SA_diff     = QueryDoubleAttributeWithDefault("KappaDiffCoeff"    , el1, 0);
	SA_convec   = QueryDoubleAttributeWithDefault("KappaConvection"   , el1, 0);
	SA_vA       = QueryDoubleAttributeWithDefault("KappaReaccCoeff"   , el1, 0);
	SA_source   = QueryDoubleAttributeWithDefault("KappaSource"       , el1, 0);
	SA_DMsource = QueryDoubleAttributeWithDefault("KappaDMSource"     , el1, 0);
	SA_gasH2    = QueryDoubleAttributeWithDefault("KappaGasH2"        , el1, 0);
	SA_gasHI    = QueryDoubleAttributeWithDefault("KappaGasHI"        , el1, 0);
	SA_gasHII   = QueryDoubleAttributeWithDefault("KappaGasHII"       , el1, 0);
	SA_MagField = QueryDoubleAttributeWithDefault("KappaMagneticField", el1, 0);
	SA_ISRFStar = QueryDoubleAttributeWithDefault("KappaISRFStarlight", el1, 0);
	SA_ISRFDust = QueryDoubleAttributeWithDefault("KappaISRFDust"     , el1, 0);
         
	SA_cut_diff     = QueryDoubleAttributeWithDefault("CutDiffCoeff"    , el1, 100);
	SA_cut_convec   = QueryDoubleAttributeWithDefault("CutConvection"   , el1, 100);
	SA_cut_vA       = QueryDoubleAttributeWithDefault("CutReaccCoeff"   , el1, 100);
	SA_cut_source   = QueryDoubleAttributeWithDefault("CutSource"       , el1, 1e16);
	SA_cut_DMsource = QueryDoubleAttributeWithDefault("CutDMSource"     , el1, 1e16);
	SA_cut_gasH2    = QueryDoubleAttributeWithDefault("CutGasH2"        , el1, 1e16);
	SA_cut_gasHI    = QueryDoubleAttributeWithDefault("CutGasHI"        , el1, 1e16);
	SA_cut_gasHII   = QueryDoubleAttributeWithDefault("CutGasHII"       , el1, 1e16);
	SA_cut_MagField = QueryDoubleAttributeWithDefault("CutMagneticField", el1, 1e16);
	SA_cut_ISRF     = QueryDoubleAttributeWithDefault("CutISRF"         , el1, 1e16);
         
	if (SA_type == "Daniele" || SA_type == "Blasi" || SA_type == "BlasiModel") {
            
	  num_arms  = QueryIntAttribute("NumArms", el1);
	  spiral_width = QueryDoubleAttributeWithDefault("SpiralWidth", el1, 0.4);
	  cout << "You specified a spiral arm structure! " << endl;
            
	  for (int i_arms = 0; i_arms < num_arms; i_arms++)
            {
	      ostringstream Kstring, r0string, theta0string;
	      Kstring  << "K_"   << i_arms+1;
	      r0string << "r0_"  << i_arms+1;
	      theta0string << "theta0_" << i_arms+1;
	      arms_Kvec.push_back( QueryDoubleAttribute(Kstring.str().c_str(), el1)  );
	      arms_r0vec.push_back( QueryDoubleAttribute(r0string.str().c_str(), el1)  );
	      arms_theta0vec.push_back( QueryDoubleAttribute(theta0string.str().c_str(), el1)  );
	      cout << "Arm number " << i_arms << " K = " << arms_Kvec.back() << " r0 = " << arms_r0vec.back() << " theta0 = " << arms_theta0vec.back() << endl;
            }
	}
      }
      
    // Flags fk 130701
    el = doc->FirstChildElement("Flags");
      
    DIFFUSION  = (el) ? QueryIntAttribute("diff_flag", el) : true;
    CONVECTION = (el) ? QueryIntAttribute("conv_flag", el) : CONVECTION;
    REACC      = (el) ? QueryIntAttribute("reacc_flag", el) : REACC;
    ELOSS      = (el) ? QueryIntAttribute("eloss_flag", el) : true;
    RDECAY     = (el) ? QueryIntAttribute("decay_flag", el) : true;
    SPALL      = (el) ? QueryIntAttribute("spall_flag", el) : true;
    if(!DIFFUSION)
      {
	D0 = 0.; // for D0=0 there is no diffusion
	REACC = 0; // diffusive Reacceleration makes no sense without Diffusion
      }
      
    // Cosmic Rays
    el = doc->FirstChildElement("CR");
    if (!el) {
      cerr << "You must set CR properties." << endl;
      exit(NOCR);
    }
    else {
      el->QueryStringAttribute("type", &gmod);
      if (gmod == "PowerLawBreak") Spectrum = PowerLawBreak;
      else if (gmod == "ExpSuppressed") Spectrum = ExpSuppressed;
      else {
	if (feedback >0) cout << "You did not choose any spectral shape. It will be set to PowerLawBreak." << endl;
	Spectrum = PowerLawBreak;
      }
      
      // DG30.09.2013 NEW implementation of injection slopes for the nuclei and electrons with arbitrary number of breaks
      // Nuclear slopes and break positions (arbitrary number of breaks)
      UseInjectionIndexAllNuclei = false;
      if (prop_nuclei==true) {
	el1 = el->FirstChildElement("InjectionIndexAllNuclei");
	if (!el1){
	  if (feedback >0) cout << "You did not specify an InjectionIndexAllNuclei, using .source.param values." << endl;
	}
	else if (el1) {
	  UseInjectionIndexAllNuclei = true;
	  cout << "Using the same injection index for all species from xml file " << endl;
	  inp_break_positions = ReadInjectionStructure("rho_", el1);
	  inp_inj_indexes = ReadInjectionStructure("alpha_", el1);
	  cutoff_rig = QueryDoubleAttribute("CutoffRig",el1);
	}
      }
      
      //Electron slopes and break positions (arbitrary number of breaks)
      if (prop_lep==true) {
	el1 = el->FirstChildElement("InjectionIndexElectrons");
	if (!el1) {
	  cerr << "You must set electron injection properties." << endl;
	  exit(NOCR);
	}
	else {
	  inp_break_el_positions = ReadInjectionStructure("rho_", el1);
	  inp_inj_el_indexes = ReadInjectionStructure("alpha_", el1);
	  cout << "+++ Electron indexes +++" << endl;
	  for (int i=0; i<inp_inj_el_indexes.size(); i++)
	   cout << inp_inj_el_indexes[i] << " "<<endl;
          cout << endl;
	  cutoff_rig_el = QueryDoubleAttributeWithDefault("CutoffRigEl",el1, 1.e12);
	  cout << "electron cutoff at " << cutoff_rig_el << endl;
	}
      }
          
      //ExtraComponent slopes and break positions (arbitrary number of breaks)
      if (prop_extracomp==true) {
	el1 = el->FirstChildElement("InjectionIndexExtraComponent");
	if (!el1) {
	  cerr << "You must set extra component injection properties." << endl;
	  exit(NOCR);
	}
	else {
	  inp_break_extra_positions = ReadInjectionStructure("rho_", el1);
	  inp_inj_extra_indexes = ReadInjectionStructure("alpha_", el1);
	  cout << "+++ Extra component indexes +++" << endl;
	  for (int i=0; i<inp_inj_extra_indexes.size(); i++)
	   cout << inp_inj_extra_indexes[i] << " "<<endl;
          cout << endl;
	  cutoff_rig_extra = QueryDoubleAttributeWithDefault("CutoffRigExtra",el1, 1.e12);
          cout << "extra component cutoff at " << cutoff_rig_extra << endl;
	}
      }
      
      //Normalization energies and corresponding fluxes
      sp_ref_rig_norm = QueryDoubleAttribute("ProtNormEn_GeV", el);
      sp_ref_rig_el   = QueryDoubleAttribute("ElNormEn_GeV", el);
      sp_ref_rig_el_extra = QueryDoubleAttribute("ElNormEnExtra_GeV", el);
      cout << "norm energy of the extra component = " << sp_ref_rig_el_extra << endl;
      spect_norm    = QueryDoubleAttribute("ProtNormFlux", el);
      spect_norm_el = QueryDoubleAttribute("ElNormFlux", el);
      spect_norm_el_extra = QueryDoubleAttribute("ElNormFluxExtra", el);
      cout << "norm energy of the extra component = " << spect_norm_el_extra << endl;
         
      //Add the possibility to put a exp cutoff again and reorder this stuff (DG) -----------------------
      //sp_ref_rig = QueryDoubleAttribute("NucRigBreak", el);
      //sp_ref_rig_exp_el = QueryDoubleAttribute("ElRigExp", el);;
      //sp_ref_rig_exp = QueryDoubleAttribute("NucRigExp", el);
      //exp_cut_index_el = QueryDoubleAttribute("ExpCutIndexEl", el);
      //cutoff_rig = QueryDoubleAttribute("ElCutoffRig", el);
      //spect_ind_el_low = QueryDoubleAttribute("ElSpect_index_below_Rig", el);
      //spect_ind_el = QueryDoubleAttribute("ElSpect_index_above_Rig", el);
      //sp_ref_rig_break_el = QueryDoubleAttribute("ElBreakRig", el);
      //-------------------------------------------------------------------------------------------------
         
      //This is for the extra component
      //spect_ind_el_low_extra = QueryDoubleAttribute("ElSpect_index_below_RigExtra", el);
      //	spect_ind_el_extra = QueryDoubleAttribute("ElSpect_index_above_RigExtra", el);
      //	sp_ref_rig_break_el_extra = QueryDoubleAttribute("ElBreakRigExtra", el);
      //	cutoff_rig_extra = QueryDoubleAttribute("ElCutoffRigExtra", el);

      sourcedata = string(inputfilename, 0, inputfilename.size()-4);
      sourcedata += ".source.param";
      //cout << "Check sourcedata: " << sourcedata << endl;
      
      //antiproton_cs = QueryIntAttribute("ApCs", el);
      //scaling = ( el->FirstChildElement("Scaling") );
    }
    
    el = doc->FirstChildElement("DarkMatter");
    //if (!el) cerr << "You did not set DM properties." << endl;
    if (el) {
      el1 = el->FirstChildElement("MovingClump");
         
      if (el1) MOVING_CLUMP = true;
      else MOVING_CLUMP = false;
      
      if (MOVING_CLUMP) {
	clump_x0 = QueryDoubleAttribute("clump_x0", el1);
	clump_y0 = QueryDoubleAttribute("clump_y0", el1);
	clump_z0 = QueryDoubleAttribute("clump_z0", el1);
	clump_vx = QueryDoubleAttribute("clump_vx", el1);
	clump_vy = QueryDoubleAttribute("clump_vy", el1);
	clump_vz = QueryDoubleAttribute("clump_vz", el1);
	clump_norm = QueryDoubleAttribute("clump_norm", el1);
	clump_size = QueryDoubleAttributeWithDefault("clump_size", el1, 0.1);
	clump_inj = QueryDoubleAttribute("clump_inj", el1);
	clump_cutoff = QueryDoubleAttribute("clump_cutoff", el1);
	clump_deltat = QueryDoubleAttribute("clump_deltat",el1);
	analytical_refinement = ( el1->FirstChildElement("AnalyticalRefinement") );
	stop_after_timestep = QueryDoubleAttributeWithDefault("StopAfterTimestep",el1,-1);
            
	cout << "I read all the MOVING CLUMP parameters " << endl;
	//if (gmod == "Delta") DMs = Delta;
      }
         
      propDM = true;
      prop_DMap = ( el->FirstChildElement("PropDMAntiProton") );
      prop_DMdeuteron = ( el->FirstChildElement("PropDMAntiDeuteron") );
      prop_DMel = ( el->FirstChildElement("PropDMLepton") );
         
      mx = QueryDoubleAttribute("Mass", el);
      taudec = QueryDoubleAttribute("LifeTime", el);
      sigmav = QueryDoubleAttribute("SigmaV", el);
      rhos = QueryDoubleAttribute("SSDensity", el);
      EkDelta = QueryDoubleAttribute("EkDelta", el);
      dmmode = QueryIntAttribute("Channel", el);
         
      el->QueryStringAttribute("Profile", &gmod);
      if (gmod == "Iso") dmprof = ISO;
      else if (gmod == "NFW") dmprof = NFW;
      else if (gmod == "Kra") dmprof = Kra;
      else if (gmod == "Moore") dmprof = Moore;
      else if (gmod == "Einasto") dmprof = Einasto;
      else {
	cerr << "Your choice of DM density profile was not implemented. You chose: " << gmod << endl;
	exit(NODMPROFILE);
      }
      el->QueryStringAttribute("Reaction", &gmod);
      if (gmod == "Annihilation") DMr = Annihilation;
      else if (gmod == "Decay") DMr = Decay;
      else {
	cerr << "Your choice of DM reaction was not implemented. You chose: " << gmod << endl;
	exit(NODMREACTION);
      }
         
      el->QueryStringAttribute("Model", &gmod);
         
      if (gmod == "Delta") {DMs = Delta; cout << "delta injection! " <<endl;}
      else if (gmod == "DarkSUSY") DMs = DarkSUSY;
      else if (gmod == "EWCorrections") {
#ifdef HAVE_ROOT
            
	DMs = EWCorrections;
#else
	cerr << "You must compile with ROOT support to use EW Corrections." << endl;
	exit(NO_ROOT_EWCORRECTIONS);
#endif
      }
      else if (gmod == "SelfTable") {
	DMs = SelfTable;
	if (prop_DMap){cout<<"Using user defined table for DM ap"<<endl;
	  //el->QueryStringAttribute("AntiprotonDatafile", &MySelfTableDMap);cout<<MySelfTableDMap.c_str()<<endl;}
	  MySelfTableDMap = QueryStringAttribute("AntiprotonDatafile", el);cout<<MySelfTableDMap.c_str()<<endl;}
	if (prop_DMdeuteron)
	  //el->QueryStringAttribute("AntideuteronDatafile", &MySelfTableDMdbar);
	  MySelfTableDMdbar = QueryStringAttribute("AntideuteronDatafile", el);
	if (prop_DMel){cout<<"Using user defined table for DM leptons"<<endl;
	  //el->QueryStringAttribute("LeptonDatafile", &MySelfTableDMel);cout<<MySelfTableDMel.c_str()<<endl;}
	  MySelfTableDMel = QueryStringAttribute("LeptonDatafile", el);cout<<MySelfTableDMel.c_str()<<endl;}
            
      }
      else {
	cerr << "Your choice of DM model was not implemented. You chose: " << gmod << endl;
	exit(NODMMODEL);
      }
         
    }
    else {
      propDM = false;
      prop_DMap = false;
      prop_DMdeuteron = false;
      prop_DMel = false;
    }
      
    prop_ap = (prop_secap || prop_DMap);
    prop_deuteron = (prop_secdeuteron || prop_DMdeuteron);
    //prop_el = (prop_secel || prop_DMel || prop_TPP);
      
    // Simon: DMCMC ************************************
    el = doc->FirstChildElement("DMCMC");
    write_flag=0;
    fit_mod=0;
    if (el){
         
      el1 = el->FirstChildElement("run_id");
      if (el1) {
	el1->QueryStringAttribute("type", &gmod);
	run_id = gmod;
      }
         
      write_flag = QueryIntAttribute("write_flag", el);
      fit_mod = QueryIntAttribute("fit_mod", el);
    }
    // end DMCMC ****************************************    
  }
  else {
    cerr << "Problems with the input file " << inputfilename << endl;
    exit(PROBINPUT);
  }
   
  delete doc;
   
  cout << "Input file read successfully!" << endl;
   
  return 0;
}

void Input::Print() {
   
  cout << endl;
  cout << "*******************************************************" << endl;
  cout << "*******************************************************" << endl;
   
  cout << "Settings are " << endl;
  cout << "Delta E / E = " << Ekfact << endl;
  cout << "Ekmax = " << Ekmax << endl;
  cout << "Ekmin = " << Ekmin << endl;
  cout << "Dimension of R grid = " << numr << endl;
  cout << "Dimension of X grid = " << numx << endl;
  cout << "Dimension of Y grid = " << numy << endl;
  cout << "Dimension of Z grid = " << numz << endl;
  cout << "Rmin = " <<  Rmin << endl;
  cout << "Rmax = " <<  Rmax << endl;
  cout << "Zmax = " <<  zmax << endl;
  cout << "Grid type = " << gridtype << endl;
   
  // Limits for the nuclear chain
  cout << "Maximal propagated charge = " << Zmax << endl;
  cout << "Minimal propagated charge = " << Zmin << endl;
  cout << "Whether to propagate antiprotons " << prop_ap << endl;
  cout << "Whether to propagate antideuterons " << prop_deuteron << endl;
  cout << "Whether to propagate electrons and positrons " <<  prop_lep << endl;
  cout << "Whether to include K-capture " << Kcapture << endl;
   
   
  // I/O
  cout << "Fullstore " << fullstore << endl;
  cout << "Partial store " << partialstore << endl;
  cout << "Run name " << filename << endl;
   
  // Algorithm
  cout << "Nrept = " << Nrept << endl;
  cout << "Dt factor = " << dtfactor << endl;
  cout << "Dt min = " << dtmin << endl;
  cout << "Dt max = " << dtmax << endl;
   
  cout << "ADI? " << ADI << endl;
  cout << "OpSplit? " << OpSplit << endl;
  cout << "Tol = " << tol << endl;
  cout << "Alpha = " << alpha << endl;
   
  // Galaxy
  cout << "GasType = " << gas_model << endl;
  cout << "SNRType = " << SNR_model << endl;
  cout << "xco_modes = " << xco_mode << endl;
  cout << "Ringmin = " << ringmin << endl;
  cout << "Ringmax = " << ringmax << endl;
   
   
  //Flags fk 130701
  cout << "Transport equation with:" << endl;
  cout << "DIFF? " << DIFFUSION << endl;
  cout << "CONV? " << CONVECTION << endl;
  cout << "REACC? " << REACC << endl;
  cout << "ELOSS? " << ELOSS << endl;
  cout << "RDECAY? " << RDECAY << endl;
  cout << "SPALL? " << SPALL << endl;
   
   
  // Reacceleration
  cout << "vA = " << vAlfven << endl;
  cout << "ReaccType = " << diff_reacc << endl;
   
  // Convection
  cout << "v0 = " << v0 << endl;
  cout << "dvdz = " << dvdz << endl;
   
  // Diffusion
  cout << "DPerpType = " << set_profile << endl;
  cout << "D0 = " << D0 << endl;
  cout << "zt = " << zt << endl;
  cout << "etaT = " << etaT << endl;
  cout << "Delta = " << delta << endl;
  cout << "index_radial = " << index_radial << endl;

  if (VariableDelta == 1) {
    cout << "deltaA = " << delta_A << endl;		
    cout << "deltaB = " << delta_B << endl;
  }

  cout << "DPar = " << Dpar << endl;
  cout << "etaTpar = " << etaTpar << endl;
  cout << "DeltaPar = " << DeltaPar << endl;
   
  cout << "DPerp = " << Dperp << endl;
  cout << "etaTperp = " << etaTperp << endl;
  cout << "DeltaPerp = " << DeltaPerp << endl;
   
  // Cross sections
  cout << "SpallXSec = " << spallationxsec << endl;
  cout << "LeptonYields = " <<ly << endl;
  cout << "AntiprotonYiels = " <<apy <<endl;
   
  // CR normalization
  cout << "sp_ref_rig_norm = " << sp_ref_rig_norm << endl;
  cout << "sp_ref_rig_el = " << sp_ref_rig_el << endl;
  cout << "spect_norm = " << spect_norm << endl;
  cout << "spect_norm_el = " << spect_norm_el << endl;
   
  /*cout << "Nuclear injection slopes " << endl;
    cout << alpha_0 << " up to " << rho_0 << " GeV, then " << endl;
    cout << alpha_1 << " up to " << rho_1 << " GeV, then " << endl;
    cout << alpha_2 << " up to " << rho_2 << " GeV, then " << endl;
    cout << alpha_3 << " up to Emax" << endl;
  */
   
  /*cout << "Electron injection slopes " << endl;
    cout << alpha_el_0 << " up to " << rho_el_0 << " GeV, then " << endl;
    cout << alpha_el_1 << " up to " << rho_el_1 << " GeV, then " << endl;
    cout << alpha_el_2 << " up to " << rho_el_2 << " GeV, then " << endl;
    cout << alpha_el_3 << " up to Emax" << endl;*/
   
  // Parameter space
  //cout << "sp_ref_rig = " << sp_ref_rig << endl;
  // Leptons
  //cout << "cutoff_rig = " << cutoff_rig << endl;
  //cout << "spect_ind_el_low = " << spect_ind_el_low << endl;
  //cout << "spect_ind_el = " << spect_ind_el << endl;
  //cout << "sp_ref_rig_break_el = " << sp_ref_rig_break_el << endl;
   
  // Antiprotons
  cout << "antiproton_cs = " << antiproton_cs << endl;
  cout << "Scaling for ap cross sections = " << scaling << endl;
   
  // DM
  cout << "Propagate DM? " << propDM << endl;
  cout << "Mass = " << mx << endl;
  cout << "taudec = " << taudec << endl;
  cout << "sigmav = " << sigmav << endl;
  cout << "DMreaction = " << DMr << endl;
  cout << "DMspectrum = " << DMs << endl;
  cout << "rhos = " << rhos << endl;
  cout << "E Delta = " << EkDelta << endl;
  //cout << "MySelfTable for antiprotons = " << MySelfTableDMap << endl;
  //cout << "MySelfTable for leptons     = " << MySelfTableDMel << endl;
   
  cout << "Dm mode = " << dmmode << endl;
   
  cout << "DMprofile = " << dmprof << endl;
   
  cout << endl;
  cout << "*******************************************************" << endl;
  cout << "*******************************************************" << endl;
   
  return ;
}

/*
  Input::Input(char** argv, int argc) :
  Ekfact(atof(argv[2])),
  #ifndef CLUMPS
  numr((atoi(argv[3])%2 == 1) ? atoi(argv[3]): atoi(argv[3])+1),
  numz((atoi(argv[4])%2 == 1) ? atoi(argv[4]): atoi(argv[4])+1),
  D0(atof(argv[5])*1.e28/kpc/kpc*Myr),
  zt(atof(argv[6])),
  zmax(atof(argv[7])),
  delta(atof(argv[8])),
  index_radial(atof(argv[9])),
  #ifdef HAVE_DS
  mx(atof(argv[10])),
  taudec(atof(argv[11])),
  #else
  ab_C(atof(argv[10])),
  ab_N(atof(argv[11])),
  #endif
  #ifdef REAC
  vAlfven(atof(argv[12])*km*Myr/kpc),
  #endif
  v0(atof(argv[13])*km*Myr/kpc),
  dvdz(atof(argv[14])*km*Myr/kpc),
  etaT(atof(argv[15])),
  #ifdef HAVE_DS
  sigmav(atof(argv[16])),
  dmmode(atoi(argv[17])),
  dmprof(DMprofile(atoi(argv[18]))),
  #endif
  #else
  numx((atoi(argv[3])%2 == 1) ? atoi(argv[3]): atoi(argv[3])+1),
  numy((atoi(argv[4])%2 == 1) ? atoi(argv[4]): atoi(argv[4])+1),
  numz((atoi(argv[5])%2 == 1) ? atoi(argv[5]): atoi(argv[5])+1),
  D0(atof(argv[6])*1.e28/kpc/kpc*Myr),
  zt(atof(argv[7])),
  zmax(atof(argv[8])),
  delta(atof(argv[9])),
  index_radial(atof(argv[10])),
  mx(atof(argv[11])),
  taudec(atof(argv[12])),
  #ifdef REAC
  vAlfven(atof(argv[13])*km*Myr/kpc),
  #endif
  v0(atof(argv[14])*km*Myr/kpc),
  dvdz(atof(argv[15])*km*Myr/kpc),
  etaT(atof(argv[16])),
  sigmav(atof(argv[17])),
  dmmode(atoi(argv[18])),
  dmprof(DMprofile(atoi(argv[19]))),
  #endif
  filename(argv[1]) {
 
  if (!strcmp(argv[argc-1], "const")) set_profile = Constant;
  else if (!strcmp(argv[argc-1], "exp")) set_profile = Exp;
  else if (!strcmp(argv[argc-1], "expr")) set_profile = Expr;
  else if (!strcmp(argv[argc-1], "radial")) set_profile = Radial;
  else if (!strcmp(argv[argc-1], "exp_radial")) set_profile = ExpRadial;
  else if (!strcmp(argv[argc-1], "blasi")) set_profile = Blasi;
  else {
  cerr << "What DPerpType did you choose?" << endl;
  exit(-2);
  }
  }
 
  Input::Input(char* run, double& Ekfact_, int& Nr, int& Nz, double& D0_, double& zt_, double& zmax_, double& delta_, double& index_radial_, double& ab_C_, double& ab_N_, double& vA_, double& vC_, double& dvdz_conv_, double etaT_, char* profile) :
  Ekfact(Ekfact_),
  numr((Nr%2 == 1) ? Nr: Nr+1),
  numz((Nz%2 == 1) ? Nz: Nz+1),
  D0(D0_*1.e28/kpc/kpc*Myr),
  zt(zt_),
  zmax(zmax_),
  delta(delta_),
  index_radial(index_radial_),
  #ifndef HAVE_DS
  ab_C(ab_C_),
  ab_N(ab_N_),
  #endif
  #ifdef REAC
  vAlfven(vA_*km*Myr/kpc),
  #endif
  v0(vC_*km*Myr/kpc),
  dvdz(dvdz_conv_*km*Myr/kpc),
  etaT(etaT_),
  filename(run) {
 
  int ipos = filename.find(' ');
  if (ipos < string::npos) filename.erase(ipos);
 
  if (!strcmp(profile, "const")) set_profile = Constant;
  else if (!strcmp(profile, "exp")) set_profile = Exp;
  else if (!strcmp(profile, "expr")) set_profile = Expr;
  else if (!strcmp(profile, "radial")) set_profile = Radial;
  else if (!strcmp(profile, "exp_radial")) set_profile = ExpRadial;
  else if (!strcmp(profile, "blasi")) set_profile = Blasi;
  else {
  cerr << "What DPerpType did you choose?" << endl;
  exit(-2);
  }
  }
*/
