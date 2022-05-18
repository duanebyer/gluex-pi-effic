#include "DSelector_omega_pi_effic.h"

#include <algorithm>
#include <cctype>
#include <map>
#include <set>
#include <sstream>
#include <string>

#include <TDirectory.h>
#include <TMath.h>
#include <TPRegexp.h>
#include <TString.h>

// Length units in cm, time units in ns.
Double_t const MASS_PI0 = 0.1349770; // GeV
Double_t const SPEED_OF_LIGHT = 29.9792458; // cm/ns
Double_t const PI = TMath::Pi();
Double_t const DEG_TO_RAD = TMath::DegToRad();
Double_t const RAD_TO_DEG = TMath::RadToDeg();

// Ranges on histograms.
// TODO: Adjust these back to what they were later.
Double_t const MASS_MIN = 0.0;//0.3;
Double_t const MASS_MAX = 3.5;//1.3;
UInt_t const BINS_1D = 400;//600
UInt_t const BINS_2D = 100;

struct ParticleId {
	Particle_t p; 
	Int_t id;
};
struct UniquenessTracker {
	std::set<std::map<Particle_t, std::set<Int_t> > > dUsedSoFar;
	Bool_t Check(std::vector<ParticleId> const& locParticleIds) {
		// TODO: Does this uniqueness tracking actually work?
		std::map<Particle_t, std::set<Int_t> > locUsedThisCombo;
		for (ParticleId locParticleId : locParticleIds) {
			locUsedThisCombo[locParticleId.p].insert(locParticleId.id);
		}
		auto locInsertResult = dUsedSoFar.insert(locUsedThisCombo);
		return locInsertResult.second;
	}
};

Double_t NormalizeAngle(Double_t locAngle) {
	while (locAngle > PI) {
		locAngle -= 2.*PI;
	}
	while (locAngle <= -PI) {
		locAngle += 2.*PI;
	}
	return locAngle;
}

std::string Trim(std::string locStr) {
	auto check = [](char ch) {
		return !std::isspace(ch);
	};
	locStr.erase(
		locStr.begin(),
		std::find_if(locStr.begin(), locStr.end(), check));
	locStr.erase(
		std::find_if(locStr.rbegin(), locStr.rend(), check).base(),
		locStr.end());
	return locStr;
}

Bool_t CheckCut(Bool_t locCut[CUT_COUNT], UInt_t locVariant) {
	for (UInt_t loc_i = 0; loc_i < CUT_COUNT; ++loc_i) {
		if (CUT_VARIANT[locVariant].dEnabled[loc_i] && !locCut[loc_i]) {
			return false;
		}
	}
	return true;
}

Bool_t CheckCut(Bool_t locCut[CUT_COUNT], Bool_t locCutFound[CUT_FOUND_COUNT], UInt_t locVariant) {
	if (!CheckCut(locCut, locVariant)) {
		return false;
	}
	for (UInt_t loc_i = 0; loc_i < CUT_FOUND_COUNT; ++loc_i) {
		if (CUT_VARIANT[locVariant].dEnabledFound[loc_i] && !locCutFound[loc_i]) {
			return false;
		}
	}
	return true;
}

// Attempts to read a string from the RCDB.
enum class RCDBType {
	RCND,
	CCDB,
};
Bool_t QueryRCDB(
		RCDBType type,
		UInt_t runNumber, std::string const& key,
		std::string* result) {
	int error;
	std::ostringstream command;
	UInt_t discardLines;
	if (type == RCDBType::RCND) {
		command << "rcnd " << runNumber << " " << key;
		discardLines = 0;
	} else if (type == RCDBType::CCDB) {
		command << "ccdb dump " << key << " -r " << runNumber;
		discardLines = 1;
	} else {
		return false;
	}
	std::FILE* file = gSystem->OpenPipe(command.str().c_str(), "r");
	char buffer[1024];
	for (UInt_t i = 0; i < discardLines; ++i) {
		std::fgets(buffer, sizeof(buffer) / sizeof(char), file);
	}
	std::fgets(buffer, sizeof(buffer) / sizeof(char), file);
	error = std::ferror(file);
	gSystem->ClosePipe(file);
	if (error == 0) {
		*result = Trim(buffer);
		return true;
	} else {
		return false;
	}
}

// Read the polarization from the RCDB.
Bool_t QueryPolarization(UInt_t runNumber, Polarization* polarization) {
	std::string polStr;
	if (!QueryRCDB(RCDBType::RCND, runNumber, "polarization_direction", &polStr)) {
		return false;
	}
	if (polStr == "PARA") {
		*polarization = Polarization::PARA;
		return true;
	} else if (polStr == "PERP") {
		*polarization = Polarization::PERP;
		return true;
	} else if (polStr == "N/A") {
		*polarization = Polarization::AMO;
		return true;
	} else {
		return false;
	}
}

// Read the beam bunch spacing from the RCDB.
Bool_t QueryBeamBunchSpacing(UInt_t runNumber, Double_t* spacing) {
	std::string spacingStr;
	if (!QueryRCDB(RCDBType::CCDB, runNumber, "PHOTON_BEAM/RF/beam_period", &spacingStr)) {
		return false;
	}
	std::istringstream spacingStrStream(spacingStr);
	spacingStrStream >> *spacing;
	return (Bool_t) spacingStrStream;
}

// Read the accidental scaling factor.
Bool_t QueryAccidScalingFactor(UInt_t runNumber, AccidScalingFactor* scaling) {
	std::string scalingStr;
	if (!QueryRCDB(RCDBType::CCDB, runNumber, "ANALYSIS/accidental_scaling_factor", &scalingStr)) {
		return false;
	}
	std::istringstream scalingStrStream(scalingStr);
	Double_t blank;
	scalingStrStream
		>> scaling->hodoscopeHiFactor >> blank
		>> scaling->hodoscopeLoFactor >> blank
		>> scaling->microscopeFactor >> blank
		>> scaling->TAGMEnergyBoundHi >> scaling->TAGMEnergyBoundLo;
	return (Bool_t) scalingStrStream;
}

void DSelector_omega_pi_effic::Init(TTree* locTree) {
	// The Init() function is called when the selector needs to initialize a
	// new tree or chain.  Typically here the branch addresses and branch
	// pointers of the tree will be set.  Init() will be called many times when
	// running on PROOF (once per file to be processed).

	// Set output file name (can be overriden by user in PROOF).
	dOutputFileName = "omega_pi_effic.root";
	// If blank, doesn't output tree.
	dOutputTreeFileName = "";
	// Output flat tree (one combo per tree entry), "" for none.
	dFlatTreeFileName = "";
	// Tree name (if blank, chooses a default).
	dFlatTreeName = "";
	// Don't save default branches, reduce disk footprint.
	dSaveDefaultFlatBranches = false;
	// False: save particles as TLorentzVector objects.
	// True: save as four doubles instead.
	dSaveTLorentzVectorsAsFundamentaFlatTree = false;

	dInputTreeName = locTree->GetName();
	Info("Init", "Processing tree: %s", dInputTreeName.Data());
	// If nothing found, these are the default parameters.
	dNumBeamBunches = 1;
	dMassConstraint = -1;
	dUnusedTracks = 0;
	dExtraChargedTracks = 3;
	dFitType = 4;
	TPRegexp locRegex("^(pi0pipmisspim|pi0pimmisspip)_(?:(_B\\d)|(_M\\d)|(_U\\d)|(_T\\d)|(_F\\d))+(?:_.*)?_Tree$");
	TObjArray* locSubstrArr = locRegex.MatchS(dInputTreeName);
	if (locSubstrArr->GetEntries() == 0) {
		Warning("Init", "Unable to extract tree parameters, using defaults");
	} else {
		TString locPiSignStr = ((TObjString*) locSubstrArr->At(1))->GetString();
		// TODO: Actually store the pion signs and validate them against the
		// data?
		if (locPiSignStr == "pi0pipmisspim") {
			Info("Init", "Pi1:+, Pi2:-");
		} else if (locPiSignStr == "pi0pimmisspip") {
			Info("Init", "Pi1:-, Pi2:+");
		} else {
			Warning("Init", "Error extracting pi sign");
		}
		for (Int_t loc_i = 2; loc_i < locSubstrArr->GetEntries(); ++loc_i) {
			TString locParamStr = ((TObjString*) locSubstrArr->At(loc_i))->GetString();
			if (locParamStr(1) == 'B') {
				dNumBeamBunches = locParamStr(2) - '0';
			} else if (locParamStr(1) == 'M') {
				// TODO: Technically, this is buggy because 'M' can be specified
				// multiple times.
				dMassConstraint = locParamStr(2) - '0';
			} else if (locParamStr(1) == 'U') {
				dUnusedTracks = locParamStr(2) - '0';
			} else if (locParamStr(1) == 'T') {
				dExtraChargedTracks = locParamStr(2) - '0';
			} else if (locParamStr(1) == 'F') {
				dFitType = locParamStr(2) - '0';
			}
		}
		Info("Init", "B:%d, M:%d, U:%d, T:%d, F:%d",
			dNumBeamBunches,
			dMassConstraint,
			dUnusedTracks,
			dExtraChargedTracks,
			dFitType);
	}

	dUseDefaults = false;
	char const* chUseDefaults = gSystem->Getenv("USE_DEFAULTS");
	if (chUseDefaults != nullptr) {
		dUseDefaults = true;
	}
	Info("Init", "Using defaults? %d", dUseDefaults);

	dAccidMethod = AccidMethod::NONE;
	char const* chAccidMethod = gSystem->Getenv("ACCID_METHOD");
	if (chAccidMethod != nullptr) {
		TString envAccidMethod(chAccidMethod);
		envAccidMethod.ToUpper();
		if (envAccidMethod == "NONE") {
			dAccidMethod = AccidMethod::NONE;
		} else if (envAccidMethod == "SUBTRACT_ALL") {
			dAccidMethod = AccidMethod::SUBTRACT_ALL;
		} else if (envAccidMethod == "SUBTRACT_ALL_BUT_NEAREST") {
			dAccidMethod = AccidMethod::SUBTRACT_ALL_BUT_NEAREST;
		} else if (envAccidMethod == "SRC2021") {
			dAccidMethod = AccidMethod::SRC2021;
		} else {
			Fatal("Init", "Unrecognized ACCID_METHOD value %s", envAccidMethod.Data());
		}
		Info("Init", "Accidental subtraction method: %s", envAccidMethod.Data());
	} else {
		Fatal("Init", "No ACCID_METHOD specified");
	}

	//Because this function gets called for each TTree in the TChain, we must
	//be careful: We need to re-initialize the tree interface & branch
	//wrappers, but don't want to recreate histograms.

	// Return early if this is the first TTree in the TChain.
	Bool_t locInitializedPriorFlag = dInitializedFlag;
	DSelector::Init(locTree);
	if (locInitializedPriorFlag) {
		return;
	}

	Get_ComboWrappers();
	dPreviousRunNumber = 0;

	// Analysis actions, executed in order of added to dAnalysisActions.
	dAnalysisActions.push_back(new DHistogramAction_ParticleID(dComboWrapper, false));
	//dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 0.5, KPlus, SYS_BCAL));
	// TODO: PIDFOM not included in original analysis.
	dAnalysisActions.push_back(new DHistogramAction_PIDFOM(dComboWrapper));
	//dAnalysisActions.push_back(new DCutAction_PIDFOM(dComboWrapper, KPlus, 0.1));
	//dAnalysisActions.push_back(new DCutAction_EachPIDFOM(dComboWrapper, 0.1));
	//dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, false, Lambda, 1000, 1.0, 1.2, "Lambda"));
	//dAnalysisActions.push_back(new DHistogramAction_MissingMassSq(dComboWrapper, false, 1000, -0.1, 0.1));
	dAnalysisActions.push_back(new DHistogramAction_KinFitResults(dComboWrapper));
	//dAnalysisActions.push_back(new DCutAction_MissingMassSq(dComboWrapper, false, -0.03, 0.02));
	//dAnalysisActions.push_back(new DCutAction_ShowerQuality(dComboWrapper, SYS_FCAL, 0.5));
	dAnalysisActions.push_back(new DHistogramAction_BeamEnergy(dComboWrapper, false));
	//dAnalysisActions.push_back(new DCutAction_BeamEnergy(dComboWrapper, false, 8.2, 8.8));
	dAnalysisActions.push_back(new DHistogramAction_ParticleComboKinematics(dComboWrapper, false));

	// TODO: Consider putting cuts as analysis actions here.

	// Example of a manual action.
	std::deque<Particle_t> etaDecay;
	etaDecay.push_back(PiPlus);
	etaDecay.push_back(PiMinus);
	// Put neutral pion or two gammas here?
	etaDecay.push_back(Gamma);
	etaDecay.push_back(Gamma);
	dAnalyzeCutActions = new DHistogramAction_AnalyzeCutActions(
		dAnalysisActions,
		dComboWrapper,
		false, 0,
		etaDecay,
		1000, 0.9, 2.4,
		"CutActionEffect" );

	Initialize_Actions();
	dAnalyzeCutActions->Initialize();

	// CREATE HISTOGRAMS
	TDirectory* locDir = gDirectory;
	for (UInt_t locCutVariant = 0; locCutVariant < CUT_VARIANT_COUNT; ++locCutVariant) {
		TDirectory* locCutDir = locDir->mkdir(Form("Cut%s", CUT_VARIANT[locCutVariant].dName));
		locCutDir->cd();
		dHist_KinFitCL[locCutVariant] = new TH1D(
			"KinFitCL",
			";Kin. Fit CL",
			BINS_1D, 0., 1.);
		dHist_MissingMass[locCutVariant] = new TH1D(
			"MissingMass",
			";Missing Mass (GeV/c^{2})",
			BINS_1D, -0.25, 0.25);
		dHist_MissingMass_Found[locCutVariant] = new TH1D(
			"MissingMass_Found",
			"Pion Found;Missing Mass (GeV/c^{2})",
			BINS_1D, 0.0, 0.25);
		dHist_MissingMass_Missing[locCutVariant] = new TH1D(
			"MissingMass_Missing",
			"Pion Missing;Missing Mass (GeV/c^{2})",
			BINS_1D, 0.0, 0.25);
		dHist_MissingMass_Measured[locCutVariant] = new TH1D(
			"MissingMass_Measured",
			";Missing Mass (GeV/c^{2})",
			BINS_1D, -0.25, 0.25);
		dHist_MissingMass_Measured_Found[locCutVariant] = new TH1D(
			"MissingMass_Measured_Found",
			"Pion Found;Missing Mass (GeV/c^{2})",
			BINS_1D, 0.0, 0.25);
		dHist_MissingMass_Measured_Missing[locCutVariant] = new TH1D(
			"MissingMass_Measured_Missing",
			"Pion Missing;Missing Mass (GeV/c^{2})",
			BINS_1D, 0.0, 0.25);
		dHist_BeamEnergy[locCutVariant] = new TH1D(
			"BeamEnergy",
			";Beam Energy (GeV)",
			BINS_1D, 0., 12.);
		dHist_FoundPi2[locCutVariant] = new TH1D(
			"FoundPi2",
			";Number of pions found",
			4, 0., 4.);
		dHist_Pi0MassSq_Measured_Found[locCutVariant] = new TH1D(
			"Pi0MassSq_Measured_Found",
			";Neutral Pion Mass Sq (GeV/c^{2})^{2}",
			BINS_1D, 0.09, 0.18);
		dHist_Pi0MassSq_Measured_Missing[locCutVariant] = new TH1D(
			"Pi0MassSq_Measured_Missing",
			";Neutral Pion Mass Sq (GeV/c^{2})^{2}",
			BINS_1D, 0.09, 0.18);
		dHist_Pi0MassSq_Found[locCutVariant] = new TH1D(
			"Pi0MassSq_Found",
			";Neutral Pion Mass Sq (GeV/c^{2})^{2}",
			BINS_1D, 0.09, 0.18);
		dHist_Pi0MassSq_Missing[locCutVariant] = new TH1D(
			"Pi0MassSq_Missing",
			";Neutral Pion Mass Sq (GeV/c^{2})^{2}",
			BINS_1D, 0.09, 0.18);
		dHist_ProtonMomentum[locCutVariant] = new TH3D(
			"ProtonMomentum",
			"Proton;Momentum (GeV/c^{2});Polar Angle (deg.);Azimuthal Angle (deg.)",
			BINS_2D, 0., 10.,
			BINS_2D, 0., 180.,
			BINS_2D, -180., 180.);
		dHist_Pi0Momentum[locCutVariant] = new TH3D(
			"Pi0Momentum",
			"Neutral pion;Momentum (GeV/c^{2});Polar Angle (deg.);Azimuthal Angle (deg.)",
			BINS_2D, 0., 10.,
			BINS_2D, 0., 180.,
			BINS_2D, -180., 180.);
		dHist_Pi1Momentum[locCutVariant] = new TH3D(
			"Pi1Momentum",
			"Pion;Momentum (GeV/c^{2});Polar Angle (deg.);Azimuthal Angle (deg.)",
			BINS_2D, 0., 10.,
			BINS_2D, 0., 180.,
			BINS_2D, -180., 180.);
		dHist_Pi2Momentum_Kinfit[locCutVariant] = new TH3D(
			"Pi2Momentum_Kinfit",
			"Pion;Momentum (GeV/c^{2});Polar Angle (deg.);Azimuthal Angle (deg.)",
			BINS_2D, 0., 10.,
			BINS_2D, 0., 180.,
			BINS_2D, -180., 180.);
		dHist_Pi2Momentum_Found[locCutVariant] = new TH3D(
			"Pi2Momentum_Found",
			"Pion Found;Momentum (GeV/c^{2});Polar Angle (deg.);Azimuthal Angle (deg.)",
			BINS_2D, 0., 10.,
			BINS_2D, 0., 180.,
			BINS_2D, -180., 180.);
		dHist_Pi2Momentum_Missing[locCutVariant] = new TH3D(
			"Pi2Momentum_Missing",
			"Pion Missing;Momentum (GeV/c^{2});Polar Angle (deg.);Azimuthal Angle (deg.)",
			BINS_2D, 0., 10.,
			BINS_2D, 0., 180.,
			BINS_2D, -180., 180.);
		dHist_Proton3PiPhi_Kinfit[locCutVariant] = new TH2D(
			"Proton3PiPhi_Kinfit",
			";Proton Azimuthal Angle (deg.);3 Pions Azimuthal Angle (deg.)",
			BINS_2D, -180., 180.,
			BINS_2D, -180., 180.);
		dHist_Proton3PiCollinearity_Kinfit[locCutVariant] = new TH1D(
			"Proton3PiCollinearity_Kinfit",
			";Proton Azimuthal Angle - 3 Pions Azimuthal Angle (deg.)",
			BINS_1D, -20., 20.);
		dHist_Proton3PiPhi_Found[locCutVariant] = new TH2D(
			"Proton3PiPhi_Found",
			";Proton Azimuthal Angle (deg.);3 Pions Azimuthal Angle (deg.)",
			BINS_2D, -180., 180.,
			BINS_2D, -180., 180.);
		dHist_Proton3PiCollinearity_Found[locCutVariant] = new TH1D(
			"Proton3PiCollinearity_Found",
			";Proton Azimuthal Angle - 3 Pions Azimuthal Angle (deg.)",
			BINS_1D, -180., 180.);
		dHist_Proton3PiCollinearity_Mixed[locCutVariant] = new TH1D(
			"Proton3PiCollinearity_Mixed",
			";Proton Azimuthal Angle - 3 Pions Azimuthal Angle (deg.)",
			BINS_1D, -20., 20.);
		// Diagnostic histograms.
		dHist_ProtonMomentum_MMOP_Found[locCutVariant] = new TH2D(
			"ProtonMomentum_MMOP_Found",
			"Pion Found;Proton Momentum (GeV/c^{2});MMOP (GeV/c^{2})",
			BINS_2D, 0., 10.,
			BINS_2D, MASS_MIN, MASS_MAX);
		dHist_ProtonMomentum_MMOP_Missing[locCutVariant] = new TH2D(
			"ProtonMomentum_MMOP_Missing",
			"Pion Missing;Proton Momentum (GeV/c^{2});MMOP (GeV/c^{2})",
			BINS_2D, 0., 10.,
			BINS_2D, MASS_MIN, MASS_MAX);
		dHist_Proton3PiCollinearity_Mixed_MMOP_Found[locCutVariant] = new TH2D(
			"Proton3PiCollinearity_Mixed_MMOP_Found",
			"Pion Found;Proton Azimuthal Angle - 3 Pions Azimuthal Angle (deg.);MMOP (GeV/c^{2})",
			BINS_2D, -20., 20.,
			BINS_2D, MASS_MIN, MASS_MAX);
		dHist_Proton3PiCollinearity_Mixed_MMOP_Missing[locCutVariant] = new TH2D(
			"Proton3PiCollinearity_Mixed_MMOP_Missing",
			"Pion Missing;Proton Azimuthal Angle - 3 Pions Azimuthal Angle (deg.);MMOP (GeV/c^{2})",
			BINS_2D, -20., 20.,
			BINS_2D, MASS_MIN, MASS_MAX);
		dHist_Proton3PiCollinearity_MMOP_Found[locCutVariant] = new TH2D(
			"Proton3PiCollinearity_MMOP_Found",
			"Pion Found;Proton Azimuthal Angle - 3 Pions Azimuthal Angle (deg.);MMOP (GeV/c^{2})",
			BINS_2D, -180., 180.,
			BINS_2D, MASS_MIN, MASS_MAX);
		dHist_KinFitCL_MMOP_Found[locCutVariant] = new TH2D(
			"KinFitCL_MMOP_Found",
			"Pion Found;Kin. Fit CL;MMOP (GeV/c^{2})",
			BINS_2D, 0., 1.,
			BINS_2D, MASS_MIN, MASS_MAX);
		dHist_KinFitCL_MMOP_Missing[locCutVariant] = new TH2D(
			"KinFitCL_MMOP_Missing",
			"Pion Missing;Kin. Fit CL;MMOP (GeV/c^{2})",
			BINS_2D, 0., 1.,
			BINS_2D, MASS_MIN, MASS_MAX);
		dHist_MMOP_Measured_Found[locCutVariant] = new TH1D(
			"MMOP_Measured_Found",
			"Pion Found;MMOP (GeV/c^{2})",
			BINS_1D, MASS_MIN, MASS_MAX);
		dHist_MMOP_Measured_Missing[locCutVariant] = new TH1D(
			"MMOP_Measured_Missing",
			"Pion Missing;MMOP (GeV/c^{2})",
			BINS_1D, MASS_MIN, MASS_MAX);
		dHist_MMOP_Found[locCutVariant] = new TH1D(
			"MMOP_Found",
			"Pion Found;MMOP (GeV/c^{2})",
			BINS_1D, MASS_MIN, MASS_MAX);
		dHist_MMOP_Missing[locCutVariant] = new TH1D(
			"MMOP_Missing",
			"Pion Missing;MMOP (GeV/c^{2})",
			BINS_1D, MASS_MIN, MASS_MAX);
		dHist_M3PI_Measured_Found[locCutVariant] = new TH1D(
			"M3PI_Measured_Found",
			"Pion Found;M_{#pi^{0}#pi^{+}#pi^{-}} (GeV/c^{2})",
			BINS_1D, MASS_MIN, MASS_MAX);
		dHist_M3PI_Kinfit_Found[locCutVariant] = new TH1D(
			"M3PI_Found",
			"Pion Found;M_{#pi^{0}#pi^{+}#pi^{-}} (GeV/c^{2})",
			BINS_1D, MASS_MIN, MASS_MAX);
		dHist_M3PI_Kinfit_Missing[locCutVariant] = new TH1D(
			"M3PI_Missing",
			"Pion Found;M_{#pi^{0}#pi^{+}#pi^{-}} (GeV/c^{2})",
			BINS_1D, MASS_MIN, MASS_MAX);
	}
	locDir->cd();

	/************************** EXAMPLE USER INITIALIZATION: CUSTOM OUTPUT BRANCHES - MAIN TREE *************************/

	//EXAMPLE MAIN TREE CUSTOM BRANCHES (OUTPUT ROOT FILE NAME MUST FIRST BE GIVEN!!!! (ABOVE: TOP)):
	//The type for the branch must be included in the brackets
	//1st function argument is the name of the branch
	//2nd function argument is the name of the branch that contains the size of the array (for fundamentals only)
	/*
	dTreeInterface->Create_Branch_Fundamental<Int_t>("my_int"); //fundamental = char, int, float, double, etc.
	dTreeInterface->Create_Branch_FundamentalArray<Int_t>("my_int_array", "my_int");
	dTreeInterface->Create_Branch_FundamentalArray<Float_t>("my_combo_array", "NumCombos");
	dTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("my_p4");
	dTreeInterface->Create_Branch_ClonesArray<TLorentzVector>("my_p4_array");
	*/

	if (dFlatTreeInterface != nullptr) {
		dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("accidweight");
	}

	//EXAMPLE FLAT TREE CUSTOM BRANCHES (OUTPUT ROOT FILE NAME MUST FIRST BE GIVEN!!!! (ABOVE: TOP)):
	//The type for the branch must be included in the brackets
	//1st function argument is the name of the branch
	//2nd function argument is the name of the branch that contains the size of the array (for fundamentals only)
	/*
	dFlatTreeInterface->Create_Branch_Fundamental<Int_t>("flat_my_int"); //fundamental = char, int, float, double, etc.
	dFlatTreeInterface->Create_Branch_FundamentalArray<Int_t>("flat_my_int_array", "flat_my_int");
	dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("flat_my_p4");
	dFlatTreeInterface->Create_Branch_ClonesArray<TLorentzVector>("flat_my_p4_array");
	*/

	/************************************* ADVANCED EXAMPLE: CHOOSE BRANCHES TO READ ************************************/

	//TO SAVE PROCESSING TIME
		//If you know you don't need all of the branches/data, but just a subset of it, you can speed things up
		//By default, for each event, the data is retrieved for all branches
		//If you know you only need data for some branches, you can skip grabbing data from the branches you don't need
		//Do this by doing something similar to the commented code below

	//dTreeInterface->Clear_GetEntryBranches(); //now get none
	//dTreeInterface->Register_GetEntryBranch("Proton__P4"); //manually set the branches you want

	dIsMC = (dTreeInterface->Get_Branch("MCWeight") != NULL);
}

Bool_t DSelector_omega_pi_effic::Process(Long64_t locEntry) {
	// The Process() function is called for each entry in the tree. The entry
	// argument specifies which entry in the currently loaded tree is to be
	// processed.
	//
	// This function should contain the "body" of the analysis. It can contain
	// simple or elaborate selection criteria, run algorithms on the data of
	// the event and typically fill histograms.
	//
	// The processing can be stopped by calling Abort().  Use fStatus to set
	// the return value of TTree::Process().  The return value is currently not
	// used.

	// Gets the data from the tree for the entry.
	DSelector::Process(locEntry);

	// Get polarization orientation from the run number.
	UInt_t locRunNumber = Get_RunNumber();
	if(locRunNumber != dPreviousRunNumber) {
		dPreviousRunNumber = locRunNumber;
		if (!QueryPolarization(locRunNumber, &dPolarization)) {
			dPolarization = Polarization::AMO;
			if (!dUseDefaults) {
				Fatal("Process", "Couldn't query RCDB for polarization");
			}
		}
		if (!QueryBeamBunchSpacing(locRunNumber, &dBeamBunchSpacing)) {
			dBeamBunchSpacing = 4.008016032;
			if (!dUseDefaults) {
				Fatal("Process", "Couldn't query RCDB for beam bunch spacing");
			}
		}
		if (!QueryAccidScalingFactor(locRunNumber, &dAccidScalingFactor)) {
			dAccidScalingFactor = AccidScalingFactor();
			if (!dUseDefaults) {
				Fatal("Process", "Couldn't query RCDB for accidental scaling factor");
			}
		}
	}

	Reset_Actions_NewEvent();
	// Manual actions must have Reset_NewEvent() called.
	dAnalyzeCutActions->Reset_NewEvent();

	// UNIQUENESS TRACKING
	// Sometimes, some content is the exact same between one combo and the next
	// e.g. maybe two combos have different beam particles, but the same data
	// for the final-state When histogramming, you don't want to double-count
	// when this happens: artificially inflates your signal (or background) So,
	// for each quantity you histogram, keep track of what particles you used
	// (for a given combo) Then for each combo, just compare to what you used
	// before, and make sure it's unique

	// Create uniqueness trackers.
	UniquenessTracker locTracker_BeamEnergy[CUT_VARIANT_COUNT];
	UniquenessTracker locTracker_MissingMass[CUT_VARIANT_COUNT];
	UniquenessTracker locTracker_MissingMass_Found[CUT_VARIANT_COUNT];
	UniquenessTracker locTracker_MissingMass_Missing[CUT_VARIANT_COUNT];
	UniquenessTracker locTracker_Proton[CUT_VARIANT_COUNT];
	UniquenessTracker locTracker_Pi0[CUT_VARIANT_COUNT];
	UniquenessTracker locTracker_Pi0_Found[CUT_VARIANT_COUNT];
	UniquenessTracker locTracker_Pi0_Missing[CUT_VARIANT_COUNT];
	UniquenessTracker locTracker_Pi1[CUT_VARIANT_COUNT];
	UniquenessTracker locTracker_Pi2_Kinfit[CUT_VARIANT_COUNT];
	UniquenessTracker locTracker_Pi2_Found[CUT_VARIANT_COUNT];
	UniquenessTracker locTracker_Pi2_Missing[CUT_VARIANT_COUNT];
	UniquenessTracker locTracker_Proton3Pi_Kinfit[CUT_VARIANT_COUNT];
	UniquenessTracker locTracker_Proton3Pi_Kinfit_Found[CUT_VARIANT_COUNT];
	UniquenessTracker locTracker_Proton3Pi_Kinfit_Missing[CUT_VARIANT_COUNT];
	UniquenessTracker locTracker_Proton3Pi_Found[CUT_VARIANT_COUNT];
	UniquenessTracker locTracker_MMOP_Found[CUT_VARIANT_COUNT];
	UniquenessTracker locTracker_MMOP_Missing[CUT_VARIANT_COUNT];
	UniquenessTracker locTracker_MMOP[CUT_VARIANT_COUNT];
	UniquenessTracker locTracker_M3PI_Found[CUT_VARIANT_COUNT];
	UniquenessTracker locTracker_M3PI_Kinfit_Found[CUT_VARIANT_COUNT];
	UniquenessTracker locTracker_M3PI_Kinfit_Missing[CUT_VARIANT_COUNT];

	// CUSTOM OUTPUT BRANCHES
	//Int_t locMyInt = 7;
	//dTreeInterface->Fill_Fundamental<Int_t>("my_int", locMyInt);
	//TLorentzVector locMyP4(4.0, 3.0, 2.0, 1.0);
	//dTreeInterface->Fill_TObject<TLorentzVector>("my_p4", locMyP4);
	//for(int loc_i = 0; loc_i < locMyInt; ++loc_i) {
	//	dTreeInterface->Fill_Fundamental<Int_t>("my_int_array", 3*loc_i, loc_i); //2nd argument = value, 3rd = array index
	//}

	// NUMBER OF CHARGED TRACKS
	set<Int_t> locChargedTrackIDs;
	for (UInt_t loc_j = 0; loc_j < Get_NumChargedHypos(); ++loc_j) {
		dChargedHypoWrapper->Set_ArrayIndex(loc_j);
		locChargedTrackIDs.insert(dChargedHypoWrapper->Get_TrackID());
	}
	int locChargedTrackCount = locChargedTrackIDs.size();

	// COMBO LOOP
	for(UInt_t loc_i = 0; loc_i < Get_NumCombos(); ++loc_i) {
		dComboWrapper->Set_ComboIndex(loc_i);
		if(dComboWrapper->Get_IsComboCut()) {
			continue;
		}

		// LOAD PARTICLE INFO
		Int_t locBeamID = dComboBeamWrapper->Get_BeamID();
		Int_t locPi1TrackID = dPi1Wrapper->Get_TrackID();
		Int_t locProtonTrackID = dProtonWrapper->Get_TrackID();
		Int_t locPhoton1NeutralID = dPhoton1Wrapper->Get_NeutralID();
		Int_t locPhoton2NeutralID = dPhoton2Wrapper->Get_NeutralID();

		// The wrappers give the kinematic fit if it was performed, otherwise
		// give the measured momenta.
		TLorentzVector locBeamP4 = dComboBeamWrapper->Get_P4();
		TLorentzVector locPi1P4 = dPi1Wrapper->Get_P4();
		TLorentzVector locProtonP4 = dProtonWrapper->Get_P4();
		//TLorentzVector locPi2P4 = dPi2Wrapper->Get_P4();
		TLorentzVector locPhoton1P4 = dPhoton1Wrapper->Get_P4();
		TLorentzVector locPhoton2P4 = dPhoton2Wrapper->Get_P4();
		TLorentzVector locPi2P4_Kinfit = dPi2Wrapper->Get_P4();
		TLorentzVector locPi0P4 = locPhoton1P4 + locPhoton2P4;
		TLorentzVector loc3PiP4_Kinfit = locPi0P4 + locPi1P4 + locPi2P4_Kinfit; 

		// Measured momenta.
		TLorentzVector locBeamP4_Measured = dComboBeamWrapper->Get_P4_Measured();
		TLorentzVector locPi1P4_Measured = dPi1Wrapper->Get_P4_Measured();
		TLorentzVector locProtonP4_Measured = dProtonWrapper->Get_P4_Measured();
		TLorentzVector locPhoton1P4_Measured = dPhoton1Wrapper->Get_P4_Measured();
		TLorentzVector locPhoton2P4_Measured = dPhoton2Wrapper->Get_P4_Measured();
		TLorentzVector locPi0P4_Measured = locPhoton1P4_Measured + locPhoton2P4_Measured;
		TLorentzVector loc3PiP4_Mixed = locPi0P4_Measured + locPi1P4_Measured + locPi2P4_Kinfit; 

		// Combine 4-vectors.
		TLorentzVector locMissingP4 = locBeamP4 + dTargetP4;
		locMissingP4 -= locProtonP4 + locPi0P4 + locPi1P4;
		TLorentzVector locMissingP4_Measured = locBeamP4_Measured + dTargetP4;
		locMissingP4_Measured -= locProtonP4_Measured + locPi0P4_Measured + locPi1P4_Measured;

		// Kinematic quantities.
		Double_t locMissingMass = locMissingP4.M();
		Double_t locMissingMass_Measured = locMissingP4_Measured.M();
		Double_t locPi0MassSq_Measured = locPi0P4_Measured.M2();
		Double_t locUnusedShowerEnergy = dComboWrapper->Get_Energy_UnusedShowers();
		Double_t locKinFitCL = dComboWrapper->Get_ConfidenceLevel_KinFit("");

		Double_t locProton3PiCollinearity_Kinfit = NormalizeAngle(
			locProtonP4.Phi() - loc3PiP4_Kinfit.Phi() + PI);
		Double_t locProton3PiCollinearity_Mixed = NormalizeAngle(
			locProtonP4_Measured.Phi() - loc3PiP4_Mixed.Phi() + PI);

		Double_t locMMOP = (locBeamP4 + dTargetP4 - locProtonP4).M();
		Double_t locMMOP_Measured = (locBeamP4_Measured + dTargetP4 - locProtonP4_Measured).M();
		Double_t locM3PI_Kinfit = (locPi0P4 + locPi1P4 + locPi2P4_Kinfit).M();

		// RF timing information.
		TLorentzVector locBeamX4_Measured = dComboBeamWrapper->Get_X4_Measured();
		Double_t locRFTime_Measured = dComboWrapper->Get_RFTime();
		Double_t locBeamDeltaT_Measured = locBeamX4_Measured.T()
			- (locRFTime_Measured + (locBeamX4_Measured.Z() - dTargetCenter.Z())/SPEED_OF_LIGHT);
		Int_t locRelBeamBucket = TMath::Nint(locBeamDeltaT_Measured/dBeamBunchSpacing);
		Double_t locAccidWeight;
		if (dAccidMethod == AccidMethod::NONE) {
			// No accidental subtraction.
			locAccidWeight = 1.;
		} else if (dAccidMethod == AccidMethod::SUBTRACT_ALL) {
			// Accidental subtraction.
			if (locRelBeamBucket == 0) {
				locAccidWeight = 1.;
			} else if (TMath::Abs(locRelBeamBucket) <= dNumBeamBunches) {
				locAccidWeight = -1./(2*dNumBeamBunches);
				locAccidWeight *= dAccidScalingFactor.GetFactor(locBeamP4.E());
			} else {
				locAccidWeight = 0.;
			}
		} else if (dAccidMethod == AccidMethod::SUBTRACT_ALL_BUT_NEAREST) {
			// Accidental subtraction skipping the adjacent peaks.
			if (locRelBeamBucket == 0) {
				locAccidWeight = 1.;
			} else if (TMath::Abs(locRelBeamBucket) == 1) {
				locAccidWeight = 0.;
			} else if (TMath::Abs(locRelBeamBucket) <= dNumBeamBunches) {
				locAccidWeight = -1./(2*(dNumBeamBunches - 1));
				locAccidWeight *= dAccidScalingFactor.GetFactor(locBeamP4.E());
			} else {
				locAccidWeight = 0.;
			}
		} else if (dAccidMethod == AccidMethod::SRC2021) {
			// Can't trust the RCDB, use hard-coded values for now.
			locAccidWeight = 1.;
			dBeamBunchSpacing = 4.008016032;
			locRelBeamBucket = TMath::Nint(locBeamDeltaT_Measured/dBeamBunchSpacing);
			if (locRelBeamBucket == 0) {
				locAccidWeight = 1.;
			} else if (TMath::Abs(locRelBeamBucket) <= dNumBeamBunches) {
				locAccidWeight = -1./(2*dNumBeamBunches);
			} else {
				locAccidWeight = 0.;
			}
		} else {
			Fatal("Process", "Unrecognized accidental subtraction method %d", dAccidMethod);
		}

		//Double_t locBunchPeriod = dAnalysisUtilities.Get_BeamBunchPeriod(Get_RunNumber());
		//Double_t locDeltaT_RF = dAnalysisUtilities.Get_DeltaT_RF(Get_RunNumber(), locBeamX4_Measured, dComboWrapper);
		//Int_t locRelBeamBucket = dAnalysisUtilities.Get_RelativeBeamBucket(Get_RunNumber(), locBeamX4_Measured, dComboWrapper); // 0 for in-time events, non-zero integer for out-of-time photons
		//Int_t locNumOutOfTimeBunchesInTree = XXX; //YOU need to specify this number
		//Number of out-of-time beam bunches in tree (on a single side, so that total number out-of-time bunches accepted is 2 times this number for left + right bunches) 

		//Bool_t locSkipNearestOutOfTimeBunch = true; // True: skip events from nearest out-of-time bunch on either side (recommended).
		//Int_t locNumOutOfTimeBunchesToUse = locSkipNearestOutOfTimeBunch ? locNumOutOfTimeBunchesInTree-1:locNumOutOfTimeBunchesInTree; 
		//Double_t locAccidentalScalingFactor = dAnalysisUtilities.Get_AccidentalScalingFactor(Get_RunNumber(), locBeamP4.E(), dIsMC); // Ideal value would be 1, but deviations require added factor, which is different for data and MC.
		//Double_t locAccidentalScalingFactorError = dAnalysisUtilities.Get_AccidentalScalingFactorError(Get_RunNumber(), locBeamP4.E()); // Ideal value would be 1, but deviations observed, need added factor.
		//Double_t locHistAccidWeightFactor = locRelBeamBucket==0 ? 1 : -locAccidentalScalingFactor/(2*locNumOutOfTimeBunchesToUse) ; // Weight by 1 for in-time events, ScalingFactor*(1/NBunches) for out-of-time
		//if(locSkipNearestOutOfTimeBunch && abs(locRelBeamBucket)==1) { // Skip nearest out-of-time bunch: tails of in-time distribution also leak in
		//	dComboWrapper->Set_IsComboCut(true); 
		//	continue; 
		//} 

		// ANALYSIS ACTIONS
		dAnalyzeCutActions->Perform_Action();
		if(!Execute_Actions()) {
			continue;
		}

		// CUTS
		Bool_t locCut[] = {
			locUnusedShowerEnergy <= 1.,
			locChargedTrackCount <= 3,
			locMissingMass_Measured <= 0.25,
			locKinFitCL >= 0.1,
			locPi0MassSq_Measured >= TMath::Sq(MASS_PI0 - 0.015)
				|| locPi0MassSq_Measured <= TMath::Sq(MASS_PI0 + 0.015),
			locProtonP4_Measured.P() >= 0.4,
			locProtonP4_Measured.P() >= 0.5,
			locProtonP4_Measured.P() >= 0.6,
			locProtonP4_Measured.P() >= 0.7,
			locProtonP4_Measured.P() >= 0.8,
			locKinFitCL >= 0.2,
			locKinFitCL >= 0.3,
			locKinFitCL >= 0.4,
			locProtonP4_Measured.Theta()*RAD_TO_DEG >= 20. + 10. * locProtonP4_Measured.P(),
			TMath::Abs(locProton3PiCollinearity_Mixed) <= 1.*DEG_TO_RAD,
		};
		static_assert(
			sizeof(locCut) / sizeof(locCut[0]) == CUT_COUNT,
			"Conflict in CUT_COUNT");
		// TODO: Since we're trying many different cuts, for now we don't set
		// the flag. These cuts are applied to histograms only.
		/*
		if (
				locUnusedShowerEnergy > 1.
				|| locChargedTrackCount > 3
				|| locMissingMass_Measured > 0.25
				|| locKinFitCL < 0.1
				|| locPi0MassSq_Measured < TMath::Sq(MASS_PI0 - 0.015)
				|| locPi0MassSq_Measured > TMath::Sq(MASS_PI0 + 0.015)) {
			dComboWrapper->Set_IsComboCut(true);
			continue;
		}
		*/

		// CHARGED TRACK LOOP
		UInt_t locFoundPi2[CUT_COUNT] = { 0 };
		for (UInt_t loc_j = 0; loc_j < Get_NumChargedHypos(); ++loc_j) {
			dChargedHypoWrapper->Set_ArrayIndex(loc_j);
			// Filter so that we consider all missing pion tracks only.
			Int_t locPi2TrackID = dChargedHypoWrapper->Get_TrackID();
			if (dChargedHypoWrapper->Get_PID() != dPi2PID) {
				continue;
			}
			if (locPi2TrackID == locPi1TrackID || locPi2TrackID == locProtonTrackID) {
				continue;
			}
			// CUTS
			// TODO: If the track reconstructed kinematics is very different
			// from the missing momentum, then skip it. This is slightly
			// different than the cut done in the 2017 study.
			TLorentzVector locPi2P4 = dChargedHypoWrapper->Get_P4();
			TLorentzVector locPi2P4_Measured = dChargedHypoWrapper->Get_P4_Measured();
			TLorentzVector loc3PiP4 = locPi0P4 + locPi1P4 + locPi2P4;
			TLorentzVector loc3PiP4_Measured = locPi0P4_Measured + locPi1P4_Measured + locPi2P4_Measured;
			Double_t locDeltaP = locPi2P4_Measured.P() - locMissingP4_Measured.P();
			Double_t locDeltaTheta = locPi2P4_Measured.Theta() - locMissingP4_Measured.Theta();
			Double_t locDeltaPhi = locPi2P4_Measured.Phi() - locMissingP4_Measured.Phi();

			// KINEMATIC QUANTITIES
			Double_t locM3PI_Measured = (locPi0P4_Measured + locPi1P4_Measured + locPi2P4_Measured).M();
			Double_t locProton3PiCollinearity = NormalizeAngle(
				locProtonP4.Phi() - loc3PiP4.Phi() + PI);
			Double_t locProton3PiCollinearity_Measured = NormalizeAngle(
				locProtonP4_Measured.Phi() - loc3PiP4_Measured.Phi() + PI);

			Bool_t locCutFound[] = {
				TMath::Abs(locDeltaP) <= 3.,
				TMath::Abs(NormalizeAngle(locDeltaTheta)) <= 30.*DEG_TO_RAD,
				TMath::Abs(NormalizeAngle(locDeltaPhi)) <= 90.*DEG_TO_RAD,
				TMath::Abs(locProton3PiCollinearity_Measured) <= 1.*DEG_TO_RAD,
			};
			static_assert(
				sizeof(locCutFound) / sizeof(locCutFound[0]) == CUT_FOUND_COUNT,
				"Conflict in CUT_FOUND_COUNT");

			// FILL HISTOGRAMS
			for (UInt_t locCutVariant = 0; locCutVariant < CUT_VARIANT_COUNT; ++locCutVariant) {
				if (!CheckCut(locCut, locCutFound, locCutVariant)) {
					continue;
				}
				locFoundPi2[locCutVariant] += 1;
				// TODO: The uniqueness tracking for a lot of these is
				// questionable: should the restriction of "finding" a pion be
				// enough to include the found pion in the uniqueness tracking
				// or not?
				if (locTracker_MissingMass_Found[locCutVariant].Check({
						{ Unknown, locBeamID },
						{ dPi1PID, locPi1TrackID },
						{ Proton, locProtonTrackID },
						{ Gamma, locPhoton1NeutralID },
						{ Gamma, locPhoton2NeutralID }})) {
					dHist_MissingMass_Found[locCutVariant]->Fill(
						locMissingMass,
						locAccidWeight);
					dHist_MissingMass_Measured_Found[locCutVariant]->Fill(
						locMissingMass_Measured,
						locAccidWeight);
				}
				if (locTracker_Pi0_Found[locCutVariant].Check({
						{ Gamma, locPhoton1NeutralID },
						{ Gamma, locPhoton2NeutralID }})) {
					dHist_Pi0MassSq_Measured_Found[locCutVariant]->Fill(
						locPi0P4_Measured.M(),
						locAccidWeight);
					dHist_Pi0MassSq_Found[locCutVariant]->Fill(
						locPi0P4.M(),
						locAccidWeight);
				}
				if (locTracker_Pi2_Found[locCutVariant].Check({
						{ dPi2PID, locPi2TrackID }})) {
					// TODO: For non-measured, is this the right way to do
					// uniqueness tracking?
					dHist_Pi2Momentum_Found[locCutVariant]->Fill(
						locPi2P4.P(),
						locPi2P4.Theta()*RAD_TO_DEG,
						locPi2P4.Phi()*RAD_TO_DEG,
						locAccidWeight);
				}
				if (locTracker_Proton3Pi_Found[locCutVariant].Check({
						{ Proton, locProtonTrackID },
						{ Gamma, locPhoton1NeutralID },
						{ Gamma, locPhoton2NeutralID },
						{ dPi1PID, locPi1TrackID },
						{ dPi2PID, locPi2TrackID }})) {
					dHist_Proton3PiCollinearity_MMOP_Found[locCutVariant]->Fill(
						locProton3PiCollinearity*RAD_TO_DEG,
						locMMOP,
						locAccidWeight);
					dHist_Proton3PiPhi_Found[locCutVariant]->Fill(
						locProtonP4.Phi()*RAD_TO_DEG,
						loc3PiP4.Phi()*RAD_TO_DEG,
						locAccidWeight);
					dHist_Proton3PiCollinearity_Found[locCutVariant]->Fill(
						locProton3PiCollinearity*RAD_TO_DEG,
						locAccidWeight);
				}
				if (locTracker_Proton3Pi_Kinfit_Found[locCutVariant].Check({
						{ Proton, locProtonTrackID },
						{ Gamma, locPhoton1NeutralID },
						{ Gamma, locPhoton2NeutralID },
						{ dPi1PID, locPi1TrackID }})) {
					dHist_Proton3PiCollinearity_Mixed_MMOP_Found[locCutVariant]->Fill(
						locProton3PiCollinearity_Mixed*RAD_TO_DEG,
						locMMOP,
						locAccidWeight);
				}
				if (locTracker_MMOP_Found[locCutVariant].Check({
						{ Unknown, locBeamID },
						{ Proton, locProtonTrackID }})) {
					dHist_ProtonMomentum_MMOP_Found[locCutVariant]->Fill(
						locProtonP4.P(),
						locMMOP,
						locAccidWeight);
					dHist_KinFitCL_MMOP_Found[locCutVariant]->Fill(
						locKinFitCL,
						locMMOP,
						locAccidWeight);
					dHist_MMOP_Measured_Found[locCutVariant]->Fill(
						locMMOP_Measured,
						locAccidWeight);
					dHist_MMOP_Found[locCutVariant]->Fill(
						locMMOP,
						locAccidWeight);
				}
				if (locTracker_M3PI_Found[locCutVariant].Check({
						{ dPi1PID, locPi1TrackID },
						{ dPi2PID, locPi2TrackID },
						{ Gamma, locPhoton1NeutralID },
						{ Gamma, locPhoton2NeutralID }})) {
					dHist_M3PI_Measured_Found[locCutVariant]->Fill(
						locM3PI_Measured,
						locAccidWeight);
				}
				if (locTracker_M3PI_Kinfit_Found[locCutVariant].Check({
						{ dPi1PID, locPi1TrackID },
						{ Gamma, locPhoton1NeutralID },
						{ Gamma, locPhoton2NeutralID }})) {
					dHist_M3PI_Kinfit_Found[locCutVariant]->Fill(
						locM3PI_Kinfit,
						locAccidWeight);
				}
			}
		}

		// FILL CUSTOM BRANCHES
		/*
		TLorentzVector locMyComboP4(8.0, 7.0, 6.0, 5.0);
		//for arrays below: 2nd argument is value, 3rd is array index
		//NOTE: By filling here, AFTER the cuts above, some indices won't be updated (and will be whatever they were from the last event)
			//So, when you draw the branch, be sure to cut on "IsComboCut" to avoid these.
		dTreeInterface->Fill_Fundamental<Float_t>("my_combo_array", -2*loc_i, loc_i);
		dTreeInterface->Fill_TObject<TLorentzVector>("my_p4_array", locMyComboP4, loc_i);
		*/

		// FILL HISTOGRAMS
		for (UInt_t locCutVariant = 0; locCutVariant < CUT_VARIANT_COUNT; ++locCutVariant) {
			if (!CheckCut(locCut, locCutVariant)) {
				continue;
			}
			dHist_KinFitCL[locCutVariant]->Fill(
				locKinFitCL,
				locAccidWeight);
			if (locTracker_BeamEnergy[locCutVariant].Check({
					{ Unknown, locBeamID }})) {
				dHist_BeamEnergy[locCutVariant]->Fill(
					locBeamP4.E(),
					locAccidWeight);
			}
			if (locTracker_MissingMass[locCutVariant].Check({
					{ Unknown, locBeamID },
					{ dPi1PID, locPi1TrackID },
					{ Proton, locProtonTrackID },
					{ Gamma, locPhoton1NeutralID },
					{ Gamma, locPhoton2NeutralID }})) {
				dHist_MissingMass[locCutVariant]->Fill(
					locMissingMass,
					locAccidWeight);
				dHist_MissingMass_Measured[locCutVariant]->Fill(
					locMissingMass_Measured,
					locAccidWeight);
			}
			dHist_FoundPi2[locCutVariant]->Fill(
				locFoundPi2[locCutVariant],
				locAccidWeight);
			if (locTracker_Proton[locCutVariant].Check({
					{ Proton, locProtonTrackID }})) {
				dHist_ProtonMomentum[locCutVariant]->Fill(
					locProtonP4.P(),
					locProtonP4.Theta()*RAD_TO_DEG,
					locProtonP4.Phi()*RAD_TO_DEG,
					locAccidWeight);
			}
			if (locTracker_Pi0[locCutVariant].Check({
					{ Gamma, locPhoton1NeutralID },
					{ Gamma, locPhoton2NeutralID }})) {
				dHist_Pi0Momentum[locCutVariant]->Fill(
					locPi0P4.P(),
					locPi0P4.Theta()*RAD_TO_DEG,
					locPi0P4.Phi()*RAD_TO_DEG,
					locAccidWeight);
			}
			if (locTracker_Pi1[locCutVariant].Check({
					{ dPi1PID, locPi1TrackID }})) {
				dHist_Pi1Momentum[locCutVariant]->Fill(
					locPi1P4.P(),
					locPi1P4.Theta()*RAD_TO_DEG,
					locPi1P4.Phi()*RAD_TO_DEG,
					locAccidWeight);
			}
			// TODO: What uniqueness tracking should be done here (true for all
			// `Kinfit` histograms)? Since a kinematic fit involves all
			// particles, should every particle be included in the tracker?
			if (locTracker_Pi2_Kinfit[locCutVariant].Check({})) {
				dHist_Pi2Momentum_Kinfit[locCutVariant]->Fill(
					locPi2P4_Kinfit.P(),
					locPi2P4_Kinfit.Theta()*RAD_TO_DEG,
					locPi2P4_Kinfit.Phi()*RAD_TO_DEG,
					locAccidWeight);
			}
			if (locTracker_Proton3Pi_Kinfit[locCutVariant].Check({
					{ Proton, locProtonTrackID },
					{ Gamma, locPhoton1NeutralID },
					{ Gamma, locPhoton2NeutralID },
					{ dPi1PID, locPi1TrackID }})) {
				dHist_Proton3PiPhi_Kinfit[locCutVariant]->Fill(
					locProtonP4.Phi()*RAD_TO_DEG,
					loc3PiP4_Kinfit.Phi()*RAD_TO_DEG,
					locAccidWeight);
				dHist_Proton3PiCollinearity_Kinfit[locCutVariant]->Fill(
					locProton3PiCollinearity_Kinfit*RAD_TO_DEG,
					locAccidWeight);
				dHist_Proton3PiCollinearity_Mixed[locCutVariant]->Fill(
					locProton3PiCollinearity_Mixed*RAD_TO_DEG,
					locAccidWeight);
			}
			if (locTracker_MMOP[locCutVariant].Check({
					{ Unknown, locBeamID },
					{ Proton, locProtonTrackID }})) {
			}
			if (locFoundPi2[locCutVariant] == 0) {
				if (locTracker_MissingMass_Missing[locCutVariant].Check({
						{ Unknown, locBeamID },
						{ dPi1PID, locPi1TrackID },
						{ Proton, locProtonTrackID },
						{ Gamma, locPhoton1NeutralID },
						{ Gamma, locPhoton2NeutralID }})) {
					dHist_MissingMass_Missing[locCutVariant]->Fill(
						locMissingMass,
						locAccidWeight);
					dHist_MissingMass_Measured_Missing[locCutVariant]->Fill(
						locMissingMass_Measured,
						locAccidWeight);
				}
				// Second pion was not found.
				if (locTracker_Pi0_Missing[locCutVariant].Check({
						{ Gamma, locPhoton1NeutralID },
						{ Gamma, locPhoton2NeutralID }})) {
					dHist_Pi0MassSq_Measured_Missing[locCutVariant]->Fill(
						locPi0P4_Measured.M(),
						locAccidWeight);
					dHist_Pi0MassSq_Missing[locCutVariant]->Fill(
						locPi0P4.M(),
						locAccidWeight);
				}
				if (locTracker_Pi2_Missing[locCutVariant].Check({
						{ Unknown, locBeamID },
						{ dPi1PID, locPi1TrackID },
						{ Proton, locProtonTrackID },
						{ Gamma, locPhoton1NeutralID },
						{ Gamma, locPhoton2NeutralID }})) {
					dHist_Pi2Momentum_Missing[locCutVariant]->Fill(
						locMissingP4.P(),
						locMissingP4.Theta()*RAD_TO_DEG,
						locMissingP4.Phi()*RAD_TO_DEG,
						locAccidWeight);
				}
				if (locTracker_Proton3Pi_Kinfit_Missing[locCutVariant].Check({
						{ Proton, locProtonTrackID },
						{ Gamma, locPhoton1NeutralID },
						{ Gamma, locPhoton2NeutralID },
						{ dPi1PID, locPi1TrackID }})) {
					dHist_Proton3PiCollinearity_Mixed_MMOP_Missing[locCutVariant]->Fill(
						locProton3PiCollinearity_Mixed*RAD_TO_DEG,
						locMMOP,
						locAccidWeight);
				}
				if (locTracker_MMOP_Missing[locCutVariant].Check({
						{ Unknown, locBeamID },
						{ Proton, locProtonTrackID }})) {
					dHist_ProtonMomentum_MMOP_Missing[locCutVariant]->Fill(
						locProtonP4.P(),
						locMMOP,
						locAccidWeight);
					dHist_KinFitCL_MMOP_Missing[locCutVariant]->Fill(
						locKinFitCL,
						locMMOP,
						locAccidWeight);
					dHist_MMOP_Measured_Missing[locCutVariant]->Fill(
						locMMOP_Measured,
						locAccidWeight);
					dHist_MMOP_Missing[locCutVariant]->Fill(
						locMMOP,
						locAccidWeight);
				}
				if (locTracker_M3PI_Kinfit_Missing[locCutVariant].Check({
						{ dPi1PID, locPi1TrackID },
						{ Gamma, locPhoton1NeutralID },
						{ Gamma, locPhoton2NeutralID }})) {
					dHist_M3PI_Kinfit_Missing[locCutVariant]->Fill(
						locM3PI_Kinfit,
						locAccidWeight);
				}
			}
		}

		// FILL FLAT TREE
		if (dFlatTreeInterface != nullptr) {
			dFlatTreeInterface->Fill_Fundamental<Double_t>("accidweight", locAccidWeight);
		}

		/*
		//FILL ANY CUSTOM BRANCHES FIRST!!
		Int_t locMyInt_Flat = 7;
		dFlatTreeInterface->Fill_Fundamental<Int_t>("flat_my_int", locMyInt_Flat);

		TLorentzVector locMyP4_Flat(4.0, 3.0, 2.0, 1.0);
		dFlatTreeInterface->Fill_TObject<TLorentzVector>("flat_my_p4", locMyP4_Flat);

		for(int loc_j = 0; loc_j < locMyInt_Flat; ++loc_j)
		{
			dFlatTreeInterface->Fill_Fundamental<Int_t>("flat_my_int_array", 3*loc_j, loc_j); //2nd argument = value, 3rd = array index
			TLorentzVector locMyComboP4_Flat(8.0, 7.0, 6.0, 5.0);
			dFlatTreeInterface->Fill_TObject<TLorentzVector>("flat_my_p4_array", locMyComboP4_Flat, loc_j);
		}
		*/
		if (dFlatTreeInterface != nullptr) {
			Fill_FlatTree();
		}
	}

	// Record number of combos surviving actions.
	Fill_NumCombosSurvivedHists();

	/******************************************* LOOP OVER THROWN DATA (OPTIONAL) ***************************************/
/*
	//Thrown beam: just use directly
	if(dThrownBeam != NULL)
		Double_t locEnergy = dThrownBeam->Get_P4().E();

	//Loop over throwns
	for(UInt_t loc_i = 0; loc_i < Get_NumThrown(); ++loc_i)
	{
		//Set branch array indices corresponding to this particle
		dThrownWrapper->Set_ArrayIndex(loc_i);

		//Do stuff with the wrapper here ...
	}
*/
	/****************************************** LOOP OVER OTHER ARRAYS (OPTIONAL) ***************************************/
/*
	//Loop over beam particles (note, only those appearing in combos are present)
	for(UInt_t loc_i = 0; loc_i < Get_NumBeam(); ++loc_i)
	{
		//Set branch array indices corresponding to this particle
		dBeamWrapper->Set_ArrayIndex(loc_i);

		//Do stuff with the wrapper here ...
	}

	//Loop over charged track hypotheses
	for(UInt_t loc_i = 0; loc_i < Get_NumChargedHypos(); ++loc_i)
	{
		//Set branch array indices corresponding to this particle
		dChargedHypoWrapper->Set_ArrayIndex(loc_i);

		//Do stuff with the wrapper here ...
	}

	//Loop over neutral particle hypotheses
	for(UInt_t loc_i = 0; loc_i < Get_NumNeutralHypos(); ++loc_i)
	{
		//Set branch array indices corresponding to this particle
		dNeutralHypoWrapper->Set_ArrayIndex(loc_i);

		//Do stuff with the wrapper here ...
	}
*/

	/************************************ EXAMPLE: FILL CLONE OF TTREE HERE WITH CUTS APPLIED ************************************/
/*
	Bool_t locIsEventCut = true;
	for(UInt_t loc_i = 0; loc_i < Get_NumCombos(); ++loc_i) {
		//Set branch array indices for combo and all combo particles
		dComboWrapper->Set_ComboIndex(loc_i);
		// Is used to indicate when combos have been cut
		if(dComboWrapper->Get_IsComboCut())
			continue;
		locIsEventCut = false; // At least one combo succeeded
		break;
	}
	if(!locIsEventCut && dOutputTreeFileName != "")
		Fill_OutputTree();
*/

	return kTRUE;
}

void DSelector_omega_pi_effic::Finalize(void) {
	// Save anything to output here that you do not want to be in the default
	// DSelector output ROOT file.

	// Otherwise, don't do anything else (especially if you are using PROOF).
	// If you are using PROOF, this function is called on each thread, so
	// anything you do will not have the combined information from the various
	// threads.  Besides, it is best-practice to do post-processing (e.g.
	// fitting) separately, in case there is a problem.

	// Call this after everything else. Saves results to the output file.
	DSelector::Finalize();
}

