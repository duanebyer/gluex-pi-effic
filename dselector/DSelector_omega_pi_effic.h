#ifndef DSelector_omega_pi_effic_h
#define DSelector_omega_pi_effic_h

#include <iostream>

#include <DSelector/DSelector.h>
#include <DSelector/DHistogramActions.h>
#include <DSelector/DCutActions.h>

#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>

UInt_t const CUT_COUNT = 15;
UInt_t const CUT_FOUND_COUNT = 4;
struct Cut {
	char const* dName;
	Bool_t dEnabled[CUT_COUNT];
	Bool_t dEnabledFound[CUT_FOUND_COUNT];
};
// We include the following variants:
//  * All cuts applied
//  * No cuts applied
//  * All cuts applied individually
//  * All but one cut applied
Cut const CUT_VARIANT[] = {
	{ "Std",          { 1,1,1,1,1,0,0,0,0,0,0,0,0,0,0 }, { 1,1,1,0 } },
	{ "Std+Diag",     { 1,1,1,1,1,0,0,0,0,0,0,0,0,1,0 }, { 1,1,1,0 } },
	{ "Std+Coll1",    { 1,1,1,1,1,0,0,0,0,0,0,0,0,0,1 }, { 1,1,1,0 } },
	{ "Std+Coll2",    { 1,1,1,1,1,0,0,0,0,0,0,0,0,0,0 }, { 1,1,1,1 } },
	{ "Std+Pmom>0.4", { 1,1,1,1,1,1,0,0,0,0,0,0,0,0,0 }, { 1,1,1,0 } },
	{ "Std+Pmom>0.5", { 1,1,1,1,1,0,1,0,0,0,0,0,0,0,0 }, { 1,1,1,0 } },
	{ "Std+Pmom>0.6", { 1,1,1,1,1,0,0,1,0,0,0,0,0,0,0 }, { 1,1,1,0 } },
	{ "Std+Pmom>0.7", { 1,1,1,1,1,0,0,0,1,0,0,0,0,0,0 }, { 1,1,1,0 } },
	{ "Std+Pmom>0.8", { 1,1,1,1,1,0,0,0,0,1,0,0,0,0,0 }, { 1,1,1,0 } },
	{ "Std+CL>0.2",   { 1,1,1,1,1,0,0,0,0,0,1,0,0,0,0 }, { 1,1,1,0 } },
	{ "Std+CL>0.3",   { 1,1,1,1,1,0,0,0,0,0,0,1,0,0,0 }, { 1,1,1,0 } },
	{ "Std+CL>0.4",   { 1,1,1,1,1,0,0,0,0,0,0,0,1,0,0 }, { 1,1,1,0 } },
	{ "None",         { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 }, { 0,0,0,0 } },
};
UInt_t const CUT_VARIANT_COUNT = sizeof(CUT_VARIANT) / sizeof(CUT_VARIANT[0]);

enum class AccidMethod {
	NONE,
	SUBTRACT_ALL,
	SUBTRACT_ALL_BUT_NEAREST,
	SRC2021,
};

enum class Polarization {
	PARA,
	PERP,
	AMO,
};

struct AccidScalingFactor {
	Double_t hodoscopeHiFactor;
	Double_t hodoscopeLoFactor;
	Double_t microscopeFactor;
	Double_t TAGMEnergyBoundHi;
	Double_t TAGMEnergyBoundLo;
	AccidScalingFactor() :
		hodoscopeHiFactor(1.),
		hodoscopeLoFactor(1.),
		microscopeFactor(1.),
		TAGMEnergyBoundHi(8.9131731046),
		TAGMEnergyBoundLo(7.8887606296) { }
	Double_t GetFactor(Double_t beamEnergy) const {
		if (beamEnergy > TAGMEnergyBoundHi) {
			return hodoscopeHiFactor;
		} else if (beamEnergy > TAGMEnergyBoundLo) {
			return microscopeFactor;
		} else {
			return hodoscopeLoFactor;
		}
	}
};

class DSelector_omega_pi_effic : public DSelector {
	public:

		DSelector_omega_pi_effic(TTree* locTree = NULL) : DSelector(locTree){}
		virtual ~DSelector_omega_pi_effic(){}

		void Init(TTree *tree);
		Bool_t Process(Long64_t entry);

	private:

		void Get_ComboWrappers(void);
		void Finalize(void);

		TString dInputTreeName;
		Int_t dNumBeamBunches;     // param 'B'
		Int_t dFitType;            // param 'F'
		Int_t dUnusedTracks;       // param 'U'
		Int_t dExtraChargedTracks; // param 'T'
		Int_t dMassConstraint;     // param 'M'

		Bool_t dUseDefaults;
		AccidMethod dAccidMethod;
		UInt_t dPreviousRunNumber;
		Polarization dPolarization;
		Double_t dBeamBunchSpacing;
		AccidScalingFactor dAccidScalingFactor;

		Bool_t dIsMC;

		// ANALYZE CUT ACTIONS
		// // Automatically makes mass histograms where one cut is missing
		DHistogramAction_AnalyzeCutActions* dAnalyzeCutActions;

		//CREATE REACTION-SPECIFIC PARTICLE ARRAYS

		//Step 0
		DParticleComboStep* dStep0Wrapper;
		DBeamParticle* dComboBeamWrapper;
		DChargedTrackHypothesis* dPi1Wrapper;
		DChargedTrackHypothesis* dProtonWrapper;
		DChargedTrackHypothesis* dPi2Wrapper;
		Particle_t dPi1PID;
		Particle_t dPi2PID;

		//Step 1
		DParticleComboStep* dStep1Wrapper;
		DNeutralParticleHypothesis* dPhoton1Wrapper;
		DNeutralParticleHypothesis* dPhoton2Wrapper;

		// Histograms.
		// Some general naming conventions:
		//  * By default, quantities are computed with momenta from the
		//    kinematic fit.
		//  * `Measured` refers to quantities computed with momenta from tracks.
		//  * `KinFit` refers to when the missing pion momentum is also computed
		//    from the kinematic fit.
		TH1D* dHist_KinFitCL[CUT_VARIANT_COUNT];
		TH1D* dHist_MissingMass[CUT_VARIANT_COUNT];
		TH1D* dHist_MissingMass_Found[CUT_VARIANT_COUNT];
		TH1D* dHist_MissingMass_Missing[CUT_VARIANT_COUNT];
		TH1D* dHist_MissingMass_Measured[CUT_VARIANT_COUNT];
		TH1D* dHist_MissingMass_Measured_Found[CUT_VARIANT_COUNT];
		TH1D* dHist_MissingMass_Measured_Missing[CUT_VARIANT_COUNT];
		TH1D* dHist_BeamEnergy[CUT_VARIANT_COUNT];
		TH1D* dHist_FoundPi2[CUT_VARIANT_COUNT];

		TH3D* dHist_ProtonMomentum[CUT_VARIANT_COUNT];
		TH3D* dHist_Pi0Momentum[CUT_VARIANT_COUNT];
		TH3D* dHist_Pi1Momentum[CUT_VARIANT_COUNT];
		TH3D* dHist_Pi2Momentum_Kinfit[CUT_VARIANT_COUNT];
		TH3D* dHist_Pi2Momentum_Found[CUT_VARIANT_COUNT];
		TH3D* dHist_Pi2Momentum_Missing[CUT_VARIANT_COUNT];
		// Coplanarity of proton and pion.
		TH2D* dHist_Proton3PiPhi_Kinfit[CUT_VARIANT_COUNT];
		TH1D* dHist_Proton3PiCollinearity_Kinfit[CUT_VARIANT_COUNT];
		TH2D* dHist_Proton3PiPhi_Found[CUT_VARIANT_COUNT];
		TH1D* dHist_Proton3PiCollinearity_Found[CUT_VARIANT_COUNT];
		TH1D* dHist_Proton3PiCollinearity_Mixed[CUT_VARIANT_COUNT];
		// MMOP diagnostics.
		TH2D* dHist_ProtonMomentum_MMOP_Found[CUT_VARIANT_COUNT];
		TH2D* dHist_ProtonMomentum_MMOP_Missing[CUT_VARIANT_COUNT];
		TH2D* dHist_Proton3PiCollinearity_Mixed_MMOP_Found[CUT_VARIANT_COUNT];
		TH2D* dHist_Proton3PiCollinearity_Mixed_MMOP_Missing[CUT_VARIANT_COUNT];
		TH2D* dHist_Proton3PiCollinearity_MMOP_Found[CUT_VARIANT_COUNT];
		TH2D* dHist_KinFitCL_MMOP_Missing[CUT_VARIANT_COUNT];
		TH2D* dHist_KinFitCL_MMOP_Found[CUT_VARIANT_COUNT];
		// Invariant mass of neutral pion.
		TH1D* dHist_Pi0MassSq_Measured_Found[CUT_VARIANT_COUNT];
		TH1D* dHist_Pi0MassSq_Measured_Missing[CUT_VARIANT_COUNT];
		TH1D* dHist_Pi0MassSq_Found[CUT_VARIANT_COUNT];
		TH1D* dHist_Pi0MassSq_Missing[CUT_VARIANT_COUNT];
		// Missing mass off proton.
		TH1D* dHist_MMOP_Measured_Found[CUT_VARIANT_COUNT];
		TH1D* dHist_MMOP_Measured_Missing[CUT_VARIANT_COUNT];
		TH1D* dHist_MMOP_Found[CUT_VARIANT_COUNT];
		TH1D* dHist_MMOP_Missing[CUT_VARIANT_COUNT];
		// Mass of three pions.
		TH1D* dHist_M3PI_Measured_Found[CUT_VARIANT_COUNT];
		TH1D* dHist_M3PI_Kinfit_Found[CUT_VARIANT_COUNT];
		TH1D* dHist_M3PI_Kinfit_Missing[CUT_VARIANT_COUNT];

	ClassDef(DSelector_omega_pi_effic, 0);
};

void DSelector_omega_pi_effic::Get_ComboWrappers() {
	//Step 0
	dStep0Wrapper = dComboWrapper->Get_ParticleComboStep(0);
	dComboBeamWrapper = static_cast<DBeamParticle*>(dStep0Wrapper->Get_InitialParticle());
	dPi1Wrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(1));
	dProtonWrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(2));
	dPi2Wrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(3));
	dPi1PID = dPi1Wrapper->Get_PID();
	dPi2PID = dPi2Wrapper->Get_PID();

	//Step 1
	dStep1Wrapper = dComboWrapper->Get_ParticleComboStep(1);
	dPhoton1Wrapper = static_cast<DNeutralParticleHypothesis*>(dStep1Wrapper->Get_FinalParticle(0));
	dPhoton2Wrapper = static_cast<DNeutralParticleHypothesis*>(dStep1Wrapper->Get_FinalParticle(1));
}

#endif

