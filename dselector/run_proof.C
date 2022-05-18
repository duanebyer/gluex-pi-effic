R__ADD_INCLUDE_PATH($HALLD_RECON_HOME/$BMS_OSNAME/include/)
R__ADD_INCLUDE_PATH($ROOT_ANALYSIS_HOME/$BMS_OSNAME/include/)
R__LOAD_LIBRARY($ROOT_ANALYSIS_HOME/$BMS_OSNAME/lib/libDSelector.so)

void run_proof(
		char const* in_name,
		char const* out_name="omega_pi_effic.root",
		char const* tree_name="pi0pipmisspim__B1_T1_U1_M7_Effic_Tree",
		UInt_t num_threads=4) {
	gROOT->ProcessLine(".x $(ROOT_ANALYSIS_HOME)/scripts/Load_DSelector.C");
	std::cout << "File: " << in_name << "\n";
	std::cout << "Tree name: " << tree_name << "\n";
	TChain* chain = new TChain(tree_name);
	chain->Add(in_name);
	std::cout << "# of trees: " << chain->GetNtrees() << "\n";
	std::cout << "# of entries: " << chain->GetEntries() << "\n";
	//chain->Process("DSelector_omega_pi_effic.C++");
	DPROOFLiteManager::Process_Chain(
		chain,
		"DSelector_omega_pi_effic.C+g",
		num_threads,
		out_name,
		"", "");
}

