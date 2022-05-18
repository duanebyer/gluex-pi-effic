R__ADD_INCLUDE_PATH($HALLD_RECON_HOME/$BMS_OSNAME/include/)
R__ADD_INCLUDE_PATH($ROOT_ANALYSIS_HOME/$BMS_OSNAME/include/)
R__LOAD_LIBRARY($ROOT_ANALYSIS_HOME/$BMS_OSNAME/lib/libDSelector.so)

void run_simple(char const* config_name) {
	ifstream config(config_name);
	std::string tree_name;
	config >> tree_name;
	TChain* chain = new TChain(tree_name.c_str());
	while (config) {
		std::string next_file;
		config >> next_file;
		if (next_file != "") {
			chain->Add(next_file.c_str());
		}
	}
	config.close();
	gROOT->ProcessLine(".x $ROOT_ANALYSIS_HOME/scripts/Load_DSelector.C");
	chain->Process("DSelector_omega_pi_effic.C++g");
}

