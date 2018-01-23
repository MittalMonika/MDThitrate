
#include <iostream>
#include <string>
#include "TChain.h"
#include "MuTrigNtExampleAna/ExampleLooper.h"
#include "PhaseIIMuonTriggerNtuple/ChainHelper.h"
#include "PhaseIIMuonTriggerNtuple/string_utils.h"


void printHelp() {
    std::cout << "  Options:\n"
              << "  -n number of events to process\n"
              << "     defaults: -1 (all events)\n"

              << "  -k number of events to skip\n"
              << "     defaults: 0\n"

              << "  -d debug printout level\n"
              << "     defaults: 0 (quiet)\n"

              << "  -i input (file, list, or dir)\n"
              << "     defaults: ''\n"

              << "  -o output file \n"
              << "     defaults: ''\n"
 
              << "  -s sample name, for naming files\n"
              << "     defaults: ntuple sample name\n"


              << "  -h print this help" << std::endl;
}

int main(int argc, char **argv) {
    int nEvt = -1;
    int nSkip = 0;
    int dbg = 0;
    std::string sample;
    std::string dataset  = "zb";
    std::string input, output;
    std::cout << "run_example_looper\n\n";

    /** Read inputs to program */
    for (int i = 1; i < argc; ++i) {
        if      (strcmp(argv[i], "-n") == 0) nEvt = atoi(argv[++i]);
        else if (strcmp(argv[i], "-k") == 0) nSkip = atoi(argv[++i]);
        else if (strcmp(argv[i], "-d") == 0) dbg = atoi(argv[++i]);
        else if (strcmp(argv[i], "-i") == 0) input = argv[++i];
	else if (strcmp(argv[i], "-o") == 0) output = argv[++i]; 
        else if (strcmp(argv[i], "-s") == 0) sample = argv[++i];
	else if (strcmp(argv[i], "-c") == 0) dataset = argv[++i];
        else {
            printHelp();
            return 0;
        }
    }

    if (input.empty()) {
        std::cout << "You must specify an input\n";
        return 1;
    }

    if (dbg) {
        std::cout << "Being called as: " << TrigNtup::utils::commandLineArguments(argc, argv) << "\n";
    }

    std::cout << "flags:\n"
              << "  sample  " << sample << "\n"
              << "  nEvt    " << nEvt   << "\n"
              << "  nSkip   " << nSkip  << "\n"
              << "  dbg     " << dbg    << "\n"
              << "  input   " << input  << "\n"
              << "  output   " << output  << "\n"
	      << " dataset "<< dataset << "\n"
              << std::endl;

    bool verbose = dbg > 0;
    // Build the input chain
    TChain chain {"MuonTriggerNt"};
    ChainHelper::addInput(&chain, input, verbose);
    Long64_t nEntries = chain.GetEntries();
    chain.ls();

    // Build the TSelector
    //ExampleLooper looper {&chain};
    ExampleLooper looper (&chain, output, dataset.c_str());
    looper.setDebug(dbg);
    //looper.setSampleName(ChainHelper::sampleName(input, verbose));

    // Run the job
    if (nEvt < 0) nEvt = nEntries;
    std::cout << "\n"
              << "Total entries:   " << nEntries << "\n"
              << "Process entries: " << nEvt     << "\n";
    if (nEvt > 0) chain.Process(&looper, sample.c_str(), nEvt, nSkip);

    std::cout << "\nfinished\n";

    return 0;
}

