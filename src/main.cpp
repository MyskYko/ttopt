#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <cassert>

extern void ReadBlifHeader(std::ifstream &f, std::string &modulename, std::vector<std::string> &inputs, std::vector<std::string> &outputs);
extern void ReadSim(std::string filename, int nInputs, char **pBPats, int &nBpatterns);
extern int ReadBlifFuncs(std::ifstream &f, int nGroupSize, std::vector<std::string> &LUTInputs, std::vector<std::string> &LUTOutputs, std::vector<std::vector<int> > &onsets);

extern void GeneratePla(std::string filename, std::vector<std::vector<int> > onsets, std::vector<char *> pBPats, int nBPats, int rarity);
extern void ReadPla(std::string filename, std::vector<std::vector<std::string> > &onsets);

void RunEspresso(std::vector<std::vector<int> > onsets, std::vector<char *> pBPats, int nBPats, int rarity, std::vector<std::vector<std::string> > &optimized_onsets) {
  std::string planame = "test.pla";
  GeneratePla(planame, onsets, pBPats, nBPats, rarity);
  std::string planame2 = planame + ".esp.pla";
  std::string cmd = "espresso " + planame + " > " + planame2;
  int r = std::system(cmd.c_str());
  assert(r == 0);
  ReadPla(planame2, optimized_onsets);
}

int main(int argc, char **argv) {
  std::string ifname = argv[1];
  std::string simname = argv[2];
  std::string ofname = ifname + ".opt.blif";
  int nGroupSize  = 3;
  int rarity = 1;

  std::ifstream f(ifname);
  std::ofstream of(ofname);

  std::string modulename;
  std::vector<std::string> inputs, outputs;
  ReadBlifHeader(f, modulename, inputs, outputs);
  int nInputs = inputs.size();
  std::map<std::string, int> input2index;
  for(uint i = 0; i < inputs.size(); i++) {
    input2index[inputs[i]] = i;
  }

  of << ".model " << modulename << std::endl;
  of << ".inputs";
  for(auto input: inputs) {
    of << " " << input;
  }
  of << std::endl;
  of << ".outputs";
  for(auto output: outputs) {
    of << " " << output;
  }
  of << std::endl;

  char **pBPats = (char **)malloc(sizeof(char *) * nInputs);
  int nBPats;
  ReadSim(simname, nInputs, pBPats, nBPats);

  std::vector<std::string> LUTInputs;
  std::vector<std::string> LUTOutputs;
  std::vector<std::vector<int> > onsets;
  while(ReadBlifFuncs(f, nGroupSize, LUTInputs, LUTOutputs, onsets)) {
    std::vector<char *> pBPatsSubset;
    for(uint i = 0; i < LUTInputs.size(); i++) {
      pBPatsSubset.push_back(pBPats[input2index[LUTInputs[i]]]);
    }

    std::vector<std::vector<std::string> > optimized_onsets;
    RunEspresso(onsets, pBPatsSubset, nBPats, rarity, optimized_onsets);
    
    for(int j = 0; j < nGroupSize; j++) {
      of << ".names";
      for(auto input: LUTInputs) {
        of << " " << input;
      }
      of << " " << LUTOutputs[j] << std::endl;
      for(auto onset: optimized_onsets[j]) {
        of << onset << " 1" << std::endl;
      }
    }
  }

  of << ".end" << std::endl;

  for(int i = 0; i < nInputs; i++) {
    free(pBPats[i]);
  }
  free(pBPats);

  return 0;
}
