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
extern void ReadSim(std::string filename, int nInputs, std::vector<char *> &vpBPats, int &nBpatterns);
extern int ReadBlifFuncs(std::ifstream &f, int nGroupSize, std::vector<std::string> &LUTInputs, std::vector<std::string> &LUTOutputs, std::vector<std::vector<int> > &onsets);
extern void GeneratePla(std::string filename, std::vector<std::vector<int> > const &onsets, std::vector<char *> const &vpBPats, int nBPats, int rarity);
extern void ReadPla(std::string filename, std::vector<std::vector<std::string> > &onsets);

extern void TTTest(std::vector<std::vector<int> > const &onsets, std::vector<char *> const &vpBPats, int nBPats, int rarity, std::vector<std::string> const &inputs, std::vector<std::string> const &outputs, std::ofstream &f);

extern void BddTest(std::vector<std::vector<int> > const &onsets, std::vector<std::string> const &inputs, std::vector<std::string> const &outputs, std::ofstream &f);

void RunEspresso(std::vector<std::vector<int> > const &onsets, std::vector<char *> const &vpBPats, int nBPats, int rarity, std::vector<std::string> const &inputs, std::vector<std::string> const &outputs, std::ofstream &f) {
  std::string planame = "test.pla";
  GeneratePla(planame, onsets, vpBPats, nBPats, rarity);
  std::string planame2 = planame + ".esp.pla";
  std::string cmd = "espresso " + planame + " > " + planame2;
  int r = std::system(cmd.c_str());
  assert(r == 0);
  std::vector<std::vector<std::string> > onsets2;
  ReadPla(planame2, onsets2);
  for(uint i = 0; i < outputs.size(); i++) {
    f << ".names";
    for(auto input: inputs) {
      f << " " << input;
    }
    f << " " << outputs[i] << std::endl;
    for(auto pat: onsets2[i]) {
      f << pat << " 1" << std::endl;
    }
  }
}

int main(int argc, char **argv) {
  std::string ifname = argv[1];
  std::string ofname = argv[2];
  std::string simname = argv[3];
  int nGroupSize = std::stoi(argv[4]);
  int rarity = std::stoi(argv[5]);

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

  std::vector<char *> vpBPats;
  int nBPats;
  if(!simname.empty()) {
    ReadSim(simname, nInputs, vpBPats, nBPats);
  }

  std::vector<std::string> LUTInputs;
  std::vector<std::string> LUTOutputs;
  std::vector<std::vector<int> > onsets;
  while(ReadBlifFuncs(f, nGroupSize, LUTInputs, LUTOutputs, onsets)) {
    std::vector<char *> vpBPatsSubset(LUTInputs.size());
    if(!simname.empty()) {
      for(uint i = 0; i < LUTInputs.size(); i++) {
        vpBPatsSubset[i] = vpBPats[input2index[LUTInputs[i]]];
      }
    }

#ifdef REORDER_WITH_CUDD
    BddTest(onsets, LUTInputs, LUTOutputs, of);
#else
    TTTest(onsets, vpBPatsSubset, nBPats, rarity, LUTInputs, LUTOutputs, of);
#endif
    //RunEspresso(onsets, vpBPatsSubset, nBPats, rarity, LUTInputs, LUTOutputs, of);
  }

  of << ".end" << std::endl;

  for(auto pBPats: vpBPats) {
    free(pBPats);
  }
  vpBPats.clear();

  return 0;
}
