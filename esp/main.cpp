#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <cassert>

extern void ReadSim(char *filename, int nInputs, char **pValues, int &nBpatterns);
extern void GeneratePLA(std::string filename, std::vector<std::vector<int> > onset, std::vector<char *> pValues, int nBPatterns, int rarity);

int main(int argc, char **argv) {
  char *filename = argv[1];
  char *simname = argv[2];
  std::string filename2 = filename;
  filename2 += ".opt.blif";

  int nGroupOutputs  = 3;
  
  std::ifstream f(filename);
  std::ofstream f2(filename2);
  std::string str;
  
  std::vector<std::string> inputs;
  int nInputs;
  while(std::getline(f, str)) {
    f2 << str << std::endl;
    if(str.length() > 7 && str.substr(0, 7) == ".inputs") {
      break;
    }
  }
  {
    std::stringstream ss(str);
    std::string input;
    std::getline(ss, input, ' ');
    while(std::getline(ss, input, ' ')) {
      inputs.push_back(input);
    }
  }
  nInputs = inputs.size();

  char **pValues = (char **)malloc(sizeof(char *) * nInputs);
  int nBPatterns;
  ReadSim(simname, nInputs, pValues, nBPatterns);

  std::map<std::string, int> input2index;
  for(int i = 0; i < inputs.size(); i++) {
    input2index[inputs[i]] = i;
  }
  
  std::getline(f, str);
  while(f) {
    if(str.length() < 6 || str.substr(0, 6) != ".names") {
      f2 << str << std::endl;
      std::getline(f, str);
      continue;
    }
    std::vector<std::string> LUTInputs;
    std::vector<std::string> LUTOutputs;
    std::vector<std::vector<int> > onsets;
    for(int j = 0; j < nGroupOutputs; j++) {
      std::vector<std::string> LUTInputs_;
      std::string LUTOutput;
      {
        std::cout << str << std::endl;
        std::stringstream ss(str);
        std::string input;
        std::getline(ss, input, ' ');
        while(std::getline(ss, input, ' ')) {
          LUTInputs_.push_back(input);
        }
        LUTOutput = LUTInputs_.back();
        LUTInputs_.pop_back();
      }
      LUTOutputs.push_back(LUTOutput);
      if(j == 0) {
        LUTInputs = LUTInputs_;
      } else {
        assert(LUTInputs == LUTInputs_);
      }
      std::vector<int> onset;
      while(std::getline(f, str)) {
        if(str.length() == 0 || (str[0] != '0' && str[0] != '1')) {
          break;
        }
        int value = 0;
        for(int i = 0; i < LUTInputs.size(); i++) {
          value = (str[i] - '0') + (value << 1);
        }
        onset.push_back(value);
      }
      onsets.push_back(onset);
      if(j == nGroupOutputs - 1) {
        break;
      }
      while(str.length() < 6 || str.substr(0, 6) != ".names") {
        f2 << str << std::endl;    
        std::getline(f, str);
      }
    }
    std::vector<char *> pValuesSubset;
    for(int i = 0; i < LUTInputs.size(); i++) {
      pValuesSubset.push_back(pValues[input2index[LUTInputs[i]]]);
    }
      
    std::string planame = "test.pla";
    int rarity = 1;
    GeneratePLA(planame, onsets, pValuesSubset, nBPatterns, rarity);

    std::string planame2 = planame + ".esp.pla";
    std::string cmd = "espresso " + planame + " > " + planame2;
    std::system(cmd.c_str());

    std::vector<std::vector<std::string> > optimized_onsets(nGroupOutputs);
    std::ifstream f3(planame2);
    std::string str2;
    while(std::getline(f3, str2)) {
      if(str2.length() > 0 && (str2[0] == '0' || str2[0] == '1' || str2[0] == '-')) {
        int pos = str2.find(' ');
        std::string value = str2.substr(0, pos);
        for(int j = 0; j < nGroupOutputs; j++) {
          if(str2[pos+1+j] == '1') {
            optimized_onsets[j].push_back(value);
          }
        }
      }
    }
    f3.close();
    
    for(int j = 0; j < nGroupOutputs; j++) {
      f2 << ".names";
      for(auto input: LUTInputs) {
        f2 << " " << input;
      }
      f2 << " " << LUTOutputs[j] << std::endl;
      for(auto onset: optimized_onsets[j]) {
        f2 << onset << " 1" << std::endl;
      }
    }
  }
    
  for(int i = 0; i < nInputs; i++) {
    free(pValues[i]);
  }
  free(pValues);
  

  return 0;
  /*
  for(int i = 0; i < 4; i++) {
    int c = pValues[1][i];
    for(int j = 0; j < 8; j++) {
      if(c & 1) {
        printf("1\n");
      } else {
        printf("0\n");
      }
      c = c >> 1;
    }
  }


  return 0;
  */
}
