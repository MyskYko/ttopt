#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <iomanip>

extern std::string BinaryToString(int bin, int size);

void GeneratePla(std::string filename, std::vector<std::vector<int> > const &onsets, std::vector<char *> const &pBPats, int nBPats, int rarity) {
  int LUTSize = pBPats.size();
  int nOutputs = onsets.size();
  std::vector<int> careset;
  if(rarity == 0) {
    for(int pat = 0; pat < (1 << LUTSize); pat++) {
      careset.push_back(pat);
    }
  } else {
    std::vector<int> count(1 << LUTSize);
    for(int i = 0; i < nBPats; i++) {
      for(int j = 0; j < 8; j++) {
        std::string str;
        for(auto pBPat: pBPats) {
          str += ((pBPat[i] >> j) & 1) + '0';
        }
        int pat = std::stoi(str, 0, 2);
        count[pat]++;
        if(count[pat] == rarity) {
          careset.push_back(pat);
        }
      }
    }
  }
  
  std::vector<int> count(1 << LUTSize);
  for(int i = 0; i < nOutputs; i++) {
    int pow10 = 1;
    for(int j = i; j < nOutputs-1; j++) {
      pow10 = pow10 * 10;
    }
    for(int pat: onsets[i]) {
      count[pat] += pow10;
    }
  }

  std::ofstream f(filename);
  f << ".i " << LUTSize << std::endl;
  f << ".o " << nOutputs << std::endl;
  f << ".type fr" << std::endl;
  for(int pat: careset) {
    f << BinaryToString(pat, LUTSize) << " " << std::setfill('0') << std::setw(nOutputs) << count[pat] << std::endl;
  }

  f.close();
}

void ReadPla(std::string filename, std::vector<std::vector<std::string> > &onsets) {
  onsets.clear();
  std::ifstream f(filename);
  std::string str;
  while(std::getline(f, str)) {
    if(str.length() > 2 && str.substr(0, 2) == ".i") {
      break;
    }
  }
  int nInputs = std::stoi(str.substr(3));
  while(std::getline(f, str)) {
    if(str.length() > 2 && str.substr(0, 2) == ".o") {
      break;
    }
  }
  int nOutputs = std::stoi(str.substr(3));
  onsets.resize(nOutputs);
  while(std::getline(f, str)) {
    if(str.length() > 0 && (str[0] == '0' || str[0] == '1' || str[0] == '-')) {
      std::string pat = str.substr(0, nInputs);
      for(int j = 0; j < nOutputs; j++) {
        if(str[nInputs+1+j] == '1') {
          onsets[j].push_back(pat);
        }
      }
    }
  }
  f.close();
}
