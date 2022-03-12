#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <list>
#include <cassert>

extern std::string BinaryToString(int bin, int size);

class TT {
public:
  int nInputs;
  int nSize;
  int nOutputs;
  std::vector<uint> t;

  static const uint ones[];

  TT(std::vector<std::vector<int> > const &onsets, int nInputs): nInputs(nInputs) {
    assert(nInputs >= 5);
    nSize = 1 << (nInputs - 5);
    nOutputs = onsets.size();
    t.resize(nSize * nOutputs);
    for(int i = 0; i < nOutputs; i++) {
      for(int pat: onsets[i]) {
        int index = pat / 32;
        int pos = pat % 32;
        t[index + nSize * i] |= 1 << pos;
      }
    }
  }
  
  void GeneratePla(std::string filename) {
    std::ofstream f(filename);
    f << ".i " << nInputs << std::endl;
    f << ".o " << nOutputs << std::endl;
    for(int index = 0; index < nSize; index++) {
      for(int pos = 0; pos < 32; pos++) {
        int pat = (index << 5) + pos;
        f << BinaryToString(pat, nInputs) << " ";
        for(int i = 0; i < nOutputs; i++) {
          f << (char)('0' + ((t[index + nSize * i] >> pos) & 1));
        }
        f << std::endl;
      }
    }
  }

  int CountBDDNodesOneGetValue(std::vector<int> &vNodes, int node, int input) {
    assert(nInputs - input <= 5);
    int index = node >> (input + 5 - nInputs);
    if(nInputs - input == 5) {
      return t[index];
    }
    int pos = (node % (1 << (input + 5 - nInputs))) << (nInputs - input);
    return (t[index] >> pos) & ones[nInputs - input];
  }

  void CountBDDNodesOne(std::vector<int> &vNodes, int node, int input) {
    if(nInputs - input > 5) {
      int nScope = 1 << (nInputs - input - 5);
      bool fZero = true;
      bool fOne = true;
      uint one = ones[5];
      std::list<int> nodes;
      for(int node2: vNodes) {
        nodes.push_back(node2 << 1);
        nodes.push_back((node2 << 1) ^ 1);
      }
      for(int i = 0; i < nScope; i++) {
        fZero &= t[nScope * node + i] == 0;
        fOne &= t[nScope * node + i] == one;
        for(auto it = nodes.begin(); it != nodes.end(); it++) {
          int node2 = *it >> 1;
          if(*it & 1) {
            if((t[nScope * node + i] ^ one) != t[nScope * node2 + i]) {
              it = nodes.erase(it);
              it--;
            }
          } else {
            if(t[nScope * node + i] != t[nScope * node2 + i]) {
              it = nodes.erase(it);
              it--;
            }
          }
        }
        if(!fZero && !fOne && nodes.empty()) {
          break;
        }
      }
      if(!fZero && !fOne && nodes.empty()) {
        vNodes.push_back(node);
      }
    } else {
      uint one = ones[nInputs - input];
      uint value = CountBDDNodesOneGetValue(vNodes, node, input);
      if(value == 0 || value == one) {
        return;
      }
      for(int node2: vNodes) {
        uint value2 = CountBDDNodesOneGetValue(vNodes, node2, input);
        if(value2 == value || (value2 ^ one) == value) {
          return;
        }
      }
      vNodes.push_back(node);
    }
  }
  
  int CountBDDNodes() {
    std::vector<std::vector<int> > vvNodes(nInputs);
    for(int i = 0; i < nOutputs; i++) {
      CountBDDNodesOne(vvNodes[0], i, 0);
    }
    for(int i = 1; i < nInputs; i++) {
      for(int node: vvNodes[i-1]) {
        CountBDDNodesOne(vvNodes[i], node << 1, i);
        CountBDDNodesOne(vvNodes[i], (node << 1) + 1, i);
      }
    }
    int count = 0;
    for(auto vNodes: vvNodes) {
      count += vNodes.size();
    }
    return count;
  }
};

const uint TT::ones[] = {0,
                0x00000003,
                0x0000000f,
                0x000000ff,
                0x0000ffff,
                0xffffffff};

void TTTest(std::vector<std::vector<int> > const &onsets, std::vector<char *> const &pBPats, int nBPats, int rarity, std::vector<std::vector<std::string> > &optimized_onsets) {
  int nInputs = pBPats.size();
  TT tt(onsets, nInputs);
  //tt.GeneratePla("test2.pla");
  std::cout << tt.CountBDDNodes() << std::endl;
}
