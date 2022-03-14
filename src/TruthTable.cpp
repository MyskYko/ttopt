#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <list>
#include <algorithm>
#include <numeric>
#include <cassert>

extern std::string BinaryToString(int bin, int size);

class TT {
public:
  int nInputs;
  int nSize;
  int nOutputs;
  std::vector<uint> t;
  std::vector<int> vLevels;

  static const uint ones[];
  static const uint sieve[];

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
    vLevels.resize(nInputs);
    std::iota(vLevels.begin(), vLevels.end(), 0);
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

  int BDDGetValue(int node, int lev) {
    assert(nInputs - lev <= 5);
    int index = node >> (lev + 5 - nInputs);
    int pos = (node % (1 << (lev + 5 - nInputs))) << (nInputs - lev);
    return (t[index] >> pos) & ones[nInputs - lev];
  }

  int CountBDDNodesOne(std::vector<int> &vNodes, int node, int lev) {
    if(nInputs - lev > 5) {
      int nScope = 1 << (nInputs - lev - 5);
      bool fZero = true;
      bool fOne = true;
      uint one = ones[5];
      // vertical procedure begin
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
            if(~t[nScope * node + i] != t[nScope * node2 + i]) {
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
      if(fZero) {
        return -1;
      }
      if(fOne) {
        return -2;
      }
      if(!nodes.empty()) {
        return nodes.front();
      }
      // vertical procedure end
    } else {
      uint one = ones[nInputs - lev];
      uint value = BDDGetValue(node, lev);
      if(value == 0) {
        return -1;
      }
      if(value == one) {
        return -2;
      }
      for(int node2: vNodes) {
        uint value2 = BDDGetValue(node2, lev);
        if(value2 == value) {
          return node2 << 1;
        }
        if((value2 ^ one) == value) {
          return (node2 << 1) ^ 1;
        }
      }
    }
    vNodes.push_back(node);
    return node << 1;
  }
  
  int CountBDDNodes() {
    std::vector<std::vector<int> > vvNodes(nInputs);
    std::vector<std::vector<int> > vvNodesRedundant(nInputs);
    for(int i = 0; i < nOutputs; i++) {
      CountBDDNodesOne(vvNodes[0], i, 0);
    }
    for(int i = 1; i < nInputs; i++) {
      for(int node: vvNodes[i-1]) {
        int cof0 = CountBDDNodesOne(vvNodes[i], node << 1, i);
        int cof1 = CountBDDNodesOne(vvNodes[i], (node << 1) ^ 1, i);
        if(cof0 == cof1) {
          vvNodesRedundant[i-1].push_back(node);
        }
      }
    }
    int count = 1; // const node
    for(int i = 0; i < nInputs; i++) {
      count += vvNodes[i].size() - vvNodesRedundant[i].size();
    }
    return count;
  }

  bool Imply(int node1, int node2, int lev) {
    if(nInputs - lev > 5) {
      int nScope = 1 << (nInputs - lev - 5);
      for(int i = 0; i < nScope; i++) {
        if(t[nScope * node1 + i] & ~t[nScope * node2 + i]) {
          return 0;
        }
      }
      return 1;
    } else {
      uint value1 = BDDGetValue(node1, lev);
      uint value2 = BDDGetValue(node2, lev);
      return !(value1 & (value2 ^ ones[nInputs - lev]));
    }
  }

  int GenerateBDDBlifRec(std::vector<std::vector<int> > &vvNodes, std::vector<std::vector<int> > &vvNodeIDs, int node, int lev, int &nNodes, std::ofstream &f, std::string const &prefix) {
    if(nInputs - lev > 5) {
      int nScope = 1 << (nInputs - lev - 5);
      bool fZero = true;
      bool fOne = true;
      uint one = ones[5];
      // vertical procedure begin
      std::list<int> nodes;
      for(int node2: vvNodes[lev]) {
        nodes.push_back(node2 << 1);
        nodes.push_back((node2 << 1) ^ 1);
      }
      for(int i = 0; i < nScope; i++) {
        fZero &= t[nScope * node + i] == 0;
        fOne &= t[nScope * node + i] == one;
        for(auto it = nodes.begin(); it != nodes.end(); it++) {
          int node2 = *it >> 1;
          if(*it & 1) {
            if(~t[nScope * node + i] != t[nScope * node2 + i]) {
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
      if(fZero) {
        return 0;
      }
      if(fOne) {
        return 1;
      }
      if(!nodes.empty()) {
        auto it = std::lower_bound(vvNodes[lev].begin(), vvNodes[lev].end(), nodes.front() >> 1);
        int i = it - vvNodes[lev].begin();
        return (vvNodeIDs[lev][i] << 1) ^ (nodes.front() & 1);
      }
      // vertical procedure end
    } else {
      uint one = ones[nInputs - lev];
      uint value = BDDGetValue(node, lev);
      if(value == 0) {
        return 0;
      }
      if(value == one) {
        return 1;
      }
      for(uint i = 0; i < vvNodes[lev].size(); i++) {
        int node2 = vvNodes[lev][i];
        uint value2 = BDDGetValue(node2, lev);
        if(value2 == value) {
          return vvNodeIDs[lev][i] << 1;
        }
        if((value2 ^ one) == value) {
          return (vvNodeIDs[lev][i] << 1) ^ 1;
        }
      }
    }
    int cof0 = GenerateBDDBlifRec(vvNodes, vvNodeIDs, node << 1, lev + 1, nNodes, f, prefix);
    int cof1 = GenerateBDDBlifRec(vvNodes, vvNodeIDs, (node << 1) ^ 1, lev + 1, nNodes, f, prefix);
    if(cof0 == cof1) {
      return cof0;
    }
    int cof0id = cof0 >> 1;
    int cof1id = cof1 >> 1;
    bool cof0c = cof0 & 1;
    bool cof1c = cof1 & 1;
    bool imp01 = Imply(node << 1, (node << 1) ^ 1, lev + 1);
    bool imp10 = Imply((node << 1) ^ 1, node << 1, lev + 1);
    f << ".names " << prefix << "v" << lev << " " << prefix << "n" << cof0id << " " << prefix << "n" << cof1id << " " << prefix << "n" << nNodes << std::endl;
    f << (imp01? "-" : "0") << !cof0c << "- 1" << std::endl;
    f << (imp10? "--" : "1-") << !cof1c << " 1" << std::endl;
    vvNodes[lev].push_back(node);
    vvNodeIDs[lev].push_back(nNodes);
    return (nNodes++) << 1;
  }

  void GenerateBDDBlif(std::vector<std::string> const &inputs, std::vector<std::string> const &outputs, std::ofstream &f) {
    std::string prefix = outputs.front();
    int nNodes = 1; // const node
    std::vector<std::vector<int> > vvNodes(nInputs);
    std::vector<std::vector<int> > vvNodeIDs(nInputs);
    std::vector<int> vOutputs;
    f << ".names " << prefix << "n0" << std::endl;
    for(int i = 0; i < nInputs; i++) {
      f << ".names " << inputs[i] << " " << prefix << "v" << vLevels[i] << std::endl;
      f << "1 1" << std::endl;
    }
    for(int i = 0; i < nOutputs; i++) {
      int node = GenerateBDDBlifRec(vvNodes, vvNodeIDs, i, 0, nNodes, f, prefix);
      int id = node >> 1;
      bool c = node & 1;
      f << ".names " << prefix << "n" << id << " " << outputs[i] << std::endl;
      f << !c << " 1" << std::endl;
    }
  }

  void SwapLevel(int lev) {
    auto it0 = std::find(vLevels.begin(), vLevels.end(), lev);
    auto it1 = std::find(vLevels.begin(), vLevels.end(), lev + 1);
    std::swap(*it0, *it1);
    if(nInputs - lev > 6) {
      int nScope = 1 << (nInputs - lev - 5 - 2);
      for(int i = nScope; i < nSize * nOutputs; i += (nScope << 2)) {
        for(int j = 0; j < nScope; j++) {
          std::swap(t[i + j], t[i + nScope + j]);
        }
      }
    } else if(nInputs - lev == 6) {
      for(int i = 0; i < nSize * nOutputs; i += 2) {
        t[i+1] ^= t[i] >> 16;
        t[i] ^= t[i+1] << 16;
        t[i+1] ^= t[i] >> 16;
      }
    } else {
      for(int i = 0; i < nSize * nOutputs; i++) {
        int d = (nInputs - lev - 2);
        int shamt = 1 << d;
        t[i] ^= (t[i] >> shamt) & sieve[d];
        t[i] ^= (t[i] & sieve[d]) << shamt;
        t[i] ^= (t[i] >> shamt) & sieve[d];
      }
    }
  }
};

const uint TT::ones[] = {0x00000001,
                         0x00000003,
                         0x0000000f,
                         0x000000ff,
                         0x0000ffff,
                         0xffffffff};

const uint TT::sieve[] = {0x22222222,
                          0x0c0c0c0c,
                          0x00f000f0,
                          0x0000ff00};

void TTTest(std::vector<std::vector<int> > const &onsets, std::vector<char *> const &pBPats, int nBPats, int rarity, std::vector<std::string> const &inputs, std::vector<std::string> const &outputs, std::ofstream &f) {
  int nInputs = inputs.size();
  TT tt(onsets, nInputs);
  //tt.GeneratePla("test2.pla");
  std::cout << tt.CountBDDNodes() << std::endl;
  for(int i = 0; i < 100; i++) {
    tt.SwapLevel(rand() % (nInputs-1));
  }
  std::cout << tt.CountBDDNodes() << std::endl;
  tt.GenerateBDDBlif(inputs, outputs, f);
}
