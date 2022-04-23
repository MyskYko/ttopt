#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <random>
#include <cassert>
#include <map>
#include <cuddObj.hh>

class Bdd {
public:
  typedef uint64_t word;
  const int ww = 64; // word width
  const int lww = 6; // log word width

  int nInputs;
  int nSize;
  int nTotalSize;
  int nOutputs;

  std::vector<word> t;

  CUDD::Cudd *man;
  std::vector<CUDD::BDD> vNodes;

  static const word ones[];

  Bdd(std::vector<std::vector<int> > const &onsets, int nInputs): nInputs(nInputs) {
    man = new CUDD::Cudd(nInputs, 0);
    nOutputs = onsets.size();
    if(nInputs >= lww) {
      nSize = 1 << (nInputs - lww);
      nTotalSize = nSize * nOutputs;
      t.resize(nTotalSize);
      for(int i = 0; i < nOutputs; i++) {
        for(int pat: onsets[i]) {
          int index = pat / ww;
          int pos = pat % ww;
          t[nSize * i + index] |= 1ull << pos;
        }
      }
    } else {
      nSize = 0;
      nTotalSize = ((1 << nInputs) * nOutputs + ww - 1) / ww;
      t.resize(nTotalSize);
      for(int i = 0; i < nOutputs; i++) {
        int padding = i * (1 << nInputs);
        for(int pat: onsets[i]) {
          int pos = (padding + pat) % ww;
          t[padding / ww] |= 1ull << pos;
        }
      }
    }
  }

  ~Bdd() {
    vNodes.clear();
    delete man;
  }

  word GetValue(int index_lev, int lev) {
    assert(index_lev >= 0);
    assert(nInputs - lev <= lww);
    int logwidth = nInputs - lev;
    int index = index_lev >> (lww - logwidth);
    int pos = (index_lev % (1 << (lww - logwidth))) << logwidth;
    return (t[index] >> pos) & ones[logwidth];
  }

  CUDD::BDD BuildRec(int index, int var) {
    int logwidth = nInputs - var;
    if(logwidth <= lww) {
      word value = GetValue(index, var);
      if(value == 0) {
        return !man->bddOne();
      }
      if(value == ones[logwidth]) {
        return man->bddOne();
      }
    }
    auto cof0 = BuildRec(index << 1, var + 1);
    auto cof1 = BuildRec((index << 1) + 1, var + 1);
    auto v = man->ReadVars(var);
    return v.Ite(cof1, cof0);
  }

  void Build() {
    for(int i = 0; i < nOutputs; i++) {
      auto p = BuildRec(i, 0);
      vNodes.push_back(p);
    }
  }

  int GenerateBlifRec(CUDD::BDD const &x, int &nNodes, std::map<uint64_t, int> &m, std::string const &prefix, std::ofstream &f) {
    if(Cudd_IsComplement(x.getNode())) {
      return GenerateBlifRec(!x, nNodes, m, prefix, f) ^ 1;
    }
    uint64_t id = (uint64_t)x.getNode();
    if(m.count(id)) {
      return m[id] << 1;
    }
    auto cof0 = GenerateBlifRec(CUDD::BDD(*man, Cudd_E(x.getNode())), nNodes, m, prefix, f);
    auto cof1 = GenerateBlifRec(CUDD::BDD(*man, Cudd_T(x.getNode())), nNodes, m, prefix, f);
    int var = x.NodeReadIndex();
    f << ".names " << prefix << "v" << var << " " << prefix << "n" << (cof0 >> 1) << " " << prefix << "n" << (cof1 >> 1) << " " << prefix << "n" << nNodes << std::endl;
    f << "0" << !(cof0 & 1) << "- 1" << std::endl;
    f << "1-" << !(cof1 & 1) << " 1" << std::endl;
    m[id] = nNodes;
    return (nNodes++) << 1;
  }

  void GenerateBlif(std::vector<std::string> const &inputs, std::vector<std::string> const &outputs, std::ofstream &f) {
    std::string prefix = outputs.front();
    std::map<uint64_t, int> m;
    f << ".names " << prefix << "n0" << std::endl << "1" << std::endl;
    uint64_t id = (uint64_t)man->bddOne().getNode();
    m[id] = 0;
    int nNodes = 1; // const node
    for(int i = 0; i < nInputs; i++) {
      f << ".names " << inputs[i] << " " << prefix << "v" << i << std::endl;
      f << "1 1" << std::endl;
    }
    for(int i = 0; i < nOutputs; i++) {
      int node = GenerateBlifRec(vNodes[i], nNodes, m, prefix, f);
      f << ".names " << prefix << "n" << (node >> 1) << " " << outputs[i] << std::endl;
      f << !(node & 1) << " 1" << std::endl;
    }
  }
};

const Bdd::word Bdd::ones[] = {0x0000000000000001ull,
                               0x0000000000000003ull,
                               0x000000000000000full,
                               0x00000000000000ffull,
                               0x000000000000ffffull,
                               0x00000000ffffffffull,
                               0xffffffffffffffffull};

void BddTest(std::vector<std::vector<int> > const &onsets, std::vector<std::string> const &inputs, std::vector<std::string> const &outputs, std::ofstream &f) {
  int nInputs = inputs.size();
  auto m = Bdd(onsets, nInputs);
  m.Build();

  m.man->SetMaxGrowth(100);
  m.man->SetSiftMaxSwap(1000000000);
  m.man->SetSiftMaxVar(1000000000);

  m.man->ReduceHeap(CUDD_REORDER_SIFT);
  int best = m.man->SharingSize(m.vNodes);
  int perm[nInputs];
  for(int j = 0; j < nInputs; j++) {
    perm[m.man->ReadPerm(j)] = j;
  }
  std::mt19937 rng;
  for(int i = 0; i < 20; i++) {
    int vLevelsNew[nInputs];
    std::iota(vLevelsNew, vLevelsNew + nInputs, 0);
    std::shuffle(vLevelsNew, vLevelsNew + nInputs, rng);
    m.man->ShuffleHeap(vLevelsNew);
    m.man->ReduceHeap(CUDD_REORDER_SIFT);
    if(best > m.man->SharingSize(m.vNodes)){
      best = m.man->SharingSize(m.vNodes);
      for(int j = 0; j < nInputs; j++) {
        perm[m.man->ReadPerm(j)] = j;
      }
    }
  }
  // m.man->ShuffleHeap(perm);

  // m.GenerateBlif(inputs, outputs, f);
}
