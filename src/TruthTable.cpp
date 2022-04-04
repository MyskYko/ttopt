#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <list>
#include <algorithm>
#include <numeric>
#include <random>
#include <cassert>
#include <map>
#include <bitset>
#include <unordered_map>

extern std::string BinaryToString(int bin, int size);

template <class T>
inline void hash_combine(std::size_t & seed, const T & v)
{
  std::hash<T> hasher;
  seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

namespace std
{
  template<typename S, typename T> struct hash<pair<S, T>>
  {
    inline size_t operator()(const pair<S, T> & v) const
    {
      size_t seed = 0;
      ::hash_combine(seed, v.first);
      ::hash_combine(seed, v.second);
      return seed;
    }
  };
}

class TruthTable {
public:
  typedef uint64_t word;
  const int ww = 64; // word width
  const int lww = 6; // log word width
  typedef std::bitset<64> bsw;

  int nInputs;
  int nSize;
  int nTotalSize;
  int nOutputs;
  std::vector<word> t;

  std::vector<std::vector<int> > vvIndices;
  std::vector<std::vector<int> > vvRedundantIndices;
  std::vector<int> vLevels;

  std::vector<std::vector<word> > savedt;
  std::vector<std::vector<std::vector<int> > > vvIndicesSaved;
  std::vector<std::vector<std::vector<int> > > vvRedundantIndicesSaved;
  std::vector<std::vector<int> > vLevelsSaved;

  std::mt19937 rng;
  static const word ones[];
  static const word swapmask[];

  TruthTable(std::vector<std::vector<int> > const &onsets, int nInputs): nInputs(nInputs) {
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
    vLevels.resize(nInputs);
    std::iota(vLevels.begin(), vLevels.end(), 0);
  }

  void GeneratePla(std::string filename) {
    std::ofstream f(filename);
    f << ".i " << nInputs << std::endl;
    f << ".o " << nOutputs << std::endl;
    if(nSize) {
      for(int index = 0; index < nSize; index++) {
        for(int pos = 0; pos < ww; pos++) {
          int pat = (index << lww) + pos;
          f << BinaryToString(pat, nInputs) << " ";
          for(int i = 0; i < nOutputs; i++) {
            f << ((t[nSize * i + index] >> pos) & 1);
          }
          f << std::endl;
        }
      }
    } else {
      for(int pos = 0; pos < (1 << nInputs); pos++) {
        f << BinaryToString(pos, nInputs) << " ";
        for(int i = 0; i < nOutputs; i++) {
          int padding = i * (1 << nInputs);
          f << ((t[padding / ww] >> (pos + padding % ww)) & 1);
        }
        f << std::endl;
      }
    }
  }

  virtual void Save(uint i) {
    if(savedt.size() < i + 1) {
      savedt.resize(i + 1);
      vLevelsSaved.resize(i + 1);
    }
    savedt[i] = t;
    vLevelsSaved[i] = vLevels;
  }

  virtual void Load(uint i) {
    assert(i < savedt.size());
    t = savedt[i];
    vLevels = vLevelsSaved[i];
  }

  virtual void SaveIndices(uint i) {
    if(vvIndicesSaved.size() < i + 1) {
      vvIndicesSaved.resize(i + 1);
      vvRedundantIndicesSaved.resize(i + 1);
    }
    vvIndicesSaved[i] = vvIndices;
    vvRedundantIndicesSaved[i] = vvRedundantIndices;
  }

  virtual void LoadIndices(uint i) {
    vvIndices = vvIndicesSaved[i];
    vvRedundantIndices = vvRedundantIndicesSaved[i];
  }

  word GetValue(int index_lev, int lev) {
    assert(index_lev >= 0);
    assert(nInputs - lev <= lww);
    int logwidth = nInputs - lev;
    int index = index_lev >> (lww - logwidth);
    int pos = (index_lev % (1 << (lww - logwidth))) << logwidth;
    return (t[index] >> pos) & ones[logwidth];
  }

  int IsEq(int index1, int index2, int lev, bool fCompl = false) {
    assert(index1 >= 0);
    assert(index2 >= 0);
    int logwidth = nInputs - lev;
    bool fEq = true;
    if(logwidth > lww) {
      int nScopeSize = 1 << (logwidth - lww);
      for(int i = 0; i < nScopeSize && (fEq || fCompl); i++) {
        fEq &= t[nScopeSize * index1 + i] == t[nScopeSize * index2 + i];
        fCompl &= t[nScopeSize * index1 + i] == ~t[nScopeSize * index2 + i];
      }
    } else {
      word value = GetValue(index1, lev) ^ GetValue(index2, lev);
      fEq &= !value;
      fCompl &= !(value ^ ones[logwidth]);
    }
    return 2 * fCompl + fEq;
  }

  bool Imply(int index1, int index2, int lev) {
    assert(index1 >= 0);
    assert(index2 >= 0);
    int logwidth = nInputs - lev;
    if(logwidth > lww) {
      int nScopeSize = 1 << (logwidth - lww);
      for(int i = 0; i < nScopeSize; i++) {
        if(t[nScopeSize * index1 + i] & ~t[nScopeSize * index2 + i]) {
          return false;
        }
      }
      return true;
    }
    return !(GetValue(index1, lev) & (GetValue(index2, lev) ^ ones[logwidth]));
  }

  int BDDNodeCountLevel(int lev) {
    return vvIndices[lev].size() - vvRedundantIndices[lev].size();
  }

  int BDDNodeCount() {
    int count = 1; // const node
    for(int i = 0; i < nInputs; i++) {
      count += BDDNodeCountLevel(i);
    }
    return count;
  }

  int BDDFind(int index, int lev) {
    int logwidth = nInputs - lev;
    if(logwidth > lww) {
      int nScopeSize = 1 << (logwidth - lww);
      bool fZero = true;
      bool fOne = true;
      for(int i = 0; i < nScopeSize && (fZero || fOne); i++) {
        word value = t[nScopeSize * index + i];
        fZero &= !value;
        fOne &= !(~value);
      }
      if(fZero || fOne) {
        return -2 ^ fOne;
      }
      for(uint j = 0; j < vvIndices[lev].size(); j++) {
        int index2 = vvIndices[lev][j];
        bool fEq = true;
        bool fCompl = true;
        for(int i = 0; i < nScopeSize && (fEq || fCompl); i++) {
          fEq &= t[nScopeSize * index + i] == t[nScopeSize * index2 + i];
          fCompl &= t[nScopeSize * index + i] == ~t[nScopeSize * index2 + i];
        }
        if(fEq || fCompl) {
          return (j << 1) ^ fCompl;
        }
      }
    } else {
      word value = GetValue(index, lev);
      if(!value) {
        return -2;
      }
      if(!(value ^ ones[logwidth])) {
        return -1;
      }
      for(uint j = 0; j < vvIndices[lev].size(); j++) {
        int index2 = vvIndices[lev][j];
        word value2 = value ^ GetValue(index2, lev);
        if(!(value2)) {
          return j << 1;
        }
        if(!(value2 ^ ones[logwidth])) {
          return (j << 1) ^ 1;
        }
      }
    }
    return -3;
  }

  virtual int BDDBuildOne(int index, int lev) {
    int r = BDDFind(index, lev);
    if(r >= -2) {
      return r;
    }
    vvIndices[lev].push_back(index);
    return (vvIndices[lev].size() - 1) << 1;
  }

  virtual void BDDBuildStartup() {
    vvIndices.clear();
    vvIndices.resize(nInputs);
    vvRedundantIndices.clear();
    vvRedundantIndices.resize(nInputs);
    for(int i = 0; i < nOutputs; i++) {
      BDDBuildOne(i, 0);
    }
  }

  virtual void BDDBuildLevel(int lev) {
    for(int index: vvIndices[lev-1]) {
      int cof0 = BDDBuildOne(index << 1, lev);
      int cof1 = BDDBuildOne((index << 1) ^ 1, lev);
      if(cof0 == cof1) {
        vvRedundantIndices[lev-1].push_back(index);
      }
    }
  }

  virtual int BDDBuild() {
    BDDBuildStartup();
    for(int i = 1; i < nInputs; i++) {
      BDDBuildLevel(i);
    }
    return BDDNodeCount();
  }

  virtual int BDDRebuild(int lev) {
    vvIndices[lev].clear();
    vvIndices[lev+1].clear();
    for(int i = lev; i < lev + 2; i++) {
      if(!i) {
        for(int j = 0; j < nOutputs; j++) {
          BDDBuildOne(j, 0);
        }
      } else {
        vvRedundantIndices[i-1].clear();
        BDDBuildLevel(i);
      }
    }
    if(lev < nInputs - 2) {
      vvRedundantIndices[lev+1].clear();
      for(int index: vvIndices[lev+1]) {
        if(IsEq(index << 1, (index << 1) ^ 1, lev + 2)) {
          vvRedundantIndices[lev+1].push_back(index);
        }
      }
    }
    return BDDNodeCount();
  }

  virtual void Swap(int lev) {
    assert(lev < nInputs - 1);
    auto it0 = std::find(vLevels.begin(), vLevels.end(), lev);
    auto it1 = std::find(vLevels.begin(), vLevels.end(), lev + 1);
    std::swap(*it0, *it1);
    if(nInputs - lev - 1 > lww) {
      int nScopeSize = 1 << (nInputs - lev - 2 - lww);
      for(int i = nScopeSize; i < nTotalSize; i += (nScopeSize << 2)) {
        for(int j = 0; j < nScopeSize; j++) {
          std::swap(t[i + j], t[i + nScopeSize + j]);
        }
      }
    } else if(nInputs - lev - 1 == lww) {
      for(int i = 0; i < nTotalSize; i += 2) {
        t[i+1] ^= t[i] >> (ww / 2);
        t[i] ^= t[i+1] << (ww / 2);
        t[i+1] ^= t[i] >> (ww / 2);
      }
    } else {
      for(int i = 0; i < nTotalSize; i++) {
        int d = nInputs - lev - 2;
        int shamt = 1 << d;
        t[i] ^= (t[i] >> shamt) & swapmask[d];
        t[i] ^= (t[i] & swapmask[d]) << shamt;
        t[i] ^= (t[i] >> shamt) & swapmask[d];
      }
    }
  }

  void SwapIndex(int &index, int d) {
    if((index >> d) % 4 == 1) {
      index += 1 << d;
    } else if((index >> d) % 4 == 2) {
      index -= 1 << d;
    }
  }

  virtual int BDDSwap(int lev) {
    Swap(lev);
    for(int i = lev + 2; i < nInputs; i++) {
      for(uint j = 0; j < vvIndices[i].size(); j++) {
        SwapIndex(vvIndices[i][j], i - (lev + 2));
      }
    }
    // swapping vvRedundantIndices is unnecessary for node counting
    return BDDRebuild(lev);
  }

  int SiftReo() {
    int best = BDDBuild();
    Save(0);
    SaveIndices(0);
    std::vector<int> vars(nInputs);
    std::iota(vars.begin(), vars.end(), 0);
    std::sort(vars.begin(), vars.end(), [&](int i1, int i2) {return BDDNodeCountLevel(vLevels[i1]) > BDDNodeCountLevel(vLevels[i2]);});
    bool turn = true;
    for(int var: vars) {
      bool updated = false;
      int lev = vLevels[var];
      for(int i = lev; i < nInputs - 1; i++) {
        int count = BDDSwap(i);
        if(best > count) {
          best = count;
          updated = true;
          Save(turn);
          SaveIndices(turn);
        }
      }
      if(lev) {
        Load(!turn);
        LoadIndices(!turn);
        for(int i = lev - 1; i >= 0; i--) {
          int count = BDDSwap(i);
          if(best > count) {
            best = count;
            updated = true;
            Save(turn);
            SaveIndices(turn);
          }
        }
      }
      turn ^= updated;
      Load(!turn);
      LoadIndices(!turn);
    }
    return best;
  }

  void Reo(std::vector<int> vLevelsNew) {
    for(int i = 0; i < nInputs; i++) {
      int var = std::find(vLevelsNew.begin(), vLevelsNew.end(), i) - vLevelsNew.begin();
      int lev = vLevels[var];
      if(lev < i) {
        for(int j = lev; j < i; j++) {
          Swap(j);
        }
      } else if(lev > i) {
        for(int j = lev - 1; j >= i; j--) {
          Swap(j);
        }
      }
    }
    assert(vLevels == vLevelsNew);
  }

  int RandomSiftReo(int nRound) {
    int best = SiftReo();
    Save(2);
    for(int i = 0; i < nRound; i++) {
      std::vector<int> vLevelsNew(nInputs);
      std::iota(vLevelsNew.begin(), vLevelsNew.end(), 0);
      std::shuffle(vLevelsNew.begin(), vLevelsNew.end(), rng);
      Reo(vLevelsNew);
      int r = SiftReo();
      if(best > r) {
        best = r;
        Save(2);
      }
    }
    Load(2);
    return best;
  }

  int BDDGenerateBlifRec(std::vector<std::vector<int> > &vvNodes, int &nNodes, int index, int lev, std::ofstream &f, std::string const &prefix) {
    int r = BDDFind(index, lev);
    if(r >= 0) {
      return (vvNodes[lev][r >> 1] << 1) ^ (r & 1);
    }
    if(r >= -2) {
      return r + 2;
    }
    int cof0 = BDDGenerateBlifRec(vvNodes, nNodes, index << 1, lev + 1, f, prefix);
    int cof1 = BDDGenerateBlifRec(vvNodes, nNodes, (index << 1) ^ 1, lev + 1, f, prefix);
    if(cof0 == cof1) {
      return cof0;
    }
    f << ".names " << prefix << "v" << lev << " " << prefix << "n" << (cof0 >> 1) << " " << prefix << "n" << (cof1 >> 1) << " " << prefix << "n" << nNodes << std::endl;
    f << (Imply(index << 1, (index << 1) ^ 1, lev + 1)? "-" : "0") << !(cof0 & 1) << "- 1" << std::endl;
    f << (Imply((index << 1) ^ 1, index << 1, lev + 1)? "--" : "1-") << !(cof1 & 1) << " 1" << std::endl;
    vvIndices[lev].push_back(index);
    vvNodes[lev].push_back(nNodes);
    return (nNodes++) << 1;
  }

  virtual void BDDGenerateBlif(std::vector<std::string> const &inputs, std::vector<std::string> const &outputs, std::ofstream &f) {
    std::string prefix = outputs.front();
    int nNodes = 1; // const node
    vvIndices.clear();
    vvIndices.resize(nInputs);
    std::vector<std::vector<int> > vvNodes(nInputs);
    std::vector<int> vOutputs;
    f << ".names " << prefix << "n0" << std::endl;
    for(int i = 0; i < nInputs; i++) {
      f << ".names " << inputs[i] << " " << prefix << "v" << vLevels[i] << std::endl;
      f << "1 1" << std::endl;
    }
    for(int i = 0; i < nOutputs; i++) {
      int node = BDDGenerateBlifRec(vvNodes, nNodes, i, 0, f, prefix);
      f << ".names " << prefix << "n" << (node >> 1) << " " << outputs[i] << std::endl;
      f << !(node & 1) << " 1" << std::endl;
    }
  }
};

const TruthTable::word TruthTable::ones[] = {0x0000000000000001ull,
                                             0x0000000000000003ull,
                                             0x000000000000000full,
                                             0x00000000000000ffull,
                                             0x000000000000ffffull,
                                             0x00000000ffffffffull,
                                             0xffffffffffffffffull};

const TruthTable::word TruthTable::swapmask[] = {0x2222222222222222ull,
                                                 0x0c0c0c0c0c0c0c0cull,
                                                 0x00f000f000f000f0ull,
                                                 0x0000ff000000ff00ull,
                                                 0x00000000ffff0000ull};

class TruthTableReo : public TruthTable {
public:
  bool fBuilt = false;
  std::vector<std::vector<int> > vvChildren;
  std::vector<std::vector<std::vector<int> > > vvChildrenSaved;

  TruthTableReo(std::vector<std::vector<int> > const &onsets, int nInputs): TruthTable(onsets, nInputs) {}

  void Save(uint i) override {
    if(vLevelsSaved.size() < i + 1) {
      vLevelsSaved.resize(i + 1);
    }
    vLevelsSaved[i] = vLevels;
  }

  void Load(uint i) override {
    assert(i < vLevelsSaved.size());
    vLevels = vLevelsSaved[i];
  }

  void SaveIndices(uint i) override {
    TruthTable::SaveIndices(i);
    if(vvChildrenSaved.size() < i + 1) {
      vvChildrenSaved.resize(i + 1);
    }
    vvChildrenSaved[i] = vvChildren;
  }

  void LoadIndices(uint i) override {
    TruthTable::LoadIndices(i);
    vvChildren = vvChildrenSaved[i];
  }

  void BDDBuildStartup() override {
    vvChildren.clear();
    vvChildren.resize(nInputs);
    TruthTable::BDDBuildStartup();
  }

  void BDDBuildLevel(int lev) override {
    for(int index: vvIndices[lev-1]) {
      int cof0 = BDDBuildOne(index << 1, lev);
      int cof1 = BDDBuildOne((index << 1) ^ 1, lev);
      vvChildren[lev-1].push_back(cof0);
      vvChildren[lev-1].push_back(cof1);
      if(cof0 == cof1) {
        vvRedundantIndices[lev-1].push_back(index);
      }
    }
  }

  int BDDBuild() override {
    if(fBuilt) {
      return BDDNodeCount();
    }
    fBuilt = true;
    BDDBuildStartup();
    for(int i = 1; i < nInputs + 1; i++) {
      BDDBuildLevel(i);
    }
    return BDDNodeCount();
  }

  int BDDRebuildOne(int index, int cof0, int cof1, int lev, std::unordered_map<std::pair<int, int>, int> &unique, std::vector<int> &vChildrenLow) {
    if(cof0 < 0 && cof0 == cof1) {
      return cof0;
    }
    bool fCompl = cof0 & 1;
    if(fCompl) {
      cof0 ^= 1;
      cof1 ^= 1;
    }
    if(unique.count({cof0, cof1})) {
      return (unique[{cof0, cof1}] << 1) ^ fCompl;
    }
    vvIndices[lev].push_back(index);
    unique[{cof0, cof1}] = vvIndices[lev].size() - 1;
    vChildrenLow.push_back(cof0);
    vChildrenLow.push_back(cof1);
    if(cof0 == cof1) {
      vvRedundantIndices[lev].push_back(index);
    }
    return ((vvIndices[lev].size() - 1) << 1) ^ fCompl;
  }

  int BDDRebuild(int lev) override {
    vvRedundantIndices[lev].clear();
    vvRedundantIndices[lev+1].clear();
    std::vector<int> vChildrenHigh;
    std::vector<int> vChildrenLow;
    std::unordered_map<std::pair<int, int>, int> unique;
    unique.reserve(2 * vvIndices[lev+1].size());
    vvIndices[lev+1].clear();
    for(uint i = 0; i < vvIndices[lev].size(); i++) {
      int index = vvIndices[lev][i];
      int cof0index = vvChildren[lev][i+i] >> 1;
      int cof1index = vvChildren[lev][i+i+1] >> 1;
      bool cof0c = vvChildren[lev][i+i] & 1;
      bool cof1c = vvChildren[lev][i+i+1] & 1;
      int cof00, cof01, cof10, cof11;
      if(cof0index < 0) {
        cof00 = -2 ^ cof0c;
        cof01 = -2 ^ cof0c;
      } else {
        cof00 = vvChildren[lev+1][cof0index+cof0index] ^ cof0c;
        cof01 = vvChildren[lev+1][cof0index+cof0index+1] ^ cof0c;
      }
      if(cof1index < 0) {
        cof10 = -2 ^ cof1c;
        cof11 = -2 ^ cof1c;
      } else {
        cof10 = vvChildren[lev+1][cof1index+cof1index] ^ cof1c;
        cof11 = vvChildren[lev+1][cof1index+cof1index+1] ^ cof1c;
      }
      int newcof0 = BDDRebuildOne(index << 1, cof00, cof10, lev + 1, unique, vChildrenLow);
      int newcof1 = BDDRebuildOne((index << 1) ^ 1, cof01, cof11, lev + 1, unique, vChildrenLow);
      vChildrenHigh.push_back(newcof0);
      vChildrenHigh.push_back(newcof1);
      if(newcof0 == newcof1) {
        vvRedundantIndices[lev].push_back(index);
      }
    }
    vvChildren[lev] = vChildrenHigh;
    vvChildren[lev+1] = vChildrenLow;
    return BDDNodeCount();
  }

  void Swap(int lev) override {
    assert(lev < nInputs - 1);
    auto it0 = std::find(vLevels.begin(), vLevels.end(), lev);
    auto it1 = std::find(vLevels.begin(), vLevels.end(), lev + 1);
    std::swap(*it0, *it1);
    BDDRebuild(lev);
  }

  int BDDSwap(int lev) override {
    Swap(lev);
    return BDDNodeCount();
  }

  void BDDGenerateBlif(std::vector<std::string> const &inputs, std::vector<std::string> const &outputs, std::ofstream &f) override {
    abort();
  }
};

class TruthTableRewrite : public TruthTable {
public:
  TruthTableRewrite(std::vector<std::vector<int> > const &onsets, int nInputs): TruthTable(onsets, nInputs) {}

  void SetValue(int index_lev, int lev, word value) {
    assert(index_lev >= 0);
    assert(nInputs - lev <= lww);
    int logwidth = nInputs - lev;
    int index = index_lev >> (lww - logwidth);
    int pos = (index_lev % (1 << (lww - logwidth))) << logwidth;
    t[index] &= ~(ones[logwidth] << pos);
    t[index] ^= value << pos;
  }

  void CopyFunc(int index1, int index2, int lev, bool fCompl) {
    assert(index1 >= 0);
    int logwidth = nInputs - lev;
    if(logwidth > lww) {
      int nScopeSize = 1 << (logwidth - lww);
      if(!fCompl) {
        if(index2 < 0) {
          for(int i = 0; i < nScopeSize; i++) {
            t[nScopeSize * index1 + i] = 0;
          }
        } else {
          for(int i = 0; i < nScopeSize; i++) {
            t[nScopeSize * index1 + i] = t[nScopeSize * index2 + i];
          }
        }
      } else {
        if(index2 < 0) {
          for(int i = 0; i < nScopeSize; i++) {
            t[nScopeSize * index1 + i] = ones[lww];
          }
        } else {
          for(int i = 0; i < nScopeSize; i++) {
            t[nScopeSize * index1 + i] = ~t[nScopeSize * index2 + i];
          }
        }
      }
    } else {
      word value = 0;
      if(index2 >= 0) {
        value = GetValue(index2, lev);
      }
      if(fCompl) {
        value ^= ones[logwidth];
      }
      SetValue(index1, lev, value);
    }
  }

  void ShiftToMajority(int index, int lev) {
    assert(index >= 0);
    int logwidth = nInputs - lev;
    int count = 0;
    if(logwidth > lww) {
      int nScopeSize = 1 << (logwidth - lww);
      for(int i = 0; i < nScopeSize; i++) {
        count += bsw(t[nScopeSize * index + i]).count();
      }
    } else {
      count = bsw(GetValue(index, lev)).count();
    }
    bool majority = count > (1 << (logwidth - 1));
    CopyFunc(index, -1, lev, majority);
  }
};

class TruthTableCare : public TruthTableRewrite {
public:
  std::vector<word> originalt;
  std::vector<word> caret;
  std::vector<word> care;

  std::vector<std::vector<std::pair<int, int> > > vvMergedIndices;

  std::vector<std::vector<word> > savedcare;
  std::vector<std::vector<std::vector<std::pair<int, int> > > > vvMergedIndicesSaved;

  TruthTableCare(std::vector<std::vector<int> > const &onsets, int nInputs, std::vector<char *> const &pBPats, int nBPats, int rarity): TruthTableRewrite(onsets, nInputs) {
    if(nSize) {
      care.resize(nSize);
    } else {
      care.resize(1);
    }
    std::vector<int> count(1 << nInputs);
    for(int i = 0; i < nBPats; i++) {
      for(int j = 0; j < 8; j++) {
        int pat = 0;
        for(auto pBPat: pBPats) {
          pat <<= 1;
          pat |= (pBPat[i] >> j) & 1;
        }
        count[pat]++;
        if(count[pat] == rarity) {
          int index = pat / ww;
          int pos = pat % ww;
          care[index] |= 1ull << pos;
        }
      }
    }
  }

  void Save(uint i) override {
    TruthTable::Save(i);
    if(savedcare.size() < i + 1) {
      savedcare.resize(i + 1);
    }
    savedcare[i] = care;
  }

  void Load(uint i) override {
    TruthTable::Load(i);
    care = savedcare[i];
  }

  void SaveIndices(uint i) override {
    TruthTable::SaveIndices(i);
    if(vvMergedIndicesSaved.size() < i + 1) {
      vvMergedIndicesSaved.resize(i + 1);
    }
    vvMergedIndicesSaved[i] = vvMergedIndices;
  }

  void LoadIndices(uint i) override {
    TruthTable::LoadIndices(i);
    vvMergedIndices = vvMergedIndicesSaved[i];
  }

  void Swap(int lev) override {
    TruthTable::Swap(lev);
    if(nInputs - lev - 1 > lww) {
      int nScopeSize = 1 << (nInputs - lev - 2 - lww);
      for(int i = nScopeSize; i < nSize; i += (nScopeSize << 2)) {
        for(int j = 0; j < nScopeSize; j++) {
          std::swap(care[i + j], care[i + nScopeSize + j]);
        }
      }
    } else if(nInputs - lev - 1 == lww) {
      for(int i = 0; i < nSize; i += 2) {
        care[i+1] ^= care[i] >> (ww / 2);
        care[i] ^= care[i+1] << (ww / 2);
        care[i+1] ^= care[i] >> (ww / 2);
      }
    } else {
      for(int i = 0; i < nSize || (i == 0 && !nSize); i++) {
        int d = nInputs - lev - 2;
        int shamt = 1 << d;
        care[i] ^= (care[i] >> shamt) & swapmask[d];
        care[i] ^= (care[i] & swapmask[d]) << shamt;
        care[i] ^= (care[i] >> shamt) & swapmask[d];
      }
    }
  }

  void GeneratePlaCare(std::string filename) {
    std::ofstream f(filename);
    f << ".i " << nInputs << std::endl;
    f << ".o 1" << std::endl;
    if(nSize) {
      for(int index = 0; index < nSize; index++) {
        for(int pos = 0; pos < ww; pos++) {
          int pat = (index << lww) + pos;
          f << BinaryToString(pat, nInputs) << " ";
          f << ((care[index] >> pos) & 1);
          f << std::endl;
        }
      }
    } else {
      for(int pos = 0; pos < (1 << nInputs); pos++) {
        f << BinaryToString(pos, nInputs) << " ";
        f << ((care[0] >> pos) & 1);
        f << std::endl;
      }
    }
  }

  void GeneratePlaMasked(std::string filename) {
    std::ofstream f(filename);
    f << ".i " << nInputs << std::endl;
    f << ".o " << nOutputs << std::endl;
    if(nSize) {
      for(int index = 0; index < nSize; index++) {
        for(int pos = 0; pos < ww; pos++) {
          int pat = (index << lww) + pos;
          f << BinaryToString(pat, nInputs) << " ";
          for(int i = 0; i < nOutputs; i++) {
            if((care[index] >> pos) & 1) {
              f << ((t[nSize * i + index] >> pos) & 1);
            } else {
              f << ((originalt[nSize * i + index] >> pos) & 1);
            }
          }
          f << std::endl;
        }
      }
    } else {
      for(int pos = 0; pos < (1 << nInputs); pos++) {
        f << BinaryToString(pos, nInputs) << " ";
        for(int i = 0; i < nOutputs; i++) {
          int padding = i * (1 << nInputs);
          if((care[0] >> pos) & 1) {
            f << ((t[padding / ww] >> (pos + padding % ww)) & 1);
          } else {
            f << ((originalt[padding / ww] >> (pos + padding % ww)) & 1);
          }
        }
        f << std::endl;
      }
    }
  }

  void RestoreCare() {
    caret.clear();
    if(nSize) {
      for(int i = 0; i < nOutputs; i++) {
        caret.insert(caret.end(), care.begin(), care.end());
      }
    } else {
      caret.resize(nTotalSize);
      for(int i = 0; i < nOutputs; i++) {
        int padding = i * (1 << nInputs);
        caret[padding / ww] |= care[0] << (padding % ww);
      }
    }
  }

  word GetCare(int index_lev, int lev) {
    assert(index_lev >= 0);
    assert(nInputs - lev <= lww);
    int logwidth = nInputs - lev;
    int index = index_lev >> (lww - logwidth);
    int pos = (index_lev % (1 << (lww - logwidth))) << logwidth;
    return (caret[index] >> pos) & ones[logwidth];
  }

  void CopyFuncMasked(int index1, int index2, int lev, bool fCompl) {
    assert(index1 >= 0);
    assert(index2 >= 0);
    int logwidth = nInputs - lev;
    if(logwidth > lww) {
      int nScopeSize = 1 << (logwidth - lww);
      for(int i = 0; i < nScopeSize; i++) {
        word value = t[nScopeSize * index2 + i];
        if(fCompl) {
          value = ~value;
        }
        word cvalue = caret[nScopeSize * index2 + i];
        t[nScopeSize * index1 + i] &= ~cvalue;
        t[nScopeSize * index1 + i] |= cvalue & value;
      }
    } else {
      word one = ones[logwidth];
      word value1 = GetValue(index1, lev);
      word value2 = GetValue(index2, lev);
      if(fCompl) {
        value2 ^= one;
      }
      word cvalue = GetCare(index2, lev);
      value1 &= cvalue ^ one;
      value1 |= cvalue & value2;
      SetValue(index1, lev, value1);
    }
  }

  bool IsDC(int index, int lev) {
    if(nInputs - lev > lww) {
      int nScopeSize = 1 << (nInputs - lev - lww);
      for(int i = 0; i < nScopeSize; i++) {
        if(caret[nScopeSize * index + i]) {
          return false;
        }
      }
    } else if(GetCare(index, lev)) {
      return false;
    }
    return true;
  }

  int Include(int index1, int index2, int lev, bool fCompl) {
    assert(index1 >= 0);
    assert(index2 >= 0);
    int logwidth = nInputs - lev;
    bool fEq = true;
    if(logwidth > lww) {
      int nScopeSize = 1 << (logwidth - lww);
      for(int i = 0; i < nScopeSize && (fEq || fCompl); i++) {
        word cvalue = caret[nScopeSize * index2 + i];
        if(~caret[nScopeSize * index1 + i] & cvalue) {
          return 0;
        }
        word value = t[nScopeSize * index1 + i] ^ t[nScopeSize * index2 + i];
        fEq &= !(value & cvalue);
        fCompl &= !(~value & cvalue);
      }
    } else {
      word cvalue = GetCare(index2, lev);
      if((GetCare(index1, lev) ^ ones[logwidth]) & cvalue) {
        return 0;
      }
      word value = GetValue(index1, lev) ^ GetValue(index2, lev);
      fEq &= !(value & cvalue);
      fCompl &= !((value ^ ones[logwidth]) & cvalue);
    }
    return 2 * fCompl + fEq;
  }

  int Intersect(int index1, int index2, int lev, bool fCompl, bool fEq = true) {
    assert(index1 >= 0);
    assert(index2 >= 0);
    int logwidth = nInputs - lev;
    if(logwidth > lww) {
      int nScopeSize = 1 << (logwidth - lww);
      for(int i = 0; i < nScopeSize && (fEq || fCompl); i++) {
        word value = t[nScopeSize * index1 + i] ^ t[nScopeSize * index2 + i];
        word cvalue = caret[nScopeSize * index1 + i] & caret[nScopeSize * index2 + i];
        fEq &= !(value & cvalue);
        fCompl &= !(~value & cvalue);
      }
    } else {
      word value = GetValue(index1, lev) ^ GetValue(index2, lev);
      word cvalue = GetCare(index1, lev) & GetCare(index2, lev);
      fEq &= !(value & cvalue);
      fCompl &= !((value ^ ones[logwidth]) & cvalue);
    }
    return 2 * fCompl + fEq;
  }

  void MergeCare(int index1, int index2, int lev) {
    assert(index1 >= 0);
    assert(index2 >= 0);
    int logwidth = nInputs - lev;
    if(logwidth > lww) {
      int nScopeSize = 1 << (logwidth - lww);
      for(int i = 0; i < nScopeSize; i++) {
        caret[nScopeSize * index1 + i] |= caret[nScopeSize * index2 + i];
      }
    } else {
      word value = GetCare(index2, lev);
      int index = index1 >> (lww - logwidth);
      int pos = (index1 % (1 << (lww - logwidth))) << logwidth;
      caret[index] |= value << pos;
    }
  }

  void Merge(int index1, int index2, int lev, bool fCompl) {
    MergeCare(index1, index2, lev);
    vvMergedIndices[lev].push_back({(index1 << 1) ^ fCompl, index2});
  }

  int BDDBuildOne(int index, int lev) override {
    int r = BDDFind(index, lev);
    if(r >= -2) {
      if(r >= 0) {
        Merge(vvIndices[lev][r >> 1], index, lev, r & 1);
      }
      return r;
    }
    vvIndices[lev].push_back(index);
    return (vvIndices[lev].size() - 1) << 1;
  }

  void CompleteMerge() {
    for(int i = nInputs - 1; i >= 0; i--) {
      for(auto it = vvMergedIndices[i].rbegin(); it != vvMergedIndices[i].rend(); it++) {
        CopyFunc((*it).second, (*it).first >> 1, i, (*it).first & 1);
      }
    }
  }

  void BDDBuildStartup() override {
    RestoreCare();
    vvIndices.clear();
    vvIndices.resize(nInputs);
    vvRedundantIndices.clear();
    vvRedundantIndices.resize(nInputs);
    vvMergedIndices.clear();
    vvMergedIndices.resize(nInputs);
    for(int i = 0; i < nOutputs; i++) {
      if(!IsDC(i, 0)) {
        BDDBuildOne(i, 0);
      }
    }
  }

  virtual void BDDRebuildByMerge(int lev) {
    for(auto &p: vvMergedIndices[lev]) {
      MergeCare(p.first >> 1, p.second, lev);
    }
  }

  int BDDRebuild(int lev) override {
    RestoreCare();
    for(int i = lev; i < nInputs; i++) {
      vvIndices[i].clear();
      vvMergedIndices[i].clear();
      if(i) {
        vvRedundantIndices[i-1].clear();
      }
    }
    for(int i = 0; i < lev; i++) {
      BDDRebuildByMerge(i);
    }
    for(int i = lev; i < nInputs; i++) {
      if(!i) {
        for(int j = 0; j < nOutputs; j++) {
          if(!IsDC(j, 0)) {
            BDDBuildOne(j, 0);
          }
        }
      } else {
        BDDBuildLevel(i);
      }
    }
    return BDDNodeCount();
  }

  int BDDSwap(int lev) override {
    Swap(lev);
    return BDDRebuild(lev);
  }

  void OptimizationStartup() {
    originalt = t;
    RestoreCare();
    vvIndices.clear();
    vvIndices.resize(nInputs);
    vvMergedIndices.clear();
    vvMergedIndices.resize(nInputs);
    for(int i = 0; i < nOutputs; i++) {
      if(!IsDC(i, 0)) {
        BDDBuildOne(i, 0);
      } else {
        ShiftToMajority(i, 0);
      }
    }
  }

  virtual void Optimize() {
    OptimizationStartup();
    for(int i = 1; i < nInputs; i++) {
      for(int index: vvIndices[i-1]) {
        BDDBuildOne(index << 1, i);
        BDDBuildOne((index << 1) ^ 1, i);
      }
    }
    CompleteMerge();
  }
};

class TruthTableCareReduce : public TruthTableCare {
public:
  std::vector<std::vector<int> > vvChildren;
  std::vector<std::vector<std::pair<int, int> > > vvSkippedIndices;
  std::pair<int, int> empty = {-1, -1};

  std::vector<std::vector<std::vector<int> > > vvChildrenSaved;
  std::vector<std::vector<std::vector<std::pair<int, int> > > > vvSkippedIndicesSaved;

  TruthTableCareReduce(std::vector<std::vector<int> > const &onsets, int nInputs, std::vector<char *> const &pBPats, int nBPats, int rarity): TruthTableCare(onsets, nInputs, pBPats, nBPats, rarity) {}

  void SaveIndices(uint i) override {
    TruthTableCare::SaveIndices(i);
    if(vvChildrenSaved.size() < i + 1) {
      vvChildrenSaved.resize(i + 1);
      vvSkippedIndicesSaved.resize(i + 1);
    }
    vvChildrenSaved[i] = vvChildren;
    vvSkippedIndicesSaved[i] = vvSkippedIndices;
  }

  void LoadIndices(uint i) override {
    TruthTableCare::LoadIndices(i);
    vvChildren = vvChildrenSaved[i];
    vvSkippedIndices = vvSkippedIndicesSaved[i];
  }

  void BDDBuildStartup() override {
    vvChildren.clear();
    vvChildren.resize(nInputs);
    vvSkippedIndices.clear();
    vvSkippedIndices.resize(nInputs);
    vvSkippedIndices[nInputs-1].push_back(empty);
    TruthTableCare::BDDBuildStartup();
  }

  void BDDBuildLevel(int lev) override {
    for(int index: vvIndices[lev-1]) {
      int cof0 = BDDBuildOne(index << 1, lev);
      int cof1 = BDDBuildOne((index << 1) ^ 1, lev);
      vvChildren[lev-1].push_back(cof0);
      vvChildren[lev-1].push_back(cof1);
    }
  }

  int BDDBuild() override {
    BDDBuildStartup();
    for(int i = 1; i < nInputs; i++) {
      BDDBuildLevel(i);
    }
    BDDReduce(nInputs - 2);
    return BDDNodeCount();
  }

  int BDDRebuild(int lev) override {
    for(int i = lev; i < nInputs; i++) {
      if(i) {
        vvChildren[i-1].clear();
      }
    }
    TruthTableCare::BDDRebuild(lev);
    BDDReduce(nInputs - 2);
    return BDDNodeCount();
  }

  void BDDReduce(int lev) {
    if(lev > nInputs - 2) {
      lev = nInputs - 2;
    }
    for(int i = lev; i >= 0; i--) {
      vvRedundantIndices[i].clear();
      vvSkippedIndices[i].clear();
      vvSkippedIndices[i].resize(vvIndices[i].size(), empty);
      std::unordered_map<std::pair<std::pair<int, int>, std::pair<int, int> >, int> unique;
      unique.reserve(2 * vvIndices[i].size());
      for(uint j = 0; j < vvIndices[i].size(); j++) {
        std::pair<int, int> cof0, cof1;
        int cof0index = vvChildren[i][j+j] >> 1;
        int cof1index = vvChildren[i][j+j+1] >> 1;
        bool cof0c = vvChildren[i][j+j] & 1;
        bool cof1c = vvChildren[i][j+j+1] & 1;
        if(cof0index < 0) {
          cof0 = {nInputs, cof0c};
        } else if(vvSkippedIndices[i+1][cof0index] != empty) {
          cof0 = vvSkippedIndices[i+1][cof0index];
          cof0.second ^= cof0c;
        } else {
          cof0 = {i+1, vvChildren[i][j+j]};
        }
        if(cof1index < 0) {
          cof1 = {nInputs, cof1c};
        } else if(vvSkippedIndices[i+1][cof1index] != empty) {
          cof1 = vvSkippedIndices[i+1][cof1index];
          cof1.second ^= cof1c;
        } else {
          cof1 = {i+1, vvChildren[i][j+j+1]};
        }
        if(cof0 == cof1) {
          vvSkippedIndices[i][j] = cof0;
          vvRedundantIndices[i].push_back(vvIndices[i][j]);
          continue;
        }
        bool fCompl = cof0.second & 1;
        if(fCompl) {
          cof0.second ^= 1;
          cof1.second ^= 1;
        }
        if(unique.count({cof0, cof1})) {
          vvSkippedIndices[i][j] = {i, unique[{cof0, cof1}] ^ fCompl};
          vvRedundantIndices[i].push_back(vvIndices[i][j]);
          continue;
        }
        unique[{cof0, cof1}] = (j << 1) ^ fCompl;
      }
    }
  }
};

class TruthTableOSDM : public TruthTableCareReduce {
public:
  TruthTableOSDM(std::vector<std::vector<int> > const &onsets, int nInputs, std::vector<char *> const &pBPats, int nBPats, int rarity): TruthTableCareReduce(onsets, nInputs, pBPats, nBPats, rarity) {}

  void BDDBuildLevel(int lev) override {
    for(int index: vvIndices[lev-1]) {
      int cof0index = index << 1;
      int cof1index = cof0index ^ 1;
      int cof0, cof1;
      if(IsDC(cof0index, lev)) {
        cof1 = BDDBuildOne(cof1index, lev);
        cof0 = cof1;
      } else if(IsDC(cof1index, lev)) {
        cof0 = BDDBuildOne(cof0index, lev);
        cof1 = cof0;
      } else {
        cof0 = BDDBuildOne(cof0index, lev);
        cof1 = BDDBuildOne(cof1index, lev);
      }
      vvChildren[lev-1].push_back(cof0);
      vvChildren[lev-1].push_back(cof1);
    }
  }

  int BDDBuildOneNoCare(int index, int lev) {
    int r = BDDFind(index, lev);
    if(r >= 0 && vvIndices[lev][r >> 1] != index) {
      vvMergedIndices[lev].push_back({(vvIndices[lev][r >> 1] << 1) ^ (r & 1), index});
    }
    return r;
  }

  void BDDBuildLevelNoCare(int lev) {
    for(int index: vvIndices[lev-1]) {
      int cof0index = index << 1;
      int cof1index = cof0index ^ 1;
      int cof0, cof1;
      if(IsDC(cof0index, lev)) {
        cof1 = BDDBuildOneNoCare(cof1index, lev);
        cof0 = cof1;
      } else if(IsDC(cof1index, lev)) {
        cof0 = BDDBuildOneNoCare(cof0index, lev);
        cof1 = cof0;
      } else {
        cof0 = BDDBuildOneNoCare(cof0index, lev);
        cof1 = BDDBuildOneNoCare(cof1index, lev);
      }
      vvChildren[lev-1].push_back(cof0);
      vvChildren[lev-1].push_back(cof1);
    }
  }

  int BDDRebuild(int lev) override {
    RestoreCare();
    for(int i = lev; i < lev + 2; i++) {
      vvIndices[i].clear();
      vvMergedIndices[i].clear();
      if(i) {
        vvChildren[i-1].clear();
      }
    }
    for(int i = 0; i < lev; i++) {
      BDDRebuildByMerge(i);
    }
    for(int i = lev; i < lev + 2; i++) {
      if(!i) {
        for(int j = 0; j < nOutputs; j++) {
          if(!IsDC(j, 0)) {
            BDDBuildOne(j, 0);
          }
        }
      } else {
        BDDBuildLevel(i);
      }
    }
    if(lev < nInputs - 2) {
      vvMergedIndices[lev+2].clear();
      vvChildren[lev+1].clear();
      BDDBuildLevelNoCare(lev + 2);
    }
    BDDReduce(lev + 1);
    return BDDNodeCount();
  }

  int BDDSwap(int lev) override {
    for(int i = lev + 2; i < nInputs; i++) {
      for(auto &p: vvMergedIndices[i]) {
        if(p.first >= 0) {
          SwapIndex(p.first, i - (lev + 2) + 1);
        }
        SwapIndex(p.second, i - (lev + 2));
      }
    }
    return TruthTable::BDDSwap(lev);
  }

  void Optimize() override {
    OptimizationStartup();
    for(int i = 1; i < nInputs; i++) {
      for(int index: vvIndices[i-1]) {
        int cof0index = index << 1;
        int cof1index = cof0index ^ 1;
        if(IsDC(cof0index, i)) {
          vvMergedIndices[i].push_back({cof1index << 1, cof0index});
          BDDBuildOne(cof1index, i);
        } else if(IsDC(cof1index, i)) {
          vvMergedIndices[i].push_back({cof0index << 1, cof1index});
          BDDBuildOne(cof0index, i);
        } else {
          BDDBuildOne(cof0index, i);
          BDDBuildOne(cof1index, i);
        }
      }
    }
    CompleteMerge();
  }
};

class TruthTableOSM : public TruthTableCareReduce {
public:
  bool fComplOSM;

  TruthTableOSM(std::vector<std::vector<int> > const &onsets, int nInputs, std::vector<char *> const &pBPats, int nBPats, int rarity, bool fComplOSM = true): TruthTableCareReduce(onsets, nInputs, pBPats, nBPats, rarity), fComplOSM(fComplOSM) {}

  void BDDBuildLevel(int lev) override {
    for(int index: vvIndices[lev-1]) {
      int cof0index = index << 1;
      int cof1index = cof0index ^ 1;
      int cof0, cof1;
      if(int r = Include(cof0index, cof1index, lev, fComplOSM)) {
        cof0 = BDDBuildOne(cof0index, lev);
        cof1 = cof0 ^ !(r & 1);
      } else if(int r = Include(cof1index, cof0index, lev, fComplOSM)) {
        cof1 = BDDBuildOne(cof1index, lev);
        cof0 = cof1 ^ !(r & 1);
      } else {
        cof0 = BDDBuildOne(cof0index, lev);
        cof1 = BDDBuildOne(cof1index, lev);
      }
      vvChildren[lev-1].push_back(cof0);
      vvChildren[lev-1].push_back(cof1);
    }
  }

  void Optimize() override {
    OptimizationStartup();
    for(int i = 1; i < nInputs; i++) {
      for(int index: vvIndices[i-1]) {
        int cof0index = index << 1;
        int cof1index = cof0index ^ 1;
        if(int r = Include(cof0index, cof1index, i, fComplOSM)) {
          vvMergedIndices[i].push_back({(cof0index << 1) ^ !(r & 1), cof1index});
          BDDBuildOne(cof0index, i);
        } else if(int r = Include(cof1index, cof0index, i, fComplOSM)) {
          vvMergedIndices[i].push_back({(cof1index << 1) ^ !(r & 1), cof0index});
          BDDBuildOne(cof1index, i);
        } else {
          BDDBuildOne(cof0index, i);
          BDDBuildOne(cof1index, i);
        }
      }
    }
    CompleteMerge();
  }
};

class TruthTableTSM : public TruthTableCareReduce {
public:
  bool fComplTSM;

  TruthTableTSM(std::vector<std::vector<int> > const &onsets, int nInputs, std::vector<char *> const &pBPats, int nBPats, int rarity, bool fComplTSM = true): TruthTableCareReduce(onsets, nInputs, pBPats, nBPats, rarity), fComplTSM(fComplTSM) {}

  void BDDBuildLevel(int lev) override {
    for(int index: vvIndices[lev-1]) {
      int cof0index = index << 1;
      int cof1index = cof0index ^ 1;
      int cof0, cof1;
      if(int r = Intersect(cof0index, cof1index, lev, fComplTSM)) {
        CopyFuncMasked(cof0index, cof1index, lev, !(r & 1));
        Merge(cof0index, cof1index, lev, !(r & 1));
        cof0 = BDDBuildOne(cof0index, lev);
        cof1 = cof0 ^ !(r & 1);
      } else {
        cof0 = BDDBuildOne(cof0index, lev);
        cof1 = BDDBuildOne(cof1index, lev);
      }
      vvChildren[lev-1].push_back(cof0);
      vvChildren[lev-1].push_back(cof1);
    }
  }

  int BDDBuild() override {
    TruthTable::Save(3);
    int r = TruthTableCareReduce::BDDBuild();
    TruthTable::Load(3);
    return r;
  }

  void BDDRebuildByMerge(int lev) override {
    for(auto &p: vvMergedIndices[lev]) {
      CopyFuncMasked(p.first >> 1, p.second, lev, p.first & 1);
      MergeCare(p.first >> 1, p.second, lev);
    }
  }

  int BDDRebuild(int lev) override {
    TruthTable::Save(3);
    int r = TruthTableCareReduce::BDDRebuild(lev);
    TruthTable::Load(3);
    return r;
  }

  void Optimize() override {
    OptimizationStartup();
    for(int i = 1; i < nInputs; i++) {
      for(int index: vvIndices[i-1]) {
        int cof0index = index << 1;
        int cof1index = cof0index ^ 1;
        if(int r = Intersect(cof0index, cof1index, i, fComplTSM)) {
          CopyFuncMasked(cof0index, cof1index, i, !(r & 1));
          Merge(cof0index, cof1index, i, !(r & 1));
          BDDBuildOne(cof0index, i);
        } else {
          BDDBuildOne(cof0index, i);
          BDDBuildOne(cof1index, i);
        }
      }
    }
    CompleteMerge();
  }
};

class TruthTableLevelTSM : public TruthTableCare {
public:
  TruthTableLevelTSM(std::vector<std::vector<int> > const &onsets, int nInputs, std::vector<char *> const &pBPats, int nBPats, int rarity): TruthTableCare(onsets, nInputs, pBPats, nBPats, rarity) {}

  int BDDFindTSM(int index, int lev) {
    int logwidth = nInputs - lev;
    if(logwidth > lww) {
      int nScopeSize = 1 << (logwidth - lww);
      bool fZero = true;
      bool fOne = true;
      for(int i = 0; i < nScopeSize && (fZero || fOne); i++) {
        word value = t[nScopeSize * index + i];
        word cvalue = caret[nScopeSize * index + i];
        fZero &= !(value & cvalue);
        fOne &= !(~value & cvalue);
      }
      if(fZero || fOne) {
        return -2 ^ fOne;
      }
      for(int index2: vvIndices[lev]) {
        bool fEq = true;
        bool fCompl = true;
        for(int i = 0; i < nScopeSize && (fEq || fCompl); i++) {
          word value = t[nScopeSize * index + i] ^ t[nScopeSize * index2 + i];
          word cvalue = caret[nScopeSize * index + i] & caret[nScopeSize * index2 + i];
          fEq &= !(value & cvalue);
          fCompl &= !(~value & cvalue);
        }
        if(fEq || fCompl) {
          return (index2 << 1) ^ !fEq;
        }
      }
    } else {
      word one = ones[logwidth];
      word value = GetValue(index, lev);
      word cvalue = GetCare(index, lev);
      if(!(value & cvalue)) {
        return -2;
      }
      if(!((value ^ one) & cvalue)) {
        return -1;
      }
      for(int index2: vvIndices[lev]) {
        word value2 = value ^ GetValue(index2, lev);
        word cvalue2 = cvalue & GetCare(index2, lev);
        if(!(value2 & cvalue2)) {
          return index2 << 1;
        }
        if(!((value2 ^ one) & cvalue2)) {
          return (index2 << 1) ^ 1;
        }
      }
    }
    return -3;
  }

  int BDDBuildOne(int index, int lev) override {
    int r = BDDFindTSM(index, lev);
    if(r >= -2) {
      if(r >= 0) {
        CopyFuncMasked(r >> 1, index, lev, r & 1);
        Merge(r >> 1, index, lev, r & 1);
      } else {
        vvMergedIndices[lev].push_back({r, index});
      }
      return r;
    }
    vvIndices[lev].push_back(index);
    return index << 1;
  }

  int BDDBuild() override {
    TruthTable::Save(3);
    int r = TruthTable::BDDBuild();
    TruthTable::Load(3);
    return r;
  }

  void BDDRebuildByMerge(int lev) override {
    for(auto &p: vvMergedIndices[lev]) {
      if(p.first >= 0) {
        CopyFuncMasked(p.first >> 1, p.second, lev, p.first & 1);
        MergeCare(p.first >> 1, p.second, lev);
      }
    }
  }

  int BDDRebuild(int lev) override {
    TruthTable::Save(3);
    int r = TruthTableCare::BDDRebuild(lev);
    TruthTable::Load(3);
    return r;
  }
};

void TTTest(std::vector<std::vector<int> > const &onsets, std::vector<char *> const &pBPats, int nBPats, int rarity, std::vector<std::string> const &inputs, std::vector<std::string> const &outputs, std::ofstream &f) {
  int nInputs = inputs.size();
  if(!rarity) {
    TruthTable tt(onsets, nInputs);
    tt.RandomSiftReo(20);
    tt.BDDGenerateBlif(inputs, outputs, f);
  } else {
    TruthTableLevelTSM tt(onsets, nInputs, pBPats, nBPats, rarity);
    tt.RandomSiftReo(20);
    tt.Optimize();
    tt.BDDGenerateBlif(inputs, outputs, f);
  }

  // TruthTableReo tt(onsets, nInputs);
  // tt.RandomSiftReo(20);
  // TruthTable tt2(onsets, nInputs);
  // tt2.Reo(tt.vLevels);
  // tt2.BDDGenerateBlif(inputs, outputs, f);

  // std::vector<int> vLevels;
  // {
  //   TruthTable tt(onsets, nInputs);
  //   tt.RandomSiftReo(20);
  //   vLevels = tt.vLevels;
  // }
  // TruthTableOSM tt(onsets, nInputs, pBPats, nBPats, rarity, false);
  // tt.Reo(vLevels);
  // tt.Optimize();
  // tt.BDDGenerateBlif(inputs, outputs, f);

  // TruthTableCare tt(onsets, nInputs, pBPats, nBPats, rarity);
  // TruthTableCareReduce tt(onsets, nInputs, pBPats, nBPats, rarity);
  // TruthTableOSDM tt(onsets, nInputs, pBPats, nBPats, rarity);
  // TruthTableOSM tt(onsets, nInputs, pBPats, nBPats, rarity);
  // TruthTableTSM tt(onsets, nInputs, pBPats, nBPats, rarity);
  // TruthTableLevelTSM tt(onsets, nInputs, pBPats, nBPats, rarity);
  // tt.RandomSiftReo(20);
  // tt.Optimize();
  // tt.BDDGenerateBlif(inputs, outputs, f);

  // TruthTableOSM tt1(onsets, nInputs, pBPats, nBPats, rarity, false);
  // int r1 = tt1.RandomSiftReo(20);
  // TruthTableOSM tt2(onsets, nInputs, pBPats, nBPats, rarity, true);
  // int r2 = tt2.RandomSiftReo(20);
  // TruthTableCare &tt = (r1 < r2)? tt1: tt2;
  // tt.Optimize();
  // tt.BDDGenerateBlif(inputs, outputs, f);

  // std::cout << tt.BDDBuild() << std::endl;
  // tt.GeneratePlaCare("care.pla");
  // tt.GeneratePlaMasked("test.pla");
  // tt.BDDGenerateBlif(inputs, outputs, f);
}
