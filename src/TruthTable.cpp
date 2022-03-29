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
  std::vector<int> vLevels;

  std::vector<std::vector<word> > savedt;
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

  virtual void Save(uint i) {
    if(savedt.size() < i + 1) {
      savedt.resize(i + 1);
      vLevelsSaved.resize(i + 1);
    }
    savedt[i] = t;
    vLevelsSaved[i] = vLevels;
  }

  virtual void Load(uint i) {
    t = savedt[i];
    vLevels = vLevelsSaved[i];
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

  word GetValue(int index_lev, int lev) {
    assert(index_lev >= 0);
    assert(nInputs - lev <= lww);
    int logwidth = nInputs - lev;
    int index = index_lev >> (lww - logwidth);
    int pos = (index_lev % (1 << (lww - logwidth))) << logwidth;
    return (t[index] >> pos) & ones[logwidth];
  }

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
      for(int i = 0; i < nScopeSize; i++) {
        if(index2 < 0) {
          t[nScopeSize * index1 + i] = 0;
        } else {
          t[nScopeSize * index1 + i] = t[nScopeSize * index2 + i];
        }
        if(fCompl) {
          t[nScopeSize * index1 + i] = ~t[nScopeSize * index1 + i];
        }
      }
    } else {
      word value;
      if(index2 < 0) {
        value = 0;
      } else {
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
      for(int index2: vvIndices[lev]) {
        bool fEq = true;
        bool fCompl = true;
        for(int i = 0; i < nScopeSize && (fEq || fCompl); i++) {
          word value = t[nScopeSize * index + i];
          word value2 = t[nScopeSize * index2 + i];
          fEq &= value == value2;
          fCompl &= value == ~value2;
        }
        if(fEq || fCompl) {
          return (index2 << 1) ^ fCompl;
        }
      }
    } else {
      word one = ones[logwidth];
      word value = GetValue(index, lev);
      if(!value) {
        return -2;
      }
      if(!(value ^ one)) {
        return -1;
      }
      for(int index2: vvIndices[lev]) {
        word value2 = value ^ GetValue(index2, lev);
        if(!(value2)) {
          return index2 << 1;
        }
        if(!(value2 ^ one)) {
          return (index2 << 1) ^ 1;
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
    return index << 1;
  }

  void BDDRemoveRedundantIndices(std::vector<std::vector<int> > const &vvIndicesRedundant) {
    for(int i = 0; i < nInputs; i++) {
      auto it = vvIndicesRedundant[i].begin();
      std::vector<int> vIndicesNew;
      for(int j: vvIndices[i]) {
        if(it == vvIndicesRedundant[i].end() || j != *it) {
          vIndicesNew.push_back(j);
        } else {
          it++;
        }
      }
      vvIndices[i] = vIndicesNew;
    }
  }

  int BDDNodeCount() {
    int count = 1; // const node
    for(int i = 0; i < nInputs; i++) {
      count += vvIndices[i].size();
    }
    return count;
  }

  virtual void BDDBuildStartup() {
    vvIndices.clear();
    vvIndices.resize(nInputs);
    for(int i = 0; i < nOutputs; i++) {
      BDDBuildOne(i, 0);
    }
  }

  virtual int BDDBuild() {
    BDDBuildStartup();
    std::vector<std::vector<int> > vvIndicesRedundant(nInputs);
    for(int i = 1; i < nInputs; i++) {
      for(int index: vvIndices[i-1]) {
        int cof0 = BDDBuildOne(index << 1, i);
        int cof1 = BDDBuildOne((index << 1) ^ 1, i);
        if(cof0 == cof1) {
          vvIndicesRedundant[i-1].push_back(index);
        }
      }
    }
    BDDRemoveRedundantIndices(vvIndicesRedundant);
    return BDDNodeCount();
  }

  bool Imply(int index1, int index2, int lev) {
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

  int BDDGenerateBlifRec(std::vector<std::vector<int> > &vvNodes, int &nNodes, int index, int lev, std::ofstream &f, std::string const &prefix) {
    int r = BDDFind(index, lev);
    if(r >= 0) {
      auto it = std::lower_bound(vvIndices[lev].begin(), vvIndices[lev].end(), r >> 1);
      int i = it - vvIndices[lev].begin();
      return (vvNodes[lev][i] << 1) ^ (r & 1);
    }
    if(r >= -2) {
      return r + 2;
    }
    int cof0 = BDDGenerateBlifRec(vvNodes, nNodes, index << 1, lev + 1, f, prefix);
    int cof1 = BDDGenerateBlifRec(vvNodes, nNodes, (index << 1) ^ 1, lev + 1, f, prefix);
    if(cof0 == cof1) {
      return cof0;
    }
    int cof0id = cof0 >> 1;
    int cof1id = cof1 >> 1;
    bool cof0c = cof0 & 1;
    bool cof1c = cof1 & 1;
    bool imp01 = Imply(index << 1, (index << 1) ^ 1, lev + 1);
    bool imp10 = Imply((index << 1) ^ 1, index << 1, lev + 1);
    f << ".names " << prefix << "v" << lev << " " << prefix << "n" << cof0id << " " << prefix << "n" << cof1id << " " << prefix << "n" << nNodes << std::endl;
    f << (imp01? "-" : "0") << !cof0c << "- 1" << std::endl;
    f << (imp10? "--" : "1-") << !cof1c << " 1" << std::endl;
    vvIndices[lev].push_back(index);
    vvNodes[lev].push_back(nNodes);
    return (nNodes++) << 1;
  }

  void BDDGenerateBlif(std::vector<std::string> const &inputs, std::vector<std::string> const &outputs, std::ofstream &f) {
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
      int id = node >> 1;
      bool c = node & 1;
      f << ".names " << prefix << "n" << id << " " << outputs[i] << std::endl;
      f << !c << " 1" << std::endl;
    }
  }

  virtual void SwapLevel(int lev) {
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

  int SiftReo() {
    int best = BDDBuild();
    Save(0);
    std::list<int> vars(nInputs);
    std::iota(vars.begin(), vars.end(), 0);
    while(!vars.empty()) {
      int maxvar = -1;
      uint maxnodes = 0;
      std::list<int>::iterator maxit;
      for(auto it = vars.begin(); it != vars.end(); it++) {
        if(vvIndices[vLevels[maxvar]].size() > maxnodes) {
          maxnodes = vvIndices[vLevels[maxvar]].size();
          maxvar = *it;
          maxit = it;
        }
      }
      if(maxvar == -1) {
        break;
      }
      vars.erase(maxit);
      Save(1);
      for(int i = vLevels[maxvar]; i < nInputs - 1; i++) {
        SwapLevel(i);
        int count = BDDBuild();
        if(best > count) {
          best = count;
          Save(0);
        }
      }
      Load(1);
      for(int i = vLevels[maxvar]-1; i >= 0; i--) {
        SwapLevel(i);
        int count = BDDBuild();
        if(best > count) {
          best = count;
          Save(0);
        }
      }
      Load(0);
    }
    return best;
  }

  void Reo(std::vector<int> vLevelsNew) {
    for(int i = 0; i < nInputs; i++) {
      int var = std::find(vLevelsNew.begin(), vLevelsNew.end(), i) - vLevelsNew.begin();
      int lev = vLevels[var];
      if(lev < i) {
        for(int j = lev; j < i; j++) {
          SwapLevel(j);
        }
      } else if(lev > i) {
        for(int j = lev-1; j >= i; j--) {
          SwapLevel(j);
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

  void ShowIndices() {
    for(uint i = 0; i < vvIndices.size(); i++) {
      std::cout << "var " << i << ":" << std::endl;
      for(int index: vvIndices[i]) {
        std::cout << index << ", ";
      }
      std::cout << std::endl;
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

class TruthTableCare : public TruthTable{
public:
  std::vector<word> originalt;
  std::vector<word> caret;
  std::vector<word> care;

  std::vector<std::vector<std::pair<int, int> > > vvIndicesMerged;

  std::vector<std::vector<word> > savedcare;

  TruthTableCare(std::vector<std::vector<int> > const &onsets, int nInputs, std::vector<char *> const &pBPats, int nBPats, int rarity): TruthTable(onsets, nInputs) {
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
    RestoreCare();
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
    RestoreCare();
  }

  void SwapLevel(int lev) override {
    TruthTable::SwapLevel(lev);
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
    RestoreCare();
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

  word GetCare(int index_lev, int lev) {
    assert(index_lev >= 0);
    assert(nInputs - lev <= lww);
    int logwidth = nInputs - lev;
    int index = index_lev >> (lww - logwidth);
    int pos = (index_lev % (1 << (lww - logwidth))) << logwidth;
    return (caret[index] >> pos) & ones[logwidth];
  }

  bool IsDC(int index, int lev) {
    if(nInputs - lev > lww) {
      int nScopeSize = 1 << (nInputs - lev - lww);
      for(int i = 0; i < nScopeSize; i++) {
        word value = caret[nScopeSize * index + i];
        if(value != 0) {
          return false;
        }
      }
      return true;
    }
    word value = GetCare(index, lev);
    if(value != 0) {
      return false;
    }
    return true;
  }

  void MergeCare(int index1, int index2, int lev) {
    assert(index2 >= 0);
    if(index1 < 0) {
      return;
    }
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

  int BDDBuildOne(int index, int lev) override {
    int r = BDDFind(index, lev);
    if(r >= -2) {
      if(r >= 0) {
        MergeCare(r >> 1, index, lev);
      }
      if(!vvIndicesMerged.empty()) {
        vvIndicesMerged[lev].push_back({r, index});
      }
      return r;
    }
    vvIndices[lev].push_back(index);
    return index << 1;
  }

  void BDDRemoveRedundantIndicesFromChildren(std::vector<std::vector<int> > const &vvChildren) {
    std::vector<std::vector<int> > vvIndicesRedundant(nInputs);
    std::map<int, std::pair<int, int> > skipped;
    for(int i = nInputs - 2; i >= 0; i--) {
      std::map<int, std::pair<int, int> > nextskipped;
      std::unordered_map<std::pair<std::pair<int, int>, std::pair<int, int> >, int> unique;
      unique.reserve(2 * vvIndices[i].size());
      for(uint j = 0; j < vvIndices[i].size(); j++) {
        std::pair<int, int> cof0, cof1;
        int cof0index = vvChildren[i][j+j] >> 1;
        int cof1index = vvChildren[i][j+j+1] >> 1;
        bool cof0c = vvChildren[i][j+j] & 1;
        bool cof1c = vvChildren[i][j+j+1] & 1;
        if(cof0index < 0) {
          cof0 = {nInputs, vvChildren[i][j+j]};
        } else if(skipped.count(cof0index)) {
          cof0 = skipped[cof0index];
          cof0.second ^= cof0c;
        } else {
          cof0 = {i+1, vvChildren[i][j+j]};
        }
        if(cof1index < 0) {
          cof1 = {nInputs, vvChildren[i][j+j+1]};
        } else if(skipped.count(cof1index)) {
          cof1 = skipped[cof1index];
          cof1.second ^= cof1c;
        } else {
          cof1 = {i+1, vvChildren[i][j+j+1]};
        }
        if(cof0 == cof1) {
          nextskipped[vvIndices[i][j]] = cof0;
          vvIndicesRedundant[i].push_back(vvIndices[i][j]);
          continue;
        }
        bool fCompl = cof0.second & 1;
        if(fCompl) {
          cof0.second ^= 1;
          cof1.second ^= 1;
        }
        if(unique.count({cof0, cof1})) {
          nextskipped[vvIndices[i][j]] = {i, unique[{cof0, cof1}] ^ fCompl};
          vvIndicesRedundant[i].push_back(vvIndices[i][j]);
          continue;
        }
        unique[{cof0, cof1}] = (vvIndices[i][j] << 1) ^ fCompl;
      }
      skipped = nextskipped;
    }
    BDDRemoveRedundantIndices(vvIndicesRedundant);
  }

  void Merge() {
    for(int i = nInputs - 1; i >= 0; i--) {
      for(auto it = vvIndicesMerged[i].rbegin(); it != vvIndicesMerged[i].rend(); it++) {
        CopyFunc((*it).second, (*it).first >> 1, i, (*it).first & 1);
      }
    }
    vvIndicesMerged.clear();
  }

  void BDDBuildStartup() override {
    vvIndices.clear();
    vvIndices.resize(nInputs);
    for(int i = 0; i < nOutputs; i++) {
      if(!IsDC(i, 0)) {
        BDDBuildOne(i, 0);
      }
    }
  }

  void OptimizationStartup() {
    originalt = t;
    vvIndices.clear();
    vvIndices.resize(nInputs);
    vvIndicesMerged.resize(nInputs);
    for(int i = 0; i < nOutputs; i++) {
      if(!IsDC(i, 0)) {
        BDDBuildOne(i, 0);
      } else {
        ShiftToMajority(i, 0);
      }
    }
  }

  int Include(int index1, int index2, int lev, bool fCompl) {
    assert(index1 >= 0);
    assert(index2 >= 0);
    int logwidth = nInputs - lev;
    bool fEq = true;
    if(logwidth > lww) {
      int nScopeSize = 1 << (logwidth - lww);
      for(int i = 0; i < nScopeSize && (fEq || fCompl); i++) {
        word value = t[nScopeSize * index1 + i] ^ t[nScopeSize * index2 + i];
        word cvalue = caret[nScopeSize * index2 + i];
        fEq &= !(value & cvalue);
        fCompl &= !(~value & cvalue);
      }
    } else {
      word value = GetValue(index1, lev) ^ GetValue(index2, lev);
      word cvalue = GetCare(index2, lev);
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

  virtual void Optimize() = 0;
};

class TruthTableOSDM : public TruthTableCare{
public:
  TruthTableOSDM(std::vector<std::vector<int> > const &onsets, int nInputs, std::vector<char *> const &pBPats, int nBPats, int rarity): TruthTableCare(onsets, nInputs, pBPats, nBPats, rarity) {}

  int BDDBuild() override {
    BDDBuildStartup();
    std::vector<std::vector<int> > vvChildren(nInputs);
    for(int i = 1; i < nInputs; i++) {
      for(int index: vvIndices[i-1]) {
        int cof0index = index << 1;
        int cof1index = cof0index ^ 1;
        int cof0, cof1;
        if(IsDC(cof0index, i)) {
          cof1 = BDDBuildOne(cof1index, i);
          cof0 = cof1;
        } else if(IsDC(cof1index, i)) {
          cof0 = BDDBuildOne(cof0index, i);
          cof1 = cof0;
        } else {
          cof0 = BDDBuildOne(cof0index, i);
          cof1 = BDDBuildOne(cof1index, i);
        }
        vvChildren[i-1].push_back(cof0);
        vvChildren[i-1].push_back(cof1);
      }
    }
    BDDRemoveRedundantIndicesFromChildren(vvChildren);
    RestoreCare();
    return BDDNodeCount();
  }

  void Optimize() override {
    OptimizationStartup();
    for(int i = 1; i < nInputs; i++) {
      for(int index: vvIndices[i-1]) {
        int cof0index = index << 1;
        int cof1index = cof0index ^ 1;
        if(IsDC(cof0index, i)) {
          vvIndicesMerged[i].push_back({cof1index << 1, cof0index});
          BDDBuildOne(cof1index, i);
        } else if(IsDC(cof1index, i)) {
          vvIndicesMerged[i].push_back({cof0index << 1, cof1index});
          BDDBuildOne(cof0index, i);
        } else {
          BDDBuildOne(cof0index, i);
          BDDBuildOne(cof1index, i);
        }
      }
    }
    Merge();
  }
};

class TruthTableOSM : public TruthTableCare{
public:
  bool fComplOSM;

  TruthTableOSM(std::vector<std::vector<int> > const &onsets, int nInputs, std::vector<char *> const &pBPats, int nBPats, int rarity, bool fComplOSM = true): TruthTableCare(onsets, nInputs, pBPats, nBPats, rarity), fComplOSM(fComplOSM) {}

  int BDDBuild() override {
    BDDBuildStartup();
    std::vector<std::vector<int> > vvChildren(nInputs);
    for(int i = 1; i < nInputs; i++) {
      for(int index: vvIndices[i-1]) {
        int cof0index = index << 1;
        int cof1index = cof0index ^ 1;
        int cof0, cof1;
        if(int r = Include(cof0index, cof1index, i, fComplOSM)) {
          MergeCare(cof0index, cof1index, i);
          cof0 = BDDBuildOne(cof0index, i);
          cof1 = cof0 ^ !(r & 1);
        } else if(int r = Include(cof1index, cof0index, i, fComplOSM)) {
          MergeCare(cof1index, cof0index, i);
          cof1 = BDDBuildOne(cof1index, i);
          cof0 = cof1 ^ !(r & 1);
        } else {
          cof0 = BDDBuildOne(cof0index, i);
          cof1 = BDDBuildOne(cof1index, i);
        }
        vvChildren[i-1].push_back(cof0);
        vvChildren[i-1].push_back(cof1);
      }
    }
    BDDRemoveRedundantIndicesFromChildren(vvChildren);
    RestoreCare();
    return BDDNodeCount();
  }

  void Optimize() override {
    OptimizationStartup();
    for(int i = 1; i < nInputs; i++) {
      for(int index: vvIndices[i-1]) {
        int cof0index = index << 1;
        int cof1index = cof0index ^ 1;
        if(int r = Include(cof0index, cof1index, i, fComplOSM)) {
          MergeCare(cof0index, cof1index, i);
          vvIndicesMerged[i].push_back({(cof0index << 1) ^ !(r & 1), cof1index});
          BDDBuildOne(cof0index, i);
        } else if(int r = Include(cof1index, cof0index, i, fComplOSM)) {
          MergeCare(cof1index, cof0index, i);
          vvIndicesMerged[i].push_back({(cof1index << 1) ^ !(r & 1), cof0index});
          BDDBuildOne(cof1index, i);
        } else {
          BDDBuildOne(cof0index, i);
          BDDBuildOne(cof1index, i);
        }
      }
    }
    Merge();
  }
};

class TruthTableTSM : public TruthTableCare{
public:
  bool fComplTSM;

  TruthTableTSM(std::vector<std::vector<int> > const &onsets, int nInputs, std::vector<char *> const &pBPats, int nBPats, int rarity, bool fComplTSM = true): TruthTableCare(onsets, nInputs, pBPats, nBPats, rarity), fComplTSM(fComplTSM) {}

  int BDDBuild() {
    Save(3);
    BDDBuildStartup();
    std::vector<std::vector<int> > vvChildren(nInputs);
    for(int i = 1; i < nInputs; i++) {
      for(int index: vvIndices[i-1]) {
        int cof0index = index << 1;
        int cof1index = cof0index ^ 1;
        int cof0, cof1;
        if(int r = Intersect(cof0index, cof1index, i, fComplTSM)) {
          CopyFuncMasked(cof0index, cof1index, i, !(r & 1));
          MergeCare(cof0index, cof1index, i);
          cof0 = BDDBuildOne(cof0index, i);
          cof1 = cof0 ^ !(r & 1);
        } else {
          cof0 = BDDBuildOne(cof0index, i);
          cof1 = BDDBuildOne(cof1index, i);
        }
        vvChildren[i-1].push_back(cof0);
        vvChildren[i-1].push_back(cof1);
      }
    }
    BDDRemoveRedundantIndicesFromChildren(vvChildren);
    Load(3);
    return BDDNodeCount();
  }

  void Optimize() override {
    OptimizationStartup();
    for(int i = 1; i < nInputs; i++) {
      for(int index: vvIndices[i-1]) {
        int cof0index = index << 1;
        int cof1index = cof0index ^ 1;
        if(int r = Intersect(cof0index, cof1index, i, fComplTSM)) {
          CopyFuncMasked(cof0index, cof1index, i, !(r & 1));
          MergeCare(cof0index, cof1index, i);
          vvIndicesMerged[i].push_back({(cof0index << 1) ^ !(r & 1), cof1index});
          BDDBuildOne(cof0index, i);
        } else {
          BDDBuildOne(cof0index, i);
          BDDBuildOne(cof1index, i);
        }
      }
    }
    Merge();
  }
};

class TruthTableTSMNew : public TruthTableCare{
public:
  bool fComplOSM, fComplTSM;

  TruthTableTSMNew(std::vector<std::vector<int> > const &onsets, int nInputs, std::vector<char *> const &pBPats, int nBPats, int rarity, bool fComplOSM = true, bool fComplTSM = true): TruthTableCare(onsets, nInputs, pBPats, nBPats, rarity), fComplOSM(fComplOSM), fComplTSM(fComplTSM) {}

  int BDDBuild() override {
    Save(3);
    BDDBuildStartup();
    std::vector<std::vector<int> > vvChildren(nInputs);
    for(int i = 1; i < nInputs; i++) {
      for(int index: vvIndices[i-1]) {
        int cof0index = index << 1;
        int cof1index = cof0index ^ 1;
        int cof0, cof1;
        if(int r = Include(cof0index, cof1index, i, fComplOSM)) {
          MergeCare(cof0index, cof1index, i);
          cof0 = BDDBuildOne(cof0index, i);
          cof1 = cof0 ^ !(r & 1);
        } else if(int r = Include(cof1index, cof0index, i, fComplOSM)) {
          MergeCare(cof1index, cof0index, i);
          cof1 = BDDBuildOne(cof1index, i);
          cof0 = cof1 ^ !(r & 1);
        } else {
          cof0 = BDDBuildOne(cof0index, i);
          cof1 = BDDBuildOne(cof1index, i);
        }
        vvChildren[i-1].push_back(cof0);
        vvChildren[i-1].push_back(cof1);
      }
      std::map<int, int> m;
      for(uint j = 0; j < vvIndices[i-1].size(); j++) {
        int cof0 = vvChildren[i-1][j+j];
        while(m.count(cof0 >> 1) && cof0 != (m[cof0 >> 1] ^ (cof0 & 1))) {
          cof0 = m[cof0 >> 1] ^ (cof0 & 1);
        }
        int cof1 = vvChildren[i-1][j+j+1];
        while(m.count(cof1 >> 1) && cof1 != (m[cof1 >> 1] ^ (cof1 & 1))) {
          cof1 = m[cof1 >> 1] ^ (cof1 & 1);
        }
        int cof0index = cof0 >> 1;
        int cof1index = cof1 >> 1;
        if(cof0index < 0 || cof1index < 0 || cof0index == cof1index) {
          continue;
        }
        bool fComplCof = (cof0 & 1) ^ (cof1 & 1);
        if(int r = Intersect(cof0index, cof1index, i, fComplCof || fComplTSM, !fComplCof || fComplTSM)) {
          auto it = std::find(vvIndices[i].begin(), vvIndices[i].end(), cof0index);
          assert(it != vvIndices[i].end());
          vvIndices[i].erase(it);
          it = std::find(vvIndices[i].begin(), vvIndices[i].end(), cof1index);
          assert(it != vvIndices[i].end());
          vvIndices[i].erase(it);
          CopyFuncMasked(cof0index, cof1index, i, !(r & 1));
          MergeCare(cof0index, cof1index, i);
          m[cof0index] = BDDBuildOne(cof0index, i);
          m[cof1index] = m[cof0index] ^ !(r & 1);
        }
      }
      for(uint j = 0; j < vvChildren[i-1].size(); j++) {
        int cof = vvChildren[i-1][j];
        while(m.count(cof >> 1) && cof != (m[cof >> 1] ^ (cof & 1))) {
          cof = m[cof >> 1] ^ (cof & 1);
        }
        vvChildren[i-1][j] = cof;
      }
    }
    BDDRemoveRedundantIndicesFromChildren(vvChildren);
    Load(3);
    return BDDNodeCount();
  }

  void Optimize() override {
    OptimizationStartup();
    for(int i = 1; i < nInputs; i++) {
      std::vector<int> vChildren;
      for(int index: vvIndices[i-1]) {
        int cof0index = index << 1;
        int cof1index = cof0index ^ 1;
        int cof0, cof1;
        if(int r = Include(cof0index, cof1index, i, fComplOSM)) {
          MergeCare(cof0index, cof1index, i);
          vvIndicesMerged[i].push_back({(cof0index << 1) ^ !(r & 1), cof1index});
          cof0 = BDDBuildOne(cof0index, i);
          cof1 = cof0 ^ !(r & 1);
        } else if(int r = Include(cof1index, cof0index, i, fComplOSM)) {
          MergeCare(cof1index, cof0index, i);
          vvIndicesMerged[i].push_back({(cof1index << 1) ^ !(r & 1), cof0index});
          cof1 = BDDBuildOne(cof1index, i);
          cof0 = cof1 ^ !(r & 1);
        } else {
          cof0 = BDDBuildOne(cof0index, i);
          cof1 = BDDBuildOne(cof1index, i);
        }
        vChildren.push_back(cof0);
        vChildren.push_back(cof1);
      }
      std::map<int, int> m;
      for(uint j = 0; j < vvIndices[i-1].size(); j++) {
        int cof0 = vChildren[j+j];
        while(m.count(cof0 >> 1) && cof0 != (m[cof0 >> 1] ^ (cof0 & 1))) {
          cof0 = m[cof0 >> 1] ^ (cof0 & 1);
        }
        int cof1 = vChildren[j+j+1];
        while(m.count(cof1 >> 1) && cof1 != (m[cof1 >> 1] ^ (cof1 & 1))) {
          cof1 = m[cof1 >> 1] ^ (cof1 & 1);
        }
        int cof0index = cof0 >> 1;
        int cof1index = cof1 >> 1;
        if(cof0index < 0 || cof1index < 0 || cof0index == cof1index) {
          continue;
        }
        bool fComplCof = (cof0 & 1) ^ (cof1 & 1);
        if(int r = Intersect(cof0index, cof1index, i, fComplCof || fComplTSM, !fComplCof || fComplTSM)) {
          auto it = std::find(vvIndices[i].begin(), vvIndices[i].end(), cof0index);
          assert(it != vvIndices[i].end());
          vvIndices[i].erase(it);
          it = std::find(vvIndices[i].begin(), vvIndices[i].end(), cof1index);
          assert(it != vvIndices[i].end());
          vvIndices[i].erase(it);
          CopyFuncMasked(cof0index, cof1index, i, !(r & 1));
          MergeCare(cof0index, cof1index, i);
          vvIndicesMerged[i].push_back({(cof0index << 1) ^ !(r & 1), cof1index});
          m[cof0index] = BDDBuildOne(cof0index, i);
          m[cof1index] = m[cof0index] ^ !(r & 1);
        }
      }
    }
    Merge();
  }
};

class TruthTableLevelTSM : public TruthTableCare{
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
          return (index2 << 1) ^ fCompl;
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
        MergeCare(r >> 1, index, lev);
      }
      if(!vvIndicesMerged.empty()) {
        vvIndicesMerged[lev].push_back({r, index});
      }
      return r;
    }
    vvIndices[lev].push_back(index);
    return index << 1;
  }

  int BDDBuild() override {
    Save(3);
    TruthTable::BDDBuild();
    Load(3);
    return BDDNodeCount();
  }

  void Optimize() override {
    OptimizationStartup();
    for(int i = 1; i < nInputs; i++) {
      for(int index: vvIndices[i-1]) {
        BDDBuildOne(index << 1, i);
        BDDBuildOne((index << 1) ^ 1, i);
      }
    }
    Merge();
  }
};

void TTTest(std::vector<std::vector<int> > const &onsets, std::vector<char *> const &pBPats, int nBPats, int rarity, std::vector<std::string> const &inputs, std::vector<std::string> const &outputs, std::ofstream &f) {
  int nInputs = inputs.size();
  // TruthTable tt(onsets, nInputs);
  // // tt.SiftReo();
  // tt.RandomSiftReo(20);
  // tt.BDDGenerateBlif(inputs, outputs, f);

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

  TruthTableLevelTSM tt(onsets, nInputs, pBPats, nBPats, rarity);
  tt.RandomSiftReo(20);
  tt.Optimize();
  tt.BDDGenerateBlif(inputs, outputs, f);

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
