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

extern std::string BinaryToString(int bin, int size);

class TT {
public:
  typedef uint word;
  const int ww = 32; // word width
  const int lww = 5; // log word width

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

  TT(std::vector<std::vector<int> > const &onsets, int nInputs): nInputs(nInputs) {
    nOutputs = onsets.size();
    if(nInputs >= lww) {
      nSize = 1 << (nInputs - lww);
      nTotalSize = nSize * nOutputs;
      t.resize(nTotalSize);
      for(int i = 0; i < nOutputs; i++) {
        for(int pat: onsets[i]) {
          int index = pat / ww;
          int pos = pat % ww;
          t[nSize * i + index] |= 1 << pos;
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
          t[padding / ww] |= 1 << pos;
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
        count += std::bitset<32>(t[nScopeSize * index + i]).count();
      }
    } else {
      count = std::bitset<32>(GetValue(index, lev)).count();
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
      word one = ones[lww];
      for(int i = 0; i < nScopeSize && (fZero || fOne); i++) {
        word value = t[nScopeSize * index + i];
        fZero &= value == 0;
        fOne &= value == one;
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
      if(value == 0) {
        return -2;
      }
      if(value == one) {
        return -1;
      }
      for(int index2: vvIndices[lev]) {
        word value2 = GetValue(index2, lev);
        if(value2 == value) {
          return index2 << 1;
        }
        if((value2 ^ one) == value) {
          return (index2 << 1) ^ 1;
        }
      }
    }
    return -3;
  }

  int BDDCountNodesOne(int index, int lev) {
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
  
  virtual int BDDCountNodes() {
    vvIndices.clear();
    vvIndices.resize(nInputs);
    std::vector<std::vector<int> > vvIndicesRedundant(nInputs);
    for(int i = 0; i < nOutputs; i++) {
      BDDCountNodesOne(i, 0);
    }
    for(int i = 1; i < nInputs; i++) {
      for(int index: vvIndices[i-1]) {
        int cof0 = BDDCountNodesOne(index << 1, i);
        int cof1 = BDDCountNodesOne((index << 1) ^ 1, i);
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
    int best = BDDCountNodes();
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
        int count = BDDCountNodes();
        if(best > count) {
          best = count;
          Save(0);
        }
      }
      Load(1);
      for(int i = vLevels[maxvar]-1; i >= 0; i--) {
        SwapLevel(i);
        int count = BDDCountNodes();
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

const TT::word TT::ones[] = {0x00000001,
                             0x00000003,
                             0x0000000f,
                             0x000000ff,
                             0x0000ffff,
                             0xffffffff};

const TT::word TT::swapmask[] = {0x22222222,
                                 0x0c0c0c0c,
                                 0x00f000f0,
                                 0x0000ff00};

class TTCare : public TT{
public:
  std::vector<word> originalt;
  std::vector<word> caret;
  std::vector<word> care;

  std::vector<std::vector<std::pair<int, int> > > vvIndicesMerged;

  std::vector<std::vector<word> > savedcare;

  TTCare(std::vector<std::vector<int> > const &onsets, int nInputs, std::vector<char *> const &pBPats, int nBPats, int rarity): TT(onsets, nInputs) {
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
          care[index] |= 1 << pos;
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
    TT::Save(i);
    if(savedcare.size() < i + 1) {
      savedcare.resize(i + 1);
    }
    savedcare[i] = care;
  }

  void Load(uint i) override {
    TT::Load(i);
    care = savedcare[i];
    RestoreCare();
  }

  void SwapLevel(int lev) override {
    TT::SwapLevel(lev);
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

  int BDDCountNodesCareOne(int index, int lev) {
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
      std::map<std::pair<std::pair<int, int>, std::pair<int, int> >, int> unique;
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

  void BDDCountNodesStartup() {
    vvIndices.clear();
    vvIndices.resize(nInputs);
    for(int i = 0; i < nOutputs; i++) {
      if(!IsDC(i, 0)) {
        BDDCountNodesCareOne(i, 0);
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
        BDDCountNodesCareOne(i, 0);
      } else {
        ShiftToMajority(i, 0);
      }
    }
  }

  int BDDCountNodesOSDM() {
    BDDCountNodesStartup();
    std::vector<std::vector<int> > vvChildren(nInputs);
    for(int i = 1; i < nInputs; i++) {
      for(int index: vvIndices[i-1]) {
        int cof0index = index << 1;
        int cof1index = cof0index ^ 1;
        int cof0, cof1;
        if(IsDC(cof0index, i)) {
          cof1 = BDDCountNodesCareOne(cof1index, i);
          cof0 = cof1;
        } else if(IsDC(cof1index, i)) {
          cof0 = BDDCountNodesCareOne(cof0index, i);
          cof1 = cof0;
        } else {
          cof0 = BDDCountNodesCareOne(cof0index, i);
          cof1 = BDDCountNodesCareOne(cof1index, i);
        }
        vvChildren[i-1].push_back(cof0);
        vvChildren[i-1].push_back(cof1);
      }
    }
    BDDRemoveRedundantIndicesFromChildren(vvChildren);
    RestoreCare();
    return BDDNodeCount();
  }

  void OSDM() {
    OptimizationStartup();
    for(int i = 1; i < nInputs; i++) {
      for(int index: vvIndices[i-1]) {
        int cof0index = index << 1;
        int cof1index = cof0index ^ 1;
        if(IsDC(cof0index, i)) {
          vvIndicesMerged[i].push_back({cof1index << 1, cof0index});
          BDDCountNodesCareOne(cof1index, i);
        } else if(IsDC(cof1index, i)) {
          vvIndicesMerged[i].push_back({cof0index << 1, cof1index});
          BDDCountNodesCareOne(cof0index, i);
        } else {
          BDDCountNodesCareOne(cof0index, i);
          BDDCountNodesCareOne(cof1index, i);
        }
      }
    }
    Merge();
  }

  int Include(int index1, int index2, int lev, bool fCompl) {
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

  int BDDCountNodesOSM(bool fCompl) {
    BDDCountNodesStartup();
    std::vector<std::vector<int> > vvChildren(nInputs);
    for(int i = 1; i < nInputs; i++) {
      for(int index: vvIndices[i-1]) {
        int cof0index = index << 1;
        int cof1index = cof0index ^ 1;
        int cof0, cof1;
        if(int r = Include(cof0index, cof1index, i, fCompl)) {
          MergeCare(cof0index, cof1index, i);
          cof0 = BDDCountNodesCareOne(cof0index, i);
          cof1 = cof0 ^ !(r & 1);
        } else if(int r = Include(cof1index, cof0index, i, fCompl)) {
          MergeCare(cof1index, cof0index, i);
          cof1 = BDDCountNodesCareOne(cof1index, i);
          cof0 = cof1 ^ !(r & 1);
        } else {
          cof0 = BDDCountNodesCareOne(cof0index, i);
          cof1 = BDDCountNodesCareOne(cof1index, i);
        }
        vvChildren[i-1].push_back(cof0);
        vvChildren[i-1].push_back(cof1);
      }
    }
    BDDRemoveRedundantIndicesFromChildren(vvChildren);
    RestoreCare();
    return BDDNodeCount();
  }

  void OSM(bool fCompl) {
    OptimizationStartup();
    for(int i = 1; i < nInputs; i++) {
      for(int index: vvIndices[i-1]) {
        int cof0index = index << 1;
        int cof1index = cof0index ^ 1;
        if(int r = Include(cof0index, cof1index, i, fCompl)) {
          MergeCare(cof0index, cof1index, i);
          vvIndicesMerged[i].push_back({(cof0index << 1) ^ !(r & 1), cof1index});
          BDDCountNodesCareOne(cof0index, i);
        } else if(int r = Include(cof1index, cof0index, i, fCompl)) {
          MergeCare(cof1index, cof0index, i);
          vvIndicesMerged[i].push_back({(cof1index << 1) ^ !(r & 1), cof0index});
          BDDCountNodesCareOne(cof1index, i);
        } else {
          BDDCountNodesCareOne(cof0index, i);
          BDDCountNodesCareOne(cof1index, i);
        }
      }
    }
    Merge();
  }

  int Intersect(int index1, int index2, int lev, bool fCompl, bool fEq = true) {
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

  int BDDCountNodesTSM(bool fCompl) {
    Save(3);
    BDDCountNodesStartup();
    std::vector<std::vector<int> > vvChildren(nInputs);
    for(int i = 1; i < nInputs; i++) {
      for(int index: vvIndices[i-1]) {
        int cof0index = index << 1;
        int cof1index = cof0index ^ 1;
        int cof0, cof1;
        if(int r = Intersect(cof0index, cof1index, i, fCompl)) {
          CopyFuncMasked(cof0index, cof1index, i, !(r & 1));
          MergeCare(cof0index, cof1index, i);
          cof0 = BDDCountNodesCareOne(cof0index, i);
          cof1 = cof0 ^ !(r & 1);
        } else {
          cof0 = BDDCountNodesCareOne(cof0index, i);
          cof1 = BDDCountNodesCareOne(cof1index, i);
        }
        vvChildren[i-1].push_back(cof0);
        vvChildren[i-1].push_back(cof1);
      }
    }
    BDDRemoveRedundantIndicesFromChildren(vvChildren);
    Load(3);
    return BDDNodeCount();
  }

  void TSM(bool fCompl) {
    OptimizationStartup();
    for(int i = 1; i < nInputs; i++) {
      for(int index: vvIndices[i-1]) {
        int cof0index = index << 1;
        int cof1index = cof0index ^ 1;
        if(int r = Intersect(cof0index, cof1index, i, fCompl)) {
          CopyFuncMasked(cof0index, cof1index, i, !(r & 1));
          MergeCare(cof0index, cof1index, i);
          vvIndicesMerged[i].push_back({(cof0index << 1) ^ !(r & 1), cof1index});
          BDDCountNodesCareOne(cof0index, i);
        } else {
          BDDCountNodesCareOne(cof0index, i);
          BDDCountNodesCareOne(cof1index, i);
        }
      }
    }
    Merge();
  }

  void TSMNew(bool fCompl) {
    originalt = t;
    std::vector<std::vector<std::pair<int, int> > > merged(nInputs);
    vvIndices.clear();
    vvIndices.resize(nInputs);
    for(int i = 0; i < nOutputs; i++) {
      if(IsDC(i, 0)) {
        if(nSize) {
          for(int j = 0; j < nSize; j++) {
            t[j + nSize * i] = 0;
          }
        } else {
          SetValue(i, 0, 0);
        }
        continue;
      }
      int r = BDDCountNodesOne(i, 0);
      if(r != i << 1) {
        MergeCare(r >> 1, i, 0);
        merged[0].push_back({r, i});
      }
    }
    for(int i = 1; i < nInputs; i++) {
      std::vector<int> vChildren;
      for(int index: vvIndices[i-1]) {
        if(int r = Include(index << 1, (index << 1) ^ 1, i, fCompl)) {
          MergeCare(index << 1, (index << 1) ^ 1, i);
          int cof0 = BDDCountNodesOne(index << 1, i);
          if(cof0 != index << 2) {
            MergeCare(cof0 >> 1, index << 1, i);
            merged[i].push_back({cof0, index << 1});
          }
          merged[i].push_back({(index << 2) ^ 1 ^ (r & 1), (index << 1) ^ 1});
          vChildren.push_back(cof0);
          vChildren.push_back(cof0 ^ 1 ^ (r & 1));
        } else if(int r = Include((index << 1) ^ 1, index << 1, i, fCompl)) {
          MergeCare((index << 1) ^ 1, index << 1, i);
          int cof1 = BDDCountNodesOne((index << 1) ^ 1, i);
          if(cof1 != ((index << 2) ^ 2)) {
            MergeCare(cof1 >> 1, (index << 1) ^ 1, i);
            merged[i].push_back({cof1, (index << 1) ^ 1});
          }
          merged[i].push_back({(index << 2) ^ 2 ^ 1 ^ (r & 1), index << 1});
          vChildren.push_back(cof1 ^ 1 ^ (r & 1));
          vChildren.push_back(cof1);
        } else {
          int cof0 = BDDCountNodesOne(index << 1, i);
          if(cof0 != index << 2) {
            MergeCare(cof0 >> 1, index << 1, i);
            merged[i].push_back({cof0, index << 1});
          }
          int cof1 = BDDCountNodesOne((index << 1) ^ 1, i);
          if(cof1 != ((index << 2) ^ 2)) {
            MergeCare(cof1 >> 1, (index << 1) ^ 1, i);
            merged[i].push_back({cof1, (index << 1) ^ 1});
          }
          vChildren.push_back(cof0);
          vChildren.push_back(cof1);
        }
      }
      vvIndices[i].clear();
      std::map<int, int> m;
      for(uint j = 0; j < vvIndices[i-1].size(); j++) {
        int cof0old = vChildren[j+j];
        if(m.count(cof0old >> 1)) {
          cof0old = m[cof0old >> 1] ^ (cof0old & 1);
        }
        int cof1old = vChildren[j+j+1];
        if(m.count(cof1old >> 1)) {
          cof1old = m[cof1old >> 1] ^ (cof1old & 1);
        }
        int cof0index = cof0old >> 1;
        int cof1index = cof1old >> 1;
        bool cof0c = cof0old & 1;
        bool cof1c = cof1old & 1;
        if(cof0index < 0 || cof1index < 0) {
          if(cof0index >= 0) {
            int cof0 = BDDCountNodesOne(cof0index, i);
            if(cof0 != cof0index << 1) {
              MergeCare(cof0 >> 1, cof0index, i);
              merged[i].insert(merged[i].begin(), {cof0, cof0index});
            }
            m[cof0index] = cof0;
          }
          if(cof1index >= 0) {
            int cof1 = BDDCountNodesOne(cof1index, i);
            if(cof1 != cof1index << 1) {
              MergeCare(cof1 >> 1, cof1index, i);
              merged[i].insert(merged[i].begin(), {cof1, cof1index});
            }
            m[cof1index] = cof1;
          }
          continue;
        }
        if(cof0index == cof1index) {
          int cof0 = BDDCountNodesOne(cof0index, i);
          if(cof0 != cof0index << 1) {
            MergeCare(cof0 >> 1, cof0index, i);
            merged[i].insert(merged[i].begin(), {cof0, cof0index});
          }
          m[cof0index] = cof0;
          continue;
        }
        bool fComplCof = cof0c ^ cof1c;
        if(int r = Intersect(cof0index, cof1index, i, fComplCof || fCompl, !fComplCof || fCompl)) {
          CopyFuncMasked(cof0index, cof1index, i, !(r & 1));
          MergeCare(cof0index, cof1index, i);
          merged[i].insert(merged[i].begin(), {(cof0index << 1) ^ 1 ^ (r & 1), cof1index});
          int cof0 = BDDCountNodesOne(cof0index, i);
          if(cof0 != cof0index << 1) {
            MergeCare(cof0 >> 1, cof0index, i);
            merged[i].insert(merged[i].begin(), {cof0, cof0index});
          }
          m[cof0index] = cof0;
        } else {
          int cof0 = BDDCountNodesOne(cof0index, i);
          if(cof0 != cof0index << 1) {
            MergeCare(cof0 >> 1, cof0index, i);
            merged[i].insert(merged[i].begin(), {cof0, cof0index});
          }
          m[cof0index] = cof0;
          int cof1 = BDDCountNodesOne(cof1index, i);
          if(cof1 != cof1index << 1) {
            MergeCare(cof1 >> 1, cof1index, i);
            merged[i].insert(merged[i].begin(), {cof1, cof1index});
          }
          m[cof1index] = cof1;
        }
      }
    }
    for(int i = nInputs - 1; i >= 0; i--) {
      for(auto p: merged[i]) {
        CopyFunc(p.second, p.first >> 1, i, p.first & 1);
      }
    }
  }
};

class TTOSDM : public TTCare{
public:
  TTOSDM(std::vector<std::vector<int> > const &onsets, int nInputs, std::vector<char *> const &pBPats, int nBPats, int rarity): TTCare(onsets, nInputs, pBPats, nBPats, rarity) {}

  int BDDCountNodes() override {
    return BDDCountNodesOSDM();
  }
};

class TTOSM : public TTCare{
public:
  TTOSM(std::vector<std::vector<int> > const &onsets, int nInputs, std::vector<char *> const &pBPats, int nBPats, int rarity): TTCare(onsets, nInputs, pBPats, nBPats, rarity) {}

  int BDDCountNodes() override {
    return BDDCountNodesOSM(false);
  }
};

class TTOSMC : public TTCare{
public:
  TTOSMC(std::vector<std::vector<int> > const &onsets, int nInputs, std::vector<char *> const &pBPats, int nBPats, int rarity): TTCare(onsets, nInputs, pBPats, nBPats, rarity) {}

  int BDDCountNodes() override {
    return BDDCountNodesOSM(true);
  }
};

class TTTSM : public TTCare{
public:
  TTTSM(std::vector<std::vector<int> > const &onsets, int nInputs, std::vector<char *> const &pBPats, int nBPats, int rarity): TTCare(onsets, nInputs, pBPats, nBPats, rarity) {}

  int BDDCountNodes() override {
    return BDDCountNodesTSM(true);
  }
};

void TTTest(std::vector<std::vector<int> > const &onsets, std::vector<char *> const &pBPats, int nBPats, int rarity, std::vector<std::string> const &inputs, std::vector<std::string> const &outputs, std::ofstream &f) {
  int nInputs = inputs.size();
  // TT tt(onsets, nInputs);
  // std::cout << tt.BDDCountNodes() << std::endl;
  // // tt.SiftReo();
  // tt.RandomSiftReo(20);
  // std::cout << tt.BDDCountNodes() << std::endl;
  // tt.GeneratePla("test.pla");

  // TTCare tt(onsets, nInputs, pBPats, nBPats, rarity);
  // std::cout << tt.BDDCountNodes() << std::endl;
  // tt.RandomSiftReo(20);
  // std::cout << tt.BDDCountNodes() << std::endl;
  // tt.OSM(true);
  // std::cout << tt.BDDCountNodes() << std::endl;
  // std::cout << tt.BDDCountNodesOSM(false) << std::endl;
  // std::cout << tt.BDDCountNodesOSM(true) << std::endl;
  // tt.TSM(true);
  // std::cout << tt.BDDCountNodes() << std::endl;
  // tt.GeneratePlaCare("care.pla");
  // tt.GeneratePlaMasked("test.pla");
  // tt.BDDGenerateBlif(inputs, outputs, f);

  TTTSM tt(onsets, nInputs, pBPats, nBPats, rarity);
  tt.RandomSiftReo(20);
  tt.TSM(true);
  tt.BDDGenerateBlif(inputs, outputs, f);

  // TTCare tt(onsets, nInputs, pBPats, nBPats, rarity);
  // tt.RandomSiftReo(20);
  // int r1 = tt.BDDCountNodesOSM(false);
  // int r2 = tt.BDDCountNodesOSM(true);
  // tt.OSM(r2 <= r1);
  // tt.BDDGenerateBlif(inputs, outputs, f);
  
  // TTOSM tt1(onsets, nInputs, pBPats, nBPats, rarity);
  // int r1 = tt1.RandomSiftReo(20);
  // TTOSMC tt2(onsets, nInputs, pBPats, nBPats, rarity);
  // int r2 = tt2.RandomSiftReo(20);
  // if(r1 < r2){
  //   tt1.OSM(false);
  //   tt1.BDDGenerateBlif(inputs, outputs, f);
  // } else {
  //   tt2.OSM(true);
  //   tt2.BDDGenerateBlif(inputs, outputs, f);
  // }
}
