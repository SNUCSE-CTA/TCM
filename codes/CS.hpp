#ifndef _CS_H
#define _CS_H

#include <algorithm>
#include <cstring>
#include <limits>
#include <queue>
#include <vector>

#include "../checktime.hpp"
#include "graph.hpp"

using namespace std;

long long numMatches = 0;
long long numLocalMatched = 1;
long long numBacktrack = 0;
long long numFiltered = 0;

bool debugCheck = false;
struct CSNode {
  int numCandidates;
  // int* candidates; // use G.verticesByLabel
  int* mark;
  int* parentCount;
  int* childCount;
  int* validParentCand;
  int* validChildCand;
  int* validNeighborCand;
  bool* NLFCheck;

  int** E1min;  // E1min[e, v]
  int** E1max;
  int** E2min;
  int** E2max;
  int** E1minEdge;
  int** E1maxEdge;
  int** E2minEdge;
  int** E2maxEdge;
  bool** E1minLazy;
  bool** E1maxLazy;
  bool** E2minLazy;
  bool** E2maxLazy;
};

struct CSEdge {
  int qeID, deID;

  CSEdge(const int _qeID = -1, const int _deID = -1)
      : qeID(_qeID), deID(_deID) {}

  CSEdge(const CSEdge& other) : qeID(other.qeID), deID(other.deID) {}
};

class CSStructure {
 public:
  DataGraph G;
  QueryGraph Q;

  int numNodes;
  CSNode* CSNodes;

  int rootVertex;

  // priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int,
  // int>>> Qmin; priority_queue<pair<int, int>> Qmax;
  queue<pair<int, int>> Qmin, Qmax;
  priority_queue<pair<int, int>, vector<pair<int, int>>,
                 greater<pair<int, int>>>
      QEmin;
  priority_queue<pair<int, int>> QEmax;
  vector<vector<bool>> QEminVisited, QEmaxVisited;
  vector<vector<int>> QEminVertices, QEmaxVertices;

  void allocate();
  void init();

  vector<vector<bool>> E1CheckPrev, E2CheckPrev;
  vector<vector<bool>> E1CheckNow, E2CheckNow;

  bool computeNLF(int u, int v);
  bool isTCMatchable(int qeID, int deID);
  bool isTCMatchableE1(int qeID, int deID);
  bool isTCMatchableE2(int qeID, int deID);

  // TODO: inline function?
  vector<CSEdge> findCSChangedEdges(const Edge& e);
  vector<CSEdge> findNLFChangedEdges(const Edge& e);

  vector<CSEdge> insertionUpdateE(vector<CSEdge>& CSEdges);
  void insertionUpdateETopDown(int qeID, int deID, vector<CSEdge>& TCMCSEdges);
  void insertionUpdateEBottomUp(int qeID, int deID, vector<CSEdge>& TCMCSEdges);

  void insertionUpdateD(vector<CSEdge>& CSEdges, vector<CSEdge>& NLFEdges);
  void insertionUpdateDTopDown(int u, int v, int uc, int vc);
  void insertionUpdateDBottomUp(int up, int vp, int u, int v);

  vector<CSEdge> deletionUpdateE(vector<CSEdge>& CSEdges);
  void deletionUpdateETopDown(int qeID, int deID, vector<CSEdge>& TCMCSEdges);
  void deletionUpdateEBottomUp(int qeID, int deID, vector<CSEdge>& TCMCSEdges);

  void deletionUpdateD(vector<CSEdge>& CSEdges, vector<CSEdge>& NLFEdges);
  void deletionUpdateDTopDown(int u, int v, int uc, int vc);
  void deletionUpdateDBottomUp(int up, int vp, int u, int v);

  // vector<vector<pair<int, int>>> Cm; // (Vertex, Time)
  int** CmVertices;
  int** CmVerticesDuplicateCount;
  int** EuSizes;
  int** CmMinNeighbor;
  int* CmSizes;
  vector<vector<int>> ECm;
  bool* isProduct;

  bool* visitedQuery;
  bool* visitedData;
  bool* checkData;
  int* extendibleQuery;  // number of matched neighbor vertices
  int* vertexMatch;
  int* edgeMatch;
  int* minEdge;

  bool* isIsolated;
  int* isolatedVertices;  // [?] pair<pair<int, int>, int>* or compare function
  long long* matchedCountDP;
  long long* storedMatches;

  uint64_t** edgeFailingSet;
  uint64_t** curEdgeFailingSet;
  int depth;
  int numMatchedEdges;

  vector<int> computeCm(int u, int lvl);
  void computeECm(int qeID, bool& existBefore, bool& existAfter);
  void updateVertexMatch(int u, int v, int lvl);
  void restoreVertexMatch(int u, int v, int lvl);
  void updateEdgeMatch(int qeID, int deID, bool dummy);
  void restoreEdgeMatch(int qeID, int deID, bool dummy);
  // void backtrackIsolated(int s, int e, int idx, int cnt, int start);
  void computeNumMatches(int qeID);
  bool backtrack(int lvl, vector<int>& checkEdges, int idx);
  void findMatches(vector<CSEdge>& CSEdges);

  void printCS() {
    for (int i = 0; i < numNodes; i++) {
      cout << i << " node" << endl;
      cout << "mark: ";
      for (int j = 0; j < CSNodes[i].numCandidates; j++) {
        cout << CSNodes[i].mark[j] << " ";
      }
      cout << endl;
    }
  }
};

bool CSStructure::computeNLF(int u, int v) {
  for (int i = 0; i < NLFArraySize; i++) {
    if ((Q.queryNLFBit[u * NLFArraySize + i] &
         G.dataNLFBit[v * NLFArraySize + i]) !=
        Q.queryNLFBit[u * NLFArraySize + i]) {
      return false;
    }
  }
  return true;
}

bool CSStructure::isTCMatchable(int qeID, int deID) {
  auto& qe = Q.edges[qeID];
  auto& de = G.stream[deID];
  int u = qe.u, up = qe.v;
  int v = de.u, vp = de.v;
  if (Q.dagOrderInv[u] < Q.dagOrderInv[up]) {
    swap(u, up);
    swap(v, vp);
  }
  int vIndex = G.vertexIDByLabel[v], vpIndex = G.vertexIDByLabel[vp];

  for (auto beID : Q.beforeEdges[qeID]) {
    for (auto loc : Q.minPathTreeRevEdgeLoc[up][beID]) {
      // if(CSNodes[up].E1min[loc] == NULL) continue;
      if (CSNodes[up].E1min[loc][vpIndex] == -1) return false;
      if (de.time < G.stream[CSNodes[up].E1min[loc][vpIndex]].time)
        return false;
    }

    for (auto loc : Q.minPathTreeEdgeLoc[u][beID]) {
      // if(CSNodes[u].E2min[loc] == NULL) continue;
      if (CSNodes[u].E2min[loc][vIndex] == -1) return false;
      if (de.time < G.stream[CSNodes[u].E2min[loc][vIndex]].time) return false;
    }
  }

  for (auto aeID : Q.afterEdges[qeID]) {
    for (auto loc : Q.maxPathTreeRevEdgeLoc[up][aeID]) {
      // if(CSNodes[up].E1max[loc] == NULL) continue;
      if (CSNodes[up].E1max[loc][vpIndex] == -1) return false;
      if (de.time > G.stream[CSNodes[up].E1max[loc][vpIndex]].time)
        return false;
    }

    for (auto loc : Q.maxPathTreeEdgeLoc[u][aeID]) {
      // if(CSNodes[u].E2max[loc] == NULL) continue;
      if (CSNodes[u].E2max[loc][vIndex] == -1) return false;
      if (de.time > G.stream[CSNodes[u].E2max[loc][vIndex]].time) return false;
    }
  }

  return true;
}

bool CSStructure::isTCMatchableE1(int qeID, int deID) {
  auto& qe = Q.edges[qeID];
  auto& de = G.stream[deID];
  int u = qe.u, up = qe.v;
  int v = de.u, vp = de.v;
  if (Q.dagOrderInv[u] < Q.dagOrderInv[up]) {
    swap(u, up);
    swap(v, vp);
  }
  int vIndex = G.vertexIDByLabel[v], vpIndex = G.vertexIDByLabel[vp];

  for (auto beID : Q.beforeEdges[qeID]) {
    for (auto loc : Q.minPathTreeRevEdgeLoc[up][beID]) {
      // if(CSNodes[up].E1min[loc] == NULL) continue;
      if (CSNodes[up].E1min[loc][vpIndex] == -1) return false;
      if (de.time < G.stream[CSNodes[up].E1min[loc][vpIndex]].time)
        return false;
    }
  }

  for (auto aeID : Q.afterEdges[qeID]) {
    for (auto loc : Q.maxPathTreeRevEdgeLoc[up][aeID]) {
      // if(CSNodes[up].E1max[loc] == NULL) continue;
      if (CSNodes[up].E1max[loc][vpIndex] == -1) return false;
      if (de.time > G.stream[CSNodes[up].E1max[loc][vpIndex]].time)
        return false;
    }
  }

  return true;
}

bool CSStructure::isTCMatchableE2(int qeID, int deID) {
  auto& qe = Q.edges[qeID];
  auto& de = G.stream[deID];
  int u = qe.u, up = qe.v;
  int v = de.u, vp = de.v;
  if (Q.dagOrderInv[u] < Q.dagOrderInv[up]) {
    swap(u, up);
    swap(v, vp);
  }
  int vIndex = G.vertexIDByLabel[v], vpIndex = G.vertexIDByLabel[vp];

  for (auto beID : Q.beforeEdges[qeID]) {
    for (auto loc : Q.minPathTreeEdgeLoc[u][beID]) {
      // if(CSNodes[u].E2min[loc] == NULL) continue;
      if (CSNodes[u].E2min[loc][vIndex] == -1) return false;
      if (de.time < G.stream[CSNodes[u].E2min[loc][vIndex]].time) return false;
    }
  }

  for (auto aeID : Q.afterEdges[qeID]) {
    for (auto loc : Q.maxPathTreeEdgeLoc[u][aeID]) {
      // if(CSNodes[u].E2max[loc] == NULL) continue;
      if (CSNodes[u].E2max[loc][vIndex] == -1) return false;
      if (de.time > G.stream[CSNodes[u].E2max[loc][vIndex]].time) return false;
    }
  }

  return true;
}

void CSStructure::allocate() {
  // NLF
  NLFBitSize = NLFHash.size() * 4;
  NLFArraySize = (NLFBitSize + 63) >> 6;

  Q.queryNLF = new int[Q.numVertices * NLFHash.size()]();
  Q.queryNLFBit = new uint64_t[Q.numVertices * NLFArraySize]();

  G.dataNLF = new int[G.numVertices * NLFHash.size()]();
  G.dataNLFBit = new uint64_t[G.numVertices * NLFArraySize]();

  // Query NEC
  // Q.isLeafNode = new bool[Q.numVertices]();
  // Q.NECMapping = new int[Q.numVertices];
  // Q.NECInverse = new int[Q.numVertices];
  // Q.NECSize = new int[Q.numVertices];

  // CS
  numNodes = Q.numVertices;
  CSNodes = new CSNode[numNodes];
  QEminVisited.resize(numNodes);
  QEmaxVisited.resize(numNodes);
  QEminVertices.resize(numNodes);
  QEmaxVertices.resize(numNodes);

  for (int i = 0; i < numNodes; i++) {
    CSNodes[i].numCandidates = G.numVerticesByLabel[Q.vLabels[i]];
    CSNodes[i].mark = new int[CSNodes[i].numCandidates]();
    CSNodes[i].parentCount = new int[CSNodes[i].numCandidates]();
    CSNodes[i].childCount = new int[CSNodes[i].numCandidates]();
    CSNodes[i].validParentCand =
        new int[CSNodes[i].numCandidates * Q.dagParent[i].size()]();
    CSNodes[i].validChildCand =
        new int[CSNodes[i].numCandidates * Q.dagChild[i].size()]();
    CSNodes[i].validNeighborCand =
        new int[CSNodes[i].numCandidates * Q.neighbors[i].size()]();
    CSNodes[i].NLFCheck = new bool[CSNodes[i].numCandidates]();

    CSNodes[i].E1min = new int*[Q.minPathTreeRev[i].size()];
    CSNodes[i].E1minEdge = new int*[Q.minPathTreeRev[i].size()];
    CSNodes[i].E1minLazy = new bool*[Q.minPathTreeRev[i].size()];
    for (int j = 0; j < Q.minPathTreeRev[i].size(); j++) {
      CSNodes[i].E1min[j] = new int[CSNodes[i].numCandidates];
      CSNodes[i].E1minEdge[j] = new int[CSNodes[i].numCandidates];
      CSNodes[i].E1minLazy[j] = new bool[CSNodes[i].numCandidates]();
      memset(CSNodes[i].E1min[j], -1, sizeof(int) * CSNodes[i].numCandidates);
      memset(CSNodes[i].E1minEdge[j], -1,
             sizeof(int) * CSNodes[i].numCandidates);
    }
    CSNodes[i].E1max = new int*[Q.maxPathTreeRev[i].size()];
    CSNodes[i].E1maxEdge = new int*[Q.maxPathTreeRev[i].size()];
    CSNodes[i].E1maxLazy = new bool*[Q.maxPathTreeRev[i].size()];
    for (int j = 0; j < Q.maxPathTreeRev[i].size(); j++) {
      CSNodes[i].E1max[j] = new int[CSNodes[i].numCandidates];
      CSNodes[i].E1maxEdge[j] = new int[CSNodes[i].numCandidates];
      CSNodes[i].E1maxLazy[j] = new bool[CSNodes[i].numCandidates]();
      memset(CSNodes[i].E1max[j], -1, sizeof(int) * CSNodes[i].numCandidates);
      memset(CSNodes[i].E1maxEdge[j], -1,
             sizeof(int) * CSNodes[i].numCandidates);
    }

    CSNodes[i].E2min = new int*[Q.minPathTree[i].size()];
    CSNodes[i].E2minEdge = new int*[Q.minPathTree[i].size()];
    CSNodes[i].E2minLazy = new bool*[Q.minPathTree[i].size()];
    for (int j = 0; j < Q.minPathTree[i].size(); j++) {
      CSNodes[i].E2min[j] = new int[CSNodes[i].numCandidates];
      CSNodes[i].E2minEdge[j] = new int[CSNodes[i].numCandidates];
      CSNodes[i].E2minLazy[j] = new bool[CSNodes[i].numCandidates]();
      memset(CSNodes[i].E2min[j], -1, sizeof(int) * CSNodes[i].numCandidates);
      memset(CSNodes[i].E2minEdge[j], -1,
             sizeof(int) * CSNodes[i].numCandidates);
    }
    CSNodes[i].E2max = new int*[Q.maxPathTree[i].size()];
    CSNodes[i].E2maxEdge = new int*[Q.maxPathTree[i].size()];
    CSNodes[i].E2maxLazy = new bool*[Q.maxPathTree[i].size()];
    for (int j = 0; j < Q.maxPathTree[i].size(); j++) {
      CSNodes[i].E2max[j] = new int[CSNodes[i].numCandidates];
      CSNodes[i].E2maxEdge[j] = new int[CSNodes[i].numCandidates];
      CSNodes[i].E2maxLazy[j] = new bool[CSNodes[i].numCandidates]();
      memset(CSNodes[i].E2max[j], -1, sizeof(int) * CSNodes[i].numCandidates);
      memset(CSNodes[i].E2maxEdge[j], -1,
             sizeof(int) * CSNodes[i].numCandidates);
    }

    QEminVisited[i].resize(CSNodes[i].numCandidates, false);
    QEmaxVisited[i].resize(CSNodes[i].numCandidates, false);
  }

  E1CheckPrev.resize(Q.numEdges);
  E2CheckPrev.resize(Q.numEdges);
  E1CheckNow.resize(Q.numEdges);
  E2CheckNow.resize(Q.numEdges);
  for (int i = 0; i < Q.numEdges; i++) {
    E1CheckPrev[i].resize(G.numEdges);
    E2CheckPrev[i].resize(G.numEdges);
    E1CheckNow[i].resize(G.numEdges);
    E2CheckNow[i].resize(G.numEdges);
  }

  // Backtrack
  visitedQuery = new bool[Q.numVertices]();
  visitedData = new bool[G.numVertices]();
  checkData = new bool[G.numVertices]();
  extendibleQuery = new int[Q.numVertices]();
  vertexMatch = new int[Q.numVertices]();
  memset(vertexMatch, -1, sizeof(int) * Q.numVertices);
  edgeMatch = new int[Q.numEdges]();
  memset(edgeMatch, -1, sizeof(int) * Q.numEdges);
  minEdge = new int[Q.numEdges + 1]();

  CmVertices = new int*[Q.numVertices];
  CmVerticesDuplicateCount = new int*[Q.numVertices];
  for (int i = 0; i < Q.numVertices; i++) {
    CmVertices[i] = new int[G.numVertices]();
    CmVerticesDuplicateCount[i] = new int[G.numVertices]();
  }
  EuSizes = new int*[Q.numVertices];
  CmMinNeighbor = new int*[Q.numVertices];
  for (int i = 0; i < Q.numVertices; i++) {
    EuSizes[i] = new int[Q.numVertices]();
    CmMinNeighbor[i] = new int[Q.numVertices]();
  }
  CmSizes = new int[Q.numVertices]();
  ECm.resize(Q.numEdges);
  isProduct = new bool[Q.numEdges]();

  isIsolated = new bool[Q.numVertices]();
  isolatedVertices = new int[Q.numVertices]();
  matchedCountDP = new long long[Q.numVertices + 1]();
  storedMatches = new long long[G.numEdges]();

  edgeFailingSet = new uint64_t*[Q.numEdges + 1]();
  for (int i = 0; i < Q.numEdges + 1; i++) {
    edgeFailingSet[i] = new uint64_t[(Q.numEdges + 63) >> 6]();
  }
  curEdgeFailingSet = new uint64_t*[Q.numVertices + 1]();
  for (int i = 0; i < Q.numVertices + 1; i++) {
    curEdgeFailingSet[i] = new uint64_t[(Q.numEdges + 63) >> 6]();
  }
  numMatchedEdges = 0;
  depth = 0;
}

void CSStructure::init() {}

vector<CSEdge> CSStructure::findCSChangedEdges(const Edge& de) {
  vector<CSEdge> CSEdges;
  for (auto qe : Q.edges) {
    if (de.isMatch(qe)) {
      CSEdges.emplace_back(qe.id, de.id);
    }
  }

  return CSEdges;
}

vector<CSEdge> CSStructure::findNLFChangedEdges(const Edge& de) {
  vector<CSEdge> NLFEdges;
  for (auto qe : Q.edges) {
    if (de.isNLFMatch(qe)) {
      NLFEdges.emplace_back(qe.id, de.id);
    }
  }

  return NLFEdges;
}

vector<CSEdge> CSStructure::insertionUpdateE(vector<CSEdge>& CSEdges) {
  vector<CSEdge> TCMCSEdges;

  for (auto CSEdge : CSEdges) {
    if (!E1CheckNow[CSEdge.qeID][CSEdge.deID])
      E1CheckNow[CSEdge.qeID][CSEdge.deID] =
          isTCMatchableE1(CSEdge.qeID, CSEdge.deID);
    if (E1CheckNow[CSEdge.qeID][CSEdge.deID])
      insertionUpdateETopDown(CSEdge.qeID, CSEdge.deID, TCMCSEdges);
    if (!E2CheckNow[CSEdge.qeID][CSEdge.deID])
      E2CheckNow[CSEdge.qeID][CSEdge.deID] =
          isTCMatchableE2(CSEdge.qeID, CSEdge.deID);
    if (E2CheckNow[CSEdge.qeID][CSEdge.deID])
      insertionUpdateEBottomUp(CSEdge.qeID, CSEdge.deID, TCMCSEdges);
  }

  while (!QEmin.empty()) {
    int u = Q.dagOrder[QEmin.top().first];
    int v = QEmin.top().second;
    int vIndex = G.vertexIDByLabel[v];
    QEmin.pop();
    QEminVisited[u][vIndex] = false;

    if (!QEminVertices[u].empty()) {
      for (auto v2 : QEminVertices[u]) {
        for (auto qe : Q.dagChildInEdges[u])
          for (auto de : G.inEdges[v2][Q.edgeHash[qe.id]])
            if (de.isMatch(qe) && !E1CheckNow[qe.id][de.id])
              E1CheckNow[qe.id][de.id] = isTCMatchableE1(qe.id, de.id);

        for (auto qe : Q.dagChildOutEdges[u])
          for (auto de : G.outEdges[v2][Q.edgeHash[qe.id]])
            if (de.isMatch(qe) && !E1CheckNow[qe.id][de.id])
              E1CheckNow[qe.id][de.id] = isTCMatchableE1(qe.id, de.id);
      }
      QEminVertices[u].clear();
    }

    for (auto qe : Q.dagChildInEdges[u])
      for (auto de : G.inEdges[v][Q.edgeHash[qe.id]])
        if (de.isMatch(qe) && E1CheckNow[qe.id][de.id])
          insertionUpdateETopDown(qe.id, de.id, TCMCSEdges);

    for (auto qe : Q.dagChildOutEdges[u])
      for (auto de : G.outEdges[v][Q.edgeHash[qe.id]])
        if (de.isMatch(qe) && E1CheckNow[qe.id][de.id])
          insertionUpdateETopDown(qe.id, de.id, TCMCSEdges);
  }

  while (!QEmax.empty()) {
    int u = Q.dagOrder[QEmax.top().first];
    int v = QEmax.top().second;
    int vIndex = G.vertexIDByLabel[v];
    QEmax.pop();
    QEmaxVisited[u][vIndex] = false;

    if (!QEmaxVertices[u].empty()) {
      for (auto v2 : QEmaxVertices[u]) {
        for (auto qe : Q.dagParentInEdges[u])
          for (auto de : G.inEdges[v2][Q.edgeHash[qe.id]])
            if (de.isMatch(qe) && !E2CheckNow[qe.id][de.id])
              E2CheckNow[qe.id][de.id] = isTCMatchableE2(qe.id, de.id);

        for (auto qe : Q.dagParentOutEdges[u])
          for (auto de : G.outEdges[v2][Q.edgeHash[qe.id]])
            if (de.isMatch(qe) && !E2CheckNow[qe.id][de.id])
              E2CheckNow[qe.id][de.id] = isTCMatchableE2(qe.id, de.id);
      }
      QEmaxVertices[u].clear();
    }

    for (auto qe : Q.dagParentInEdges[u])
      for (auto de : G.inEdges[v][Q.edgeHash[qe.id]])
        if (de.isMatch(qe) && E2CheckNow[qe.id][de.id])
          insertionUpdateEBottomUp(qe.id, de.id, TCMCSEdges);

    for (auto qe : Q.dagParentOutEdges[u])
      for (auto de : G.outEdges[v][Q.edgeHash[qe.id]])
        if (de.isMatch(qe) && E2CheckNow[qe.id][de.id])
          insertionUpdateEBottomUp(qe.id, de.id, TCMCSEdges);
  }

  return TCMCSEdges;
}

void CSStructure::insertionUpdateETopDown(int qeID, int deID,
                                          vector<CSEdge>& TCMCSEdges) {
  Edge& qe = Q.edges[qeID];
  Edge& de = G.stream[deID];

  int u1 = qe.u, u2 = qe.v, v1 = de.u, v2 = de.v;
  if (Q.dagOrderInv[u1] > Q.dagOrderInv[u2]) {
    swap(u1, u2);
    swap(v1, v2);
  }

  int v1Index = G.vertexIDByLabel[v1];
  int v2Index = G.vertexIDByLabel[v2];

  int u1Index = Q.dagIndex[u2][u1];
  bool updated = false;
  int start, end;

  // E1min update
  start = Q.minPathTreeRevIndex[u2][u1Index];
  end = Q.minPathTreeRevIndex[u2][u1Index + 1];
  for (int u2Pos = start; u2Pos < end; u2Pos++) {
    int u1Pos = Q.minPathTreeRevPos[u2][u2Pos];
    if (u1Pos != -1) {
      if (CSNodes[u1].E1min[u1Pos][v1Index] == -1) continue;
      if (CSNodes[u2].E1min[u2Pos][v2Index] == -1 ||
          G.stream[CSNodes[u2].E1min[u2Pos][v2Index]].time >
              G.stream[CSNodes[u1].E1min[u1Pos][v1Index]].time) {
        CSNodes[u2].E1min[u2Pos][v2Index] = CSNodes[u1].E1min[u1Pos][v1Index];
        CSNodes[u2].E1minEdge[u2Pos][v2Index] = de.id;
        updated = true;
      }
    } else {
      if (CSNodes[u2].E1min[u2Pos][v2Index] == -1 ||
          G.stream[CSNodes[u2].E1min[u2Pos][v2Index]].time >
              G.stream[de.id].time) {
        CSNodes[u2].E1min[u2Pos][v2Index] = de.id;
        CSNodes[u2].E1minEdge[u2Pos][v2Index] = de.id;
        updated = true;
      }
    }
  }

  // E1max update
  start = Q.maxPathTreeRevIndex[u2][u1Index];
  end = Q.maxPathTreeRevIndex[u2][u1Index + 1];
  for (int u2Pos = start; u2Pos < end; u2Pos++) {
    int u1Pos = Q.maxPathTreeRevPos[u2][u2Pos];
    if (u1Pos != -1) {
      if (CSNodes[u1].E1max[u1Pos][v1Index] == -1) continue;
      if (CSNodes[u2].E1max[u2Pos][v2Index] == -1 ||
          G.stream[CSNodes[u2].E1max[u2Pos][v2Index]].time <
              G.stream[CSNodes[u1].E1max[u1Pos][v1Index]].time) {
        CSNodes[u2].E1max[u2Pos][v2Index] = CSNodes[u1].E1max[u1Pos][v1Index];
        CSNodes[u2].E1maxEdge[u2Pos][v2Index] = de.id;
        updated = true;
      }
    } else {
      if (CSNodes[u2].E1max[u2Pos][v2Index] == -1 ||
          G.stream[CSNodes[u2].E1max[u2Pos][v2Index]].time <
              G.stream[de.id].time) {
        CSNodes[u2].E1max[u2Pos][v2Index] = de.id;
        CSNodes[u2].E1maxEdge[u2Pos][v2Index] = de.id;
        updated = true;
      }
    }
  }

  if (updated && !QEminVisited[u2][v2Index]) {
    QEmin.emplace(Q.dagOrderInv[u2], v2);
    QEminVisited[u2][v2Index] = true;
    QEminVertices[u2].push_back(v2);
  }

  if (!E1CheckPrev[qeID][deID] && E2CheckPrev[qeID][deID])
    TCMCSEdges.emplace_back(qeID, deID);
  E1CheckPrev[qeID][deID] = true;
}

void CSStructure::insertionUpdateEBottomUp(int qeID, int deID,
                                           vector<CSEdge>& TCMCSEdges) {
  Edge& qe = Q.edges[qeID];
  Edge& de = G.stream[deID];

  int u1 = qe.u, u2 = qe.v, v1 = de.u, v2 = de.v;
  if (Q.dagOrderInv[u1] < Q.dagOrderInv[u2]) {
    swap(u1, u2);
    swap(v1, v2);
  }

  int v1Index = G.vertexIDByLabel[v1];
  int v2Index = G.vertexIDByLabel[v2];

  int u1Index = Q.dagIndex[u2][u1];
  bool updated = false;
  int start, end;

  // E2min update
  start = Q.minPathTreeIndex[u2][u1Index];
  end = Q.minPathTreeIndex[u2][u1Index + 1];
  for (int u2Pos = start; u2Pos < end; u2Pos++) {
    int u1Pos = Q.minPathTreePos[u2][u2Pos];
    if (u1Pos != -1) {
      if (CSNodes[u1].E2min[u1Pos][v1Index] == -1) continue;
      if (CSNodes[u2].E2min[u2Pos][v2Index] == -1 ||
          G.stream[CSNodes[u2].E2min[u2Pos][v2Index]].time >
              G.stream[CSNodes[u1].E2min[u1Pos][v1Index]].time) {
        CSNodes[u2].E2min[u2Pos][v2Index] = CSNodes[u1].E2min[u1Pos][v1Index];
        CSNodes[u2].E2minEdge[u2Pos][v2Index] = de.id;
        updated = true;
      }
    } else {
      if (CSNodes[u2].E2min[u2Pos][v2Index] == -1 ||
          G.stream[CSNodes[u2].E2min[u2Pos][v2Index]].time >
              G.stream[de.id].time) {
        CSNodes[u2].E2min[u2Pos][v2Index] = de.id;
        CSNodes[u2].E2minEdge[u2Pos][v2Index] = de.id;
        updated = true;
      }
    }
  }

  // E2max update
  start = Q.maxPathTreeIndex[u2][u1Index];
  end = Q.maxPathTreeIndex[u2][u1Index + 1];
  for (int u2Pos = start; u2Pos < end; u2Pos++) {
    int u1Pos = Q.maxPathTreePos[u2][u2Pos];
    if (u1Pos != -1) {
      if (CSNodes[u1].E2max[u1Pos][v1Index] == -1) continue;
      if (CSNodes[u2].E2max[u2Pos][v2Index] == -1 ||
          G.stream[CSNodes[u2].E2max[u2Pos][v2Index]].time <
              G.stream[CSNodes[u1].E2max[u1Pos][v1Index]].time) {
        CSNodes[u2].E2max[u2Pos][v2Index] = CSNodes[u1].E2max[u1Pos][v1Index];
        CSNodes[u2].E2maxEdge[u2Pos][v2Index] = de.id;
        updated = true;
      }
    } else {
      if (CSNodes[u2].E2max[u2Pos][v2Index] == -1 ||
          G.stream[CSNodes[u2].E2max[u2Pos][v2Index]].time <
              G.stream[de.id].time) {
        CSNodes[u2].E2max[u2Pos][v2Index] = de.id;
        CSNodes[u2].E2maxEdge[u2Pos][v2Index] = de.id;
        updated = true;
      }
    }
  }

  if (updated && !QEmaxVisited[u2][v2Index]) {
    QEmax.emplace(Q.dagOrderInv[u2], v2);
    QEmaxVisited[u2][v2Index] = true;
    QEmaxVertices[u2].push_back(v2);
  }

  if (E1CheckPrev[qeID][deID] && !E2CheckPrev[qeID][deID])
    TCMCSEdges.emplace_back(qeID, deID);
  E2CheckPrev[qeID][deID] = true;
}

vector<CSEdge> CSStructure::deletionUpdateE(vector<CSEdge>& CSEdges) {
  vector<CSEdge> TCMCSEdges;

  for (auto CSEdge : CSEdges) {
    E1CheckNow[CSEdge.qeID][CSEdge.deID] =
        E2CheckNow[CSEdge.qeID][CSEdge.deID] = false;
    if (E1CheckPrev[CSEdge.qeID][CSEdge.deID])
      deletionUpdateETopDown(CSEdge.qeID, CSEdge.deID, TCMCSEdges);
    if (E2CheckPrev[CSEdge.qeID][CSEdge.deID])
      deletionUpdateEBottomUp(CSEdge.qeID, CSEdge.deID, TCMCSEdges);
    if (E1CheckPrev[CSEdge.qeID][CSEdge.deID] &&
        E2CheckPrev[CSEdge.qeID][CSEdge.deID]) {
      TCMCSEdges.emplace_back(CSEdge.qeID, CSEdge.deID);
    }
    E1CheckPrev[CSEdge.qeID][CSEdge.deID] =
        E2CheckPrev[CSEdge.qeID][CSEdge.deID] = false;
  }

  while (!QEmin.empty()) {
    int u = Q.dagOrder[QEmin.top().first];
    int v = QEmin.top().second;
    int vIndex = G.vertexIDByLabel[v];
    QEmin.pop();
    QEminVisited[u][vIndex] = false;

    if (!QEminVertices[u].empty()) {
      for (auto v2 : QEminVertices[u]) {
        for (auto qe : Q.dagChildInEdges[u])
          for (auto de : G.inEdges[v2][Q.edgeHash[qe.id]])
            if (de.isMatch(qe) && E1CheckNow[qe.id][de.id])
              E1CheckNow[qe.id][de.id] = isTCMatchableE1(qe.id, de.id);

        for (auto qe : Q.dagChildOutEdges[u])
          for (auto de : G.outEdges[v2][Q.edgeHash[qe.id]])
            if (de.isMatch(qe) && E1CheckNow[qe.id][de.id])
              E1CheckNow[qe.id][de.id] = isTCMatchableE1(qe.id, de.id);
      }
      QEminVertices[u].clear();
    }

    for (auto qe : Q.dagChildInEdges[u])
      for (auto de : G.inEdges[v][Q.edgeHash[qe.id]])
        if (de.isMatch(qe) && E1CheckPrev[qe.id][de.id])
          deletionUpdateETopDown(qe.id, de.id, TCMCSEdges);

    for (auto qe : Q.dagChildOutEdges[u])
      for (auto de : G.outEdges[v][Q.edgeHash[qe.id]])
        if (de.isMatch(qe) && E1CheckPrev[qe.id][de.id])
          deletionUpdateETopDown(qe.id, de.id, TCMCSEdges);
  }

  while (!QEmax.empty()) {
    int u = Q.dagOrder[QEmax.top().first];
    int v = QEmax.top().second;
    int vIndex = G.vertexIDByLabel[v];
    QEmax.pop();
    QEmaxVisited[u][vIndex] = false;

    if (!QEmaxVertices[u].empty()) {
      for (auto v2 : QEmaxVertices[u]) {
        for (auto qe : Q.dagParentInEdges[u])
          for (auto de : G.inEdges[v2][Q.edgeHash[qe.id]])
            if (de.isMatch(qe) && E2CheckNow[qe.id][de.id])
              E2CheckNow[qe.id][de.id] = isTCMatchableE2(qe.id, de.id);

        for (auto qe : Q.dagParentOutEdges[u])
          for (auto de : G.outEdges[v2][Q.edgeHash[qe.id]])
            if (de.isMatch(qe) && E2CheckNow[qe.id][de.id])
              E2CheckNow[qe.id][de.id] = isTCMatchableE2(qe.id, de.id);
      }
      QEmaxVertices[u].clear();
    }

    for (auto qe : Q.dagParentInEdges[u])
      for (auto de : G.inEdges[v][Q.edgeHash[qe.id]])
        if (de.isMatch(qe) && E2CheckPrev[qe.id][de.id])
          deletionUpdateEBottomUp(qe.id, de.id, TCMCSEdges);

    for (auto qe : Q.dagParentOutEdges[u])
      for (auto de : G.outEdges[v][Q.edgeHash[qe.id]])
        if (de.isMatch(qe) && E2CheckPrev[qe.id][de.id])
          deletionUpdateEBottomUp(qe.id, de.id, TCMCSEdges);
  }

  return TCMCSEdges;
}

void CSStructure::deletionUpdateETopDown(int qeID, int deID,
                                         vector<CSEdge>& TCMCSEdges) {
  Edge& qe = Q.edges[qeID];
  Edge& de = G.stream[deID];

  int u1 = qe.u, u2 = qe.v, v1 = de.u, v2 = de.v;
  if (Q.dagOrderInv[u1] > Q.dagOrderInv[u2]) {
    swap(u1, u2);
    swap(v1, v2);
  }

  bool v2Out = (de.u == v2) ? true : false;
  auto& neighborEdges = (v2Out) ? G.outEdges[v2][Q.edgeHash[qe.id]]
                                : G.inEdges[v2][Q.edgeHash[qe.id]];

  int v1Index = G.vertexIDByLabel[v1];
  int v2Index = G.vertexIDByLabel[v2];

  int u1Index = Q.dagIndex[u2][u1];
  bool updated = false;
  int start, end;

  // E1min update
  start = Q.minPathTreeRevIndex[u2][u1Index];
  end = Q.minPathTreeRevIndex[u2][u1Index + 1];
  for (int u2Pos = start; u2Pos < end; u2Pos++) {
    int u1Pos = Q.minPathTreeRevPos[u2][u2Pos];
    if (u1Pos != -1) {
      if (CSNodes[u2].E1min[u2Pos][v2Index] == -1) continue;
      if (neighborEdges.size() == 0) {
        CSNodes[u2].E1min[u2Pos][v2Index] = -1;
        CSNodes[u2].E1minEdge[u2Pos][v2Index] = -1;
        updated = true;
        CSNodes[u2].E1minLazy[u2Pos][v2Index] = false;
        continue;
      }
      if (!CSNodes[u2].E1minLazy[u2Pos][v2Index]) {
        if (CSNodes[u2].E1minEdge[u2Pos][v2Index] != de.id) continue;
        if (CSNodes[u2].E1min[u2Pos][v2Index] ==
                CSNodes[u1].E1min[u1Pos][v1Index] &&
            E1CheckNow[qeID][deID])
          continue;
      }

      // if D_1[u2][v2] == 0 then pass
      if (CSNodes[u2].mark[v2Index] < 1) {
        CSNodes[u2].E1minLazy[u2Pos][v2Index] = true;
        continue;
      }
      CSNodes[u2].E1minLazy[u2Pos][v2Index] = false;

      int prevVal = CSNodes[u2].E1min[u2Pos][v2Index];
      CSNodes[u2].E1min[u2Pos][v2Index] = -1;
      CSNodes[u2].E1minEdge[u2Pos][v2Index] = -1;

      for (auto dne : neighborEdges) {
        if (!dne.isMatch(qe) || !E1CheckNow[qeID][dne.id]) continue;
        int vn = (v2Out) ? dne.v : dne.u;
        int vnIndex = G.vertexIDByLabel[vn];

        if (CSNodes[u1].E1min[u1Pos][vnIndex] == -1) continue;
        if (CSNodes[u2].E1min[u2Pos][v2Index] == -1 ||
            G.stream[CSNodes[u2].E1min[u2Pos][v2Index]].time >
                G.stream[CSNodes[u1].E1min[u1Pos][vnIndex]].time) {
          CSNodes[u2].E1min[u2Pos][v2Index] = CSNodes[u1].E1min[u1Pos][vnIndex];
          CSNodes[u2].E1minEdge[u2Pos][v2Index] = dne.id;
        }
      }
      if (prevVal != CSNodes[u2].E1min[u2Pos][v2Index]) updated = true;
    } else {
      if (neighborEdges.size() == 0) {
        CSNodes[u2].E1min[u2Pos][v2Index] = -1;
        CSNodes[u2].E1minEdge[u2Pos][v2Index] = -1;
        updated = true;
        CSNodes[u2].E1minLazy[u2Pos][v2Index] = false;
        continue;
      }
      if (!CSNodes[u2].E1minLazy[u2Pos][v2Index]) {
        if (CSNodes[u2].E1min[u2Pos][v2Index] != de.id) continue;
        if (E1CheckNow[qeID][deID]) continue;
      }

      // if D_1[u2][v2] == 0 then pass
      if (CSNodes[u2].mark[v2Index] < 1) {
        CSNodes[u2].E1minLazy[u2Pos][v2Index] = true;
        continue;
      }
      CSNodes[u2].E1minLazy[u2Pos][v2Index] = false;

      int prevVal = CSNodes[u2].E1min[u2Pos][v2Index];
      CSNodes[u2].E1min[u2Pos][v2Index] = -1;
      CSNodes[u2].E1minEdge[u2Pos][v2Index] = -1;

      for (auto dne : neighborEdges) {
        if (!dne.isMatch(qe) || !E1CheckNow[qeID][dne.id]) continue;
        if (CSNodes[u2].E1min[u2Pos][v2Index] == -1 ||
            G.stream[CSNodes[u2].E1min[u2Pos][v2Index]].time >
                G.stream[dne.id].time) {
          CSNodes[u2].E1min[u2Pos][v2Index] = dne.id;
          CSNodes[u2].E1minEdge[u2Pos][v2Index] = dne.id;
          break;
        }
      }
      if (prevVal != CSNodes[u2].E1min[u2Pos][v2Index]) updated = true;
    }
  }

  // E1max update
  start = Q.maxPathTreeRevIndex[u2][u1Index];
  end = Q.maxPathTreeRevIndex[u2][u1Index + 1];
  for (int u2Pos = start; u2Pos < end; u2Pos++) {
    int u1Pos = Q.maxPathTreeRevPos[u2][u2Pos];
    if (u1Pos != -1) {
      if (CSNodes[u2].E1max[u2Pos][v2Index] == -1) continue;
      if (neighborEdges.size() == 0) {
        CSNodes[u2].E1max[u2Pos][v2Index] = -1;
        CSNodes[u2].E1maxEdge[u2Pos][v2Index] = -1;
        updated = true;
        CSNodes[u2].E1maxLazy[u2Pos][v2Index] = false;
        continue;
      }
      if (!CSNodes[u2].E1maxLazy[u2Pos][v2Index]) {
        if (CSNodes[u2].E1maxEdge[u2Pos][v2Index] != de.id) continue;
        if (CSNodes[u2].E1max[u2Pos][v2Index] ==
                CSNodes[u1].E1max[u1Pos][v1Index] &&
            E1CheckNow[qeID][deID])
          continue;
      }

      // if D_1[u2][v2] == 0 then pass
      if (CSNodes[u2].mark[v2Index] < 1) {
        CSNodes[u2].E1maxLazy[u2Pos][v2Index] = true;
        continue;
      }
      CSNodes[u2].E1maxLazy[u2Pos][v2Index] = false;

      int prevVal = CSNodes[u2].E1max[u2Pos][v2Index];
      CSNodes[u2].E1max[u2Pos][v2Index] = -1;
      CSNodes[u2].E1maxEdge[u2Pos][v2Index] = -1;

      for (auto dne : neighborEdges) {
        if (!dne.isMatch(qe) || !E1CheckNow[qeID][dne.id]) continue;
        int vn = (v2Out) ? dne.v : dne.u;
        int vnIndex = G.vertexIDByLabel[vn];

        if (CSNodes[u1].E1max[u1Pos][vnIndex] == -1) continue;
        if (CSNodes[u2].E1max[u2Pos][v2Index] == -1 ||
            G.stream[CSNodes[u2].E1max[u2Pos][v2Index]].time <
                G.stream[CSNodes[u1].E1max[u1Pos][vnIndex]].time) {
          CSNodes[u2].E1max[u2Pos][v2Index] = CSNodes[u1].E1max[u1Pos][vnIndex];
          CSNodes[u2].E1maxEdge[u2Pos][v2Index] = dne.id;
        }
      }
      if (prevVal != CSNodes[u2].E1max[u2Pos][v2Index]) updated = true;
    } else {
      if (neighborEdges.size() == 0) {
        CSNodes[u2].E1max[u2Pos][v2Index] = -1;
        CSNodes[u2].E1maxEdge[u2Pos][v2Index] = -1;
        updated = true;
        CSNodes[u2].E1maxLazy[u2Pos][v2Index] = false;
        continue;
      }
      if (!CSNodes[u2].E1maxLazy[u2Pos][v2Index]) {
        if (CSNodes[u2].E1max[u2Pos][v2Index] != de.id) continue;
        if (E1CheckNow[qeID][deID]) continue;
      }

      // if D_1[u2][v2] == 0 then pass
      if (CSNodes[u2].mark[v2Index] < 1) {
        CSNodes[u2].E1maxLazy[u2Pos][v2Index] = true;
        continue;
      }
      CSNodes[u2].E1maxLazy[u2Pos][v2Index] = false;

      int prevVal = CSNodes[u2].E1max[u2Pos][v2Index];
      CSNodes[u2].E1max[u2Pos][v2Index] = -1;
      CSNodes[u2].E1maxEdge[u2Pos][v2Index] = -1;

      for (auto rit = neighborEdges.rbegin(); rit != neighborEdges.rend();
           rit++) {
        auto dne = *rit;
        if (!dne.isMatch(qe) || !E1CheckNow[qeID][dne.id]) continue;
        if (CSNodes[u2].E1max[u2Pos][v2Index] == -1 ||
            G.stream[CSNodes[u2].E1max[u2Pos][v2Index]].time <
                G.stream[dne.id].time) {
          CSNodes[u2].E1max[u2Pos][v2Index] = dne.id;
          CSNodes[u2].E1maxEdge[u2Pos][v2Index] = dne.id;
          break;
        }
      }
      if (prevVal != CSNodes[u2].E1max[u2Pos][v2Index]) updated = true;
    }
  }

  if (updated && !QEminVisited[u2][v2Index]) {
    QEmin.emplace(Q.dagOrderInv[u2], v2);
    QEminVisited[u2][v2Index] = true;
    QEminVertices[u2].push_back(v2);
  }

  E1CheckPrev[qeID][deID] = E1CheckNow[qeID][deID];
  if (!E1CheckPrev[qeID][deID] && E2CheckPrev[qeID][deID])
    TCMCSEdges.emplace_back(qeID, deID);
}

void CSStructure::deletionUpdateEBottomUp(int qeID, int deID,
                                          vector<CSEdge>& TCMCSEdges) {
  Edge& qe = Q.edges[qeID];
  Edge& de = G.stream[deID];

  int u1 = qe.u, u2 = qe.v, v1 = de.u, v2 = de.v;
  if (Q.dagOrderInv[u1] < Q.dagOrderInv[u2]) {
    swap(u1, u2);
    swap(v1, v2);
  }

  bool v2Out = (de.u == v2) ? true : false;
  auto& neighborEdges = (v2Out) ? G.outEdges[v2][Q.edgeHash[qe.id]]
                                : G.inEdges[v2][Q.edgeHash[qe.id]];

  int v1Index = G.vertexIDByLabel[v1];
  int v2Index = G.vertexIDByLabel[v2];

  int u1Index = Q.dagIndex[u2][u1];
  bool updated = false;
  int start, end;

  // E2min update
  start = Q.minPathTreeIndex[u2][u1Index];
  end = Q.minPathTreeIndex[u2][u1Index + 1];
  for (int u2Pos = start; u2Pos < end; u2Pos++) {
    int u1Pos = Q.minPathTreePos[u2][u2Pos];
    if (u1Pos != -1) {
      if (CSNodes[u2].E2min[u2Pos][v2Index] == -1) continue;
      if (neighborEdges.size() == 0) {
        CSNodes[u2].E2min[u2Pos][v2Index] = -1;
        CSNodes[u2].E2minEdge[u2Pos][v2Index] = -1;
        updated = true;
        CSNodes[u2].E2minLazy[u2Pos][v2Index] = false;
        continue;
      }
      if (!CSNodes[u2].E2minLazy[u2Pos][v2Index]) {
        if (CSNodes[u2].E2minEdge[u2Pos][v2Index] != de.id) continue;
        if (CSNodes[u2].E2min[u2Pos][v2Index] ==
                CSNodes[u1].E2min[u1Pos][v1Index] &&
            E2CheckNow[qeID][deID])
          continue;
      }

      // if D_2[u2][v2] == 0 then pass
      if (CSNodes[u2].mark[v2Index] < 2) {
        CSNodes[u2].E2minLazy[u2Pos][v2Index] = true;
        continue;
      }
      CSNodes[u2].E2minLazy[u2Pos][v2Index] = false;

      int prevVal = CSNodes[u2].E2min[u2Pos][v2Index];
      CSNodes[u2].E2min[u2Pos][v2Index] = -1;
      CSNodes[u2].E2minEdge[u2Pos][v2Index] = -1;

      for (auto dne : neighborEdges) {
        if (!dne.isMatch(qe) || !E2CheckNow[qeID][dne.id]) continue;
        int vn = (v2Out) ? dne.v : dne.u;
        int vnIndex = G.vertexIDByLabel[vn];

        if (CSNodes[u1].E2min[u1Pos][vnIndex] == -1) continue;
        if (CSNodes[u2].E2min[u2Pos][v2Index] == -1 ||
            G.stream[CSNodes[u2].E2min[u2Pos][v2Index]].time >
                G.stream[CSNodes[u1].E2min[u1Pos][vnIndex]].time) {
          CSNodes[u2].E2min[u2Pos][v2Index] = CSNodes[u1].E2min[u1Pos][vnIndex];
          CSNodes[u2].E2minEdge[u2Pos][v2Index] = dne.id;
        }
      }
      if (prevVal != CSNodes[u2].E2min[u2Pos][v2Index]) updated = true;
    } else {
      if (neighborEdges.size() == 0) {
        CSNodes[u2].E2min[u2Pos][v2Index] = -1;
        CSNodes[u2].E2minEdge[u2Pos][v2Index] = -1;
        updated = true;
        CSNodes[u2].E2minLazy[u2Pos][v2Index] = false;
        continue;
      }
      if (!CSNodes[u2].E2minLazy[u2Pos][v2Index]) {
        if (CSNodes[u2].E2min[u2Pos][v2Index] != de.id) continue;
        if (E2CheckNow[qeID][deID]) continue;
      }

      // if D_2[u2][v2] == 0 then pass
      if (CSNodes[u2].mark[v2Index] < 2) {
        CSNodes[u2].E2minLazy[u2Pos][v2Index] = true;
        continue;
      }
      CSNodes[u2].E2minLazy[u2Pos][v2Index] = false;

      int prevVal = CSNodes[u2].E2min[u2Pos][v2Index];
      CSNodes[u2].E2min[u2Pos][v2Index] = -1;
      CSNodes[u2].E2minEdge[u2Pos][v2Index] = -1;

      for (auto dne : neighborEdges) {
        if (!dne.isMatch(qe) || !E2CheckNow[qeID][dne.id]) continue;
        if (CSNodes[u2].E2min[u2Pos][v2Index] == -1 ||
            G.stream[CSNodes[u2].E2min[u2Pos][v2Index]].time >
                G.stream[dne.id].time) {
          CSNodes[u2].E2min[u2Pos][v2Index] = dne.id;
          CSNodes[u2].E2minEdge[u2Pos][v2Index] = dne.id;
          break;
        }
      }
      if (prevVal != CSNodes[u2].E2min[u2Pos][v2Index]) updated = true;
    }
  }

  // E2max update
  start = Q.maxPathTreeIndex[u2][u1Index];
  end = Q.maxPathTreeIndex[u2][u1Index + 1];
  for (int u2Pos = start; u2Pos < end; u2Pos++) {
    int u1Pos = Q.maxPathTreePos[u2][u2Pos];
    if (u1Pos != -1) {
      if (CSNodes[u2].E2max[u2Pos][v2Index] == -1) continue;
      if (neighborEdges.size() == 0) {
        CSNodes[u2].E2max[u2Pos][v2Index] = -1;
        CSNodes[u2].E2maxEdge[u2Pos][v2Index] = -1;
        updated = true;
        CSNodes[u2].E2maxLazy[u2Pos][v2Index] = false;
        continue;
      }
      if (!CSNodes[u2].E2maxLazy[u2Pos][v2Index]) {
        if (CSNodes[u2].E2maxEdge[u2Pos][v2Index] != de.id) continue;
        if (CSNodes[u2].E2max[u2Pos][v2Index] ==
                CSNodes[u1].E2max[u1Pos][v1Index] &&
            E2CheckNow[qeID][deID])
          continue;
      }

      // if D_2[u2][v2] == 0 then pass
      if (CSNodes[u2].mark[v2Index] < 2) {
        CSNodes[u2].E2maxLazy[u2Pos][v2Index] = true;
        continue;
      }
      CSNodes[u2].E2maxLazy[u2Pos][v2Index] = false;

      int prevVal = CSNodes[u2].E2max[u2Pos][v2Index];
      CSNodes[u2].E2max[u2Pos][v2Index] = -1;
      CSNodes[u2].E2maxEdge[u2Pos][v2Index] = -1;

      for (auto dne : neighborEdges) {
        if (!dne.isMatch(qe) || !E2CheckNow[qeID][dne.id]) continue;
        int vn = (v2Out) ? dne.v : dne.u;
        int vnIndex = G.vertexIDByLabel[vn];

        if (CSNodes[u1].E2max[u1Pos][vnIndex] == -1) continue;
        if (CSNodes[u2].E2max[u2Pos][v2Index] == -1 ||
            G.stream[CSNodes[u2].E2max[u2Pos][v2Index]].time <
                G.stream[CSNodes[u1].E2max[u1Pos][vnIndex]].time) {
          CSNodes[u2].E2max[u2Pos][v2Index] = CSNodes[u1].E2max[u1Pos][vnIndex];
          CSNodes[u2].E2maxEdge[u2Pos][v2Index] = dne.id;
        }
      }
      if (prevVal != CSNodes[u2].E2max[u2Pos][v2Index]) updated = true;
    } else {
      if (neighborEdges.size() == 0) {
        CSNodes[u2].E2max[u2Pos][v2Index] = -1;
        CSNodes[u2].E2maxEdge[u2Pos][v2Index] = -1;
        updated = true;
        CSNodes[u2].E2maxLazy[u2Pos][v2Index] = false;
        continue;
      }
      if (!CSNodes[u2].E2maxLazy[u2Pos][v2Index]) {
        if (CSNodes[u2].E2max[u2Pos][v2Index] != de.id) continue;
        if (E2CheckNow[qeID][deID]) continue;
      }

      // if D_2[u2][v2] == 0 then pass
      if (CSNodes[u2].mark[v2Index] < 2) {
        CSNodes[u2].E2maxLazy[u2Pos][v2Index] = true;
        continue;
      }
      CSNodes[u2].E2maxLazy[u2Pos][v2Index] = false;

      int prevVal = CSNodes[u2].E2max[u2Pos][v2Index];
      CSNodes[u2].E2max[u2Pos][v2Index] = -1;
      CSNodes[u2].E2maxEdge[u2Pos][v2Index] = -1;

      for (auto rit = neighborEdges.rbegin(); rit != neighborEdges.rend();
           rit++) {
        auto dne = *rit;
        if (!dne.isMatch(qe) || !E2CheckNow[qeID][dne.id]) continue;
        if (CSNodes[u2].E2max[u2Pos][v2Index] == -1 ||
            G.stream[CSNodes[u2].E2max[u2Pos][v2Index]].time <
                G.stream[dne.id].time) {
          CSNodes[u2].E2max[u2Pos][v2Index] = dne.id;
          CSNodes[u2].E2maxEdge[u2Pos][v2Index] = dne.id;
          break;
        }
      }
      if (prevVal != CSNodes[u2].E2max[u2Pos][v2Index]) updated = true;
    }
  }

  if (updated && !QEmaxVisited[u2][v2Index]) {
    QEmax.emplace(Q.dagOrderInv[u2], v2);
    QEmaxVisited[u2][v2Index] = true;
    QEmaxVertices[u2].push_back(v2);
  }

  E2CheckPrev[qeID][deID] = E2CheckNow[qeID][deID];
  if (E1CheckPrev[qeID][deID] && !E2CheckPrev[qeID][deID])
    TCMCSEdges.emplace_back(qeID, deID);
}

void CSStructure::insertionUpdateD(vector<CSEdge>& CSEdges,
                                   vector<CSEdge>& NLFEdges) {
  for (auto CSEdge : NLFEdges) {
    Edge& qe = Q.edges[CSEdge.qeID];
    Edge& de = G.stream[CSEdge.deID];

    int u1 = qe.u, u2 = qe.v, v1 = de.u, v2 = de.v;
    if (Q.dagOrderInv[u1] > Q.dagOrderInv[u2]) {
      swap(u1, u2);
      swap(v1, v2);
    }

    int v1Index = G.vertexIDByLabel[v1];
    int v2Index = G.vertexIDByLabel[v2];

    // NLF
    if (CSNodes[u1].mark[v1Index] == 0 &&
        CSNodes[u1].parentCount[v1Index] == Q.dagParent[u1].size() &&
        !CSNodes[u1].NLFCheck[v1Index]) {
      bool newNLFCheck = computeNLF(u1, v1);
      if (newNLFCheck) {
        CSNodes[u1].NLFCheck[v1Index] = true;
        Qmin.emplace(Q.dagOrderInv[u1], v1);
        if (CSNodes[u1].childCount[v1Index] == Q.dagChild[u1].size()) {
          Qmax.emplace(Q.dagOrderInv[u1], v1);
        }
      }
    }
    CSNodes[u1].NLFCheck[v1Index] = computeNLF(u1, v1);

    if (CSNodes[u2].mark[v2Index] == 0 &&
        CSNodes[u2].parentCount[v2Index] == Q.dagParent[u2].size() &&
        !CSNodes[u2].NLFCheck[v2Index]) {
      bool newNLFCheck = computeNLF(u2, v2);
      if (newNLFCheck) {
        CSNodes[u2].NLFCheck[v2Index] = true;
        Qmin.emplace(Q.dagOrderInv[u2], v2);
        if (CSNodes[u2].childCount[v2Index] == Q.dagChild[u2].size()) {
          Qmax.emplace(Q.dagOrderInv[u2], v2);
        }
      }
    }
    CSNodes[u2].NLFCheck[v2Index] = computeNLF(u2, v2);
  }

  for (auto CSEdge : CSEdges) {
    Edge& qe = Q.edges[CSEdge.qeID];
    Edge& de = G.stream[CSEdge.deID];

    int u1 = qe.u, u2 = qe.v, v1 = de.u, v2 = de.v;
    if (Q.dagOrderInv[u1] > Q.dagOrderInv[u2]) {
      swap(u1, u2);
      swap(v1, v2);
    }

    int v1Index = G.vertexIDByLabel[v1];
    int v2Index = G.vertexIDByLabel[v2];

    if (CSNodes[u1].mark[v1Index] >= 1) insertionUpdateDTopDown(u1, v1, u2, v2);
    if (CSNodes[u1].mark[v1Index] == 2) {
      int neighborIndex = Q.neighborIndex[u2][u1];
      int neighborCandIndex = v2Index * Q.neighbors[u2].size() + neighborIndex;
      CSNodes[u2].validNeighborCand[neighborCandIndex]++;
    }
    if (CSNodes[u2].mark[v2Index] == 2) {
      insertionUpdateDBottomUp(u1, v1, u2, v2);
      int neighborIndex = Q.neighborIndex[u1][u2];
      int neighborCandIndex = v1Index * Q.neighbors[u1].size() + neighborIndex;
      CSNodes[u1].validNeighborCand[neighborCandIndex]++;
    }
  }

  while (!Qmin.empty()) {
    int u = Q.dagOrder[Qmin.front().first];
    int v = Qmin.front().second;
    int vIndex = G.vertexIDByLabel[v];
    Qmin.pop();

    CSNodes[u].mark[vIndex] = 1;

    for (auto qe : Q.dagChildOutEdges[u]) {
      int uc = qe.v;
      for (auto de : G.outEdges[v][Q.edgeHash[qe.id]]) {
        if (de.isMatch(qe) && E1CheckPrev[qe.id][de.id] &&
            E2CheckPrev[qe.id][de.id]) {
          insertionUpdateDTopDown(u, v, uc, de.v);
        }
      }
    }

    for (auto qe : Q.dagChildInEdges[u]) {
      int uc = qe.u;
      for (auto de : G.inEdges[v][Q.edgeHash[qe.id]]) {
        if (de.isMatch(qe) && E1CheckPrev[qe.id][de.id] &&
            E2CheckPrev[qe.id][de.id]) {
          insertionUpdateDTopDown(u, v, uc, de.u);
        }
      }
    }

    if (CSNodes[u].childCount[vIndex] == Q.dagChild[u].size())
      Qmax.emplace(Q.dagOrderInv[u], v);
  }

  while (!Qmax.empty()) {
    int u = Q.dagOrder[Qmax.front().first];
    int v = Qmax.front().second;
    int vIndex = G.vertexIDByLabel[v];
    Qmax.pop();

    CSNodes[u].mark[vIndex] = 2;

    for (auto qe : Q.dagParentOutEdges[u]) {
      int up = qe.v;
      for (auto de : G.outEdges[v][Q.edgeHash[qe.id]]) {
        if (de.isMatch(qe) && E1CheckPrev[qe.id][de.id] &&
            E2CheckPrev[qe.id][de.id]) {
          insertionUpdateDBottomUp(up, de.v, u, v);
        }
      }
    }

    for (auto qe : Q.dagParentInEdges[u]) {
      int up = qe.u;
      for (auto de : G.inEdges[v][Q.edgeHash[qe.id]]) {
        if (de.isMatch(qe) && E1CheckPrev[qe.id][de.id] &&
            E2CheckPrev[qe.id][de.id]) {
          insertionUpdateDBottomUp(up, de.u, u, v);
        }
      }
    }

    for (auto qe : Q.outEdges[u]) {
      int un = qe.v;
      for (auto de : G.outEdges[v][Q.edgeHash[qe.id]]) {
        if (de.isMatch(qe) && E1CheckPrev[qe.id][de.id] &&
            E2CheckPrev[qe.id][de.id]) {
          int vn = de.v;
          int vnIndex = G.vertexIDByLabel[vn];
          int neighborIndex = Q.neighborIndex[un][u];
          int neighborCandIndex =
              vnIndex * Q.neighbors[un].size() + neighborIndex;
          CSNodes[un].validNeighborCand[neighborCandIndex]++;
        }
      }
    }

    for (auto qe : Q.inEdges[u]) {
      int un = qe.u;
      for (auto de : G.inEdges[v][Q.edgeHash[qe.id]]) {
        if (de.isMatch(qe) && E1CheckPrev[qe.id][de.id] &&
            E2CheckPrev[qe.id][de.id]) {
          int vn = de.u;
          int vnIndex = G.vertexIDByLabel[vn];
          int neighborIndex = Q.neighborIndex[un][u];
          int neighborCandIndex =
              vnIndex * Q.neighbors[un].size() + neighborIndex;
          CSNodes[un].validNeighborCand[neighborCandIndex]++;
        }
      }
    }
  }
}

void CSStructure::insertionUpdateDTopDown(int u, int v, int uc, int vc) {
  int vcIndex = G.vertexIDByLabel[vc];
  int parentIndex = Q.dagIndex[uc][u];
  int parentCandIndex = vcIndex * Q.dagParent[uc].size() + parentIndex;

  if (CSNodes[uc].validParentCand[parentCandIndex] == 0) {
    CSNodes[uc].parentCount[vcIndex]++;

    if (CSNodes[uc].parentCount[vcIndex] == Q.dagParent[uc].size()) {
      if (CSNodes[uc].NLFCheck[vcIndex]) Qmin.emplace(Q.dagOrderInv[uc], vc);
    }
  }
  CSNodes[uc].validParentCand[parentCandIndex]++;
}

void CSStructure::insertionUpdateDBottomUp(int up, int vp, int u, int v) {
  int vpIndex = G.vertexIDByLabel[vp];
  int childIndex = Q.dagIndex[up][u];
  int childCandIndex = vpIndex * Q.dagChild[up].size() + childIndex;

  if (CSNodes[up].validChildCand[childCandIndex] == 0) {
    CSNodes[up].childCount[vpIndex]++;
    if (CSNodes[up].mark[vpIndex] == 1 &&
        CSNodes[up].childCount[vpIndex] == Q.dagChild[up].size()) {
      Qmax.emplace(Q.dagOrderInv[up], vp);
    }
  }
  CSNodes[up].validChildCand[childCandIndex]++;
}

void CSStructure::deletionUpdateD(vector<CSEdge>& CSEdges,
                                  vector<CSEdge>& NLFEdges) {
  for (auto CSEdge : NLFEdges) {
    Edge& qe = Q.edges[CSEdge.qeID];
    Edge& de = G.stream[CSEdge.deID];

    int u1 = qe.u, u2 = qe.v, v1 = de.u, v2 = de.v;
    if (Q.dagOrderInv[u1] > Q.dagOrderInv[u2]) {
      swap(u1, u2);
      swap(v1, v2);
    }

    int v1Index = G.vertexIDByLabel[v1];
    int v2Index = G.vertexIDByLabel[v2];

    // NLF
    if (CSNodes[u1].mark[v1Index] >= 1 && CSNodes[u1].NLFCheck[v1Index]) {
      bool newNLFCheck = computeNLF(u1, v1);
      if (!newNLFCheck) {
        CSNodes[u1].NLFCheck[v1Index] = false;
        Qmin.emplace(Q.dagOrderInv[u1], v1);
      }
    }
    CSNodes[u1].NLFCheck[v1Index] = computeNLF(u1, v1);

    if (CSNodes[u2].mark[v2Index] >= 1 && CSNodes[u2].NLFCheck[v2Index]) {
      bool newNLFCheck = computeNLF(u2, v2);
      if (!newNLFCheck) {
        CSNodes[u2].NLFCheck[v2Index] = false;
        Qmin.emplace(Q.dagOrderInv[u2], v2);
      }
    }
    CSNodes[u2].NLFCheck[v2Index] = computeNLF(u2, v2);
  }

  for (auto CSEdge : CSEdges) {
    Edge& qe = Q.edges[CSEdge.qeID];
    Edge& de = G.stream[CSEdge.deID];

    int u1 = qe.u, u2 = qe.v, v1 = de.u, v2 = de.v;
    if (Q.dagOrderInv[u1] > Q.dagOrderInv[u2]) {
      swap(u1, u2);
      swap(v1, v2);
    }

    int v1Index = G.vertexIDByLabel[v1];
    int v2Index = G.vertexIDByLabel[v2];

    if (CSNodes[u1].mark[v1Index] >= 1) deletionUpdateDTopDown(u1, v1, u2, v2);
    if (CSNodes[u1].mark[v1Index] == 2) {
      int neighborIndex = Q.neighborIndex[u2][u1];
      int neighborCandIndex = v2Index * Q.neighbors[u2].size() + neighborIndex;
      CSNodes[u2].validNeighborCand[neighborCandIndex]--;
    }
    if (CSNodes[u2].mark[v2Index] == 2) {
      deletionUpdateDBottomUp(u1, v1, u2, v2);
      int neighborIndex = Q.neighborIndex[u1][u2];
      int neighborCandIndex = v1Index * Q.neighbors[u1].size() + neighborIndex;
      CSNodes[u1].validNeighborCand[neighborCandIndex]--;
    }
  }

  while (!Qmin.empty()) {
    int u = Q.dagOrder[Qmin.front().first];
    int v = Qmin.front().second;
    int vIndex = G.vertexIDByLabel[v];
    Qmin.pop();

    if (CSNodes[u].mark[vIndex] == 2) {
      for (auto qe : Q.dagParentOutEdges[u]) {
        int up = qe.v;
        for (auto de : G.outEdges[v][Q.edgeHash[qe.id]]) {
          if (de.isMatch(qe) && E1CheckPrev[qe.id][de.id] &&
              E2CheckPrev[qe.id][de.id]) {
            deletionUpdateDBottomUp(up, de.v, u, v);
          }
        }
      }

      for (auto qe : Q.dagParentInEdges[u]) {
        int up = qe.u;
        for (auto de : G.inEdges[v][Q.edgeHash[qe.id]]) {
          if (de.isMatch(qe) && E1CheckPrev[qe.id][de.id] &&
              E2CheckPrev[qe.id][de.id]) {
            deletionUpdateDBottomUp(up, de.u, u, v);
          }
        }
      }

      for (auto qe : Q.outEdges[u]) {
        int un = qe.v;
        for (auto de : G.outEdges[v][Q.edgeHash[qe.id]]) {
          if (de.isMatch(qe) && E1CheckPrev[qe.id][de.id] &&
              E2CheckPrev[qe.id][de.id]) {
            int vn = de.v;
            int vnIndex = G.vertexIDByLabel[vn];
            int neighborIndex = Q.neighborIndex[un][u];
            int neighborCandIndex =
                vnIndex * Q.neighbors[un].size() + neighborIndex;
            CSNodes[un].validNeighborCand[neighborCandIndex]--;
          }
        }
      }

      for (auto qe : Q.inEdges[u]) {
        int un = qe.u;
        for (auto de : G.inEdges[v][Q.edgeHash[qe.id]]) {
          if (de.isMatch(qe) && E1CheckPrev[qe.id][de.id] &&
              E2CheckPrev[qe.id][de.id]) {
            int vn = de.u;
            int vnIndex = G.vertexIDByLabel[vn];
            int neighborIndex = Q.neighborIndex[un][u];
            int neighborCandIndex =
                vnIndex * Q.neighbors[un].size() + neighborIndex;
            CSNodes[un].validNeighborCand[neighborCandIndex]--;
          }
        }
      }
    }

    for (auto qe : Q.dagChildOutEdges[u]) {
      int uc = qe.v;
      for (auto de : G.outEdges[v][Q.edgeHash[qe.id]]) {
        if (de.isMatch(qe) && E1CheckPrev[qe.id][de.id] &&
            E2CheckPrev[qe.id][de.id]) {
          deletionUpdateDTopDown(u, v, uc, de.v);
        }
      }
    }

    for (auto qe : Q.dagChildInEdges[u]) {
      int uc = qe.u;
      for (auto de : G.inEdges[v][Q.edgeHash[qe.id]]) {
        if (de.isMatch(qe) && E1CheckPrev[qe.id][de.id] &&
            E2CheckPrev[qe.id][de.id]) {
          deletionUpdateDTopDown(u, v, uc, de.u);
        }
      }
    }

    CSNodes[u].mark[vIndex] = 0;
  }

  while (!Qmax.empty()) {
    int u = Q.dagOrder[Qmax.front().first];
    int v = Qmax.front().second;
    int vIndex = G.vertexIDByLabel[v];
    Qmax.pop();
    if (CSNodes[u].mark[vIndex] != 2) continue;

    for (auto qe : Q.dagParentOutEdges[u]) {
      int up = qe.v;
      for (auto de : G.outEdges[v][Q.edgeHash[qe.id]]) {
        if (de.isMatch(qe) && E1CheckPrev[qe.id][de.id] &&
            E2CheckPrev[qe.id][de.id]) {
          deletionUpdateDBottomUp(up, de.v, u, v);
        }
      }
    }

    for (auto qe : Q.dagParentInEdges[u]) {
      int up = qe.u;
      for (auto de : G.inEdges[v][Q.edgeHash[qe.id]]) {
        if (de.isMatch(qe) && E1CheckPrev[qe.id][de.id] &&
            E2CheckPrev[qe.id][de.id]) {
          deletionUpdateDBottomUp(up, de.u, u, v);
        }
      }
    }

    for (auto qe : Q.outEdges[u]) {
      int un = qe.v;
      for (auto de : G.outEdges[v][Q.edgeHash[qe.id]]) {
        if (de.isMatch(qe) && E1CheckPrev[qe.id][de.id] &&
            E2CheckPrev[qe.id][de.id]) {
          int vn = de.v;
          int vnIndex = G.vertexIDByLabel[vn];
          int neighborIndex = Q.neighborIndex[un][u];
          int neighborCandIndex =
              vnIndex * Q.neighbors[un].size() + neighborIndex;
          CSNodes[un].validNeighborCand[neighborCandIndex]--;
        }
      }
    }

    for (auto qe : Q.inEdges[u]) {
      int un = qe.u;
      for (auto de : G.inEdges[v][Q.edgeHash[qe.id]]) {
        if (de.isMatch(qe) && E1CheckPrev[qe.id][de.id] &&
            E2CheckPrev[qe.id][de.id]) {
          int vn = de.u;
          int vnIndex = G.vertexIDByLabel[vn];
          int neighborIndex = Q.neighborIndex[un][u];
          int neighborCandIndex =
              vnIndex * Q.neighbors[un].size() + neighborIndex;
          CSNodes[un].validNeighborCand[neighborCandIndex]--;
        }
      }
    }

    CSNodes[u].mark[vIndex] = 1;
  }
}

void CSStructure::deletionUpdateDTopDown(int u, int v, int uc, int vc) {
  int vcIndex = G.vertexIDByLabel[vc];
  int parentIndex = Q.dagIndex[uc][u];
  int parentCandIndex = vcIndex * Q.dagParent[uc].size() + parentIndex;

  if (CSNodes[uc].validParentCand[parentCandIndex] == 1) {
    if (CSNodes[uc].mark[vcIndex] >= 1 &&
        CSNodes[uc].parentCount[vcIndex] == Q.dagParent[uc].size()) {
      if (CSNodes[uc].NLFCheck[vcIndex]) Qmin.emplace(Q.dagOrderInv[uc], vc);
    }
    CSNodes[uc].parentCount[vcIndex]--;
  }
  CSNodes[uc].validParentCand[parentCandIndex]--;
}

void CSStructure::deletionUpdateDBottomUp(int up, int vp, int u, int v) {
  int vpIndex = G.vertexIDByLabel[vp];
  int childIndex = Q.dagIndex[up][u];
  int childCandIndex = vpIndex * Q.dagChild[up].size() + childIndex;

  if (CSNodes[up].validChildCand[childCandIndex] == 1) {
    if (CSNodes[up].mark[vpIndex] == 2 &&
        CSNodes[up].childCount[vpIndex] == Q.dagChild[up].size()) {
      if (CSNodes[up].NLFCheck[vpIndex]) Qmax.emplace(Q.dagOrderInv[up], vp);
    }
    CSNodes[up].childCount[vpIndex]--;
  }
  CSNodes[up].validChildCand[childCandIndex]--;
}

vector<int> CSStructure::computeCm(int u, int lvl) {
  vector<int> ret;
  int ub = CmMinNeighbor[lvl][u];
  int vb = vertexMatch[ub];
  CmSizes[u] = 0;

  for (auto qe : Q.outEdges[u]) {
    if (qe.v != ub) {
      if (vertexMatch[qe.v] != -1) ret.push_back(qe.id);
      continue;
    }
    ret.push_back(qe.id);
    swap(ret[0], ret[ret.size() - 1]);

    int size = 0;

    mTime lowTime(0, 0);
    for (auto beID : Q.beforeEdges[qe.id]) {
      int deID = edgeMatch[beID];
      if (deID != -1 && lowTime < G.stream[deID].time) {
        lowTime = G.stream[deID].time;
      }
    }

    mTime highTime(numeric_limits<int>::max(), numeric_limits<int>::max());
    for (auto aeID : Q.afterEdges[qe.id]) {
      int deID = edgeMatch[aeID];
      if (deID != -1 && highTime > G.stream[deID].time) {
        highTime = G.stream[deID].time;
      }
    }
    if (lowTime > highTime) continue;

    std::deque<Edge>::iterator lo, hi;
    lo = lower_bound(G.inEdges[vb][Q.edgeHash[qe.id]].begin(),
                     G.inEdges[vb][Q.edgeHash[qe.id]].end(),
                     lowTime + mTime(0, 1));
    hi = upper_bound(G.inEdges[vb][Q.edgeHash[qe.id]].begin(),
                     G.inEdges[vb][Q.edgeHash[qe.id]].end(), highTime);

    for (auto it = lo; it != hi; it++) {
      auto de = *it;
      if (!de.isMatch(qe) || !E1CheckPrev[qe.id][de.id] ||
          !E2CheckPrev[qe.id][de.id])
        continue;

      int v = de.u;
      if (visitedData[v] || checkData[v]) continue;
      int vIndex = G.vertexIDByLabel[v];
      if (CSNodes[u].mark[vIndex] != 2) continue;

      CmVertices[u][size] = v;
      size++;
      checkData[v] = true;
    }
    CmSizes[u] = size;
    for (int i = 0; i < size; i++) checkData[CmVertices[u][i]] = false;
  }

  for (auto qe : Q.inEdges[u]) {
    if (qe.u != ub) {
      if (vertexMatch[qe.u] != -1) ret.push_back(qe.id);
      continue;
    }
    ret.push_back(qe.id);
    swap(ret[0], ret[ret.size() - 1]);

    int size = 0;

    mTime lowTime(0, 0);
    for (auto beID : Q.beforeEdges[qe.id]) {
      int deID = edgeMatch[beID];
      if (deID != -1 && lowTime < G.stream[deID].time) {
        lowTime = G.stream[deID].time;
      }
    }

    mTime highTime(numeric_limits<int>::max(), numeric_limits<int>::max());
    for (auto aeID : Q.afterEdges[qe.id]) {
      int deID = edgeMatch[aeID];
      if (deID != -1 && highTime > G.stream[deID].time) {
        highTime = G.stream[deID].time;
      }
    }
    if (lowTime > highTime) continue;

    std::deque<Edge>::iterator lo, hi;
    lo = lower_bound(G.outEdges[vb][Q.edgeHash[qe.id]].begin(),
                     G.outEdges[vb][Q.edgeHash[qe.id]].end(),
                     lowTime + mTime(0, 1));
    hi = upper_bound(G.outEdges[vb][Q.edgeHash[qe.id]].begin(),
                     G.outEdges[vb][Q.edgeHash[qe.id]].end(), highTime);

    for (auto it = lo; it != hi; it++) {
      auto de = *it;
      if (!de.isMatch(qe) || !E1CheckPrev[qe.id][de.id] ||
          !E2CheckPrev[qe.id][de.id])
        continue;

      int v = de.v;
      if (visitedData[v] || checkData[v]) continue;
      int vIndex = G.vertexIDByLabel[v];
      if (CSNodes[u].mark[vIndex] != 2) continue;

      CmVertices[u][size] = v;
      size++;
      checkData[v] = true;
    }
    CmSizes[u] = size;
    for (int i = 0; i < size; i++) checkData[CmVertices[u][i]] = false;
  }

  return ret;
}

void CSStructure::computeECm(int qeID, bool& existBefore, bool& existAfter) {
  ECm[qeID].clear();

  mTime lowTime(0, 0);
  for (auto beID : Q.beforeEdges[qeID]) {
    int deID = edgeMatch[beID];
    if (deID != -1 && lowTime < G.stream[deID].time) {
      lowTime = G.stream[deID].time;
    }
    if (deID == -1) existBefore = true;
  }

  mTime highTime(numeric_limits<int>::max(), numeric_limits<int>::max());
  for (auto aeID : Q.afterEdges[qeID]) {
    int deID = edgeMatch[aeID];
    if (deID != -1 && highTime > G.stream[deID].time) {
      highTime = G.stream[deID].time;
    }
    if (deID == -1) existAfter = true;
  }
  if (lowTime > highTime) return;

  auto qe = Q.edges[qeID];
  int u1 = qe.u, v1 = vertexMatch[u1];
  int u2 = qe.v, v2 = vertexMatch[u2];
  std::deque<Edge>::iterator lo, hi;
  lo = lower_bound(G.outEdges[v1][Q.edgeHash[qe.id]].begin(),
                   G.outEdges[v1][Q.edgeHash[qe.id]].end(),
                   lowTime + mTime(0, 1));
  hi = upper_bound(G.outEdges[v1][Q.edgeHash[qe.id]].begin(),
                   G.outEdges[v1][Q.edgeHash[qe.id]].end(), highTime);

  for (auto it = lo; it != hi; it++) {
    auto de = *it;
    if (de.v != v2 || !de.isMatch(qe) || !E1CheckPrev[qe.id][de.id] ||
        !E2CheckPrev[qe.id][de.id])
      continue;
    ECm[qeID].push_back(de.id);
  }
}

void CSStructure::updateVertexMatch(int u, int v, int lvl) {
  vertexMatch[u] = v;
  visitedQuery[u] = true;
  visitedData[v] = true;
  // cout << "lvl: " << lvl << ", " << u << "->" << v << " match" << endl;

  if (lvl > 0) {
    for (int i = 0; i < Q.numVertices; i++) {
      EuSizes[lvl][i] = EuSizes[lvl - 1][i];
      CmMinNeighbor[lvl][i] = CmMinNeighbor[lvl - 1][i];
    }
  }

  int vIndex = G.vertexIDByLabel[v];

  int iterCount;
  for (iterCount = 0; iterCount < Q.dagParent[u].size(); iterCount++) {
    int u2 = Q.dagParent[u][iterCount];
    if (visitedQuery[u2]) continue;

    int candIndex = vIndex * Q.neighbors[u].size() + Q.neighborIndex[u][u2];
    if (extendibleQuery[u2] == 0 ||
        EuSizes[lvl][u2] > CSNodes[u].validNeighborCand[candIndex]) {
      EuSizes[lvl][u2] = CSNodes[u].validNeighborCand[candIndex];
      CmMinNeighbor[lvl][u2] = u;
    }

    extendibleQuery[u2]++;
    if (extendibleQuery[u2] == Q.degree[u2]) {
      isIsolated[u2] = true;
    }
  }

  for (iterCount = 0; iterCount < Q.dagChild[u].size(); iterCount++) {
    int u2 = Q.dagChild[u][iterCount];
    if (visitedQuery[u2]) continue;

    int candIndex = vIndex * Q.neighbors[u].size() + Q.neighborIndex[u][u2];
    if (extendibleQuery[u2] == 0 ||
        EuSizes[lvl][u2] > CSNodes[u].validNeighborCand[candIndex]) {
      EuSizes[lvl][u2] = CSNodes[u].validNeighborCand[candIndex];
      CmMinNeighbor[lvl][u2] = u;
    }

    extendibleQuery[u2]++;
    if (extendibleQuery[u2] == Q.degree[u2]) {
      isIsolated[u2] = true;
    }
  }
}

void CSStructure::restoreVertexMatch(int u, int v, int lvl) {
  visitedQuery[u] = false;
  visitedData[v] = false;
  vertexMatch[u] = -1;
  // cout << "lvl: " << lvl << ", " << u << "->" << v << " unmatch" << endl;

  int iterCount;
  for (iterCount = 0; iterCount < Q.dagParent[u].size(); iterCount++) {
    int u2 = Q.dagParent[u][iterCount];
    if (visitedQuery[u2]) continue;

    if (extendibleQuery[u2] == Q.degree[u2]) {
      isIsolated[u2] = false;
    }
    extendibleQuery[u2]--;
  }
  for (iterCount = 0; iterCount < Q.dagChild[u].size(); iterCount++) {
    int u2 = Q.dagChild[u][iterCount];
    if (visitedQuery[u2]) continue;

    if (extendibleQuery[u2] == Q.degree[u2]) {
      isIsolated[u2] = false;
    }
    extendibleQuery[u2]--;
  }
}

void CSStructure::updateEdgeMatch(int qeID, int deID, bool dummy = false) {
  edgeMatch[qeID] = deID;
  depth += 1;
  if (!dummy) {
    numMatchedEdges += 1;
    if (numMatchedEdges == 1)
      minEdge[numMatchedEdges] = deID;
    else {
      minEdge[numMatchedEdges] = minEdge[numMatchedEdges - 1];
      if (G.stream[minEdge[numMatchedEdges]].time > G.stream[deID].time)
        minEdge[numMatchedEdges] = deID;
    }
  }
  // cout << "edge: " << qeID << "(" << Q.edges[qeID].u << "->" <<
  // Q.edges[qeID].v << ") -> " << deID << "(" << G.stream[deID].u << "->" <<
  // G.stream[deID].v << ") match" << endl;
}

void CSStructure::restoreEdgeMatch(int qeID, int deID, bool dummy = false) {
  edgeMatch[qeID] = -1;
  depth -= 1;
  if (!dummy) numMatchedEdges -= 1;
  // cout << "edge: " << qeID << "(" << Q.edges[qeID].u << "->" <<
  // Q.edges[qeID].v << ") -> " << deID << "(" << G.stream[deID].u << "->" <<
  // G.stream[deID].v << ") unmatch" << endl;
}

void CSStructure::computeNumMatches(int qeID) {
  if (qeID == Q.numEdges) {
    numMatches += 1;
    storedMatches[minEdge[numMatchedEdges]] += 1;
    return;
  }

  if (isProduct[qeID]) {
    for (auto deID : ECm[qeID]) {
      updateEdgeMatch(qeID, deID);
      computeNumMatches(qeID + 1);
      restoreEdgeMatch(qeID, deID);
    }
  } else {
    computeNumMatches(qeID + 1);
  }
}

bool CSStructure::backtrack(int lvl, vector<int>& checkEdges, int idx = 0) {
  if (lvl == Q.numVertices) {
    computeNumMatches(0);
    return true;
  }

  bool ret = false;
  // choose next vertex according to matching order
  if (idx == checkEdges.size()) {
    int u = -1;
    for (int i = 0; i < Q.numVertices; i++) {
      // 1. if u is already matched / not extendible, then pass
      if (visitedQuery[i] || !extendibleQuery[i]) continue;

      // 2. choose u whose E(u) is minimum
      if (u == -1 || EuSizes[lvl - 1][i] < EuSizes[lvl - 1][u]) u = i;

      // 3. choose u which has the maximum number of non-matched neighbor when
      // |Cm[u]| is same
      if (EuSizes[lvl - 1][i] == EuSizes[lvl - 1][u])
        if (Q.degree[i] - extendibleQuery[i] > Q.degree[u] - extendibleQuery[u])
          u = i;
    }
    auto newCheckEdges = computeCm(u, lvl - 1);

    // |ECm((umin, u))|=0 and newCheckEdges[0] = (umin, u)
    if (CmSizes[u] == 0) {
      int eID = newCheckEdges[0];
      for (int i = 0; i < ((Q.numEdges + 63) >> 6); i++) {
        edgeFailingSet[depth][i] = Q.relatedEdges[eID][i];
      }
    } else {
      for (int i = 0; i < ((Q.numEdges + 63) >> 6); i++) {
        curEdgeFailingSet[lvl][i] = 0ULL;
      }
      for (int i = 0; i < CmSizes[u]; i++) {
        int v = CmVertices[u][i];
        updateVertexMatch(u, v, lvl);
        ret |= backtrack(lvl, newCheckEdges, 0);
        restoreVertexMatch(u, v, lvl);
        for (int i = 0; i < ((Q.numEdges + 63) >> 6); i++) {
          curEdgeFailingSet[lvl][i] |= edgeFailingSet[depth][i];
        }
      }
      for (int i = 0; i < ((Q.numEdges + 63) >> 6); i++) {
        edgeFailingSet[depth][i] = curEdgeFailingSet[lvl][i];
      }
    }
  } else {
    int qeID = checkEdges[idx];
    bool existBefore = false, existAfter = false;
    computeECm(qeID, existBefore, existAfter);
    isProduct[qeID] = (existBefore || existAfter) ? false : true;
    if (ECm[qeID].size() == 0) {
      for (int i = 0; i < ((Q.numEdges + 63) >> 6); i++) {
        edgeFailingSet[depth][i] = Q.relatedEdges[qeID][i];
      }
      return false;
    }

    for (int i = 0; i < ((Q.numEdges + 63) >> 6); i++) {
      edgeFailingSet[depth][i] = 0ULL;
    }

    // isProduct[qeID] = false;
    if (!existBefore && !existAfter)  // isProduct[qeID] == true
    {
      updateEdgeMatch(qeID, ECm[qeID][0], true);
      if (idx != checkEdges.size() - 1)
        ret |= backtrack(lvl, checkEdges, idx + 1);
      else
        ret |= backtrack(lvl + 1, checkEdges, idx + 1);
      restoreEdgeMatch(qeID, ECm[qeID][0], true);
      if (!ret) {
        for (int i = 0; i < ((Q.numEdges + 63) >> 6); i++) {
          edgeFailingSet[depth][i] = edgeFailingSet[depth + 1][i];
        }
        if (edgeFailingSet[depth][qeID >> 6] & (1ULL << (qeID & 0x3F))) {
          for (int i = 0; i < ((Q.numEdges + 63) >> 6); i++) {
            edgeFailingSet[depth][i] |= Q.relatedEdges[qeID][i];
          }
        }
      }
    } else if (!existBefore && existAfter) {
      for (auto deID : ECm[qeID]) {
        bool result;
        updateEdgeMatch(qeID, deID);
        if (idx != checkEdges.size() - 1)
          result = backtrack(lvl, checkEdges, idx + 1);
        else
          result = backtrack(lvl + 1, checkEdges, idx + 1);
        restoreEdgeMatch(qeID, deID);
        ret |= result;
        if (!ret) {
          for (int i = 0; i < ((Q.numEdges + 63) >> 6); i++) {
            edgeFailingSet[depth][i] = edgeFailingSet[depth + 1][i];
          }
          if (edgeFailingSet[depth][qeID >> 6] & (1ULL << (qeID & 0x3F))) {
            for (int i = 0; i < ((Q.numEdges + 63) >> 6); i++) {
              edgeFailingSet[depth][i] |= Q.relatedEdges[qeID][i];
            }
          }
        }
        if (!result) return ret;
      }
    } else if (existBefore && !existAfter) {
      for (int i = ECm[qeID].size() - 1; i >= 0; i--) {
        auto deID = ECm[qeID][i];
        bool result;
        updateEdgeMatch(qeID, deID);
        if (idx != checkEdges.size() - 1)
          result = backtrack(lvl, checkEdges, idx + 1);
        else
          result = backtrack(lvl + 1, checkEdges, idx + 1);
        restoreEdgeMatch(qeID, deID);
        ret |= result;
        if (!ret) {
          for (int i = 0; i < ((Q.numEdges + 63) >> 6); i++) {
            edgeFailingSet[depth][i] = edgeFailingSet[depth + 1][i];
          }
          if (edgeFailingSet[depth][qeID >> 6] & (1ULL << (qeID & 0x3F))) {
            for (int i = 0; i < ((Q.numEdges + 63) >> 6); i++) {
              edgeFailingSet[depth][i] |= Q.relatedEdges[qeID][i];
            }
          }
        }
        if (!result) return ret;
      }
    } else {
      for (auto deID : ECm[qeID]) {
        bool result;
        updateEdgeMatch(qeID, deID);
        if (idx != checkEdges.size() - 1)
          result = backtrack(lvl, checkEdges, idx + 1);
        else
          result = backtrack(lvl + 1, checkEdges, idx + 1);
        restoreEdgeMatch(qeID, deID);
        ret |= result;
        if (!ret) {
          if (!(edgeFailingSet[depth + 1][qeID >> 6] &
                (1ULL << (qeID & 0x3F)))) {
            for (int i = 0; i < ((Q.numEdges + 63) >> 6); i++) {
              edgeFailingSet[depth][i] = edgeFailingSet[depth + 1][i];
            }
            return ret;
          } else {
            for (int i = 0; i < ((Q.numEdges + 63) >> 6); i++) {
              edgeFailingSet[depth][i] |= edgeFailingSet[depth + 1][i];
            }
          }
          if (edgeFailingSet[depth][qeID >> 6] & (1ULL << (qeID & 0x3F))) {
            for (int i = 0; i < ((Q.numEdges + 63) >> 6); i++) {
              edgeFailingSet[depth][i] |= Q.relatedEdges[qeID][i];
            }
          }
        }
        if (ret && !result) return ret;
      }
    }
  }

  return ret;
}

void CSStructure::findMatches(vector<CSEdge>& CSEdges) {
  vector<int> idx;
  for (int i = 0; i < CSEdges.size(); i++) {
    auto CSEdge = CSEdges[i];
    Edge& qe = Q.edges[CSEdge.qeID];
    Edge& de = G.stream[CSEdge.deID];

    int u1 = qe.u, u2 = qe.v, v1 = de.u, v2 = de.v;

    int v1Index = G.vertexIDByLabel[v1];
    int v2Index = G.vertexIDByLabel[v2];
    if (!E1CheckPrev[qe.id][de.id] || !E2CheckPrev[qe.id][de.id]) continue;
    if (CSNodes[u1].mark[v1Index] != 2 || CSNodes[u2].mark[v2Index] != 2)
      continue;
    // if(Q.NECMapping[CSEdge.u1] != CSEdge.u1 || Q.NECMapping[CSEdge.u2] !=
    // CSEdge.u2) continue;

    idx.push_back(i);
  }
  if (idx.empty()) return;

  int prevu1 = -1, prevu2 = -1, prevv1 = -1, prevv2 = -1;
  for (int i = 0; i < idx.size(); i++) {
    auto CSEdge = CSEdges[idx[i]];
    Edge& qe = Q.edges[CSEdge.qeID];
    Edge& de = G.stream[CSEdge.deID];
    ECm[qe.id].clear();
    ECm[qe.id].push_back(de.id);
    isProduct[qe.id] = false;

    int u1 = qe.u, u2 = qe.v, v1 = de.u, v2 = de.v;

    updateVertexMatch(u1, v1, 0);
    updateVertexMatch(u2, v2, 1);
    updateEdgeMatch(qe.id, de.id);
    vector<int> vi;
    backtrack(2, vi);
    restoreEdgeMatch(qe.id, de.id);
    restoreVertexMatch(u2, v2, 1);
    restoreVertexMatch(u1, v1, 0);
  }
  numBacktrack += idx.size();
  numFiltered += CSEdges.size() - idx.size();
}

#endif
