#ifndef _GRAPH_H
#define _GRAPH_H

#include <algorithm>
#include <cstring>
#include <deque>
#include <iostream>
#include <map>
#include <tuple>
#include <vector>

#include "../util.hpp"

#define SORT_THRESHOLD 16
using namespace std;

int maxDataVertex = -1, maxQueryVertex = -1;
map<tuple<int, int, int>, int> NLFHash;
map<tuple<int, int, int>, int> queryEdgeHash;
int NLFBitSize, NLFArraySize;

struct Edge;
struct mTime {
  int sec, msec;

  mTime(const int _sec = -1, const int _msec = -1) : sec(_sec), msec(_msec) {}
  mTime(const mTime& other) : sec(other.sec), msec(other.msec) {}

  bool operator<(const mTime& other) const {
    if (sec < other.sec) return true;
    if (sec > other.sec) return false;
    return msec < other.msec;
  }

  bool operator>(const mTime& other) const { return other < *this; }

  bool operator==(const mTime& other) const {
    if (sec != other.sec) return false;
    if (msec != other.msec) return false;

    return true;
  }

  bool operator!=(const mTime& other) const { return !(*this == other); }

  mTime& operator=(const mTime& other) {
    sec = other.sec;
    msec = other.msec;
    return *this;
  }

  mTime& operator+=(const mTime& rhs) {
    this->msec += rhs.msec;
    this->sec += rhs.sec;
    if (this->msec >= 1e9) {
      this->msec -= 1e9;
      this->sec++;
    }
    return *this;
  }

  mTime& operator-=(const mTime& rhs) {
    this->msec -= rhs.msec;
    this->sec -= rhs.sec;
    if (this->msec < 0) {
      this->msec += 1e9;
      this->sec--;
    }
    return *this;
  }
};

inline mTime operator+(mTime lhs, const mTime& rhs) {
  lhs += rhs;
  return lhs;
}

inline mTime operator-(mTime lhs, const mTime& rhs) {
  lhs -= rhs;
  return lhs;
}

struct Edge {
  int id;
  int u, v;
  int uLabel, vLabel;
  int eLabel;
  mTime time;

  Edge(const int _id = -1, const int _u = -1, const int _v = -1,
       const int _uLabel = -1, const int _vLabel = -1, const int _eLabel = -1,
       const int _sec = -1, const int _msec = -1)
      : id(_id),
        u(_u),
        v(_v),
        uLabel(_uLabel),
        vLabel(_vLabel),
        eLabel(_eLabel),
        time(_sec, _msec) {}

  Edge(const Edge& other)
      : id(other.id),
        u(other.u),
        v(other.v),
        uLabel(other.uLabel),
        vLabel(other.vLabel),
        eLabel(other.eLabel),
        time(other.time) {}

  bool operator<(const Edge& other) const { return time < other.time; }

  bool operator<(const mTime& other) const { return time < other; }

  bool operator==(const Edge& other) const {
    if (id != other.id) return false;
    if (u != other.u) return false;
    if (v != other.v) return false;
    if (uLabel != other.uLabel) return false;
    if (vLabel != other.vLabel) return false;
    if (eLabel != other.eLabel) return false;
    if (time != other.time) return false;

    return true;
  }

  bool operator!=(const Edge& other) const { return !(*this == other); }

  Edge& operator=(const Edge& other) {
    this->id = other.id;
    this->u = other.u;
    this->v = other.v;
    this->uLabel = other.uLabel;
    this->vLabel = other.vLabel;
    this->eLabel = other.eLabel;
    this->time = other.time;

    return *this;
  }

  bool isMatch(const Edge& qEdge) const {
    if (uLabel != qEdge.uLabel) return false;
    if (vLabel != qEdge.vLabel) return false;
    if (eLabel != qEdge.eLabel) return false;

    return true;
  }

  bool isNLFMatch(const Edge& qEdge) const {
    if (uLabel != qEdge.uLabel) return false;
    if (vLabel != qEdge.vLabel) return false;
    if (eLabel != qEdge.eLabel) return false;

    return true;
  }
};

bool operator<(const mTime& lhs, const Edge& rhs) { return lhs < rhs.time; }
class DataGraph {
 public:
  size_t numVertices;
  size_t numEdges;
  vector<int> vLabels;

  vector<Edge> stream;
  vector<vector<int>> streamHash;
  typedef deque<Edge> adjListType;
  vector<vector<adjListType>> outEdges;
  vector<vector<adjListType>> inEdges;

  DataGraph(const size_t _numVertices = 0) {
    numVertices = _numVertices;
    vLabels.resize(numVertices);
    outEdges.resize(numVertices);
    inEdges.resize(numVertices);
    for (int i = 0; i < numVertices; i++) {
      outEdges[i].resize(queryEdgeHash.size());
      inEdges[i].resize(queryEdgeHash.size());
    }
  }

  void setNumVertices(size_t _numVertices) {
    numVertices = _numVertices;
    vLabels.resize(numVertices);
    outEdges.resize(numVertices);
    inEdges.resize(numVertices);
    for (int i = 0; i < numVertices; i++) {
      outEdges[i].resize(queryEdgeHash.size());
      inEdges[i].resize(queryEdgeHash.size());
    }
  }

  void setNumEdges(size_t _numEdges) {
    numEdges = _numEdges;
    streamHash.resize(numEdges);
  }

  void insertEdge(const Edge& e) {
    for (auto hash : streamHash[e.id]) {
      outEdges[e.u][hash].push_back(e);
      inEdges[e.v][hash].push_back(e);
    }
  }

  void deleteEdge(const Edge& e) {
    for (auto hash : streamHash[e.id]) {
      outEdges[e.u][hash].pop_front();
      inEdges[e.v][hash].pop_front();
    }
  }

  vector<int> numVerticesByLabel;
  vector<int> vertexIDByLabel;
  vector<vector<int>> verticesByLabel;

  int windowSize;
  mTime avgSpan;
  void calTimeSpan(int wSize);

  int* dataNLF;
  uint64_t* dataNLFBit;
  void updateNLF(Edge& e, bool insert);
  void increaseNLF(int u, int idx);
  void decreaseNLF(int u, int idx);
};

void readDataGraph(DataGraph& G, char* dataGraphPath) {
  int numVertices = 0, numEdges = 0;
  int leftVertex, rightVertex, vLabel, eLabel, stime, mtime, graphID, temp;
  char ch;
  int maxVLabel = -1;
  map<tuple<int, int, int>, int>::iterator it;

  char* inFileBufferPtr = initFile(dataGraphPath);
  while (*inFileBufferPtr) {
    ch = parseChar(&inFileBufferPtr);
    switch (ch) {
      case 'v':
        leftVertex = parseInt(&inFileBufferPtr);
        vLabel = parseInt(&inFileBufferPtr);
        maxVLabel = max(maxVLabel, vLabel);

        numVertices++;

        break;
      case 'e':
        leftVertex = parseInt(&inFileBufferPtr);
        rightVertex = parseInt(&inFileBufferPtr);
        eLabel = parseInt(&inFileBufferPtr);
        stime = parseInt(&inFileBufferPtr);
        mtime = 0;

        numEdges++;

        break;
      case 't':
        ch = parseChar(&inFileBufferPtr);  // #
        graphID = parseInt(&inFileBufferPtr);

        break;
    }
  }
  clearFile();

  G.setNumVertices(numVertices);
  G.setNumEdges(numEdges);
  numVertices = numEdges = 0;
  inFileBufferPtr = initFile(dataGraphPath);

  while (*inFileBufferPtr) {
    ch = parseChar(&inFileBufferPtr);
    switch (ch) {
      case 'v':
        leftVertex = parseInt(&inFileBufferPtr);
        vLabel = parseInt(&inFileBufferPtr);
        G.vLabels[leftVertex] = vLabel;
        maxVLabel = max(maxVLabel, vLabel);

        break;
      case 'e':
        leftVertex = parseInt(&inFileBufferPtr);
        rightVertex = parseInt(&inFileBufferPtr);
        eLabel = parseInt(&inFileBufferPtr);
        stime = parseInt(&inFileBufferPtr);
        mtime = 0;

        G.stream.emplace_back(numEdges, leftVertex, rightVertex,
                              G.vLabels[leftVertex], G.vLabels[rightVertex],
                              eLabel, stime, mtime);

        it = queryEdgeHash.find(
            make_tuple(G.vLabels[leftVertex], G.vLabels[rightVertex], eLabel));
        if (it != queryEdgeHash.end())
          G.streamHash[numEdges].push_back(it->second);

        numEdges++;

        break;
      case 't':
        ch = parseChar(&inFileBufferPtr);  // #
        graphID = parseInt(&inFileBufferPtr);

        break;
    }
  }
  clearFile();

  G.numVerticesByLabel.resize(maxVLabel + 1, 0);
  G.vertexIDByLabel.resize(G.numVertices);
  G.verticesByLabel.resize(maxVLabel + 1);

  for (int i = 0; i < G.numVertices; i++) {
    G.vertexIDByLabel[i] = G.numVerticesByLabel[G.vLabels[i]]++;
    G.verticesByLabel[G.vLabels[i]].push_back(i);
  }
}

void DataGraph::calTimeSpan(int wSize) {
  windowSize = wSize;

  long long int span = stream[stream.size() - 1].time.sec - stream[0].time.sec;
  long long int winTimes = 1 + stream.size() / windowSize;
  long long int _avgSpan = 1 + span / winTimes;
  avgSpan.sec = _avgSpan;
  avgSpan.msec = 0;
}

void DataGraph::updateNLF(Edge& e, bool insert) {
  if (insert) {
    auto it = NLFHash.find(make_tuple(vLabels[e.v], e.eLabel, 1));
    if (it != NLFHash.end()) increaseNLF(e.u, it->second);

    it = NLFHash.find(make_tuple(vLabels[e.u], e.eLabel, 0));
    if (it != NLFHash.end()) increaseNLF(e.v, it->second);
  } else {
    auto it = NLFHash.find(make_tuple(vLabels[e.v], e.eLabel, 1));
    if (it != NLFHash.end()) decreaseNLF(e.u, it->second);

    it = NLFHash.find(make_tuple(vLabels[e.u], e.eLabel, 0));
    if (it != NLFHash.end()) decreaseNLF(e.v, it->second);
  }
}

void DataGraph::increaseNLF(int u, int idx) {
  if (dataNLF[u * NLFHash.size() + idx] < 4) {
    int b = idx * 4 + dataNLF[u * NLFHash.size() + idx];
    dataNLFBit[u * NLFArraySize + (b >> 6)] |= (1ULL << (b & 0x3F));
  }
  dataNLF[u * NLFHash.size() + idx]++;
}

void DataGraph::decreaseNLF(int u, int idx) {
  dataNLF[u * NLFHash.size() + idx]--;
  if (dataNLF[u * NLFHash.size() + idx] < 4) {
    int b = idx * 4 + dataNLF[u * NLFHash.size() + idx];
    dataNLFBit[u * NLFArraySize + (b >> 6)] ^= (1ULL << (b & 0x3F));
  }
}

class QueryGraph {
 public:
  size_t numVertices;
  size_t numEdges;
  vector<int> vLabels;
  vector<int> degree;

  typedef deque<Edge> adjListType;
  vector<adjListType> outEdges;
  vector<adjListType> inEdges;
  vector<Edge> edges;
  vector<int> edgeHash;
  vector<vector<int>> beforeEdges, afterEdges;
  uint64_t** relatedEdges;

  QueryGraph(const size_t _numVertices = 0)
      : numVertices(_numVertices),
        outEdges(_numVertices),
        inEdges(_numVertices) {}

  void setNumVertices(size_t _numVertices) {
    numVertices = _numVertices;
    vLabels.resize(numVertices);
    outEdges.resize(numVertices);
    inEdges.resize(numVertices);
    degree.resize(_numVertices, 0);
  }

  void setNumEdges(size_t _numEdges) {
    numEdges = _numEdges;
    edgeHash.resize(_numEdges);
    beforeEdges.resize(_numEdges);
    afterEdges.resize(_numEdges);
    relatedEdges = new uint64_t*[_numEdges];
    for (int i = 0; i < _numEdges; i++) {
      relatedEdges[i] = new uint64_t[(_numEdges + 63) >> 6]();
    }
  }

  void insertEdge(const Edge& e) {
    outEdges[e.u].push_back(e);
    inEdges[e.v].push_back(e);
  }

  void deleteEdge(const Edge& e) {
    outEdges[e.u].pop_front();
    inEdges[e.v].pop_front();
  }

  vector<int> dagOrder;
  vector<int> dagOrderInv;
  vector<int> dagDepth;
  vector<vector<int>> dagParent;
  vector<vector<int>> dagChild;
  vector<vector<int>> dagIndex;
  vector<adjListType> dagParentOutEdges;
  vector<adjListType> dagParentInEdges;
  vector<adjListType> dagChildOutEdges;
  vector<adjListType> dagChildInEdges;
  vector<vector<bool>> dagAncestor;
  vector<vector<bool>> dagDescendant;
  vector<int> candScore;
  vector<bool> isCandidate;
  vector<int> candidateOrder;
  vector<bool> pathVisited;

  vector<vector<int>> neighbors;
  vector<vector<int>> neighborIndex;

  pair<int, vector<int>> findPath(int u);
  pair<int, int> buildDAGWithGivenRoot(int rootVertex);
  int buildDAG();

  vector<vector<int>> minPathTree, maxPathTree;
  vector<vector<int>> minPathTreeIndex, maxPathTreeIndex;
  vector<vector<int>> minPathTreePos, maxPathTreePos;
  vector<vector<vector<int>>> minPathTreeEdgeLoc, maxPathTreeEdgeLoc;
  vector<vector<int>> minPathTreeRev, maxPathTreeRev;
  vector<vector<int>> minPathTreeRevIndex, maxPathTreeRevIndex;
  vector<vector<int>> minPathTreeRevPos, maxPathTreeRevPos;
  vector<vector<vector<int>>> minPathTreeRevEdgeLoc, maxPathTreeRevEdgeLoc;

  void buildPathTree();

  int* queryNLF;
  uint64_t* queryNLFBit;
  void buildNLF();

  // bool* isLeafNode;
  // int* NECMapping;
  // int* NECInverse;
  // int* NECSize;

  // void constructLeafNEC();
};

void readQueryGraph(QueryGraph& Q, char* queryGraphPath) {
  int numVertices = 0, numEdges = 0;
  int leftVertex, rightVertex, vLabel, eLabel, graphID, temp;
  int leftEdge, rightEdge;
  char ch;

  char* inFileBufferPtr = initFile(queryGraphPath);

  while (*inFileBufferPtr) {
    ch = parseChar(&inFileBufferPtr);
    switch (ch) {
      case 'v':
        leftVertex = parseInt(&inFileBufferPtr);
        vLabel = parseInt(&inFileBufferPtr);

        numVertices++;
        break;
      case 'e':
        leftVertex = parseInt(&inFileBufferPtr);
        rightVertex = parseInt(&inFileBufferPtr);
        eLabel = parseInt(&inFileBufferPtr);

        numEdges++;
        break;
      case 'b':
        leftEdge = parseInt(&inFileBufferPtr);
        rightEdge = parseInt(&inFileBufferPtr);

        break;
      case 't':
        ch = parseChar(&inFileBufferPtr);  // #
        ch = parseChar(&inFileBufferPtr);  // s
        graphID = parseInt(&inFileBufferPtr);

        break;
    }
  }

  Q.setNumVertices(numVertices);
  Q.setNumEdges(numEdges);
  numVertices = numEdges = 0;

  clearFile();
  inFileBufferPtr = initFile(queryGraphPath);

  while (*inFileBufferPtr) {
    ch = parseChar(&inFileBufferPtr);
    switch (ch) {
      case 'v':
        leftVertex = parseInt(&inFileBufferPtr);
        Q.vLabels[leftVertex] = parseInt(&inFileBufferPtr);

        numVertices++;
        break;
      case 'e':
        leftVertex = parseInt(&inFileBufferPtr);
        rightVertex = parseInt(&inFileBufferPtr);
        eLabel = parseInt(&inFileBufferPtr);

        if (queryEdgeHash.find(make_tuple(Q.vLabels[leftVertex],
                                          Q.vLabels[rightVertex], eLabel)) ==
            queryEdgeHash.end()) {
          int size = queryEdgeHash.size();
          queryEdgeHash[make_tuple(Q.vLabels[leftVertex],
                                   Q.vLabels[rightVertex], eLabel)] = size;
        }
        if (NLFHash.find(make_tuple(Q.vLabels[rightVertex], eLabel, 1)) ==
            NLFHash.end()) {
          int size = NLFHash.size();
          NLFHash[make_tuple(Q.vLabels[rightVertex], eLabel, 1)] = size;
        }
        if (NLFHash.find(make_tuple(Q.vLabels[leftVertex], eLabel, 0)) ==
            NLFHash.end()) {
          int size = NLFHash.size();
          NLFHash[make_tuple(Q.vLabels[leftVertex], eLabel, 0)] = size;
        }

        Q.edges.emplace_back(numEdges, leftVertex, rightVertex,
                             Q.vLabels[leftVertex], Q.vLabels[rightVertex],
                             eLabel);
        Q.edgeHash[numEdges] = queryEdgeHash[make_tuple(
            Q.vLabels[leftVertex], Q.vLabels[rightVertex], eLabel)];
        Q.outEdges[leftVertex].emplace_back(numEdges, leftVertex, rightVertex,
                                            Q.vLabels[leftVertex],
                                            Q.vLabels[rightVertex], eLabel);
        Q.inEdges[rightVertex].emplace_back(numEdges, leftVertex, rightVertex,
                                            Q.vLabels[leftVertex],
                                            Q.vLabels[rightVertex], eLabel);
        Q.degree[leftVertex]++;
        Q.degree[rightVertex]++;

        numEdges++;
        break;
      case 'b':
        leftEdge = parseInt(&inFileBufferPtr);
        rightEdge = parseInt(&inFileBufferPtr);
        Q.beforeEdges[rightEdge].push_back(leftEdge);
        Q.afterEdges[leftEdge].push_back(rightEdge);
        Q.relatedEdges[rightEdge][leftEdge >> 6] |= (1ULL << (leftEdge & 0x3F));
        Q.relatedEdges[leftEdge][rightEdge >> 6] |=
            (1ULL << (rightEdge & 0x3F));

        break;
      case 't':
        ch = parseChar(&inFileBufferPtr);  // #
        ch = parseChar(&inFileBufferPtr);  // s
        graphID = parseInt(&inFileBufferPtr);

        break;
    }
  }

  clearFile();
}

pair<int, vector<int>> QueryGraph::findPath(int u) {
  int retScore = 0;
  vector<int> retPath;
  pathVisited[u] = true;
  for (auto e : outEdges[u]) {
    int v = e.v;
    if (pathVisited[v]) continue;
    int score = 0;
    vector<int> path;
    tie(score, path) = findPath(v);
    for (auto beID : beforeEdges[e.id]) {
      auto be = edges[beID];
      if (!pathVisited[be.u] || !pathVisited[be.v]) continue;
      score++;
    }
    for (auto aeID : afterEdges[e.id]) {
      auto ae = edges[aeID];
      if (!pathVisited[ae.u] || !pathVisited[ae.v]) continue;
      score++;
    }
    if (score > retScore) {
      retScore = score;
      retPath = path;
    }
  }
  for (auto e : inEdges[u]) {
    int v = e.u;
    if (pathVisited[v]) continue;
    int score = 0;
    vector<int> path;
    tie(score, path) = findPath(v);
    for (auto beID : beforeEdges[e.id]) {
      auto be = edges[beID];
      if (!pathVisited[be.u] || !pathVisited[be.v]) continue;
      score++;
    }
    for (auto aeID : afterEdges[e.id]) {
      auto ae = edges[aeID];
      if (!pathVisited[ae.u] || !pathVisited[ae.v]) continue;
      score++;
    }
    if (score > retScore) {
      retScore = score;
      retPath = path;
    }
  }

  pathVisited[u] = false;
  retPath.push_back(u);
  return make_pair(retScore, retPath);
}

// Returns score of the DAG
pair<int, int> QueryGraph::buildDAGWithGivenRoot(int rootVertex) {
  dagOrder.clear();
  dagOrderInv.clear();
  dagDepth.clear();
  dagParent.clear();
  dagChild.clear();
  dagIndex.clear();
  dagParentOutEdges.clear();
  dagParentInEdges.clear();
  dagChildOutEdges.clear();
  dagChildInEdges.clear();
  neighbors.clear();
  neighborIndex.clear();
  dagAncestor.clear();
  dagDescendant.clear();
  candScore.clear();
  isCandidate.clear();
  candidateOrder.clear();
  pathVisited.clear();

  // bfsCheck.resize(numVertices);
  dagOrder.resize(numVertices);
  dagOrderInv.resize(numVertices, -1);
  dagDepth.resize(numVertices, 0);
  dagParent.resize(numVertices);
  dagChild.resize(numVertices);
  dagIndex.resize(numVertices);
  dagParentOutEdges.resize(numVertices);
  dagParentInEdges.resize(numVertices);
  dagChildOutEdges.resize(numVertices);
  dagChildInEdges.resize(numVertices);
  dagAncestor.resize(numVertices);
  dagDescendant.resize(numVertices);
  for (int i = 0; i < numVertices; i++) {
    dagIndex[i].resize(numVertices, -1);
    dagAncestor[i].resize(numVertices, false);
    dagAncestor[i][i] = true;
    dagDescendant[i].resize(numVertices, false);
  }
  candScore.resize(numVertices, 0);
  isCandidate.resize(numVertices);
  candidateOrder.resize(numVertices, 0);
  pathVisited.resize(numVertices);

  neighbors.resize(numVertices);
  neighborIndex.resize(numVertices);
  for (int i = 0; i < numVertices; i++)
    neighborIndex[i].resize(numVertices, -1);

  int cnt = 0, ret = 0, candidateCnt = 0;
  int maxDepth = 0;
  vector<int> path;
  int pathScore;
  tie(pathScore, path) = findPath(rootVertex);
  reverse(path.begin(), path.end());

  dagOrderInv[rootVertex] = cnt;
  dagOrder[cnt++] = rootVertex;
  dagDepth[rootVertex] = 0;
  for (int i = 0; i < cnt; i++) {
    int currVertex = dagOrder[i];
    // dagAncestor[currVertex][currVertex] = true;
    isCandidate[currVertex] = false;
    maxDepth = max(maxDepth, dagDepth[currVertex]);
    for (auto it : outEdges[currVertex]) {
      int childVertex = it.v;
      if (dagOrderInv[childVertex] != -1) continue;
      if (dagIndex[currVertex][childVertex] == -1) {
        dagIndex[childVertex][currVertex] = dagParent[childVertex].size();
        dagParent[childVertex].push_back(currVertex);
        dagDepth[childVertex] =
            max(dagDepth[childVertex], dagDepth[currVertex] + 1);

        dagIndex[currVertex][childVertex] = dagChild[currVertex].size();
        dagChild[currVertex].push_back(childVertex);

        for (int j = 0; j < numVertices; j++)
          if (dagAncestor[currVertex][j]) dagAncestor[childVertex][j] = true;
      }
      dagChildOutEdges[currVertex].push_back(it);
      dagParentInEdges[childVertex].push_back(it);
      if (!isCandidate[childVertex])
        candidateOrder[childVertex] = candidateCnt++;
      isCandidate[childVertex] = true;
    }
    for (auto it : inEdges[currVertex]) {
      int childVertex = it.u;
      if (dagOrderInv[childVertex] != -1) continue;
      if (dagIndex[currVertex][childVertex] == -1) {
        dagIndex[childVertex][currVertex] = dagParent[childVertex].size();
        dagParent[childVertex].push_back(currVertex);
        dagDepth[childVertex] =
            max(dagDepth[childVertex], dagDepth[currVertex] + 1);

        dagIndex[currVertex][childVertex] = dagChild[currVertex].size();
        dagChild[currVertex].push_back(childVertex);

        for (int j = 0; j < numVertices; j++)
          if (dagAncestor[currVertex][j]) dagAncestor[childVertex][j] = true;
      }
      dagChildInEdges[currVertex].push_back(it);
      dagParentOutEdges[childVertex].push_back(it);
      if (!isCandidate[childVertex])
        candidateOrder[childVertex] = candidateCnt++;
      isCandidate[childVertex] = true;
    }

    int maxScore = -1, nextVertex = -1;
    for (int j = 0; j < numVertices; j++) {
      if (!isCandidate[j]) continue;
      candScore[j] = 0;
      for (auto it : outEdges[j]) {
        int childVertex = it.v;
        if (dagOrderInv[childVertex] != -1) continue;
        for (auto aeID : afterEdges[it.id]) {
          auto ae = edges[aeID];
          if (dagAncestor[j][ae.u] && dagAncestor[j][ae.v]) candScore[j]++;
          if (ae.u == childVertex && dagOrderInv[ae.v] == -1) candScore[j]++;
          if (ae.v == childVertex && dagOrderInv[ae.u] == -1) candScore[j]++;
        }
        for (auto beID : beforeEdges[it.id]) {
          auto be = edges[beID];
          if (dagAncestor[j][be.u] && dagAncestor[j][be.v]) candScore[j]++;
          if (be.u == childVertex && dagOrderInv[be.v] == -1) candScore[j]++;
          if (be.v == childVertex && dagOrderInv[be.u] == -1) candScore[j]++;
        }
      }
      for (auto it : inEdges[j]) {
        int childVertex = it.u;
        if (dagOrderInv[childVertex] != -1) continue;
        for (auto aeID : afterEdges[it.id]) {
          auto ae = edges[aeID];
          if (dagOrderInv[ae.u] == -1 || dagOrderInv[ae.v] == -1) continue;
          if (dagAncestor[j][ae.u] && dagAncestor[j][ae.v]) candScore[j]++;
        }
        for (auto beID : beforeEdges[it.id]) {
          auto be = edges[beID];
          if (dagOrderInv[be.u] == -1 || dagOrderInv[be.v] == -1) continue;
          if (dagAncestor[j][be.u] && dagAncestor[j][be.u]) candScore[j]++;
        }
      }
      if (maxScore < candScore[j] ||
          (maxScore == candScore[j] &&
           candidateOrder[nextVertex] > candidateOrder[j])) {
        maxScore = candScore[j];
        nextVertex = j;
      }
    }
    if (nextVertex == -1) break;
    if (i + 1 < path.size()) {
      nextVertex = path[i + 1];
      maxScore = candScore[nextVertex];
    }

    dagOrderInv[nextVertex] = cnt;
    dagOrder[cnt++] = nextVertex;
    ret += maxScore;
  }

  for (int i = 0; i < numVertices; i++)
    for (int j = 0; j < numVertices; j++)
      if (dagAncestor[i][j]) dagDescendant[j][i] = true;

  for (int i = 0; i < numVertices; i++) {
    sort(dagChildOutEdges[i].begin(), dagChildOutEdges[i].end());
    sort(dagChildInEdges[i].begin(), dagChildInEdges[i].end());
    sort(dagParentOutEdges[i].begin(), dagParentOutEdges[i].end());
    sort(dagParentInEdges[i].begin(), dagParentInEdges[i].end());
    for (auto it : outEdges[i]) {
      int neighborVertex = it.v;

      neighborIndex[i][it.v] = neighbors[i].size();
      neighbors[i].push_back(it.v);

      neighborIndex[it.v][i] = neighbors[it.v].size();
      neighbors[it.v].push_back(i);
    }
  }
  return make_pair(ret, maxDepth);
}

// Returns root vertex
int QueryGraph::buildDAG() {
  int rootVertex = -1;
  int maxScore = 0;
  int maxDepth = 0;

  for (int i = 0; i < numVertices; i++) {
    auto ret = buildDAGWithGivenRoot(i);
    int score = ret.first;
    int depth = ret.second;
    if (rootVertex == -1 || maxScore < score ||
        (maxScore == score && maxDepth < depth)) {
      rootVertex = i;
      maxScore = score;
    }
  }

  buildDAGWithGivenRoot(rootVertex);
  buildPathTree();
  return rootVertex;
}

void QueryGraph::buildPathTree() {
  minPathTree.resize(numVertices);
  maxPathTree.resize(numVertices);
  minPathTreeIndex.resize(numVertices);
  maxPathTreeIndex.resize(numVertices);
  minPathTreePos.resize(numVertices);
  maxPathTreePos.resize(numVertices);
  minPathTreeEdgeLoc.resize(numVertices);
  maxPathTreeEdgeLoc.resize(numVertices);
  minPathTreeRev.resize(numVertices);
  maxPathTreeRev.resize(numVertices);
  minPathTreeRevIndex.resize(numVertices);
  maxPathTreeRevIndex.resize(numVertices);
  minPathTreeRevPos.resize(numVertices);
  maxPathTreeRevPos.resize(numVertices);
  minPathTreeRevEdgeLoc.resize(numVertices);
  maxPathTreeRevEdgeLoc.resize(numVertices);
  for (int i = 0; i < numVertices; i++) {
    minPathTreeEdgeLoc[i].resize(numEdges);
    maxPathTreeEdgeLoc[i].resize(numEdges);
    minPathTreeRevEdgeLoc[i].resize(numEdges);
    maxPathTreeRevEdgeLoc[i].resize(numEdges);
  }

  for (int i = numVertices - 1; i >= 0; i--) {
    int u = dagOrder[i];
    for (auto uc : dagChild[u]) {
      minPathTreeIndex[u].push_back(minPathTree[u].size());

      for (int j = 0; j < minPathTree[uc].size(); j++) {
        auto eID = minPathTree[uc][j];
        for (auto aeID : afterEdges[eID]) {
          int au = edges[aeID].u, av = edges[aeID].v;
          if (dagAncestor[u][au] && dagAncestor[u][av]) {
            minPathTreeEdgeLoc[u][eID].push_back(minPathTree[u].size());
            minPathTree[u].push_back(eID);
            minPathTreePos[u].push_back(j);
            break;
          }
        }
      }
      for (auto e : dagChildInEdges[u]) {
        if (uc != e.u) continue;
        int eID = e.id;
        for (auto aeID : afterEdges[eID]) {
          int au = edges[aeID].u, av = edges[aeID].v;
          if (dagAncestor[u][au] && dagAncestor[u][av]) {
            minPathTreeEdgeLoc[u][eID].push_back(minPathTree[u].size());
            minPathTree[u].push_back(eID);
            minPathTreePos[u].push_back(-1);
            break;
          }
        }
      }
      for (auto e : dagChildOutEdges[u]) {
        if (uc != e.v) continue;
        int eID = e.id;
        for (auto aeID : afterEdges[eID]) {
          int au = edges[aeID].u, av = edges[aeID].v;
          if (dagAncestor[u][au] && dagAncestor[u][av]) {
            minPathTreeEdgeLoc[u][eID].push_back(minPathTree[u].size());
            minPathTree[u].push_back(eID);
            minPathTreePos[u].push_back(-1);
            break;
          }
        }
      }

      maxPathTreeIndex[u].push_back(maxPathTree[u].size());

      for (int j = 0; j < maxPathTree[uc].size(); j++) {
        auto eID = maxPathTree[uc][j];
        for (auto beID : beforeEdges[eID]) {
          int bu = edges[beID].u, bv = edges[beID].v;
          if (dagAncestor[u][bu] && dagAncestor[u][bv]) {
            maxPathTreeEdgeLoc[u][eID].push_back(maxPathTree[u].size());
            maxPathTree[u].push_back(eID);
            maxPathTreePos[u].push_back(j);
            break;
          }
        }
      }
      for (auto e : dagChildInEdges[u]) {
        if (uc != e.u) continue;
        int eID = e.id;
        for (auto beID : beforeEdges[eID]) {
          int bu = edges[beID].u, bv = edges[beID].v;
          if (dagAncestor[u][bu] && dagAncestor[u][bv]) {
            maxPathTreeEdgeLoc[u][eID].push_back(maxPathTree[u].size());
            maxPathTree[u].push_back(eID);
            maxPathTreePos[u].push_back(-1);
            break;
          }
        }
      }
      for (auto e : dagChildOutEdges[u]) {
        if (uc != e.v) continue;
        int eID = e.id;
        for (auto beID : beforeEdges[eID]) {
          int bu = edges[beID].u, bv = edges[beID].v;
          if (dagAncestor[u][bu] && dagAncestor[u][bv]) {
            maxPathTreeEdgeLoc[u][eID].push_back(maxPathTree[u].size());
            maxPathTree[u].push_back(eID);
            maxPathTreePos[u].push_back(-1);
            break;
          }
        }
      }
    }
    minPathTreeIndex[u].push_back(minPathTree[u].size());
    maxPathTreeIndex[u].push_back(maxPathTree[u].size());
  }

  for (int i = 0; i < numVertices; i++) {
    int u = dagOrder[i];
    for (auto up : dagParent[u]) {
      minPathTreeRevIndex[u].push_back(minPathTreeRev[u].size());

      for (int j = 0; j < minPathTreeRev[up].size(); j++) {
        auto eID = minPathTreeRev[up][j];
        for (auto aeID : afterEdges[eID]) {
          int au = edges[aeID].u, av = edges[aeID].v;
          if (dagDescendant[u][au] && dagDescendant[u][av]) {
            minPathTreeRevEdgeLoc[u][eID].push_back(minPathTreeRev[u].size());
            minPathTreeRev[u].push_back(eID);
            minPathTreeRevPos[u].push_back(j);
            break;
          }
        }
      }
      for (auto e : dagParentInEdges[u]) {
        if (up != e.u) continue;
        int eID = e.id;
        for (auto aeID : afterEdges[eID]) {
          int au = edges[aeID].u, av = edges[aeID].v;
          if (dagDescendant[u][au] && dagDescendant[u][av]) {
            minPathTreeRevEdgeLoc[u][eID].push_back(minPathTreeRev[u].size());
            minPathTreeRev[u].push_back(eID);
            minPathTreeRevPos[u].push_back(-1);
            break;
          }
        }
      }
      for (auto e : dagParentOutEdges[u]) {
        if (up != e.v) continue;
        int eID = e.id;
        for (auto aeID : afterEdges[eID]) {
          int au = edges[aeID].u, av = edges[aeID].v;
          if (dagDescendant[u][au] && dagDescendant[u][av]) {
            minPathTreeRevEdgeLoc[u][eID].push_back(minPathTreeRev[u].size());
            minPathTreeRev[u].push_back(eID);
            minPathTreeRevPos[u].push_back(-1);
            break;
          }
        }
      }

      maxPathTreeRevIndex[u].push_back(maxPathTreeRev[u].size());

      for (int j = 0; j < maxPathTreeRev[up].size(); j++) {
        auto eID = maxPathTreeRev[up][j];
        for (auto beID : beforeEdges[eID]) {
          int bu = edges[beID].u, bv = edges[beID].v;
          if (dagDescendant[u][bu] && dagDescendant[u][bv]) {
            maxPathTreeRevEdgeLoc[u][eID].push_back(maxPathTreeRev[u].size());
            maxPathTreeRev[u].push_back(eID);
            maxPathTreeRevPos[u].push_back(j);
            break;
          }
        }
      }
      for (auto e : dagParentInEdges[u]) {
        if (up != e.u) continue;
        int eID = e.id;
        for (auto beID : beforeEdges[eID]) {
          int bu = edges[beID].u, bv = edges[beID].v;
          if (dagDescendant[u][bu] && dagDescendant[u][bv]) {
            maxPathTreeRevEdgeLoc[u][eID].push_back(maxPathTreeRev[u].size());
            maxPathTreeRev[u].push_back(eID);
            maxPathTreeRevPos[u].push_back(-1);
            break;
          }
        }
      }
      for (auto e : dagParentOutEdges[u]) {
        if (up != e.v) continue;
        int eID = e.id;
        for (auto beID : beforeEdges[eID]) {
          int bu = edges[beID].u, bv = edges[beID].v;
          if (dagDescendant[u][bu] && dagDescendant[u][bv]) {
            maxPathTreeRevEdgeLoc[u][eID].push_back(maxPathTreeRev[u].size());
            maxPathTreeRev[u].push_back(eID);
            maxPathTreeRevPos[u].push_back(-1);
            break;
          }
        }
      }
    }
    minPathTreeRevIndex[u].push_back(minPathTreeRev[u].size());
    maxPathTreeRevIndex[u].push_back(maxPathTreeRev[u].size());
  }
}

void QueryGraph::buildNLF() {
  for (auto& e : edges) {
    int idx = NLFHash[make_tuple(e.vLabel, e.eLabel, 1)];
    if (queryNLF[e.u * NLFHash.size() + idx] < 4) {
      int b = idx * 4 + queryNLF[e.u * NLFHash.size() + idx];
      queryNLFBit[e.u * NLFArraySize + (b >> 6)] |= (1ULL << (b & 0x3F));
    }
    queryNLF[e.u * NLFHash.size() + idx]++;

    idx = NLFHash[make_tuple(e.uLabel, e.eLabel, 0)];
    if (queryNLF[e.v * NLFHash.size() + idx] < 4) {
      int b = idx * 4 + queryNLF[e.v * NLFHash.size() + idx];
      queryNLFBit[e.v * NLFArraySize + (b >> 6)] |= (1ULL << (b & 0x3F));
    }
    queryNLF[e.v * NLFHash.size() + idx]++;
  }
}

#endif
