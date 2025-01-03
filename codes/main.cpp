#include <cstdlib>
#include <iostream>

#include "CS.hpp"
#include "graph.hpp"

using namespace std;

void printEdge(const Edge& e) {
  cout << e.u << "(" << e.uLabel << ")" << "->" << e.v << "(" << e.vLabel
       << "), eLabel: " << e.eLabel;
  cout << ", sec: " << e.time.sec << ", msec: " << e.time.msec;
}

int main(int argc, char* argv[]) {
  char* streamPath = argv[1];
  char* queryGraphPath = argv[2];
  int windowSize = atoi(argv[3]);

  CSStructure CS;
  readQueryGraph(CS.Q, queryGraphPath);
  readDataGraph(CS.G, streamPath);
  CS.G.calTimeSpan(windowSize);

  CS.rootVertex = CS.Q.buildDAG();
  CS.allocate();
  CS.init();

  CS.Q.buildNLF();

  int l = 0, r = 0;
  long long expiredMatches = 0;
  for (; r < CS.G.stream.size(); r++) {
    if (r + 1 == windowSize) queryProcessTimer.start();
    auto& insertEdge = CS.G.stream[r];
    // cout << "insertEdge: ";
    // printEdge(insertEdge);
    // cout << endl;

    if (CS.G.streamHash[r].size() != 0) {
      CS.G.insertEdge(insertEdge);
      CS.G.updateNLF(insertEdge, true);
      auto CSEdges = CS.findCSChangedEdges(insertEdge);
      auto NLFEdges = CS.findNLFChangedEdges(insertEdge);
      auto TCMCSEdges = CS.insertionUpdateE(CSEdges);
      CS.insertionUpdateD(TCMCSEdges, NLFEdges);
      CS.findMatches(CSEdges);
    }

    for (; l < r; l++) {
      auto& deleteEdge = CS.G.stream[l];

      if (insertEdge.time - deleteEdge.time > CS.G.avgSpan) {
        // cout << "deleteEdge: ";
        // printEdge(deleteEdge);
        // cout << endl;
        if (CS.G.streamHash[l].size() != 0) {
          auto CSEdges = CS.findCSChangedEdges(deleteEdge);
          auto NLFEdges = CS.findNLFChangedEdges(deleteEdge);
          CS.G.deleteEdge(deleteEdge);
          CS.G.updateNLF(deleteEdge, false);
          auto TCMCSEdges = CS.deletionUpdateE(CSEdges);
          CS.deletionUpdateD(TCMCSEdges, NLFEdges);
          // cout << "delete matches: " << CS.storedMatches[deleteEdge.id] <<
          // endl;
          expiredMatches += CS.storedMatches[deleteEdge.id];
        }

      } else
        break;
    }
  }
  for (; l < r; l++) {
    auto& deleteEdge = CS.G.stream[l];
    if (CS.G.streamHash[l].size() != 0) {
      auto CSEdges = CS.findCSChangedEdges(deleteEdge);
      auto NLFEdges = CS.findNLFChangedEdges(deleteEdge);
      CS.G.deleteEdge(deleteEdge);
      CS.G.updateNLF(deleteEdge, false);
      auto TCMCSEdges = CS.deletionUpdateE(CSEdges);
      CS.deletionUpdateD(TCMCSEdges, NLFEdges);
      expiredMatches += CS.storedMatches[deleteEdge.id];
    }
  }
  if (r >= windowSize) queryProcessTimer.end();

  cout << "numMatches: " << numMatches << endl;
  cout << "expiredMatches: " << expiredMatches << endl;
  cout << "query processing time: " << queryProcessTimer << "ms" << endl;

  return 0;
}