#pragma once
#include "def.h"

edgePair EdgePair(int u, int v);
void SplitInt(const string& s, vector<int>& v, const string& c);
void FindKTrussCom(PUNGraph G, int k, int q, vector<Community>& C);
void EdgeSupport(PUNGraph graph, unordered_map<edgePair, int, pair_hash>& sup, set<pair<int, edgePair>>& nonDecSup);
int TrussDecomposition(PUNGraph G, unordered_map<int, int>& nodeTrussness, unordered_map<edgePair, int, pair_hash>& edgeTrussness);
bool cmp(node a, node b);
void ShowCommunitySize(vector<pair<Community, set<int> > > R, FILE* file);
