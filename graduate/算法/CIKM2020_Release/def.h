#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <string>
#include <unordered_map>
#include <map>
#include <queue>
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <utility>
#include <iso646.h>
#include <memory>
#include <iterator>

#include "time.h"
#include "Snap.h"

#define MAX_MEMORY 20000000

using namespace std;

typedef pair<int, int> edgePair;
typedef set<edgePair> SuperEdge;
typedef vector<edgePair> Community;
typedef int attrType;
constexpr int attrTypeFlag = 1;

struct pair_hash
{
    template<class T1, class T2>
    std::size_t operator() (const std::pair<T1, T2>& p) const
    {
        auto h1 = std::hash<T1>{}(p.first);
        auto h2 = std::hash<T2>{}(p.second);
        return h1 ^ h2;
    }
};

struct SuperNode
{
    int snID;     //supernode ID
    int k;		  //supernode trussness
    vector<edgePair> SNEdgeList; //supernode include edge pair
};

struct node {
    int id;
    int comEdgeCnt;
};