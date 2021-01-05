#pragma once
#include "def.h"

class EquiTrussAttribute
{
public:
    PUNGraph G; 
    PUNGraph superGraph; 
    vector<SuperNode> superNodes; 
    SuperEdge superEdge; 
    unordered_map<edgePair, int, pair_hash> global_edgeTrussness; 
    unordered_map<int, int> global_nodeTrussness; 
    unordered_map<int, set<int> > vertexSNIndex;

    unordered_map<int, vector<attrType> > attributesOfVertex;
    unordered_map<int, string> attr_int2str;
    unordered_map<attrType, set<int> > vertices_in_attribute;
    map<attrType, vector<int> > vertices_in_supernode;
    unordered_map<int, set<attrType> > SNAttributes;
    unordered_map<int, vector<vector<edgePair> > > attrTruss;

public:
    bool LoadSuperGraph(const char* superNode_file, const char* superEdge_file, const char* nodeTruss_file, const char* edgeTruss_file);
    void ProcessEdge(edgePair edge, int trussness, int k, unordered_map<edgePair, bool, pair_hash>& processed, queue<edgePair>& Q, unordered_map<edgePair, vector<int>, pair_hash>& SNIdentifyList, int snID);
    void BuildEquiTruss(const char* superNode_file, const char* superEdge_file, const char* nodeTruss_file, const char* edgeTruss_file);
    void BuildVertexSNIndex();
    void GetSNAttributes();
    void LoadVertexAttribute(const char* fileName);
    void BuildAttrTruss(const char* output_file);
    void LoadAttrTruss(const char* attrTruss_file); 

    map<int, set<pair<set<attrType>, uint64_t> > > FP_Growth(set<int>& nodes, int FPsupport, set<int>& psiInt, set<attrType>& SATR);

    vector<pair<Community, set<int> > > CAC_MTIndexD(int k, int q, set<attrType>& SATR);
    vector<pair<Community, set<int> > > CAC_MTIndexI(int k, const int q, set<attrType>& SATR);
    void CAC_MTIndexI_DFS(int psiAttrCnt, int sortedPsiAttr[], int index, int combinePsi[], int psiLevel, bool ifSelect, int k, int q, map<string, bool>& ifExistPsi);
};

