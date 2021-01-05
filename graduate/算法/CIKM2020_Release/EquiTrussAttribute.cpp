#include "EquiTrussAttribute.h"
#include "common.h"
#include "fptree.h"

vector<int> is_meet(MAX_MEMORY, 0);

// Load supergraph from the local equitruss index file
bool EquiTrussAttribute::LoadSuperGraph(const char* superNode_file, const char* superEdge_file, const char* nodeTruss_file, const char* edgeTruss_file)
{
	FILE* Fedge = fopen(superEdge_file, "r");
	FILE* Fnode = fopen(superNode_file, "r");
	FILE* FnodeTruss = fopen(nodeTruss_file, "r");
	FILE* FedgeTruss = fopen(edgeTruss_file, "r");

	int edgeCnt = 0;
	SuperNode sv;
	while (!feof(Fnode)) {
		int id, k;
		if (edgeCnt <= 0) {
			superNodes.push_back(sv);   
			sv.SNEdgeList.clear(); 
			fscanf(Fnode, "%d\t%d\t%d\n", &id, &k, &edgeCnt);
			sv.k = k, sv.snID = id;
		}
		else {
			int u, v;
			fscanf(Fnode, "%d\t%d\n", &u, &v);
			sv.SNEdgeList.push_back(make_pair(u, v));
			edgeCnt--;
		}
	}
	superNodes.push_back(sv);
	cout << "Load supernodes completed" << endl;
	while (!feof(Fedge)) {
		int u, v;
		fscanf(Fedge, "%d\t%d\n", &u, &v);
		superEdge.insert(make_pair(u, v));
	}
	cout << "Load superedges completed" << endl;

	while (!feof(FnodeTruss)) {
		int node, truss;
		fscanf(FnodeTruss, "%d\t%d\n", &node, &truss);
		global_nodeTrussness[node] = truss;
	}
	while (!feof(FedgeTruss)) {
		int u, v, truss;
		fscanf(FedgeTruss, "%d\t%d\t%d\n", &u, &v, &truss);
		global_edgeTrussness[make_pair(u, v)] = truss;
	}
	cout << "Load trussness completed" << endl;

	fclose(Fnode);
	fclose(Fedge);
	fclose(FedgeTruss);
	fclose(FnodeTruss);
	return true;
}

// Process the edge in building the equitruss index
void EquiTrussAttribute::ProcessEdge(edgePair edge, int trussness, int k, unordered_map<edgePair, bool, pair_hash>& processed, queue<edgePair>& Q, unordered_map<edgePair, vector<int>, pair_hash>& SNIdentifyList, int snID)
{
	if (trussness == k) {	
		if (processed[edge] == false) {
			processed[edge] = true;
			Q.push(edge);
		}
	}
	else {	
		vector<int>::iterator it = find(SNIdentifyList[edge].begin(), SNIdentifyList[edge].end(), snID);
		if (it == SNIdentifyList[edge].end()) {
			SNIdentifyList[edge].push_back(snID);
		}
	}
	return;
}

// Build the equitruss index, and save the index to files
void EquiTrussAttribute::BuildEquiTruss(const char* superNode_file, const char* superEdge_file, const char* nodeTruss_file, const char* edgeTruss_file)
{
	PUNGraph SG = TUNGraph::New();
	ifstream file;
	if (superNode_file != NULL)
	{
		file.open(superNode_file);
	}
	if (superNode_file != NULL && file)
	{
		cout << "Loading from the local index file..." << endl;
		LoadSuperGraph(superNode_file, superEdge_file, nodeTruss_file, edgeTruss_file);
		vector<SuperNode>::iterator node_it = superNodes.begin();
		node_it += 1;
		for (; node_it != superNodes.end(); node_it++) {
			SG->AddNode(node_it->snID);
		}
		set<pair<int, int>>::iterator it = superEdge.begin();
		for (; it != superEdge.end(); it++) {
			SG->AddEdge(it->first, it->second);
		}
		printf("supernodes' size: %d\n", superNodes.size() - 1);
		printf("superedges' size: %d\n", superEdge.size());
		superGraph = SG;
		BuildVertexSNIndex();
		return;
	}

	clock_t startTime = clock();
	int kMax = TrussDecomposition(G, global_nodeTrussness, global_edgeTrussness);	

	PUNGraph G_ = new TUNGraph(*G);

	int maxNodeTruss = -1, minNodeTruss = 10000, maxEdgeTruss = -1, minEdgeTruss = 10000;
	for (TUNGraph::TNodeI NI = G_->BegNI(); NI < G_->EndNI(); NI++) {
		int node = NI.GetId();
		int truss = global_nodeTrussness[node];
		if (truss > maxNodeTruss) maxNodeTruss = truss;
		if (truss < minNodeTruss) minNodeTruss = truss;
	}
	for (TUNGraph::TEdgeI EI = G_->BegEI(); EI < G_->EndEI(); EI++) {
		int u = EI.GetDstNId(), v = EI.GetSrcNId();
		if (u > v) swap(u, v);
		int truss = global_edgeTrussness[make_pair(u, v)];
		if (truss > maxEdgeTruss) maxEdgeTruss = truss;
		if (truss < minEdgeTruss) minEdgeTruss = truss;
	}


	unordered_map<edgePair, bool, pair_hash> processed;		
	unordered_map<int, set<edgePair> > sameTruss;	
	queue<edgePair> Q;				
	unordered_map<edgePair, vector<int>, pair_hash> SNIdentifyList;	

	for (TUNGraph::TEdgeI EI = G_->BegEI(); EI < G_->EndEI(); EI++) {
		edgePair e = EdgePair(EI.GetSrcNId(), EI.GetDstNId());
		processed[e] = false;
		sameTruss[global_edgeTrussness[e]].insert(e);
	}

	int snID = 0;	
	SuperNode t;
	superNodes.push_back(t); 


	for (int k = 3; k <= kMax; k++) {
		while (!sameTruss[k].empty()) {
			edgePair edge = *sameTruss[k].begin();
			processed[edge] = true;
			snID++;
			SuperNode sv;	
			sv.k = k, sv.snID = snID;
			superNodes.push_back(sv);	
			SG->AddNode(snID);

			// traverse G in BFS to identify all the edges k-truss equivalent to e by exploring its incident k-triangles
			Q.push(edge);
			while (!Q.empty()) {
				edgePair edgeUV = Q.front();
				Q.pop();
				superNodes[snID].SNEdgeList.push_back(edgeUV);
				// Create a super-edge (v, u) where u is an existing super - node with u.snID = id;
				for (vector<int>::iterator it = SNIdentifyList[edgeUV].begin(); it != SNIdentifyList[edgeUV].end(); it++) {
					int su = *it;
					superEdge.insert(EdgePair(sv.snID, su));	// Add super-edge
					if (!SG->IsEdge(sv.snID, su)) {
						SG->AddEdge(sv.snID, su);
					}
				}

				int u = edgeUV.first, v = edgeUV.second;
				int temp;
				if (G_->GetNI(u).GetDeg() > G_->GetNI(v).GetDeg()) {
					temp = u, u = v, v = temp;
				}

				for (int i = 0; i < G_->GetNI(u).GetDeg(); i++) {
					int w = G_->GetNI(u).GetNbrNId(i);	// the same neighbor of u and v
					if (!G_->IsEdge(w, v))
					{
						continue;
					}
					if (global_edgeTrussness[EdgePair(w, v)] >= k && global_edgeTrussness[EdgePair(w, u)] >= k) {
						ProcessEdge(EdgePair(u, w), global_edgeTrussness[EdgePair(u, w)], k, processed, Q, SNIdentifyList, snID);
						ProcessEdge(EdgePair(v, w), global_edgeTrussness[EdgePair(v, w)], k, processed, Q, SNIdentifyList, snID);
					}
				}
				sameTruss[k].erase(edgeUV);
				G_->DelEdge(edgeUV.first, edgeUV.second);
			}
		}
	}

	G_->Clr();

	if (superNode_file != NULL)
	{
		FILE* Fnode = fopen(superNode_file, "w");
		FILE* Fedge = fopen(superEdge_file, "w");
		FILE* FnodeTruss = fopen(nodeTruss_file, "w");
		FILE* FedgeTruss = fopen(edgeTruss_file, "w");
		printf("superNode's size: %d\n", superNodes.size() - 1);
		vector<SuperNode>::iterator node_it = superNodes.begin();
		node_it += 1;
		for (; node_it != superNodes.end(); node_it++) {
			fprintf(Fnode, "%d\t%d\t%d\n", node_it->snID, node_it->k, int(node_it->SNEdgeList.size()));
			for (edgePair ep : node_it->SNEdgeList) {
				fprintf(Fnode, "%d\t%d\n", ep.first, ep.second);
			}
		}
		cout << "Save the supernodes completed" << endl;

		printf("superEdge's size:%d\n", superEdge.size());
		set<pair<int, int>>::iterator it = superEdge.begin();
		for (; it != superEdge.end(); it++) {
			fprintf(Fedge, "%d\t%d\n", it->first, it->second);
		}
		cout << "Save the superedges completed" << endl;

		unordered_map<edgePair, int, pair_hash>::iterator edgeTruss_it = global_edgeTrussness.begin();
		for (; edgeTruss_it != global_edgeTrussness.end(); edgeTruss_it++) {
			fprintf(FedgeTruss, "%d\t%d\t%d\n", edgeTruss_it->first.first, edgeTruss_it->first.second, edgeTruss_it->second);
		}
		unordered_map<int, int>::iterator nodeTruss_it = global_nodeTrussness.begin();
		for (; nodeTruss_it != global_nodeTrussness.end(); nodeTruss_it++) {
			fprintf(FnodeTruss, "%d\t%d\n", nodeTruss_it->first, nodeTruss_it->second);
		}
		cout << "Save the trussness completed" << endl;
		fclose(Fedge);
		fclose(Fnode);
		fclose(FedgeTruss);
		fclose(FnodeTruss);
	}

	superGraph = SG;
	BuildVertexSNIndex();
	return;
}

// Build the mapping between vertices and supernodes 
void EquiTrussAttribute::BuildVertexSNIndex()
{
	if (superGraph->Empty())
	{
		return;
	}
	for (vector<SuperNode>::iterator it = superNodes.begin() + 1; it != superNodes.end(); it++) {
		for (edgePair edge : it->SNEdgeList) {
			vertexSNIndex[edge.first].insert(it->snID);
			vertexSNIndex[edge.second].insert(it->snID);
		}
	}
}

// Load the attribute set of every vertex in the graph, from the local attribute data file
void EquiTrussAttribute::LoadVertexAttribute(const char* fileName)
{
	int attr_index = 1;
	int line = 0;
	FILE* attributeF = fopen(fileName, "r");
	while (!feof(attributeF)) {
		int node;
		char attr[20005];
		fscanf(attributeF, "%d\t%s\n", &node, &attr);
		if (!G->IsNode(node)) continue;

		vector<attrType> vec;
		string str = string(attr);  
		if (str.size() == 0) continue;
		if (attrTypeFlag) SplitInt(str, vec, ",");
		for (int i = 0; i < vec.size(); i++) {
			attributesOfVertex[node].push_back(vec[i]);    
			vertices_in_attribute[vec[i]].insert(node);
			for (int sn : vertexSNIndex[node]) {
				vertices_in_supernode[vec[i]].push_back(sn);
			}
		}
	}

	for (map<int, vector<int>>::iterator iter = vertices_in_supernode.begin(); iter != vertices_in_supernode.end(); iter++)
	{
		sort(iter->second.begin(), iter->second.end());
		iter->second.erase(unique(iter->second.begin(), iter->second.end()), iter->second.end());
	}

	fclose(attributeF);
	return;
}

// Build MTIndex
void EquiTrussAttribute::BuildAttrTruss(const char* output_file)
{
	ifstream file(output_file);
	if (file)
	{
		LoadAttrTruss(output_file);
		return;
	}

	clock_t tim = clock();
	FILE* fout = fopen(output_file, "w");

	int attr_max_id = -1;
	for (auto iter = vertices_in_attribute.begin(); iter != vertices_in_attribute.end(); iter++)
	{
		if (iter->first >= attr_max_id)
		{
			attr_max_id = iter->first;
		}
	}
	cout << "attr_max_id = " << attr_max_id << endl;
	for (int i = 1; i <= attr_max_id; i++) {
		if (vertices_in_attribute[i].size() <= 2)
		{
			continue;  
		}
		PUNGraph graph = TUNGraph::New();
		TIntV nodes;
		for (int node : vertices_in_attribute[i])
		{
			nodes.Add(node);
		}
		graph = TSnap::GetSubGraph(G, nodes, false);

		clock_t startTime = clock();
		unordered_map<edgePair, int, pair_hash> edgeTrussness; 
		unordered_map<int, int> nodeTrussness;	
		int kMax = TrussDecomposition(graph, nodeTrussness, edgeTrussness);	
		fprintf(fout, "%d\t%d\n", i, -1);  

		unordered_map<edgePair, bool, pair_hash> processed;		
		unordered_map<int, set<edgePair> > sameTruss;	
		queue<edgePair> Q;				
		unordered_map<edgePair, vector<int>, pair_hash> SNIdentifyList;	

		for (TUNGraph::TEdgeI EI = graph->BegEI(); EI < graph->EndEI(); EI++) {
			edgePair e = EdgePair(EI.GetSrcNId(), EI.GetDstNId());
			processed[e] = false;
			sameTruss[edgeTrussness[e]].insert(e);
		}
		int snID = 0;

		vector<SuperNode> SNode;
		SuperEdge SEdge;
		SuperNode t;
		SNode.push_back(t); 
		

		for (int k = 3; k <= kMax; k++) {
			while (!sameTruss[k].empty()) {
				edgePair edge = *sameTruss[k].begin();
				processed[edge] = true;
				snID++;
				SuperNode sv;	
				sv.k = k, sv.snID = snID;
				SNode.push_back(sv);	

				Q.push(edge);
				while (!Q.empty()) {
					edgePair edgeUV = Q.front();
					Q.pop();
					SNode[snID].SNEdgeList.push_back(edgeUV);

					for (vector<int>::iterator it = SNIdentifyList[edgeUV].begin(); it != SNIdentifyList[edgeUV].end(); it++) {
						int su = *it;
						SEdge.insert(EdgePair(sv.snID, su));	
					}

					int u = edgeUV.first, v = edgeUV.second;
					int temp;
					if (graph->GetNI(u).GetDeg() > graph->GetNI(v).GetDeg()) {
						temp = u, u = v, v = temp;
					}

					for (int i = 0; i < graph->GetNI(u).GetDeg(); i++) {
						int w = graph->GetNI(u).GetNbrNId(i);	// the same neighbor of u and v
						if (!graph->IsEdge(w, v))
						{
							continue;
						}
						if (edgeTrussness[EdgePair(w, v)] >= k && edgeTrussness[EdgePair(w, u)] >= k) {
							ProcessEdge(EdgePair(u, w), edgeTrussness[EdgePair(u, w)], k, processed, Q, SNIdentifyList, snID);
							ProcessEdge(EdgePair(v, w), edgeTrussness[EdgePair(v, w)], k, processed, Q, SNIdentifyList, snID);
						}
					}
					sameTruss[k].erase(edgeUV);
					graph->DelEdge(edgeUV.first, edgeUV.second);
				}
			}
		}

		bool FLAG = true;
		for (SuperNode sn : SNode) {
			if (FLAG) {  
				FLAG = false;
				continue;
			}
			fprintf(fout, "%d\t%d\n", sn.k, -3);  
			for (edgePair edge : sn.SNEdgeList) {
				fprintf(fout, "%d\t%d\n", edge.first, edge.second);
			}
			fprintf(fout, "%d\t%d\n", -1, -1);
		}
		edgeTrussness.clear(), nodeTrussness.clear();
		graph->Clr();
	}
	LoadAttrTruss(output_file);
	cout << "time = " << double(clock() - tim) / CLOCKS_PER_SEC << 's' << endl << endl;
	return;
}

// Get the attributes of each supernode
void EquiTrussAttribute::GetSNAttributes()
{
	for (SuperNode sn : superNodes) {
		for (edgePair edge : sn.SNEdgeList) {
			for (attrType attr : attributesOfVertex[edge.first]) {
				SNAttributes[sn.snID].insert(attr);
			}
			for (attrType attr : attributesOfVertex[edge.second]) {
				SNAttributes[sn.snID].insert(attr);
			}
		}
	}
	return;
}

map<int, set<pair<set<attrType>, uint64_t> > > EquiTrussAttribute::FP_Growth(set<int>& nodes, int FPsupport, set<int>& psiInt, set<attrType>& SATR)
{
	vector<Transaction> transactions;
	set<int>::iterator it = nodes.begin();
	for (; it != nodes.end(); it++) {
		set<attrType> vec_to_set(attributesOfVertex[*it].begin(), attributesOfVertex[*it].end());
		set<attrType> intersection; intersection.clear();
		set_intersection(vec_to_set.begin(), vec_to_set.end(), SATR.begin(), SATR.end(), inserter(intersection, intersection.begin()));
		Transaction temp = intersection;
		transactions.push_back(temp);
	}

	const uint64_t minimum_support_threshold = FPsupport;

	const FPTree fptree{ transactions, minimum_support_threshold };

	const set<Pattern> patterns = fptree_growth(fptree);

	set<Pattern>::iterator iter = patterns.begin();
	map<int, set<Pattern> > psi;

	for (; iter != patterns.end(); iter++) 
	{
		int n = iter->first.size();
		psi[n].insert(*iter);
		psiInt.insert(n);
	}
	return psi;
}

// Load MTIndex from the local index file
void EquiTrussAttribute::LoadAttrTruss(const char* attrTruss_file) 
{
	cout << "Loading from the local attribute truss index file..." << endl;
	FILE* fin = fopen(attrTruss_file, "r");
	int attr_id, curSuperNodeK;
	vector<edgePair> curSuperNode;
	while (!feof(fin)) 
	{
		int u, v;
		fscanf(fin, "%d\t%d\n", &u, &v);
		if (u != -1 && v == -1) 
		{				
			attr_id = u;
			attrTruss[attr_id] = vector<vector<edgePair> >(2000);
		}
		else if (u == -1 && v == -1) 
		{
			for (auto edge : curSuperNode)
			{
				attrTruss[attr_id][curSuperNodeK].push_back(edge);
			}
			curSuperNode.clear();
		}
		else if (u != -1 && v == -3) 
		{
			curSuperNodeK = u;
		}
		else 
		{
			curSuperNode.push_back(EdgePair(u, v));
		}
	}
	fclose(fin);
}

// Algorithm 4: CAC-MTIndexD
// Input: A graph G (Here we use the field G in class EquiTrussAttribute), an integer k, a query vertex q, a query attribute set SATR
vector<pair<Community, set<int> > > EquiTrussAttribute::CAC_MTIndexD(int k, int q, set<attrType>& SATR)
{
	// This is the result
	// In each pair of the vector R, the first element is the community, the second element is the shared attributes of this community
	vector<pair<Community, set<int> > > R;
	set<int> NbrNodes;
	map<attrType, int> attr_nbr_cnt;
	for (attrType elem : SATR)
	{
		attr_nbr_cnt[elem] = 1;
	}
	int max_value = 0;  
	for (int i = 0; i < G->GetNI(q).GetDeg(); i++)
	{
		int w = G->GetNI(q).GetNbrNId(i);
		if (global_edgeTrussness[EdgePair(w, q)] >= k)
		{
			NbrNodes.insert(w);
			for (attrType elem : attributesOfVertex[w])
			{
				if (attr_nbr_cnt[elem] != 0)
				{
					if (attr_nbr_cnt[elem] > max_value)
					{
						max_value = attr_nbr_cnt[elem];
					}
					attr_nbr_cnt[elem] += 1;
				}
			}
		}
	}
	if (max_value < k - 1)
	{
		return R;
	}

	set<int> psiInt;
	map<int, set<pair<set<attrType>, uint64_t> > > psi;
	psi = FP_Growth(NbrNodes, k - 1, psiInt, SATR);
	if (psi.size() == 0)
	{
		return R;
	}

	unordered_map<attrType, vector<edgePair>> Hs; 
	PUNGraph G_ = TUNGraph::New();
	for (attrType attr : SATR)
	{
		auto attrSuperNodes = attrTruss[attr];
		vector<edgePair> edges;
		for (int m = k; m < attrSuperNodes.size(); m++)
		{
			for (auto edge : attrSuperNodes[m])
			{
				edges.push_back(edge);
				G_->AddEdge2(edge.first, edge.second);
			}
		}
		if (!G_->IsNode(q))
		{
			edges.clear();
		}
		else if (edges.size() >= 4000)
		{
			TIntV cnCom;
			TSnap::GetNodeWcc(G_, q, cnCom);
			G_ = TSnap::GetSubGraph(G_, cnCom);
			edges.clear();
			for (TUNGraph::TEdgeI EI = G_->BegEI(); EI != G_->EndEI(); EI++)
			{
				edges.push_back(EdgePair(EI.GetSrcNId(), EI.GetDstNId()));
			}
		}
		Hs[attr] = edges;
		G_->Clr();
	}

	vector<int> psiDex;
	for (set<int>::iterator it = psiInt.begin(); it != psiInt.end(); it++)
	{
		psiDex.push_back(*it);
	}
	int psiSize = psiDex.size();
	int i = psi.size();
	while (i >= 1)
	{
		for (set<pair<set<attrType>, uint64_t> >::iterator iter = psi[psiDex[i - 1]].begin();
			iter != psi[psiDex[i - 1]].end();
			iter++)
		{
			int j = 0;
			auto global_psi = iter->first;
			vector<edgePair> G_S;
			
			unordered_map<edgePair, int, pair_hash> edgeCount;
			for (attrType attr : global_psi)
			{
				if (j == 0)
				{
					for (auto edge : Hs[attr])
					{
						edgeCount[edge] = 1;
					}
				}
				else
				{
					for (auto edge : Hs[attr])
					{
						if (edgeCount.find(edge) != edgeCount.end())
						{
							edgeCount[edge]++;
						}
					}
				}
				j++;
			}
			for (auto iter = edgeCount.begin(); iter != edgeCount.end(); iter++)
			{
				if (iter->second == global_psi.size())
				{
					G_S.push_back(iter->first);
				}
			}
			edgeCount.clear();

			PUNGraph GRAPH = TUNGraph::New();
			for (vector<edgePair>::iterator iter = G_S.begin(); iter != G_S.end(); iter++)
			{
				int node_u = iter->first;
				int node_v = iter->second;
				GRAPH->AddEdge2(node_u, node_v);
			}
			if (GRAPH->GetNodes() == 0 || !GRAPH->IsNode(q))
			{
				continue;
			}
			G_S.clear();

			vector<Community> A;
			FindKTrussCom(GRAPH, k, q, A);
			if (A.size() != 0)
			{
				for (auto com : A) 
				{
					R.push_back(make_pair(com, global_psi));
				}
			}
			GRAPH->Clr();
		}
		if (R.size() == 0)
		{
			i--;
		}
		else 
		{
			break;
		}
	}
	return R;
}

int maxMatchPsi = 0;
vector<pair<Community, set<int> > > RRR;	
unordered_map<attrType, vector<edgePair>> attrComs;
unordered_map<attrType, set<int>> attrComVertexs;

// Algorithm 3: CAC-MTIndexI
// Input: A graph G (Here we use the field G in class EquiTrussAttribute), an integer k, a query vertex q, a query attribute set SATR
vector<pair<Community, set<int> > > EquiTrussAttribute::CAC_MTIndexI(int k, const int q, set<attrType>& SATR)
{
	clock_t origin = clock();
	vector<pair<Community, set<int> > > R;
	set<int> NbrNodes;
	unordered_map<attrType, int> attr_nbr_cnt;
	for (attrType elem : SATR) {
		attr_nbr_cnt[elem] = 1;  
	}
	int max_value = 0; 
	for (int i = 0; i < G->GetNI(q).GetDeg(); i++) {
		int w = G->GetNI(q).GetNbrNId(i);
		if (global_edgeTrussness[EdgePair(w, q)] >= k) {   
			NbrNodes.insert(w);
			for (attrType elem : attributesOfVertex[w]) {
				if (attr_nbr_cnt[elem] != 0) {
					if (attr_nbr_cnt[elem] > max_value) max_value = attr_nbr_cnt[elem];  
					attr_nbr_cnt[elem] += 1;
				}
			}
		}
	}

	if (max_value < k - 1) {
		return R;
	}

	clock_t FP_Growth_Time = clock();
	set<int> psiInt;
	map<int, set<pair<set<attrType>, uint64_t> > > psi;
	psi = FP_Growth(NbrNodes, k - 1, psiInt, SATR);
	if (psi.size() == 0) return R;    
	vector<int> psiDex;
	for (set<int>::iterator it = psiInt.begin(); it != psiInt.end(); it++) psiDex.push_back(*it);

	attrComs.clear();
	attrComVertexs.clear();
	for (attrType attr : SATR)
	{
		auto attrSuperNodes = attrTruss[attr];
		vector<edgePair> edges;
		set<int> vertexs;
		for (int m = k; m < attrSuperNodes.size(); m++)
		{
			for (auto edge : attrSuperNodes[m])
			{
				edges.push_back(edge);
				vertexs.insert(edge.first);
				vertexs.insert(edge.second);
			}
		}
		attrComs[attr] = edges;
		attrComVertexs[attr] = vertexs;
	}

	set<attrType> SATR_temp;
	map<string, bool> ifExistPsi;
	for (map<int, set<pair<set<attrType>, uint64_t> > >::iterator it1 = psi.begin(); it1 != psi.end(); it1++) {
		for (set<pair<set<attrType>, uint64_t> >::iterator it2 = it1->second.begin(); it2 != it1->second.end(); it2++) {
			string skr = "";
			set<int> tmp;
			for (attrType elem : it2->first) {
				SATR_temp.insert(elem);
				tmp.insert(elem);
			}
			for (int attr : tmp) {
				skr += to_string(attr) + "_";
			}
			ifExistPsi[skr] = true;
		}
	}
	node attr[1005];
	int psiAttrCnt = 0;
	for (attrType at : SATR_temp) {
		attr[psiAttrCnt].id = at;
		attr[psiAttrCnt++].comEdgeCnt = attrComs[at].size();   
	}
	sort(attr, attr + psiAttrCnt, cmp); 
	int sortedPsiAttr[1005];
	for (int i = 0; i < psiAttrCnt; i++) {
		sortedPsiAttr[i] = attr[i].id;
	}
	RRR.clear();
	maxMatchPsi = 0; 
	int combinePsi[1005];
	CAC_MTIndexI_DFS(psiAttrCnt, sortedPsiAttr, 0, combinePsi, 0, false, k, q, ifExistPsi);
	return RRR;
}

void EquiTrussAttribute::CAC_MTIndexI_DFS(int psiAttrCnt, int sortedPsiAttr[], int index, int combinePsi[], int psiLevel, bool ifSelect, int k, int q, map<string, bool>& ifExistPsi)
{
	if (psiLevel >= maxMatchPsi && psiLevel && ifSelect) {   			
		PUNGraph subgraph = TUNGraph::New();
		if (attrComs[combinePsi[0]].size() == 0) return;

		string skr = "";
		set<int> tmp;
		for (int i = 0; i < psiLevel; i++) tmp.insert(combinePsi[i]);
		for (int attr : tmp) {
			skr += to_string(attr) + "_";
		}
		if (!ifExistPsi[skr]) return;  

		for (edgePair edge : attrComs[combinePsi[0]]) {
			int u = edge.first, v = edge.second;
			if (!subgraph->IsNode(u)) subgraph->AddNode(u);
			if (!subgraph->IsNode(v)) subgraph->AddNode(v);
			subgraph->AddEdge(u, v);
		}

		for (int x = 1; x < psiLevel; x++) {
			if (attrComs[combinePsi[x]].size() == 0) return;
			is_meet.assign(MAX_MEMORY, 0);
			for (int ver : attrComVertexs[combinePsi[x]]) {
				is_meet[ver] = 1;
			}
			for (TUNGraph::TNodeI NI = subgraph->BegNI(); NI < subgraph->EndNI(); NI++) {
				int node = NI.GetId();
				if (is_meet[node] != 1) subgraph->DelNode(node);
			}
			if (!subgraph->IsNode(q)) return;  
		}
		for (TUNGraph::TNodeI NI = subgraph->BegNI(); NI < subgraph->EndNI(); NI++) {
			int node = NI.GetId();
			if (NI.GetDeg() < k - 1) subgraph->DelNode(node);
		}
		if (!subgraph->IsNode(q)) return;   

		queue<int> Q;
		Q.push(q);
		PUNGraph GRAPH_q = TUNGraph::New();
		GRAPH_q->AddNode(q);
		while (!Q.empty()) {
			int cur = Q.front();  Q.pop();
			for (int i = 0; i < subgraph->GetNI(cur).GetDeg(); i++) {
				int nbr = subgraph->GetNI(cur).GetNbrNId(i);
				if (!GRAPH_q->IsNode(nbr)) {
					Q.push(nbr);
					GRAPH_q->AddNode(nbr);
				}
				GRAPH_q->AddEdge(cur, nbr);
			}
		}

		vector<Community> A;
		FindKTrussCom(GRAPH_q, k, q, A);
		if (A.size() != 0) {		
			set<int> attrs;
			for (int str = 0; str < psiLevel; str++) {
				attrs.insert(combinePsi[str]);
			}
			if (psiLevel > maxMatchPsi) RRR.clear();    
			for (Community com : A)
			{
				RRR.push_back(make_pair(com, attrs));
			}
			maxMatchPsi = psiLevel;
			A.clear();
		}
	}
	if (index >= psiAttrCnt) return;
	combinePsi[psiLevel] = sortedPsiAttr[index];
	CAC_MTIndexI_DFS(psiAttrCnt, sortedPsiAttr, index + 1, combinePsi, psiLevel + 1, true, k, q, ifExistPsi); 
	CAC_MTIndexI_DFS(psiAttrCnt, sortedPsiAttr, index + 1, combinePsi, psiLevel, false, k, q, ifExistPsi);
	return;
}
