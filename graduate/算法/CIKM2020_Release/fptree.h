#include "def.h"

typedef int addAttrType;				// 属性关键字类型


using Item = addAttrType;
using Transaction = std::set<Item>;
using TransformedPrefixPath = std::pair<std::set<Item>, uint64_t>;
using Pattern = std::pair<std::set<Item>, uint64_t>;


struct FPNode {
	const Item item;
	uint64_t frequency;
	std::shared_ptr<FPNode> node_link;
	std::shared_ptr<FPNode> parent;
	std::vector<std::shared_ptr<FPNode>> children;

	FPNode(const Item&, const std::shared_ptr<FPNode>&);
};

struct FPTree {
	std::shared_ptr<FPNode> root;
	std::map<Item, std::shared_ptr<FPNode>> header_table;
	uint64_t minimum_support_threshold;

	FPTree(const std::vector<Transaction>&, uint64_t);

	bool empty() const;
};


std::set<Pattern> fptree_growth(const FPTree&);
