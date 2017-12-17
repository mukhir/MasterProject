#ifndef _PAULY_ITERATIVE_SIMPLIFICATION_
#define _PAULY_ITERATIVE_SIMPLIFICATION_

#include "Defs.hpp"

#include <queue>


class HierClustering
{
public:
    HierClustering();
    ~HierClustering();
        
    void init(const std::string filename);
    
    uint32_t contract(void* buffer, uint32_t finalCount);
    
private:
    
	static const int N_MAX = 40;
	float threshold;
	
	size_t _allCount;
	size_t _currPoints;
	splatHierClustering* _splatHierClustering;
	std::list<hCluster> _hCluster;
	
	std::list<int> f;
	
	
    //======================================================================
    //                  Private functions
    //======================================================================
    
    void _readFileIntoMemory(const std::string filename);
	void _initializeCluster();
	void _splitRecursive(hCluster& hc);
	
	void _combineSplats(uint32_t finalCount);
    
    uint32_t _writeOutputToBuffer(void* buffer);
    
	std::priority_queue< hCluster*, std::vector<hCluster*> , clusterComparison > _pQueueCluster;
	
};



#endif