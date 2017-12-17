#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <set>

#include "Defs.hpp"


using namespace vmml;

static const int HashTable[8][8] = {{0,3,6,1,4,7,2,5},
									{1,4,7,2,5,0,3,6},
									{2,5,0,3,6,1,4,7},
									{3,6,1,4,7,2,5,0},
									{4,7,2,5,0,3,6,1},
									{5,0,3,6,1,4,7,2},
									{6,1,4,7,2,5,0,3},
									{7,2,5,0,3,6,1,4}};

static const float PI = 3.1415f;

class kClustering
{
public:
	kClustering();
	~kClustering();
	
	void init(const std::string filename);
	void kCentroids(uint32_t k);
	uint32_t combine(void* buffer, uint32_t finalCount);
		
private:
	
	void _readFileIntoMemory(const std::string filename);
    void _findMinMax();
	
    void _gridifyPoints();
    void _computeNeighs();
	
	void _groupPoints();
	void _calculateOverlap(const uint32_t i, const int t);
	void _calculateDeviation(const uint32_t i);
	void _samplePoints();
	
	void _addCentroids(uint32_t k);
	void _removeCentroids(uint32_t k);
	void _makeClusters();
	
    void _freeSomeSpace();
	
	float _distance(const int i, const int j);
	void _closestCentroidSearch(const int i);

	inline
    void _computeNeighs(splatKClustering& s);
    
    inline
    bool _isNeighbor(const splatPoint& sp1, const splatPoint& sp2);
	
	void _combineSplats(const uint32_t i);
	
	uint32_t _writeOutputToBuffer(void* buffer);
	
	//======================================================================
    //              Private variables
    //======================================================================
   
	static const int K = 150;
		
	int maxPointsIndex;
	int dumpPointsIndex;
	uint32_t count[8];
	
	// Initial splats and their count
    splatPoint* _allSplats;
    gridCellEntropy* _allGridCells;
	group* _allGroups;
    uint32_t    _allCount;
    uint32_t    _gridCount;
	uint32_t    _filledGridCount;
    // Total points remaining after contraction
    uint32_t    _currPoints;
    
    std::vector<splatKClustering> _splatsKClustering;
    
    float   _minX, _minY, _minZ;
    float   _maxX, _maxY, _maxZ;
	
	std::priority_queue< splatKClustering*, std::vector<splatKClustering*> , overlapComparisonIncreasing > _overlapQueueDecreasing;
	std::priority_queue< splatKClustering*, std::vector<splatKClustering*> , overlapComparisonDecreasing > _overlapQueueIncreasing;
	
	std::priority_queue< splatKClustering*, std::vector<splatKClustering*> , distanceComparisonIncreasing > _distanceQueueIncreasing;
	std::priority_queue< splatKClustering*, std::vector<splatKClustering*> , deviationComparisonIncreasing > _deviationQueueIncreasing;
	
	std::priority_queue< gridCellEntropy*, std::vector<gridCellEntropy*> , gridComparison > _pQueueGrid;
	
	};	

