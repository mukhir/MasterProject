#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <set>

#include "Defs.hpp"



using namespace vmml;


class Entropy
{
public:
	Entropy();
	~Entropy();
	
	void init(const std::string filename);
	uint32_t combine(void* buffer, uint32_t finalCount);
		
private:
	
	void _readFileIntoMemory(const std::string filename);
    void _findMinMax();
    void _gridifyPoints();
    void _computeNeighs();
    void _freeSomeSpace();

	inline
    void _computeNeighs(splatEntropy& s);
    
    inline
    bool _isNeighbor(const splatPoint& sp1, const splatPoint& sp2);
	
	void _computeEntropy();
	float _computeSplatEntropy(const uint32_t i);
	
	void _combine(uint32_t finalCount);
	uint32_t combineSplats(const uint32_t i);
	
	uint32_t _writeOutputToBuffer(void* buffer);
	
	//======================================================================
    //              Private variables
    //======================================================================
    
	static const int K = 200;
	
	// Initial splats and their count
    splatPoint* _allSplats;
    gridCellEntropy* _allGridCells;
    uint32_t    _allCount;
    uint32_t    _gridCount;
    // Total points remaining after contraction
    uint32_t    _currPoints;
    
    std::vector<splatEntropy> _splatsEntropy;
    
    float   _minX, _minY, _minZ;
    float   _maxX, _maxY, _maxZ;
	
	std::priority_queue< gridCellEntropy*, std::vector<gridCellEntropy*> , gridComparison > _pQueueGrid;
	std::priority_queue< splatEntropy*, std::vector<splatEntropy*> , splatComparison > splatQueue;
	
	};	

