#ifndef _PAULY_ITERATIVE_SIMPLIFICATION_
#define _PAULY_ITERATIVE_SIMPLIFICATION_

#include "Defs.hpp"

#include <queue>


class PaulySimp
{
public:
    PaulySimp();
    ~PaulySimp();
        
    void init(const std::string filename);
    
    uint32_t contract(void* buffer, uint32_t finalCount);
    
private:
    static const int K = 200;
    //======================================================================
    //                  Private functions
    //======================================================================
    
    void _readFileIntoMemory(const std::string filename);
    void _findMinMax();
    void _gridifyPoints();
    void _computeNeighs();
    void _freeSomeSpace();
    
    inline
    void _computeNeighs(splatPauly& s);
    
    inline
    bool _isNeighbor(const splatPoint& sp1, const splatPoint& sp2);
    
    
    void _computeQ();
    
    void _computeInitialContErr();
    
    void _createPriorityQ();
    
    void _contract(uint32_t finalCount);
    
    uint32_t _writeOutputToBuffer(void* buffer);
    
    inline
    void _computeContErrAllNeigh(int i);
    
    inline
    vmml::Vector4<float> _fitPlane(const splatPoint& sp1, const splatPoint& sp2);
    
    inline
    vmml::Matrix4<float> _getMatrix(const vmml::Vector4<float>& v);
    
    inline
    vmml::Matrix4<float> _Q(const splatPoint& sp1, const std::vector<int>& v);
    
    // Contraction error of contracting splatPauly "s" with its ith neighbor
    inline
    float _contractionError(const splatPauly& s, const int i);
    
    // Combine splatPauly "s" with its ith neighbor
    // Only new position is calculated and returned
    inline
    vmml::Vector3<float> _combineSplats(const splatPauly& s, const int i);
    
    // Combined splat in "s"
    // New position, normal and radius are calculated and stored in "s"
    inline
    void _combineSplatsEnd(splatPauly& s, const int i);
    
    // Combine neighbors of "s" with its ith neighbor
    // The combined neighbors are copied back to s.neigh
    inline
    void _combineNeighs(splatPauly& s, const int i);
    
    // Contraction of ith splat in "_splatsPauly" with its jth neighbor
    inline
    bool _contractPair(const int i, const int j);
    
    //======================================================================
    //              Private variables
    //======================================================================
    // Initial splats and their count
    splatPoint* _allSplats;
    gridCellPauly* _allGridCells;
    uint32_t    _allCount;
    uint32_t    _gridCount;
    // Total points remaining after contraction
    uint32_t    _currPoints;
    
    std::vector<splatPauly> _splatsPauly;
    
    float   _minX, _minY, _minZ;
    float   _maxX, _maxY, _maxZ;
    
    std::priority_queue< splatPauly*, std::vector<splatPauly*> , splatComparison > _pQueue;
};



#endif