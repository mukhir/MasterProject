#ifndef _DEFS_ENTROPY_
#define _DEFS_ENTROPY_

#include <vector>
#include <queue>
#include "vmmlib/vector3.h"

using namespace std;

struct splatPoint
{
    float x, y, z, rad;
    float nx, ny, nz;
    float r, g, b;
    
    void operator = (const splatPoint& sp)
    {
        x = sp.x; y = sp.y; z = sp.z; rad = sp.rad;
        nx = sp.nx; ny = sp.ny; nz = sp.nz;
        r = sp.r; g = sp.g; b = sp.b;
    }
    
    void operator = (float val)
    {
        x = y = z = rad = val;
        nx = ny = nz = val;
        r = g = b = val;
    }
    
    splatPoint operator * (float val)
    {
        splatPoint sp;
        sp = *this;
        sp.x *= val; sp.y *= val; sp.z *= val; sp.rad *= val;
        sp.nx *= val; sp.ny *= val; sp.nz *= val;
        sp.r *= val; sp.g *= val; sp.b *= val;
        return sp;
    }
    
    void operator *= (float val)
    {
        x *= val; y *= val; z *= val; rad *= val;
        nx *= val; ny *= val; nz *= val;
        r *= val; g *= val; b *= val;
    }
    
    void operator += (const splatPoint& sp)
    {
        x += sp.x; y += sp.y; z += sp.z; rad += sp.rad;
        nx += sp.nx; ny += sp.ny; nz += sp.nz;
        r += sp.r; g += sp.g; b += sp.b;
    }  
};

struct splatKClustering
{
    // Its own index in the global splat array
    uint32_t index;
    // splat
    splatPoint sp;
    // All neighbors(only indices)
    std::vector<int> neigh;
    // Grid indices
    int gi, gj, gk;
    // if valid == false then the splat has been already merged and destroyed
    bool valid;
    // entropy
    float overlap;
	// level 
    uint8_t level;
	// to maintain priority
	bool queue_again;
	// group where it belongs
	int group;
	// cumulative distance of points in the group around it
	float distance;
	// centroid to which this point is closest
	uint32_t centroidIndex;
	// yes for centroid
	bool isCentroid;
	// list of cluster points
	std::list<int> clusterPointIndex;
	// normal deviation with neighbours
	float normalDeviation;
    
    void operator = (const splatKClustering& s)
    {
        index = s.index;
        sp = s.sp;
        neigh = s.neigh;
        gi = s.gi;
        gj = s.gj;
        gk = s.gk;
        valid = s.valid;
        level = s.level;
    }
};

class overlapComparisonIncreasing
	{
	public:
		bool operator()(const splatKClustering* const s1, const splatKClustering* const s2)
		{
			return (s1->overlap < s2->overlap);
		}
	};

class overlapComparisonDecreasing
	{
	public:
		bool operator()(const splatKClustering* const s1, const splatKClustering* const s2)
		{
			return (s1->overlap > s2->overlap);
		}
	};

class distanceComparisonIncreasing
	{
	public:
		bool operator()(const splatKClustering* const s1, const splatKClustering* const s2)
		{
			return (s1->distance < s2->distance);
		}
	};

class deviationComparisonIncreasing
	{
	public:
		bool operator()(const splatKClustering* const s1, const splatKClustering* const s2)
		{
			return (s1->normalDeviation > s2->normalDeviation);
		}
	};

struct gridCellEntropy
{
    std::vector< uint32_t > pointIndices;
	float gEntropy;
	
	std::priority_queue< splatKClustering*, std::vector<splatKClustering*> , overlapComparisonIncreasing > _pQueueSplat;
	
};

struct group
{
	std::vector< uint32_t > pointIndices;
	int index;
};

class gridComparison
{
public:
	bool operator()(const gridCellEntropy* const g1, const gridCellEntropy* const g2)
	{
		return (g1->gEntropy > g2->gEntropy);
	}
};



#endif

