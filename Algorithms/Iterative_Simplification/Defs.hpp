#ifndef _DEFS_PAULY_
#define _DEFS_PAULY_

#include <vector>
#include <vmmlib/vector4.h>
#include <vmmlib/matrix4.h>


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

struct splatPauly
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
    // Matrix Q
    vmml::Matrix4<float> Q;
    // Index of neighbor with minimum contraction error
    // This is the relative position of neigh in its own neigh list
    // and not index in the global array
    uint32_t cindex;
    // corresponding contraction error
    float cerr;
    uint8_t level;
	
	bool queue_again;
    
    void operator = (const splatPauly& s)
    {
        index = s.index;
        sp = s.sp;
        neigh = s.neigh;
        gi = s.gi;
        gj = s.gj;
        gk = s.gk;
        valid = s.valid;
        Q = s.Q;
        cindex = s.cindex;
        cerr = s.cerr;
        level = s.level;
    }
};

class splatComparison
{
public:
    bool operator()(const splatPauly* const s1, const splatPauly* const s2)
    {
        return (s1->cerr > s2->cerr);
    }
};

struct gridCellPauly
{
    std::vector< uint32_t > pointIndices;
};

#endif

