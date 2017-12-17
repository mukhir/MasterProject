     #ifndef _DEFS_PAULY_
#define _DEFS_PAULY_

#include <list>
#include <queue>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <vmmlib/vector3.h>
#include <vmmlib/vector4.h>
#include <vmmlib/matrix3.h>
#include<Eigen/Eigen/Dense>


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

struct splatHierClustering
{
    // Its own index in the global splat array
    uint32_t index;
    // splat
    splatPoint sp;
    // if valid == false then the splat has been already merged and destroyed
    bool valid;
    
};

struct hCluster
{
	std::list<uint32_t> pointIndices;
	Eigen::Matrix3f cov;
	vmml::Vector3<float> centroid;
	Eigen::Vector3f eigenvalues;
	Eigen::Matrix3f eigenvectors;
	vmml::Vector3<float> v2;
	float variation;
};

class clusterComparison
{
public:
    bool operator()(const hCluster* const h1, const hCluster* const h2)
    {
        return (h1->variation > h2->variation);
    }
};



#endif

