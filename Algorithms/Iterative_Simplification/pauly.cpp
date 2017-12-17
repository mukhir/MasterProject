#include <iostream>
#include <fstream>

#include <set>

#include "pauly.hpp"

using namespace std;

inline float reverseFloat (char *c) 
{
    float i;
    char *p = (char *)&i;
	
    p[0] = c[3];
    p[1] = c[2];
    p[2] = c[1];
    p[3] = c[0];
	
    return i;
}

PaulySimp::PaulySimp()
{}

PaulySimp::~PaulySimp()
{
    _splatsPauly.clear();
}

void PaulySimp::init(const std::string filename)
{
    _readFileIntoMemory(filename);
    _findMinMax();
    _gridifyPoints();
    _computeNeighs();
    
    _computeQ();
    _computeInitialContErr();
    _createPriorityQ();
    
    _freeSomeSpace();
    cout << "Finished preprocessing" << endl;
}

uint32_t PaulySimp::contract(void* buffer, uint32_t finalCount)
{
    if(finalCount >= _allCount)
    {
        cout << "Contraction not needed" << endl;
        cout << "Returning............." << endl;
        return _allCount;
    }
     
    _contract(finalCount);
    
    cout << "Finished contracting " << endl;
     
    return _writeOutputToBuffer(buffer);
}


void PaulySimp::_readFileIntoMemory(const std::string filename)
{
    size_t size;

    ifstream rfile(filename.c_str());
    rfile.seekg(0, ios::end);
    size = rfile.tellg() / sizeof(splatPoint);
    rfile.seekg(0, ios::beg);
	
	//size /=5;

    _allSplats = new splatPoint[size];
    _allCount = size;

    //rfile.read((char*)&_allSplats[0], size*sizeof(splatPoint));
	
	splatPoint sp;
	
	for(size_t i=0; i<size; i++)
    {
		
		rfile.read((char*)&sp, sizeof(splatPoint));
		
		sp.x = reverseFloat((char*)&sp.x);
		sp.y = reverseFloat((char*)&sp.y);
		sp.z = reverseFloat((char*)&sp.z);
		sp.rad = reverseFloat((char*)&sp.rad);
		
		sp.nx = reverseFloat((char*)&sp.nx);
		sp.ny = reverseFloat((char*)&sp.ny);
		sp.nz = reverseFloat((char*)&sp.nz);
		
		_allSplats[i] = sp;
	
    }

    rfile.close();
    
    _splatsPauly.reserve(size);
    _currPoints = size;
}

void PaulySimp::_findMinMax()
{
    _minX = _maxX = _allSplats[0].x;
    _minY = _maxY = _allSplats[1].y;
    _minZ = _maxZ = _allSplats[2].z;
    
    for(int i = 1; i < _allCount; i++)
    {
        _minX = fmin(_minX, _allSplats[i].x);
        _minY = fmin(_minY, _allSplats[i].y);
        _minZ = fmin(_minZ, _allSplats[i].z);
        
        _maxX = fmax(_maxX, _allSplats[i].x);
        _maxY = fmax(_maxY, _allSplats[i].y);
        _maxZ = fmax(_maxZ, _allSplats[i].z);
    }
    
    // Extend the min and max a bit
    // ** Never trust floating point operations completely when you need more precision **
    _minX -= 0.001 * (_maxX - _minX);
    _minY -= 0.001 * (_maxY - _minY);
    _minZ -= 0.001 * (_maxZ - _minZ);
    
    _maxX += 0.001 * (_maxX - _minX);
    _maxY += 0.001 * (_maxY - _minY);
    _maxZ += 0.001 * (_maxZ - _minZ);
    
    cout << "Total points " << _allCount << endl;
    cout << "X min and max " << _minX << " " << _maxX << endl;
    cout << "Y min and max " << _minY << " " << _maxY << endl;
    cout << "Z min and max " << _minZ << " " << _maxZ << endl;
}


void PaulySimp::_gridifyPoints()
{
	float floatK;
	floatK = (float)K;
    //_gridCount = K*K*K;
	_gridCount = floatK*floatK*floatK;
    _allGridCells = new gridCellPauly[_gridCount];
    
    float slopeX = float(floatK-1)/(_maxX - _minX);
    float slopeY = float(floatK-1)/(_maxY - _minY);
    float slopeZ = float(floatK-1)/(_maxZ - _minZ);
    
    int gx, gy, gz, index;
    int i;

    for(i = 0; i < _allCount; i++)
    {       
        gx = slopeX * (_allSplats[i].x - _minX);
        gy = slopeY * (_allSplats[i].y - _minY);
        gz = slopeZ * (_allSplats[i].z - _minZ);
    
        if(gx >= floatK || gx < 0 || gy >= floatK || gy < 0 || gz >= floatK || gz < 0)
        {
            continue;
        }   
        index = gz * floatK * floatK + gx * floatK + gy;
        _allGridCells[index].pointIndices.push_back(i);
        
        splatPauly s;
        s.index = i;
        s.sp = _allSplats[i];
        s.gi = gx;
        s.gj = gy;
        s.gk = gz;
        s.valid = true;
        s.level = 0;
		s.queue_again = false;
        _splatsPauly.push_back(s);
        //_splatsPauly[i].index = i;
        
        
        //if(i > 0)
        //    cout << _splatsPauly[i].index << " ";
    }
}

void PaulySimp::_computeNeighs()
{
    for(int i = 0; i < _allCount; i++)
        _computeNeighs(_splatsPauly[i]);
}

void PaulySimp::_computeNeighs(splatPauly& s)
{
	int floatK = (float)K;
	
    int xl = s.gi - 1; if(xl < 0) xl = 0;
    int yl = s.gj - 1; if(yl < 0) yl = 0;
    int zl = s.gk - 1; if(zl < 0) zl = 0;
 
    int xh = s.gi + 1; if(xh >= floatK) xh = K-1;
    int yh = s.gj + 1; if(yh >= floatK) yh = K-1;
    int zh = s.gk + 1; if(zh >= floatK) zh = K-1;
    
    //cout << xl << " " << xh << " * ";
    for(int i = xl; i <= xh; i++)
        for(int j = yl; j <= yh; j++)
            for(int k = zl; k <= zh; k++)
            {
                int index = k * floatK * floatK + i * floatK + j;
                
                //cout << _allGridCells[index].pointIndices.size() << " ";
                for(int p = 0; p < _allGridCells[index].pointIndices.size(); p++)
                {
                    int nsp = _allGridCells[index].pointIndices[p];
                    if( nsp == s.index)
                        continue;
                        
                    if(_isNeighbor(s.sp, _allSplats[nsp]))
                        s.neigh.push_back(nsp);
                }
            }
      /* if(s.neigh.size()<0)
	cout << s.neigh.size() << " "<<s.index<<endl;*/
}

void PaulySimp::_freeSomeSpace()
{
    if(_allSplats)
        delete[] _allSplats;
        
    if(_allGridCells)
        delete[] _allGridCells;
}

bool PaulySimp::_isNeighbor(const splatPoint& sp1, const splatPoint& sp2)
{
    float dist = sqrt((sp1.x - sp2.x)*(sp1.x - sp2.x) + (sp1.y - sp2.y)*(sp1.y - sp2.y) +
            (sp1.z - sp2.z)*(sp1.z - sp2.z));
            
    float r1r2 = sp1.rad + sp2.rad;
    
    return (r1r2 >= dist);
}



void PaulySimp::_computeQ()
{
    for(int i = 0; i < _allCount; i++)
        _splatsPauly[i].Q = _Q(_splatsPauly[i].sp, _splatsPauly[i].neigh);
}

vmml::Matrix4<float> PaulySimp::_Q(const splatPoint& sp1, const std::vector<int>& v)
{
    vmml::Matrix4<float> mat = _getMatrix(vmml::Vector4<float>(0.0, 0.0, 0.0, 0.0));
    for(int i = 0; i < v.size(); i++)
    {
        mat += _getMatrix(_fitPlane(sp1, _allSplats[v[i]]));
    }
    
    return mat;
}

// Creates matrix from the vector obtained from above
vmml::Matrix4<float> PaulySimp::_getMatrix(const vmml::Vector4<float>& v)
{
    vmml::Matrix4<float> mat;
    mat.m00 = v.x * v.x;    mat.m01 = v.x * v.y;    mat.m02 = v.x * v.z;    mat.m03 = v.x * v.w;
    mat.m10 = v.y * v.x;    mat.m11 = v.y * v.y;    mat.m12 = v.y * v.z;    mat.m13 = v.y * v.w;
    mat.m20 = v.z * v.x;    mat.m21 = v.z * v.y;    mat.m22 = v.z * v.z;    mat.m23 = v.z * v.w;
    mat.m30 = v.w * v.x;    mat.m31 = v.w * v.y;    mat.m32 = v.w * v.z;    mat.m33 = v.w * v.w;
    
    return mat;
}

// Fits a 2-point plane with the given normal 
// Stores a, b, c and d of plane
vmml::Vector4<float> PaulySimp::_fitPlane(const splatPoint& sp1, const splatPoint& sp2)
{
    vmml::Vector4<float> plane;
	vmml::Vector3<float> p1, p2, edge, vectorInPlane, planeNormal, pointNormal;
	
	p1 = vmml::Vector3<float>(sp1.x, sp1.y, sp1.z);
	p2 = vmml::Vector3<float>(sp2.x, sp2.y, sp2.z);
	pointNormal = vmml::Vector3<float>(sp1.nx, sp1.ny, sp1.nz);
	
	edge = p1 - p2;
	vectorInPlane = edge.cross(pointNormal);
	
	planeNormal = edge.cross(vectorInPlane);
	planeNormal.normalize();
	
    plane.x = planeNormal.x;
    plane.y = planeNormal.y;
    plane.z = planeNormal.z;
    plane.w = -(sp1.x * planeNormal.x + sp1.y * planeNormal.y + sp1.z * planeNormal.z);
    
    return plane;
}


void PaulySimp::_computeInitialContErr()
{
    for(int i = 0; i < _allCount; i++)
    {
        _computeContErrAllNeigh(i);
    }
}

void PaulySimp::_computeContErrAllNeigh(int i)
{
    splatPauly& sp = _splatsPauly[i];
     sp.cindex = -1;
     sp.cerr = 10e+6;
     
     for(int j = 0; j < sp.neigh.size(); j++)
     {
         if(!_splatsPauly[sp.neigh[j]].valid)
             continue;
         float err = _contractionError(sp, j);
         if(sp.cerr > err)
         {
             sp.cerr = err;
             sp.cindex = j;
         }
     }
}


float PaulySimp::_contractionError(const splatPauly& s, const int i)
{
	float err;
	
    vmml::Vector4<float> v = vmml::Vector4<float>(_combineSplats(s, i), 1.0);
    
    vmml::Matrix4<float> mat = s.Q +  _splatsPauly[s.neigh[i]].Q;
	
   /* vmml::Vector4<float> result = mat * v;
	if(v.dot(result) == 0.0)
		cout<<v.dot(result)<<endl;
	 return fabs(v.dot(result));*/
	
	err = fabs(mat.m00*v.x*v.x + 2.0*mat.m01*v.x*v.y + 2.0*mat.m02*v.x*v.z  + 2.0*mat.m03*v.x + mat.m11*v.y*v.y + 2.0*mat.m12*v.y*v.z + 2.0*mat.m13*v.y + mat.m22*v.z*v.z + 2.0*mat.m23*v.z + mat.m33);
	
	/*if(err==0.0)
		cout<<err<<endl;*/ 
	return err;
}


void PaulySimp::_createPriorityQ()
{
    for(int i = 0; i < _allCount; i++)
    {
        _pQueue.push(&_splatsPauly[i]);
    }
	
	
}

void PaulySimp::_contract(uint32_t finalCount)
{
    bool valid = true;
    while(_currPoints > finalCount && !_pQueue.empty())
    {
        // An invalidated splat would not be deleted immediately.
        // Rather the deletion is amortized so as to avoid searching
        // for it in the priority queue. Such a splat would eventually
        // be deleted when it comes to the root of heap.
        splatPauly* s = _pQueue.top();
        
        // TODO : if some index in neighbor is/was deleted, then update "cindex"
        
        if(s->valid && s->neigh.size() <= s->cindex)
        {
            _computeContErrAllNeigh(s->index);
        }
        
        if(!s->valid || s->neigh.size() == 0 || s->cindex == -1)//|| s->neigh.size() <= s->cindex)
        {
            _pQueue.pop();
            continue;
        }
		
		while(s->queue_again == true)
		{
			_pQueue.pop();
			_pQueue.push(s);
			s->queue_again = false;
			s = _pQueue.top();
			//cout<<"queue again "<<_currPoints<<endl;
		}
    
        
        //cout << s->index << " " << s->neigh.size() << " " << s->cerr << " " << s->cindex << endl;
     
  
        valid = _contractPair(s->index, s->cindex);
        //cout << s->index << " pushed " << (int)valid << " " << s->cindex << endl;
        if(valid)
        {
            // Contraction succeeded
            _currPoints--;
        }
        else
        {
            // Mark that neighbor invalid and remove it
            _splatsPauly[s->neigh[s->cindex]].valid = false;
            s->neigh.erase(s->neigh.begin() + s->cindex);
        }  
        
        _pQueue.pop();
        _pQueue.push(s);
        
        
    }
}

uint32_t PaulySimp::_writeOutputToBuffer(void* buffer)
{
    splatPoint* data = (splatPoint*)(buffer);
    uint32_t count = 0;
	
	splatPoint sp;
    
    for(int i=0; i<_allCount; i++)
    {
        if(_splatsPauly[i].valid)
        {
			sp  = _splatsPauly[i].sp;
			
			sp.x = reverseFloat((char*)&sp.x);
			sp.y = reverseFloat((char*)&sp.y);
			sp.z = reverseFloat((char*)&sp.z);
			sp.rad = reverseFloat((char*)&sp.rad);
			
			sp.nx = reverseFloat((char*)&sp.nx);
			sp.ny = reverseFloat((char*)&sp.ny);
			sp.nz = reverseFloat((char*)&sp.nz);
			
            data[count++] = sp;
        }
    }
    
    return count;
}


vmml::Vector3<float> PaulySimp::_combineSplats(const splatPauly& s, const int i)
{
    const splatPauly& s1 = _splatsPauly[s.neigh[i]];
    
    vmml::Vector3<float> v;
    v.x = s.sp.x * s.sp.rad + s1.sp.x * s1.sp.rad;
    v.y = s.sp.y * s.sp.rad + s1.sp.y * s1.sp.rad;
    v.z = s.sp.z * s.sp.rad + s1.sp.z * s1.sp.rad;

    float den = s.sp.rad + s1.sp.rad;
    v /= den;
    
    return v;
}

void PaulySimp::_combineSplatsEnd(splatPauly& s, const int i)
{
    splatPauly& ns = _splatsPauly[s.neigh[i]];
    
    vmml::Vector3<float> v = _combineSplats(s, i);
    vmml::Vector3<float> v1(s.sp.x, s.sp.y, s.sp.z);
    vmml::Vector3<float> v2(ns.sp.x, ns.sp.y, ns.sp.z);
    
    float rad = fmax((v-v1).length() + s.sp.rad, (v-v2).length() + ns.sp.rad);
    
    vmml::Vector3<float> normal(s.sp.rad * s.sp.nx + ns.sp.rad * ns.sp.nx,
                                s.sp.rad * s.sp.ny + ns.sp.rad * ns.sp.ny,
                                s.sp.rad * s.sp.nz + ns.sp.rad * ns.sp.nz);
    normal.normalize();
    
    s.sp.x = v.x; s.sp.y = v.y; s.sp.z = v.z;
    s.sp.rad = rad;
    s.sp.nx = normal.x; s.sp.ny = normal.y; s.sp.nz = normal.z; 
}


void PaulySimp::_combineNeighs(splatPauly& s, const int i)
{
    set<int> neighUnion;
    set<int>::iterator it;
    
    splatPauly& sn = _splatsPauly[s.neigh[i]];
        
    for(int j=0; j<s.neigh.size(); j++)
		neighUnion.insert(s.neigh[j]);
        
    for(int j=0; j<sn.neigh.size(); j++)
		neighUnion.insert(sn.neigh[j]); 
		
	
    it = neighUnion.find(s.neigh[i]);
    if(it != neighUnion.end())
        neighUnion.erase(it);
	
	it = neighUnion.find(s.index);
    if(it != neighUnion.end())
        neighUnion.erase(it);
    //neighUnion.erase(neighUnion.find(s.index));
        
    // Copy the neighbors from set to s.neighs
    s.neigh.clear();
    it = neighUnion.begin();
    while(it != neighUnion.end())
    {
        if(_splatsPauly[*it].valid)
            s.neigh.push_back(*it);
        it++;
    }
	
}


bool PaulySimp::_contractPair(const int i, const int j)
{
    if( i < 0 || j < 0)
        return false;

    splatPauly& s1 = _splatsPauly[i];
    splatPauly& s2 = _splatsPauly[s1.neigh[j]];
    int oldId = s1.neigh[j];
    
    set<int> neighs;
    set<int>::iterator it;
    
    if(!s1.valid || !s2.valid)
        // TODO : return an objection
        return false;
        
    _combineSplatsEnd(s1, j);
    _combineNeighs(s1, j); 

    s1.cerr = 10e+6;
    s2.valid = false;
	
	s1.Q = s1.Q + s2.Q;

   //Calculate next minimal value for error 
   // and update neighbor set        
    for(int l=0; l<s1.neigh.size(); l++)
    {
        int n = s1.neigh[l];
        splatPauly& sn = _splatsPauly[n];
        
        // All invalid neighbors have already been removed during neighbor merging
        //if(!sn.valid)
        //    continue;
            
        // Update contraction error
        float err1 = _contractionError(s1, l);
        if(err1 < s1.cerr)
        {
            s1.cindex = l;
            s1.cerr = err1;
        }  
		
    }
           
    return true;
}