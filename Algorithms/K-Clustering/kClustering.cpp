#include "kClustering.h"
#include <float.h>

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

kClustering::kClustering()
{}

kClustering::~kClustering()
{
    _splatsKClustering.clear();
}

void kClustering::init(const std::string filename)
{
	_readFileIntoMemory(filename);
    _findMinMax();
	
     _gridifyPoints();
    _computeNeighs();
	
	_groupPoints();
	_samplePoints();
	
	//cout << "Finished preprocessing" << endl;
}
void kClustering::kCentroids(uint32_t k)
{
	//cout<<"kCentroids"<<endl;
	//cout<<count[maxPointsIndex]<<" "<<k<<endl;
	
	cout<<count[maxPointsIndex]<<endl;
	
	if(count[maxPointsIndex]<k)
		_addCentroids(k);
	else
		_removeCentroids(k);
	
	_makeClusters();
	
	_freeSomeSpace();
	
			
}
uint32_t kClustering::combine(void* buffer, uint32_t finalCount)
{
    if(finalCount >= _allCount)
	 {
	 cout << "Contraction not needed" << endl;
	 cout << "Returning............." << endl;
	 return _allCount;
	 }
	 
	for(uint32_t i=0; i<_allCount; i++)
		if(_splatsKClustering[i].isCentroid)
			_combineSplats(i);

	
    cout << "Finished combining " << endl;
	
	return _writeOutputToBuffer(buffer);
	
}

void kClustering::_readFileIntoMemory(const std::string filename)
{
	size_t size;
	
    ifstream rfile(filename.c_str());
    rfile.seekg(0, ios::end);
    size = rfile.tellg() / sizeof(splatPoint);
    rfile.seekg(0, ios::beg);

    _allSplats = new splatPoint[size];
    _allCount = size;
	
#ifdef __LITTLE_ENDIAN__
	rfile.read((char*)&_allSplats[0], size*sizeof(splatPoint));
#else ifdef __BIG_ENDIAN__
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
#endif
	
    rfile.close();
    
    _splatsKClustering.reserve(size);
    _currPoints = size;
}

void kClustering::_findMinMax()
{
    _minX = _maxX = _allSplats[0].x;
    _minY = _maxY = _allSplats[1].y;
    _minZ = _maxZ = _allSplats[2].z;
    
    for(uint32_t i = 1; i < _allCount; i++)
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
    
    //cout << "Total points " << _allCount << endl;
    //cout << "X min and max " << _minX << " " << _maxX << endl;
   // cout << "Y min and max " << _minY << " " << _maxY << endl;
    //cout << "Z min and max " << _minZ << " " << _maxZ << endl;
}

void kClustering::_gridifyPoints()
{
	float floatK;
	floatK = (float)K;
    //_gridCount = K*K*K;
	_gridCount = floatK*floatK*floatK;
    _allGridCells = new gridCellEntropy[_gridCount];
	_filledGridCount = 0;
    
    float slopeX = float(floatK-1)/(_maxX - _minX);
    float slopeY = float(floatK-1)/(_maxY - _minY);
    float slopeZ = float(floatK-1)/(_maxZ - _minZ);
    
    int gx, gy, gz, index;
    uint32_t i;
	
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
        
        splatKClustering s;
        s.index = i;
        s.sp = _allSplats[i];
        s.gi = gx;
        s.gj = gy;
        s.gk = gz;
		s.overlap = 0.0;
		s.valid = true;
        s.level = 0;
		s.queue_again = false;
		s.centroidIndex = -1;
		s.isCentroid = 0;
        _splatsKClustering.push_back(s);
   
    }
	
	for(i = 0; i < _gridCount; i++)
		if(_allGridCells[i].pointIndices.size()>0)
			_filledGridCount++;
	
	//cout<<"_filledGridCount "<<_filledGridCount<<endl;
}

void kClustering::_computeNeighs()
{
    for(uint32_t i = 0; i < _allCount; i++)
        _computeNeighs(_splatsKClustering[i]);
}

void kClustering::_computeNeighs(splatKClustering& s)
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
                for(uint32_t p = 0; p < _allGridCells[index].pointIndices.size(); p++)
                {
                    uint32_t nsp = _allGridCells[index].pointIndices[p];
                    if( nsp == s.index)
                        continue;
					
                    if(_isNeighbor(s.sp, _allSplats[nsp]))
                        s.neigh.push_back(nsp);
                }
            }
	/* if(s.neigh.size()<0)
	 cout << s.neigh.size() << " "<<s.index<<endl;*/
}

bool kClustering::_isNeighbor(const splatPoint& sp1, const splatPoint& sp2)
{
    float dist = sqrt((sp1.x - sp2.x)*(sp1.x - sp2.x) + (sp1.y - sp2.y)*(sp1.y - sp2.y) +
					  (sp1.z - sp2.z)*(sp1.z - sp2.z));
	
    float r1r2 = sp1.rad + sp2.rad;
    
    return (r1r2 >= dist);
}

void kClustering::_groupPoints()
{
	int index, x, y;
	uint32_t max;
	
	count[0] = count[1] = count[2] = count[3] = count[4] = count[5] = count[6] = count[7] = 0;
	max = 0;
	
	for(uint32_t i = 0; i < _allCount; i++)
    {
		x = abs(int(_allSplats[i].x/_allSplats[i].rad));
		y = abs(int(_allSplats[i].y/_allSplats[i].rad));
		
        index = HashTable[x%8][y%8];
		
		count[index]++;
		
		_splatsKClustering[i].group = index;
		_splatsKClustering[i].isCentroid = 1;
    } 
	
	for(int j=0; j<8; j++)
	{
		if(max<count[j])
		{
			max = count[j];
			maxPointsIndex = j;
		}
		//cout<<count[j]<<endl;
	}
	
	if(maxPointsIndex == 0)
		dumpPointsIndex = 1;
	else
		dumpPointsIndex = 0;
		
	//cout<<maxPointsIndex<<endl;
}

void kClustering::_samplePoints()
{
	
	for(uint32_t i = 0; i < _allCount; i++)
		if(_splatsKClustering[i].group == maxPointsIndex)
			_calculateOverlap(i, 0);
	
	while(!_overlapQueueDecreasing.empty())
	{
		splatKClustering * s = _overlapQueueDecreasing.top();
		
		while(s->queue_again == true)
		{
			_overlapQueueDecreasing.pop();
			_calculateOverlap(s->index, 0);
			s->queue_again = false;
			//cout<<s->overlap<<endl;
			if(!_overlapQueueDecreasing.empty())
			{
				s = _overlapQueueDecreasing.top();
				//cout<<_overlapQueue.size()<<endl;
			}
		}
		
		if(!_overlapQueueDecreasing.empty())
		{
			s->group = dumpPointsIndex;
			s->isCentroid = 0;
			count[dumpPointsIndex]++;
			count[maxPointsIndex]--;
			
		
			for(uint32_t i = 0; i<s->neigh.size(); i++)
				if(_splatsKClustering[s->neigh[i]].group == maxPointsIndex)
					_splatsKClustering[s->neigh[i]].queue_again = true;
			
			_overlapQueueDecreasing.pop();
		}
	}
	
	int n;
	
	for(uint32_t i = 0; i < _allCount; i++)
	{
		splatKClustering& sk = _splatsKClustering[i];
		n = 0;
			
		if(sk.group != maxPointsIndex)
		{
			for(uint32_t k = 0; k<sk.neigh.size(); k++)
				if(_splatsKClustering[sk.neigh[k]].group == maxPointsIndex)
					n++;
			
			if(n == 0)
			{
				count[sk.group]--;
				sk.group = maxPointsIndex;
				sk.isCentroid = 1;
				count[maxPointsIndex]++;
			}
		}

	}
	
	//cout<<count[maxPointsIndex]<<endl;

}

void kClustering::_calculateOverlap(const uint32_t i, const int t)
{
	splatKClustering& sk = _splatsKClustering[i];
	float R, r, d;
	float _overlap = 0;
	
	R = sk.sp.rad;
	
	for(uint32_t j=0; j<sk.neigh.size(); j++)
	{
		splatKClustering& sn = _splatsKClustering[sk.neigh[j]];
		if(sn.group ==  maxPointsIndex)
		{	
			r = sn.sp.rad;
			d = sqrt((sk.sp.x - sn.sp.x)*(sk.sp.x - sn.sp.x) + (sk.sp.y - sn.sp.y)*(sk.sp.y - sn.sp.y) +
					 (sk.sp.z - sn.sp.z)*(sk.sp.z - sn.sp.z));
			_overlap += abs(PI*(R+r-d)*(R+r-d)*(d*d + 2.0f*d*r - 3.0f*r*r + 2.0f*d*R +6.0f*r*R - 3.0f*R*R)/(12.0f*d));
		}
	}
	
	sk.overlap = _overlap;
	
	if(sk.overlap>0.0f)
	{
		if(t == 0)
			_overlapQueueDecreasing.push(&_splatsKClustering[i]);
		if(t == 1)
			_overlapQueueIncreasing.push(&_splatsKClustering[i]);
	}
	
}

void kClustering::_freeSomeSpace()
{
    if(_allSplats)
        delete[] _allSplats;
	
	 if(_allGridCells)
	 delete[] _allGridCells;
}

void kClustering::_addCentroids(uint32_t k)
{
	//cout<<"_addCentroids"<<endl;
	for(uint32_t i=0; i<_allCount; i++)
	{
		splatKClustering& s = _splatsKClustering[i];
		if(s.group != maxPointsIndex)
			_calculateOverlap(i, 1);
	}
	
	while(count[maxPointsIndex] < k)
	{
		splatKClustering * s = _overlapQueueIncreasing.top();
		
		while(s->queue_again == true)
		{
			_overlapQueueIncreasing.pop();
			_calculateOverlap(s->index, 1);
			s->queue_again = false;
			
			s = _overlapQueueIncreasing.top();
				
		}
		
		count[s->group]--;
		s->group = maxPointsIndex;
		s->isCentroid = 1;
		count[maxPointsIndex]++;
		
		for(uint32_t i = 0; i<s->neigh.size(); i++)
			if(_splatsKClustering[s->neigh[i]].group == maxPointsIndex)
				_splatsKClustering[s->neigh[i]].queue_again = true;
		
		
		_overlapQueueIncreasing.pop();
	}
	
	//cout<<count[maxPointsIndex]<<endl;
}

void kClustering::_removeCentroids(uint32_t k)
{
	//cout<<"_removeCentroids"<<endl;
	
	for(uint32_t i=0; i<_allCount; i++)
	{
		splatKClustering& s = _splatsKClustering[i];
		if(s.group == maxPointsIndex)
			_calculateDeviation(i);
	}
	
	//cout<<"_deviationQueueIncreasing.size()"<<_deviationQueueIncreasing.size()<<endl;
	
	while(count[maxPointsIndex] > k)
	{
		splatKClustering * s = _deviationQueueIncreasing.top();
		
		count[s->group]--;
		s->group = dumpPointsIndex;
		s->isCentroid = 0;
		count[dumpPointsIndex]++;
		
		//cout<<count[maxPointsIndex]<<endl;
		
		_deviationQueueIncreasing.pop();
	}
	
	//cout<<count[maxPointsIndex]<<endl;
}

void kClustering::_calculateDeviation(const uint32_t i)
{
	splatKClustering& sk = _splatsKClustering[i];
	
	float cosTheta, denominator;
	sk.normalDeviation = 0;
	
	for(uint32_t j=0; j<sk.neigh.size(); j++)
	{
		splatKClustering& sn = _splatsKClustering[sk.neigh[j]];
		
		//cosTheta = A.B/AB
		
		denominator = (sqrt(sk.sp.nx * sk.sp.nx + sk.sp.ny * sk.sp.ny + sk.sp.nz * sk.sp.nz))*
						(sqrt(sn.sp.nx * sn.sp.nx + sn.sp.ny * sn.sp.ny + sn.sp.nz * sn.sp.nz));
		
		cosTheta = fabs((sk.sp.nx * sn.sp.nx + sk.sp.ny * sn.sp.ny + sk.sp.nz * sn.sp.nz)/denominator);
		
		if(cosTheta > 1.0f)
			cosTheta = 1.0f;
		
		sk.normalDeviation += fabs(acos(cosTheta));
	}
	
	if(sk.neigh.size() > 0)
	{
		sk.normalDeviation /= sk.neigh.size();
		_deviationQueueIncreasing.push(&_splatsKClustering[i]);
	}
	
}

void kClustering::_makeClusters()
{
	float d,dist;
	
	for(uint32_t i=0; i<_allCount; i++)
	{
		splatKClustering& s = _splatsKClustering[i];
		d = dist = FLT_MAX;
		
		if(s.group != maxPointsIndex)
		{
			for(uint32_t j = 0; j<s.neigh.size(); j++)
			{
				if(_splatsKClustering[s.neigh[j]].group == maxPointsIndex)
				{
					
					d = _distance(i, s.neigh[j]);
					if(d<dist)
					{
						dist = d;
						s.centroidIndex = s.neigh[j];
					}
					
				}
			}
			
			//cout<<i<<endl;
			
			if(s.centroidIndex == -1)
				_closestCentroidSearch(i);
			else
				_splatsKClustering[s.centroidIndex].clusterPointIndex.push_back(i);
		}
		//cout<<i<<endl;
		
	}
	cout<<"makeclustersdone"<<endl;
}

void kClustering::_closestCentroidSearch(const int pointIndex)
{ 
	float d,dist;
	int t = 0;
	int floatK = (float)K;
	
	splatKClustering& s = _splatsKClustering[pointIndex];
	
	d = dist = FLT_MAX;
	
	while(s.centroidIndex == -1)
	{
		int xl = s.gi - t; if(xl < 0) xl = 0;
		int yl = s.gj - t; if(yl < 0) yl = 0;
		int zl = s.gk - t; if(zl < 0) zl = 0;
		
		int xh = s.gi + t; if(xh >= floatK) xh = K-1;
		int yh = s.gj + t; if(yh >= floatK) yh = K-1;
		int zh = s.gk + t; if(zh >= floatK) zh = K-1;
		
		for(int i = xl; i <= xh; i++)
			for(int j = yl; j <= yh; j++)
				for(int k = zl; k <= zh; k++)
				{
					if(i == xl || i == xh || i == yl || i == yh || i == zl || i == yh)
					{
						int index = k * floatK * floatK + i * floatK + j;
					
						for(uint32_t p = 0; p < _allGridCells[index].pointIndices.size(); p++)
						{
							uint32_t nsp = _allGridCells[index].pointIndices[p];
							if( nsp == s.index)
								continue;
						
							if(_splatsKClustering[nsp].group == maxPointsIndex)
							{
								d = _distance(pointIndex, nsp);
								if(d<dist)
								{
									dist = d;
									s.centroidIndex = nsp;
								}
							}
						}
					}
				}
		
		if(s.centroidIndex != -1)
			_splatsKClustering[s.centroidIndex].clusterPointIndex.push_back(pointIndex);
		t++;
	}
}

float kClustering::_distance(const int i, const int j)
{
	float d;

	splatKClustering& s1 = _splatsKClustering[i];
	splatKClustering& s2 = _splatsKClustering[j];
	
	d = sqrt((s1.sp.x - s2.sp.x)*(s1.sp.x - s2.sp.x) + (s1.sp.y - s2.sp.y)*(s1.sp.y - s2.sp.y) +
				(s1.sp.z - s2.sp.z)*(s1.sp.z - s2.sp.z));
	
	return d;
}

void kClustering::_combineSplats(const uint32_t i)
{
	splatKClustering& s = _splatsKClustering[i];
	
	float rad, den;
	vmml::Vector3<float> p; //position
	vmml::Vector3<float> p1;
	vmml::Vector3<float> n; //normal
	
	std::list<int>::iterator it;
	
	rad = s.sp.rad;
	p.x = s.sp.x*rad;  p.y = s.sp.y*rad;  p.z = s.sp.z*rad;
	n.x = s.sp.nx*rad;  n.y = s.sp.ny*rad;  n.z = s.sp.nz*rad;
	den = rad;

		// New position and normal
	
	for(it = s.clusterPointIndex.begin(); it != s.clusterPointIndex.end(); it++)
	{
		p.x += _splatsKClustering[*it].sp.x * _splatsKClustering[*it].sp.rad;  
		p.y += _splatsKClustering[*it].sp.y * _splatsKClustering[*it].sp.rad;   
		p.z += _splatsKClustering[*it].sp.z * _splatsKClustering[*it].sp.rad; 
		den += _splatsKClustering[*it].sp.rad;
		
		n.x += _splatsKClustering[*it].sp.nx * _splatsKClustering[*it].sp.rad;
		n.y += _splatsKClustering[*it].sp.ny * _splatsKClustering[*it].sp.rad;
		n.z += _splatsKClustering[*it].sp.nz * _splatsKClustering[*it].sp.rad;
		
	}
	
	p /= den;
	n.normalize();
	
	s.sp.x = p.x;  s.sp.y = p.y;  s.sp.z = p.z;
	s.sp.nx = n.x;  s.sp.ny = n.y;  s.sp.nz = n.z;
		
		// New radius
		
	for(it = s.clusterPointIndex.begin(); it != s.clusterPointIndex.end(); it++)
	{
		p1.x = _splatsKClustering[*it].sp.x;  p1.y = _splatsKClustering[*it].sp.y;  p1.z = _splatsKClustering[*it].sp.z;
		rad = fmax(rad, (p-p1).length() + _splatsKClustering[*it].sp.rad);
	}
	
	s.sp.rad = rad;
		
}


uint32_t kClustering::_writeOutputToBuffer(void* buffer)
{
    splatPoint* data = (splatPoint*)(buffer);
    uint32_t count = 0;
	
	splatPoint sp;
  
	for(uint32_t i=0; i<_allCount; i++)
    {
		if(_splatsKClustering[i].group == maxPointsIndex)
        {
#ifdef __LITTLE_ENDIAN__
			data[count++] = _splatsKClustering[i].sp;
#else ifdef __BIG_ENDIAN__
			sp  = _splatsKClustering[i].sp;
			
			sp.x = reverseFloat((char*)&sp.x);
			sp.y = reverseFloat((char*)&sp.y);
			sp.z = reverseFloat((char*)&sp.z);
			sp.rad = reverseFloat((char*)&sp.rad);
			
			sp.nx = reverseFloat((char*)&sp.nx);
			sp.ny = reverseFloat((char*)&sp.ny);
			sp.nz = reverseFloat((char*)&sp.nz);
			//cout<<count<<" "<<i<<endl;
            data[count++] = sp;
#endif
			
        }
    }
	
	//cout<<count<<endl;

    return count;
}

int main (int argc, char* const argv[]) {
    
	/*	#ifdef __LITTLE_ENDIAN__
	 cout << "Little endian" << endl;
	 #else ifdef __BIG_ENDIAN__
	 cout << "Big Endian" << endl;
	 #endif*/

	int outputCount = 32000;

    splatPoint* sArr = new splatPoint[outputCount];
     
    kClustering kSimp;
    kSimp.init("armadillo.points");
	kSimp.kCentroids(outputCount);
	kSimp.combine(sArr, outputCount);
    
    std::ofstream wfile("output");
    wfile.write((char*)sArr, outputCount*sizeof(splatPoint));
    wfile.close();
    
    delete[] sArr;
    
    return 0;
}
