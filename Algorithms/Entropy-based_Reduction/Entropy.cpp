#include "Entropy.h"
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

Entropy::Entropy()
{}

Entropy::~Entropy()
{
    _splatsEntropy.clear();
}

void Entropy::init(const std::string filename)
{
	_readFileIntoMemory(filename);
    _findMinMax();
    _gridifyPoints();
    _computeNeighs();
	
	_computeEntropy();
	
	_freeSomeSpace();
	cout << "Finished preprocessing" << endl;
}

uint32_t Entropy::combine(void* buffer, uint32_t finalCount)
{
    if(finalCount >= _allCount)
	 {
	 cout << "Contraction not needed" << endl;
	 cout << "Returning............." << endl;
	 return _allCount;
	 }
	 
	 _combine(finalCount);
	
	
    cout << "Finished combining " << endl;
	
	return _writeOutputToBuffer(buffer);
	
}

void Entropy::_readFileIntoMemory(const std::string filename)
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
    
    _splatsEntropy.reserve(size);
    _currPoints = size;
}

void Entropy::_findMinMax()
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
    
    cout << "Total points " << _allCount << endl;
    cout << "X min and max " << _minX << " " << _maxX << endl;
    cout << "Y min and max " << _minY << " " << _maxY << endl;
    cout << "Z min and max " << _minZ << " " << _maxZ << endl;
}

void Entropy::_gridifyPoints()
{
	float floatK;
	floatK = (float)K;
    //_gridCount = K*K*K;
	_gridCount = floatK*floatK*floatK;
    _allGridCells = new gridCellEntropy[_gridCount];
    
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
        
        splatEntropy s;
        s.index = i;
        s.sp = _allSplats[i];
        s.gi = gx;
        s.gj = gy;
        s.gk = gz;
        s.valid = true;
        s.level = 0;
        _splatsEntropy.push_back(s);
        //_splatsPauly[i].index = i;
        
        
        //if(i > 0)
        //    cout << _splatsPauly[i].index << " ";
    }
}

void Entropy::_computeNeighs()
{
    for(uint32_t i = 0; i < _allCount; i++)
        _computeNeighs(_splatsEntropy[i]);
}

void Entropy::_computeNeighs(splatEntropy& s)
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

void Entropy::_freeSomeSpace()
{
    if(_allSplats)
        delete[] _allSplats;
	
   /* if(_allGridCells)
        delete[] _allGridCells;*/
}

bool Entropy::_isNeighbor(const splatPoint& sp1, const splatPoint& sp2)
{
    float dist = sqrt((sp1.x - sp2.x)*(sp1.x - sp2.x) + (sp1.y - sp2.y)*(sp1.y - sp2.y) +
					  (sp1.z - sp2.z)*(sp1.z - sp2.z));
	
    float r1r2 = sp1.rad + sp2.rad;
    
    return (r1r2 >= dist);
}

void Entropy::_computeEntropy()
{
	float entropy;
	for(uint32_t i=0; i<_gridCount; i++)
	{
		_allGridCells[i].gEntropy = 0.0f;
		if(_allGridCells[i].pointIndices.size()>0)
		{
			for(uint32_t s=0; s<_allGridCells[i].pointIndices.size(); s++)
			{
				entropy = _computeSplatEntropy(_allGridCells[i].pointIndices[s]);
				_allGridCells[i].gEntropy += entropy;
				
				// only keeping points with some entropy/neighbours
				if(entropy != 0.0f)
				{
					_allGridCells[i]._pQueueSplat.push(&_splatsEntropy[_allGridCells[i].pointIndices[s]]);
					splatQueue.push(&_splatsEntropy[_allGridCells[i].pointIndices[s]]);
					_splatsEntropy[_allGridCells[i].pointIndices[s]].queue_again = false;
				}
			}
			if(_allGridCells[i].gEntropy != 0.0f)
				_pQueueGrid.push(&_allGridCells[i]);   // only pushing grids with points
		}
	}
}

float Entropy::_computeSplatEntropy(const uint32_t i)
{
	splatEntropy& se = _splatsEntropy[i];
	float cosTheta, denominator;
	
	se.sEntropy = 0.0f;
	
	for(uint32_t j=0; j<se.neigh.size(); j++)
	{
		splatEntropy& sn = _splatsEntropy[se.neigh[j]];
		if(sn.valid)
		{	
			//Calculate entropy of each point cosTheta = A.B/AB
			denominator = (sqrt(se.sp.nx * se.sp.nx + se.sp.ny * se.sp.ny + se.sp.nz * se.sp.nz))*
						(sqrt(sn.sp.nx * sn.sp.nx + sn.sp.ny * sn.sp.ny + sn.sp.nz * sn.sp.nz));
			cosTheta = fabs((se.sp.nx * sn.sp.nx + se.sp.ny * sn.sp.ny + se.sp.nz * sn.sp.nz)/denominator);
		
			if(cosTheta > 1.0f)
				cosTheta = 1.0f;
		
			se.sEntropy += (1.0f + float(se.level))/(1.0f + cosTheta);
		}
	}
	
	if(se.neigh.size() > 0)
		se.sEntropy /= se.neigh.size();
	
	//cout<<se.sEntropy<<endl;
	
	return se.sEntropy;
}

void Entropy::_combine(uint32_t finalCount)
{
	
	//while(_currPoints > 134000 && !_pQueueGrid.empty())
	while(_currPoints > finalCount && !splatQueue.empty())
	{
		splatEntropy * s = splatQueue.top();
		
		if(s->sEntropy <= 0.0f || !s->valid) //pop and invalid points
		{
			//cout<<"invalid "<<_currPoints<<endl;
			splatQueue.pop();
			continue;
		}
		//cout<<s->index<<" "<<s->sEntropy<<endl;
		
		while(s->queue_again == true)
		{
			splatQueue.pop();
			splatQueue.push(s);
			s->queue_again = false;
			s = splatQueue.top();
			//cout<<"queue again "<<_currPoints<<endl;
		}
		
		_currPoints -= combineSplats(s->index);  // not only one point, equal to neigbours size
		
		//cout<<"combined"<<endl;
		
		splatQueue.pop();
		splatQueue.push(s);
		
		//cout<<_currPoints<<endl;
	
	}
}

uint32_t Entropy::combineSplats(const uint32_t i)
{
	splatEntropy& s = _splatsEntropy[i];
	
	uint32_t points;
	
	float rad, den;
	vmml::Vector3<float> p; //position
	vmml::Vector3<float> p1;
	vmml::Vector3<float> n; //normal
	
	std::vector<int> neigh;
	
	rad = s.sp.rad;
	p.x = s.sp.x*rad;  p.y = s.sp.y*rad;  p.z = s.sp.z*rad;
	n.x = s.sp.nx*rad;  n.y = s.sp.ny*rad;  n.z = s.sp.nz*rad;
	den = rad;
	
	s.level++;
	points = 0;
	//cout<<points<<endl;

		// New position and normal
	
	for(uint32_t j=0; j<s.neigh.size(); j++)
	{
		if(_splatsEntropy[s.neigh[j]].valid)
		{
			p.x += _splatsEntropy[s.neigh[j]].sp.x * _splatsEntropy[s.neigh[j]].sp.rad;  
			p.y += _splatsEntropy[s.neigh[j]].sp.y * _splatsEntropy[s.neigh[j]].sp.rad;   
			p.z += _splatsEntropy[s.neigh[j]].sp.z * _splatsEntropy[s.neigh[j]].sp.rad; 
			den += _splatsEntropy[s.neigh[j]].sp.rad;
		
			n.x += _splatsEntropy[s.neigh[j]].sp.nx * _splatsEntropy[s.neigh[j]].sp.rad;
			n.y += _splatsEntropy[s.neigh[j]].sp.ny * _splatsEntropy[s.neigh[j]].sp.rad;
			n.z += _splatsEntropy[s.neigh[j]].sp.nz * _splatsEntropy[s.neigh[j]].sp.rad;
		
			_splatsEntropy[s.neigh[j]].valid = false;
			points++;
	
			//cout<<"position and normal"<<endl;
		}
	}
	
	p /= den;
	n.normalize();
	
	s.sp.x = p.x;  s.sp.y = p.y;  s.sp.z = p.z;
	s.sp.nx = n.x;  s.sp.ny = n.y;  s.sp.nz = n.z;
		
		// New radius
		
	for(uint32_t j=0; j<s.neigh.size(); j++)
	{
		//if(_splatsEntropy[s.neigh[j]].valid)
		//{
			p1.x = _splatsEntropy[s.neigh[j]].sp.x;  p1.y = _splatsEntropy[s.neigh[j]].sp.y;  p1.z = _splatsEntropy[s.neigh[j]].sp.z;
			rad = fmax(rad, (p-p1).length() + _splatsEntropy[s.neigh[j]].sp.rad);
		//}
	}
	s.sp.rad = rad;
	//cout<<s.sp.rad<<endl; 
	
		// New neighbours, add only valid points, dont add to itself
	
	neigh.clear();
	
	for(uint32_t j=0; j<s.neigh.size(); j++)
		for(uint32_t k=0; k<_splatsEntropy[s.neigh[j]].neigh.size(); k++)
			if((_splatsEntropy[_splatsEntropy[s.neigh[j]].neigh[k]].valid) && (_splatsEntropy[s.neigh[j]].neigh[k] != s.index))  //add only valid point and dont add to itself
			{
				neigh.push_back(_splatsEntropy[s.neigh[j]].neigh[k]);
				//cout<<"neighbour"<<neigh.size()<<endl;
			}
	
	
	s.neigh.clear();
	s.neigh.insert(s.neigh.begin(), neigh.begin(), neigh.end());
	
	//cout<<s.neigh.size()<<endl;
	
		// Update neighbours neighbour
	/*for(uint32_t j=0; j<s.neigh.size(); j++)
	{
		splatEntropy& sn = _splatsEntropy[s.neigh[j]];
		neigh.clear();
		neigh.push_back(i);  // Add itself to its neighbour
		
		for(uint32_t k=0; k< sn.neigh.size(); k++)
			if(_splatsEntropy[sn.neigh[k]].valid)
				neigh.push_back(sn.neigh[k]);
		
		sn.neigh.clear();
		sn.neigh.insert(sn.neigh.begin(), neigh.begin(), neigh.end());
		_computeSplatEntropy(s.neigh[j]);
		sn.queue_again = true;
		//cout<<"Update neighbours neighbour"<<j<<" "<<s.neigh.size()<<endl;
	}*/
		
	neigh.clear();
		// New entropy
	
	_computeSplatEntropy(i);
		
	return points;
	
}


uint32_t Entropy::_writeOutputToBuffer(void* buffer)
{
    splatPoint* data = (splatPoint*)(buffer);
    uint32_t count = 0;
	
	splatPoint sp;
  
    for(uint32_t i=0; i<_allCount; i++)
    {
       if(_splatsEntropy[i].valid)
		//if((_allSplats[i].z > 430.0) && (_allSplats[i].x<75.0) && (_allSplats[i].x> -75.0) &&(_allSplats[i].y > 0))
        {
#ifdef __LITTLE_ENDIAN__
			data[count++] = _splatsPauly[i].sp;
#else ifdef __BIG_ENDIAN__
			sp  = _splatsEntropy[i].sp;
			
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
  
    return count;
}

int main (int argc, char* const argv[]) {
    
	#ifdef __LITTLE_ENDIAN__
	 cout << "Little endian" << endl;
	 #else ifdef __BIG_ENDIAN__
	 cout << "Big Endian" << endl;
	 #endif
	
	int outputCount = 32000;
    splatPoint* sArr = new splatPoint[outputCount];
    
    Entropy entropySimp;
    entropySimp.init("armadillo.points");
	entropySimp.combine(sArr, outputCount);
    
    std::ofstream wfile("output");
    wfile.write((char*)sArr, outputCount*sizeof(splatPoint));
    wfile.close();
    
    delete[] sArr;
    
    return 0;
}
