

#include <set>

#include "HierarchicalClustering.hpp"


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

HierClustering::HierClustering()
{}

HierClustering::~HierClustering()
{
    delete[] _splatHierClustering;
	_hCluster.clear();
}

void HierClustering::init(const std::string filename)
{
    _readFileIntoMemory(filename);
	
	threshold = 0.16;
	
	std::list<hCluster>::iterator it = _hCluster.begin();
	
	for(; it!=_hCluster.end(); it++)
		_splitRecursive(*it);
	
	cout <<_hCluster.size()<<" "<<_pQueueCluster.size()<<endl;
	
    cout << "Finished preprocessing" << endl;
}

uint32_t HierClustering::contract(void* buffer, uint32_t finalCount)
{
    if(finalCount >= _allCount)
    {
        cout << "Contraction not needed" << endl;
        cout << "Returning............." << endl;
        return _allCount;
    }
   
	_combineSplats(finalCount);
    
    cout << "Finished contracting " << endl;
     
    return _writeOutputToBuffer(buffer);
}


void HierClustering::_readFileIntoMemory(const std::string filename)
{
    size_t size;

    ifstream rfile(filename.c_str());
    rfile.seekg(0, ios::end);
    size = rfile.tellg() / sizeof(splatPoint);
    rfile.seekg(0, ios::beg);
	
	//size /=5;

    _allCount = size;
	_currPoints = _allCount;

    //rfile.read((char*)&_allSplats[0], size*sizeof(splatPoint));
	
	_splatHierClustering = new splatHierClustering[size];
	
	hCluster hc;
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
		
		_splatHierClustering[i].sp = sp;
		_splatHierClustering[i].index = i;
		_splatHierClustering[i].valid = true;
		
		hc.pointIndices.push_back(i);
	
    }
	
	_hCluster.push_back(hc);

    rfile.close();

}

void HierClustering::_splitRecursive(hCluster& hc)
{
	char ch;
	
	hc.centroid.x = .0f;
	hc.centroid.y = .0f;
	hc.centroid.z = .0f;
	
	//cout<<hc.pointIndices.size()<<endl;
	
	std::list<uint32_t>::iterator it = hc.pointIndices.begin();
	
	for(; it!=hc.pointIndices.end(); it++)
	{
		hc.centroid.x += _splatHierClustering[*it].sp.x;
		hc.centroid.y += _splatHierClustering[*it].sp.y;
		hc.centroid.z += _splatHierClustering[*it].sp.z;
	}
	
	hc.centroid /= (float)hc.pointIndices.size();
	
	//cout<<"centroids "<<hc.centroid.x<<" "<<hc.centroid.y<<" "<<hc.centroid.z<<endl;
	
	it = hc.pointIndices.begin();
	
	for(; it!=hc.pointIndices.end(); it++)
	{
		hc.cov(0,0) += ((_splatHierClustering[*it].sp.x - hc.centroid.x)*(_splatHierClustering[*it].sp.x - hc.centroid.x));
		hc.cov(0,1) += ((_splatHierClustering[*it].sp.x - hc.centroid.x)*(_splatHierClustering[*it].sp.y - hc.centroid.y));
		hc.cov(0,2) += ((_splatHierClustering[*it].sp.x - hc.centroid.x)*(_splatHierClustering[*it].sp.z - hc.centroid.z));
		
		hc.cov(1,1) += ((_splatHierClustering[*it].sp.y - hc.centroid.y)*(_splatHierClustering[*it].sp.y - hc.centroid.y));
		hc.cov(1,2) += ((_splatHierClustering[*it].sp.y - hc.centroid.y)*(_splatHierClustering[*it].sp.z - hc.centroid.z));
		
		hc.cov(2,2) += ((_splatHierClustering[*it].sp.z - hc.centroid.z)*(_splatHierClustering[*it].sp.z - hc.centroid.z));
		
	}
	
	hc.cov(1,0) = hc.cov(0,1);
	hc.cov(2,0) = hc.cov(0,2);
	hc.cov(2,1) = hc.cov(1,2);
	
	
	/*cout<<endl<<hc.cov(0,0)<<" "<<hc.cov(0,1)<<" "<<hc.cov(0,2)<<endl;
	cout<<hc.cov(1,0)<<" "<<hc.cov(1,1)<<" "<<hc.cov(1,2)<<endl;
	cout<<hc.cov(2,0)<<" "<<hc.cov(2,1)<<" "<<hc.cov(2,2)<<endl<<endl;*/
	
	
	
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eigensolver(hc.cov);
	
	hc.eigenvalues = eigensolver.eigenvalues();
	hc.eigenvectors = eigensolver.eigenvectors();
	
	//cout<<"values "<<hc.eigenvalues<<endl;
	
	//cout<<"vectors "<<hc.eigenvectors<<endl;

	
	hc.variation = fabs(hc.eigenvalues(0))/(fabs(hc.eigenvalues(0)) + fabs(hc.eigenvalues(1)) + fabs(hc.eigenvalues(2)));
	 
	if((hc.pointIndices.size()>N_MAX) || (hc.variation>threshold))
	{
		vmml::Vector4<float> plane;
		
		plane.x = hc.eigenvectors(0, 2);
		plane.y = hc.eigenvectors(1, 2);
		plane.z = hc.eigenvectors(2, 2);
		plane.w = -(hc.centroid.x * plane.x + hc.centroid.y * plane.y + hc.centroid.z * plane.z);
		
		//cout<<"plane "<<plane.x<<" "<<plane.y<<" "<<plane.z<<" "<<plane.w<<endl;
		
		
		float d;
		
		it = hc.pointIndices.begin();
		
		hCluster hc1, hc2;
		hc1.pointIndices.clear();
		hc2.pointIndices.clear();
		
		
		for(; it!=hc.pointIndices.end(); it++)
		{
			d = _splatHierClustering[*it].sp.x*plane.x + _splatHierClustering[*it].sp.y*plane.y + _splatHierClustering[*it].sp.z*plane.z + plane.w;
			
			if(d> .0f)
				hc1.pointIndices.push_back(*it);
			else
				hc2.pointIndices.push_back(*it);
		
		}
		
		//_hCluster.erase(*hc);
		
		if((hc1.pointIndices.size()>0) && (hc2.pointIndices.size()>0))
		{
			_hCluster.push_back(hc1);
			_hCluster.push_back(hc2);
			
			hc.pointIndices.clear();
		}
		
		/*if((hc1.pointIndices.size() ==0) || (hc2.pointIndices.size() ==0) )
			cout<<hc1.pointIndices.size()<<" "<<hc2.pointIndices.size()<<endl;*/
		
		//cin>>ch;
		
	}
	
	if(hc.pointIndices.size()>1)
		_pQueueCluster.push(& hc);
	
}

void HierClustering::_combineSplats(uint32_t finalCount)
{
	float rad, den;
	vmml::Vector3<float> p;
	vmml::Vector3<float> p1;
	vmml::Vector3<float> n;
	
	//cout<<_pQueueCluster.size()<<endl;
	
	while(_currPoints > finalCount && !_pQueueCluster.empty())
	{
		hCluster* h = _pQueueCluster.top();
		/*cout<<h->pointIndices.size()<<" ";
		cout<<h->variation<<endl;*/
		
		p.x = p.z = p.y = rad = den = .0f;
		
		std::list<uint32_t>::iterator it = h->pointIndices.begin();
		
		for(; it!=h->pointIndices.end(); it++)
		{
			p.x += (_splatHierClustering[*it].sp.x*_splatHierClustering[*it].sp.rad);
			p.y += (_splatHierClustering[*it].sp.y*_splatHierClustering[*it].sp.rad);
			p.z += (_splatHierClustering[*it].sp.z*_splatHierClustering[*it].sp.rad);
			den += _splatHierClustering[*it].sp.rad;
			
			n.x += (_splatHierClustering[*it].sp.nx*_splatHierClustering[*it].sp.rad);
			n.y += (_splatHierClustering[*it].sp.ny*_splatHierClustering[*it].sp.rad);
			n.z += (_splatHierClustering[*it].sp.nz*_splatHierClustering[*it].sp.rad);
			
			_splatHierClustering[*it].valid = false;
		}
		
		p /= den;
		n.normalize();
		
		
		for(it = h->pointIndices.begin(); it!=h->pointIndices.end(); it++)
		{
			p1.x = _splatHierClustering[*it].sp.x; p1.y = _splatHierClustering[*it].sp.y; p1.z = _splatHierClustering[*it].sp.z; 
			rad = fmax(rad, (p-p1).length() + _splatHierClustering[*it].sp.rad);
		}
		
		
		it = h->pointIndices.begin();
		
		_splatHierClustering[*it].sp.x = p.x;
		_splatHierClustering[*it].sp.y = p.y;
		_splatHierClustering[*it].sp.z = p.z;
		
		_splatHierClustering[*it].sp.nx = n.x;
		_splatHierClustering[*it].sp.ny = n.y;
		_splatHierClustering[*it].sp.nz = n.z;
		
		_splatHierClustering[*it].sp.rad = rad;
		_splatHierClustering[*it].valid = true;
		
		_currPoints -= (h->pointIndices.size() - 1);
		
		//cout<<_currPoints<<endl;
		
		_pQueueCluster.pop();
		
		//cout<<_currPoints<<" "<<_pQueueCluster.size()<<endl;
	}
	
	cout<<"_combineSplats"<<endl;
}
     
uint32_t HierClustering::_writeOutputToBuffer(void* buffer)
{
    splatPoint* data = (splatPoint*)(buffer);
    uint32_t count = 0;
	
	splatPoint sp;
    
    for(int i=0; i<_allCount; i++)
    {
        if(_splatHierClustering[i].valid)
        {
			sp  = _splatHierClustering[i].sp;
			
			sp.x = reverseFloat((char*)&sp.x);
			sp.y = reverseFloat((char*)&sp.y);
			sp.z = reverseFloat((char*)&sp.z);
			sp.rad = reverseFloat((char*)&sp.rad);
			
			sp.nx = reverseFloat((char*)&sp.nx);
			sp.ny = reverseFloat((char*)&sp.ny);
			sp.nz = reverseFloat((char*)&sp.nz);
			
            data[count++] = sp;
			
			//cout<<count<<" "<<i<<endl;
        }
    }
	
	//cout<<count<<endl;
    
    return count;
}

