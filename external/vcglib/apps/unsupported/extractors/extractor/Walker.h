#ifndef __VCGTEST_WALKER
#define __VCGTEST_WALKER

#include "Definitions.h"
#include "Volume.h"

// La classe Walker implementa la politica di visita del volume; conoscendo l'ordine di visita del volume
// è conveniente che il Walker stesso si faccia carico del caching dei dati utilizzati durante l'esecuzione 
// degli algoritmi MarchingCubes ed ExtendedMarchingCubes, in particolare il calcolo del volume ai vertici
// delle celle e delle intersezioni della superficie con le celle. In questo esempio il volume da processare
// viene suddiviso in fette; in questo modo se il volume ha dimensione h*l*w (rispettivamente altezza,
// larghezza e profondità), lo spazio richiesto per il caching dei vertici già allocati passa da O(h*l*w)
// a O(h*l). 
class Walker
{
private:
	typedef int VertexIndex;

public:
	Walker(const BoundingBox &bbox, const vcg::Point3i &resolution)
	{
		_bbox				= bbox;
		_resolution = resolution;
		_cell_size.X() = _bbox.DimX()/_resolution.X();
		_cell_size.Y() = _bbox.DimY()/_resolution.Y();
		_cell_size.Z() = _bbox.DimZ()/_resolution.Z();
		_slice_dimension = resolution.X()*resolution.Z();

		_x_cs = new VertexIndex[ _slice_dimension ];
		_y_cs = new VertexIndex[ _slice_dimension ];
		_z_cs = new VertexIndex[ _slice_dimension ];
		_x_ns = new VertexIndex[ _slice_dimension ];
		_z_ns = new VertexIndex[ _slice_dimension ];
		_v_cs = new float[_slice_dimension];
		_v_ns = new float[_slice_dimension];

	};

	~Walker()
	{}

	template<class EXTRACTOR_TYPE>
	void BuildMesh(Mesh &mesh, Volume &volume, EXTRACTOR_TYPE &extractor)
	{
		_volume = &volume;
		_mesh		= &mesh;
		_mesh->Clear();
		vcg::Point3i p1, p2;

		Begin();
		extractor.Initialize();
		for (int j=_bbox.min.Y(); j<_bbox.max.Y()-_cell_size.Y(); j+=_cell_size.Y())
		{
			for (int i=_bbox.min.X(); i<_bbox.max.X()-_cell_size.X(); i+=_cell_size.X())
			{
				for (int k=_bbox.min.Z(); k<_bbox.max.Z()-_cell_size.Z(); k+=_cell_size.Z())
				{
					p1.X()=i;									p1.Y()=j;									p1.Z()=k;
					p2.X()=i+_cell_size.X();	p2.Y()=j+_cell_size.Y();	p2.Z()=k+_cell_size.Z();
					extractor.ProcessCell(p1, p2);
				}
			}
			NextSlice();
		}
		extractor.Finalize();
		_volume = NULL;
		_mesh		= NULL;
	};

	float V(int pi, int pj, int pk)
	{
		int i = pi - _bbox.min.X();
		int k = pk - _bbox.min.Z();
		return (pj==_current_slice) ? _v_cs[i+k*_resolution.X()] : _v_ns[i+k*_resolution.X()];
	}

	bool Exist(const vcg::Point3i &p1, const vcg::Point3i &p2, VertexPointer &v)
	{ 
		int i_idx = p1.X()-_bbox.min.X();
		int k_idx = p2.Z()-_bbox.min.Z();
		int index = i_idx+k_idx*_resolution.X();
		if (p1.X()!=p2.X())					//intersezione della superficie con un Xedge
			return (p1.Y()==_current_slice)? _x_cs[index]!=-1 : _x_ns[index]!=-1;
		else if (p1.Y()!=p2.Y())		//intersezione della superficie con un Yedge
			return _y_cs[index]!=-1;
		else if (p1.Z()!=p2.Z())		//intersezione della superficie con un Zedge
			return (p1.Y()==_current_slice)? _z_cs[index]!=-1 : _z_ns[index]!=-1;

		assert(false); // impossibile: i due punti non erano allineati rispetto a nessuna direzione
		return false;
	}

	void GetXIntercept(const vcg::Point3i &p1, const vcg::Point3i &p2, VertexPointer &v) 
	{ 
		int i = p1.X() - _bbox.min.X();
		int z = p1.Z() - _bbox.min.Z();
		VertexIndex index = i+z*_resolution.X();
		VertexIndex pos;
		if (p1.Y()==_current_slice)
		{
			if ((pos=_x_cs[index])==-1)
			{
				_x_cs[index] = (VertexIndex) _mesh->vert.size();
				pos = _x_cs[index];
				Allocator::AddVertices( *_mesh, 1 );
				v = &_mesh->vert[pos];
				_volume->GetXIntercept(p1, p2, v);
				return;
			}
		}
		if (p1.Y()==_current_slice+_cell_size.Y())
		{
			if ((pos=_x_ns[index])==-1)
			{
				_x_ns[index] = (VertexIndex) _mesh->vert.size();
				pos = _x_ns[index];
				Allocator::AddVertices( *_mesh, 1 );
				v = &_mesh->vert[pos];
				_volume->GetXIntercept(p1, p2, v);
				return;
			}
		}
		v = &_mesh->vert[pos];
	}
	void GetYIntercept(const vcg::Point3i &p1, const vcg::Point3i &p2, VertexPointer &v) 
	{
		int i = p1.X() - _bbox.min.X();
		int z = p1.Z() - _bbox.min.Z();
		VertexIndex index = i+z*_resolution.X();
		VertexIndex pos;
		if ((pos=_y_cs[index])==-1)
		{
			_y_cs[index] = (VertexIndex) _mesh->vert.size();
			pos = _y_cs[index];
			Allocator::AddVertices( *_mesh, 1);
			v = &_mesh->vert[ pos ];
			_volume->GetYIntercept(p1, p2, v);
		}
		v = &_mesh->vert[pos];
	}
	void GetZIntercept(const vcg::Point3i &p1, const vcg::Point3i &p2, VertexPointer &v) 
	{
		int i = p1.X() - _bbox.min.X();
		int z = p1.Z() - _bbox.min.Z();
		VertexIndex index = i+z*_resolution.X();
		VertexIndex pos;
		if (p1.Y()==_current_slice)
		{
			if ((pos=_z_cs[index])==-1)
			{
				_z_cs[index] = (VertexIndex) _mesh->vert.size();
				pos = _z_cs[index];
				Allocator::AddVertices( *_mesh, 1 );
				v = &_mesh->vert[pos];
				_volume->GetZIntercept(p1, p2, v);
				return;
			}
		}
		if (p1.Y()==_current_slice+_cell_size.Y())
		{
			if ((pos=_z_ns[index])==-1)
			{
				_z_ns[index] = (VertexIndex) _mesh->vert.size();
				pos = _z_ns[index];
				Allocator::AddVertices( *_mesh, 1 );
				v = &_mesh->vert[pos];
				_volume->GetZIntercept(p1, p2, v);
				return;
			}
		}
		v = &_mesh->vert[pos];
	}

protected:
	BoundingBox		_bbox;
	vcg::Point3i	_resolution;
	vcg::Point3i	_cell_size;

	int _slice_dimension;
	int	_current_slice;
  
	float *_v_cs; // il valore del campo campionato nella fetta di volumecorrente 
	float *_v_ns; // il valore del campo campionato nella prossima fetta di volume
	
	VertexIndex *_x_cs; // indici dell'intersezioni della superficie lungo gli Xedge della fetta corrente
	VertexIndex	*_y_cs; // indici dell'intersezioni della superficie lungo gli Yedge della fetta corrente
	VertexIndex *_z_cs; // indici dell'intersezioni della superficie lungo gli Zedge della fetta corrente
	VertexIndex *_x_ns; // indici dell'intersezioni della superficie lungo gli Xedge della prossima fetta 
	VertexIndex *_z_ns; // indici dell'intersezioni della superficie lungo gli Zedge della prossima fetta 

	Mesh		*_mesh;
	Volume	*_volume;

	void NextSlice() 
	{
		memset(_x_cs, -1, _slice_dimension*sizeof(VertexIndex));
		memset(_y_cs,	-1, _slice_dimension*sizeof(VertexIndex));
		memset(_z_cs, -1, _slice_dimension*sizeof(VertexIndex));

		std::swap(_x_cs, _x_ns);
		std::swap(_z_cs, _z_ns);		
		std::swap(_v_cs, _v_ns);
		
		_current_slice += _cell_size.Y();
		int j						= _current_slice + _cell_size.Y();
		int k_idx, i_idx, index;
		for (int i=_bbox.min.X(); i<_bbox.max.X(); i+=_cell_size.X())
		{
			i_idx = i-_bbox.min.X();
			for (int k=_bbox.min.Z(); k<_bbox.max.Z(); k+=_cell_size.Z())
			{
				k_idx = k-_bbox.min.Z();
				index = i_idx+k_idx*_resolution.X();
				_v_ns[ index ] = _volume->V(i, j, k);
			}
		}
	}

	void Begin()
	{
		_current_slice = _bbox.min.Y();

		memset(_x_cs, -1, _slice_dimension*sizeof(VertexIndex));
		memset(_y_cs, -1, _slice_dimension*sizeof(VertexIndex));
		memset(_z_cs, -1, _slice_dimension*sizeof(VertexIndex));
		memset(_x_ns, -1, _slice_dimension*sizeof(VertexIndex));
		memset(_z_ns, -1, _slice_dimension*sizeof(VertexIndex));
		
		int index;
		int j = _current_slice;
		int i_idx, k_idx;
		for (int i=_bbox.min.X(); i<_bbox.max.X(); i+=_cell_size.X())
		{
			i_idx = i-_bbox.min.X();
			for (int k=_bbox.min.Z(); k<_bbox.max.Z(); k+=_cell_size.Z())
			{
				k_idx = k-_bbox.min.Z();
				index = i_idx+k_idx*_resolution.X();
				_v_cs[index] = _volume->V(i, j, k);
				_v_ns[index] = _volume->V(i, j+_cell_size.Y(), k);
			}
		}
	}
};
#endif // __VCGTEST_WALKER