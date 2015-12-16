////////////////////////////////////////////////////////////////////////////////////
// Copyright © Charalambos "Charis" Poullis, charalambos@poullis.org    		  //
// This work can only be used under an exclusive license of the author.           //
////////////////////////////////////////////////////////////////////////////////////


#ifndef __SURFACE_POINT_CPP__
#define __SURFACE_POINT_CPP__

#include "SurfacePoint.h"

SurfacePoint::SurfacePoint(Vector3f const &_point, Vector2i const &_index)	{
	point = _point;
	index = _index;
}

SurfacePoint::~SurfacePoint()	{

}

#endif
