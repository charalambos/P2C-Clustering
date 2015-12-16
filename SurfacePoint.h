////////////////////////////////////////////////////////////////////////////////////
// Copyright © Charalambos "Charis" Poullis, charalambos@poullis.org    		  //
// This work can only be used under an exclusive license of the author.           //
////////////////////////////////////////////////////////////////////////////////////

#ifndef __SURFACE_POINT_H__
#define __SURFACE_POINT_H__

#include "Vector.h"
#include "DepthMapPoint.h"

class SurfacePoint : public DepthMapPoint	{
	public:
		SurfacePoint(Vector3f const &_point, Vector2i const &_index);
		~SurfacePoint();

	private:
};


#endif
