////////////////////////////////////////////////////////////////////////////////////
// Copyright © Charalambos "Charis" Poullis, charalambos@poullis.org    		  //
// This work can only be used under an exclusive license of the author.           //
////////////////////////////////////////////////////////////////////////////////////

#ifndef __DEPTH_MAP_POINT_H__
#define __DEPTH_MAP_POINT_H__

#include "Vector.h"

class DepthMapPoint	{
	public:
		DepthMapPoint();
		DepthMapPoint(Vector3f const &_point, Vector2i const &_index);
		~DepthMapPoint();

		Vector3f getPoint() const;

		Vector2i getIndex() const;

		void setPoint(Vector3f const &_point);

		void setIndex(Vector2i const &_index);

		void set(Vector3f const &_point, Vector2i const &_index);

		///Both components must be the same to return i.e. the point as well as the index
		friend bool operator== ( DepthMapPoint const &p1,  DepthMapPoint const &p2)	{
			if (p1.point != p2.point)	return false;
			if (p1.index != p2.index) return false;
			return true;
		}

	protected:
		Vector3f point;
		Vector2i index;
};


#endif
