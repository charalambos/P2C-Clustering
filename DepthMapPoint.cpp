////////////////////////////////////////////////////////////////////////////////////
// Copyright © Charalambos "Charis" Poullis, charalambos@poullis.org    		  //
// This work can only be used under an exclusive license of the author.           //
////////////////////////////////////////////////////////////////////////////////////


#ifndef __DEPTH_MAP_POINT_CPP__
#define __DEPTH_MAP_POINT_CPP__

#include "SurfacePoint.h"

DepthMapPoint::DepthMapPoint()	{

}

DepthMapPoint::DepthMapPoint(Vector3f const &_point, Vector2i const &_index)	{
	point = _point;
	index = _index;
}

DepthMapPoint::~DepthMapPoint()	{

}

Vector3f DepthMapPoint::getPoint() const	{
	return point;
}

Vector2i DepthMapPoint::getIndex() const	{
	return index;
}

void DepthMapPoint::setPoint(Vector3f const &_point)	{
	point = _point;
	return;
}

void DepthMapPoint::setIndex(Vector2i const &_index)	{
	index = _index;
	return;
}

void DepthMapPoint::set(Vector3f const &_point, Vector2i const &_index)	{
	point = _point;
	index = _index;
	return;
}

#endif
