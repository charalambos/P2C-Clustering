////////////////////////////////////////////////////////////////////////////////////
// Copyright © Charalambos "Charis" Poullis, charalambos@poullis.org    	  //
// This work can only be used under an exclusive license of the author.           //
////////////////////////////////////////////////////////////////////////////////////
#ifndef __SURFACE_H__
#define __SURFACE_H__

#include "Image.h"
#include "SurfacePoint.h"
#include "GeometricObject.h"


///Forward declaration
class Patch;
class InformationManager;

class Surface	{
	public:
		///Constructor
		Surface(InformationManager *_information_manager);
		///Destructor
		~Surface();

		///Add a patch to the surface
		void add(Patch *patch);

		///Remove a patch from the surface
		void remove(Patch *patch);

		///The distribution of the 4D plane coefficients
		GaussianDistribution6f *getSurfaceDescriptorsDistribution() const;

		///Returns the list of patches which form the surface
		std::vector<int> getPatches();

		///Determines whether two surfaces should be merged
		friend bool similar(Surface *s1, Surface *s2);

		///Returns the similarity measured as the Bhattacharya distance
		friend float similarity(Surface *s1, Surface *s2);

		///Returns the bounding box of the surface
		BoundingBox *getBoundingBox();

		///Finalize the surface information
		void setInfo(std::vector<Vector2i> const &_point_indices, std::vector<Vector3f> const &_points);

		///Returns the point indices of the surface points to the geo_box's xyz map
		std::vector<Vector2i> getPointIndices();

		///Returns the points belonging to this surface
		std::vector<Vector3f> getPoints();
		std::vector<SurfacePoint> getSurfacePoints();

		///Returns the object resulting from the triangulation of the surface
		GeometricObject *triangulate();

	private:
        ///The bounding box for the surface
		BoundingBox *bbox;
		///The set of patches which constitute the surface
		std::vector<int> patches;
		///The points contained in the patch
		std::vector<SurfacePoint> points;
		///The distribution of the surface descriptors
		GaussianDistribution6f *surface_descriptors_distr;
		///The information manager
		InformationManager *information_manager;

		///static counter
		static int global_surface_counter;
};


#endif
