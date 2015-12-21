////////////////////////////////////////////////////////////////////////////////////
// Copyright © Charalambos "Charis" Poullis, charalambos@poullis.org    	  //
// This work can only be used under an exclusive license of the author.           //
////////////////////////////////////////////////////////////////////////////////////
#ifndef __PATCH_H__
#define __PATCH_H__

#include "Utilities.h"
#include "BoundingBox.h"
#include "GeometryProcessing.h"
#include "GaussianDistribution.h"
#include "GeospatialBoundingBox.h"

#define MIN_STABLE	3

///Forward declaration
class Surface;

class Patch	{
	public:
		///Constructor
		Patch(int _cluster_id);
		///Destructor
		~Patch();

		///Registers the geospatial bounding box the surface belongs to
		void belongsTo(GeospatialBoundingBox *_geo_box);

		///Get the object the patch belongs to
		GeometricObject *getObject();

		///Returns the point indices
		std::vector<Vector2i> getPointIndices();

		///Returns the points
		std::vector<Vector3f> getPoints();

		///Get the patch's cluster id
		int getClusterId();

		///Registers the surface the patch belongs to
		void belongsTo(Surface *_surface);

		///Returns the belonging surface
		Surface *belongsTo();

		///Returns the centroid of the inlier points
		Vector3f getInlierCentroid();

		///Returns the bounding box of the patch
		BoundingBox *getBoundingBox();

		///These function only affects the first distribution function to update the region with a new point and update the distributions
		void update(Vector3f const &new_point, Vector3f const &new_normal, Vector2i new_point_index, float h_var, float n_var);

		///Function to check if a point is likely to be part of this region
		bool isLikelyToBePartOf(float candidate_elevation, Vector3f const &candidate_normal, float h_var, float n_var);

		///Returns if the patch is stable
		bool isStable();

		///Returns the Gaussian Distribution
		GaussianDistribution6f *getPatchDistribution();

	private:
		///The distribution of the patch descriptors
		GaussianDistribution6f *patch_descriptors_distr;
		///The fitting errors
		std::vector<float> fitting_errors;
		///The point indices of the patch
		std::vector<Vector2i> point_indices;
		///The 3D points
		std::vector<Vector3f> points;
		///The 3D normals
		std::vector<Vector3f> normals;

		///The geospatial bounding box the linear surface belongs to
		GeospatialBoundingBox *geo_box;
		///The id of the corresponding cluster
		int cluster_id;
		///The surface this patch belongs to
		Surface *surface;
		///The bounding box of the patch
		BoundingBox *bbox;

		///The minimum and maximmum allowed probability which determines if a sample should be part of the distribution or not
		float min_probability, max_allowed_probability;

		///Defines if the patch is stable
		bool stable;

		///The sum of all normals
		Vector3f sum_of_normals;
};



#endif

