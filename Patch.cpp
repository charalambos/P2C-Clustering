////////////////////////////////////////////////////////////////////////////////////
// Copyright © Charalambos "Charis" Poullis, charalambos@poullis.org    	  //
// This work can only be used under an exclusive license of the author.           //
////////////////////////////////////////////////////////////////////////////////////
#ifndef __PATCH_CPP__
#define __PATCH_CPP__

#include "Patch.h"
#include "Surface.h"

#define	KAPPA	1.0f

Patch::Patch(int _cluster_id)	{
	patch_descriptors_distr = new GaussianDistribution6f();
	geo_box = 0x00;
	cluster_id = _cluster_id;
	surface = 0x00;
	bbox = new BoundingBox();
	stable = false;

	min_probability = 0.0f;
	max_allowed_probability = 0.0f;

	sum_of_normals = Vector3f(0.0f,0.0f,0.0f);

}

Patch::~Patch()	{
	if (patch_descriptors_distr)	delete patch_descriptors_distr;
	geo_box = 0x00;
	if (bbox)	delete bbox;
}

void Patch::belongsTo(GeospatialBoundingBox *_geo_box)	{
	geo_box = _geo_box;
	return;
}

GeometricObject *Patch::getObject()	{
	return geo_box->getObject();
}

BoundingBox *Patch::getBoundingBox()	{
	return bbox;
}

std::vector<Vector2i> Patch::getPointIndices()	{
	return point_indices;
}

int Patch::getClusterId()	{
	return cluster_id;
}

void Patch::belongsTo(Surface *_surface)	{
	surface = _surface;
	return;
}

Surface *Patch::belongsTo()	{
	return surface;
}

void Patch::update(Vector3f const &new_point, Vector3f const &new_normal, Vector2i new_point_index, float h_var, float n_var)	{
	///Add the point and the point index
	points.push_back(new_point);
	if (!stable && points.size() >= MIN_STABLE)	{
		stable = true;
	}

	///Add the point index of the point
	point_indices.push_back(new_point_index);

	///Add the normal
	normals.push_back(new_normal);

	///This factor controls the variance at the beginning stages so that it doesn't become zero
	float factor = 1.0f/max(1.0f,(float) log(points.size()*points.size()));

	///Form the sample for the new point
	Vector6f sample = Vector6f(new_normal(0),new_normal(1),new_normal(2),new_point(2), h_var, n_var);

	///Fit a GD to the samples again
	patch_descriptors_distr->update(sample);

	///Get the means
	Vector6f means = patch_descriptors_distr->getMeans();
	///Normalize the first 3 components of the means
	Vector3f norm_normal = Vector3f(means(0),means(1),means(2));
	norm_normal.Normalize();

	means = Vector6f(norm_normal(0),norm_normal(1),norm_normal(2),means(3), means(4), means(5));
	patch_descriptors_distr->setMeans(means);

	///The minimum probability will be that of the means
	min_probability = patch_descriptors_distr->prob(means);
	///The maximum allowed probability will be that of the mean normal + or -  a factor = K*variance
	Matrix6f cov_matrix = patch_descriptors_distr->getCovarianceMatrix();
	Vector6f variances = Vector6f(cov_matrix(0,0),cov_matrix(1,1),cov_matrix(2,2),cov_matrix(3,3),cov_matrix(4,4),cov_matrix(5,5));
	variances(0) = max(factor,variances(0));
	variances(1) = max(factor,variances(1));
	variances(2) = max(factor,variances(2));
	variances(3) = max(factor,variances(3));
	variances(4) = max(factor,variances(4));
	variances(5) = max(factor,variances(5));
	max_allowed_probability = patch_descriptors_distr->prob(means + KAPPA*variances);

	Vector3f min_pt = bbox->getMinPt();
	Vector3f max_pt = bbox->getMaxPt();
	///Check if max or min
	if (max_pt(0) < new_point(0))	max_pt(0) = new_point(0);
	if (max_pt(1) < new_point(1))	max_pt(1) = new_point(1);
	if (max_pt(2) < new_point(2))	max_pt(2) = new_point(2);

	if (min_pt(0) > new_point(0))	min_pt(0) = new_point(0);
	if (min_pt(1) > new_point(1))	min_pt(1) = new_point(1);
	if (min_pt(2) > new_point(2))	min_pt(2) = new_point(2);

	bbox->setMinPt(min_pt);
	bbox->setMaxPt(max_pt);

	return;
}

bool Patch::isLikelyToBePartOf(float candidate_elevation, Vector3f const &candidate_normal, float h_var, float n_var)	{
	///If this is the first point of the region then add it
	if ( points.size() == 0 )	return true;

	///Otherwise we need to compute the probability of the new point's sample(normal and rotated depth) of being part of this patch

	///Form the sample
	Vector6f sample = Vector6f(candidate_normal(0), candidate_normal(1), candidate_normal(2), candidate_elevation, h_var, n_var);

	///Compute the probability of the candidate point's information
	float prob =  patch_descriptors_distr->prob(sample);

	///If the probability is within the acceptable range then add the point to the patch i.e. return true
	if (prob < max_allowed_probability)	{
		return true;
	}

	///Otherwise its not part of the patch
	return false;
}

bool Patch::isStable()	{
	return stable;
}

std::vector<Vector3f> Patch::getPoints()	{
	return points;
}

GaussianDistribution6f * Patch::getPatchDistribution()	{
	return patch_descriptors_distr;
}

#endif

