////////////////////////////////////////////////////////////////////////////////////
// Copyright © Charalambos "Charis" Poullis, charalambos@poullis.org    	  //
// This work can only be used under an exclusive license of the author.           //
////////////////////////////////////////////////////////////////////////////////////
#ifndef __SURFACE_CPP__
#define __SURFACE_CPP__

#include "Patch.h"
#include "Surface.h"
#include "InformationManager.h"

int Surface::global_surface_counter = 0;

Surface::Surface(InformationManager *_information_manager)	{
	patches.clear();
	surface_descriptors_distr = new GaussianDistribution6f();
	bbox = new BoundingBox();
	information_manager = _information_manager;
}

Surface::~Surface()	{
	patches.clear();
	if (surface_descriptors_distr)	delete surface_descriptors_distr;
	if (bbox)	delete bbox;
	information_manager = 0x00;
}

void Surface::add(Patch *patch)	{

	///Add the patch to the list of patches
	patches.push_back(patch->getClusterId());

	///Register the surface to the patch
	patch->belongsTo(this);

	///Add the patch's descriptors to the surface descriptors
	std::vector<Vector6f> patch_samples= patch->getPatchDistribution()->getSamples();
	if (patch_samples.size() == 0)	{
		std::cout << "Patch descriptors are zero!!!" << std::endl;
		exit(-1);
	}

	///Check for boundary points i.e. the bounding box
	Vector3f min_pt = bbox->getMinPt();
	Vector3f max_pt = bbox->getMaxPt();
	///Get the patch's bbox
	if (min_pt(0) > patch->getBoundingBox()->getMinPt()(0))	min_pt(0) = patch->getBoundingBox()->getMinPt()(0);
	if (min_pt(1) > patch->getBoundingBox()->getMinPt()(1))	min_pt(1) = patch->getBoundingBox()->getMinPt()(1);
	if (min_pt(2) > patch->getBoundingBox()->getMinPt()(2))	min_pt(2) = patch->getBoundingBox()->getMinPt()(2);

	if (max_pt(0) < patch->getBoundingBox()->getMaxPt()(0))	max_pt(0) = patch->getBoundingBox()->getMaxPt()(0);
	if (max_pt(1) < patch->getBoundingBox()->getMaxPt()(1))	max_pt(1) = patch->getBoundingBox()->getMaxPt()(1);
	if (max_pt(2) < patch->getBoundingBox()->getMaxPt()(2))	max_pt(2) = patch->getBoundingBox()->getMaxPt()(2);
	///Set the points
	bbox->setMinPt(min_pt);
	bbox->setMaxPt(max_pt);

	///Update the Gaussian distribution of the surface descriptors
	surface_descriptors_distr->update(patch_samples);

	///Ensure that the first 3 components of the descriptor remain normalized
	Vector6f means = surface_descriptors_distr->getMeans();
	///Get just the normal
	Vector3f normal = Vector3f(means(0), means(1), means(2));
	normal.Normalize();
	///Reset the means
	surface_descriptors_distr->setMeans(Vector6f(normal(0),normal(1),normal(2), means(3),means(4),means(5)));

	///Sort the ids
	std::sort(patches.begin(),patches.end(),compare_func<int>);


	return;
}

void Surface::remove(Patch *patch)	{
	///Remove the patch
	std::vector<int>::iterator p_itor = patches.begin();
	for (int i=0;i<patches.size();i++,p_itor++)	{
		if (patches[i] == patch->getClusterId())	{
			patches.erase(p_itor);
			break;
		}
	}
	///Request an update on the neighbourhood information from the information manager
	information_manager->updateNeighbours(this);
	return;
}

bool similar(Surface *s1, Surface *s2)	{

	///Compare them to each other by computing the bhattacharya distance
	float similarity = fabs(bhattacharya_distance(s1->getSurfaceDescriptorsDistribution(), s2->getSurfaceDescriptorsDistribution()));

	if (similarity <0.125f)	{
		return true;
	}

	return false;
}

float similarity(Surface *s1, Surface *s2)	{

	///Compare them to each other by computing the bhattacharya distance
	return(fabs(bhattacharya_distance(s1->getSurfaceDescriptorsDistribution(), s2->getSurfaceDescriptorsDistribution())));
}


GaussianDistribution6f *Surface::getSurfaceDescriptorsDistribution() const	{
	return surface_descriptors_distr;
}

std::vector<int> Surface::getPatches()	{
	return patches;
}


BoundingBox *Surface::getBoundingBox()	{
	return bbox;
}

void Surface::setInfo(std::vector<Vector2i> const &_point_indices, std::vector<Vector3f> const &_points)	{
	points.clear();
	for (int i=0;i<_point_indices.size();i++)	{
		points.push_back(SurfacePoint(_points[i], _point_indices[i]));
	}
	return;
}

std::vector<Vector2i> Surface::getPointIndices()	{
	std::vector<Vector2i> point_indices;
	for (int i=0;i<points.size();i++)	point_indices.push_back(points[i].getIndex());
	return point_indices;
}

///Returns the points belonging to this surface
std::vector<Vector3f> Surface::getPoints()	{
	std::vector<Vector3f> surface_points;
	for (int i=0;i<points.size();i++)	surface_points.push_back(points[i].getPoint());
	return surface_points;
}

std::vector<SurfacePoint> Surface::getSurfacePoints()	{
	return points;
}

GeometricObject *Surface::triangulate()	{
	///Find the minimum index
	Vector2i min_pt_index = Vector2i(INT_MAX, INT_MAX);
	Vector2i max_pt_index = Vector2i(-INT_MAX, -INT_MAX);
	for (int i=0;i<points.size();i++)	{
		if (points[i].getIndex()(0) < min_pt_index(0))	min_pt_index(0) = points[i].getIndex()(0);
		if (points[i].getIndex()(1) < min_pt_index(1))	min_pt_index(1) = points[i].getIndex()(1);

		if (points[i].getIndex()(0) > max_pt_index(0))	max_pt_index(0) = points[i].getIndex()(0);
		if (points[i].getIndex()(1) > max_pt_index(1))	max_pt_index(1) = points[i].getIndex()(1);
	}

	///Find the image size occupied by this surface
	int sizex = 1 + max_pt_index(0) - min_pt_index(0);
	int sizey = 1 + max_pt_index(1) - min_pt_index(1);

	///Allocate an array of this size
	Image *surface_map = new Image(sizex, sizey, 0.0f,0.0f,0.0f,1.0f);
	///Put the 3D values corresponding to the point indices
	for (int i=0;i<points.size();i++)	{
		surface_map->setPixel(points[i].getIndex()(0) - min_pt_index(0), points[i].getIndex()(1) - min_pt_index(1), vector2color3<float>(points[i].getPoint()));
	}

	///Triangulate:create a connected mesh of the same size as the surface_map
	std::vector<Face *> new_faces;
	std::vector<Vector3f> new_vertices;
	int vertex_count = 0 ;
	Image *already_added = new Image(sizex,sizey, -1.0f,-1.0f,-1.0f,1.0f);
	for (int y=0;y<sizey-1;y++)	{
		for (int x=0;x<sizex-1;x++)	{
			bool stop_it = false;
			unsigned int offset_x=1;
			unsigned int offset_y=1;
			if (surface_map->getPixel(x,y) == Color(0.0f,0.0f,0.0f)) continue;
			while (surface_map->getPixel(x+offset_x,y) == Color(0.0f,0.0f,0.0f))	{
				offset_x++;
				if (x+offset_x>=sizex-1) {
					stop_it = true;
					break;
				}
			}
			if (stop_it) continue;
			while (surface_map->getPixel(x,y+offset_y) == Color(0.0f,0.0f,0.0f))	{
				offset_y++;
				if (y+offset_y>=sizey-1) {
					stop_it = true;
					break;
				}
			}
			if (stop_it) continue;
			if (y+offset_y >= sizey-1 || x+offset_x >= sizex-1) continue;
			if (surface_map->getPixel(x+offset_x, y+offset_y) == Color(0.0f,0.0f,0.0f)) continue;
			std::vector<int> vertex_indices;
			int index1,index2,index3,index4;
			//add the points
			if (already_added->getPixel(x,y) == Color(-1.0f,-1.0f,-1.0f))	{
				Vector3f point1 = color2vector3<float>(surface_map->getPixel(x,y));
				new_vertices.push_back(point1);
				index1 = vertex_count;
				already_added->setPixel(x,y, Color(float(vertex_count),float(vertex_count),float(vertex_count)));
				vertex_count++;
			}
			else	{
				index1 = int(already_added->getPixel(x,y).r());
			}
			if (already_added->getPixel(x+offset_x,y) == Color(-1.0f,-1.0f,-1.0f))	{
				Vector3f point2 = color2vector3<float>(surface_map->getPixel(x+offset_x,y));
				new_vertices.push_back(point2);
				index2 = vertex_count;
				already_added->setPixel(x+offset_x,y,Color(float(vertex_count),float(vertex_count),float(vertex_count)));
				vertex_count++;
			}
			else	{
				index2 = int(already_added->getPixel(x+offset_x,y).r());
			}
			if (already_added->getPixel(x,y+offset_y) == Color(-1.0f,-1.0f,-1.0f))	{
				Vector3f point3 = color2vector3<float>(surface_map->getPixel(x,y+offset_y));
				new_vertices.push_back(point3);
				index3 = vertex_count;
				already_added->setPixel(x,y+offset_y,Color(float(vertex_count),float(vertex_count),float(vertex_count)));
				vertex_count++;
			}
			else	{
				index3 = int(already_added->getPixel(x,y+offset_y).r());
			}
			if (already_added->getPixel(x+offset_x,y+offset_y) == Color(-1.0f,-1.0f,-1.0f))	{
				Vector3f point4 = color2vector3<float>(surface_map->getPixel(x+offset_x,y+offset_y));
				new_vertices.push_back(point4);
				index4 = vertex_count;
				already_added->setPixel(x+offset_x, y+offset_y, Color(float(vertex_count),float(vertex_count),float(vertex_count)));
				vertex_count++;
			}
			else	{
				index4 = int(already_added->getPixel(x+offset_x,y+offset_y).r());
			}
			//Create 2 faces for each of the triangles form
			Face *new_face=0x00;
			//printf("%d %d %d %d\n",index1,index2,index3,index4);
			vertex_indices.push_back(index1);
			vertex_indices.push_back(index3);
			vertex_indices.push_back(index2);
			//create a new face
			new_face = new Face();
			new_face->setVertices(vertex_indices);
			new_faces.push_back(new_face);

			vertex_indices.clear();
			vertex_indices.push_back(index2);
			vertex_indices.push_back(index3);
			vertex_indices.push_back(index4);
			//create a new face
			new_face = new Face();
			new_face->setVertices(vertex_indices);
			new_faces.push_back(new_face);
		}
	}

	//add a new object
	std::vector<Vector3f> normals;
	std::vector<Vector2f> tex_coords;
	std::vector<Edge *> edges;

	GeometricObject *new_object = new GeometricObject(new_vertices, normals,tex_coords,new_faces,edges);

	delete already_added;
	delete surface_map;
	return new_object;
}


#endif
