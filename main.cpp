////////////////////////////////////////////////////////////////////////////////////
// Copyright © Charalambos "Charis" Poullis, charalambos@poullis.org    	      //
// This work can only be used under an exclusive license of the author.           //
////////////////////////////////////////////////////////////////////////////////////

#include <sys/stat.h>

#include "Vector.h"
#include "BoundingBox.h"
#include "SurfacePoint.h"
#include "GeometryExporter.h"
#include "InformationManager.h"
#include "GeospatialBoundingBox.h"

#define DEFAULT_RES		        2048
#define EXPORT_PATCH            1
#define EXPORT_SURFACE          1
#define UNSTABLE_SURFACE	    10

///The data file name
std::string file_name;

///Uses probabilistic region grow to cluster points into patches
void probabilisticRegionGrow(InformationManager *information_manager, GeospatialBoundingBox *geo_box, Image *patch_index_map);
///Computes all neighbourhood information
void deriveNeighbourhoodInformation(InformationManager *information_manager, GeospatialBoundingBox *geo_box, Image *patch_index_map);
///Initializes the structures relating to the surfaces
void initializeSurfaces(InformationManager *information_manager);
///Merges similar (bhatacharrya distance) surfaces together
void mergeSimilarSurfaces(InformationManager *information_manager);
///Merges small surfaces to neighbouring surfaces
void mergeSmallSurfaces(InformationManager *information_manager, GeospatialBoundingBox *geo_box);
///Finalizes the surface information after all the editing
void finalizeSurfaceInformation(GeospatialBoundingBox *geo_box, InformationManager *information_manager);
///Builds a surface index map
Image *buildSurfaceIndexMap(int id, InformationManager *information_manager, GeospatialBoundingBox *geo_bbox, std::string const &step_label);
///Export all the surfaces as geometry
int exportSurfaces(int id, InformationManager *information_manager, GeospatialBoundingBox *geo_bbox, std::string const &step_label);
///Export all the patches as geometry
int exportPatches(int id, InformationManager *information_manager, GeospatialBoundingBox *geo_bbox, std::string const &step_label);
///Returns the points in a geo box
void getPoints(GeospatialBoundingBox *geo_box, std::vector<Vector2i> const &point_indices, std::vector<Vector3f> &points);


void probabilisticRegionGrow(InformationManager *information_manager, GeospatialBoundingBox *geo_box, Image *patch_index_map)	{
	///The search area used for growing
	int search_area = 1;
	///The region index
	int patch_index = 0;
	///Get the xyz map of the geo box
	Image *xyz_map = geo_box->getXYZMap();

	///Get the normal map based on the XYZ map of the geo_box
	Image *normal_map = GeometryProcessing::computeNormalMap(xyz_map, true);

	///Get the height variation between the points
	Image *height_and_dot_variation_map = geo_box->getHeightAndNormalVariationMap();


	///A map to hold the already processed points in the map
	Image *already_processed = new Image( xyz_map->getWidth(),xyz_map->getHeight(), 0.0f, 0.0f, 0.0f, 1.0f);

	int patch_image_index=1;
	Image *patch_image = new Image(patch_index_map->getWidth(),patch_index_map->getHeight());

	///Begin growing regions
	for ( int y=0;y<xyz_map->getHeight();y++ )	{
		for ( int x=0;x<xyz_map->getWidth();x++ )	{
			///if this point is valid
			if ( xyz_map->getPixel(x,y) == Color(0.0f,0.0f,0.0f))	continue;
			///if its not processed
			if (already_processed->getPixel(x,y) == Color(1.0f,1.0f,1.0f))	continue;
			///get the xyz value
			Vector3f point = color2vector3<float>(xyz_map->getPixel(x,y));
			///get the normal
			Vector3f normal = color2vector3<float>(normal_map->getPixel(x,y));;
			///get the point index
			Vector2i point_index = Vector2i(x,y);
			///get the height and normal variation
			float h_var = color2vector3<float>(height_and_dot_variation_map->getPixel(x,y))(0);
			float n_var = color2vector3<float>(height_and_dot_variation_map->getPixel(x,y))(1);
			///create a region with this point added
			Patch *a_patch = new Patch(patch_index);
			///add the initial point, normal and index to the region
			a_patch->update(point, normal, point_index, h_var, n_var );
			///add the point to the region map
			patch_index_map->setPixel(x,y,Color(float(patch_index)));

			///create a processing list with this point on it
			std::vector<Vector2i> processing_list, marked_points;
			marked_points.push_back(point_index);
			processing_list.push_back ( point_index );
			///mark it as processed
			already_processed->setPixel(x,y, Color(1.0f,1.0f,1.0f));

			///start the processing
			while ( processing_list.size() )	{
				///get the current point
				Vector2i current_point_index = processing_list[processing_list.size()-1];
				///remove it from the list
				processing_list.pop_back();

				///check the local neighbourhood
				for (int i=current_point_index(1) - search_area;i<=current_point_index(1) + search_area;i++)	{
					for (int j=current_point_index(0) - search_area;j<=current_point_index(0) + search_area;j++)	{
						///if its the same
						if ( i==y && j==x )	continue;
						///if its inbound
						if (outOfBounds(xyz_map, j, i))	continue;
						///if its not processed
						if (already_processed->getPixel(j,i) == Color(1.0f,1.0f,1.0f))	continue;
						///if its a valid point
						if (xyz_map->getPixel(j,i) == Color(0.0f,0.0f,0.0f))	continue;

						///get the saliency and normal of this point
						Vector3f candidate_point = color2vector3<float>(xyz_map->getPixel(j,i));
						///get the height and normal variation
						float candidate_h_var = color2vector3<float>(height_and_dot_variation_map->getPixel(j,i))(0);
						float candidate_n_var = color2vector3<float>(height_and_dot_variation_map->getPixel(j,i))(1);

						Vector3f candidate_normal = color2vector3<float>(normal_map->getPixel(j,i));
						///check is the distance of this position from the plane
						///and the normal difference is similar to the other points
						if (a_patch->isLikelyToBePartOf(candidate_point(2), candidate_normal, candidate_h_var, candidate_n_var))	{
							Vector2i candidate_point_index = Vector2i(j,i);
							///add it to the list
							processing_list.push_back(candidate_point_index);
							///add it to the region
							a_patch->update(candidate_point,candidate_normal,candidate_point_index, candidate_h_var, candidate_n_var);
							///mark it as processed
							already_processed->setPixel(j,i,Color(1.0f,1.0f,1.0f));
							///add it to the region map
							patch_index_map->setPixel(j,i,Color(float(patch_index)));
							///Save out an image of the region with the new point. Used for animations
							if (0)	{//DEBUG) {
								patch_image->setPixel(j,i,Color(patch_image_index));
								patch_image->saveImage(_format("region_patch_%.06d.png",patch_image_index++));
							}

							marked_points.push_back(candidate_point_index);
						}
					}
				}
			}

			///add the region to the resulting regions if it's stable otherwise just skip it
			if (a_patch->isStable())	{
				information_manager->addPatch(a_patch);
				patch_index++;
			}
			else	{
				///Unmark any points that were marked
				for (int i=0;i<marked_points.size();i++)	{
					patch_index_map->setPixel(marked_points[i](0), marked_points[i](1), Color(-1.0f));
				}
			}
		}
	}

	///Clean up
	delete already_processed;
	delete patch_image;

	return;
}

void deriveNeighbourhoodInformation(InformationManager *information_manager, GeospatialBoundingBox *geo_box, Image *patch_index_map)	{

	///Get the xyz_mapcorresponding to the geospatial bounding box
	Image *xyz_map = geo_box->getXYZMap();

	///Allocate the proper memory for the neighbourhood information
	information_manager->initializePatchNeighbourInfo();//clusters.size());

	///Go through the patch index map
	for (int y=0;y<patch_index_map->getHeight();y++)	{
		for (int x=0;x<patch_index_map->getWidth();x++)	{
			///Get the patch index of the current pixel
			int current_patch_index = _round(patch_index_map->getPixel(x,y).r());
			///if it's invalid then continue
			if (current_patch_index == -1)	continue;
			///Check the immediate neighbourhood of this pixel
			for (int i=-1;i<=1;i++)	{
				for (int j=-1;j<=1;j++)	{
					///Check for out of bounds
					if (outOfBounds(patch_index_map, x+j, y+i))	continue;
					///Get the neighbour's patch index
					int neighbour_patch_index = _round(patch_index_map->getPixel(x+j,y+i).r());
					///Make sure they are not in the same patch
					if (current_patch_index == neighbour_patch_index)	continue;
					///If it's invalid then continue
					if (neighbour_patch_index == -1)	continue;
					///Otherwise add the neighbouring relationship
					information_manager->addNeighbourPatchIndex(current_patch_index, neighbour_patch_index);
				}
			}
		}
	}

	///Sort the patch neighbours based on their patch indices. This is required for computing the difference later.
// 	information_manager->sortPatchNeighbours();
/*
	///Debugging
	for (int i=0;i<information_manager->getNumberOfPatches();i++)	{
		std::cout << "Patch " << i << " neighbours: ";
		for (unsigned int j=0;j<information_manager->getPatchNeighboursOf(i).size();j++)	{
			std::cout << information_manager->getPatchNeighboursOf(i)[j];
			if (j + 1 != information_manager->getNumberOfPatches()) std::cout << ", ";
		}
		std::cout << std::endl;
	}
*/
	return;
}

void initializeSurfaces(InformationManager *information_manager)	{

	///Go through the patches and initialize a surface for each patch
	for (int i=0;i<information_manager->getNumberOfPatches();i++)	{
		Patch *patch = information_manager->getPatch(i);
		///Create a new surface object
		Surface *surface = new Surface(information_manager);
		///Add the surface
		information_manager->addSurface(surface);
		///Add the patch to a new surface
		std::vector<int> patches;
		patches.push_back(i);
		information_manager->addPatchesToSurface(surface, patches);
	}

	///Call this to make a stored copy to the neighbouring surface pointers.
	information_manager->reviewSurfaceNeighbours();

	return;
}


void mergeSimilarSurfaces(InformationManager *information_manager)	{


	int variable_size = information_manager->getNumberOfSurfaces();
	bool changes = true;
	int iterations = 1;
	while (changes)	{
// 		{
// 			boost::mutex::scoped_lock lock(thread_mutex);
			std::cout << _format("Surfaces at iteration %d: %d",iterations++, information_manager->getNumberOfSurfaces()) << std::endl;
			std::cout << "Iteration " << iterations << ": " << timestamp() << std::endl;
// 		}

		///Change the flag to false
		changes = false;
		///Go through the surfaces and compare them to each other
		for (int i=0;i<variable_size;i++)	{
			//std::cout << "Variable size: " << variable_size << std::endl;
			///Get the current surface
			Surface *current_surface = information_manager->getSurface(i);
			///Get the neighbouring surfaces of the current surface
			std::vector<Surface *> surface_neighbours = information_manager->getSurfaceNeighboursOf(i);
			///The surfaces that will be merged
			std::vector<Surface *> similar_surfaces;
			///and their associated indices
			std::vector<int> similar_surfaces_indices;
			///Go through the neighbours and check if they are similar enough for merging
			for (int j=0;j<surface_neighbours.size();j++)	{
				///Get the surface this patch belongs to
				Surface *neighbouring_surface = surface_neighbours[j];
				///Compare the two surfaces for similirity
				if (similar(current_surface,neighbouring_surface))	{
					///If they are similar then add it to the list of similar surfaces
					similar_surfaces.push_back(neighbouring_surface);
					///Find the correct index of the surface in the surfaces array stored in the information manager
					int neighbouring_surface_index = information_manager->getSurfaceIndex(neighbouring_surface);
					if (neighbouring_surface_index < 0) std::cout << "problem" << std::endl;
					///Add the index to the indices array
					similar_surfaces_indices.push_back(neighbouring_surface_index);
				}
			}
			///If there are other similar surfaces then merge them and remove them from the list
			if (similar_surfaces.size()!=0)	{
/*std::cout << "Before" << std::endl;
for (int k=0;k<information_manager->getSurfaceNeighbours().size();k++)	{
	std::cout << "Surface " << k << std::endl;
	SURFACE_NEIGHBOUR_INFO neighbours = information_manager->getSurfaceNeighbours()[k];
	std::cout << "\tNeighbours: ";
	for (int l=0;l<neighbours.neighbour_patch_indices.size();l++)	{
		std::cout << neighbours.neighbour_patch_indices[l] << " ";
	}
	std::cout <<std::endl;
	std::cout << "\tPatches: ";
	for (int l=0;l<information_manager->getSurface(k)->getPatches().size();l++)	{
		std::cout << information_manager->getSurface(k)->getPatches()[l] << " ";
	}
	std::cout <<std::endl;
}*/
				///Merge the similar surfaces and the current surface and create a new one
				similar_surfaces.push_back(current_surface);
				similar_surfaces_indices.push_back(i);
				///The merge function already adds the new surface to the list contained in the information manager
				Surface *new_surface = information_manager->merge(similar_surfaces);
				///Review the neighbourhood information for the surfaces i.e. any surfaces which where neighbours to the
				///a surface which was removed (similar_surfaces) will be updated to the newly created surface
				information_manager->reviewSurfaceNeighbours(similar_surfaces, new_surface);

				///Adjust the size
				variable_size += 1 - similar_surfaces.size();
				///Adjust the pointer
				i--;
				///Mark the flags
				changes = true;


// std::cout << "Removing: ";
// for (int k=0;k<similar_surfaces_indices.size();k++)	{
// 	std::cout << similar_surfaces_indices[k] << " ";
// }
// std::cout << std::endl;


				///Erase the old surfaces from the list. Remember that the current surface is also on the list
				remove(information_manager->getSurfaces(), similar_surfaces);
				///NOTE: This MUST be in increasing order
				std::sort(similar_surfaces_indices.begin(),similar_surfaces_indices.end(),compare_func<int>);
				remove(information_manager->getSurfaceNeighbours(), similar_surfaces_indices);

// std::cout << "After" << std::endl;
// for (int k=0;k<information_manager->getSurfaceNeighbours().size();k++)	{
// 	std::cout << "Surface " << k << std::endl;
// 	SURFACE_NEIGHBOUR_INFO neighbours = information_manager->getSurfaceNeighbours()[k];
// 	std::cout << "\tNeighbours: ";
// 	for (int l=0;l<neighbours.neighbour_patch_indices.size();l++)	{
// 		std::cout << neighbours.neighbour_patch_indices[l] << " ";
// 	}
// 	std::cout <<std::endl;
// 	std::cout << "\tPatches: ";
// 	for (int l=0;l<information_manager->getSurface(k)->getPatches().size();l++)	{
// 		std::cout << information_manager->getSurface(k)->getPatches()[l] << " ";
// 	}
// 	std::cout <<std::endl;
// }

			}
// std::cout << "surfaces af: " << information_manager->getSurfaces().size() << std::endl;
			//if (changes)	break;
		}
	}

	return;
}

void mergeSmallSurfaces(InformationManager *information_manager, GeospatialBoundingBox *geo_box)	{


	int variable_size = information_manager->getNumberOfSurfaces();
	bool changes = true;
	int iterations = 1;
	while (changes)	{
// 		{
// 			boost::mutex::scoped_lock lock(thread_mutex);
			std::cout << _format("Surfaces at iteration %d: %d",iterations++, information_manager->getNumberOfSurfaces()) << std::endl;
// 		}

		///Change the flag to false
		changes = false;
		///Go through the surfaces and compare them to each other
		for (int i=0;i<variable_size;i++)	{
			///std::cout << "Variable size: " << variable_size << " i: "<< i << std::endl;
			///Get the current surface
			Surface *current_surface = information_manager->getSurface(i);
			///Get the number of points contained in this surface
			std::vector<Vector2i> point_indices;
			information_manager->getPointIndicesForSurface(i, point_indices);
			int number_of_surface_points = point_indices.size();
			///If the number of points is big enough then move on
			if (number_of_surface_points > UNSTABLE_SURFACE/*median_number_of_surface_points*/)	continue;
			///Otherwise find the closest neighbour and add it to it
			///Get the neighbouring surfaces of the current surface
			std::vector<Surface *> surface_neighbours = information_manager->getSurfaceNeighboursOf(i);
			///The surfaces that will be merged
			std::vector<Surface *> similar_surfaces;
			///and their associated indices
			std::vector<int> similar_surfaces_indices;
			///The index of the surface with the lowest bhattacharya distance
			Surface *best_neighbour = 0x00;
			float best_similarity = FLT_MAX;
			///Go through the neighbours and check which neighbour has the smallest bhatacharrya distance
			for (int j=0;j<surface_neighbours.size();j++)	{
				///Get the surface this patch belongs to
				Surface *neighbouring_surface = surface_neighbours[j];
				///Get the similarity with this neighbouring surface
				float sim = similarity(current_surface, neighbouring_surface);
				///Check if it's the smallest so far
				if (sim < best_similarity)	{
					best_similarity = sim;
					best_neighbour = neighbouring_surface;
				}
			}

			///If there is at least one neighbour which can be merged with the current surface
			if (best_neighbour != 0x00)	{
				///If they are similar then add it to the list of similar surfaces
				similar_surfaces.push_back(best_neighbour);
				///Find the correct index of the surface in the surfaces array stored in the information manager
				int neighbouring_surface_index = information_manager->getSurfaceIndex(best_neighbour);
				if (neighbouring_surface_index < 0) std::cout << "problem" << std::endl;
				///Add the index to the indices array
				similar_surfaces_indices.push_back(neighbouring_surface_index);
				//std::cout << "Merging surface id " << i << " with " << number_of_surface_points << " with  surface id " << neighbouring_surface_index << " with " << best_neighbour->getPointIndices().size() << std::endl;

				///Merge the similar surfaces and the current surface and create a new one
				similar_surfaces.push_back(current_surface);
				similar_surfaces_indices.push_back(i);
				///The merge function already adds the new surface to the list contained in the information manager
				Surface *new_surface = information_manager->merge(similar_surfaces);
				///Review the neighbourhood information for the surfaces i.e. any surfaces which where neighbours to the
				///a surface which was removed (similar_surfaces) will be updated to the newly created surface
				information_manager->reviewSurfaceNeighbours(similar_surfaces, new_surface);

				///Adjust the size
				variable_size += 1 - similar_surfaces.size();
				///Adjust the pointer
				i--;
				///Mark the flags
				changes = true;

				///Erase the old surfaces from the list. Remember that the current surface is also on the list
				remove(information_manager->getSurfaces(), similar_surfaces);
				///NOTE: This MUST be in increasing order
				std::sort(similar_surfaces_indices.begin(),similar_surfaces_indices.end(),compare_func<int>);
				remove(information_manager->getSurfaceNeighbours(), similar_surfaces_indices);

				///Update the information of the added surface
				///Get all the point indices of the surface
				std::vector<Vector2i> point_indices;
				information_manager->getPointIndicesForSurface(variable_size-1, point_indices);
				///Get all the 3D points of the surface
				std::vector<Vector3f> points;
				getPoints(geo_box, point_indices, points);
				///Set the point indices and points to the surface
				information_manager->getSurface(variable_size-1)->setInfo(point_indices, points);
			}
		}
	}

	return;
}

void getPoints(GeospatialBoundingBox *geo_box, std::vector<Vector2i> const &point_indices, std::vector<Vector3f> &points)	{
	for (int i=0;i<point_indices.size();i++)	{
		points.push_back(color2vector3<float>(geo_box->getXYZMap()->getPixel(point_indices[i](0), point_indices[i](1))));
	}

	return;
}

void finalizeSurfaceInformation(GeospatialBoundingBox *geo_box, InformationManager *information_manager)	{
	///Go through all the surfaces and set their final points and indices
	for (int i=0;i<information_manager->getSurfaces().size();i++)	{
		///Get all the point indices of the surface
		std::vector<Vector2i> point_indices;
		information_manager->getPointIndicesForSurface(i, point_indices);
		///Get all the 3D points of the surface
		std::vector<Vector3f> points;
		getPoints(geo_box, point_indices, points);
		///Set the point indices and points to the surface
		information_manager->getSurface(i)->setInfo(point_indices, points);
	}

	return;
}

Image *buildSurfaceIndexMap(int id, InformationManager *information_manager, GeospatialBoundingBox *geo_bbox, std::string const &step_label)	{

	///Create a surface index map
	Image *surface_index_map = new Image(geo_bbox->getXYZMap()->getWidth(), geo_bbox->getXYZMap()->getHeight(), 0.0f,0.0f,0.0f,1.0f);

	///Go through the surfaces
	for (int i=0;i<information_manager->getSurfaces().size();i++)	{
		///Get the surface points for each surface
		std::vector<SurfacePoint> surface_points = information_manager->getSurface(i)->getSurfacePoints();
		for (int j=0;j<surface_points.size();j++)	{
			///Mark the index in the index map
			surface_index_map->setPixel(surface_points[j].getIndex()(0), surface_points[j].getIndex()(1), Color(float(i)));
		}
	}

	surface_index_map->saveImage(_format("%d/surface_index_map_%s.pfm", id, step_label.c_str()).c_str());

	return surface_index_map;
}

int exportPatches(int id, InformationManager *information_manager, GeospatialBoundingBox *geo_bbox, std::string const &step_label)	{
	///Create a color coded map for the patches
	Image *color_coded_patch_map = new Image(geo_bbox->getXYZMap()->getWidth(), geo_bbox->getXYZMap()->getHeight(), 0.0f,0.0f,0.0f,1.0f);
 	///Export each patch separately
	std::cout << "Exporting patch map..." << std::endl;


	int total_patch_objects = 0;
	for (int i=0;i<information_manager->getPatches().size();i++)	{
		///Set the right color
	float r= ((double) rand() / (RAND_MAX));
	float g= ((double) rand() / (RAND_MAX));
	float b= ((double) rand() / (RAND_MAX));

		///Get the point indices for each patch
		if (!information_manager->getPatch(i)) {
			std::cout << "NULL patch " << i << std::endl;
			continue;
		}
		std::vector<Vector2i> point_indices = information_manager->getPatch(i)->getPointIndices();
		///Get the vertices
		std::vector<Vector3f> patch_vertices;
		for (int k=0;k<point_indices.size();k++)	{
			///Mark the point indices in the color coded map
			color_coded_patch_map->setPixel(point_indices[k](0), point_indices[k](1), Color(r,g,b));
		}

		total_patch_objects++;
	}

	///Save out the color coded patch map
	color_coded_patch_map->saveImage(_format("%d/color_coded_patch_map_%s.png", id, step_label.c_str()).c_str());

	std::cout << "...done." << std::endl;

	///Clean up
	delete color_coded_patch_map;

	return total_patch_objects;
}

int exportSurfaces(int id, InformationManager *information_manager, GeospatialBoundingBox *geo_bbox, std::string const &step_label)	{
	///Create a color coded map for the patches
	Image *color_coded_surface_map = new Image(geo_bbox->getXYZMap()->getWidth(), geo_bbox->getXYZMap()->getHeight(), 0.0f,0.0f,0.0f,1.0f);
	///Export each surface separately
	std::cout << "Exporting surface map..." << std::endl;

	int total_surface_objects = 0;
	for (int i=0;i<information_manager->getSurfaces().size();i++)	{
		///Set the right color
        float r= ((double) rand() / (RAND_MAX));
        float g= ((double) rand() / (RAND_MAX));
        float b= ((double) rand() / (RAND_MAX));


		///Get the point indices for each patch
		std::vector<Vector2i> point_indices = information_manager->getSurface(i)->getPointIndices();
		for (int k=0;k<point_indices.size();k++)	{
			///Mark the point indices in the color coded map
			color_coded_surface_map->setPixel(point_indices[k](0), point_indices[k](1), Color(r,g,b));
		}

		total_surface_objects++;
	}

	///Save out the color coded patch map
	color_coded_surface_map->saveImage(_format("%d/color_coded_surface_map_%s.png", id, step_label.c_str()).c_str());
	///Clean up
	delete color_coded_surface_map;

	return total_surface_objects;
}

int main(int argc, char *argv[])    {
    int i =0; ///Can be used to convert the code to process a series of geoboxes.

    if (argc != 2)  {
        std::cout << "Number of arguments should be exactly 2. Quiting..." << std::endl;
        std::cout << "Usage: P2C-Clustering structured_geometry.pfm" << std::endl;
        return 0;
    }
    else    {
        file_name = std::string(argv[1]);
        std::cout << "Processing file: " << file_name << std::endl;
    }

    Image *xyz_map = new Image();
    if (!xyz_map->loadImage(file_name)) {
        std::cerr << "Problem loading input file: " << file_name << std::endl;
        return 0;
    }
    xyz_map->saveImage("input_structured_geometry.pfm");

    ///Create a geometry exporter object
    GeometryExporter *ge = new GeometryExporter();

    ///Initialize a geo bounding box with the data
    GeospatialBoundingBox *geo_bbox = new GeospatialBoundingBox();

    ///Initialize the data required for further processing
    geo_bbox->initialize(xyz_map);

    ///Export the information files for the geobox
    geo_bbox->exportInformationFile();

    if (geo_bbox->getNumberOfPoints() == 0)	    {
        std::cerr << "There are 0 points in the geometry file." << std::endl;
        return 0;
    }

    ///initialize random seed:
    srand (time(NULL));

    ///Create a folder for the results of this geo box
    mkdir(_format("%d",i).c_str(),0777);

    ///Create an information manager object
    InformationManager *information_manager = new InformationManager();

    //////////////////////////////////////////////////////////////////////////////////////////
    ///PHASE ONE: POINT CLUSTERING															//
    ///Using the probabilistic region growing algorithm points are clustered into patches	//
    //////////////////////////////////////////////////////////////////////////////////////////
    ///Perform probabilistic region growing to get the points grouped into initial patches
    Image *patch_index_map = new Image(geo_bbox->getXYZMap()->getWidth(),geo_bbox->getXYZMap()->getHeight(), -1.0f,-1.0f,-1.0f,1.0f);
    std::cout << "Starting point clustering at " << timestamp() << std::endl;
    probabilisticRegionGrow(information_manager, geo_bbox, patch_index_map);
    std::cout << "Finished point clustering at " << timestamp() << std::endl;

    if (EXPORT_PATCH) {
        patch_index_map->saveImage(_format("%s_patch_map.pfm",geo_bbox->getFileName().c_str()));
        int total_exported_patch_objects = exportPatches(i, information_manager, geo_bbox,"A");
        std::cout << _format("Total patches: %d", information_manager->getNumberOfPatches()) << std::endl;
        std::cout << _format("Total patch objects: %d", total_exported_patch_objects) << std::endl;
    }
    std::cout << _format("Number of patches: %d", information_manager->getNumberOfPatches()) << std::endl;


    //////////////////////////////////////////////////////////////////////////////////////////
    ///PHASE TWO: PATCH CLUSTERING															//
    ///Using the bhattacharrya distance as a metric merge neighbouring similar patches		//
    //////////////////////////////////////////////////////////////////////////////////////////
    ///Find neighbourhood information between the patches
    std::cout << "Starting neighbourhood processing at " << timestamp() << std::endl;
    deriveNeighbourhoodInformation(information_manager, geo_bbox, patch_index_map);
    std::cout << "Finishing neighbourhood processing at " << timestamp() << std::endl;
    delete patch_index_map;
    ///Initialize a surface for each patch
    initializeSurfaces(information_manager);
    std::cout << _format("Initial number of surfaces: %d",information_manager->getNumberOfSurfaces()) << std::endl;
    ///Merge similar SURFACES
    std::cout << "Starting patch clustering at " << timestamp() << std::endl;
    mergeSimilarSurfaces(information_manager);
    std::cout << _format("Number of surfaces after merging: %d",information_manager->getSurfaces().size()) << std::endl;
    std::cout << "Finished patch clustering at " << timestamp() << std::endl;
    ///Merge small SURFACES
    mergeSmallSurfaces(information_manager, geo_bbox);
    std::cout << _format("Final surfaces after merging of small surfaces: %d",information_manager->getSurfaces().size()) << std::endl;
    ///Finalize the surface information i.e. compute the total point indices and points and save it to the InformationManager
    finalizeSurfaceInformation(geo_bbox, information_manager);
    ///Create the surface index map
    Image *surface_index_map = buildSurfaceIndexMap(i, information_manager, geo_bbox, "A");
    if (EXPORT_SURFACE)	{
        int total_exported_surface_objects = exportSurfaces(i, information_manager, geo_bbox, "A");
        std::cout << _format("Total surfaces: %d", information_manager->getNumberOfSurfaces()) << std::endl;
        std::cout << _format("Total surface objects: %d", total_exported_surface_objects) << std::endl;
    }
    surface_index_map->saveImage(_format("%d/surface_index_map.pfm",i));
    surface_index_map->normalize();
    surface_index_map->saveImage(_format("%d/surface_index_map.png",i));
    delete surface_index_map;


    delete information_manager;


    ///Export the entire mesh object
    ///Get the triangulated object corresponding to this geospatial bounding box
    GeometricObject *object = geo_bbox->getObject();
    object->translate(Vector3f(geo_bbox->getCentroid()));
    ge->exportToOBJ(geo_bbox->getFileName().c_str(),object);
    std::cout << _format("Exported geospatial bounding box: %s", _format("%s.obj",geo_bbox->getFileName().c_str()).c_str()) << std::endl;
    delete object;

    geo_bbox->cleanUp();
    delete geo_bbox;
    delete ge;

    return 1;
}



