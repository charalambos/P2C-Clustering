#ifndef __INFORMATION_MANAGER_H__
#define __INFORMATION_MANAGER_H__

#include "Patch.h"
#include "Surface.h"
#include "Utilities.h"

#include <vector>
#include <iterator>

typedef struct PATCH_NEIGHBOUR_INFO	{
	std::vector<int> neighbour_patches;
	std::vector<int> neighbour_point_indices;
} PATCH_NEIGHBOUR_INFO;

typedef struct SURFACE_NEIGHBOUR_INFO	{
	std::vector<Surface *> neighbour_surfaces;
	std::vector<int> neighbour_patch_indices;
} SURFACE_NEIGHBOUR_INFO;

class InformationManager {
	public:
		InformationManager();
		~InformationManager();

		///Returns all the patches
		std::vector<Patch *> &getPatches();

		///Returns all the surfaces
		std::vector<Surface *> &getSurfaces();

		///Returns the patch at the given index
		Patch *getPatch(int index);

		///Returns the surface at the given index
		Surface *getSurface(int index);

		///Returns the surface which contains the patch in question
		Surface *getBelongingSurfaceOf(int index);

		///Returns a list of the neighbours of the patch in question
		std::vector<int> getPatchNeighboursOf(int index);

		///Returns the patches belonging to a surface
		std::vector<int> getPatchesOfSurface(int index);

		///Returns a list of the neighbours of the surface in question
		std::vector<Surface *> getSurfaceNeighboursOf(int index);

		///Initializes the patch neighbour information
		void initializePatchNeighbourInfo();

		///Returns true if a point index exists in the neighbour information of the patch at index
		bool existsNeighbourPointIndex(int patch_index, int point_index);

		///Adds a point index to the patch at index. Checks for duplicates before adding.
		///Returns true if the addition was successfull otherwise if it's already there
		///if returns false
		bool addNeihgbourPointIndex(int patch_index, int point_index);

		///Returns the number of point indices for a patch
		int getNumberOfNeighbourPointIndicesFor(int patch_index);

		///Returns the point at the point index of the patch at index
		int getNeighbourPointIndex(int patch_index, int point_index);

		///Returns a list of the neighbouring point indices of the patch at index
		std::vector<int> &getNeighbourPointIndices(int patch_index);

		///Adds a neighbour patch. Checks for duplicates before adding
		bool addNeighbourPatchIndex(int patch_index, int neighbour_patch_index);

		///Returns true if a patch index exists in the neighbour information of the patch at index
		bool existsNeighbourPatchIndex(int patch_index, int neighbour_patch_index);

		///Returns the number of patch neighbours
		int getNumberOfPatches();

		///Returns the number of surfaces
		int getNumberOfSurfaces();

		///Sorts the patch neighbours based on the patch indices
		void sortPatchNeighbours();

		///Adds a patch
		int addPatch(Patch *patch);

		///Adds a surface
		void addSurface(Surface *surface);

		///Adds a patch to a surface
		void addPatchesToSurface(Surface *surface, std::vector<int> const &patch_indices);

		///Computes neighbour information for the surface
		void updateNeighbours(Surface *surface);

		///Reviews surface neighbours. This should be called once all the information for each surface has been stored.
		void reviewSurfaceNeighbours(int surface_index);
		void reviewSurfaceNeighbours();
		void reviewSurfaceNeighbours(std::vector<Surface *> const &merged_surfaces, Surface *new_surface);
		void reviewSurfaceNeighbours(std::vector<Surface *> unwanted_surfaces);
		void reviewSurfaceNeighbours(std::vector<Patch *> unwanted_patches);



		///Returns the point indices for a surface
		void getPointIndicesForSurface(int surface_index, std::vector<Vector2i> &point_indices);

		///Merges multiple surfaces into a new one and returns that
		Surface *merge(std::vector<Surface *> const &similar_surfaces);

		///Returns the surface neighbours
		std::vector<SURFACE_NEIGHBOUR_INFO> &getSurfaceNeighbours();

		///Returns the patch neighbours
		std::vector<PATCH_NEIGHBOUR_INFO> &getPatchNeighbours(); 

		///Returns the index of the surface in the surfaces array
		int getSurfaceIndex(Surface *surface);

		///Set the height and dot variation map
		void setHeightAndDotVariationMap(Image *_height_and_dot_variation_map);

		///Remove a surface
		void removeSurface(int pos);

		///Remove a patch
		void removePatch(int pos);

	protected:
		///The patches
		std::vector<Patch *> patches;
		///The neighbourhood information for the patches
		std::vector<PATCH_NEIGHBOUR_INFO> patch_neighbours;
		///The surfaces
		std::vector<Surface *> surfaces;
		///The neighbourhood information for the surfaces
		std::vector<SURFACE_NEIGHBOUR_INFO> surface_neighbours;
};


#endif
