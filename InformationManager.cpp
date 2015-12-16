#ifndef __INFORMATION_MANAGER_CPP__
#define __INFORMATION_MANAGER_CPP__


#include "InformationManager.h"

InformationManager::InformationManager()	{
	patches.clear();
	surfaces.clear();
}

InformationManager::~InformationManager()	{
	for (int i=0;i<patches.size();i++)	{
		delete patches[i];
	}
	for (int i=0;i<surfaces.size();i++)	{
		delete surfaces[i];
	}
	patches.clear();
	surfaces.clear();
	for (int i=0;i<patch_neighbours.size();i++)	{
		patch_neighbours[i].neighbour_patches.clear();
		patch_neighbours[i].neighbour_point_indices.clear();
	}
	for (int i=0;i<surface_neighbours.size();i++)	{
		surface_neighbours[i].neighbour_surfaces.clear();
		surface_neighbours[i].neighbour_patch_indices.clear();
	}
	patch_neighbours.clear();
	surface_neighbours.clear();
}

std::vector<Patch *> &InformationManager::getPatches()	{
	return patches;
}

std::vector<Surface *> &InformationManager::getSurfaces()	{
	return surfaces;
}

Patch *InformationManager::getPatch(int index)	{
	if (index >= 0  && index < patches.size())	{
		return patches[index];
	}
	return 0x00;
}

Surface *InformationManager::getSurface(int index)	{
	if (index >= 0 && index < surfaces.size())	{
		return surfaces[index];
	}
	return 0x00;
}

Surface *InformationManager::getBelongingSurfaceOf(int index)	{
	Patch *patch = getPatch(index);
	if (patch)	return patch->belongsTo();
	return 0x00;
}

std::vector<int> InformationManager::getPatchNeighboursOf(int index)	{
	return patch_neighbours[index].neighbour_patches;
}

std::vector<int> InformationManager::getPatchesOfSurface(int index)	{
	return surfaces[index]->getPatches();
}

std::vector<Surface *> InformationManager::getSurfaceNeighboursOf(int index)	{
	return surface_neighbours[index].neighbour_surfaces;
}

void InformationManager::initializePatchNeighbourInfo()	{
	patch_neighbours.clear();
	patch_neighbours.resize(patches.size());
	for (int i=0;i<patches.size();i++)	{
		patch_neighbours[i].neighbour_patches= std::vector<int>();
		patch_neighbours[i].neighbour_point_indices = std::vector<int>();
	}
	return;
}

bool InformationManager::existsNeighbourPointIndex(int patch_index, int point_index)	{
	std::vector<int>::const_iterator it =  std::find(patch_neighbours[patch_index].neighbour_point_indices.begin(),
							 patch_neighbours[patch_index].neighbour_point_indices.end(),
							 point_index);
	if (it == patch_neighbours[patch_index].neighbour_point_indices.end())	{
		return false;
	}
	return true;
}

bool InformationManager::addNeihgbourPointIndex(int patch_index, int point_index)	{
	if (!existsNeighbourPointIndex(patch_index,point_index))	{
		patch_neighbours[patch_index].neighbour_point_indices.push_back(point_index);
		return true;
	}
	return false;
}

int InformationManager::getNumberOfNeighbourPointIndicesFor(int patch_index)	{
	return patch_neighbours[patch_index].neighbour_point_indices.size();
}

std::vector<int> &InformationManager::getNeighbourPointIndices(int patch_index)	{
	return patch_neighbours[patch_index].neighbour_point_indices;
}

int InformationManager::getNeighbourPointIndex(int patch_index, int point_index)	{
	return patch_neighbours[patch_index].neighbour_point_indices[point_index];
}

bool InformationManager::addNeighbourPatchIndex(int patch_index, int neighbour_patch_index)	{
	if (!existsNeighbourPatchIndex(patch_index,neighbour_patch_index))	{
		if (neighbour_patch_index  >= patches.size())	{
			std::cout << "here 3: " << neighbour_patch_index << " " << patch_index << std::endl;
		}
		patch_neighbours[patch_index].neighbour_patches.push_back(neighbour_patch_index);
		return true;
	}
	return false;
}

bool InformationManager::existsNeighbourPatchIndex(int patch_index, int neighbour_patch_index)	{
	std::vector<int>::iterator it =  std::find(patch_neighbours[patch_index].neighbour_patches.begin(),
							 patch_neighbours[patch_index].neighbour_patches.end(),
							 neighbour_patch_index);
	if (it == patch_neighbours[patch_index].neighbour_patches.end())	{
		return false;
	}
	return true;
}

int InformationManager::getNumberOfPatches()	{
	return patches.size();
}

int InformationManager::getNumberOfSurfaces()	{
	return surfaces.size();
}

void InformationManager::sortPatchNeighbours()	{

	///Sort the neighbour indices
	for (int i=0;i<patch_neighbours.size();i++)	{
		std::sort(patch_neighbours[i].neighbour_patches.begin(),patch_neighbours[i].neighbour_patches.end(),compare_func<int>);
	}
	return;
}

int InformationManager::addPatch(Patch *patch)	{
	patches.push_back(patch);
	return patches.size()-1;
}

void InformationManager::addSurface(Surface *surface)	{
	surfaces.push_back(surface);
	SURFACE_NEIGHBOUR_INFO surface_neighbour_info;
	surface_neighbour_info.neighbour_surfaces = std::vector<Surface *>();
	surface_neighbour_info.neighbour_patch_indices = std::vector<int>();
	
	surface_neighbours.push_back(surface_neighbour_info);
// 	int current_size = surface_neighbours.size();
// 	surface_neighbours.resize(current_size + 1);
// 	surface_neighbours[current_size - 1].neighbour_surfaces = std::vector<Surface *>();
// 	surface_neighbours[current_size - 1].neighbour_patch_indices = std::vector<int>();
	
	return;
}

void InformationManager::addPatchesToSurface(Surface *surface, std::vector<int> const &patch_indices)	{

	for (int p=0;p<patch_indices.size();p++)	{
		///Add the patch to the surface
		surface->add(patches[patch_indices[p]]);
	}
		
	///Compute surface neighbour information
	updateNeighbours(surface);

	return;
}

void InformationManager::updateNeighbours(Surface *surface)	{

	///Find the surface's index
	int surface_index = -1;
	for (int i=0;i<surfaces.size();i++)	{
		if (surfaces[i] == surface)	{
			surface_index = i;
			break;
		}
	}
	if (surface_index == -1)	{
		std::cout << "here 2" << std::endl;
	}

	///Get all the patch indices which belong to this surface. the indices are already presorted in the surface
	std::vector<int> surface_patch_indices = surface->getPatches();


// std::cout << "p: ";
// for (int i=0;i<surface_patch_indices.size();i++)	{
// 	std::cout << surface_patch_indices[i] << " ";
// }
// std::cout << std::endl;

	///Group all the neighbours of each of the belonging patches
	std::vector<int> total_neighbouring_patch_indices;
	for (int i=0;i<surface_patch_indices.size();i++)	{
		///Go through all the neighbouring patches of the current patch
		for (int j=0;j<patch_neighbours[surface_patch_indices[i]].neighbour_patches.size();j++)	{
			///Check if it's already part of the total list
			std::vector<int>::const_iterator it = std::find(total_neighbouring_patch_indices.begin(),
					total_neighbouring_patch_indices.end(),
					patch_neighbours[surface_patch_indices[i]].neighbour_patches[j]);
			///If it's not there then add it
			if (it == total_neighbouring_patch_indices.end())	{
				total_neighbouring_patch_indices.push_back(patch_neighbours[surface_patch_indices[i]].neighbour_patches[j]);
			}
		}
	}
	
	///Sort the ids of the total neighbouring patch indices
	std::sort(total_neighbouring_patch_indices.begin(),total_neighbouring_patch_indices.end(),compare_func<int>);	
// std::cout << "s: ";
// for (int i=0;i<total_neighbouring_patch_indices.size();i++)	{
// 	std::cout << total_neighbouring_patch_indices[i] << " ";
// }
// std::cout << std::endl;

	///At this point the list contains the total number of neighbouring patches. This includes
	///patches which are already part of the surface. In order to determine the neighbouring 
	///patches which are NOT part of the surface i.e. they belong to another surface we have 
	///to compute the difference between the contained patches and the total neighbouring patches.
	///i.e. total_neighbouring_patch_indices - surface_patch_indices

	///Create a container for the result which will have the maximum size of differences
	int max_diffs = total_neighbouring_patch_indices.size() + surface_patch_indices.size();
	std::vector<int> diffs(max_diffs);
	///Compute the difference of the two sets
	std::vector<int>::iterator it = std::set_difference (total_neighbouring_patch_indices.begin(), total_neighbouring_patch_indices.end(), surface_patch_indices.begin(), surface_patch_indices.end(), diffs.begin());
	///Check how many differences exist
	int num_diffs = it - diffs.begin();	
	
// std::cout << "Surface neighbours: ";
	///Get the neighbourhood information for the SURFACE
	///1. Keep the neighbouring patch indices of the SURFACE
	///2. Keep the neighbouring surfaces of the SURFACE
	surface_neighbours[surface_index].neighbour_surfaces.clear();
	surface_neighbours[surface_index].neighbour_patch_indices.clear();
	for (int i=0;i<num_diffs;i++)	{
		///The neighbouring patch index
		surface_neighbours[surface_index].neighbour_patch_indices.push_back(diffs[i]);
		///If the neighouring surface is not part of the neighbour_surfaces then add it
		std::vector<Surface *>::const_iterator it = std::find(	surface_neighbours[surface_index].neighbour_surfaces.begin(), 
									surface_neighbours[surface_index].neighbour_surfaces.end(),
									patches[diffs[i]]->belongsTo());
		if (it == surface_neighbours[surface_index].neighbour_surfaces.end())	{
			surface_neighbours[surface_index].neighbour_surfaces.push_back(patches[diffs[i]]->belongsTo());
		}
// 		std::cout << diffs[i] << ", ";
	}
// 	std::cout << std::endl;

	return;
}


void InformationManager::reviewSurfaceNeighbours()	{
	///Go through all the surfaces neighbour information
	for (int i=0;i<surface_neighbours.size();i++)	{
		reviewSurfaceNeighbours(i);
	}

	///Debugging
// 	for (int i=0;i<surface_neighbours.size();i++)	{
// 		std:cout << "Surface " << surfaces[i] << ": ";	
// 		for (int j=0;j<surface_neighbours[i].neighbour_surfaces.size();j++)	{
// 			std::cout << surface_neighbours[i].neighbour_surfaces[j] << ", ";
// 		}
// 		std::cout << std::endl;
// 	}

	return;
}

void InformationManager::reviewSurfaceNeighbours(int surface_index)	{

	///Store the pointers to the neighbouring SURFACES
	surface_neighbours[surface_index].neighbour_surfaces.clear();
	for (int j=0;j<surface_neighbours[surface_index].neighbour_patch_indices.size();j++)	{
		///Get the surface which contains this patch
		Surface *surface = patches[surface_neighbours[surface_index].neighbour_patch_indices[j]]->belongsTo();
		///Make sure it's not already on the list
		std::vector<Surface *>::const_iterator it = std::find(surface_neighbours[surface_index].neighbour_surfaces.begin(),
				surface_neighbours[surface_index].neighbour_surfaces.end(),
				surface);
		///If it's not there then add it
		if (it == surface_neighbours[surface_index].neighbour_surfaces.end())	{
			surface_neighbours[surface_index].neighbour_surfaces.push_back(surface);
		}
	}

	return;
}

void InformationManager::reviewSurfaceNeighbours(std::vector<Surface *> const &merged_surfaces, Surface *new_surface)	{					
	///Go through all the surfaces neighbour information
	for (int i=0;i<surface_neighbours.size();i++)	{
		///first check if the surface does not belong to the merged surfaces
		bool same_surface = false;
		for (int j=0;j<merged_surfaces.size();j++)	{
			if (merged_surfaces[j] == surfaces[i])	{
				same_surface = true;
				break;
			}
		}
		if (same_surface) continue;
	
		///Otherwise check if the neighbours belong to the merged surfaces
		bool had_merged_surface = false;
		
		///Go through the surface neighbours and check if they belong in the merged surfaces
		int variable_size = surface_neighbours[i].neighbour_surfaces.size();
		for (int j=0;j<variable_size;j++)	{
			for (int k=0;k<merged_surfaces.size();k++)	{
				///If the neighbouring surface belongs to the merged surfaces which will be removed then remove it and add the new surface as a neighbour.
				if(surface_neighbours[i].neighbour_surfaces[j] == merged_surfaces[k])	{
					std::vector<Surface *>::iterator sitor = surface_neighbours[i].neighbour_surfaces.begin() + j;
					surface_neighbours[i].neighbour_surfaces.erase(sitor);
					had_merged_surface = true;
					j--;
					variable_size--;
					break;
				}
			}
		}
		
		if (had_merged_surface && new_surface != 0x00)	surface_neighbours[i].neighbour_surfaces.push_back(new_surface);
	}

	///Debugging
// 	for (int i=0;i<surface_neighbours.size();i++)	{
// 		std:cout << "Surface " << surfaces[i] << ": ";	
// 		for (int j=0;j<surface_neighbours[i].neighbour_surfaces.size();j++)	{
// 			std::cout << surface_neighbours[i].neighbour_surfaces[j] << ", ";
// 		}
// 		std::cout << std::endl;
// 	}

	return;
}

void InformationManager::reviewSurfaceNeighbours(std::vector<Surface *> unwanted_surfaces)	{
	///Get all the patch indices contained in the unwanted surfaces
	std::vector<int> unwanted_patch_indices;
	for (int  i=0;i<unwanted_surfaces.size();i++)	{
		///Get the patch indices
		std::vector<int> patch_indices = unwanted_surfaces[i]->getPatches();
		///Add them to the total list
		unwanted_patch_indices.insert(unwanted_patch_indices.end(),patch_indices.begin(), patch_indices.end());
	}

	///At this point unwanted_patch_indices contains all the patches contained in the unwanted_surfaces

	///First go through the patches and remove all the unwanted ones
	for (int i=0;i<unwanted_patch_indices.size();i++)	{
		std::vector<Patch *>::iterator p_itor = patches.begin();
		int variable_size = patches.size();
		for (int j=0;j<variable_size;j++, p_itor++)	{
			if (patches[j]->getClusterId() == unwanted_patch_indices[i])	{
				delete patches[j];
				patches.erase(p_itor);
				break;
			}
		}
	}

	///Then remove the unwanted surfaces
	std::vector<Surface *>::iterator s_itor = surfaces.begin();
	int variable_size = surfaces.size();
	for (int i=0;i<variable_size;i++,s_itor++)	{
		for (int j=0;j<unwanted_surfaces.size();j++)	{
				if (surfaces[i] == unwanted_surfaces[j])	{
					delete surfaces[i];
					surfaces.erase(s_itor);
					variable_size--;
					s_itor = surfaces.begin() + i -1;
					i-=1;
					break;
				}
		}
	}

	///Then go through the surfaces and check if any of them contained the surfaces as neighbours
	for (int i=0;i<surfaces.size();i++)	{
		std::vector<int> patch_indices = surfaces[i]->getPatches();

	}

	return;
}

void InformationManager::reviewSurfaceNeighbours(std::vector<Patch *> unwanted_patches)	{

	return;
}


void InformationManager::getPointIndicesForSurface(int surface_index, std::vector<Vector2i> &point_indices)	{
	std::vector<int> patch_indices = surfaces[surface_index]->getPatches();
	for (int i=0;i<patch_indices.size();i++)	{
		std::vector<Vector2i> patch_point_indices = patches[patch_indices[i]]->getPointIndices();
		for (int j=0;j<patch_point_indices.size();j++)	{
			point_indices.push_back(patch_point_indices[j]);
		}
	}
	
	return;
}

Surface *InformationManager::merge(std::vector<Surface *> const &similar_surfaces)	{
	///Add all the data of the similar surfaces
	std::vector<int> total_patches;
	for (int i=0;i<similar_surfaces.size();i++)	{
		///Get the patches of the surface
		std::vector<int> patches = similar_surfaces[i]->getPatches();
		for (int j=0;j<patches.size();j++)	{
			std::vector<int>::const_iterator it =  std::find(total_patches.begin(),
								 	 total_patches.end(),
								 	 patches[j]);	
			///and add it to a global list if it's not already there.
			if (it == total_patches.end())	{
				total_patches.push_back(patches[j]);
			}
		}
	}
	///sort the patches based on their id i.e. ascending order
	std::sort(total_patches.begin(), total_patches.end(), compare_func<int>);

	///Create a new surface
	Surface *new_surface = new Surface(this);
	///Add the surface to the list of surfaces
	addSurface(new_surface);
	///Add them to the new surface
	addPatchesToSurface(new_surface,total_patches);
	
	return new_surface;
}

std::vector<SURFACE_NEIGHBOUR_INFO> &InformationManager::getSurfaceNeighbours()	{
	return surface_neighbours;
}

std::vector<PATCH_NEIGHBOUR_INFO> &InformationManager::getPatchNeighbours()	{
	return patch_neighbours;
}

int InformationManager::getSurfaceIndex(Surface *surface)	{
	for (int i=0;i<surfaces.size();i++)	{
		if (surfaces[i] == surface)	{
			return i;
		}
	}
	std::cout << "Negative surface index" << std::endl;
	return -1;
}

void InformationManager::removeSurface(int pos)	{
	std::vector<Surface *>::iterator s_itor = surfaces.begin() + pos;
	surfaces.erase(s_itor);
	std::vector<SURFACE_NEIGHBOUR_INFO>::iterator sn_itor = surface_neighbours.begin() + pos;
	surface_neighbours.erase(sn_itor);
	return;
}

void InformationManager::removePatch(int pos)	{
	std::vector<Patch *>::iterator p_itor = patches.begin() + pos;
	patches.erase(p_itor);
	return;
}

#endif
