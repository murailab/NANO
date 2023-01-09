//**********************************************************
//Copyright Tabish Syed

//This program is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//This program is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with this program.  If not, see <http://www.gnu.org/licenses/>.
//**********************************************************
//
// Created by tabish on 2023-01-10.
//

#ifndef NANOARCHITECTURE_UTILS_H
#define NANOARCHITECTURE_UTILS_H
#include <string>

#include "filesystem.h"

#define BBOX_COL_START   0
#define BBOX_COL_END     1
#define BBOX_ROW_START   2
#define BBOX_ROW_END     3
#define BBOX_SLICE_START   4
#define BBOX_SLICE_END     5
#define VOLUME        6
#define CENTROID_COL   7
#define CENTROID_ROW   8
#define CENTROID_SLICE   9
#define CLOSEST_COL    10
#define CLOSEST_ROW    11
#define CLOSEST_SLICE    12
#define PSD_ASTRO    13
#define SURFACE_COL    14
#define SURFACE_ROW    15
#define SURFACE_SLICE    16
#define MITO_ASTRO   17
#define MITO_ID      18
#define CLUSTER_ID   19

namespace fs = std::filesystem;


int makeIsotropic(const std::string& filename, float isoSpacing = 4.1341146);

size_t csvFileRowCount(std::string metricsFileName);

size_t tic();

#endif //NANOARCHITECTURE_UTILS_H
