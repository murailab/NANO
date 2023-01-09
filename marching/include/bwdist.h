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

#ifndef NANOARCHITECTURE_BWDIST_H
#define NANOARCHITECTURE_BWDIST_H

#include <itkLogger.h>
#include "itkCommandLineArgumentParser.h"
#include "filesystem.h"

namespace fs = std::filesystem;

int mitoDistanceFunction(const itk::CommandLineArgumentParser::Pointer &parser,
                         const itk::Logger::Pointer &logger);

#endif //NANOARCHITECTURE_BWDIST_H
