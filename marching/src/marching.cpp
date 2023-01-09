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

#include "marching.h"
#include "bwdist.h"
#include "distances.h"

#include <string>


int marching_main(const itk::CommandLineArgumentParser::Pointer &parser,
                  const itk::Logger::Pointer &logger){
    int mode = 0;
    parser->GetCommandLineArgument("-mode",mode);
    switch (mode){
        case 1:{
            logger->Info(" Mode = 1 : Generate distance metrics table \n");
            return generateMetricFile(parser, logger);
        }
        case 0: {
            logger->Info("Mode = 0 : Running distance from mito with depth based speed\n");
            return mitoDistanceFunction(parser, logger);
        }
        default: {
            logger->Critical(" Mode = "+std::to_string(mode) + " Unknown \n");
        }
    }
    return EXIT_SUCCESS;
}



