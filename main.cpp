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

#include <iostream>
#include <string>
#include <cstdlib>
#include <itkLogger.h>
#include <itkLoggerBase.h>
#include <itkStdStreamLogOutput.h>

#include "itkCommandLineArgumentParser.h"
#include "utils.h"
#include "filesystem.h"
#include "marching.h"



namespace fs = std::filesystem;

std::string main_helpstring() {
    std::ostringstream ss;
    ss << "Tools For Astrocyte Analysis:\n";
    ss << "$ nano -input <input_file_path> -tool -mode 0\n\n";
    ss << "Options:: \n";
    ss << "========= \n";

    ss << "\t -input\n";
    ss << "\t\t path to input file\n";

    ss << "\t -output\n";
    ss << "\t\t path to output file\n";

    ss << "\t -tool\n";
    ss << "\t\t  If present  will run tool modes instead of submodules\n";

    ss << "\t -module\n";
    ss << "\t\t -module 0 : Computation related to distances  \n";
    ss << "\t\t -module 0 -mode 0 : compute distance from mitochondria to the astrocytic surface and beyond\n";
    ss << "\t\t -module 0 -mode 1 : generate table of various metrics for psd regions. \n";

    ss << "\t -tool -mode n\n";
    ss << "\t\t Stand alone scripts/tools modes. Following tools are available:\n";
    ss << "\t\t -tool -mode 0 : Resample input image to isotropic grid (default)\n";


    ss << "\n\n";
    ss << "Examples : \n";
    ss << "Generate isotropic image\n";
    ss <<  "$./nano -tool -mode 0  -input /<path>/astrocyte.tif\n";

    ss << "Generate Distance Metrics\n";
    ss << "\t Step 0: Compute Distance Functions\n";
    ss <<  "\t$./nano -module 0 -mode 0  -input ./data\n";
    ss << "\t Step 1: Generate table of distance metrics for each PSD\n";
    ss <<  "$./nano -module 0 -mode 1  -input ./data\n";

    std::string helpString( ss.str() );
    return helpString;
}

int main(int argc, char* argv[]){
    //itk::MultiThreaderBase::SetGlobalDefaultNumberOfThreads(1);
    itk::Logger::Pointer logger = itk::Logger::New();
    itk::StdStreamLogOutput::Pointer itkcout = itk::StdStreamLogOutput::New();
    itkcout->SetStream(std::cout);

    itk::StdStreamLogOutput::Pointer itkfout = itk::StdStreamLogOutput::New();
    std::string fname;
    for(int i =0; i < argc; ++i){
        fname +=  std::string(argv[i]);
    }
    std::replace(fname.begin(),fname.end(),'/','#');

    std::ofstream logFile;
    logFile.open("/tmp/" + fname + ".log");
    if(logFile.is_open()) itkfout->SetStream(logFile);
    logger->SetLevelForFlushing(itk::LoggerBaseEnums::PriorityLevel::DEBUG);

    logger->SetPriorityLevel(itk::LoggerBaseEnums::PriorityLevel::DEBUG);

    logger->AddLogOutput(itkcout);
    if(logFile.is_open()) logger->AddLogOutput(itkfout);


    std::string humanReadableFormat = "[%b-%d-%Y, %H:%M:%S]";
    logger->SetHumanReadableFormat(humanReadableFormat);
    logger->SetTimeStampFormat(itk::LoggerBaseEnums::TimeStampFormat::HUMANREADABLE);

    using CommandLineParser = itk::CommandLineArgumentParser;
    CommandLineParser::Pointer parser = CommandLineParser::New();

    parser->SetCommandLineArguments( argc, argv );
    parser->SetProgramHelpText(main_helpstring());

    if (parser->CheckForRequiredArguments() ==
        itk::CommandLineArgumentParser::HELPREQUESTED) {
      return EXIT_SUCCESS;
    }
    if (parser->CheckForRequiredArguments() ==
        itk::CommandLineArgumentParser::FAILED) {
      std::cout << parser->GetProgramHelpText() << std::endl;
      return EXIT_FAILURE;
    }

    bool toolMode = parser->ArgumentExists("-tool");
    int mode = 0,module = 0;
    bool moduleExists = parser->GetCommandLineArgument("-module",module);
    if (!moduleExists && !toolMode){
        logger->Critical( "Please specify -tool or -module option\n");
        return EXIT_FAILURE;
    }
    parser->GetCommandLineArgument("-mode",mode);
    std::string enableDisable  = toolMode? "ENABLED" : "DISABLED";
    logger->Info( "Tools " + enableDisable + "\n");
    logger->Info("Mode Set to "+ std::to_string(mode)  + "\n");

    if(toolMode){
        std::string inputFileName = "data/astrocyte.tif";
        std::string outputFileName;
        if(!parser->GetCommandLineArgument("-input", inputFileName)){
            logger->Warning("Using Default input file : " + inputFileName + "\n");
        }else{
            logger->Info("Input File : " + inputFileName + "\n");
        }
        switch (mode){
            case(0):{
               logger->Info("Running Make isotropic Tool\n");
            }
            default: {
                float spacing = 4.1341146f;
                parser->GetCommandLineArgument("-spacing", spacing);
                logger->Info("Set output spacing to " + std::to_string(spacing) + "\n");
                makeIsotropic(inputFileName, spacing);
            }
        }
    }else {
        switch (module){
            case 0:{
                logger->Info(" Running Module 1 (Distance computation) \n");
                marching_main(parser,logger);
                break;
            }
            default:{
                logger->Critical("Unknown Module\n");
            }
        }
    }
    logFile.close();
    return EXIT_SUCCESS;
}


