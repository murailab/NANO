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

#include <string>

#include <itkFastMarchingImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImage.h>
#include <itkCastImageFilter.h>
#include <itkChangeInformationImageFilter.h>
#include <itkTIFFImageIO.h>

#include "bwdist.h"

namespace {
    //Fix this hideousness!
    //these global variables are redefined in other files e.g distances.cpp..
    using InternalPixelType = float;
    constexpr unsigned int Dimension = 3;
    using InternalImageType = itk::Image<InternalPixelType, Dimension>;
    using OutputPixelType = float;
    using OutputImageType = itk::Image<OutputPixelType, Dimension>;
    using MaskPixelType = uint8_t;
    using MaskImageType = itk::Image<MaskPixelType, Dimension>;

    MaskImageType::SpacingType spacing;
}


int runOutside(const itk::CommandLineArgumentParser::Pointer &parser,
               const itk::Logger::Pointer &logger){
    std::string root;
    parser->GetCommandLineArgument("-input",root);

    std::string seedImageName = root + "astrocyte.tif";

    logger->Info("Set astrocyte Image File to " + seedImageName + "\n");
    if(!fs::exists(seedImageName)){
        logger->Critical("File " + seedImageName + " not Found!\n");
        return EXIT_FAILURE;
    }
    using MaskReaderType = itk::ImageFileReader< MaskImageType >;
    
    MaskReaderType::Pointer astroReader = MaskReaderType::New();
    
    
    astroReader->SetFileName( seedImageName );    

    using InformationChangeFilterType = itk::ChangeInformationImageFilter< MaskImageType >;
    InformationChangeFilterType::Pointer infoChangeFilter = InformationChangeFilterType::New();
    infoChangeFilter->SetInput( astroReader->GetOutput() );


    infoChangeFilter->SetOutputSpacing( spacing );
    infoChangeFilter->ChangeSpacingOn();
    MaskImageType::Pointer seedImage = infoChangeFilter->GetOutput();
    seedImage->Update();
    
    MaskImageType::RegionType region = seedImage->GetLargestPossibleRegion();

    InternalImageType::SizeType size = region.GetSize();
    
    using FastMarchingFilterType = itk::FastMarchingImageFilter< InternalImageType, InternalImageType >;
    FastMarchingFilterType::Pointer  fastMarching = FastMarchingFilterType::New();
    
    using NodeContainer = FastMarchingFilterType::NodeContainer;
    using NodeType = FastMarchingFilterType::NodeType;
    NodeContainer::Pointer seeds = NodeContainer::New();
    
    using IteratorType = itk::ImageRegionIteratorWithIndex< MaskImageType >;
    IteratorType  it( seedImage, seedImage->GetLargestPossibleRegion() );
    int i = 0;
    float initialDistance = 0.0;
    it.GoToBegin();
    seeds->Initialize();
    while( !it.IsAtEnd() ){
        if(it.Get() > 0){ 
            NodeType node;
            node.SetValue(initialDistance);
            node.SetIndex(it.GetIndex());
            seeds->InsertElement(i++, node);
        }
        ++it;
    }
   
    fastMarching->SetTrialPoints(seeds);
    fastMarching->SetSpeedConstant(1.0);
    fastMarching->SetOutputSize(size);
    fastMarching->SetOutputSpacing(seedImage->GetSpacing());

    
    using CastOutputFilterType = itk::CastImageFilter< InternalImageType, OutputImageType >;
    CastOutputFilterType::Pointer outputCaster = CastOutputFilterType::New();
    outputCaster->SetInput(fastMarching->GetOutput());

    std::string outputfilename = fs::path(root) / "astroToa.tif";
    using OutputWriterType = itk::ImageFileWriter<  OutputImageType  >;
    OutputWriterType::Pointer writer = OutputWriterType::New();
    writer->SetInput(outputCaster->GetOutput());
    writer->SetFileName(outputfilename);
    try {
        writer->Update();
        logger->Info("Wrote time of arrival to " + outputfilename + "\n");
    }catch (itk::ExceptionObject &e){
        logger->Critical("Failed to write time of arrival image\n");
        return EXIT_FAILURE;
    }
    return  EXIT_SUCCESS;
}

int runInside(const itk::CommandLineArgumentParser::Pointer &parser,
               const itk::Logger::Pointer &logger){
    std::string root;
    parser->GetCommandLineArgument("-input",root);
    
    std::string speedImageName = root + "astrocyte.tif";
    std::string seedImageName = root + "mito.tif";

    logger->Info("Set astrocyte Image File to " + speedImageName + "\n");
    if(!fs::exists(speedImageName)){
        logger->Critical("File " + speedImageName + " not Found!\n");
        return EXIT_FAILURE;
    }

    logger->Info("Set Mitochondria Image File to " + seedImageName + "\n");
    if(!fs::exists(seedImageName)){
        logger->Critical("File " + seedImageName + " not Found!\n");
        return EXIT_FAILURE;
    }

    using MaskReaderType = itk::ImageFileReader< MaskImageType >;
    
    MaskReaderType::Pointer mitoReader = MaskReaderType::New();    
    MaskReaderType::Pointer astroReader = MaskReaderType::New();
    
    
    astroReader->SetFileName( speedImageName );    
    mitoReader->SetFileName( seedImageName );

    using InformationChangeFilterType = itk::ChangeInformationImageFilter< MaskImageType >;
    InformationChangeFilterType::Pointer infoChangeFilter = InformationChangeFilterType::New();
    infoChangeFilter->SetInput( astroReader->GetOutput() );

    infoChangeFilter->SetOutputSpacing( spacing );
    infoChangeFilter->ChangeSpacingOn();
    

    using ThresholdingFilterType = itk::BinaryThresholdImageFilter< MaskImageType, MaskImageType >;
    ThresholdingFilterType::Pointer thresholder = ThresholdingFilterType::New();
    
    thresholder->SetLowerThreshold( 1 );
    thresholder->SetUpperThreshold( 100 );
    thresholder->SetOutsideValue(  0  );
    thresholder->SetInsideValue(  255 );
    thresholder->SetInput(mitoReader->GetOutput());
  
    
    MaskImageType::Pointer seedImage = thresholder->GetOutput();
    MaskImageType::Pointer speedImage = infoChangeFilter->GetOutput();
    speedImage->Update();
    seedImage->Update();

    using FastMarchingFilterType = itk::FastMarchingImageFilter< InternalImageType, MaskImageType >;
    FastMarchingFilterType::Pointer  fastMarching = FastMarchingFilterType::New();
    
    using NodeContainer = FastMarchingFilterType::NodeContainer;
    using NodeType = FastMarchingFilterType::NodeType;
    NodeContainer::Pointer seeds = NodeContainer::New();
    
    using IteratorType = itk::ImageRegionIteratorWithIndex< MaskImageType >;
    IteratorType  it( seedImage, seedImage->GetLargestPossibleRegion() );
    int i = 0;
    float initialDistance = 0.0;
    it.GoToBegin();
    seeds->Initialize();
    while( !it.IsAtEnd() ){
        if(it.Get() > 0){ 
            NodeType node;
            node.SetValue(initialDistance);
            node.SetIndex(it.GetIndex());
            seeds->InsertElement(i++, node);
        }
        ++it;
    }
    
    fastMarching->SetTrialPoints(seeds);
    fastMarching->SetInput(speedImage);
    
    using CastOutputFilterType = itk::CastImageFilter< InternalImageType, OutputImageType >;
    CastOutputFilterType::Pointer outputCaster = CastOutputFilterType::New();
    outputCaster->SetInput(fastMarching->GetOutput());
    
    std::string outputfilename = fs::path(root) / "toa.tif";
    using OutputWriterType = itk::ImageFileWriter<  OutputImageType  >;
    OutputWriterType::Pointer writer = OutputWriterType::New();
    using TIFFIOType = itk::TIFFImageIO ;
    TIFFIOType::Pointer tiffIO = TIFFIOType::New();
    writer->SetImageIO(tiffIO);
    writer->SetInput(outputCaster->GetOutput());
    writer->SetFileName(outputfilename);
    try {
        writer->Update();
        logger->Info("Wrote time of arrival to " + outputfilename + "\n");
    }catch (itk::ExceptionObject &e){
       logger->Critical("Failed to write time of arrival image\n");
       return EXIT_FAILURE;
    }
    return  EXIT_SUCCESS;
}

int mitoDistanceFunction(const itk::CommandLineArgumentParser::Pointer &parser,
                         const itk::Logger::Pointer &logger){

    spacing[0] = 4.1341146f;
    spacing[1] = 4.1341146f;
    spacing[2] = 8;
    std::string root = "./data/";
    parser->GetCommandLineArgument("-input",root);
    if( !fs::exists(root) ){
        logger->Critical("Could not find folder " + root + "\n");
        return EXIT_FAILURE;
    }
    std::string astroToaImageName = fs::path(root) / "astroToa.tif";

    if( !fs::exists(astroToaImageName) ){
        logger->Info("running outside\n");
        runOutside(parser, logger);
    }else{
      logger->Info(astroToaImageName + " exists : Skipping Outside\n");
    }

    runInside(parser, logger);
    return EXIT_SUCCESS;
}
