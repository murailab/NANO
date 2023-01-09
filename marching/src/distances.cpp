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

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkVariableSizeMatrix.h>
#include <itkMyCSVNumericObjectFileWriter.h>
#include <itkCSVArray2DFileReader.h>
#include <itkChangeInformationImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkMultiplyImageFilter.h>
#include <itkMinimumMaximumImageCalculator.h>

#include "itkLabelGeometryImageFilter.h"
#include "distances.h"
#include "utils.h"

namespace {
    //Fix this hideousness!
    //these global variables are redefined in other files e.g bwdist.cpp
    constexpr unsigned int Dimension = 3;
    using OutputPixelType = float;
    using LabelPixelType = uint16_t;
    using MitoPixelType = uint8_t;
    using OutputImageType = itk::Image<OutputPixelType, Dimension>;
    using LabelImageType = itk::Image<LabelPixelType, Dimension>;
    using MitoImageType = itk::Image<MitoPixelType, Dimension>;
    using DataMatrixType = itk::VariableSizeMatrix<OutputPixelType>;
    using MaskPixelType = uint8_t;
    using MaskImageType = itk::Image<MaskPixelType, Dimension>;

    MaskImageType::SpacingType spacing;
}

void printMatrix(DataMatrixType* matrix) {
    for (size_t i = 0; i < matrix->Rows(); ++i) {
        for (size_t j = 0; j < matrix->Cols(); ++j) {
            std::cout << (*matrix)(i, j) << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

LabelImageType::SizeValueType psdRegionCount(std::string root){
    std::cout<<"Counting number of psd regions"<<std::endl;
    std::string psdFileName = root + "psd16.tif";

    if(!fs::exists(psdFileName)){
        std::cerr<<"File "<< psdFileName<<" Not Found \n";
        return -1;//return a very large number
    }

    using LabelImageReaderType = itk::ImageFileReader< LabelImageType >;
    LabelImageReaderType::Pointer psdImageReader = LabelImageReaderType::New();
    psdImageReader->SetFileName(psdFileName);

    LabelImageType::Pointer psdLabelImage = psdImageReader->GetOutput();
    psdLabelImage->Update();

    using LabelGeometryImageFilterType = itk::LabelGeometryImageFilter< LabelImageType >;
    LabelGeometryImageFilterType::Pointer labelGeometryImageFilter = LabelGeometryImageFilterType::New();
    labelGeometryImageFilter->SetInput( psdLabelImage );
    labelGeometryImageFilter->Update();

    LabelImageType::SizeValueType psdCount = labelGeometryImageFilter->GetNumberOfObjects();
    return psdCount-1;
}


int computePsdMetrics(const itk::CommandLineArgumentParser::Pointer &parser,
                      const itk::Logger::Pointer &logger, DataMatrixType* csvData){
    /* calulates all metrics for regions outside the astrocye
     * i.e metrics for psd regions
     */
    std::string root;
    parser->GetCommandLineArgument("-input", root);
    logger->Info("Starting PSD metrics computation for :" + root + "\n");

    std::string psdFileName = root + "psd16.tif";
    std::string toaFileName = root + "astroToa.tif";
    if(!fs::exists(psdFileName)){
        logger->Error("File " + psdFileName + " Not Found \n");
        return EXIT_FAILURE;
    }
    if(!fs::exists(toaFileName)){
        logger->Critical("File" + toaFileName + " not found\n");
        logger->Critical("Run distance function computation outside astrocyte to generate the file\n");
        return EXIT_FAILURE;
    }


    using LabelImageReaderType = itk::ImageFileReader< LabelImageType >;
    using ToaImageReaderType = itk::ImageFileReader< OutputImageType >;

    LabelImageReaderType::Pointer psdImageReader = LabelImageReaderType::New();
    psdImageReader->SetFileName(psdFileName);

    using InformationChangeFilterType = itk::ChangeInformationImageFilter< LabelImageType >;
    InformationChangeFilterType::Pointer infoChangeFilter = InformationChangeFilterType::New();
    infoChangeFilter->SetInput( psdImageReader->GetOutput() );

    infoChangeFilter->SetOutputSpacing( spacing );
    infoChangeFilter->ChangeSpacingOn();

    ToaImageReaderType::Pointer toaReader = ToaImageReaderType::New();
    toaReader->SetFileName(toaFileName);

    LabelImageType::Pointer psdLabelImage = infoChangeFilter->GetOutput();
    psdLabelImage->Update();

    using InformationChangeFilterTypeFloat = itk::ChangeInformationImageFilter< OutputImageType >;
    InformationChangeFilterTypeFloat::Pointer infoChangeFilterToa = InformationChangeFilterTypeFloat::New();
    infoChangeFilterToa->SetInput( toaReader->GetOutput() );

    infoChangeFilterToa->SetOutputSpacing( spacing );
    infoChangeFilterToa->ChangeSpacingOn();

    OutputImageType::Pointer arrivalImage = infoChangeFilterToa->GetOutput();

    using BinaryThresholdingFilterType = itk::BinaryThresholdImageFilter< LabelImageType, OutputImageType >;
    BinaryThresholdingFilterType::Pointer thresholder = BinaryThresholdingFilterType::New();
    thresholder->SetLowerThreshold( 1 );
    thresholder->SetUpperThreshold( 4098 );
    thresholder->SetOutsideValue( 1e10);
    thresholder->SetInsideValue(  1 );
    thresholder->SetInput(psdLabelImage);
    logger->Debug("Generating PSD mask from labelled PSD image\n");

    using MultiplyImageFilterType = itk::MultiplyImageFilter<BinaryThresholdingFilterType::OutputImageType ,OutputImageType,OutputImageType>;
    MultiplyImageFilterType::Pointer multiplyFilter = MultiplyImageFilterType::New();
    multiplyFilter->SetInput1(thresholder->GetOutput());
    multiplyFilter->SetInput2(arrivalImage);

    logger->Debug("Generate PSD arrival image\n");
    OutputImageType::Pointer psdToaImage = multiplyFilter->GetOutput();
    psdToaImage->Update();
    thresholder->ReleaseDataFlagOn();


    using LabelGeometryImageFilterType = itk::LabelGeometryImageFilter< LabelImageType >;
    LabelGeometryImageFilterType::Pointer labelGeometryImageFilter = LabelGeometryImageFilterType::New();
    labelGeometryImageFilter->SetInput( psdLabelImage );
    labelGeometryImageFilter->Update();
    logger->Debug("Region props computation completed...\n");
    LabelGeometryImageFilterType::LabelsType allLabels = labelGeometryImageFilter->GetLabels();

    using MinMaxCalculatorType = itk::MinimumMaximumImageCalculator< OutputImageType >;
    MinMaxCalculatorType::Pointer minMaxFilter = MinMaxCalculatorType::New();
    minMaxFilter->SetImage( psdToaImage );
    psdLabelImage->ReleaseDataFlagOn();
    if(parser->ArgumentExists("-debug")) {
        logger->Debug("Saving psdToa image\n");
        using WriterType = itk::ImageFileWriter<OutputImageType>;
        WriterType::Pointer writerpsdtoa = WriterType::New();
        writerpsdtoa->SetFileName(root + "psdToa.tif");
        writerpsdtoa->SetInput(psdToaImage);
        try {
            writerpsdtoa->Update();
            logger->Debug("Wrote psdToa file to disk\n");
        }catch(itk::ExceptionObject &e) {
            logger->Critical("Failed to write psdToa image\n");
        }
    }

    std::stringstream ss;
    //Compute outside Metrics.....
    for( size_t psd = 1; psd < allLabels.size(); psd++ ){
        logger->Debug("Geometry analysis of psd " + std::to_string(psd)+"/"+std::to_string(allLabels.size()) + "\n");
        LabelGeometryImageFilterType::LabelPixelType labelValue = allLabels[psd];
        ss.str(""); // initiallize debug string stream

        // Bounding Box of psd region (Bbox xstart,xend, ystart, yend, zstart,zend)
        auto boundingBox = labelGeometryImageFilter->GetBoundingBox(labelValue);
        (*csvData)(psd-1, BBOX_COL_START) =   boundingBox[0];
        (*csvData)(psd-1, BBOX_COL_END)   =   boundingBox[1];
        (*csvData)(psd-1, BBOX_ROW_START) =   boundingBox[2];
        (*csvData)(psd-1, BBOX_ROW_END)   =   boundingBox[3];
        (*csvData)(psd-1, BBOX_SLICE_START) =   boundingBox[4];
        (*csvData)(psd-1, BBOX_SLICE_END)   =   boundingBox[5];

        // Volume of psd region (volume)
        auto volume = labelGeometryImageFilter->GetVolume(labelValue);
        (*csvData)(psd-1,VOLUME) = volume;

        // Location of centroid of psd region (centroid x,y,z)
        auto centroid = labelGeometryImageFilter->GetCentroid(labelValue);
        (*csvData)(psd-1, CENTROID_COL) = centroid[0];
        (*csvData)(psd-1, CENTROID_ROW) = centroid[1];
        (*csvData)(psd-1, CENTROID_SLICE) = centroid[2];

        // Location of closest psd point (closest x,y,z)
        auto psdRegion = labelGeometryImageFilter->GetRegion(labelValue);

        minMaxFilter->SetRegion(psdRegion);
        minMaxFilter->ComputeMinimum();
        auto closestPoint = minMaxFilter->GetIndexOfMinimum();
        (*csvData)(psd-1, CLOSEST_COL) = closestPoint[0];
        (*csvData)(psd-1, CLOSEST_ROW) = closestPoint[1];
        (*csvData)(psd-1, CLOSEST_SLICE) = closestPoint[2];

        // Distance of closest point (psd->astro)
        (*csvData)(psd-1, PSD_ASTRO) = psdToaImage->GetPixel(minMaxFilter->GetIndexOfMinimum());

        ss<< " @" << psdRegion.GetIndex() << " of size " << psdRegion.GetSize();
        ss << " " <<  (*csvData)(psd-1, PSD_ASTRO) << " nm away \n";
        //ss << "coz value of region @ " << minMaxFilter->GetIndexOfMinimum() << " is ";
        //ss << psdToaImage->GetPixel(minMaxFilter->GetIndexOfMinimum()) <<std::endl;
        logger->Debug(ss.str());
    }

    return EXIT_SUCCESS;
}

int estimateAstocyteSurfacePoint(const itk::CommandLineArgumentParser::Pointer &parser,
                                 const itk::Logger::Pointer &logger, DataMatrixType* csvData){
    logger->Info("Estimating Astro surface points\n");

    std::string root;
    parser->GetCommandLineArgument("-input", root);
    std::string maskImageName = root + "astrocyte.tif";

    if(!fs::exists(maskImageName)){
        logger->Error("File "+maskImageName + " not Found \n");
        return EXIT_FAILURE;
    }

    using MaskReaderType = itk::ImageFileReader< MaskImageType >;
    MaskReaderType::Pointer astroReader = MaskReaderType::New();
    astroReader->SetFileName( maskImageName );

    using InformationChangeFilterType = itk::ChangeInformationImageFilter< MaskImageType >;
    InformationChangeFilterType::Pointer infoChangeFilter = InformationChangeFilterType::New();
    infoChangeFilter->SetInput( astroReader->GetOutput() );

    infoChangeFilter->SetOutputSpacing( spacing );
    infoChangeFilter->ChangeSpacingOn();

    MaskImageType::Pointer astrocyte = infoChangeFilter->GetOutput();
    astrocyte->Update();
    MaskImageType::SizeType size = astrocyte->GetLargestPossibleRegion().GetSize();
    MaskImageType::IndexType bound;
    bound[0] = size[0];bound[1] = size[1]; bound[2] = size[2];
    std::stringstream ss;
    for(size_t psd = 0 ; psd < csvData->Rows(); ++psd ){
        ss.str("");
        MaskImageType::PointType psdPoint;
        psdPoint[0]=  (*csvData)(psd, CLOSEST_COL);
        psdPoint[1]=  (*csvData)(psd, CLOSEST_ROW);
        psdPoint[2]=  (*csvData)(psd, CLOSEST_SLICE);
        ss<<"Searching surface point for psd "<<psd+1<<"/"<<csvData->Rows()<<" @"<<psdPoint<<std::endl;
        logger->Debug(ss.str());ss.str("");

        bool found = false;
        int radius = std::min(std::min(spacing[0], spacing[1]), spacing[2]) + 1;
        int oldRadius = 0;
        while(!found){
            MaskImageType::IndexType index;
            for(int x = -(radius/spacing[0]); x <= (radius/spacing[0]) ; ++x ){
                for(int y = -(radius/spacing[1]); y <= (radius/spacing[1]) ; ++y ){
                    for(int z = -(radius/spacing[2]); z <= (radius/spacing[2]) ; ++z ){
                        //don't check points checked in previous iterations...
                        if( ((x*spacing[0]*spacing[0]*x + y*spacing[1]*spacing[1]*y + z*spacing[2]*spacing[2]*z) < radius*radius)
                            && ((x*spacing[0]*spacing[0]*x + y*spacing[1]*spacing[1]*y + z*spacing[2]*spacing[2]*z) >= oldRadius*oldRadius) ){
                            index[0] = std::max(0.0, std::floor(psdPoint[0] + x));
                            index[1] = std::max(0.0, std::floor(psdPoint[1] + y));
                            index[2] = std::max(0.0, std::floor(psdPoint[2] + z));
                            index[0] = std::min(index[0], bound[0]);
                            index[1] = std::min(index[1], bound[1]);
                            index[2] = std::min(index[2], bound[2]);
                            if( astrocyte->GetPixel(index) > 0 ){
                                (*csvData)(psd, SURFACE_COL) = index[0];
                                (*csvData)(psd, SURFACE_ROW) = index[1];
                                (*csvData)(psd, SURFACE_SLICE) = index[2];
                                found = true;
                                break;
                            }
                        }
                    }
                }
            }
            if (radius > 1000){
                logger->Warning("Giving Up search for astrocyte point for psd "+  std::to_string(psd) + "\n");
                break;
            }
            if(!found){
                oldRadius = radius;
                radius += 1;
                ss << "Increasing search radius for psd " << psd << " to (" << radius << ")\n";
            }else{
                ss<<"Found astrocyte in sphere of radius ("<<radius<<") from psd "<<psd<<std::endl;
                (*csvData)(psd,PSD_ASTRO) = radius;
            }
            //logger->Debug(ss.str());
        }
    }
    return EXIT_SUCCESS;
}


OutputImageType::IndexType getAstroIndexFromToa(OutputImageType::Pointer arrivalImage, OutputImageType::IndexType point){
    /*Return index of minimum valued pixel in a nbd*/
    using MinMaxCalculatorType = itk::MinimumMaximumImageCalculator< OutputImageType >;
    MinMaxCalculatorType::Pointer minMaxFilter = MinMaxCalculatorType::New();
    minMaxFilter->SetImage( arrivalImage );

    OutputImageType::SizeType size = arrivalImage->GetLargestPossibleRegion().GetSize();
    OutputImageType::IndexType bound;
    for(size_t d = 0; d < Dimension; ++d)
        bound[d] = size[d] - 1;
    int width = 1;
    bool flag = true;
    while(flag) {
        OutputImageType::RegionType neighbourhood;
        OutputImageType::SizeType nbdSize;
        OutputImageType::IndexType index;
        for (size_t i = 0; i < Dimension; ++i) {
            index[i] = std::max(0.0, std::floor(point[i] - width/2));
            index[i] = std::min(bound[i] - width, index[i]);
            nbdSize[i] = width;
        }
        neighbourhood.SetIndex(index);
        neighbourhood.SetSize(nbdSize);

        minMaxFilter->SetRegion(neighbourhood);
        minMaxFilter->Compute();
        flag = false;
        //GetMinimum behaves is weird!
        //if (minMaxFilter->GetMinimum() > 999999999 && width < 9){
        if((arrivalImage->GetPixel(minMaxFilter->GetIndexOfMinimum()) > 99999999) && (width < 19) ){//9
            width += 2;
            flag = true;
        }
    }
    //std::cout<<" in box of width = " << width;
    OutputImageType::IndexType location = minMaxFilter->GetIndexOfMinimum();
    return location;
}

int mitoToSurfaceDistance(const itk::CommandLineArgumentParser::Pointer &parser,
                          const itk::Logger::Pointer &logger, DataMatrixType *csvData){
    logger->Info("Starting mito metric analysis\n");

    std::string root;
    parser->GetCommandLineArgument("-input", root);
    std::string mitoFileName = root + "mito.tif";
    std::string toaFileName = root + "toa.tif";

    if(!fs::exists(mitoFileName)){
       logger->Error("File "+mitoFileName+" not found \n");
        return EXIT_FAILURE;
    }

    if(!fs::exists(toaFileName)){
        logger->Error("File "+toaFileName+" not found \n");
        return EXIT_FAILURE;
    }

    using MitoImageReaderType = itk::ImageFileReader< MitoImageType >;
    MitoImageReaderType::Pointer mitoImageReader = MitoImageReaderType::New();
    mitoImageReader->SetFileName(mitoFileName);

    using InformationChangeFilterType = itk::ChangeInformationImageFilter< MitoImageType >;
    InformationChangeFilterType::Pointer infoChangeFilterMito = InformationChangeFilterType::New();
    infoChangeFilterMito->SetInput( mitoImageReader->GetOutput() );

    infoChangeFilterMito->SetOutputSpacing( spacing );
    infoChangeFilterMito->ChangeSpacingOn();


    using ToaImageReaderType = itk::ImageFileReader< OutputImageType >;
    ToaImageReaderType::Pointer toaReader = ToaImageReaderType::New();
    toaReader->SetFileName(toaFileName);

    using InformationChangeFilterTypeToa = itk::ChangeInformationImageFilter< OutputImageType >;
    InformationChangeFilterTypeToa::Pointer infoChangeFilterToa = InformationChangeFilterTypeToa::New();
    infoChangeFilterToa->SetInput( toaReader->GetOutput() );

    infoChangeFilterToa->SetOutputSpacing( spacing );
    infoChangeFilterToa->ChangeSpacingOn();

    MitoImageType::Pointer mitoImage = infoChangeFilterMito->GetOutput();
    mitoImage->Update();

    OutputImageType::Pointer arrivalImage = infoChangeFilterToa->GetOutput();
    arrivalImage->Update();
    std::stringstream ss;
    for(size_t r = 0; r < csvData->Rows(); ++r){
        ss.str("");
        OutputImageType::IndexType approxLandingPoint;
        approxLandingPoint[0] = ((*csvData)(r, SURFACE_COL));
        approxLandingPoint[1] = ((*csvData)(r, SURFACE_ROW));
        approxLandingPoint[2] = ((*csvData)(r, SURFACE_SLICE));
        ss<<"Landing Pt -> "<<approxLandingPoint;
        auto index = getAstroIndexFromToa(arrivalImage, approxLandingPoint);
        (*csvData)(r,MITO_ASTRO) = arrivalImage->GetPixel(index);
        ss << " is @ distance = " << (*csvData)(r,MITO_ASTRO)<< std::endl;
        logger->Debug(ss.str());
    }
    return 0;
}


int weightedMetricAnalysis(const itk::CommandLineArgumentParser::Pointer &parser,
                           const itk::Logger::Pointer &logger, DataMatrixType* csvData){
    /*Perform weight Metric analysis using arrival function outside astrocyte
     */
    bool success = true;
    if(parser->ArgumentExists("-debug")) printMatrix(csvData);
    //---------------------------------------------------------------
    if(computePsdMetrics(parser, logger, csvData) == EXIT_FAILURE){
        logger->Warning("ComputepsdMetrics failed!\n");
        success = false;
    }
    logger->Debug("Finished Computing PSD metrics\n");
    if(parser->ArgumentExists("-debug")) printMatrix(csvData);
    //-------------------------------------------------------------

    if(estimateAstocyteSurfacePoint(parser, logger, csvData) == EXIT_FAILURE){
        logger->Warning("estimateAstrocyteSurfacePoint failed!\n");
        success = false;
    }
    logger->Debug("Finished Computing Astrocytic Surface Points\n");
    if(parser->ArgumentExists("-debug")) printMatrix(csvData);
    //-------------------------------------------------------------

    if(mitoToSurfaceDistance(parser, logger, csvData) == EXIT_FAILURE){
        logger->Warning("mitoToSurfaceDistance failed!\n");
        success = false;
    }
    logger->Debug("Finished Computing Mito Distances\n");
    if(parser->ArgumentExists("-debug")) printMatrix(csvData);
    //------------------------------------------------------------

    return success ? EXIT_SUCCESS: EXIT_FAILURE;
}

int generateMetricFile(const itk::CommandLineArgumentParser::Pointer &parser,
                       const itk::Logger::Pointer &logger){

    std::string root;
    parser->GetCommandLineArgument("-input", root);

    logger->Info("Set Input folder to " + root + "\n");
    if(!fs::exists(root)){
        logger->Critical("Folder " + root + " not Found!\n");
        return EXIT_FAILURE;
    }
    spacing[0] = 4.1341146;
    spacing[1] = 4.1341146;
    spacing[2] = 8.0;
    logger->Info("Set dataset spacing : (" +std::to_string(spacing[0])+","+std::to_string(spacing[1])+","+
                 std::to_string(spacing[2])+")\n");
    using CSVWriterType =  itk::MyCSVNumericObjectFileWriter< OutputPixelType>;
    std::string baseName = "metrics";

    std::string metricsFilename;
    bool oldFile = false;
    if( !fs::exists(root + baseName +".csv")){
        metricsFilename = root + baseName +  ".csv";
        logger->Info("Writing metric file to: " + metricsFilename + "\n");
    }else{
        metricsFilename = root + baseName +"_"+ std::to_string(tic())+ ".csv";
        oldFile = true;
        logger->Warning(baseName +".csv already exists! writing new file to"+metricsFilename+"\n");
    }
    CSVWriterType::Pointer csvWriter = CSVWriterType::New();
    csvWriter->SetFileName(metricsFilename);
    csvWriter->ColumnHeadersPushBack("PSD label");//0

    csvWriter->ColumnHeadersPushBack("Bbox_col_start");//1
    csvWriter->ColumnHeadersPushBack("Bbox_col_end");//2
    csvWriter->ColumnHeadersPushBack("Bbox_row_start");//3
    csvWriter->ColumnHeadersPushBack("Bbox_row_end");//4
    csvWriter->ColumnHeadersPushBack("Bbox_slice_start");//5
    csvWriter->ColumnHeadersPushBack("Bbox_slice_end");//6

    csvWriter->ColumnHeadersPushBack("volume");//7

    csvWriter->ColumnHeadersPushBack("Centroid_col");//8
    csvWriter->ColumnHeadersPushBack("Centroid_row");//9
    csvWriter->ColumnHeadersPushBack("Centroid_slice");//10

    csvWriter->ColumnHeadersPushBack("Closest_col");//11
    csvWriter->ColumnHeadersPushBack("Closest_row");//12
    csvWriter->ColumnHeadersPushBack("Closest_slice");//13

    csvWriter->ColumnHeadersPushBack("psd->astro");//14

    csvWriter->ColumnHeadersPushBack("surface_col");//15
    csvWriter->ColumnHeadersPushBack("surface_row");//16
    csvWriter->ColumnHeadersPushBack("surface_slice");//17

    csvWriter->ColumnHeadersPushBack("mito->astro");//18

    //csvWriter->ColumnHeadersPushBack("mito_ID");//19

    LabelImageType::SizeValueType psdCount;
    if(oldFile){
        psdCount = csvFileRowCount(root + baseName +  ".csv");
        logger->Debug("Number of rows(psd) from old metric file = "+std::to_string(psdCount)+ "\n");
    }
    else{
        psdCount = psdRegionCount(root);
        logger->Debug("Number of psds from psd file = "+std::to_string(psdCount) + "\n");
    }

    for(size_t psd = 0; psd < psdCount; ++psd){
        std::string rowName = "psd " + std::to_string(psd+1);
        csvWriter->RowHeadersPushBack(rowName);
    }

    const unsigned int CSVFieldCount = csvWriter->GetColumnHeaderSize()-1;
    const unsigned int CSVRowCount = csvWriter->GetRowHeaderSize();
    DataMatrixType* csvData = new DataMatrixType;
    csvData->SetSize( CSVRowCount, CSVFieldCount );
    logger->Debug("Created ["+std::to_string(CSVRowCount)+"x"+std::to_string(CSVFieldCount)+"] data matrix\n");
    std::stringstream ss;
    //fill from previous data...
//    if (oldFile) {
//        logger->Warning("Old file exists.. filling data from old data\n");
//        using ReaderType = itk::CSVArray2DFileReader<OutputPixelType>;
//        ReaderType::Pointer reader = ReaderType::New();
//        reader->SetFileName(root + baseName + ".csv");
//        reader->SetFieldDelimiterCharacter(',');
//        reader->HasColumnHeadersOn();
//        reader->HasRowHeadersOn();
//        reader->Parse();
//        ReaderType::Array2DDataObjectPointer metrics = reader->GetOutput();
//        ss.str("\n");
//        for (size_t i = 0; i < CSVRowCount; ++i){
//            for (size_t j = 0; j < CSVFieldCount; ++j){
//                ss << (*metrics)(i,j) << "\t";
//                (*csvData)(i,j) = (*metrics)(i,j);
//            }
//            ss << "\n";
//        }
//        ss << "\n";
//        logger->Debug(ss.str());
//    }

    weightedMetricAnalysis(parser, logger, csvData);

    logger->Debug("Finished Computations\n");
    if(parser->ArgumentExists("-debug")) printMatrix(csvData);

    csvWriter->SetInput(csvData);
    csvWriter->Update();
    delete csvData;
    return EXIT_SUCCESS;
}

