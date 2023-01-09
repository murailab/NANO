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

#include <chrono>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkChangeInformationImageFilter.h>
#include <itkIdentityTransform.h>
#include <itkResampleImageFilter.h>
#include <itkCSVArray2DFileReader.h>

#include "itkMinNeighborInterpolateImageFunction.h"
#include "utils.h"


int makeIsotropic(const std::string& filename, float isoSpacing){
    size_t position = filename.find_last_of('.');
    std::string namePrefix = filename.substr(0,position);
    std::string nameSuffix = filename.substr(position, filename.size()-position);

    constexpr unsigned int Dimension = 3;
    using InputPixelType = float;
    using InternalPixelType = float;
    using OutputPixelType = float;

    using InputImageType = itk::Image<InputPixelType, Dimension>;
    using InternalImageType = itk::Image<InternalPixelType, Dimension>;

    using ReaderType = itk::ImageFileReader<InputImageType>;
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(filename);

    typename InputImageType::SpacingType  inputSpacing;
    inputSpacing[0] = 4.1341146f;
    inputSpacing[1] = 4.1341146f;
    inputSpacing[2] = 8;

    using InformationChangeFilterTypeFloat = itk::ChangeInformationImageFilter<InputImageType >;
    InformationChangeFilterTypeFloat::Pointer changeSpacing = InformationChangeFilterTypeFloat::New();
    changeSpacing->SetInput(reader->GetOutput());
    changeSpacing->SetOutputSpacing(inputSpacing);
    changeSpacing->ChangeSpacingOn();
    changeSpacing->Update();
    typename InputImageType::Pointer inputImage = changeSpacing->GetOutput();

    using OutputImageType = itk::Image<OutputPixelType, Dimension>;

    using ResampleFilterType = itk::ResampleImageFilter<InternalImageType, OutputImageType>;
    typename ResampleFilterType::Pointer resampler = ResampleFilterType::New();

    using TransformType = itk::IdentityTransform<double, Dimension>;
    typename TransformType::Pointer transform = TransformType::New();
    transform->SetIdentity();

    resampler->SetTransform(transform);

    using InterpolatorType = itk::MinNeighborInterpolateImageFunction<InternalImageType, double>;
    typename InterpolatorType::Pointer interpolator = InterpolatorType::New();
    InterpolatorType::SizeType radius;
    radius.Fill(1);
    interpolator->SetRadius(radius);
    resampler->SetInterpolator(interpolator);

    OutputImageType::SpacingType outputSpacing;
    outputSpacing[0] = isoSpacing;
    outputSpacing[1] = isoSpacing;
    outputSpacing[2] = isoSpacing;
    resampler->SetOutputSpacing(outputSpacing);
    resampler->SetOutputOrigin(inputImage->GetOrigin());
    resampler->SetOutputDirection(inputImage->GetDirection());


    typename InputImageType::SizeType inputSize = inputImage->GetLargestPossibleRegion().GetSize();

    using SizeValueType = InputImageType::SizeType::SizeValueType;

    const double dx = inputSize[0] * inputSpacing[0] / isoSpacing;
    const double dy = inputSize[1] * inputSpacing[1] / isoSpacing;
    const double dz = (inputSize[2] - 1) * inputSpacing[2] / isoSpacing;

    InputImageType::SizeType size;

    size[0] = static_cast<SizeValueType>(dx);
    size[1] = static_cast<SizeValueType>(dy);
    size[2] = static_cast<SizeValueType>(dz);

    resampler->SetSize(size);
    resampler->SetInput(inputImage);
    resampler->Update();

    using FileWriterType = itk::ImageFileWriter<OutputImageType>;
    typename FileWriterType::Pointer writer = FileWriterType::New();
    writer->SetFileName(namePrefix+ "-isotropic" + nameSuffix );
    writer->SetInput( resampler->GetOutput() );
    try
    {
        writer->Update();
        std::cout << "wrote file" << writer->GetFileName() << std::endl;
    }catch (itk::ExceptionObject &e){
        std::cout << e.what() << std::endl;
    }

    return EXIT_SUCCESS;
}

size_t csvFileRowCount(std::string metricsFileName){
    using ReaderType = itk::CSVArray2DFileReader<double>;
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( metricsFileName  );
    reader->SetFieldDelimiterCharacter( ',' );
    reader->HasColumnHeadersOn();
    reader->HasRowHeadersOn();
    reader->Parse();

    ReaderType::Array2DDataObjectPointer metrics = reader->GetOutput();
    auto surfacex = metrics->GetColumn(0);
    return surfacex.size();
}

size_t tic(){
    const auto p1 = std::chrono::system_clock::now();
    return std::chrono::duration_cast<std::chrono::seconds>(
            p1.time_since_epoch()).count();
}