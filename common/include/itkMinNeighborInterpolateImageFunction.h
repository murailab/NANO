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

#ifndef NANOARCHITECTURE_ITKMINNEIGHBORINTERPOLATEIMAGEFUNCTION_H
#define NANOARCHITECTURE_ITKMINNEIGHBORINTERPOLATEIMAGEFUNCTION_H
#include "itkInterpolateImageFunction.h"
#include "itkConstNeighborhoodIterator.h"

namespace itk
{
    template <typename TInputImage, typename TCoordRep = double>
    class ITK_TEMPLATE_EXPORT MinNeighborInterpolateImageFunction
            : public InterpolateImageFunction<TInputImage, TCoordRep>
    {
    public:
        ITK_DISALLOW_COPY_AND_ASSIGN(MinNeighborInterpolateImageFunction);

        using Self = MinNeighborInterpolateImageFunction;
        using Superclass = InterpolateImageFunction<TInputImage, TCoordRep>;
        using Pointer = SmartPointer<Self>;
        using ConstPointer = SmartPointer<const Self>;

        itkTypeMacro(MinNeighborInterpolateImageFunction, InterpolateImageFunction);

        itkNewMacro(Self);

        using OutputType = typename Superclass::OutputType;

        using InputImageType = typename Superclass::InputImageType;

        static constexpr unsigned int ImageDimension = Superclass::ImageDimension;

        using IndexType = typename Superclass::IndexType;

        using SizeType = typename Superclass::SizeType;

        using ContinuousIndexType = typename Superclass::ContinuousIndexType;

        itkSetMacro(Radius, SizeType);

        SizeType
        GetRadius() const override
        {
            return m_Radius;
        }

        OutputType
        EvaluateAtContinuousIndex(const ContinuousIndexType & index) const override
        {
            IndexType nindex;
            OutputType outputValue = std::numeric_limits<OutputType>::max();
            this->ConvertContinuousIndexToNearestIndex(index, nindex);

            ConstNeighborhoodIterator<TInputImage> iterator(m_Radius, this->GetInputImage(), this->GetInputImage()->GetRequestedRegion());
            iterator.SetLocation(nindex);
            for(size_t i = 0; i < iterator.Size(); ++i){
                auto currentValue = static_cast<OutputType>(iterator.GetPixel(i));
                if (outputValue > currentValue){
                    outputValue = currentValue;
                }
            }
            return outputValue;
        }

    protected:
        MinNeighborInterpolateImageFunction(){
            this->m_Radius.Fill(1);
        }
        ~MinNeighborInterpolateImageFunction() override = default;
        void
        PrintSelf(std::ostream & os, Indent indent) const override {
            Superclass::PrintSelf(os, indent);
        }
    private:
        SizeType m_Radius;
    };
} // end namespace itk
#endif //NANOARCHITECTURE_ITKMINNEIGHBORINTERPOLATEIMAGEFUNCTION_H
