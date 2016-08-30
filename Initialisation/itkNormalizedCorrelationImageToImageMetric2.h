/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef itkNormalizedCorrelationImageToImageMetric2_h
#define itkNormalizedCorrelationImageToImageMetric2_h

#include "itkImageToImageMetric.h"
#include "itkPoint.h"
#include <dlib/matrix.h>

namespace itk
{
/** \class NormalizedCorrelationImageToImageMetric
 * \brief Computes similarity between two images to be registered
 *
 * This metric computes the correlation between pixels in the fixed image
 * and pixels in the moving image. The spatial correspondance between
 * fixed and moving image is established through a Transform. Pixel values are
 * taken from the fixed image, their positions are mapped to the moving
 * image and result in general in non-grid position on it. Values at these
 * non-grid position of the moving image are interpolated using a user-selected
 * Interpolator. The correlation is normalized by the autocorrelations of both
 * the fixed and moving images.
 *
 * A more negative metric value indicates a greater degree of correlation
 * between the fixed and moving image. This makes the metric simpler to use
 * with optimizers that strive to minimize their cost function by default.
 *
 * \ingroup RegistrationMetrics
 * \ingroup ITKRegistrationCommon
 */
template< typename TFixedImage, typename TMovingImage >
class NormalizedCorrelationImageToImageMetric2:
  public ImageToImageMetric< TFixedImage, TMovingImage >
{
public:

  /** Standard class typedefs. */
  typedef NormalizedCorrelationImageToImageMetric2         Self;
  typedef ImageToImageMetric< TFixedImage, TMovingImage > Superclass;
  typedef SmartPointer< Self >                            Pointer;
  typedef SmartPointer< const Self >                      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(NormalizedCorrelationImageToImageMetric, Object);

  /** Types transferred from the base class */
  typedef typename Superclass::RealType                RealType;
  typedef typename Superclass::TransformType           TransformType;
  typedef typename Superclass::TransformPointer        TransformPointer;
  typedef typename Superclass::TransformParametersType TransformParametersType;
  typedef typename Superclass::TransformJacobianType   TransformJacobianType;
  typedef typename Superclass::GradientPixelType       GradientPixelType;
  typedef typename Superclass::OutputPointType         OutputPointType;
  typedef typename Superclass::InputPointType          InputPointType;

  typedef typename Superclass::MeasureType             MeasureType;
  typedef typename Superclass::DerivativeType          DerivativeType;
  typedef typename Superclass::FixedImageType          FixedImageType;
  typedef typename Superclass::MovingImageType         MovingImageType;
  typedef typename Superclass::FixedImageConstPointer  FixedImageConstPointer;
  typedef typename Superclass::MovingImageConstPointer MovingImageConstPointer;

  /** Get the derivatives of the match measure. */
  void GetDerivative(const TransformParametersType & parameters,
                     DerivativeType & Derivative) const ITK_OVERRIDE;

  /**  Get the value for single valued optimizers. */
  MeasureType GetValue(const TransformParametersType & parameters) const ITK_OVERRIDE;

  /**  Get value and derivatives for multiple valued optimizers. */
  void GetValueAndDerivative(const TransformParametersType & parameters,
                             MeasureType & Value, DerivativeType & Derivative) const ITK_OVERRIDE;
    

  double operator()(const dlib::matrix<double> &params)const
    {
        //int nb_param = m_Transform->
        typename TransformType::ParametersType parameters;
        parameters[0]=params(0);
        parameters[1]=params(1);
        parameters[2]=params(2);
        parameters[3]=params(3);
        parameters[4]=params(4);
        parameters[5]=params(5);
        
        FixedImageConstPointer fixedImage = this->m_FixedImage;
        
        if( !fixedImage )
        {
            itkExceptionMacro(<< "Fixed image has not been assigned");
        }
        
        typedef  itk::ImageRegionConstIteratorWithIndex<FixedImageType> FixedIteratorType;
        
        FixedIteratorType ti( fixedImage, this->GetFixedImageRegion() );
        
        typename FixedImageType::IndexType index;
        
        MeasureType measure;
        
        this->m_NumberOfPixelsCounted = 0;
        
        this->SetTransformParameters(parameters);
        
        typedef  typename NumericTraits<MeasureType>::AccumulateType AccumulateType;
        
        AccumulateType sff = NumericTraits<AccumulateType>::ZeroValue();
        AccumulateType smm = NumericTraits<AccumulateType>::ZeroValue();
        AccumulateType sfm = NumericTraits<AccumulateType>::ZeroValue();
        AccumulateType sf  = NumericTraits<AccumulateType>::ZeroValue();
        AccumulateType sm  = NumericTraits<AccumulateType>::ZeroValue();
        
        while( !ti.IsAtEnd() )
        {
            index = ti.GetIndex();
            
            InputPointType inputPoint;
            fixedImage->TransformIndexToPhysicalPoint(index, inputPoint);
            
            if( this->m_FixedImageMask && !this->m_FixedImageMask->IsInside(inputPoint) )
            {
                ++ti;
                continue;
            }
            
            OutputPointType transformedPoint = this->m_Transform->TransformPoint(inputPoint);
            
            if( this->m_MovingImageMask && !this->m_MovingImageMask->IsInside(transformedPoint) )
            {
                ++ti;
                continue;
            }
            
            if( this->m_Interpolator->IsInsideBuffer(transformedPoint) )
            {
                const RealType movingValue  = this->m_Interpolator->Evaluate(transformedPoint);
                const RealType fixedValue   = ti.Get();
                sff += fixedValue  * fixedValue;
                smm += movingValue * movingValue;
                sfm += fixedValue  * movingValue;
                if( this->m_SubtractMean )
                {
                    sf += fixedValue;
                    sm += movingValue;
                }
                this->m_NumberOfPixelsCounted++;
            }
            
            ++ti;
        }
        
        if( this->m_SubtractMean && this->m_NumberOfPixelsCounted > 0 )
        {
            sff -= ( sf * sf / this->m_NumberOfPixelsCounted );
            smm -= ( sm * sm / this->m_NumberOfPixelsCounted );
            sfm -= ( sf * sm / this->m_NumberOfPixelsCounted );
        }
        
        const RealType denom = -1.0 * std::sqrt(sff * smm);
        
        if( this->m_NumberOfPixelsCounted > 0 && denom != 0.0 )
        {
            measure = sfm / denom;
        }
        else
        {
            measure = NumericTraits<MeasureType>::ZeroValue();
        }
        
        return measure;

    }
    

  /** Set/Get SubtractMean boolean. If true, the sample mean is subtracted
   * from the sample values in the cross-correlation formula and
   * typically results in narrower valleys in the cost function.
   * Default value is false. */
  itkSetMacro(SubtractMean, bool);
  itkGetConstReferenceMacro(SubtractMean, bool);
  itkBooleanMacro(SubtractMean);

protected:
  NormalizedCorrelationImageToImageMetric2();
  virtual ~NormalizedCorrelationImageToImageMetric2() {}
  void PrintSelf(std::ostream & os, Indent indent) const ITK_OVERRIDE;

private:
  NormalizedCorrelationImageToImageMetric2(const Self &); //purposely not
                                                         // implemented
  void operator=(const Self &);                          //purposely not
                                                         // implemented

  bool m_SubtractMean;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkNormalizedCorrelationImageToImageMetric2.hxx"
#endif

#endif
