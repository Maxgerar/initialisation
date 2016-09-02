//
//  NCC_function.hpp
//  Initialisation
//
//  Created by Maxime Gérard on 27/08/16.
//  Copyright © 2016 Maxime Gérard. All rights reserved.
//

#ifndef NCC_function_hpp
#define NCC_function_hpp

#include <stdio.h>
#include <iostream>
#include <time.h>

#include <dlib/matrix.h>
#include "itkImage.h"
#include "itkEuler3DTransform.h"
#include "itkResampleImageFilter.h"
#include "itkImageRegionConstIteratorWithIndex.h"

#include "itkLinearInterpolateImageFunction.h"

using namespace dlib;

typedef itk::Image<double,3> ImageType;

//transform
typedef itk::Euler3DTransform<double> EulerTransformType;
typedef itk::ResampleImageFilter<ImageType, ImageType> ResamplerType;
typedef itk::LinearInterpolateImageFunction<ImageType,double> InterpolatorType;

class NCC_function
{
public:
    //constructeur
    NCC_function();
    NCC_function(ImageType::Pointer moving, ImageType::Pointer fixed);
    
    //metric computation
    double operator()(const dlib::matrix<double>& params)const
    {
        double value =0;
        
        //std::cout<<"params from operator : "<<params<<std::endl;
        
        //defining images
        
        if( !m_fixedImage )
        {
            std::cout<< "Fixed image has not been assigned"<<std::endl;
        }
        
        
        
        //image qui bouge = image US
        
        
        if(!m_movingImage)
        {
            std::cout<<"Moving image has not been assigned"<<std::endl;
        }
        
        //transform the image
        ImageType::Pointer movedImage = TransformImage(params, 1);
        InterpolatorType::Pointer Interpolator = InterpolatorType::New();
        Interpolator->SetInputImage(movedImage);
        
        
        //init timer
        std::srand(time(NULL));
        std::time_t tbegin,tend;
        double texec = 0;
        tbegin = std::time(NULL);
        
        typedef  itk::ImageRegionConstIteratorWithIndex<ImageType> FixedIteratorType;
        
        FixedIteratorType ti( m_fixedImage, m_fixedImage->GetLargestPossibleRegion() );
        
        typename ImageType::IndexType index;
        
        //std::cout<<"fixed image dim : "<<m_fixedImage->GetLargestPossibleRegion().GetSize()<<std::endl;
        //std::cout<<"moved image dim : "<<movedImage->GetLargestPossibleRegion().GetSize()<<std::endl;
        
        //MeasureType measure;
        
        int NumberOfPixelsCounted = 0;
        
        //this->SetTransformParameters(parameters);
        
        typedef  typename itk::NumericTraits<double>::AccumulateType AccumulateType;
        
        AccumulateType sff = itk::NumericTraits<AccumulateType>::ZeroValue();
        AccumulateType smm = itk::NumericTraits<AccumulateType>::ZeroValue();
        AccumulateType sfm = itk::NumericTraits<AccumulateType>::ZeroValue();
        AccumulateType sf  = itk::NumericTraits<AccumulateType>::ZeroValue();
        AccumulateType sm  = itk::NumericTraits<AccumulateType>::ZeroValue();
        
        while( !ti.IsAtEnd() )
        {
            index = ti.GetIndex();
            
            ImageType::PointType inputPoint;
            m_fixedImage->TransformIndexToPhysicalPoint(index, inputPoint);
            
            //i'm not using masks and since the iterator is on the fixedimage, there is no reason the inputPoint wouldn't be inside the fixed image
           ImageType::IndexType ind;
//            if(!m_fixedImage->TransformPhysicalPointToIndex(inputPoint, ind))
//            {
//                ++ti;
//                continue;
//            }
            
            ImageType::PointType transformedPoint = this->m_transform->TransformPoint(inputPoint);
            //ImageType::IndexType ind;
            //if the transformed point isn't in the moved image, can't compute the metric
            if(!movedImage->TransformPhysicalPointToIndex(inputPoint,ind) ) //transformed
            {
                ++ti;
                continue;
            }
            
            
            
            if( Interpolator->IsInsideBuffer(inputPoint) ) //transformed
            {
                const double movingValue  = Interpolator->Evaluate(inputPoint);//transformpoint
                const double fixedValue   = ti.Get();
                sff += fixedValue  * fixedValue;
                smm += movingValue * movingValue;
                sfm += fixedValue  * movingValue;
                if(m_SubtractMean )
                {
                    //std::cout<<"subs"<<std::endl;
                    sf += fixedValue;
                    sm += movingValue;
                }
                NumberOfPixelsCounted++;
            }
            
            ++ti;
            ++ti;
        }
        
        //std::cout<<"number of pixels counted : "<<NumberOfPixelsCounted<<std::endl;
        
        if(m_SubtractMean && NumberOfPixelsCounted > 0 )
        {
            //std::cout<<"subs"<<std::endl;
            sff -= ( sf * sf / NumberOfPixelsCounted );
            smm -= ( sm * sm / NumberOfPixelsCounted );
            sfm -= ( sf * sm / NumberOfPixelsCounted );
        }
        
        const double denom = -1.0 * std::sqrt(sff * smm);
        //std::cout<<"denom : "<<denom<<std::endl;
        //std::cout<<"sfm : "<<sfm<<std::endl;
        //std::cout<<"smm : "<<smm<<std::endl;
        //std::cout<<"sff : "<<sff<<std::endl;
        
        if( NumberOfPixelsCounted > 0 && denom != 0.0 )
        {
            value = sfm / denom;
        }
        else
        {
            value = itk::NumericTraits<double>::ZeroValue();
        }
        
        tend = std::time(NULL);
        texec = std::difftime(tend,tbegin);
        
        std::cout<<"temps de parcours en s : "<<texec<<std::endl;
        std::cout<<"NCC value : "<<value<<std::endl;
        return value;
        
    }
    
    ImageType::Pointer TransformImage(const dlib::matrix<double>&params, int ind) const;
    
    //setters
    void setMaxRot(double rot){m_maxRot =rot;}
    void setMaxTrans(double trans){m_maxTrans = trans;}
    void setRadius(double rad){m_radius = rad;}
    
    void setMovingImage(ImageType::Pointer moving){m_movingImage = moving;}
    void setFixedImage(ImageType::Pointer fixed){m_fixedImage = fixed;}
    void setSubstractMeanBool(bool subMean){m_SubtractMean=subMean;}
    
    //getters
    double GetMaxRot(){return m_maxRot;}
    double GetMaxTrans(){return m_maxTrans;}
    double GetRadius(){return m_radius;}
    
private:
    ImageType::Pointer m_movingImage;
    ImageType::Pointer m_fixedImage;
    EulerTransformType::Pointer m_transform;
    double m_maxRot;
    double m_maxTrans;
    double m_radius;
    
    bool m_SubtractMean=false;
    
};

#endif /* NCC_function_hpp */
