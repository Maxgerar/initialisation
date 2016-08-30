//
//  NCC_function.cpp
//  Initialisation
//
//  Created by Maxime Gérard on 27/08/16.
//  Copyright © 2016 Maxime Gérard. All rights reserved.
//

#include "NCC_function.hpp"

NCC_function::NCC_function()
{
    m_maxRot =0;
    m_maxTrans=0;
    m_SubtractMean = false;
    
    m_transform = EulerTransformType::New();
}

NCC_function::NCC_function(ImageType::Pointer moving, ImageType::Pointer fixed)
{
    m_movingImage = moving;
    m_fixedImage = fixed;
    m_maxTrans = 0;
    m_maxRot = 0;
    
    m_SubtractMean = false;
    
    m_transform = EulerTransformType::New();
}

ImageType::Pointer NCC_function::TransformImage(const dlib::matrix<double>&params,int ind)const
{
    ImageType::Pointer imageTransformed;
    
    //1 represents euler tsf
    if(ind==1)
    {
        //m_transform = EulerTransformType::New();
        EulerTransformType::ParametersType parameters(6);
        //mise a l'echelle des parametres
        //std::cout<<"params : "<<params<<std::endl;
        parameters[0] = params(0)*(m_maxRot)/(m_radius);
        parameters[1] = params(1)*(m_maxRot)/(m_radius);
        parameters[2] = params(2)*(m_maxRot)/(m_radius);
        parameters[3] = params(3)*(m_maxTrans)/(m_radius);
        parameters[4] = params(4)*(m_maxTrans)/(m_radius);
        parameters[5] = params(5)*(m_maxTrans)/(m_radius);
        m_transform->SetParameters(parameters);
        std::cout<<"euler tsf parameters : "<<m_transform->GetParameters()<<std::endl;
        
        typename ImageType::SizeType sizeUS = m_movingImage->GetLargestPossibleRegion().GetSize();
        typename ImageType::PointType origin = m_movingImage->GetOrigin();
        typename ImageType::SpacingType spacing = m_movingImage->GetSpacing();
        typename ImageType::PointType center;
        center[0] = origin[0]+spacing[0]*sizeUS[0]/2;
        center[1] = origin[1]+spacing[1]*sizeUS[1]/2;
        center[2] = origin[2]+spacing[2]*sizeUS[2]/2;
        
        
        EulerTransformType::ParametersType eulerFixedParameters(3);
        eulerFixedParameters[0] =center[0];
        eulerFixedParameters[1] =center[1];
        eulerFixedParameters[2] =center[2];
        
        m_transform->SetFixedParameters(eulerFixedParameters);
        //std::cout<<"tsf fixed param : "<<transform->GetFixedParameters()<<std::endl;
        
        
        
        typename ResamplerType::Pointer resamplefilter = ResamplerType::New();
        resamplefilter->SetInput(m_movingImage);
        resamplefilter->SetTransform(m_transform);
        resamplefilter->SetSize(m_movingImage->GetLargestPossibleRegion().GetSize());
        resamplefilter->SetOutputOrigin(m_transform->GetInverseTransform()->TransformPoint(m_movingImage->GetOrigin()));
        resamplefilter->SetOutputSpacing(m_movingImage->GetSpacing());
        resamplefilter->SetOutputDirection(m_transform->GetInverseMatrix()*m_movingImage->GetDirection());
//        resamplefilter->SetSize(m_fixedImage->GetLargestPossibleRegion().GetSize());
//        resamplefilter->SetOutputSpacing(m_fixedImage->GetSpacing());
//        resamplefilter->SetOutputDirection(m_fixedImage->GetDirection());
//        resamplefilter->SetOutputOrigin(m_fixedImage->GetOrigin());
       
        //resamplefilter->SetTransform(transform);
        
        try {
            resamplefilter->Update();
        } catch (itk::ExceptionObject &e) {
            std::cerr<<"error while transforming moving image"<<std::endl;
            std::cerr<<e<<std::endl;
        }
        
        imageTransformed= resamplefilter->GetOutput();
    }
    
    return imageTransformed;

}