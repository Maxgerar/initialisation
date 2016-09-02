//
//  main.cpp
//  Initialisation
//
//  Created by Maxime Gérard on 19/04/16.
//  Copyright © 2016 Maxime Gérard. All rights reserved.
//

#include <iostream>
#include <string>
#include <cmath>
#include <fstream>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkNiftiImageIO.h"
#include "itkMetaImageIO.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkGradientMagnitudeImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkAdaptiveHistogramEqualizationImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkEuler3DTransform.h"
#include "itkNearestNeighborInterpolateImageFunction.h"

#include "itkHessian3DToVesselnessMeasureImageFilter.h"
#include "itkMultiScaleHessianBasedMeasureImageFilter.h"
#include "itkHessianToObjectnessMeasureImageFilter.h"
#include "itkSymmetricSecondRankTensor.h"
#include "itkRescaleIntensityImageFilter.h"

//image processing
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryMorphologicalClosingImageFilter.h"
#include "itkShrinkImageFilter.h"

#include "itkRescaleIntensityImageFilter.h"
//classic registration

#include "itkImageRegistrationMethod.h"
//#include "itkNormalizedCorrelationImageToImageMetric.h"
#include "itkNormalizedCorrelationImageToImageMetric2.h"
#include "itkKappaStatisticImageToImageMetric.h"
#include "itkMeanReciprocalSquareDifferenceImageToImageMetric.h"
#include "itkAmoebaOptimizer.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkCommand.h"
#include "itkMeanSquaresImageToImageMetric.h"

#include "itkLBFGSBOptimizer2.h"
#include "itkRegularStepGradientDescentOptimizer.h"

#include "NCC_function.hpp"

//mutliresolution registration
#include "itkMultiResolutionImageRegistrationMethod.h"
#include "itkMultiResolutionPyramidImageFilter.h"

#include "itkCommand.h"



typedef itk::Image<double,3> ImageType;
typedef itk::Image<unsigned char,3> BinaryImageType;
typedef itk::SymmetricSecondRankTensor< double, 3 > HessianPixelType;
typedef itk::Image< HessianPixelType, 3 > HessianImageType;

typedef itk::ImageFileReader<ImageType> ReaderType;
typedef itk::ImageFileReader<BinaryImageType> BinaryReaderType;
typedef itk::ImageFileWriter<ImageType> WriterType;
typedef itk::ImageFileWriter<BinaryImageType> BinaryWriterType;

typedef itk::MinimumMaximumImageCalculator<ImageType> MinMaxCalculatorType;
typedef itk::GradientMagnitudeImageFilter<ImageType, ImageType> GradientFilterType;

typedef itk::BinaryBallStructuringElement<BinaryImageType::PixelType,3> StructuringElementType;
typedef itk::BinaryDilateImageFilter<BinaryImageType, BinaryImageType, StructuringElementType> DilateFilterType;
typedef itk::BinaryMorphologicalClosingImageFilter<BinaryImageType, BinaryImageType, StructuringElementType> ClosingFilterType;

typedef itk::Euler3DTransform<double> EulerTransformType;
typedef itk::ResampleImageFilter<ImageType, ImageType> ResampleFilterType;
//typedef itk::ShrinkImageFilter<BinaryImageType, BinaryImageType> ShrinkFilterType;
typedef itk::ShrinkImageFilter<ImageType, ImageType> ShrinkFilterType;
typedef itk::RescaleIntensityImageFilter<ImageType,ImageType> IntensityRescaleFilterType;

typedef itk::BinaryThresholdImageFilter<ImageType, ImageType> BinaryThresholdFilterType;
typedef itk::AdaptiveHistogramEqualizationImageFilter<ImageType> HistogramEqualizerType;

typedef itk::HessianToObjectnessMeasureImageFilter<HessianImageType, ImageType> ObjectnessFilterType;
typedef itk::Hessian3DToVesselnessMeasureImageFilter<double> VesselnessMeasureFilterType;
typedef itk::MultiScaleHessianBasedMeasureImageFilter<ImageType, HessianImageType,ImageType> MultiScaleEnhancementFilterType;

#include<dlib/optimization.h>
#include <dlib/matrix.h>

using namespace std;
using namespace dlib;

////classes recalage multi res
//
//template <typename TRegistration>
//class RegistrationInterfaceCommand : public itk::Command
//{
//    // Software Guide : EndCodeSnippet
//    // Software Guide : BeginLatex
//    //
//    // We then define \code{Self}, \code{Superclass}, \code{Pointer},
//    // \code{New()} and a constructor in a similar fashion to the
//    // \code{CommandIterationUpdate} class in Section
//    // \ref{sec:MonitoringImageRegistration}.
//    //
//    // Software Guide : EndLatex
//    // Software Guide : BeginCodeSnippet
//public:
//    typedef  RegistrationInterfaceCommand   Self;
//    typedef  itk::Command                   Superclass;
//    typedef  itk::SmartPointer<Self>        Pointer;
//    itkNewMacro( Self );
//protected:
//    RegistrationInterfaceCommand() {};
//    // Software Guide : EndCodeSnippet
//    // Software Guide : BeginLatex
//    //
//    // For convenience, we declare types useful for converting pointers
//    // in the \code{Execute()} method.
//    //
//    // Software Guide : EndLatex
//    // Software Guide : BeginCodeSnippet
//public:
//    typedef   TRegistration                              RegistrationType;
//    typedef   RegistrationType *                         RegistrationPointer;
//    //typedef   itk::AmoebaOptimizer   OptimizerType;
//    typedef itk::LBFGSBOptimizer OptimizerType;
//    typedef   OptimizerType *                            OptimizerPointer;
//    // Software Guide : EndCodeSnippet
//    // Software Guide : BeginLatex
//    //
//    // Two arguments are passed to the \code{Execute()} method: the first
//    // is the pointer to the object which invoked the event and the
//    // second is the event that was invoked.
//    //
//    // Software Guide : EndLatex
//    // Software Guide : BeginCodeSnippet
//    void Execute(itk::Object * object, const itk::EventObject & event)
//    {
//        // Software Guide : EndCodeSnippet
//        // Software Guide : BeginLatex
//        //
//        // First we verify if that the event invoked is of the right type.
//        // If not, we return without any further action.
//        //
//        // Software Guide : EndLatex
//        // Software Guide : BeginCodeSnippet
//        if( !(itk::IterationEvent().CheckEvent( &event )) )
//        {
//            return;
//        }
//        // Software Guide : EndCodeSnippet
//        // Software Guide : BeginLatex
//        //
//        // We then convert the input object pointer to a RegistrationPointer.
//        // Note that no error checking is done here to verify if the
//        // \code{dynamic\_cast} was successful since we know the actual object
//        // is a multi-resolution registration method.
//        //
//        // Software Guide : EndLatex
//        // Software Guide : BeginCodeSnippet
//        RegistrationPointer registration = static_cast<RegistrationPointer>( object );
//        // Software Guide : EndCodeSnippet
//        // Software Guide : BeginLatex
//        //
//        // If this is the first resolution level we set the maximum step length
//        // (representing the first step size) and the minimum step length (representing
//        // the convergence criterion) to large values.  At each subsequent resolution
//        // level, we will reduce the minimum step length by a factor of 10 in order to
//        // allow the optimizer to focus on progressively smaller regions. The maximum
//        // step length is set up to the current step length. In this way, when the
//        // optimizer is reinitialized at the beginning of the registration process for
//        // the next level, the step length will simply start with the last value used
//        // for the previous level. This will guarantee the continuity of the path
//        // taken by the optimizer through the parameter space.
//        //
//        // Software Guide : EndLatex
//        // Software Guide : BeginCodeSnippet
//        OptimizerPointer optimizer = static_cast< OptimizerPointer >(registration->GetModifiableOptimizer() );
//        std::cout << "-------------------------------------" << std::endl;
//        std::cout << "MultiResolution Level : "
//        << registration->GetCurrentLevel()  << std::endl;
//        std::cout << std::endl;
//        if ( registration->GetCurrentLevel() == 0 )
//        {
//            //a adapter a amoeba
//            //optimizer->SetFunctionConvergenceTolerance(0.1);
//            //optimizer->SetParametersConvergenceTolerance(1.0);
//            //optimizer->SetMaximumStepLength( 16.00 );
//            //optimizer->SetMinimumStepLength( 0.01 );
//        }
//        else
//        {
//            //optimizer->SetCostFunctionConvergenceFactor(optimizer->GetCostFunctionConvergenceFactor()*(1e-2));
//            //optimizer->SetFunctionConvergenceTolerance(0.1);
//            //optimizer->SetParametersConvergenceTolerance(1.0);
//            //optimizer->SetMaximumStepLength(optimizer->GetMaximumStepLength() * 0.25 );
//            //optimizer->SetMinimumStepLength(optimizer->GetMinimumStepLength() * 0.1 );
//        }
//    }
//    
//    // Software Guide : EndCodeSnippet
//    // Software Guide : BeginLatex
//    //
//    // Another version of the \code{Execute()} method accepting a \code{const}
//    // input object is also required since this method is defined as pure virtual
//    // in the base class.  This version simply returns without taking any action.
//    //
//    // Software Guide : EndLatex
//    // Software Guide : BeginCodeSnippet
//    void Execute(const itk::Object * , const itk::EventObject & )
//    { return; }
//};

//classe pour suivre le recalage
class CommandIterationUpdate : public itk::Command
{
    public :
    typedef CommandIterationUpdate Self;
    typedef itk::Command SuperClass;
    typedef itk::SmartPointer<Self> Pointer;
    itkNewMacro(Self);
    
protected:
    CommandIterationUpdate()
    {
        m_IterationNumber =0;
    }
public:
    //typedef itk::AmoebaOptimizer OptimizerType;
    typedef itk::LBFGSBOptimizer2 OptimizerType;
    //typedef itk::LBFGSBOptimizer OptimizerType;
    //typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
    typedef const OptimizerType * OptimizerPointer;
    
    void Execute(itk::Object *caller, const itk::EventObject &event) ITK_OVERRIDE
    {
        Execute( (const itk::Object *)caller, event);
    }
    
    void Execute(const itk::Object * object,
                 const itk::EventObject & event) ITK_OVERRIDE
    {
        OptimizerPointer optimizer = static_cast< OptimizerPointer >( object );
        if( ! itk::IterationEvent().CheckEvent( &event ) )
        {
            return;
        }
        std::cout << m_IterationNumber++ << "   "<<endl;
        std::cout << optimizer->GetCachedValue() << "   "<<endl;
        std::cout << optimizer->GetCachedCurrentPosition() << std::endl;
        
  
    }
private:
    unsigned long m_IterationNumber;
};


int main(int argc, const char * argv[]) {
    

    
    string filenameUS;
    string filenameIRM;
    string filenameVeins;
    string filenameCTA;
    string outputPath;
    string filenameVesselness;
    
    bool vesselness = true;
    bool initialized = false;

    /*********************
     * RUNTIME ARG
     ******************/
    
    cout<<"runtime arguments acquisition"<<endl;
    
    for(int i = 0; i < argc; ++i)
    {
        //input image US
        if(strcmp(argv[i], "-iUS")==0)
        {
            i++;
            filenameUS = argv[i];
        }
        
        //input image IRM
        if(strcmp(argv[i], "-iIRM")==0)
        {
            i++;
            filenameIRM = argv[i];
        }
        
        if(strcmp(argv[i], "-iVesselness")==0)
        {
            i++;
            filenameVesselness = argv[i];
            vesselness = false;
            cout<<"vesselness provided"<<endl;
        }
        
        if(strcmp(argv[i], "-initialized")==0)
        {
            i++;
            initialized = true;
            cout<<"intialized"<<endl;
        }
        

        if(strcmp(argv[i], "-iMaskVeins")==0)
        {
            i++;
            filenameVeins= argv[i];
            cout<<"Use of mask image of liver"<<endl;
        }
        
        if(strcmp(argv[i], "-iCTA")==0)
        {
            i++;
            filenameCTA= argv[i];
            cout<<"Use of CTA US data"<<endl;
        }
        
        if(strcmp(argv[i], "-o")==0)
        {
            i++;
            outputPath = argv[i];
        }
        
        
    }
    
    /**********************
     * VERIFICATION INPUTS
     ********************/
    
    cout<<"input error handling"<<endl;
    
    if(filenameUS == "")
    {
        cerr<<"Input US file not provided"<<endl;
        return EXIT_FAILURE;
        
    }
    
    if(filenameIRM=="")
    {
        cerr<<"input IRM file not provided"<<endl;
        return EXIT_FAILURE;
    }
    
//    if(filenameVeins=="")
//    {
//        cerr<<"input veins mask file not provided"<<endl;
//        return EXIT_FAILURE;
//    }
    
    if(filenameCTA=="")
    {
        cerr<<"input CTA file not provided"<<endl;
        return EXIT_FAILURE;
    }
    
    
    
    if(outputPath == "")
    {
        cerr<<"output path not provided"<<endl;
        return EXIT_FAILURE;
    }
    
    /**********************
     * US READING
     *********************/
    
    cout<<"Reading images"<<endl;
    
    ImageType::Pointer image_US = ImageType::New();
  
    //image US de base
    
    ReaderType::Pointer reader1 = ReaderType::New();
    itk::MetaImageIO::Pointer m_io = itk::MetaImageIO::New();
    reader1->SetImageIO(m_io);
    reader1->SetFileName(filenameUS);
    try {
        reader1->Update();
    } catch (itk::ExceptionObject &e) {
        cerr<<"Error while reading US image"<<endl;
        cerr<<e<<endl;
        EXIT_FAILURE;
    }
    
    image_US = reader1->GetOutput();
    cout<<"test lecture US"<<endl;
    cout<<"dimensions US : "<<image_US->GetLargestPossibleRegion().GetSize()<<endl;
    
    //CTA signal image US
    
    ImageType::Pointer CTA = ImageType::New();
    
    ReaderType::Pointer readerCTA = ReaderType::New();
    readerCTA->SetImageIO(m_io);
    readerCTA->SetFileName(filenameCTA);
    try {
        readerCTA->Update();
    } catch (itk::ExceptionObject &e) {
        cerr<<"Error while reading US image"<<endl;
        cerr<<e<<endl;
        EXIT_FAILURE;
    }
    
    CTA = readerCTA->GetOutput();
    cout<<"test lecture CTA"<<endl;
    cout<<"dimensions CTA : "<<CTA->GetLargestPossibleRegion().GetSize()<<endl;
    
    
    /**********************
     * MRI READING
     ***********************/
    
    ImageType::Pointer image_IRM = ImageType::New();
    
    ReaderType::Pointer reader2 = ReaderType::New();
    itk::NiftiImageIO::Pointer n_io = itk::NiftiImageIO::New();
    reader2->SetImageIO(n_io);
    reader2->SetFileName(filenameIRM);
    try {
        reader2->Update();
    } catch (itk::ExceptionObject &e) {
        cerr<<"Error while reading MRI image"<<endl;
        cerr<<e<<endl;
        EXIT_FAILURE;
    }
    
    image_IRM= reader2->GetOutput();
    cout<<"test lecture IRM"<<endl;
    cout<<"dimensions IRM : "<<image_IRM->GetLargestPossibleRegion().GetSize()<<endl;
    

    

        cout<<"reading mask sus-hep"<<endl;
//        BinaryImageType::Pointer mask_veins = BinaryImageType::New();
//        
//        BinaryReaderType::Pointer readerBin = BinaryReaderType::New();
//        readerBin->SetImageIO(n_io);
//        readerBin->SetFileName(filenameVeins);
//        try {
//            readerBin->Update();
//        } catch (itk::ExceptionObject &e) {
//            cerr<<"Error while reading MRI image"<<endl;
//            cerr<<e<<endl;
//            EXIT_FAILURE;
//        }
//        
//        mask_veins = readerBin->GetOutput();
//        cout<<"test lecture mask veins "<<endl;
//        cout<<"dimensions mask veins : "<<mask_veins->GetLargestPossibleRegion().GetSize()<<endl;
        
        ImageType::Pointer mask_veins = ImageType::New();
        ReaderType::Pointer reader_P = ReaderType::New();
        reader_P->SetImageIO(n_io);
        reader_P->SetFileName(filenameVeins);
        try {
            reader_P->Update();
        } catch (itk::ExceptionObject &e) {
            cerr<<"Error while reading MRI image"<<endl;
            cerr<<e<<endl;
            EXIT_FAILURE;
        }
        
        mask_veins = reader_P->GetOutput();
                cout<<"test lecture mask veins "<<endl;
                cout<<"dimensions mask veins : "<<mask_veins->GetLargestPossibleRegion().GetSize()<<endl;
        

    
    /**************************
     * TEST VESSELNESS PORTAL
     **************************/
    
    //vesselness sur IRM
    cout<<"vesselness image IRM"<<endl;
    ImageType::Pointer vesselImage = ImageType::New();
    
    if(vesselness)
    {
        typedef itk::SymmetricSecondRankTensor<double,3> HessianPixelType;
        typedef itk::Image<HessianPixelType,3> HessianImageType;
        typedef itk::HessianToObjectnessMeasureImageFilter<HessianImageType, ImageType> ObjectnessFilterType;
        
        ObjectnessFilterType::Pointer vesselnessFilter = ObjectnessFilterType::New();
        vesselnessFilter->SetBrightObject(true);
        vesselnessFilter->SetScaleObjectnessMeasure(false);
        vesselnessFilter->SetObjectDimension(1);
        vesselnessFilter->SetAlpha(0.5);
        vesselnessFilter->SetBeta(0.5);
        vesselnessFilter->SetGamma(10.0);
        //vesselnessFilter->GetInput()->GetN
        
        typedef itk::MultiScaleHessianBasedMeasureImageFilter<ImageType, HessianImageType,ImageType> MultiScaleEnhancementFilterType;
        
        MultiScaleEnhancementFilterType::Pointer enhancer = MultiScaleEnhancementFilterType::New();
        enhancer->SetInput(image_IRM);
        enhancer->SetHessianToMeasureFilter(vesselnessFilter);
        enhancer->SetSigmaStepMethodToLogarithmic();
        enhancer->SetSigmaMaximum(2.5); //3.0
        enhancer->SetSigmaMinimum(2.0); //2.5
        enhancer->SetNumberOfSigmaSteps(5);
        
        vesselImage = enhancer->GetOutput();
        
    }
    
    else
    {
        cout<<"reading the vesselness image"<<endl;
        ReaderType::Pointer reader_V = ReaderType::New();
        reader_V->SetImageIO(n_io);
        reader_V->SetFileName(filenameVesselness);
        try {
            reader_V->Update();
        } catch (itk::ExceptionObject &e) {
            cerr<<"Error while reading MRI image"<<endl;
            cerr<<e<<endl;
            EXIT_FAILURE;
        }
        
        vesselImage = reader_V->GetOutput();

        
    }

    
    
    /********************************
     * Rescale intensity VESSELNESS
     ************************************/
    
    IntensityRescaleFilterType::Pointer rescalerI = IntensityRescaleFilterType::New();
    rescalerI->SetInput(vesselImage);
    rescalerI->SetOutputMaximum(255);
    rescalerI->SetOutputMinimum(0);
    rescalerI->Update();
    
    ImageType::Pointer vesselImageR = rescalerI->GetOutput();
    
    WriterType::Pointer writer4 =  WriterType::New();
    writer4->SetImageIO(n_io);
    string out4 = outputPath+"/vesselnessIRMImage.nii.gz";
    writer4->SetFileName(out4);
    writer4->SetInput(vesselImageR);
    try {
        writer4->Update();
    } catch (itk::ExceptionObject &e) {
        cerr<<"error whilte writing registered image"<<endl;
        cerr<<e<<endl;
        return EXIT_FAILURE;
    }
    

    
    
    /********************
     * DOWNSAMPLE IRM
     *******************/
    ShrinkFilterType::Pointer shrinkFilter2 = ShrinkFilterType::New();
    //shrinkFilter2->SetInput(mask_veins);
    shrinkFilter2->SetInput(vesselImageR);
    
    shrinkFilter2->SetShrinkFactor(0, 3);
    shrinkFilter2->SetShrinkFactor(1, 3);
    shrinkFilter2->SetShrinkFactor(2, 3);
    try {
        shrinkFilter2->Update();
    } catch (itk::ExceptionObject &e) {
        cerr<<"error while downsampling US image"<<endl;
        cerr<<e<<endl;
        return EXIT_FAILURE;
    }
    
    //BinaryImageType::Pointer shrunk_mask_veins = shrinkFilter2->GetOutput();
    ImageType::Pointer shrunk_mask_IRM = shrinkFilter2->GetOutput();
    
    WriterType::Pointer writerIRM_shrunk = WriterType::New();
    writerIRM_shrunk->SetImageIO(n_io);
    writerIRM_shrunk->SetInput(shrunk_mask_IRM);
    string out_IRM_shrunk = outputPath+"/shrunk_MRI.nii.gz";
    writerIRM_shrunk->SetFileName(out_IRM_shrunk);
    try {
        writerIRM_shrunk->Update();
    } catch (itk::ExceptionObject &e) {
        cerr<<"error while writing shrunk MRI image"<<endl;
        cerr<<e<<endl;
    }
    
    
    //MEASURE COMPUTATION TIME CONSIDERING ONLY THE US PROCESSING + registration
    
    //to determine time of computation
    std::srand(time(NULL));
    std::time_t tbegin,tend;
    double texec = 0;
    tbegin = std::time(NULL);

    
    
    /**********************************************************
     * ORIENTATION : GROSS INITIAL ALIGNMENT A PRIORI KNOWLEDGE
     ***********************************************************/
    //aligning image origins
    
    ImageType::Pointer US_translated = ImageType::New();
    ImageType::Pointer CTA_translated = ImageType::New();
    
    ImageType::SizeType sizeIRM = image_IRM->GetLargestPossibleRegion().GetSize();
    ImageType::PointType origineIRM = image_IRM->GetOrigin();
    ImageType::SpacingType spacingIRM = image_IRM->GetSpacing();
    
    ImageType::SizeType sizeUS = image_US->GetLargestPossibleRegion().GetSize();
    ImageType::PointType origineUS = image_US->GetOrigin();
    ImageType::SpacingType spacingUS = image_US->GetSpacing();
    
    if(!initialized)
    {
        cout<<"aligning image origins"<<endl;
        
        //euler tsf
        
        //determining image centers
        //MRI
//        ImageType::SizeType sizeIRM = image_IRM->GetLargestPossibleRegion().GetSize();
//        ImageType::PointType origineIRM = image_IRM->GetOrigin();
//        ImageType::SpacingType spacingIRM = image_IRM->GetSpacing();
        ImageType::PointType centerIRM;
        
        centerIRM[0] = origineIRM[0]+spacingIRM[0]*sizeIRM[0]/2;
        centerIRM[1] = origineIRM[1]+spacingIRM[1]*sizeIRM[2]/2;
        centerIRM[2] = origineIRM[2]+spacingIRM[2]*sizeIRM[1]/2;
        
        cout<<"center MRI : "<<centerIRM<<endl;
        
//        ImageType::SizeType sizeUS = image_US->GetLargestPossibleRegion().GetSize();
//        ImageType::PointType origineUS = image_US->GetOrigin();
//        ImageType::SpacingType spacingUS = image_US->GetSpacing();
        ImageType::PointType centerUS;
        
        centerUS[0] = origineUS[0]+spacingUS[0]*sizeUS[0]/2;
        centerUS[1] = origineUS[1]+spacingUS[1]*sizeUS[1]/2;
        centerUS[2] = origineUS[2]+spacingUS[2]*sizeUS[2]/2;
        
        cout<<"center US : "<<centerUS<<endl;
        
        //ImageType::PointType tsl = origineIRM-origineUS;
        
        EulerTransformType::Pointer initial_tsf = EulerTransformType::New();
        EulerTransformType::ParametersType Parameters(6);
        Parameters[0] = 0.0;//+0.2;
        Parameters[1] = 0.0;
        Parameters[2] = 0.0;//-0.1
        //on sait que le foie est situe a droite et dans la partie superieure dans l'abdomen
        Parameters[3] = -(centerIRM[0]-centerUS[0])+50;//+30;//+20;//10
        Parameters[4] = -(centerIRM[1]-centerUS[1]);//+10;//+10; //+5,15,0 +20 ?
        Parameters[5] = -(centerIRM[2]-centerUS[2])-45;//-45;//-25 -50-90
        
        cout<<"parameters : "<<Parameters<<endl;
        
        initial_tsf->SetParameters(Parameters);
        
        EulerTransformType::ParametersType fixedParameters(3);
        fixedParameters[0] = centerUS[0];
        fixedParameters[1] = centerUS[1];
        fixedParameters[2] = centerUS[2];
        
        initial_tsf->SetFixedParameters(fixedParameters);
        
        
        ResampleFilterType::Pointer resampler = ResampleFilterType::New();
        resampler->SetInput(image_US);
        resampler->SetTransform(initial_tsf);
        resampler->SetSize(image_US->GetLargestPossibleRegion().GetSize());
        resampler->SetOutputSpacing(image_US->GetSpacing());
        resampler->SetOutputDirection(initial_tsf->GetInverseMatrix()*image_US->GetDirection());
        resampler->SetOutputOrigin(initial_tsf->GetInverseTransform()->TransformPoint(image_US->GetOrigin()));
        
        US_translated = resampler->GetOutput();
        
        
        WriterType::Pointer writer3 = WriterType::New();
        string out8 =outputPath+"/translated_US_test.nii.gz";
        writer3->SetImageIO(n_io);
        writer3->SetInput(US_translated);
        writer3->SetFileName(out8);
        try {
            writer3->Update();
        } catch (itk::ExceptionObject &e) {
            cerr<<"error whilte writing registered image"<<endl;
            cerr<<e<<endl;
            return EXIT_FAILURE;
        }
        
        //YOU GOT TO APPLY THE SAME TSF TO THE CTA DATA !!!
        
        ResampleFilterType::Pointer resampler2 = ResampleFilterType::New();
        resampler2->SetInput(CTA);
        resampler2->SetTransform(initial_tsf);
        resampler2->SetSize(CTA->GetLargestPossibleRegion().GetSize());
        resampler2->SetOutputSpacing(CTA->GetSpacing());
        resampler2->SetOutputDirection(initial_tsf->GetInverseMatrix()*CTA->GetDirection());
        resampler2->SetOutputOrigin(initial_tsf->GetInverseTransform()->TransformPoint(CTA->GetOrigin()));
        
        CTA_translated = resampler2->GetOutput();
        
        WriterType::Pointer writerCTA = WriterType::New();
        string outCTA =outputPath+"/translated_CTA_test.nii.gz";
        writerCTA->SetImageIO(n_io);
        writerCTA->SetInput(CTA_translated);
        writerCTA->SetFileName(outCTA);
        try {
            writerCTA->Update();
        } catch (itk::ExceptionObject &e) {
            cerr<<"error whilte writing registered image"<<endl;
            cerr<<e<<endl;
            return EXIT_FAILURE;
        }

        
    }
    
    else
    {
        US_translated = image_US;
        CTA_translated = CTA;
    }
    
    
    /*************************
     * THRESHOLDING
     ***********************/
    
        cout<<"binary Thresholding"<<endl;
    
        BinaryThresholdFilterType::Pointer thresholder = BinaryThresholdFilterType::New();
        thresholder->SetInput(US_translated);
        thresholder->SetLowerThreshold(1.0);
        thresholder->SetUpperThreshold(25.0);
        thresholder->SetInsideValue(255);
        thresholder->SetOutsideValue(0);
    
        ImageType::Pointer mask_US = thresholder->GetOutput();
    
        WriterType::Pointer binWriter = WriterType::New();
        binWriter->SetImageIO(n_io);
        string out_bin = outputPath+"/Thresholded_US.nii.gz";
        binWriter->SetInput(mask_US);
        binWriter->SetFileName(out_bin);
        try {
            binWriter->Update();
        } catch (itk::ExceptionObject &e) {
            cerr<<"error whilte writing registered image"<<endl;
            cerr<<e<<endl;
            return EXIT_FAILURE;
        }
    

    
  

    
    /******************************
     * Diminuation resolution US
     *****************************/
    
    cout<<"downsampling US data"<<endl;
    
    ShrinkFilterType::Pointer shrinkFilter = ShrinkFilterType::New();
    //shrinkFilter->SetInput(mask_US);
    shrinkFilter->SetInput(CTA_translated);
    int shrinkX = int(2*spacingIRM[0]/spacingUS[0]);
    int shrinkY = int(2*spacingIRM[1]/spacingUS[1]);
    int shrinkZ = int(2*spacingIRM[2]/spacingUS[2]);
    
    cout<<"shrinking factors : "<<endl;
    cout<<shrinkX<<" "<<shrinkY<<" "<<shrinkZ<<endl;
    
    shrinkFilter->SetShrinkFactor(0, shrinkX);
    shrinkFilter->SetShrinkFactor(1, shrinkY);
    shrinkFilter->SetShrinkFactor(2, shrinkZ);
    try {
        shrinkFilter->Update();
    } catch (itk::ExceptionObject &e) {
        cerr<<"error while downsampling US image"<<endl;
        cerr<<e<<endl;
        return EXIT_FAILURE;
    }
    
    //BinaryImageType::Pointer shrunk_mask_US = shrinkFilter->GetOutput();
    ImageType::Pointer shrunk_CTA_US = shrinkFilter->GetOutput();
    
    WriterType::Pointer writer_shrunk_US = WriterType::New();
    writer_shrunk_US->SetImageIO(n_io);
    string out_shrunk_US = outputPath+"/shrunk_US.nii.gz";
    writer_shrunk_US->SetFileName(out_shrunk_US);
    writer_shrunk_US->SetInput(shrunk_CTA_US);
    try {
        writer_shrunk_US->Update();
    } catch (itk::ExceptionObject &e) {
        cerr<<"error whilte writing shrunk US image"<<endl;
        cerr<<e<<endl;
        return EXIT_FAILURE;
    }
    




    

    
    
//    /*********************************
//     * REGISTRATION BASED ON VESSELS
//     **********************************/
//    
//    //declare types for all elements
//    
//    // 1. transform
//    //  The transform that will map the fixed image into the moving image.
//    //HERE, we have already declared the transformtype earlier for the gross alignment
//    //maybe we should try affine !
//    
//    //2. optimizer
//    //  An optimizer is required to explore the parameter space of the transform
//    //  in search of optimal values of the metric.
//    
//    //typedef itk::AmoebaOptimizer OptimizerType;
//    //
//    typedef itk::LBFGSBOptimizer2 OptimizerType;
//    
//    //typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
//    
//    
//    //3. metric
//    //  The metric will compare how well the two images match each other. Metric
//    //  types are usually parameterized by the image types as it can be seen in
//    //  the following type declaration.
//    
//    typedef itk::NormalizedCorrelationImageToImageMetric2<ImageType, ImageType> MetricType; //si vesselnesss et maskUS
//    
//    //typedef itk::NormalizedCorrelationImageToImageMetric<BinaryImageType, BinaryImageType> MetricType;
//    //typedef itk::MeanSquaresImageToImageMetric<BinaryImageType, BinaryImageType> MetricType;
//    
//    //typedef itk::KappaStatisticImageToImageMetric<BinaryImageType, BinaryImageType> MetricType;
//    //typedef itk::MeanReciprocalSquareDifferenceImageToImageMetric<BinaryImageType, BinaryImageType> MetricType;
//    
//    //4. interpolator
//    //  Finally, the type of the interpolator is declared. The interpolator will
//    //  evaluate the intensities of the moving image at non-grid positions.
//    
//    typedef itk::LinearInterpolateImageFunction<ImageType,double> InterpolatorType;
//    
//    //typedef itk::LinearInterpolateImageFunction<BinaryImageType,double> InterpolatorType;
//    //typedef itk::NearestNeighborInterpolateImageFunction<BinaryImageType> InterpolatorType;
//    
//    
//    //framework
//    //  The registration method type is instantiated using the types of the
//    //  fixed and moving images. This class is responsible for interconnecting
//    //  all the components that we have described so far.
//    
//    typedef itk::ImageRegistrationMethod<ImageType,ImageType> RegistrationType;
//
//    
//    //typedef itk::ImageRegistrationMethod<BinaryImageType, BinaryImageType> RegistrationType;
//    //typedef itk::ImageRegistrationMethod<ImageType, BinaryImageType> RegistrationType;
//    
//    //create components
//    
//    MetricType::Pointer metric = MetricType::New();
//    //metric->SetSubtractMean(true);
//    EulerTransformType::Pointer transform = EulerTransformType::New();
//    OptimizerType::Pointer optimizer = OptimizerType::New();
//    InterpolatorType::Pointer interpolator = InterpolatorType::New();
//    RegistrationType::Pointer registration = RegistrationType::New();
//    
//    //connect component to framework
//    registration->SetMetric(metric);
//    registration->SetOptimizer(optimizer);
//    registration->SetTransform(transform);
//    registration->SetInterpolator(interpolator);
//    
//    //set the source and target images
//    registration->SetFixedImage(vesselImageR);
//    //registration->SetFixedImage(vesselImage);
//    registration->SetMovingImage(shrunk_CTA_US);
//    
//    // here limit the region to the smallest bounding box around the vessels !-> maybe you should crop the mask_veins image !
//    registration->SetFixedImageRegion(vesselImageR->GetLargestPossibleRegion());
//    
//    //initialize transform
//    
//    typedef RegistrationType::ParametersType ParametersType;
//    ParametersType initialParameters2(transform->GetNumberOfParameters());
//    initialParameters2[0] =0.0;
//    initialParameters2[1] =0;
//    initialParameters2[2] =0;
//    initialParameters2[3] =0;
//    initialParameters2[4] =0;
//    initialParameters2[5] =0;
//
//    registration->SetInitialTransformParameters(initialParameters2);
//    
////     //Amoeba
////    OptimizerType::ParametersType simplexDelta(transform->GetNumberOfParameters());
////    simplexDelta[0] =0.4;
////    simplexDelta[1] =0.4;
////    simplexDelta[2] =0.4;
////    simplexDelta[3] = 30;
////    simplexDelta[4] = 30;
////    simplexDelta[5] = 30;
////    
////    optimizer->AutomaticInitialSimplexOff();
////    optimizer->SetInitialSimplexDelta(simplexDelta);
////    
////    //tolerance
////    optimizer->SetParametersConvergenceTolerance(0.1);
////    optimizer->SetFunctionConvergenceTolerance(0.01);
////    optimizer->SetMaximumNumberOfIterations(200);
//    
//    //LBFGSB
//    
//    OptimizerType::BoundSelectionType boundSelect(transform->GetNumberOfParameters());
//    OptimizerType::BoundValueType upperBound(transform->GetNumberOfParameters());
//    OptimizerType::BoundValueType lowerBound(transform->GetNumberOfParameters());
//    
//    boundSelect.Fill(2); //2 means lower and upper bounds for all parameters
//    
//    upperBound[0] = 0.4;
//    upperBound[1] = 0.4;
//    upperBound[2] = 0.4;
//    upperBound[3] = 30;
//    upperBound[4] = 30;
//    upperBound[5] = 30;
//    
//    lowerBound[0]= -0.4;
//    lowerBound[1]= -0.4;
//    lowerBound[2]= -0.4;
//    lowerBound[3]= -30;
//    lowerBound[4]= -30;
//    lowerBound[5]= -30;
//    
//    optimizer->SetBoundSelection(boundSelect);
//    optimizer->SetUpperBound(upperBound);
//    optimizer->SetLowerBound(lowerBound);
//    
//    optimizer->SetCostFunctionConvergenceFactor(1e+7);
//    //optimizer->SetProjectedGradientTolerance(1.0);
//    optimizer->SetMaximumNumberOfIterations(200);
//    optimizer->SetMaximumNumberOfEvaluations(200);
//    optimizer->SetMaximumNumberOfCorrections(5);
////    optimizer->SetXTol(10);
////    optimizer->SetFTol(1.0);
//    
//
//
//    
//    CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
//    optimizer->AddObserver(itk::IterationEvent(), observer);
//    
//    cout<<"registration"<<endl;
//    
//    try {
//        registration->Update();
//        cout<<"Optimizer stop condition : "<<registration->GetOptimizer()->GetStopConditionDescription()<<endl;
//
//    } catch (itk::ExceptionObject &e) {
//        cerr<<"erreur dans le recalage"<<endl;
//        cerr<<e<<endl;
//        return EXIT_FAILURE;
//    }
//
//    ParametersType finalParameters = registration->GetLastTransformParameters();
//    
//    double bestValue = optimizer->GetValue();
//    
//    cout<<"Results : "<<endl;
//    cout<<"parametres : "<<finalParameters<<endl;
//    cout<<"metric value : "<<bestValue<<endl;

    
    /*
     * TEST BOBYQA
     */
    
    cout<<"registration routine"<<endl;
    
    matrix<double> initialParameters (6,1); //values didvided by 0.3 and 30 and * radius
    initialParameters(0) = 0;
    initialParameters(1) = 0;
    initialParameters(2) = 0;
    initialParameters(3) = 0;
    initialParameters(4) = 0;
    initialParameters(5) = 0;
    
    //the bobyqa function modifies the X values if it comes two close ( within rhobegin) of the bounds
    
    matrix<double> x_lower (6,1); //-0.5 rot, -10 trans
    x_lower(0) = -1;
    x_lower(1) = -1;
    x_lower(2) = -1;
    x_lower(3) = -1;
    x_lower(4) = -1;
    x_lower(5) = -1;
    
    matrix<double> x_upper (6,1); //0.5 rot, 10 trans
    x_upper(0) = 1;
    x_upper(1) = 1;
    x_upper(2) = 1;
    x_upper(3) = 1;
    x_upper(4) = 1;
    x_upper(5) = 1;
    
    double m = 13;
    
    //rho begin
    double radius = 0.9; //0.9 //zero makes no sense !!!
    
    //rho end
    double precision = 0.01; //0.01 0.001
    
    //niter
    double nombreIteration = 200;
    
    
//    //METRIC SETTING
//    //metric
//    typedef itk::NormalizedCorrelationImageToImageMetric2<ImageType, ImageType> MetricType; //si vesselnesss et maskUS
//    //interpolator
//    typedef itk::LinearInterpolateImageFunction<ImageType,double> InterpolatorType;
//    
//    MetricType::Pointer metric = MetricType::New();
//    InterpolatorType::Pointer interpolator = InterpolatorType::New();
//    EulerTransformType::Pointer transform = EulerTransformType::New();
//    metric->SetInterpolator(interpolator);
//    metric->SetTransform(transform);
//    metric->SetMovingImage(shrunk_CTA_US);
//    metric->SetFixedImage(vesselImageR);
    
    
    
    NCC_function NCC = NCC_function(shrunk_CTA_US, shrunk_mask_IRM);
    NCC.setMaxRot(0.3);
    NCC.setMaxTrans(10);
    NCC.setRadius(radius);
    //NCC.setSubstractMeanBool(true);
    
    //double radius = 0.9
//    double maxRot = 0.2;
//    double maxTrans = 10;
    
    //cout<<"initialParameters"<<initialParameters<<endl;
    
    double best_score = find_min_bobyqa(NCC, initialParameters, m, x_lower, x_upper, radius, precision, nombreIteration);
    

    cout<<"Best score : "<<best_score<<endl;
//
//    //transformation de l'image
   EulerTransformType::Pointer finalTsf = EulerTransformType::New();
//    
    EulerTransformType::ParametersType finalParameters(6);
    finalParameters[0] = initialParameters(0)*NCC.GetMaxRot()/(NCC.GetRadius)();
    finalParameters[1] = initialParameters(1)*NCC.GetMaxRot()/(NCC.GetRadius)();
    finalParameters[2] = initialParameters(2)*NCC.GetMaxRot()/(NCC.GetRadius)();
    finalParameters[3] = initialParameters(3)*NCC.GetMaxTrans()/(NCC.GetRadius)();
    finalParameters[4] = initialParameters(4)*NCC.GetMaxTrans()/(NCC.GetRadius)();
    finalParameters[5] = initialParameters(5)*NCC.GetMaxTrans()/(NCC.GetRadius)();
  finalTsf->SetParameters(finalParameters);
//    
   cout<<"final parameters : "<<finalTsf->GetParameters()<<endl;
//    
    ImageType::SizeType sizeUS2 = US_translated->GetLargestPossibleRegion().GetSize();
    ImageType::PointType origin2 = US_translated->GetOrigin();
    ImageType::SpacingType spacing2 = US_translated->GetSpacing();
    ImageType::PointType center2;
    center2[0] = origin2[0]+spacing2[0]*sizeUS2[0]/2;
    center2[1] = origin2[1]+spacing2[1]*sizeUS2[1]/2;
    center2[2] = origin2[2]+spacing2[2]*sizeUS2[2]/2;
    
    
    EulerTransformType::ParametersType eulerFixedParameters2(3);
    eulerFixedParameters2[0] =center2[0];
    eulerFixedParameters2[1] =center2[1];
    eulerFixedParameters2[2] =center2[2];
    
    finalTsf->SetFixedParameters(eulerFixedParameters2);
    
    cout<<"Write results"<<endl;
    
    cout<<"Enregistrement dans fichier parametres : "<<endl;
    
    string outP = outputPath+"/parameters_initialisation.txt";
    ofstream fichier(outP.c_str(),ios::out | ios::trunc);
    
    if(fichier)
    {
        fichier<<"Scaling factors rotation and translation : "<<endl;
        fichier<<"Rotation : "<<NCC.GetMaxRot()<<endl;
        fichier<<"Translation : "<<NCC.GetMaxTrans()<<endl;
        fichier<<"Radius : "<< NCC.GetRadius()<<endl;
        fichier<<"Precision "<<precision<<endl;
        fichier<<"Parameters for rigid transform : "<<endl;
        fichier<<finalTsf->GetParameters()<<endl;
        fichier<<" Score for this position : "<<endl;
        fichier<<best_score<<endl;
        fichier.close();
    }
    
    else
    {
        cerr<<"Error in opening txt file for parameters"<<endl;
    }
    
//
    ResampleFilterType::Pointer resampler1 = ResampleFilterType::New();
    resampler1->SetInput(CTA_translated);
    resampler1->SetTransform(finalTsf);
    resampler1->SetSize(CTA_translated->GetLargestPossibleRegion().GetSize());
    resampler1->SetOutputSpacing(CTA_translated->GetSpacing());
    resampler1->SetOutputDirection(finalTsf->GetInverseMatrix()*CTA_translated->GetDirection());
    resampler1->SetOutputOrigin(finalTsf->GetInverseTransform()->TransformPoint(CTA_translated->GetOrigin()));
    
    ImageType::Pointer finalVesselUS = resampler1->GetOutput();
    
    
    WriterType::Pointer writer6 =  WriterType::New();
    writer6->SetImageIO(n_io);
    string out6 = outputPath+"/registered_CTA.nii.gz";
    writer6->SetFileName(out6);
    writer6->SetInput(finalVesselUS);
    try {
        writer6->Update();
    } catch (itk::ExceptionObject &e) {
        cerr<<"error whilte writing registered image"<<endl;
        cerr<<e<<endl;
        return EXIT_FAILURE;
    }
    
    ResampleFilterType::Pointer resampler3 = ResampleFilterType::New();
    resampler3->SetInput(US_translated);
    resampler3->SetTransform(finalTsf);
    resampler3->SetSize(US_translated->GetLargestPossibleRegion().GetSize());
    resampler3->SetOutputSpacing(US_translated->GetSpacing());
    resampler3->SetOutputDirection(finalTsf->GetInverseMatrix()*US_translated->GetDirection());
    resampler3->SetOutputOrigin(finalTsf->GetInverseTransform()->TransformPoint(US_translated->GetOrigin()));
    
    ImageType::Pointer finalUS = resampler3->GetOutput();
    
    
    WriterType::Pointer writer7 =  WriterType::New();
    writer7->SetImageIO(n_io);
    string out7 = outputPath+"/registered_US.nii.gz";
    writer7->SetFileName(out7);
    writer7->SetInput(finalUS);
    try {
        writer7->Update();
    } catch (itk::ExceptionObject &e) {
        cerr<<"error whilte writing registered image"<<endl;
        cerr<<e<<endl;
        return EXIT_FAILURE;
    }
    


    tend = std::time(NULL);
    texec = std::difftime(tend,tbegin);
    
    cout<<"temps en s : "<<texec<<endl;
    
    return 0;
}

//    /**********************
//     * MULTI RES REGISTRATION
//     *****************************/
//
//    //declare types for all elements
//
//    // 1. transform
//    //  The transform that will map the fixed image into the moving image.
//    //HERE, we have already declared the transformtype earlier for the gross alignment
//    //maybe we should try affine !
//
//    //2. optimizer
//    //  An optimizer is required to explore the parameter space of the transform
//    //  in search of optimal values of the metric.
//    //typedef itk::AmoebaOptimizer OptimizerType;
//    typedef itk::LBFGSBOptimizer OptimizerType;
//
//    //3. metric
//    //  The metric will compare how well the two images match each other. Metric
//    //  types are usually parameterized by the image types as it can be seen in
//    //  the following type declaration.
//
//    typedef itk::NormalizedCorrelationImageToImageMetric<BinaryImageType, BinaryImageType> MetricType;
//
//    //4. interpolator
//    //  Finally, the type of the interpolator is declared. The interpolator will
//    //  evaluate the intensities of the moving image at non-grid positions.
//
//    typedef itk::LinearInterpolateImageFunction<BinaryImageType,double> InterpolatorType;
//
//    //5. Framework
//
//    typedef itk::MultiResolutionImageRegistrationMethod<BinaryImageType, BinaryImageType> MultiResRegistrationType;
//
//    //6. Pyramids
//    //create the multi-res pyramids for the two images
//    typedef itk::MultiResolutionPyramidImageFilter<BinaryImageType, BinaryImageType> PyramidType;
//
//    //create components
//    MetricType::Pointer metric = MetricType::New();
//    EulerTransformType::Pointer transform = EulerTransformType::New();
//    OptimizerType::Pointer optimizer = OptimizerType::New();
//    InterpolatorType::Pointer interpolator = InterpolatorType::New();
//    MultiResRegistrationType::Pointer registration = MultiResRegistrationType::New();
//
//    PyramidType::Pointer fixedImagePyramid = PyramidType::New();
//    PyramidType::Pointer movingImagePyramid = PyramidType::New();
//
//    //connect component to framework
//    registration->SetMetric(metric);
//    registration->SetOptimizer(optimizer);
//    registration->SetTransform(transform);
//    registration->SetInterpolator(interpolator);
//    registration->SetFixedImagePyramid(fixedImagePyramid);
//    registration->SetMovingImagePyramid(movingImagePyramid);
//    registration->SetFixedImage(mask_veins);
//    registration->SetMovingImage(mask_US);
//    registration->SetFixedImageRegion(mask_veins->GetLargestPossibleRegion());
//
//    //initialize transform
//
//    typedef MultiResRegistrationType::ParametersType ParametersType;
//    ParametersType initialParameters(transform->GetNumberOfParameters());
//    initialParameters[0] =0;
//    initialParameters[1] =0;
//    initialParameters[2] =0;
//    initialParameters[3] =0;
//    initialParameters[4] =0;
//    initialParameters[5] =0;
//
//    registration->SetInitialTransformParameters(initialParameters);
//
//    //setting parameters for optimizer
//    //Amoeba
////    OptimizerType::ParametersType simplexDelta(transform->GetNumberOfParameters());
////    simplexDelta[0] =0.3;
////    simplexDelta[1] =0.3;
////    simplexDelta[2] =0.3;
////    simplexDelta[3] = 20;
////    simplexDelta[4] = 20;
////    simplexDelta[5] = 20;
////
////    optimizer->AutomaticInitialSimplexOff();
////    optimizer->SetInitialSimplexDelta(simplexDelta);
////
////    //tolerance
////    optimizer->SetParametersConvergenceTolerance(1.0);
////    optimizer->SetFunctionConvergenceTolerance(0.1);
////
////    optimizer->SetMaximumNumberOfIterations(200);
//
//        //LBFGSB
//
//        OptimizerType::BoundSelectionType boundSelect(transform->GetNumberOfParameters());
//        OptimizerType::BoundValueType upperBound(transform->GetNumberOfParameters());
//        OptimizerType::BoundValueType lowerBound(transform->GetNumberOfParameters());
//
//        boundSelect.Fill(2); //2 means lower and upper bounds for all parameters
//
//        upperBound[0] = 0.4;
//        upperBound[1] = 0.4;
//        upperBound[2] = 0.4;
//        upperBound[3] = 30;
//        upperBound[4] = 30;
//        upperBound[5] = 30;
//
//        lowerBound[0]= -0.4;
//        lowerBound[1]= -0.4;
//        lowerBound[2]= -0.4;
//        lowerBound[3]= -30;
//        lowerBound[4]= -30;
//        lowerBound[5]= -30;
//
//        optimizer->SetBoundSelection(boundSelect);
//        optimizer->SetUpperBound(upperBound);
//        optimizer->SetLowerBound(lowerBound);
//
//        optimizer->SetCostFunctionConvergenceFactor(1e+7);
//        //optimizer->SetProjectedGradientTolerance(1.0);
//        optimizer->SetMaximumNumberOfIterations(200);
//        optimizer->SetMaximumNumberOfEvaluations(200);
//        optimizer->SetMaximumNumberOfCorrections(5);
//
//
//
//
//    //observer on the optimisation process
//    CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
//    optimizer->AddObserver(itk::IterationEvent(), observer);
//
//    //observer for multi res reg process
//    typedef RegistrationInterfaceCommand<MultiResRegistrationType> CommandType;
//    CommandType::Pointer command = CommandType::New();
//    registration->AddObserver(itk::IterationEvent(), command);
//
//    registration->SetNumberOfLevels(3);
//
//    try {
//        registration->Update();
//        std::cout << "Optimizer stop condition: "
//        << registration->GetOptimizer()->GetStopConditionDescription()
//        << std::endl;
//    } catch (itk::ExceptionObject e) {
//        std::cout << "ExceptionObject caught !" << std::endl;
//        std::cout << e << std::endl;
//        return EXIT_FAILURE;
//    }
//
//    ParametersType finalParameters = registration->GetLastTransformParameters();
//
//    double bestValue = optimizer->GetValue();
//
//    cout<<"Results : "<<endl;
//    cout<<"parametres : "<<finalParameters<<endl;
//    cout<<"metric value : "<<bestValue<<endl;
//
//    //transformation de l'image
//    EulerTransformType::Pointer finalTsf = EulerTransformType::New();
//    finalTsf->SetParameters(finalParameters);
//
////    ResampleFilterType::Pointer resampler1 = ResampleFilterType::New();
////    resampler1->SetInput(vesselImageUS);
////    resampler1->SetTransform(finalTsf);
////    resampler1->SetSize(vesselImageUS->GetLargestPossibleRegion().GetSize());
////    resampler1->SetOutputSpacing(vesselImageUS->GetSpacing());
////    resampler1->SetOutputDirection(finalTsf->GetInverseMatrix()*vesselImageUS->GetDirection());
////    resampler1->SetOutputOrigin(finalTsf->GetInverseTransform()->TransformPoint(vesselImageUS->GetOrigin()));
////
////    ImageType::Pointer finalVesselUS = resampler1->GetOutput();
////
////
////    WriterType::Pointer writer6 =  WriterType::New();
////    writer6->SetImageIO(n_io);
////    string out6 = outputPath+"/registeredVessel_US.nii.gz";
////    writer6->SetFileName(out6);
////    writer6->SetInput(finalVesselUS);
////    try {
////        writer6->Update();
////    } catch (itk::ExceptionObject &e) {
////        cerr<<"error whilte writing registered image"<<endl;
////        cerr<<e<<endl;
////        return EXIT_FAILURE;
////    }
//
//    ResampleFilterType::Pointer resampler2 = ResampleFilterType::New();
//    resampler2->SetInput(US_translated);
//    resampler2->SetTransform(finalTsf);
//    resampler2->SetSize(US_translated->GetLargestPossibleRegion().GetSize());
//    resampler2->SetOutputSpacing(US_translated->GetSpacing());
//    resampler2->SetOutputDirection(finalTsf->GetInverseMatrix()*US_translated->GetDirection());
//    resampler2->SetOutputOrigin(finalTsf->GetInverseTransform()->TransformPoint(US_translated->GetOrigin()));
//
//    ImageType::Pointer finalUS = resampler2->GetOutput();
//
//
//    WriterType::Pointer writer7 =  WriterType::New();
//    writer7->SetImageIO(n_io);
//    string out7 = outputPath+"/registered_US_multires.nii.gz";
//    writer7->SetFileName(out7);
//    writer7->SetInput(finalUS);
//    try {
//        writer7->Update();
//    } catch (itk::ExceptionObject &e) {
//        cerr<<"error whilte writing registered image"<<endl;
//        cerr<<e<<endl;
//        return EXIT_FAILURE;
//    }



//    /*****************************
//     * Equalisation d'histogramme
//     *****************************/
//
//    HistogramEqualizerType::Pointer equalizer = HistogramEqualizerType::New();
//    equalizer->SetInput(image_US);
//    equalizer->SetRadius(1);
//    equalizer->SetAlpha(0);
//    equalizer->Update();
//
//    ImageType::Pointer US_equalized = equalizer->GetOutput();
//
//    typename WriterType::Pointer writer1 = WriterType::New();
//    string out1 = outputPath+"/EqualizedUS.nii.gz";
//    itk::NiftiImageIO::Pointer io = itk::NiftiImageIO::New();
//    writer1->SetInput(US_equalized);
//    writer1->SetImageIO(io);
//    writer1->SetFileName(out1);
//    try {
//        writer1->Update();
//    } catch (itk::ExceptionObject &e) {
//        cerr<<"error while writing image file"<<endl;
//        cerr<<e<<endl;
//        EXIT_FAILURE;
//    }
//
//    cout<<"done writing equalized image"<<endl;



//
//    /*******************
//     * TEST GRADIENT US
//     ********************/
//
//    cout<<"compute gradient image of US"<<endl;
//
//    GradientFilterType::Pointer filterG = GradientFilterType::New();
//    filterG->SetInput(image_US);
//    try {
//        filterG->Update();
//    } catch (itk::ExceptionObject &e) {
//        std::cerr<<"error while computing gradient image"<<std::endl;
//        std::cerr<<e<<std::endl;
//        EXIT_FAILURE;
//    }
//
//    ImageType::Pointer US_grad = filterG->GetOutput();
//
//
//    typename WriterType::Pointer writer1 = WriterType::New();
//    string out1 = outputPath+"/testgradUS.nii.gz";
//    itk::NiftiImageIO::Pointer io = itk::NiftiImageIO::New();
//    writer1->SetInput(US_grad);
//    writer1->SetImageIO(io);
//    writer1->SetFileName(out1);
//    try {
//        writer1->Update();
//    } catch (itk::ExceptionObject &e) {
//        cerr<<"error while writing image file"<<endl;
//        cerr<<e<<endl;
//        EXIT_FAILURE;
//    }
//
//    cout<<"done writing gradient image"<<endl;

/******************
 * Thresholding
 *******************/

//    BinaryThresholdFilterType::Pointer thresholder = BinaryThresholdFilterType::New();
//    thresholder->SetInput(image_US);
//    thresholder->SetInsideValue(255);
//    thresholder->SetOutsideValue(0);
//    thresholder->SetLowerThreshold(0);
//    thresholder->SetUpperThreshold(20);
//
//    try {
//        thresholder->Update();
//    } catch (itk::ExceptionObject &e) {
//        std::cerr<<"error while binarizing US image"<<std::endl;
//        std::cerr<<e<<std::endl;
//    }
//
//    BinaryImageType::Pointer mask = thresholder->GetOutput();
//
//    //writing mask images
//
//    BinaryWriterType::Pointer writer3 = BinaryWriterType::New();
//    std::string out3 = outputPath+"/testmaskUS2.nii.gz";
//    //itk::NiftiImageIO::Pointer io = itk::NiftiImageIO::New();
//    writer3->SetInput(mask);
//    writer3->SetImageIO(io);
//    writer3->SetFileName(out3);
//    try {
//        writer3->Update();
//    } catch (itk::ExceptionObject &e) {
//        std::cerr<<"error while writing image file"<<std::endl;
//        std::cerr<<e<<std::endl;
//        EXIT_FAILURE;
//    }
//
//    std::cout<<"done writing mask image"<<std::endl;


//
//    /**************************
//     * TEST VESSELNESS PORTAL
//     **************************/
//
//    //vesselness sur IRM
//    cout<<"vesselness image IRM"<<endl;
//
//    typedef itk::SymmetricSecondRankTensor<double,3> HessianPixelType;
//    typedef itk::Image<HessianPixelType,3> HessianImageType;
//    typedef itk::HessianToObjectnessMeasureImageFilter<HessianImageType, ImageType> ObjectnessFilterType;
//
//    ObjectnessFilterType::Pointer vesselnessFilter = ObjectnessFilterType::New();
//    vesselnessFilter->SetBrightObject(true);
//    vesselnessFilter->SetScaleObjectnessMeasure(false);
//    vesselnessFilter->SetObjectDimension(1);
//    vesselnessFilter->SetAlpha(0.5);
//    vesselnessFilter->SetBeta(0.5);
//    vesselnessFilter->SetGamma(10.0);
//    //vesselnessFilter->GetInput()->GetN
//
//    typedef itk::MultiScaleHessianBasedMeasureImageFilter<ImageType, HessianImageType,ImageType> MultiScaleEnhancementFilterType;
//
//    MultiScaleEnhancementFilterType::Pointer enhancer = MultiScaleEnhancementFilterType::New();
//    enhancer->SetInput(image_IRM);
//    enhancer->SetHessianToMeasureFilter(vesselnessFilter);
//    enhancer->SetSigmaStepMethodToLogarithmic();
//    enhancer->SetSigmaMaximum(2.5); //3.0
//    enhancer->SetSigmaMinimum(2.0); //2.5
//    enhancer->SetNumberOfSigmaSteps(5);
//
//    ImageType::Pointer vesselImage = enhancer->GetOutput();
//
//    WriterType::Pointer writer4 =  WriterType::New();
//    writer4->SetImageIO(n_io);
//    string out4 = outputPath+"/vesselnessIRMImage.nii.gz";
//    writer4->SetFileName(out4);
//    writer4->SetInput(vesselImage);
//    try {
//        writer4->Update();
//    } catch (itk::ExceptionObject &e) {
//        cerr<<"error whilte writing registered image"<<endl;
//        cerr<<e<<endl;
//        return EXIT_FAILURE;
//    }
//
//
//    //test vesselness sur US
//
//    cout<<"test vesselnessImage US"<<endl;
//
//    //should we downsample first ?
//
//    ObjectnessFilterType::Pointer vesselnessFilterUS = ObjectnessFilterType::New();
//    vesselnessFilterUS->SetBrightObject(false);
//    vesselnessFilterUS->SetScaleObjectnessMeasure(false);
//    vesselnessFilterUS->SetAlpha(0.5);
//    vesselnessFilterUS->SetBeta(0.5);
//    vesselnessFilterUS->SetGamma(1.0);
//
////    VesselnessMeasureFilterType::Pointer vesselnessFilter = VesselnessMeasureFilterType::New();
////    vesselnessFilter->SetAlpha1(1.5);
////    vesselnessFilter->SetAlpha2(1.5);//sato doesn't work here cause the structure we try to enhance are dark
//
//
//
//    MultiScaleEnhancementFilterType::Pointer enhancerUS = MultiScaleEnhancementFilterType::New();
//    enhancerUS->SetInput(US_translated);
//    enhancerUS->SetHessianToMeasureFilter(vesselnessFilterUS);
//    enhancerUS->SetSigmaStepMethodToLogarithmic();
//    enhancerUS->SetSigmaMaximum(4.0); //3.0
//    enhancerUS->SetSigmaMinimum(3.0); //2.5
//    enhancerUS->SetNumberOfSigmaSteps(3);
//
//    ImageType::Pointer vesselImageUS = enhancerUS->GetOutput();
//
//    WriterType::Pointer writer5 =  WriterType::New();
//    writer5->SetImageIO(n_io);
//    string out5 = outputPath+"/vesselnessImage_US.nii.gz";
//    writer5->SetFileName(out5);
//    writer5->SetInput(vesselImageUS);
//    try {
//        writer5->Update();
//    } catch (itk::ExceptionObject &e) {
//        cerr<<"error whilte writing registered image"<<endl;
//        cerr<<e<<endl;
//        return EXIT_FAILURE;
//    }

//    /*************************
//     * MORPHOLOGICAL CLOSING
//     *************************/
//
//    cout<<"dilatation"<<endl;
//
//    StructuringElementType kernel;
//    BinaryImageType::SizeType radius;
//    radius[0] = 3;
//    radius[1] = 3;
//    radius[2] = 3;
//    kernel.SetRadius(radius);
//    kernel.CreateStructuringElement();
//
////    DilateFilterType::Pointer dilateFilter = DilateFilterType::New();
////    dilateFilter->SetInput(mask_US);
////    dilateFilter->SetKernel(kernel);
//
//    ClosingFilterType::Pointer closingFilter = ClosingFilterType::New();
//    closingFilter->SetInput(mask_US);
//    closingFilter->SetKernel(kernel);
//    closingFilter->Update();
//
//
//    BinaryImageType::Pointer US_dilated = closingFilter->GetOutput();
//
//    BinaryWriterType::Pointer binWriter2 = BinaryWriterType::New();
//    binWriter2->SetImageIO(n_io);
//    string out_bin2 = outputPath+"/dilated_mask_US.nii.gz";
//    binWriter2->SetInput(US_dilated);
//    binWriter2->SetFileName(out_bin2);
//    try {
//        binWriter2->Update();
//    } catch (itk::ExceptionObject &e) {
//        cerr<<"error whilte writing registered image"<<endl;
//        cerr<<e<<endl;
//        return EXIT_FAILURE;
//    }
///*************************
// * THRESHOLDING
// ***********************/
//
//    cout<<"binary Thresholding"<<endl;
//
//    BinaryThresholdFilterType::Pointer thresholder = BinaryThresholdFilterType::New();
//    thresholder->SetInput(US_translated);
//    thresholder->SetLowerThreshold(1.0);
//    thresholder->SetUpperThreshold(30.0);
//    thresholder->SetInsideValue(255);
//    thresholder->SetOutsideValue(0);
//
//    BinaryImageType::Pointer mask_US = thresholder->GetOutput();
//
//    BinaryWriterType::Pointer binWriter = BinaryWriterType::New();
//    binWriter->SetImageIO(n_io);
//    string out_bin = outputPath+"/Thresholded_US.nii.gz";
//    binWriter->SetInput(mask_US);
//    binWriter->SetFileName(out_bin);
//    try {
//        binWriter->Update();
//    } catch (itk::ExceptionObject &e) {
//        cerr<<"error whilte writing registered image"<<endl;
//        cerr<<e<<endl;
//        return EXIT_FAILURE;
//    }
//

