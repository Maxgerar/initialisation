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

//classic registration

#include "itkImageRegistrationMethod.h"
#include "itkNormalizedCorrelationImageToImageMetric.h"
#include "itkKappaStatisticImageToImageMetric.h"
#include "itkMeanReciprocalSquareDifferenceImageToImageMetric.h"
#include "itkAmoebaOptimizer.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkCommand.h"

#include "itkLBFGSBOptimizer.h"

//mutliresolution registration
#include "itkMultiResolutionImageRegistrationMethod.h"
#include "itkMultiResolutionPyramidImageFilter.h"




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
typedef itk::ShrinkImageFilter<BinaryImageType, BinaryImageType> ShrinkFilterType;


typedef itk::BinaryThresholdImageFilter<ImageType, BinaryImageType> BinaryThresholdFilterType;
typedef itk::AdaptiveHistogramEqualizationImageFilter<ImageType> HistogramEqualizerType;

typedef itk::HessianToObjectnessMeasureImageFilter<HessianImageType, ImageType> ObjectnessFilterType;
typedef itk::Hessian3DToVesselnessMeasureImageFilter<double> VesselnessMeasureFilterType;
typedef itk::MultiScaleHessianBasedMeasureImageFilter<ImageType, HessianImageType,ImageType> MultiScaleEnhancementFilterType;
using namespace std;

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
//    typedef itk::GradientDescentOptimizer OptimizerType;
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
//            optimizer->SetFunctionConvergenceTolerance(1.0);
//            optimizer->SetParametersConvergenceTolerance(1.0);
//            //optimizer->SetMaximumStepLength( 16.00 );
//            //optimizer->SetMinimumStepLength( 0.01 );
//        }
//        else
//        {
//            optimizer->SetFunctionConvergenceTolerance(optimizer->GetFunctionConvergenceTolerance()*0.1);
//            optimizer->SetParametersConvergenceTolerance(optimizer->GetParametersConvergenceTolerance()*0.25);
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
    typedef itk::LBFGSBOptimizer OptimizerType;
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
//        std::cout << optimizer->GetCachedValue() << "   "<<endl;
//        std::cout << optimizer->GetCachedCurrentPosition() << std::endl;
        
        //std::cout<<optimizer->GetCurrentIteration()<<"  ";
        std::cout<<optimizer->GetCachedValue()<<endl;
        //std::cout<<optimizer->GetInfinityNormOfProjectedGradient()<<endl;
        std::cout<<optimizer->GetCachedCurrentPosition()<<endl;
  
    }
private:
    unsigned long m_IterationNumber;
};


int main(int argc, const char * argv[]) {
    
    //to determine time of computation
    std::srand(time(NULL));
    std::time_t tbegin,tend;
    double texec = 0;
    tbegin = std::time(NULL);
    
    string filenameUS;
    string filenameIRM;
    string filenameVeins;
    string outputPath;

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

        if(strcmp(argv[i], "-iMaskVeins")==0)
        {
            i++;
            filenameVeins= argv[i];
            cout<<"Use of mask image of liver"<<endl;
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
    
    if(filenameVeins=="")
    {
        cerr<<"input veins mask file not provided"<<endl;
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
    //ImageType::Pointer image_IRM = ImageType::New();
    
    
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
    
    
    //min max image US
    MinMaxCalculatorType::Pointer minMaxUS = MinMaxCalculatorType::New();
    minMaxUS->SetImage(image_US);
    minMaxUS->Compute();
    
    cout<<"intensity range US image : "<<"[ "<<minMaxUS->GetMinimum()<<","<<minMaxUS->GetMaximum()<<" ]"<<endl;
    
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
    
    
    //min max image US
    MinMaxCalculatorType::Pointer minMaxIRM = MinMaxCalculatorType::New();
    minMaxIRM->SetImage(image_IRM);
    minMaxIRM->Compute();
    
    cout<<"intensity range MRI image : "<<"[ "<<minMaxIRM->GetMinimum()<<","<<minMaxIRM->GetMaximum()<<" ]"<<endl;
    
    BinaryImageType::Pointer mask_veins = BinaryImageType::New();
    
    BinaryReaderType::Pointer readerBin = BinaryReaderType::New();
    readerBin->SetImageIO(n_io);
    readerBin->SetFileName(filenameVeins);
    try {
        readerBin->Update();
    } catch (itk::ExceptionObject &e) {
        cerr<<"Error while reading MRI image"<<endl;
        cerr<<e<<endl;
        EXIT_FAILURE;
    }
    
    mask_veins = readerBin->GetOutput();
    cout<<"test lecture mask veins "<<endl;
    cout<<"dimensions mask veins : "<<mask_veins->GetLargestPossibleRegion().GetSize()<<endl;
    
    
    
    /**********************************************************
     * ORIENTATION : GROSS INITIAL ALIGNMENT A PRIORI KNOWLEDGE
     ***********************************************************/
    //aligning image origins
    
    cout<<"aligning image origins"<<endl;
    
    //euler tsf
    
    //determining image centers
    //MRI
    ImageType::SizeType sizeIRM = image_IRM->GetLargestPossibleRegion().GetSize();
    ImageType::PointType origineIRM = image_IRM->GetOrigin();
    ImageType::SpacingType spacingIRM = image_IRM->GetSpacing();
    ImageType::PointType centerIRM;
    
    centerIRM[0] = origineIRM[0]+spacingIRM[0]*sizeIRM[0]/2;
    centerIRM[1] = origineIRM[1]+spacingIRM[1]*sizeIRM[2]/2;
    centerIRM[2] = origineIRM[2]+spacingIRM[2]*sizeIRM[1]/2;
    
    cout<<"center MRI : "<<centerIRM<<endl;
    
    ImageType::SizeType sizeUS = image_US->GetLargestPossibleRegion().GetSize();
    ImageType::PointType origineUS = image_US->GetOrigin();
    ImageType::SpacingType spacingUS = image_US->GetSpacing();
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
    Parameters[2] = 0;//-0.1
    //on sait que le foie est situe a droite et dans la partie superieure dans l'abdomen
    Parameters[3] = -(centerIRM[0]-centerUS[0])+30;//+20;//10
    Parameters[4] = -(centerIRM[1]-centerUS[1])+10;//+10; //+5,15,0 +20 ?
    Parameters[5] = -(centerIRM[2]-centerUS[2])-25;// -50-90
    
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
    
    ImageType::Pointer US_translated = resampler->GetOutput();
    
    
    
//    WriterType::Pointer writer3 = WriterType::New();
//    string out8 =outputPath+"/translated_US_test_p1.nii.gz";
//    writer3->SetImageIO(n_io);
//    writer3->SetInput(US_translated);
//    writer3->SetFileName(out8);
//    try {
//        writer3->Update();
//    } catch (itk::ExceptionObject &e) {
//        cerr<<"error whilte writing registered image"<<endl;
//        cerr<<e<<endl;
//        return EXIT_FAILURE;
//    }
    
/*************************
 * THRESHOLDING
 ***********************/
    
    cout<<"binary Thresholding"<<endl;
    
    BinaryThresholdFilterType::Pointer thresholder = BinaryThresholdFilterType::New();
    thresholder->SetInput(US_translated);
    thresholder->SetLowerThreshold(1.0);
    thresholder->SetUpperThreshold(10.0);
    thresholder->SetInsideValue(255);
    thresholder->SetOutsideValue(0);
    
    BinaryImageType::Pointer mask_US = thresholder->GetOutput();
    
    BinaryWriterType::Pointer binWriter = BinaryWriterType::New();
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
    
    /*************
     * RES IRM
     ***************/
    ShrinkFilterType::Pointer shrinkFilter2 = ShrinkFilterType::New();
    shrinkFilter2->SetInput(mask_veins);
    
    shrinkFilter2->SetShrinkFactor(0, 2);
    shrinkFilter2->SetShrinkFactor(1, 2);
    shrinkFilter2->SetShrinkFactor(2, 2);
    try {
        shrinkFilter2->Update();
    } catch (itk::ExceptionObject &e) {
        cerr<<"error while downsampling US image"<<endl;
        cerr<<e<<endl;
        return EXIT_FAILURE;
    }
    
    BinaryImageType::Pointer shrunk_mask_veins = shrinkFilter2->GetOutput();
    
    /******************************
     * Diminuation resolution US
     *****************************/
    
    ShrinkFilterType::Pointer shrinkFilter = ShrinkFilterType::New();
    shrinkFilter->SetInput(mask_US);
    
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
    
    BinaryImageType::Pointer shrunk_mask_US = shrinkFilter->GetOutput();
    
    
    
    
    
    
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
    
    
    /*********************************
     * REGISTRATION BASED ON VESSELS
     **********************************/
    
    //declare types for all elements
    
    // 1. transform
    //  The transform that will map the fixed image into the moving image.
    //HERE, we have already declared the transformtype earlier for the gross alignment
    //maybe we should try affine !
    
    //2. optimizer
    //  An optimizer is required to explore the parameter space of the transform
    //  in search of optimal values of the metric.
    
    //typedef itk::AmoebaOptimizer OptimizerType;
    typedef itk::LBFGSBOptimizer OptimizerType;
    
    
    //3. metric
    //  The metric will compare how well the two images match each other. Metric
    //  types are usually parameterized by the image types as it can be seen in
    //  the following type declaration.
    
    typedef itk::NormalizedCorrelationImageToImageMetric<BinaryImageType, BinaryImageType> MetricType;
    //typedef itk::NormalizedCorrelationImageToImageMetric<ImageType, BinaryImageType> MetricType; //si vesselnesss et maskUS
    
    //typedef itk::KappaStatisticImageToImageMetric<BinaryImageType, BinaryImageType> MetricType;
    //typedef itk::MeanReciprocalSquareDifferenceImageToImageMetric<BinaryImageType, BinaryImageType> MetricType;
    
    //4. interpolator
    //  Finally, the type of the interpolator is declared. The interpolator will
    //  evaluate the intensities of the moving image at non-grid positions.
    
    typedef itk::LinearInterpolateImageFunction<BinaryImageType,double> InterpolatorType;
    
    
    //framework
    //  The registration method type is instantiated using the types of the
    //  fixed and moving images. This class is responsible for interconnecting
    //  all the components that we have described so far.
    
    typedef itk::ImageRegistrationMethod<BinaryImageType, BinaryImageType> RegistrationType;
    //typedef itk::ImageRegistrationMethod<ImageType, BinaryImageType> RegistrationType;
    
    //create components
    
    MetricType::Pointer metric = MetricType::New();
    EulerTransformType::Pointer transform = EulerTransformType::New();
    OptimizerType::Pointer optimizer = OptimizerType::New();
    InterpolatorType::Pointer interpolator = InterpolatorType::New();
    RegistrationType::Pointer registration = RegistrationType::New();
    
    //connect component to framework
    registration->SetMetric(metric);
    registration->SetOptimizer(optimizer);
    registration->SetTransform(transform);
    registration->SetInterpolator(interpolator);
    
    //set the source and target images
    registration->SetFixedImage(mask_veins);
    //registration->SetFixedImage(vesselImage);
    registration->SetMovingImage(shrunk_mask_US);
    
    registration->SetFixedImageRegion(mask_veins->GetLargestPossibleRegion());
    
    //initialize transform
    
    typedef RegistrationType::ParametersType ParametersType;
    ParametersType initialParameters(transform->GetNumberOfParameters());
    initialParameters[0] =0;
    initialParameters[1] =0;
    initialParameters[2] =0;
    initialParameters[3] =0;
    initialParameters[4] =0;
    initialParameters[5] =0;

    registration->SetInitialTransformParameters(initialParameters);
    
//    //setting parameters for optimizer Amoeba
//    OptimizerType::ParametersType simplexDelta(transform->GetNumberOfParameters());
//    simplexDelta[0] =0.4;
//    simplexDelta[1] =0.4;
//    simplexDelta[2] =0.4;
//    simplexDelta[3] = 20;
//    simplexDelta[4] = 20;
//    simplexDelta[5] = 20;
//    
//    optimizer->AutomaticInitialSimplexOff();
//    optimizer->SetInitialSimplexDelta(simplexDelta);
//    
//    //tolerance
//    optimizer->SetParametersConvergenceTolerance(1.0);
//    optimizer->SetFunctionConvergenceTolerance(0.01);
//    optimizer->SetMaximumNumberOfIterations(200);
    
    //LBFGSB
    
    OptimizerType::BoundSelectionType boundSelect(transform->GetNumberOfParameters());
    OptimizerType::BoundValueType upperBound(transform->GetNumberOfParameters());
    OptimizerType::BoundValueType lowerBound(transform->GetNumberOfParameters());
    
    boundSelect.Fill(2); //means lower and upper bounds for all parameters
    
    upperBound[0] = 0.4;
    upperBound[1] = 0.4;
    upperBound[2] = 0.4;
    upperBound[3] = 30;
    upperBound[4] = 30;
    upperBound[5] = 30;
    
    lowerBound[0]= -0.4;
    lowerBound[1]= -0.4;
    lowerBound[2]= -0.4;
    lowerBound[3]= -30;
    lowerBound[4]= -30;
    lowerBound[5]= -30;
    
    optimizer->SetBoundSelection(boundSelect);
    optimizer->SetUpperBound(upperBound);
    optimizer->SetLowerBound(lowerBound);
    
    optimizer->SetCostFunctionConvergenceFactor(1e+5);
    //optimizer->SetProjectedGradientTolerance(1.0);
    optimizer->SetMaximumNumberOfIterations(200);
    optimizer->SetMaximumNumberOfEvaluations(200);
    optimizer->SetMaximumNumberOfCorrections(5);
    

    
    CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
    optimizer->AddObserver(itk::IterationEvent(), observer);
    
    cout<<"registration"<<endl;
    
    try {
        registration->Update();
        cout<<"Optimizer stop condition : "<<registration->GetOptimizer()->GetStopConditionDescription()<<endl;

    } catch (itk::ExceptionObject &e) {
        cerr<<"erreur dans le recalage"<<endl;
        cerr<<e<<endl;
        return EXIT_FAILURE;
    }
    
    
    ParametersType finalParameters = registration->GetLastTransformParameters();
    
    double bestValue = optimizer->GetValue();
    
    cout<<"Results : "<<endl;
    cout<<"parametres : "<<finalParameters<<endl;
    cout<<"metric value : "<<bestValue<<endl;
    
    //transformation de l'image
    EulerTransformType::Pointer finalTsf = EulerTransformType::New();
    finalTsf->SetParameters(finalParameters);
    
    cout<<"Write results"<<endl;
//
//    ResampleFilterType::Pointer resampler1 = ResampleFilterType::New();
//    resampler1->SetInput(US_translated);
//    resampler1->SetTransform(finalTsf);
//    resampler1->SetSize(US_translated->GetLargestPossibleRegion().GetSize());
//    resampler1->SetOutputSpacing(US_translated->GetSpacing());
//    resampler1->SetOutputDirection(finalTsf->GetInverseMatrix()*US_translated->GetDirection());
//    resampler1->SetOutputOrigin(finalTsf->GetInverseTransform()->TransformPoint(US_translated->GetOrigin()));
//    
//    ImageType::Pointer finalVesselUS = resampler1->GetOutput();
//    
//    
//    WriterType::Pointer writer6 =  WriterType::New();
//    writer6->SetImageIO(n_io);
//    string out6 = outputPath+"/registeredVessel_US.nii.gz";
//    writer6->SetFileName(out6);
//    writer6->SetInput(finalVesselUS);
//    try {
//        writer6->Update();
//    } catch (itk::ExceptionObject &e) {
//        cerr<<"error whilte writing registered image"<<endl;
//        cerr<<e<<endl;
//        return EXIT_FAILURE;
//    }
    
    ResampleFilterType::Pointer resampler2 = ResampleFilterType::New();
    resampler2->SetInput(US_translated);
    resampler2->SetTransform(finalTsf);
    resampler2->SetSize(US_translated->GetLargestPossibleRegion().GetSize());
    resampler2->SetOutputSpacing(US_translated->GetSpacing());
    resampler2->SetOutputDirection(finalTsf->GetInverseMatrix()*US_translated->GetDirection());
    resampler2->SetOutputOrigin(finalTsf->GetInverseTransform()->TransformPoint(US_translated->GetOrigin()));
    
    ImageType::Pointer finalUS = resampler2->GetOutput();
    
    
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
//    typedef itk::AmoebaOptimizer OptimizerType;
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
//    typedef itk::MultiResolutionImageRegistrationMethod<BinaryImageType, BianryImageType> MultiResRegistrationType;
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
//    registration->SetMovingImage();
//    registration->SetFixedImageRegion(vesselImage->GetLargestPossibleRegion());
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
//    
//    OptimizerType::ParametersType simplexDelta(transform->GetNumberOfParameters());
//    simplexDelta[0] =0.3;
//    simplexDelta[1] =0.3;
//    simplexDelta[2] =0.3;
//    simplexDelta[3] = 20;
//    simplexDelta[4] = 20;
//    simplexDelta[5] = 20;
//    
//    optimizer->AutomaticInitialSimplexOff();
//    optimizer->SetInitialSimplexDelta(simplexDelta);
//    
//    //tolerance
//    optimizer->SetParametersConvergenceTolerance(1);
//    optimizer->SetFunctionConvergenceTolerance(1);
//    
//    optimizer->SetMaximumNumberOfIterations(200);
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
//    ResampleFilterType::Pointer resampler1 = ResampleFilterType::New();
//    resampler1->SetInput(vesselImageUS);
//    resampler1->SetTransform(finalTsf);
//    resampler1->SetSize(vesselImageUS->GetLargestPossibleRegion().GetSize());
//    resampler1->SetOutputSpacing(vesselImageUS->GetSpacing());
//    resampler1->SetOutputDirection(finalTsf->GetInverseMatrix()*vesselImageUS->GetDirection());
//    resampler1->SetOutputOrigin(finalTsf->GetInverseTransform()->TransformPoint(vesselImageUS->GetOrigin()));
//    
//    ImageType::Pointer finalVesselUS = resampler1->GetOutput();
//    
//    
//    WriterType::Pointer writer6 =  WriterType::New();
//    writer6->SetImageIO(n_io);
//    string out6 = outputPath+"/registeredVessel_US.nii.gz";
//    writer6->SetFileName(out6);
//    writer6->SetInput(finalVesselUS);
//    try {
//        writer6->Update();
//    } catch (itk::ExceptionObject &e) {
//        cerr<<"error whilte writing registered image"<<endl;
//        cerr<<e<<endl;
//        return EXIT_FAILURE;
//    }
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
//    string out7 = outputPath+"/registered_US.nii.gz";
//    writer7->SetFileName(out7);
//    writer7->SetInput(finalUS);
//    try {
//        writer7->Update();
//    } catch (itk::ExceptionObject &e) {
//        cerr<<"error whilte writing registered image"<<endl;
//        cerr<<e<<endl;
//        return EXIT_FAILURE;
//    }


    tend = std::time(NULL);
    texec = std::difftime(tend,tbegin);
    
    cout<<"temps en s : "<<texec<<endl;
    
    return 0;
}


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

