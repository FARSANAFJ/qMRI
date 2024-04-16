#pragma once

#include "afxwin.h"
#include <afx.h>
#include <vector>
#include <map>
#include <cmath>
#include <iostream>
#include <sstream>
#include <algorithm>
#include "itkImageSeriesReader.h"
#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkStatisticsImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include<omp.h>



typedef double PixelType;
constexpr unsigned int Dimension = 3;
typedef itk::Image<PixelType, Dimension> ImageType;
typedef itk::ImageSeriesReader<ImageType> ReaderType;
typedef itk::Image<PixelType, 2> ImageTypeR2Star;
typedef itk::Image<PixelType, 2> ImageTypeRho;
typedef itk::Image<PixelType, 2> ImageTypeT2Star;

using CalculatorType = itk::MinimumMaximumImageCalculator<ImageType>;
using StatisticsFilterType = itk::StatisticsImageFilter<ImageType>;











class  QuantitativeMapAPI
{
public:
    QuantitativeMapAPI();
    ~QuantitativeMapAPI();



public:

    bool Initialize(const CString csDICOMDir);
   /* CStringArray& GetAttributeNameList();*/

   /* std::string GetModality();*/
    bool DicomDataReading(const char* path);
    int EchoTrainLength(const CString& dicomdirectory);

    /*void ExtractModality();*/

    void parametricmap_Computation();
    void Parametricmap_Preprocessing();

    
    
    ImageType::Pointer getvolR2Star();
    ImageType::Pointer getvolRho();
    ImageType::Pointer getvolT2Star();
    ImageType::Pointer getvolMR();

    


public:
    typedef double PixelType;

    ReaderType::Pointer reader;
    ReaderType::Pointer mreader;
    
    ImageType::Pointer volMR;
    ImageType::Pointer dicomImage;
    itk::GDCMImageIO::Pointer dicomImageIO;
    ImageType::RegionType OriginalRegion;
    ImageType::SpacingType OriginalSpacing;
    ImageType::PointType OriginalOrigin;
    ImageType::DirectionType OriginalDirection;
    ImageType::SizeType size;

    ReaderType::Pointer volreader;
    ImageType::Pointer voldicomImage;
    ImageType::RegionType volRegion;
    ImageType::SpacingType volOriginalSpacing;
    ImageType::PointType volOriginalOrigin;
    ImageType::DirectionType volOriginalDirection;
    ImageType::SizeType volsize;
    itk::GDCMImageIO::Pointer voldicomImageIO;
    std::string EchoNumber;
    int Echovalue;


private:

    ImageType::Pointer volR2Star;
    ImageType::Pointer volT2Star;
    ImageType::Pointer volRho;
  

    ImageType::Pointer R2Star;
    ImageType::Pointer T2Star;
    ImageType::Pointer Rho;
    


    CalculatorType::Pointer calculatorR2Star;
    CalculatorType::Pointer calculatorT2Star;
    CalculatorType::Pointer calculatorRho;
    StatisticsFilterType::Pointer StatisticsFilterR2Star;
    StatisticsFilterType::Pointer StatisticsFilterT2Star;
    StatisticsFilterType::Pointer StatisticsFilterRho;




    ImageTypeR2Star::RegionType Region;




    int EchoNum;
    std::vector<double>EchoTimes;
    double EchoTime;

    double t2starMin = std::numeric_limits<double>::infinity();
    double t2starMax = -std::numeric_limits<double>::infinity();
    double T2starStd;
    double T2starMean;


    double r2starMin = std::numeric_limits<double>::infinity();
    double r2starMax = -std::numeric_limits<double>::infinity();


    double R2starStd;
    double R2starMean;

    double rhoMin = std::numeric_limits<double>::infinity();
    double rhoMax = -std::numeric_limits<double>::infinity();
    double RhoStd;
    double RhoMean;

    std::string m_strModality;
    std::string m_strDICOMDir;
    CStringArray m_csAttributeNameList;


private:
   

    QuantitativeMapAPI* m_QuantitativeMapAPI;
};


