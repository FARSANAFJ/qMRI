#include "pch.h"
#include "qMRI.h"
#include<omp.h>

/**
 * Default constructor.
*
* @param       Nil
* @return      Nil
* @exception   Nil
* @see         Nil
* @since       1.0
*/
QuantitativeMapAPI::QuantitativeMapAPI()
{
    dicomImageIO = itk::GDCMImageIO::New();
    voldicomImageIO = itk::GDCMImageIO::New();
    reader = ReaderType::New();
    volreader = ReaderType::New();
    dicomImage = ImageType::New();
    voldicomImage = ImageType::New();


    volR2Star = ImageType::New();
    volT2Star = ImageType::New();
    volRho = ImageType::New();

    R2Star = ImageType::New();
    T2Star = ImageType::New();
    Rho = ImageType::New();
    calculatorR2Star = CalculatorType::New();
    calculatorT2Star = CalculatorType::New();
    calculatorRho = CalculatorType::New();
    StatisticsFilterR2Star = StatisticsFilterType::New();
    StatisticsFilterT2Star = StatisticsFilterType::New();
    StatisticsFilterRho = StatisticsFilterType::New();
}


/**
 * Default destructor.
 *
 * @param       Nil
 * @return      Nil
 * @exception   Nil
 * @see         Nil
 * @since       1.0
 */
QuantitativeMapAPI::~QuantitativeMapAPI()
{

}


bool QuantitativeMapAPI::Initialize(CString csDICOMDir)
{
    CStringA csDICOMDirA(csDICOMDir);
    if (DicomDataReading(csDICOMDirA))
    {
        m_strDICOMDir = csDICOMDirA;
      getvolMR();
      Parametricmap_Preprocessing();
      getvolT2Star();
        return true;
    }
    //T/*RACE_ERROR(ParametricMap, L"ParaMetricAPI Failed to initialze!");*/
    return false;
}


// Function to fit a polynomial using the least squares method
std::vector<double> fitPolynomial(const std::vector<double>& x, const std::vector<double>& y, int degree) {
    int n = x.size();
    int m = degree + 1;

    std::vector<double> result(m, 0.0);
    std::vector<std::vector<double>> A(m, std::vector<double>(m, 0.0));
    std::vector<double> B(m, 0.0);

    // Build the coefficient matrix A and the right-hand side vector B
    for (int i = 0; i < n; i++) {
        double xi = x[i];
        double yi = y[i];

        for (int j = 0; j < m; j++) {
            for (int k = 0; k < m; k++) {
                A[j][k] += std::pow(xi, j + k);
            }
            B[j] += std::pow(xi, j) * yi;
        }
    }

    // Solve the linear system
    for (int i = 0; i < m; i++) {
        for (int j = i + 1; j < m; j++) {
            double factor = A[j][i] / A[i][i];
            for (int k = i; k < m; k++) {
                A[j][k] -= factor * A[i][k];
            }
            B[j] -= factor * B[i];
        }
    }

    for (int i = m - 1; i >= 0; i--) {
        result[i] = B[i];
        for (int j = i + 1; j < m; j++) {
            result[i] -= A[i][j] * result[j];
        }
        result[i] /= A[i][i];
    }

    return result;
}

bool QuantitativeMapAPI::DicomDataReading(const char* dicomdirectory)
{
    reader->SetImageIO(dicomImageIO);

    // Get a list of DICOM series file names
    itk::GDCMSeriesFileNames::Pointer fileNames = itk::GDCMSeriesFileNames::New();
    fileNames->SetUseSeriesDetails(true);
    fileNames->SetDirectory(dicomdirectory);

    const itk::GDCMSeriesFileNames::SeriesUIDContainerType& seriesUIDs = fileNames->GetSeriesUIDs();
    if (seriesUIDs.empty()) {
        std::cerr << "No DICOM series found in the specified directory." << std::endl;
        return false;
    }

    // Iterate through DICOM series and process each series
    for (const auto& seriesUID : seriesUIDs) {
        const std::vector<std::string>& files = fileNames->GetFileNames(seriesUID);

        // Set the sorted file names for the reader
        reader->SetFileNames(files);

        try {
            reader->Update();
        }
        catch (itk::ExceptionObject& ex) {
            std::cerr << "Error reading DICOM series: " << ex << std::endl;
            CString csErr(ex.GetDescription());
            /* TRACE_ERROR(ParametricMap, csErr);*/
            return false;
        }
    }

    dicomImage = reader->GetOutput();

    return true;
}



int QuantitativeMapAPI::EchoTrainLength(const CString& dicomdirectory) {

    CString searchpath = dicomdirectory + _T("\\*.*");
    CFileFind finder;
    BOOL bWorking = finder.FindFile(searchpath); 
    std::string tagEchoNum;
    std::string tagEchoTime;
    std::string Instack;

    while (bWorking) {
        bWorking = finder.FindNextFile();
        if (!finder.IsDirectory())
        {

            CString filePath = finder.GetFilePath();

            dicomImageIO->GetValueFromTag("0018|0091", tagEchoNum);
            EchoNum = std::stod(tagEchoNum);

            dicomImageIO->SetFileName(CT2A(filePath));
            dicomImageIO->ReadImageInformation();
            dicomImageIO->GetValueFromTag("0020|9057", Instack);
            int InstackNum = std::stod(Instack);
            if (InstackNum == 57)
            {

                dicomImageIO->ReadImageInformation();
                dicomImageIO->GetValueFromTag("0018|0081", tagEchoTime);
                EchoTime = std::stoi(tagEchoTime);
                EchoTimes.push_back(EchoTime);

            }

        }
    }

    return EchoNum;

}

ImageType::Pointer QuantitativeMapAPI::getvolMR()
{
   
    return dicomImage;
}


ImageType::Pointer QuantitativeMapAPI::getvolR2Star()
{


    return volR2Star;

}

ImageType::Pointer QuantitativeMapAPI::getvolT2Star()
{

    return T2Star;

}

ImageType::Pointer QuantitativeMapAPI::getvolRho()
{

    return volRho;

}



void QuantitativeMapAPI::Parametricmap_Preprocessing() {

    dicomImage = reader->GetOutput();
    itk::Size<Dimension>Size{};
    int EchoNum = 5;
    Size[0] = dicomImage->GetLargestPossibleRegion().GetSize(0);
    Size[1] = dicomImage->GetLargestPossibleRegion().GetSize(1);
    Size[2] = dicomImage->GetLargestPossibleRegion().GetSize(2) / EchoNum;
    int N = Size[2];

    std::vector<PixelType> Xvalues;
    Xvalues.resize(EchoNum);

    volOriginalSpacing = dicomImage->GetSpacing();
    volOriginalSpacing[2] *= EchoNum;
    volOriginalDirection = dicomImage->GetDirection();
    volRegion.SetSize(Size);
    R2Star->SetRegions(volRegion);
    R2Star->SetSpacing(volOriginalSpacing);
    R2Star->SetDirection(volOriginalDirection);
    R2Star->Allocate();

    T2Star->SetRegions(volRegion);
    T2Star->SetSpacing(volOriginalSpacing);
    T2Star->SetDirection(volOriginalDirection);
    T2Star->Allocate();


    Rho->SetRegions(volRegion);
    Rho->SetSpacing(volOriginalSpacing);
    Rho->SetDirection(volOriginalDirection);
    Rho->Allocate();


    for (int n = 0; n < N; ++n)
    {
  #pragma omp parallel for
        for (int i = 0; i < Size[0]; ++i)
        {
 #pragma omp parallel for

            for (int j = 0; j < Size[1]; ++j)
            {

                std::fill(Xvalues.begin(), Xvalues.end(), 0);
                int cnt = 0;

                for (int k = 0 + n * EchoNum; k < EchoNum + n * EchoNum; ++k)
                {

                    double pixelvalue = dicomImage->GetPixel({ {i, j, k} });
                    double abspixelvalue = std::abs(pixelvalue);
                    double XlogValues = std::log(abspixelvalue);
                    Xvalues[cnt] = XlogValues;

                    cnt++;

                }
                EchoTimes = { 1,2,3,4,5 };
                std::vector<double> coefficients = fitPolynomial(Xvalues, EchoTimes, 1); // fitting
                std::vector<double> invcoefficients = fitPolynomial(EchoTimes, Xvalues, 1);
                double t2star = 1 / std::abs((invcoefficients[1])); // computing t2star

                if (std::isnan(t2star)) {
                    t2star = 0.0;
                }
                if (std::isinf(t2star)) {
                    t2star = 0.0;
                }
                T2Star->SetPixel({ {i, j, n} }, t2star);


                double r2star = std::abs((invcoefficients[1]));  // computing r2star
                if (std::isnan(r2star)) {
                    r2star = 0.0;
                }
                if (std::isinf(r2star)) {
                    r2star = 0.0;
                }
                R2Star->SetPixel({ {i, j, n} }, r2star);

                double rho = std::abs(invcoefficients[0]); // computing rho
                if (std::isnan(rho)) {
                    rho = 0.0;
                }
                if (std::isinf(rho)) {
                    rho = 0.0;
                }
                Rho->SetPixel({ {i, j, n} }, rho);




            }

        }
    }

    // Statistical Analysis T2Star

    calculatorT2Star->SetImage(T2Star);
    calculatorT2Star->Compute();


    t2starMin = calculatorT2Star->GetMinimum();
    t2starMax = calculatorT2Star->GetMaximum();

    StatisticsFilterT2Star->SetInput(T2Star);
    StatisticsFilterT2Star->Update();

    T2starStd = StatisticsFilterT2Star->GetSigma();
    T2starMean = StatisticsFilterT2Star->GetMean();

    // Statistical Analysis R2Star

    calculatorR2Star->SetImage(R2Star);
    calculatorR2Star->Compute();


    r2starMin = calculatorR2Star->GetMinimum();
    r2starMax = calculatorR2Star->GetMaximum();

    StatisticsFilterR2Star->SetInput(R2Star);
    StatisticsFilterR2Star->Update();

    R2starStd = StatisticsFilterR2Star->GetSigma();
    R2starMean = StatisticsFilterR2Star->GetMean();

    // Statistical Analysis Rho

    calculatorRho->SetImage(Rho);
    calculatorRho->Compute();


    rhoMin = calculatorRho->GetMinimum();
    rhoMax = calculatorRho->GetMaximum();

    StatisticsFilterRho->SetInput(Rho);
    StatisticsFilterRho->Update();

    RhoStd = StatisticsFilterRho->GetSigma();
    RhoMean = StatisticsFilterRho->GetMean();



}



void QuantitativeMapAPI::parametricmap_Computation()

// initializing

{


    dicomImage = reader->GetOutput();

    itk::Size<Dimension>Size{};
    Size[0] = dicomImage->GetLargestPossibleRegion().GetSize(0);
    Size[1] = dicomImage->GetLargestPossibleRegion().GetSize(1);
    Size[2] = dicomImage->GetLargestPossibleRegion().GetSize(2) / EchoNum;

    volOriginalSpacing = dicomImage->GetSpacing();
    volOriginalSpacing[2] *= EchoNum;
    volOriginalDirection = dicomImage->GetDirection();
    volRegion.SetSize(Size);
    volR2Star->SetRegions(volRegion);
    volR2Star->SetSpacing(volOriginalSpacing);
    volR2Star->SetDirection(volOriginalDirection);
    volR2Star->Allocate();

    volT2Star->SetRegions(volRegion);
    volT2Star->SetSpacing(volOriginalSpacing);
    volT2Star->SetDirection(volOriginalDirection);
    volT2Star->Allocate();


    volRho->SetRegions(volRegion);
    volRho->SetSpacing(volOriginalSpacing);
    volRho->SetDirection(volOriginalDirection);
    volRho->Allocate();

    // preprocessing



    // preprocessing



    std::vector<PixelType> Xvalues;
    Xvalues.resize(EchoNum);


    for (int n = 0; n < Size[2]; ++n)
    {

        for (int i = 0; i < Size[0]; ++i)
        {

            for (int j = 0; j < Size[1]; ++j)
            {
                std::fill(Xvalues.begin(), Xvalues.end(), 0);
                int cnt = 0;


                for (int k = 0 + n * EchoNum; k < EchoNum + n * EchoNum; ++k)
                {

                    double pixelvalue = dicomImage->GetPixel({ {i, j, k} });
                    double abspixelvalue = std::abs(pixelvalue);
                    double XlogValues = std::log(abspixelvalue);
                    Xvalues[cnt] = XlogValues;

                    cnt++;

                }

                std::vector<double> coefficients = fitPolynomial(Xvalues, EchoTimes, 1); // fitting
                std::vector<double> invcoefficients = fitPolynomial(EchoTimes, Xvalues, 1);


                // R2Star computation
                double r2star = std::abs(invcoefficients[1]);
                double newMaxR2Star = r2starMax + (R2starMean + 2 * R2starStd);
                double newMinR2Star = r2starMin + (R2starMean - 2 * R2starStd);
                double scaled_r2star = ((r2star - r2starMin) * (newMaxR2Star - newMinR2Star) / (r2starMax - r2starMin)) + newMinR2Star;
            }
        }
    }
}
