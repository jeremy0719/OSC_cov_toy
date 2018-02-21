///////////////////////////////////////////////////////////////////////
/// Author: P T Surukuchi
/// Date: 2017
/// TODO: Try to avoid using array, instead use a vector or TArrayD
/// Implement better exception handling
/// Implement private, public, protected functions correctly
/// Implement functions that takes a TFile and update it
////////////////////////////////////////////////////////////////////////

#ifndef TCOVARIANCEMATRIXGENERATOR_HH
#define TCOVARIANCEMATRIXGENERATOR_HH
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
using namespace std;

#include "TPrincipal.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TMatrixD.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TKey.h"
#include "TString.h"


#include "TCovUtilities.hh"


/// The number of toy histograms to draw
const int nHists=20;

/// \brief Parent class inherited by other classes to generate covariance matrix using all the input vectors/arrays
/**  More description */
class TCovMatGenerator
{
public:
    //TCovMatGenerator(){};
    /// Default Constructor, takes number of bins
    TCovMatGenerator(int);
    
    /// Constructor which takes a TString for the output file, takes number of bins
    TCovMatGenerator(int, TString);

    /// Constructor which takes a user-defined TFile for the output, takes number of bins
    TCovMatGenerator(int, TFile&);
    
    /// Getter for the number of bins #nBins.
    /** Returns the number of bins corresponding to the size of covariance matrix.
     It is a constant value defined in the TCovMatGenerator()*/
    int matSize(){return nBins;}
    
    /// Getter for the number of samples #nSamples.
    /** Returns number of samples/toys used. It is a changing number and each call to addSample() increases this value by one.*/
    int nMat(){return nSamples;}
    
    ///Function to increment the number of samples by one
    /** Returns 1 if successful or 0 */
    int appendNSample()
    {
        ++nSamples;
        return 1;
    }
    
    /// Add a sample to the the exisiting samples to be used to produce a mean covarince matrix.
    /** Uses the following formula:
     \f$ V=V+cV \f$
     Where V is already existing covariance matrix and cV is the covariance matrix generated using the argument accepted by the function. Check the function makeCovMat. */
    virtual int addSample(const vector<double>&)
    {
        std::cout<< "Inherit me!";
        return 0;
    }
    
    /// Initialize the covariance vector with the size defined by the input integer
    /** Returns: integer value, 1 for pass and 0 for fail. */
    virtual int initialize()
    {
        std::cout<< "Inherit me!";
        return 0;
    }
    /// TODO, Forgot what this is supposed to do
    virtual int finalize()
    {
        std::cout<< "Inherit me!";
        return 0;
    }
    
    /// Destructor
    //  May not be needed afterall
    virtual ~TCovMatGenerator()
    {
        delete [] energyBins;
    }
    
protected:
    
    /// Number of input samples/toys and is incremented everytime addSample() is successful
    int nSamples=0;
    /// Number of bins corresponding to the size of covariance matrix
    int nBins;
    
    /// Final covariance matrix of all the samples
    const TMatrixD* meanCovMat;
    ///  Final Correlation matrix of all samples
    TMatrixD* meanCorrMat;
    
    /// #nHists Histograms corresponsing to samples/toys
    //TH1D *hSample[nHists];
    /// #nHists Histograms corresponsing to fractional differences in samples/toys
    //TH1D *hFracSample[nHists];
    /// Histogram of mean fractional value
    TH1D *hMeanFracSample;
    /// Histograms for first #nHists covariance matrix
    //TH2D *hCovMat[nHists];
    /// Histogram of mean covariance matrix
    TH2D *hCovMatrix;
    /// Histogram of mean covariance matrix
    TH2D *hCorrMatrix;
    /// Histogram of reduced (statistics free) covariance matrix
    TH2D *hRedCovMatrix;
    /// Output ROOT file with toy, fractional variation, covariance and correlation histograms
    TFile* outputFile;
    /// sigma values
    TH1 *hSigmaValues;
    /// Array used for variable binning of energy used. TODO: Generalize this functionality, perhaps use an argument for one of the functions to do this.
    double* energyBins;
    
private:
};


/// Generate a covariance matrix using a reference vector
/**  Uses the folowing formula for generation of the matrix.
 \f$ \Sigma_{ij}= V_ {ij} = \frac{1}{M} \sum\limits^{toys}(F^{obs,toy}_{i} - F^{ref}_{i})(F^{obs,toy}_{j} - F^{ref}_{j}) \f$
 Where \f$ V_ {ij} \f$ is covariance matrix, \f$F^{ref}\f$ is the reference vector assigned in the initiailize() function and \f$F^{obs,toy}\f$ is the toy matrix.*/
class RefCovMatGenerator: public TCovMatGenerator
{
public:
    /// Default Constructor, takes the reference vector to be used for covariance matrix generation
    RefCovMatGenerator(const vector<double>&);
    
    /// Constructor, takes the reference vector to be used for covariance matrix generation and TString for root output file
    RefCovMatGenerator(const vector<double>&, TString);
    
    /// Constructor, takes the reference vector to be used for covariance matrix generation and a TFile that has already been defined by the user
    RefCovMatGenerator(const vector<double>&,TFile&);
    
    /// Add a sample to the the exisiting samples to be used to produce a mean covarince matrix.
    /** Each successful call to this function appends 1 to #nSamples. In addition, a covariance matrix is generated for each toy/sample corresponding to the inout vector w.r.t to the #refVector. For #nSamples < 20, the spectrum, fractional spectrum are drawn. Will be saved along with the other other histograms when finalize() is called.  */
    int addSample(const vector<double>&);
    
    /// A desctructor that calls finalize
    virtual ~RefCovMatGenerator()
    {
        finalize();
    }
    
private:
    
    /// Reference vector used for calculating the covariances
    vector<double> refVector;
    
    /*
    /// Covariance matrix corresponding to each of the vector input in addSample()
    TMatrixD* toyCovMat;
    
    /// aggregate of covariance matrices for taking sum of all #toyCovMat, perhaps since it is temporary can be moved to the
    TMatrixD* aggCovMat;*/
    
    /// Toy covariance matrices for the printing
    std::vector<TH2D*> hToyCovMatrix;
    std::vector<TH2D*> hRedToyCovMatrix;
    
    /// Initialize the covariance vector with the size defined by the input integer
    /** Returns: integer value, 1 for pass and 0 for fail. */
    int initialize()
    {
        std::cerr << "Oops! Wrong function chosen"<<std::endl;
        return 0;
    }
    
    /// Overloaded function which takes reference to a vector as input
    /** Returns: integer value, 1 for pass and 0 for fail. */
    int initialize(const vector<double>&);
    
    /// Calculate and plot the #meanCovMat, draw the mean fractional diff vector and other histograms.
    int finalize();
    
};


/// Generate a covariance matrix using mean of all the input vectors
/**  Uses the folowing formula for generation of the matrix.
 \f$ \Sigma_{ij}= V_ {ij} = \frac{1}{M} \sum\limits^{toys}(F^{obs,toy}_{i} - F^{mean}_{i})(F^{obs,toy}_{j} - F^{mean}_{j}) \f$
 Where \f$ V_ {ij} \f$ is covariance matrix, \f$F^{mean} =\sum\limits^{toys} F^{obs,toy} \f$ is the reference vector assigned in the initiailize() function and \f$F^{obs,toy}\f$ is the toy matrix.*/
class MeanTCovMatGenerator: public TCovMatGenerator
{
public:
    /// Default Constructor, takes number of bins
    MeanTCovMatGenerator(int);
    
    /// Constructor, takes number of bins and TString for output ROOT file
    MeanTCovMatGenerator(int, TString);
    
    /// Add a sample to the the exisiting samples to be used to produce a mean covarince matrix.
    /** Uses the following formula:
     \f$ V=V+cV \f$
     Where V is already existing covariance matrix and cV is the covariance matrix generated using the argument accepted by the function. Check the function makeCovMat. Each step along the way generates a histogram corresponsing to the covariance matrix */
    int addSample(const vector<double>&);
    
    virtual ~MeanTCovMatGenerator()
    {
        delete [] binContent;
        finalize();
    }
    
private:
    /// Array used for the value of the content of the bin
    double* binContent;
    
    /// TPrincipal object for direct calculation of the covariance matrix using Principal Components Analysis <A HREF="https://root.cern.ch/doc/master/classTPrincipal.html"> package by ROOT.
    TPrincipal* principal;
    
    /// Initialize the covariance vector with the size defined by the input integer
    /** Returns: integer value, 1 for pass and 0 for fail. */
    int initialize();
    
    /// Make the final mean covariance matrix and generate all the histograms
    int finalize();
    
};

#endif
