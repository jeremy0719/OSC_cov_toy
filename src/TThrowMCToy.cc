#include "TThrowMCToy.hh"

// Default Class Constructor
TThrowMCToy::TThrowMCToy() {
}


// Class Constructor
TThrowMCToy::TThrowMCToy(TVectorD &v_inputParameters, TMatrixDSym &m_inputCovariance) {
  // Assign true status for decomposition
  status = true;
  
  // Number of parameters (or L-E bins)
  nParameters = v_inputParameters.GetNrows();

  // Set the central values parameters vector
  v_Parameters = new TVectorD(v_inputParameters);
  
  // Set the covariance matrix
  m_Covariance = new TMatrixDSym(m_inputCovariance);
  
  // Create a Cholesky matrix (an nxn matrix)
  m_CholeskyDecomp = new TMatrixD(nParameters, nParameters);
  // Fill Cholesky matrix with zeros
  m_CholeskyDecomp->Zero();
  
  // Calculate the Cholesky Decomposition
  CholeskyDecompose((*m_CholeskyDecomp));
  
}

// Class Deconstructor
TThrowMCToy::~TThrowMCToy() {
  if (v_Parameters != NULL) delete v_Parameters;
  if (m_Covariance != NULL) delete m_Covariance;
  if (m_CholeskyDecomp != NULL) delete m_CholeskyDecomp;
}

// Performs a Cholesky decomposition
// Specifically, it decomposes m_Covariance into m_CholDec
// There is a root function that will do this... but it may or may not have a bug
// At least that is the word on the street
void TThrowMCToy::CholeskyDecompose(TMatrixD &m_CholDec) {
  
  for (int i = 0; i < nParameters; i++) {
    for (int j = 0; j < i+1; j++) {
      // Diagonal Elements
      if (i == j) {
        double element = (*m_Covariance)(i,i);
        for (int k = 0; k < i; k++) {
          element -= TMath::Power(m_CholDec(i,k), 2);
        }
  
        // If the diagonal element is less than zero,
        // the decomposition will have imaginary numbers... bad!
        if (element < 0) {

	  std::cout << "Cholesky Decomposition failed!" << std::endl;
          // Set the decomposition status
          status = false;
          // Reset matrix
          m_CholDec.Zero();
          return;
        }
        
        m_CholDec(i,i) = TMath::Sqrt(element);
      } // Diagonals
      // Off-Diagonal Elements (only lower half)
      else if (j < i) {
        double element = (*m_Covariance)(i,j);
        for (int k = 0; k < j; k++) {
          element -= m_CholDec(i,k) * m_CholDec(j,k);
        }
        element /= m_CholDec(j,j);
        
        m_CholDec(i,j) = element;
        
      } // End of Off-Diagonals
      else m_CholDec(i,j) = 0.;
      
    } // End of loop through j values
  } // End of loop through i values
  
}

// A function to create a random pair of numbers
// The pair will be independant as well as
// created in a gaussian distribution with a mean of 0 and a sigma of 1
// Uses the Box-Muller polar form
void TThrowMCToy::RandomGaussian(double *z) {
  
  double v1 = 2.0 * randomGen.Rndm() - 1.0;
  double v2 = 2.0 * randomGen.Rndm() - 1.0;
  
  double r = v1*v1 + v2*v2;
  
  while (r == 0 || r >= 1.0) {
    v1 = 2.0 * randomGen.Rndm() - 1.0;
    v2 = 2.0 * randomGen.Rndm() - 1.0;
    r = v1*v1 + v2*v2;
  }
  
  z[0] = v1 * TMath::Sqrt(-2.0 * TMath::Log(r)/r);
  z[1] = v2 * TMath::Sqrt(-2.0 * TMath::Log(r)/r);
  
}

// This function takes in a vector that will contain the thrown parameters.
// It will calculate throws based on the appropriate correlated errors.
void TThrowMCToy::ThrowExperiment(std::vector<double> &thrownParameters) {
  // Ensure vector is empty and the right size
  thrownParameters.empty();
  thrownParameters.resize(nParameters);
  
  // Check that the process is possible, i.e. the covariance matrix must be decomposable
  if (!status) {
    std::cout << "Cholesky matrix is not decomposable." << std::endl;
    std::cout << "Process not possible, returning central values." << std::endl;
    for (int i = 0 ; i < nParameters; i++) {
      thrownParameters[i] = (*v_Parameters)(i);
    }
    return;
  }
  
  // A vector of the random gaussian numbers
  // mu = 0 and sigma = 1 for each one.
  TVectorD randParameters (nParameters);
  
  // Fill the TVectorD with values.
  for (int i = 0; i < nParameters; i++) {
    // Check that the next entry is still in bounds
    if (2*i >= nParameters) break;
    // Get two random gaussian numbers
    double z[2];
    RandomGaussian(z);
    randParameters(2*i) = z[0];
    
    // Check that the next entry is still in bounds
    if (2*i+1 >= nParameters) break;
    randParameters(2*i+1) = z[1];
  }
  
  // Find the correlated adjustments
  TVectorD adjParameters = (*m_CholeskyDecomp) * randParameters;

  // Calculate Thrown Parameters
  for (int i = 0; i < nParameters; i++) {
    thrownParameters[i] = (*v_Parameters)(i) + adjParameters(i);
  }
  
}
