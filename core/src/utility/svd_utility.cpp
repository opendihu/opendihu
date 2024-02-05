#include "utility/svd_utility.h"

#ifdef HAVE_LAPACK
#include <cblas.h>
#endif

#include "easylogging++.h"

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <sstream>

using namespace std;

// takes real matrix input (rows x cols) as double[] in column major order
// performs singular-value decomposition on input utilizing LAPACKE_dgesvd
// stores the left-singular vectors (column-wise leftSingularVectors)
void SvdUtility::getSVD(double input[], int rows, int cols,
                        double leftSingVec[]) {
#ifdef HAVE_LAPACK
  int min = std::min(cols, rows);
  double *singVal = new double[min];
  double *superb = new double[min];
  double *rightSingVecT = new double[min * cols];

  double *inputCopy = new double[rows * cols];
  for (int col = 0; col < cols; ++col) {
    for (int row = 0; row < rows; ++row) {
      inputCopy[col * rows + row] = input[col * rows + row];
    }
  }

  LAPACKE_dgesvd(LAPACK_COL_MAJOR, 's', 's', rows, cols, inputCopy, rows,
                 singVal, leftSingVec, rows, rightSingVecT, min, superb);
#else
  LOG(FATAL) << "Not compiled with LAPACK, getSVD is not possible.";
#endif
}

// takes real matrix input (rows x cols) as double[] in column major order
// performs singular-value decomposition on input utilizing LAPACKE_dgesvd
// stores the left-singular vectors (column-wise leftSingularVectors), singular
// values (singularValues as vector, diagonal entries of sigma as matrix) and
// right-singular vectors (row-wise rightSingularVectorsTransposed) as double
// arrays
void SvdUtility::getSVD(double input[], int rows, int cols,
                        double leftSingVec[], double sigma[],
                        double rightSingVecT[]) {
#ifdef HAVE_LAPACK
  int min = std::min(cols, rows);
  double *singVal = new double[min];
  double *superb = new double[min];

  printMatrix("snapshots", input, rows, cols);

  double *inputCopy = new double[rows * cols];
  for (int col = 0; col < cols; ++col) {
    for (int row = 0; row < rows; ++row) {
      inputCopy[col * rows + row] = input[col * rows + row];
    }
  }

  // if error or wrong numbers by 'a' use 's' instead
  LAPACKE_dgesvd(LAPACK_COL_MAJOR, 's', 's', rows, cols, inputCopy, rows,
                 singVal, leftSingVec, rows, rightSingVecT, min, superb);

  // build diagonal matrix sigma from vector singularValues
  for (int row = 0; row < min; ++row) {
    for (int col = 0; col < min; ++col) {
      if (row == col) {
        sigma[row + col * min] = singVal[row];
      } else {
        sigma[row + col * min] = 0;
      }
    }
  }

  printMatrix("leftSingVec", leftSingVec, rows, min);
  printMatrix("sigma", sigma, min, min);
  printMatrix("rightSingVecT", rightSingVecT, min, cols);
#else
  LOG(FATAL) << "Not compiled with LAPACK, getSVD is not possible.";
#endif
}

// takes complex matrix input (rows x cols) as double _Complex[] in column major
// order performs singular-value decomposition on input utilizing LAPACKE_zgesvd
// stores the left-singular vectors (column-wise leftSingularVectors), singular
// values (singularValues as vector, diagonal entries of sigma as matrix) and
// right-singular vectors (row-wise rightSingularVectorsTransposed) as double
// _Complex arrays
void SvdUtility::getSVD(double _Complex input[], int rows, int cols,
                        double _Complex leftSingVec[], double sigma[],
                        double _Complex rightSingVecT[]) {
#ifdef HAVE_LAPACK
  int min = std::min(cols, rows);
  double *singVal = new double[min];
  double *superb = new double[min];

  double _Complex *inputCopy = new double _Complex[rows * cols];
  for (int col = 0; col < cols; ++col) {
    for (int row = 0; row < rows; ++row) {
      inputCopy[col * rows + row] = input[col * rows + row];
    }
  }

  LAPACKE_zgesvd(LAPACK_COL_MAJOR, 's', 's', rows, cols, inputCopy, rows,
                 singVal, leftSingVec, rows, rightSingVecT, min, superb);

  // cout << "info: " << info << endl << endl;

  // build diagonal matrix sigma from vector singularValues
  for (int row = 0; row < min; ++row) {
    for (int col = 0; col < min; ++col) {
      if (row == col) {
        sigma[row + col * min] = singVal[row];
      } else {
        sigma[row + col * min] = 0;
      }
    }
  }
#else
  LOG(FATAL) << "Not compiled with LAPACK, getSVD is not possible.";
#endif
}

void SvdUtility::reconstructSnapshots(int rows, int cols, double leftSingVec[],
                                      double sigma[], double rightSingVecT[],
                                      double output[]) {
#ifdef HAVE_LAPACK

  int rank = std::min(rows, cols);
  double *C = new double[rank * cols];

  // Sigma multiplied by transpose of the right singular vector
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, rank, cols, rank, 1.0,
              sigma, rank, rightSingVecT, rank, 0.0, C, rank);
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, rows, cols, rank, 1.0,
              leftSingVec, rows, C, rank, 0.0, output, rows);

  cout << "Reconstructed matrix:" << endl;
  for (int row = 0; row < rows; ++row) {
    for (int col = 0; col < cols; ++col)
      cout << output[col * rows + row] << ' ';
    cout << endl;
  }

  printMatrix("v_reconst", output, rows, cols);
#else
  LOG(FATAL)
      << "Not compiled with LAPACK, reconstructSnapshots is not possible.";
#endif
}

// takes real matrix input (rows x cols) as double[]
// prints name and input entry-wise
void SvdUtility::printMatrix(std::string name, double input[], int rows,
                             int cols) {
  cout << "    " << name << " =" << endl
       << "    \\begin{bmatrix*}[r]" << endl
       << "        " << input[0];

  for (int row = 0; row < rows; ++row) {
    for (int col = 1; col < cols; ++col) {
      cout << " & " << input[col * rows + col]; //???? col * rows + row
    }

    if (row < rows - 1) {
      cout << " \\\\" << endl << "        " << input[row + 1];
    } else {
      cout << endl << "    \\end{bmatrix*}" << endl;
    }
  }

  cout << endl;
}

// reads CSV cell by cell as vector
std::vector<double> SvdUtility::readCSV(string filename) {
  ifstream data(filename);
  string line;
  vector<double> parsedCsv;
  int i = 0;
  int j = 0;
  while (getline(data, line)) {
    stringstream lineStream(line);
    i++;
    string cell;
    j = 0;
    while (getline(lineStream, cell, ',')) {
      parsedCsv.push_back(stof(cell));
      j++;
    }
    // cout << i << endl;
    // cout << j << endl;
  }
  return parsedCsv;
}

// reads specified rows of CSV cell by cell as vector
std::vector<double> SvdUtility::readCSV(string filename, int rows) {
  ifstream data(filename);
  string line;
  vector<double> parsedCsv;
  for (int i = 0; i < rows; i++) {
    getline(data, line);
    stringstream lineStream(line);
    string cell;
    while (getline(lineStream, cell, ',')) {
      parsedCsv.push_back(stof(cell));
    }
  }
  return parsedCsv;
}

// writes vector cell by cell as CSV
void SvdUtility::writeCSV(string filename, std::vector<double> values, int m,
                          int n) {
  ofstream data;
  data.open(filename);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      data << std::to_string(values[i * n + j]) << ",";
      cout << std::to_string(values[i * n + j]) << ",";
    }
    data << "\n";
    cout << endl;
  }
  data.close();
}

// writes vector cell by cell as CSV
void SvdUtility::writeCSV(string filename, double values[], int m, int n,
                          bool columnWise) {
  ofstream data;
  data.open(filename, ios_base::app);

  if (!columnWise)
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        data << std::to_string(values[i * n + j]) << ",";
      }
      data << "\n";
    }
  else
    for (int j = 0; j < n; j++) {
      for (int i = 0; i < m; i++) {
        data << std::to_string(values[j * m + i]) << ",";
      }
      data << "\n";
    }

  data.close();
}

// returns number of rows of CSV
int SvdUtility::getCSVRowCount(string filename) {
  ifstream data(filename);
  string line;
  int i = 0;
  while (getline(data, line)) {
    stringstream lineStream(line);
    i++;
  }
  return i;
}

// returns number of columns of CSV
int SvdUtility::getCSVColumnCount(string filename) {
  ifstream data(filename);
  string line;
  vector<double> parsedCsv;
  int j = 0;
  if (getline(data, line)) {
    stringstream lineStream(line);
    string cell;
    while (getline(lineStream, cell, ',')) {
      parsedCsv.push_back(stof(cell));
      j++;
    }
  }
  return j;
}
