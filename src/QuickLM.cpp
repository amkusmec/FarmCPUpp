#include <Rcpp.h>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppParallel)]]
#include <RcppEigen.h>
#include <RcppParallel.h>

#if !defined(EIGEN_USE_MKL) // Don't use R Lapack.h if MKL is enabled
#include <R_ext/Lapack.h>
#endif

// [[Rcpp::depends(BH, bigmemory)]]
#include <bigmemory/MatrixAccessor.hpp>
#include <numeric>

// [[Rcpp:plugins(cpp11)]]

using namespace RcppParallel;

using Eigen::ColPivHouseholderQR;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::RowVectorXd;
using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::Upper;

using Rcpp::_;
using Rcpp::as;
using Rcpp::IntegerVector;
using Rcpp::List;
using Rcpp::NumericMatrix;
using Rcpp::NumericVector;
using Rcpp::wrap;
using Rcpp::XPtr;

using std::invalid_argument;
using std::numeric_limits;

typedef MatrixXd::Index Index;
typedef ColPivHouseholderQR<MatrixXd>::PermutationType Permutation;

// Parallel worker structure
struct LM : public Worker {
  // Input variables
  const VectorXd y;
  const MatrixXd X;
  VectorXi observations;
  const int npcs;
  const int nqtn;
  const Map<MatrixXd> bM;
  const int df;

  // Output structures
  RVector<double> coefficients;
  RVector<double> stderr;
  RMatrix<double> seqQTN;

  // Constructor
  LM(const VectorXd y, const MatrixXd X, const XPtr<BigMatrix> bM,
     IntegerVector observations, const int npcs, const int nqtn,
     NumericVector coefficients, NumericVector stderr,
     const int df, NumericMatrix seqQTN)
    : y(y), X(X), bM(Map<MatrixXd>((double *)bM->matrix(), bM->nrow(), bM->ncol()  )),
      observations(as<VectorXi>(observations)), npcs(npcs), nqtn(nqtn), coefficients(coefficients),
      stderr(stderr), df(df), seqQTN(seqQTN) {}

  // For extracting elements of a vector by a vector of indices
  // Modified to only accept my predetermined inputs and output
  // From: http://stackoverflow.com/questions/26267194/how-to-extract-a-subvector-of-a-eigenvector-from-a-vector-of-indices-in-eige?rq=1
  inline VectorXd extract2(const VectorXd full, const VectorXi ind) {
    int num_indices = ind.innerSize();
    VectorXd target(num_indices);
    for (int i = 0; i < num_indices; i++)
    {
      target[i] = full[ind[i]];
    }
    return target;
  }

  // Operator
  void operator()(std::size_t begin, std::size_t end) {
    // Pre-compute the constant parts of the system
    MatrixXd XtX(X.transpose() * X);
    MatrixXd Xty(X.transpose() * y);
    MatrixXd yty(y.transpose() * y);
    ColPivHouseholderQR<MatrixXd> PQR(XtX);
    MatrixXd iXtX(PQR.inverse());

    // Reserve space for the inverse of the LHS
    MatrixXd iXX = MatrixXd::Zero(X.cols() + 1, X.cols() + 1);
    Index p = iXX.cols();

    // Main loop over all SNPs
    for (size_t i = begin; i < end; ++i) {
      // Construct the new model matrix
      // Single-value accession selects a column in column-major storage
      VectorXd snpCalls(bM.col(i));
      VectorXd x(extract2(snpCalls, observations));

      // Process the marker effects
      MatrixXd xty(x.transpose() * y);
      MatrixXd xtx(x.transpose() * x);
      MatrixXd Xtx(X.transpose() * x);

      // Compute the partitioned inverses
      MatrixXd B21(Xtx.transpose() * iXtX);
      MatrixXd t2(B21 * Xtx);
      double B22 = xtx(0, 0) - t2(0, 0);
      double iB22 = 1.0 / B22;
      MatrixXd NiB22B21(-iB22 * B21);
      MatrixXd iXX11(iXtX + iB22 * B21.transpose() * B21);

      // Fill in the inverse of the LHS
      iXX.block(0, 0, p - 1, p - 1) = iXX11;
      iXX.topRightCorner(NiB22B21.cols(), 1) = NiB22B21.row(0);
      iXX.bottomLeftCorner(1, NiB22B21.cols()) = NiB22B21.row(0);
      iXX(p - 1, p - 1) = iB22;

      // Compute coefficients and standard errors
      MatrixXd rhs(Xty.rows() + xty.rows(), Xty.cols());
      rhs << Xty, xty;
      VectorXd m_coef(iXX.transpose() * rhs);
      double s2 = (yty - m_coef.transpose() * rhs).norm() / double(df);
      VectorXd se((iXX.diagonal() * s2).cwiseSqrt());

      if (B22 < 10.0e-08) {
        coefficients[i] = ::NA_REAL;
        stderr[i] = ::NA_REAL;

        if (nqtn != 0) {
          VectorXd tvalues(VectorXd::Constant(nqtn, ::NA_REAL));
          for (int j = 0; j < nqtn; j++) {
            seqQTN(i, j) = tvalues(1 + npcs + j);
          }
        }
      } else {
        coefficients[i] = m_coef(p - 1);
        stderr[i] = se(p - 1);

        if (nqtn != 0) {
          RowVectorXd tvalues = m_coef.array() / se.array();
          for (int j = 0; j < nqtn; j++) {
            seqQTN(i, j) = tvalues(1 + npcs + j);
          }
        }
      }
    }
  }
};

//' Parallelized general linear model.
//'
//' @export
// [[Rcpp::export]]
Rcpp::List QuickLM(Rcpp::NumericVector ys, Rcpp::NumericMatrix Xs, SEXP pBigMat, IntegerVector observations, int npcs, int nqtn) {
  // Check input sizes
  if (ys.length() != Xs.nrow()) throw invalid_argument("Size mismatch.");

  // Check the big.matrix type
  XPtr<BigMatrix> xpMat(pBigMat);
  if (xpMat->matrix_type() != 8) throw Rcpp::exception("big.matrix is not numeric.");

  // Allocate output data structures
  size_t n_var = xpMat->ncol();
  NumericVector coefficients(n_var);
  NumericVector stderr(n_var);
  NumericMatrix seqQTN(n_var, nqtn);

  // Degrees of freedom
  int df = Xs.nrow() - Xs.ncol() - 1;

  // Pass input and output to the worker
  LM lm(as<VectorXd>(ys), as<MatrixXd>(Xs), xpMat, observations, npcs, nqtn, coefficients, stderr, df, seqQTN);

  // parallelFor
  parallelFor(0, xpMat->ncol(), lm, 1000);

  return List::create(_["coefficients"] = coefficients,
                      _["stderr"] = stderr,
                      _["df"] = df,
                      _["seqQTN"] = seqQTN);
}
