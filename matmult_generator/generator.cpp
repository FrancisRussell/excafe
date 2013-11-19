#include <cstdlib>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include <cassert>
#include <algorithm>
#include <excafe/numeric/excafe_expression.hpp>
#include <excafe/mp/cln_conversions.hpp>
#include <excafe/cse/cse_optimiser.hpp>
#include <excafe/codegen/dynamic_cxx.hpp>
#include <excafe/util/timer.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim_all.hpp>
#include <cln/cln.h>
#include "mat_mult_code_generator.hpp"
#include "vector_entry.hpp"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <cblas.h>

#include <malloc.h>

typedef void (*mult_function_t)(const double* x, double *b, size_t count);
static const int BENCHMARK_VECTOR_COUNT = 500000;
static const int REPETITIONS = 100;

template<typename T>
class RowMajorMatrix
{
private:
  typedef T value_type;
  size_t rows;
  size_t cols;
  std::vector<value_type> data;

public:
  RowMajorMatrix() : rows(0), cols(0)
  {
  }

  RowMajorMatrix(size_t _rows, size_t _cols) :
    rows(_rows), cols(_cols), data(rows * cols)
  {
  }

  RowMajorMatrix(std::istream& stream) : rows(0), cols(0)
  {
    while(stream.good())
    {
      std::string line;
      std::getline(stream, line);
      boost::algorithm::trim_all(line);

      if (line.empty())
        break;

      std::istringstream lineStream(line);
      cols = 0;

      while(lineStream.good())
      {
        value_type value;
        lineStream >> value;
        ++cols;
        data.push_back(value);
      }

      ++rows;
    }
    assert(data.size() == rows * cols);
  }

  size_t numRows() const
  {
    return rows;
  }

  size_t numCols() const
  {
    return cols;
  }

  value_type& operator()(size_t row, size_t col)
  {
    return data[row * cols + col];
  }

  const value_type operator()(size_t row, size_t col) const
  {
    return data[row * cols + col];
  }

  const double *getData() const
  {
    return &data[0];
  }
};

static double timeGenerated(const int m, const int n, const int k, const double *mat, const double* in, double *out,
  const mult_function_t func)
{
  excafe::util::Timer timer;
  timer.start();
  for(int rep = 0; rep < REPETITIONS; ++rep)
    func(in, out, n);
  timer.stop();

  return timer.getSeconds() / REPETITIONS;
}

static double timeBLAS(const int m, const int n, const int k, const double *mat, const double* in, double *out)
{
  excafe::util::Timer timer;
  timer.start();
  for(int rep = 0; rep < REPETITIONS; ++rep)
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0, mat, k, in, n, 0.0, out, n);
  timer.stop();

  return timer.getSeconds() / REPETITIONS;
}

static double computeDelta(const double* a, const double *b, const size_t size)
{
  double delta = 0.0;

  for(size_t i = 0; i < size; ++i)
  {
    const double diff = a[i] - b[i];
    delta += diff * diff;
  }

  return std::sqrt(delta);
}

int main(int argc, char **argv)
{
  if (argc == 2)
  {
    std::fstream stream(argv[1]);
    RowMajorMatrix<double> mat(stream);
    std::cout << "Successfully read matrix of dimensions " << mat.numRows() << " x " << mat.numCols() << "." << std::endl;

    typedef excafe::ExcafeExpression<VectorEntry> poly_t;
    std::vector<poly_t> polys;

    for(size_t row = 0; row < mat.numRows(); ++row)
    {
      poly_t poly;

      for(size_t col = 0; col < mat.numCols(); ++col)
      {
        const double coeff = mat(row, col);
        cln::cl_RA coeffAsCLNRational = cln::rationalize(coeff);
        excafe::mp::Rational rational = excafe::mp::fromCLN(coeffAsCLNRational);
        poly += rational * poly_t(VectorEntry(col));
      }

      polys.push_back(poly);
    }

    std::ostringstream codeStream;
    excafe::cse::CSEOptimiser<VectorEntry> optimiser(polys.begin(), polys.end());
    MatMultCodeGenerator generator(codeStream, mat.numCols(), mat.numRows());
    generator.output(optimiser);

    std::cout << codeStream.str() << std::endl;

    excafe::codegen::DynamicCXX dynamicCXX(codeStream.str());
    dynamicCXX.compileAndLoad();

    mult_function_t generatedFunction = reinterpret_cast<mult_function_t>(dynamicCXX.getFunction("frgmv"));

    std::vector<double> in(BENCHMARK_VECTOR_COUNT * mat.numCols());
    typedef boost::uniform_real<> coeff_distribution_type;
    typedef boost::variate_generator<boost::mt19937&, coeff_distribution_type> coeff_gen_type;
    boost::mt19937 rng;
    coeff_gen_type coefficientGenerator(rng, coeff_distribution_type(-1.0, 1.0));
    std::generate(in.begin(), in.end(), coefficientGenerator);

    const size_t numOutputElements = BENCHMARK_VECTOR_COUNT * mat.numRows();
    double *outGenerated = (double*) memalign(16, sizeof(double) * numOutputElements);
    double *outReference = (double*) memalign(16, sizeof(double) * numOutputElements);
    const double generatedTime = timeGenerated(mat.numRows(), BENCHMARK_VECTOR_COUNT, mat.numCols(),
      mat.getData(), &in[0], &outGenerated[0], generatedFunction);

    std::cout << "Generated code execution time: " << generatedTime << " seconds." << std::endl;

    const double blasTime =
      timeBLAS(mat.numRows(), BENCHMARK_VECTOR_COUNT, mat.numCols(), mat.getData(), &in[0], &outReference[0]);

    std::cout << "BLAS dgemm execution time: " << blasTime << " seconds." << std::endl;
    std::cout << "Delta: " << computeDelta(&outReference[0], &outGenerated[0], numOutputElements) << std::endl;

    free(outGenerated);
    free(outReference);

    return EXIT_SUCCESS;
  }
  else
  {
    std::cerr << "usage: generator matrix_file" << std::endl;
    return EXIT_FAILURE;
  }
}
