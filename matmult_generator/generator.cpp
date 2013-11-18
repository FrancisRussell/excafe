#include <cstdlib>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include <cassert>
#include <excafe/numeric/excafe_expression.hpp>
#include <excafe/mp/cln_conversions.hpp>
#include <excafe/cse/cse_optimiser.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim_all.hpp>
#include <cln/cln.h>
#include "mat_mult_code_generator.hpp"
#include "vector_entry.hpp"

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
};

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

    excafe::cse::CSEOptimiser<VectorEntry> optimiser(polys.begin(), polys.end());
    MatMultCodeGenerator generator(std::cout, mat.numCols(), mat.numRows());

    generator.outputPrefix();
    optimiser.accept(generator);
    generator.outputPostfix();

    return EXIT_SUCCESS;
  }
  else
  {
    std::cerr << "usage: generator matrix_file" << std::endl;
    return EXIT_FAILURE;
  }
}
