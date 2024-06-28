#pragma once

#include <exception>
#include <string>

struct SingularMatrixError : std::exception {
  SingularMatrixError() {}

  char const* what() const throw() { return "Matrix is singular."; }
};