#pragma once

template <typename T>
class Matrix;

template <typename T>
class IMatrixMultiplicationStrategy {
 public:
  virtual ~IMatrixMultiplicationStrategy() = default;
  virtual Matrix<T> multiply(const Matrix<T>& a, const Matrix<T>& b) const = 0;
};
