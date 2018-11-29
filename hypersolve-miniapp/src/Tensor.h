#pragma once

#include <array>
#include "Point.h"

namespace HS {
    template<typename T, size_t rows, size_t columns>
    class Tensor {
    public:
        Tensor() {
            zero();
        }

        inline T& operator()(int row, int col) {
            return m[row][col];
        }

        inline const T& operator()(int row, int col) const {
            return m[row][col];
        }

        inline std::array<T, columns>& operator()(int row) {
            return m[row];
        }

        inline const std::array<T, columns>& operator()(int row) const {
            return m[row];
        }

        inline void assignRow(int row_index, const std::array<T, columns> &values) {
            m[row_index] = values;
        }

        template<typename S>
        inline void operator*=(const S &scalar) {
            for (size_t row = 0; row < rows; ++row) {
                for (size_t col = 0; col < columns; ++col) {
                    m[row][col] *= scalar;
                }
            }
        }

        template<typename T2>
        inline void operator+=(const Tensor<T2, rows, columns>& tensor) {
            for (size_t row = 0; row < rows; ++row) {
                for (size_t col = 0; col < columns; ++col) {
                    m[row][col] += tensor(row, col);
                }
            }
        }

        template<typename T2>
        inline Tensor<T, rows, columns> operator+(const Tensor<T2, rows, columns>& tensor) const {
            Tensor<T, rows, columns> output;
            for (size_t row = 0; row < rows; ++row) {
                for (size_t col = 0; col < columns; ++col) {
                    output(row, col) = (*this)(row, col) + tensor(row, col);
                }
            }
            return output;
        }

        template<typename T2>
        inline void operator-=(const Tensor<T2, rows, columns>& tensor) {
            for (size_t row = 0; row < rows; ++row) {
                for (size_t col = 0; col < columns; ++col) {
                    m[row][col] -= tensor(row, col);
                }
            }
        }

        template<typename T2>
        inline Tensor<T, rows, columns> operator-(const Tensor<T2, rows, columns>& tensor) const {
            Tensor<T, rows, columns> output;
            for (size_t row = 0; row < rows; ++row) {
                for (size_t col = 0; col < columns; ++col) {
                    output(row, col) = (*this)(row, col) - tensor(row, col);
                }
            }
            return output;
        }

        template<typename T2>
        inline Tensor<T, rows, columns> operator*(const T2 &scalar) const {
            Tensor<T, rows, columns> output;
            for (size_t row = 0; row < rows; ++row)
                for (size_t col = 0; col < columns; ++col)
                    output(row, col) = (*this)(row, col) * scalar;
            return output;
        }

        template<typename T2>
        inline Tensor<T, rows, columns> operator*(const Tensor<T2, rows, columns>& tensor) const {
            Tensor<T, rows, columns> output;
            for (size_t row = 0; row < rows; ++row) {
                for (size_t col = 0; col < columns; ++col) {
                    output(row, col) = 0.0;
                    for (size_t i = 0; i < columns; ++i) {
                        output(row, col) += (*this)(row, i) * tensor(i, col);
                    }
                }
            }
            return output;
        }

        template<typename T2>
        inline std::array<T, rows> operator*(const std::array<T2, columns>& vec) const {
            std::array<T, rows> output;
            for (size_t row = 0; row < rows; ++row) {
                output[row] = 0.0;
                for (size_t col = 0; col < columns; ++col) {
                    output[row] += (*this)(row, col) * vec[col];
                }
            }
            return output;
        }

        template<typename T2>
        inline HS::Point<T> operator*(const HS::Point<T2>& vec) const {
            static_assert(columns == 3, "Tensor must be 3 columns");
            static_assert(rows == 3, "Tensor must have 3 rows");
            HS::Point<T> output;
            for (size_t row = 0; row < rows; ++row) {
                output[row] = 0.0;
                for (size_t col = 0; col < columns; ++col) {
                    output[row] += (*this)(row, col) * vec[col];
                }
            }
            return output;
        }

        inline void zero() {
            for (size_t i = 0; i < rows; ++i)
                for (size_t j = 0; j < rows; ++j)
                    m[i][j] = 0.0;
        }

        inline T trace() const {
            T trace(0.0);
            for (size_t i = 0; i < rows; ++i) {
                trace += m[i][i];
            }
            return trace;
        }

        inline Tensor<T, rows, columns> transpose() const {
            Tensor<T, rows, columns> output;
            for (size_t row = 0; row < rows; ++row)
                for (size_t col = 0; col < columns; ++col)
                    output(row, col) = m[col][row];
            return output;
        };

        template<typename T2>
        inline bool operator==(const Tensor<T2, rows, columns>& input) const {
            for (size_t row = 0; row < rows; ++row)
                for (size_t col = 0; col < columns; ++col)
                    if ((*this)(row, col) != input(row, col)) return false;
            return true;
        }

        template<typename T2>
        inline bool operator!=(const Tensor<T2, rows, columns>& input) const {
            for (size_t row = 0; row < rows; ++row)
                for (size_t col = 0; col < columns; ++col)
                    if ((*this)(row, col) != input(row, col)) return true;
            return false;
        }

        inline void print() const {
            for (size_t row = 0; row < rows; ++row) {
                for (size_t col = 0; col < columns; ++col)
                    printf("%e ", m[row][col]);
                printf("\n");
            }
        }

    private:
        std::array<std::array<T, rows>, columns> m;
    };
}

template <typename ScalarType, typename TensorType, size_t rows, size_t columns>
inline HS::Tensor<TensorType, rows, columns> operator*(const ScalarType& scalar, const HS::Tensor<TensorType, rows, columns>& tensor) {
    return tensor * scalar;
}

template <typename T, size_t rows, size_t columns>
inline std::array<T, columns> operator*(const std::array<T, rows>& vec, const HS::Tensor<T, rows, columns>& tensor) {
    std::array<T, rows> output;
    output.fill(0.0);
    for (size_t row = 0; row < rows; ++row) {
        for (size_t col = 0; col < columns; ++col) {
            output[col] += vec[row] * tensor(row, col);
        }
    }
    return output;
};

template <typename T>
inline HS::Point<T> operator*(const HS::Point<T>& vec, const HS::Tensor<T, 3, 3>& tensor) {
    HS::Point<T> output;
    for (int i = 0; i < 3; ++i)
        output[i] = 0.0;
    for (int row = 0; row < 3; ++row) {
        for (int col = 0; col < 3; ++col) {
            output[col] += vec[row] * tensor(row, col);
        }
    }
    return output;
};
