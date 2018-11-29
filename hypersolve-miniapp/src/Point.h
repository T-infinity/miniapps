#pragma once

#include <array>

namespace HS {
    template<typename T>
    class Point {
    public:
        std::array<T, 3> coordinates;
        Point() = default;
        Point(const T& x, const T& y, const T& z) {
            coordinates[0] = x;
            coordinates[1] = y;
            coordinates[2] = z;
        }
        T& operator[](int index) {
            return coordinates[index];
        }
        const T& operator[](int index) const {
            return coordinates[index];
        }
        Point<T> operator+(const Point<T>& p) const {
            std::array<T, 3> new_point = {p[0], p[1], p[2]};
            new_point[0] += coordinates[0];
            new_point[1] += coordinates[1];
            new_point[2] += coordinates[2];
            return Point(new_point[0], new_point[1], new_point[2]);
        }
        Point<T> operator-(const Point<T>& p) const {
            std::array<T, 3> new_point = {p[0], p[1], p[2]};
            new_point[0] -= coordinates[0];
            new_point[1] -= coordinates[1];
            new_point[2] -= coordinates[2];
            return Point(new_point[0], new_point[1], new_point[2]);
        }
        Point<T> operator*(const T& scalar) const {
            return Point(scalar * coordinates[0],
                         scalar * coordinates[1],
                         scalar * coordinates[2]);
        }
        Point<T> operator/(const T& scalar) const {
            return Point(coordinates[0] / scalar,
                         coordinates[1] / scalar,
                         coordinates[2] / scalar);
        }
        static Point<T> cross(const Point<T> &point1, const Point<T> &point2) {
            return Point(point1[1] * point2[2] - point2[1] * point1[2],
                         point2[0] * point1[2] - point1[0] * point2[2],
                         point1[0] * point2[1] - point2[0] * point1[1]);
        }

        static T dot(const Point<T> &point1, const Point<T> &point2) {
            return point1[0] * point2[0] + point1[1] * point2[1] + point1[2] * point2[2];
        }
    };
}

template<typename T>
HS::Point<T> operator*(const T &scalar, const HS::Point<T> &p) {
    return HS::Point<T>(scalar * p[0],
                        scalar * p[1],
                        scalar * p[2]);
}

