#include <iostream>
#include "flann/flann.hpp"
#include "typeinfo"

#include "flann/util/matrix.h"

using namespace std;
namespace vip
{

template <class T>
struct ApertureConicalDistance
{
  typedef bool is_kdtree_distance;
  typedef T ElementType;
  typedef double ResultType;
  template <typename Iterator1, typename Iterator2>
  ResultType operator()(Iterator1 a, Iterator2 b, size_t size,
                        ResultType /*worst_dist*/ = -1) const
  {
    std::cout << "size: " << size << std::endl;
    ResultType result = ResultType();
    ResultType diff;
    for (size_t i = 0; i < size; ++i)
    {
      diff = *a++ - *b++;
      result += diff * diff;
    }
    return result;
  }
  template <typename U, typename V>
  inline ResultType accum_dist(const U &a, const V &b, int) const
  {
    std::cout << "a type: " << typeid(U).name() << std::endl;
    return (a - b) * (a - b);
  }
};
} // namespace flann

template <class T>
struct ApertureConicalDistance1
{
  typedef bool is_kdtree_distance;
  typedef T ElementType;
  typedef double ResultType;
  template <typename Iterator1, typename Iterator2>
  ResultType operator()(Iterator1 a, Iterator2 b, size_t size,
                        ResultType /*worst_dist*/ = -1) const
  {
    std::cout << "size: " << size << std::endl;
    ResultType result = ResultType();
    ResultType diff;
    for (size_t i = 0; i < size; ++i)
    {
      diff = *a++ - *b++;
      result += diff * diff;
    }
    return result;
  }
  template <typename U, typename V>
  inline ResultType accum_dist(const U &a, const V &b, int) const
  {
    std::cout << "a type: " << typeid(U).name() << std::endl;
    return (a - b) * (a - b);
  }
};

template <class T>
struct ApertureConicalDistance
{
  typedef bool is_kdtree_distance;
  typedef T ElementType;
  typedef double ResultType;
  // double spheretoxyz(std::vector<double> &s, std::vector<double> &xyz)
  //     const
  // {
  //   double rsin = s.at(0) * std::sin(s.at(1));
  //   xyz.push_back(rsin * std::cos(s.at(2)));    // x= r*sintheta*cospha
  //   xyz.push_back(rsin * std::sin(s.at(2)));    // y = r*sintheta*cospha
  //   xyz.push_back(s.at(0) * std::cos(s.at(1))); // z = r*costheta
  // }
  template <typename Iterator1, typename Iterator2>
  ResultType operator()(Iterator1 a, Iterator2 b, size_t size, ResultType /*worst_dist*/ = -1)const
  {
    // ResultType result = 0;
    std::cout << "size: " << size << std::endl;
    std::cout << "*a0: " << *a << std::endl;
    std::cout << "*b0: " << *b << std::endl;
    ResultType result = ResultType();
    ResultType diff;
    for (size_t i = 0; i < size; ++i)
    {
      diff = *a++ - *b++;
      result += diff * diff;
    }
 
    // std::vector<double> axyz(3), bxyz(3);
    // std::vector<double> as{*a, *(a + 1), *(a + 2)};
    // std::vector<double> bs{*b, *(b + 1), *(b + 2)};

    // spheretoxyz(as, axyz);
    // spheretoxyz(bs, bxyz);

    // double rx = *a;
    // double ry = *b;

    // double adotb = axyz[0] * bxyz[0] + axyz[1] * bxyz[1] + axyz[2] * bxyz[0];
    // if ((rx * ry) != 0)
    //   result = std::acos(adotb / (rx * ry));
    // else
    // {
    //   result = 0;
    //   std::cerr << "warninng: rx*ry is zero !" << std::endl;
    // }

    //       result = std::acos((a*b)/(a.mod*b.mod));
    return result;
  }

  template <typename U, typename V>
  inline ResultType accum_dist(const U &a, const V &b, int) const
  {
    return (a - b) * (a - b);
  }
};
struct MyPoint
{
  MyPoint(double x, double y, double z, double g, double h)
  {
    this->x = x;
    this->y = y;
    this->z = z;
    this->g = g;
    this->h = h;
  }
  double x, y, z, g, h;
};

int main()
{
  vector<MyPoint> points;
  points.push_back(MyPoint(0, 0, 0, 0, 0));
  points.push_back(MyPoint(1, 1, 1, 1, 1));
  points.push_back(MyPoint(-1, -1, -1, 1, 1));
  flann::Matrix<double> points_mat = flann::Matrix<double>(&points[0].x, points.size(), 5);

std::cout<<"points_mat0: " << *(points_mat.ptr()+5) << std::endl;
  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 5; ++j)
    {
      cout << points_mat[i][j] << "\t";
    }
    cout << endl;
  }

  int nn = 1;

  flann::Matrix<double> dataset = flann::Matrix<double>(&points[0].x, points.size(), 5);
  flann::Matrix<double> query = flann::Matrix<double>(&points[2].x, 1, 5);

  flann::Matrix<int> indices(new int[query.rows * nn], query.rows, nn);
  flann::Matrix<double> dists(new double[query.rows * nn], query.rows, nn);

  // flann::Index<vip::L2_test<double>> index(dataset, flann::KDTreeIndexParams(1));
   flann::Index<ApertureConicalDistance<double>> index(dataset, flann::KDTreeIndexParams(1)); 
  index.buildIndex();
  // do a knn search, using 128 checks

  index.knnSearch(query, indices, dists, nn, flann::SearchParams(-1));

  int i = 0;
  for (auto it = &indices[0][0]; i < indices.rows; i++, it++)
    cout << "indices1: " << *it << "\t" << std::endl;

  // delete[] dataset.ptr();
  // delete[] query.ptr();
  // delete[] indices.ptr();
  // delete[] dists.ptr();

  return 0;
}
