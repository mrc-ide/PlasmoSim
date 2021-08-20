
#pragma once

#include "misc_v10.h"

#include <vector>

//------------------------------------------------
// 2d array of integers
class array_2d_int  {
public:
  std::vector<std::vector<int>> arr;
  array_2d_int(int d1 = 0, int d2 = 0, int x = 0) {
    arr = std::vector<std::vector<int>>(d1, std::vector<int>(d2, x));
  };
  std::vector<int> & operator [](int i) {return arr[i];}
};

//------------------------------------------------
// 2d array of doubles
class array_2d_double  {
public:
  std::vector<std::vector<double>> arr;
  array_2d_double(int d1 = 0, int d2 = 0, int x = 0) {
    arr = std::vector<std::vector<double>>(d1, std::vector<double>(d2, x));
  }
  std::vector<double> & operator [](int i) {return arr[i];}
};

//------------------------------------------------
// 3d array of integers
class array_3d_int  {
public:
  std::vector<std::vector<std::vector<int>>> arr;
  array_3d_int(int d1 = 0, int d2 = 0, int d3 = 0, int x = 0) {
    arr = std::vector<std::vector<std::vector<int>>>(d1, std::vector<std::vector<int>>(d2, std::vector<int>(d3, x)));
  }
  std::vector<std::vector<int>> & operator [](int i) {return arr[i];}
};

//------------------------------------------------
// 3d array of doubles
class array_3d_double  {
public:
  std::vector<std::vector<std::vector<double>>> arr;
  array_3d_double(int d1 = 0, int d2 = 0, int d3 = 0, int x = 0) {
    arr = std::vector<std::vector<std::vector<double>>>(d1, std::vector<std::vector<double>>(d2, std::vector<double>(d3, x)));
  }
  std::vector<std::vector<double>> & operator [](int i) {return arr[i];}
};

//------------------------------------------------
// 4d array of integers
class array_4d_int  {
public:
  std::vector<std::vector<std::vector<std::vector<int>>>> arr;
  array_4d_int(int d1 = 0, int d2 = 0, int d3 = 0, int d4 = 0, int x = 0) {
    arr = std::vector<std::vector<std::vector<std::vector<int>>>>(d1, std::vector<std::vector<std::vector<int>>>(d2, std::vector<std::vector<int>>(d3, std::vector<int>(d4, x))));
  }
  std::vector<std::vector<std::vector<int>>> & operator [](int i) {return arr[i];}
};

//------------------------------------------------
// 4d array of doubles
class array_4d_double  {
public:
  std::vector<std::vector<std::vector<std::vector<double>>>> arr;
  array_4d_double(int d1 = 0, int d2 = 0, int d3 = 0, int d4 = 0, int x = 0) {
    arr = std::vector<std::vector<std::vector<std::vector<double>>>>(d1, std::vector<std::vector<std::vector<double>>>(d2, std::vector<std::vector<double>>(d3, std::vector<double>(d4, x))));
  }
  std::vector<std::vector<std::vector<double>>> & operator [](int i) {return arr[i];}
};

//------------------------------------------------
// 5d array of integers
class array_5d_int  {
public:
  std::vector<std::vector<std::vector<std::vector<std::vector<int>>>>> arr;
  array_5d_int(int d1 = 0, int d2 = 0, int d3 = 0, int d4 = 0, int d5 = 0, int x = 0) {
    arr = std::vector<std::vector<std::vector<std::vector<std::vector<int>>>>>(d1, std::vector<std::vector<std::vector<std::vector<int>>>>(d2, std::vector<std::vector<std::vector<int>>>(d3, std::vector<std::vector<int>>(d4, std::vector<int>(d5, x)))));
  }
  std::vector<std::vector<std::vector<std::vector<int>>>> & operator [](int i) {return arr[i];}
};

//------------------------------------------------
// 5d array of doubles
class array_5d_double  {
public:
  std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> arr;
  array_5d_double(int d1 = 0, int d2 = 0, int d3 = 0, int d4 = 0, int d5 = 0, int x = 0) {
    arr = std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>(d1, std::vector<std::vector<std::vector<std::vector<double>>>>(d2, std::vector<std::vector<std::vector<double>>>(d3, std::vector<std::vector<double>>(d4, std::vector<double>(d5, x)))));
  }
  std::vector<std::vector<std::vector<std::vector<double>>>> & operator [](int i) {return arr[i];}
};

