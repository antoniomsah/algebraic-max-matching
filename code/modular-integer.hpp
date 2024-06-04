#pragma once

#include <iostream>

template <const int MOD>
struct modular_int {
  long long x;

  modular_int() : x(0) {}
  modular_int(int _x) : x(_x % MOD < 0 ? _x % MOD + MOD : _x % MOD) {}

  void operator+=(modular_int rhs) {
    x += rhs.x;
    if (x >= MOD) x -= MOD;
  }

  void operator-=(modular_int rhs) {
    x -= rhs.x;
    if (x < 0) x += MOD;
  }

  void operator*=(modular_int rhs) {
    x *= rhs.x;
    x %= MOD;
  }

  void operator/=(modular_int rhs) { *this *= rhs.inv(); }

  modular_int operator+(const modular_int& rhs) {
    modular_int res = *this;
    res += rhs;
    return res;
  }

  modular_int operator-(const modular_int& rhs) {
    modular_int res = *this;
    res -= rhs;
    return res;
  }

  modular_int operator*(const modular_int& rhs) {
    modular_int res = *this;
    res *= rhs;
    return res;
  }

  modular_int operator/(const modular_int& rhs) {
    modular_int res = *this;
    res /= rhs;
    return res;
  }

  modular_int inv() { return this->pow(MOD - 2); }

  modular_int pow(int e) {
    modular_int res(1);
    for (modular_int p = *this; e > 0; e /= 2, p *= p)
      if (e % 2) res *= p;
    return res;
  }

  bool operator<(int y) const { return x < y; }
  bool operator>(int y) const { return x > y; }
  bool operator==(int y) const { return x == y; }
  bool operator!=(int y) const { return x != y; }

  bool operator<(modular_int y) const { return x < y.x; }
  bool operator>(modular_int y) const { return x > y.x; }
  bool operator==(modular_int y) const { return x == y.x; }
  bool operator!=(modular_int y) const { return x != y.x; }

  friend std::ostream& operator<<(std::ostream& os, modular_int x) {
    return (os << x.x);
  }
};