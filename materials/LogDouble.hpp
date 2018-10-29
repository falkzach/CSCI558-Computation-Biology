// The MIT License (MIT)

// Copyright (c) 2014 Oliver Serang

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

#ifndef _LOGDOUBLE_H
#define _LOGDOUBLE_H

#include <math.h>
#include <limits>
#include <iostream>

#include "ComparableMixin.hpp"

template <typename T> int sgn(T val) {
    return (val > T(0)) - (val < T(0));
}

class LogDouble : public ComparableMixin<LogDouble> {
protected:
  signed char sign;
  double log_absolute_value;
  
  static double logaddexp(double log_a, double log_b);
  static double logaddexp_first_larger(double log_a, double log_b);
  static double logsubabsexp(double log_a, double log_b);
  static double logsubexp_first_larger(double log_a, double log_b);
public:
  LogDouble();
  explicit LogDouble(double x);

  static LogDouble create_from_log_absolute_value(double log_absolute_value_param);

  // +=, -=, *=, /=
  const LogDouble & operator +=(const LogDouble & rhs);
  const LogDouble & operator -=(const LogDouble & rhs);
  const LogDouble & operator *=(const LogDouble & rhs);
  const LogDouble & operator /=(const LogDouble & rhs);
  LogDouble operator -() const;

  explicit operator double() const { return sign * exp(log_absolute_value); }

  double get_log_absolute_value() const { return log_absolute_value; }
  double get_sign() const { return sign; }

  bool operator <(LogDouble rhs) const;
  bool operator ==(LogDouble rhs) const;
  //  bool operator !=(LogDouble rhs) const { return ! (*this == rhs); }

  //  static bool is_nan(LogDouble x) { return isnan(double(x)); }
  //  static bool is_inf(LogDouble x) { return isinf(x.get_log_absolute_value()); }

  friend LogDouble exp(LogDouble rhs);
  friend std::ostream & operator <<(std::ostream & os, LogDouble rhs);

  class LogDoubleNumericException : public std::exception {
  private:
    std::string my_error_msg;
  public:
    LogDoubleNumericException(const std::string & str) throw() {
      my_error_msg = str;
    }
    virtual ~LogDoubleNumericException() throw() { }
    virtual const char* what() const throw() {
      return my_error_msg.c_str();
    }
  };

  // Note: this function is a replacement for log1p (if your compiler doesn't have it)
  static double logOnePlusX(double x) {
    if (x < -1.0)
      throw LogDouble::LogDoubleNumericException("Cannot execute log of negative value (on the reals)");
    
    if (fabs(x) > 1e-4)
      // x is far enough from zero to accurately use direct computation
      return log(1.0 + x);
    
    // Use a Taylor approximation:
    // log(1+x) = x - x^2/2 + x^3/3 - x^4/4 + ...
    // \approx x - x^2/2
    // 
    // Qualitatively, when x is so close to zero, each next term has a vanishing
    // effect.
    //
    // The error will be \leq |x^3/3| (|x/4| + |x^2/5| + |x^3/6| + ...)
    // \leq |x^3/3| (|x| + |x|^2 + |x|^3 + ...) \leq |x^3/3| 1/(1-x)
    //
    // When x \approx 0, 1/(1-x) is very close to 1.
    //
    // The absolute error is thus less than x^3, which is at or below
    // 1E-12 in the interval |x| < 1E-4.
    //
    // The relative error is thus bounded above by x^3 / log(1+x); the
    // derivative of this function has only one zero in the region |x|
    // < 1E-4, but it corresponds to a minimum rather than a maximum,
    // and thus a maximum must occur at a boundary point. Thus, the
    // maximum relative error, which is achieved at x=1E-4, is
    // \leq 1E-12/log(1+1E-4).
    
    return (-0.5 * x + 1.0) * x;
  }
};

LogDouble operator +(LogDouble lhs, LogDouble rhs);
LogDouble operator -(LogDouble lhs, LogDouble rhs);
LogDouble operator *(LogDouble lhs, LogDouble rhs);
LogDouble operator /(LogDouble lhs, LogDouble rhs);

LogDouble pow(LogDouble lhs, LogDouble rhs);

#endif
