#include "LogDouble.hpp"

LogDouble::LogDouble()
{
  log_absolute_value = std::numeric_limits<double>::quiet_NaN();
  sign = 1;
}

LogDouble::LogDouble(double x)
{
  if ( x == 0.0 )
    sign = 1;
  else
    sign = (signed char) sgn(x);

  log_absolute_value = log(fabs(x));
}

LogDouble LogDouble::create_from_log_absolute_value(double log_absolute_value_param)
{
  LogDouble result;
  result.log_absolute_value = log_absolute_value_param;
  result.sign = 1;
  return result;
}

double LogDouble::logaddexp(double log_a, double log_b)
{
  // returns the log( exp(log_a) + exp(log_b) )

  // if both are infinite, taking the difference will result in a NaN;
  // simply return infinity
  if (isinf(log_a) && isinf(log_b))
    return log_a;

  if ( log_a > log_b )
    return logaddexp_first_larger(log_a, log_b);
  else
    return logaddexp_first_larger(log_b, log_a);
}

double LogDouble::logaddexp_first_larger(double log_a, double log_b)
{
  // note: this check is for internal testing, and can be removed for greater speed
  if ( log_a < log_b )
    throw LogDoubleNumericException("logaddexp_first_larger with first argument < second");
  if ( log_a == -std::numeric_limits<double>::infinity() )
    return log_b;
  return log1p( exp(log_b-log_a) ) + log_a;
}

double LogDouble::logsubabsexp(double log_a, double log_b)
{
  // returns the log( abs( exp(log_a) - exp(log_b) ) )
  if ( log_a > log_b )
    return logsubexp_first_larger(log_a, log_b);
  else
    return logsubexp_first_larger(log_b, log_a);
}

double LogDouble::logsubexp_first_larger(double log_a, double log_b)
{
  // note: this check is for internal testing, and can be removed for greater speed
  if (log_a < log_b)
    throw LogDoubleNumericException("First argument must be larger");
  if ( log_a == -std::numeric_limits<double>::infinity() )
    return log_b;
  return log1p( -exp(log_b-log_a) ) + log_a;
}

const LogDouble & LogDouble::operator +=(const LogDouble & rhs)
{
  // if the signs are the same, simply use logaddexp
  if (sign == rhs.sign)
    log_absolute_value = logaddexp(log_absolute_value, rhs.log_absolute_value);
  else
    {
      double new_log_absolute_value = logsubabsexp(log_absolute_value, rhs.log_absolute_value);
      if ( log_absolute_value < rhs.log_absolute_value )
	// *this "loses"
	sign *= -1;
      log_absolute_value = new_log_absolute_value;
    }
  return *this;
}

const LogDouble & LogDouble::operator -=(const LogDouble & rhs)
{
  // if the signs are different, they will be the same after negation
  if (sign != rhs.sign)
    log_absolute_value = logaddexp(log_absolute_value, rhs.log_absolute_value);
  else
    {
      double new_log_absolute_value = logsubabsexp(log_absolute_value, rhs.log_absolute_value);
      if ( log_absolute_value < rhs.log_absolute_value )
	// *this "loses"
	sign *= -1;
      log_absolute_value = new_log_absolute_value;
    }
  return *this;
}
const LogDouble & LogDouble::operator *=(const LogDouble & rhs)
{
  sign *= rhs.sign;
  log_absolute_value += rhs.log_absolute_value;
  return *this;
}
const LogDouble & LogDouble::operator /=(const LogDouble & rhs)
{
  sign *= rhs.sign;
  log_absolute_value -= rhs.log_absolute_value;
  return *this;
}

LogDouble LogDouble::operator -() const
{
  LogDouble result(*this);
  result.sign = -result.sign;
  return result;
}

bool LogDouble::operator <(LogDouble rhs) const
{
  return (sign < rhs.sign) || ( sign == rhs.sign && ( (sign == 1 && log_absolute_value < rhs.log_absolute_value) || (sign == -1 && log_absolute_value > rhs.log_absolute_value) ) );
}

bool LogDouble::operator ==(LogDouble rhs) const
{
  // magnitudes must be equal, and if nonzero, signs must be equal (if
  // it's zero, disregard sign)
  return log_absolute_value == rhs.log_absolute_value && (sign == rhs.sign || (log_absolute_value == -std::numeric_limits<double>::infinity() ) );
}

// +, -, *, /
LogDouble operator +(LogDouble lhs, LogDouble rhs)
{
  lhs += rhs;
  return lhs;
}

LogDouble operator -(LogDouble lhs, LogDouble rhs)
{
  lhs -= rhs;
  return lhs;
}

LogDouble operator *(LogDouble lhs, LogDouble rhs)
{
  lhs *= rhs;
  return lhs;
}

LogDouble operator /(LogDouble lhs, LogDouble rhs)
{
  lhs /= rhs;
  return lhs;
}

LogDouble exp(LogDouble rhs)
{
  rhs.log_absolute_value = double(rhs);
  rhs.sign = 1;
  return rhs;
}

std::ostream & operator <<(std::ostream & os, LogDouble rhs)
{
  if (rhs.sign == -1)
    os << '-';
  os << "exp(" << rhs.log_absolute_value << ")~" << double(rhs);
  return os;
}

LogDouble pow(LogDouble lhs, LogDouble rhs)
{
  if (lhs.get_sign() <= 0)
    throw LogDouble::LogDoubleNumericException("Sign must be greater than zero");

  // for all x, x^0 --> 1
  if ( rhs.get_log_absolute_value() == -std::numeric_limits<double>::infinity() )
    return LogDouble(1.0);
  // for all y>0 (guaranteed by previous check), 0^y --> 0
  if ( lhs.get_log_absolute_value() == -std::numeric_limits<double>::infinity() )
    return LogDouble(0.0);
  return LogDouble::create_from_log_absolute_value(double(rhs) * lhs.get_log_absolute_value());
}
