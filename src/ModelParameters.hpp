/**
 * @file ModelParameters.hpp
 * @author The Ecomeristem Development Team
 * See the AUTHORS file
 */

/*
 * Copyright (C) 2012-2017 ULCO http://www.univ-littoral.fr
 * Copyright (C) 2005-2017 Cirad http://www.cirad.fr
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef MODELPARAMETERS_HPP
#define MODELPARAMETERS_HPP 1

#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace ecomeristem {
struct Climate {
   double Temperature;
   double Par;
   double Etp;
   double Irrigation;
   double P;

   Climate( double Temperature, double Par, double Etp, double Irrigation,
            double P ) :
      Temperature( Temperature ), Par( Par ), Etp( Etp ), Irrigation( Irrigation ),
      P( P )
   { }
};

class ModelParameters {
 public:
   ModelParameters()
   { }

   virtual ~ModelParameters()
   { }


   double get( const std::string &paramName ) const
   {
      std::map < std::string, double >::const_iterator it;
      it = mParams.find( paramName );

      if( it == mParams.end() )
         std::cout << "Warning: no value for " << paramName << std::endl;

      return ( it == mParams.end() ) ? 0 : it->second;
   }

   Climate get( double time ) const
   {
      return meteoValues[time-beginDate];
   }


   inline void set( const std::string &key, const double &value )
   {
      mParams[key] = value;
   }

   inline void clear()
   {
      mParams.clear();
   }

   
   std::vector < Climate > meteoValues;
   std::map < std::string, double > mParams;//!< Represent the parameters.
public:
    std::map < std::string, double > * getRawParameters() { return &mParams; }
    std::vector < Climate > * getMeteoValues() { return &meteoValues; }
    double beginDate;
};

}

#endif
