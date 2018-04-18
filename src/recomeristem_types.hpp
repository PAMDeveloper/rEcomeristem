/**
 * @file rcpp_samara.hpp
 * @author See the AUTHORS file
 */

/*
 * Copyright (C) 2012-2017 ULCO http://www.univ-littoral.fr
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

#define UNSAFE_RUN
#include "defines.hpp"

#include <utils/ParametersReader.hpp>
#include <utils/resultparser.h>
#include <utils/juliancalculator.h>
#include <plant/PlantModel.hpp>
#include <observer/PlantView.hpp>
#ifndef UNSAFE_RUN
#include <boost/property_tree/ptree.hpp>
#include <boost/lexical_cast.hpp>
#endif


