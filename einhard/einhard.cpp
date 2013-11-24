/**
 * @file
 *
 * Symbol declaration file for the view symbols needed by Einhard.
 *
 * Copyright 2010 Matthias Bach <marix@marix.org>
 *
 * This file is part of Einhard.
 *
 * Einhard is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Einhard is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Einhard.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "einhard.hpp"

using namespace einhard;

#define _COLOR(name, code) \
char const _##name::ANSI[] = "\33[" code "m"

_COLOR(DGray,   "01;30");
_COLOR(Black,   "00;30");
_COLOR(Red,     "01;31");
_COLOR(DRed,    "00;31");
_COLOR(Green,   "01;32");
_COLOR(DGreen,  "00;32");
_COLOR(Yellow,  "01;33");
_COLOR(Orange,  "00;33");
_COLOR(Blue,    "01;34");
_COLOR(DBlue,   "00;34");
_COLOR(Magenta, "01;35");
_COLOR(DMagenta, "00;35");
_COLOR(Cyan,    "01;36");
_COLOR(DCyan,   "00;36");
_COLOR(White,   "01;37");
_COLOR(Gray,    "00;37");
_COLOR(NoColor, "0"    );
#undef _COLOR

