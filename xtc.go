// +build xtc

/*
 * xtc.go part of goMD
 *
 *
 *
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License along
 *  with this program; if not, write to the Free Software Foundation, Inc.,
 *  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *
 */

/*To the long life of the Ven. Khenpo Phuntzok Tenzin Rinpoche*/

package main

import (
	chem "github.com/rmera/gochem"
	"github.com/rmera/gochem/traj/xtc"
)

func OpenXTC(name string) (chem.Traj, error) {
	var traj chem.Traj
	var err error
	traj, err = xtc.New(name) //now Bartender will not compile without the xdrfile libraries, which sucks.
	return traj, err
}
