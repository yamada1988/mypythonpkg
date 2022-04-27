##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2013 Stanford University and the Authors
#
# Authors: Robert McGibbon
# Contributors:
#
# MDTraj is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 2.1
# of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with MDTraj. If not, see <http://www.gnu.org/licenses/>.
##############################################################################

"""OpenMM Reporter for saving the state of a molecular dynamics simulation
through time in the GROMACS TRR format
"""

##############################################################################
# Imports
##############################################################################

from __future__ import print_function, division
from mdtraj.formats.trr import TRRTrajectoryFile
from .mybasereporter import _MyBaseReporter

##############################################################################
# Imports
##############################################################################


class TRRReporter(_MyBaseReporter):
    """TRRReporter stores a molecular dynamics trajectory in the GROMACS
    TRR format.

    Parameters
    ----------
    file : str, or TRRTrajectoryFile
        Either an open TRRTrajectoryFile object to write to, or a string
        specifying the filename of a new TRR file
    reportInterval : int
        The interval (in time steps) at which to write frames.
    coordinates : bool
        Whether to write the coordinates to the file.
    time : bool
        Whether to write the current time to the file.
    cell : bool
        Whether to write the current unitcell dimensions to the file.
    atomSubset : array_like, default=None
        Only write a subset of the atoms, with these (zero based) indices
        to the file. If None, *all* of the atoms will be written.

    Examples
    --------
    >>> simulation = Simulation(topology, system, integrator)
    >>> trr_reporter = TRRReporter('traj.trr', 100)
    >>> simulation.reporters.append(trr_reporter)
    >>> simulation.step(10000)

    >>> traj = mdtraj.trajectory.load('traj.trr')
    """
    @property
    def backend(self):
        return TRRTrajectoryFile

    def __init__(self, file, reportInterval, coordinates=True, 
                 cell=True, time=True, atomSubset=None):
        """Create a TRRReporter.
        """
        super(TRRReporter, self).__init__(file, reportInterval,
            coordinates, cell, time, potentialEnergy=False, kineticEnergy=False,
            temperature=False, velocities=False, atomSubset=atomSubset)
