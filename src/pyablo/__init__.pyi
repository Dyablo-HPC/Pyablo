""" Pyablo Module"

.. curretmodule:: pyablo
.. autosummary::
   :toctree: _autosummary

   Line
   Snapshot
   XdmfReader
"""

from typing import overload

# Type aliases
# TODO: enventually switch to use numpy dtypes ?
Vec = tuple[int, int, int]
RealArray =  list[float]
IntArray = list[int]
UIntArray = list[int]  # no unsigned int in Python ?
BoolArray = list[bool]
UInt64Array = list[int]
VecArray = list[Vec]
RealTable = list[RealArray]
VecTable = list[VecArray]
BoundingBox = tuple[Vec, Vec]


class Line:
    """Class representing a line in the simulation domain."""
    Nl: int
    start: Vec
    end: Vec
    pos: VecArray
    rho: float
    prs: float
    E: float
    vel: VecArray
    cellIds: UIntArray


class Snapshot:
    """
    Class storing info about dyablo snapshots
     
    A general note : The strategy here is to never store arrays directly in
    ram except for coordinates and connectivity for rapid lookup. 
    Future optimization might include streaming or partial buffering
    but we do this for the moment to avoid any memory explosion for big runs
    
    """

    def close(self) -> None:
        """
        Closes all the handles in the file
        This is not done in the destructor to ensure compatibility
        with pybind11 !
        """

    def print(self) -> None:
        """Pretty prints a brief summary of the snapshot"""

    def setName(self, name: str) -> None:
        """
        Sets the associated name of the snapshot

        :param name: the name to give to the snapshot
        :type name: str
        """
    
    # def setTime(self, time: float) -> None:
    #     """
    #     Sets the associated time of the snapshot
        
    #     :param time the time of the current snapshot
    #     :type time: float
    #     """
    
    def setNDim(self, ndim: int) -> None:
        """
        Defines the number of dimensions of the snapshot.
        This also sets the size of an element connectivity in memory

        :param ndim: the number of dimensions. Should be 2 or 3.
        :type ndim: int
        """
    
    
    def getCellFromPosition(self, pos: Vec) -> int:
        """
        Finds which element corresponds to a specific location

        :param pos: the position of the element to probe
        :type pos: Vec

        :return: an integer corresponding to the cell id holding position pos
        :rtype: int
        """
    
    def getCellsFromPositions(self, pos: VecArray) -> UIntArray:
        """
        Finds which elements correspond to specific locations
        This function can be especially long on large data sets !

        :param pos: a vector of positions to identify
        :type pos: VecArray

        :return: a vector of integers corre
        :rtype: UIntArray
        """

    def getCellBoundingBox(self, iCell: int) -> BoundingBox:
        """
        Returns the bounding box corresponding to a cell
        
        :param iCell: the id of the cell to probe
        :type iCell: int

        :return: a pair of Vec storing the minimum and maximum coordinates of the bounding box
        :rtype: BoundingBox
        """

    def getNCells(self) -> int:
        """
        Returns the number of cells in the domain

        :return: the number of cells in the domain
        :rtype: int
        """

    def getUniqueCells(self, pos: VecArray) -> VecArray:
        """
        Returns the positions corresponding to the center of the cells 
        traversed by the points along pos

        :param pos: a vector of positions
        :type pos: VecArray

        :return: the center of the unique cells along the trajectory pos
        :rtype: VecArray
        """
    
    def hasAttribute(self, attribute: str) -> bool:
        """
        Returns whether or not the Snapshot contains an attribute
        
        :param attribute: the name of the attribute to test
        :type attribute: str

        :return: a boolean indicating if attribute is stored in the Snapshot
        :rtype: bool
        """
    
    def getDomainBoundingBox(self) -> BoundingBox:
        """
        Returns the bounding box of the domain
        
        :note: This is a very naive implementation where we consider only a cartesian-box
        domain ! This will not work using other geometries
        
        :return: Returns the bounding box of the domain
        :rtype: BoundingBox
        """
    
    def getCellsCenter(self, iCells: UIntArray) -> VecArray:
        """
        Returns the centers of given cells
        
        :param iCells: the indices of the cells to probe
        :type iCells: UIntArray

        :return: a vector of positions corresponding to the center of the cells
        :rtype: VecArray
        """

    def getCellsSize(self, iCells: UIntArray) -> VecArray:
        """
        Returns the size of a list of cells

        :param iCells: a vector of cells to probe
        :type iCells: UIntArray

        :return: a vector of Vec indicating the size of each cell probed
        :rtype: Vec
        """
    
    def getCellsVolume(self, iCells: UIntArray) -> RealArray:
        """
        Returns the volume/surface of a cell

        :param iCells: a vector of cells to probe
        :type iCells: UIntArray

        :return: a vector of doubles indicating the volume/surface of each cell probed
        :rtype: RealArray
        """

    def probeQuantity(self, pos: VecArray, attribute: str) -> RealArray:
        """
        Probes a location for a specified attribute.

        :param pos: the position to probe
        :type pos: VecArray

        :param attribute: the attribute to probe
        :type attribute: str

        :return: the chosen attribute at position pos
        :rtype: RealArray
        """
    
    def probeDensity(self, pos: VecArray) -> RealArray:
        """
        Probes a location for density.

        :param pos: the position to probe
        :type pos: VecArray

        :return: the density at position pos
        :rtype: RealArray
        """
    
    def probePressure(self, pos: VecArray) -> RealArray:
        """
        Probes a location for pressure.

        :param pos: the position to probe
        :type pos: VecArray

        :return: the pressure at position pos
        :rtype: RealArray
        """
    
    def probeEnergy(self, pos: VecArray) -> RealArray:
        """
        Probes a location for total energy.

        :param pos: the position to probe
        :type pos: VecArray

        :return: the total energy at position pos
        :rtype: RealArray
        """

    def probeMach(self, pos: VecArray) -> RealArray:
        """
        Probes a location for mach number.

        :param pos: the position to probe
        :type pos: VecArray

        :return: the mach number at position pos
        :rtype: RealArray
        """
    
    def probeMomentum(self, pos: VecArray) -> VecArray:
        """
        Probes a location for momentum.

        :param pos: the position to probe
        :type pos: VecArray

        :return: the momentum at position pos
        :rtype: VecArray
        """
    
    def probeVelocity(self, pos: VecArray) -> VecArray:
        """
        Probes a location for Velocity.

        :param pos: the position to probe
        :type pos: VecArray

        :return: the Velocity at position pos
        :rtype: VecArray
        """
    
    def probeRank(self, pos: VecArray) -> IntArray:
        """
        Probes multiple locations for the mpi-rank

        :param pos: a vector of positions to probe
        :type pos: VecArray
        
        :return: the mpi-rank at positions pos
        :rtype: IntArray
        """

    def probeLevel(self, pos: VecArray) -> IntArray:
        """
        Probes multiple locations for the refinement level
        
        :param pos: a vector of positions to probe
        :type pos: VecArray
        
        :return: the refinement level at positions pos
        :rtype: IntArray
        """
    
    def probeOctant(self, pos: VecArray) -> IntArray:
        """
        Probes multiple locations for the octant index
        
        :param pos: a vector of positions to probe
        :type pos: VecArray

        :return: the octant indices at positions pos
        :rtype: IntArray
        """
    
    def getQuantity(self, iCells: UIntArray, attribute: str) -> RealArray:
        """
        Extracts a quantity from a list of cells
        
        :param iCells: iCells the ids of the cells to extract
        :type iCells: UIntArray
        :param attribute: attribute the name of the field to extract
        :type attribute: str

        :return: a vector of the given quantity for each cell
        :rtype: RealArray
        """
    
    def getDensity(self, iCells: UIntArray) -> RealArray:
        """
        Extracts density from a list of cells
        
        :param iCells: the ids of the cells to extract
        :type iCells: UIntArray
        
        :return: a vector of densities for each cell
        :rtype: RealArray
        """

    def getPressure(self, iCells: UIntArray) -> RealArray:
        """
        Extracts the pressure from a list of cells
        
        :param iCells: the ids of the cells to extract
        :type iCells: UIntArray
        
        :return: a vector of pressures for each cell
        :rtype: RealArray
        """
    
    def getEnergy(self, iCells: UIntArray) -> RealArray:
        """
        Extracts the total energy from a list of cells
        
        :param iCells: the ids of the cells to extract
        :type iCells: UIntArray
        
        :return: a vector of energies for each cell
        :rtype: RealArray
        """
    
    def getMach(self, iCells: UIntArray) -> RealArray:
        """
        Extracts the Mach number from a list of cells
        
        :param iCells: the ids of the cells to extract
        :type iCells: UIntArray
        
        :return: a vector of Mach number for each cell
        :rtype: RealArray
        """

    def getVelocity(self, iCells: UIntArray) -> VecArray:
        """
        Extracts the velocities from a list of cells
        
        :param iCells: the ids of the cells to extract
        :type iCells: UIntArray
        
        :return: a vector of velocities for each cell
        :rtype: VecArray
        """

    def getMomentum(self, iCells: UIntArray) -> VecArray:
        """
        Extracts the momentum from a list of cells
        
        :param iCells: the ids of the cells to extract
        :type iCells: UIntArray
        
        :return: a vector of momentum for each cell
        :rtype: VecArray
        """
    
    def getLevel(self, iCells: UIntArray) -> IntArray:
        """
        Extracts the AMR level from a list of cells
        
        :param iCells: the ids of the cells to extract
        :type iCells: UIntArray
        
        :return: a vector of levels for each cell
        :rtype: IntArray
        """
    
    def getRank(self, iCells: UIntArray) -> IntArray:
        """
        Extracts the MPI rank from a list of cells
        
        :param iCells: the ids of the cells to extract
        :type iCells: UIntArray
        
        :return: a vector of rank for each cell
        :rtype: IntArray
        """
    
    def getOctant(self, iCells: UIntArray) -> IntArray:
        """
        Extracts the mesh octant from a list of cells
        
        :param iCells: the ids of the cells to extract
        :type iCells: UIntArray
        
        :return: a vector of octant ids for each cell
        :rtype: IntArray
        """
    
    def getRefinementCriterion(self, iCells: UIntArray) -> IntArray:
        """
        Returns the value of the refinement criterion at a set of positions
        
        :param iCells: the position vector where to probe the refinement criterion
        :type iCells: UIntArray

        :return: a vector of doubleing point values corresponding to the error for refinement at each position. 
        :rtype: IntArray

        :note: Result will be 0.0f for each position at the edge of the domain
        """
    
    def getTotalMass(self) -> float:
        """
        Returns the integrated total mass over the domain
        
        :return: the total mass in the domain
        :rtype: float
        """
    
    def getTotalEnergy(self) -> float:
        """
        Returns the integrated total energy over the domain
        
        :return: the total energy in the domain
        :rtype: float
        """
    
    def getTotalKineticEnergy(self) -> float:
        """
        Returns the integrated total kinetic energy over the domain
        
        :return: the total kinetic energy in the domain
        :rtype: float
        """
    
    def getTotalInternalEnergy(self, gamma: float) -> float:
        """
        Returns the integrated total internal energy over the domain
        
        :param gamma: adiabatic index
        :type gamma: float

        :return: the total internal energy in the domain
        :rtype: float
        """
    
    def getMaxMach(self) -> float:
        """
        Returns the maximum Mach number of the domain
        
        :return: the maximum Mach number of the flow in the domain
        :rtype: float
        """
    
    def getAverageMach(self) -> float:
        """
        Returns the average Mach number of the domain
        
        :return: the average Mach number of the flow in the domain
        :rtype: float
        """

    def getTime(self) -> float:
        """
        Gets the associated time of the snapshot
        
        :return: the time of the current snapshot
        :rtype: float
        """
    
    def readAllFloat(self, attribute: str) -> RealArray:
        """
        Reads all the data corresponding to a field
        
        :param attribute: the name of the field to read
        :type attribute: str
        
        :return: an array of double correponding to the linearized array of the field
        :rtype: RealArray
        """
    
    def mortonSort2d(self, vec: RealArray, iLevel: int, bx: int, by: int) -> RealArray:
        """
        Sorts an array in 2D extracted from readAllXXXX along dimension using morton index.
        
        :param vec: the vector to sort
        :type vec: RealArray
        :param iLevel: the level at which we are
        :type iLevel: int
        :note: (only works for regular grids)
        :param bx: number of blocks along x
        :type bx: int
        :param by: number of blocks along y
        :type by: int
        
        :return: the sorted vector
        :rtype: RealArray
        """
    
    def mortonSort3d(self, vec: RealArray, iLevel: int, bx: int, by: int, bz: int) -> RealArray:
        """
        Sorts an array in 2D extracted from readAllXXXX along dimension using morton index.
        
        :param vec: the vector to sort
        :type vec: RealArray
        :param iLevel: the level at which we are
        :type iLevel: int
        :note: (only works for regular grids)
        :param bx: number of blocks along x
        :type bx: int
        :param by: number of blocks along y
        :type by: int
        :param bz: number of blocks along z
        :type bz: int

        :return: the sorted vector
        :rtype: RealArray
        """
    
    @overload
    def getSortingMask2d(self, iLevel: int, bx: int, by: int) -> UInt64Array:
        """
         Returns the sequence of morton indices to get a regular grid sorted in memory

        :param iLevel: level at which the regular grid is build
        :type iLevel: int
        :param bx: block size along x
        :type bx: int
        :param by: block size along y
        :type by: int
        
        :return: The sorted vector of morton ids
        :rtype: UInt64Array
        """
    
    @overload
    def getSortingMask2d(self, iLevel: int, bx: int, by: int,
                         coarse_res_x: int, coarse_res_y: int) -> UInt64Array:
        """
         Returns the sequence of morton indices to get a regular grid sorted in memory

        :param iLevel: level at which the regular grid is build
        :type iLevel: int
        :param bx: block size along x
        :type bx: int
        :param by: block size along y
        :type by: int       
        :param coarse_res_x: only extracts the first coarse_res_x octants along the X direction
        :type coarse_res_x: int
        :param coarse_res_y: only extracts the first coarse_res_y octants along the Y direction
        :type coarse_res_y: int

        :return: The sorted vector of morton ids
        :rtype: UInt64Array
        """
    
    @overload
    def getSortingMask3d(self, iLevel: int, bx: int, by: int, bz: int) -> UInt64Array:
        """
         Returns the sequence of morton indices to get a regular grid sorted in memory

        :param iLevel: level at which the regular grid is build
        :type iLevel: int
        :param bx: block size along x
        :type bx: int
        :param by: block size along y
        :type by: int
        :param bz: block size along z
        :type bz: int

        :return: The sorted vector of morton ids
        :rtype: UInt64Array
        """
    
    @overload
    def getSortingMask3d(self, iLevel: int, bx: int, by: int, bz: int,
                         coarse_res_x: int, coarse_res_y: int, coarse_res_z: int) -> UInt64Array:
        """
         Returns the sequence of morton indices to get a regular grid sorted in memory

        :param iLevel: level at which the regular grid is build
        :type iLevel: int
        :param bx: block size along x
        :type bx: int
        :param by: block size along y
        :type by: int
        :param bz: block size along z
        :type bz: int
        :param coarse_res_x: only extracts the first coarse_res_x octants along the X direction
        :type coarse_res_x: int
        :param coarse_res_y: only extracts the first coarse_res_y octants along the Y direction
        :type coarse_res_y: int
        :param coarse_res_z: only extracts the first coarse_res_z octants along the Z direction
        :type coarse_res_z: int

        :return: The sorted vector of morton ids
        :rtype: UInt64Array
        """
    
    def fillLine(self, line: Line) -> None:
        """
        Extracts quantities along a line in the dataset
        
        :param line: the object to fill in.
        :note: Note that Nl, start and end should already be filled in.
        :type line: Line
        """

    def fillLineUnique(self, line: Line) -> None:
        """
        Extracts quantities along a line and removes duplicate ids
        
        :param line: the object to fill in.
        :note: Note that Nl, start and end should already be filled in.
        :type line: Line
        """


class XdmfReader:
    """Class to read XDMF files."""
    def readSnapshot(self, filename: str) -> Snapshot:
        """
        Reads a snapshot from an xdmf file
        
        :param filename: the path to the file to read
        :type filename: str

        :return: A Snapshot object
        :rtype: Snapshot
        """
    
    def readTimeSeries(self, filename: str) -> list[Snapshot]:
        """
        Reads a time-series from an xdmf file

        :param filename: the path to the file to read
        :type filename: str

        :return: A TimeSeries object, a list of Snapshot.
        :rtype: list[Snapshot]
        """
