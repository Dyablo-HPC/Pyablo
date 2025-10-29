""" Pyablo Module"

.. curretmodule:: pyablo
.. autosummary::
   :toctree: _autosummary

   Line
   Snapshot
   XdmfReader
"""


class Snapshot:
    """Class representing a snapshot of simulation data."""

    def close(self) -> None:
        """Closes the snapshot and releases resources."""
    
    def setName(self, name: str) -> None:
        """Sets the name of the snapshot."""
    
    def setNDim(self, ndim: int) -> None:
        """Sets the number of dimensions of the snapshot."""
    
    def print(self) -> None:
        """Prints a summary of the snapshot."""
    
    def getCellsFromPosition(self, position: tuple[float, float, float]) -> list[int]:
        """Gets cell indices from a given position."""

    def probeDensity(self, points: list[tuple[float, float, float]]) -> list[float]:
        """Probes densities at specified points."""
    

class Line:
    """Class representing a line in the simulation domain."""
    Nl: int
    start: tuple[float, float, float]
    end: tuple[float, float, float]
    pos: list[tuple[float, float, float]]
    rho: float
    prs: float
    E: float
    vel: tuple[float, float, float]
    cellIds: list[int]


class XdmfReader:
    """Class to read XDMF files."""
    def readSnapshot(self, filename: str) -> None:
        """Reads a snapshot from an XDMF file."""
    
    def readTimeSeries(self, filename: str) -> list[Snapshot]:
        """Reads a time series of snapshots from an XDMF file."""
