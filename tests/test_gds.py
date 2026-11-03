"""Create one rectangle for a manual GDSII import test in Sonnet."""

from pathlib import Path
import sys

import gdstk


def main():
    # Coordinates are in microns; the storage precision is 1 nanometer.
    library = gdstk.Library(unit=1e-6, precision=1e-9)
    cell = library.new_cell("TOP")

    # A single 200 x 100 um metal rectangle, from (0, 0) to (200, 100).
    rectangle = gdstk.rectangle((0, 0), (200, 100), layer=1, datatype=0)
    cell.add(rectangle)

    # Save beside this script, regardless of the terminal's working directory.
    output_path = Path(__file__).resolve().parent / "test_sonnet_import.gds"
    library.write_gds(str(output_path))

    print(f"Generated: {output_path}")
    print(f"Cell: {cell.name}")
    print(f"Layer: {rectangle.layer}; datatype: {rectangle.datatype}")
    print("Rectangle: width = 200 um; height = 100 um")
    print(f"Bounding box (um): {rectangle.bounding_box()}")
    print("Units: 1 um; precision: 1 nm")


if __name__ == "__main__":
    main()
