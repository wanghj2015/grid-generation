# Grid Generation for Curvilinear Field-Line-Following Coordinates

Fortran code for generating curvilinear magnetic field-line-following coordinate systems for ionosphere-plasmasphere modeling.

## Paper

This code accompanies the paper:

> **A general curvilinear magnetic field-line-following coordinate system for ionosphere-plasmasphere modeling**  
> Houjun Wang  
> CIRES, University of Colorado Boulder

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5768675.svg)](https://doi.org/10.5281/zenodo.5768675)

## Contents

- `dipole/` - Fortran 90 modules for dipole coordinate generation
- `igrf/` - IGRF (International Geomagnetic Reference Field) coefficients
- `plot/` - Visualization scripts
- `scripts/` - Utility scripts

## Building

```bash
cd dipole
make
```

## Usage

See the dipole coordinate modules for generating field-aligned grids for ionosphere-plasmasphere simulations.
