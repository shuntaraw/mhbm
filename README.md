mHBM: Markerless Homologous Body Modeling Software
====

Markerless Homologous Body Modeling (mHBM) software is developed for creating a 3D triangular mesh
model which approximates the 3D shape of a point cloud obtained by 3D scanning (scan data) by registering a
template mesh model (generic model) onto the scan data by non-rigid mesh deformation. The registered mesh
is the homologous model which has the same vertex connectivity as the generic model and represents geometric
features of the scan data.

This software is a successor of Homologous Body Modeling (HBM) software developed by our research
group, carrying several advantages over the previous software, including registration capability without predetermined
correspondences (landmarks) by an iterative closest point (ICP) technique, accurate registration
by non-rigid mesh deformation, small mesh distortion by global numerical optimization, and fast computation
using hardware acceleration.

The software required to build the mHBM software from source code is listed below.
- Microsoft Visual Studio Express 2013 for Windows Desktop Update 3 
- Intel Math Kernel Library 11.0 
- Boost C++ Template Library 1.54.0 