// Copyright (c) 2012-2014 Shuntaro Yamazaki (shun-yamazaki (at) aist.go.jp)
// and the National Institute of Advanced Industrial Science and Technology

#pragma once

#include "constant.h"

#include "mesh.h"

namespace hbm {

/// import a landmark file in map format
/**
@verbatim
 num_landmarks
 landmark_id vertex_id
 :
@endverbatim
*/
MHBMLIB_API void ImportLandmarkMap(
    const std::string& filename, ///< file
    const CMesh& mesh, ///< source mesh
    const slib::CSparseMatrix<double>& matAllD, ///< mx3 matrix of target 3D coordinates
    slib::CSparseMatrix<double>& C, ///< source landmark matrix
    slib::CSparseMatrix<double>& D ///< target landmark matrix
);

/// import landmarks from a file in (xyz) coordinate format
/**
@verbatim
 x y z
 :
@endverbatim
*/
MHBMLIB_API slib::CSparseMatrix<double> ImportLandmarkCoordinate(
    const std::string& filename, ///< file
    const CMesh& mesh  ///< polygon mesh
);

/// export landmarks to a file.
/**
@verbatim
 x y z
 :
@endverbatim
*/
MHBMLIB_API  void ExportLandmarkCoordinate(
    const std::string& filename, ///< file
    const CMesh& mesh, ///< polygon mesh
    const slib::CSparseMatrix<double>& landmark ///< landmark matrix
);

} // namespace hbm
