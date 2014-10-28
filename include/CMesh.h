// Copyright (c) 2012-2014 Shuntaro Yamazaki (shun-yamazaki (at) aist.go.jp)
// and the National Institute of Advanced Industrial Science and Technology

#pragma once

#include "MeshBase.h"
#include "mkl_csr_driver.h"

namespace hbm {

/// vertex of a mesh
struct CCustomVertex :
        public slib::mesh::CVertexPosition,
        public slib::mesh::CVertexColor,
        public slib::mesh::CVertexNormal {
    float elasticity = 1; ///< elasticity
    float area = 1; ///<  voronoi area
};

/// face of a mesh
struct CCustomFace :
    public slib::mesh::CTriangleIndex {};

/// triangular mesh
class CMesh : public slib::mesh::CPolygonMesh<CCustomVertex, CCustomFace> {
public:
    CMesh() = default;
    CMesh(const CMesh&) = default;
    CMesh& operator=(const CMesh& m) = default;
    CMesh(CMesh&& m) : slib::mesh::CPolygonMesh<CCustomVertex, CCustomFace>(std::move(m)) {}
    CMesh& operator=(CMesh && m) {
        slib::mesh::CPolygonMesh<CCustomVertex, CCustomFace>::operator=(std::move(m));
        return *this;
    }

    /// update vertex property
    void UpdateGeometry();

    void Dump()const;
    void Read(const std::string&filename);
    void Write(const std::string&filename)const;

    /// compute a laplacian matrix to generate differential coordinate
    /// @return laplacian matrix
    slib::CSparseMatrix<double> GenerateLaplacianMatrix( const slib::CMatrix<double>& org_pos)const;

    void ConvertMeshToCoordinate( slib::CMatrix<double>& x )const ;
    void ConvertCoordinateToMesh( const slib::CMatrix<double>& x);

    void AffineTransform(const slib::CMatrix<float, 4, 4>& mat);

    std::vector<int> SelectDuplicatedFaces()const;
    std::vector<int> SelectNonManifoldVertices( )const;
    std::vector<int> SelectCollapsedFaces( )const;
    std::vector<int> SelectUnusedVertices()const;
    void DeleteVertices(
        const std::vector<int>& vid_to_delete///< list of face ids.
        );
    void DeleteFaces(
        const std::vector<int>& fid_to_delete ///< list of face ids.
        ) ;

    CMesh SubdivideModifiedButterfly() const;
    CMesh SubdivideLoop() const;
};


///// affine transform 3D coordinate matrix
///**
//@verbatim
//[x,y,z] = [R(x,y,z)+(tx,ty,tz)]
//[  :  ]   [        :          ]
//@endverbatim
//*/
//void AffineTransformCoordinateMatrix(
//    const slib::CMatrix<double, 4, 4>& trans, // affine transformation
//    slib::CMatrix<double>& coordinates // mx3 coordinate matrix
//    );

///// reorder 3D coordinate matrix from mx3 to 3mx1
///**
//@verbatim
//[x y z] -> [x]
//[: : :]    [y]
//[     ]    [z]
//[     ]    [:]
//@endverbatim
//*/
///// @return 3mx1 coordinate matrix
//slib::CMatrix<double> VectorizeCoordinates(const slib::CMatrix<double>& mat);

///// reorder 3D coordinate matrix from mx3 to 3mx1
///**
//@verbatim
//[x] -> [x y z]
//[y]    [  :  ]
//[z]
//[:]
//@endverbatim
//*/
///// @return 3mx1 coordinate matrix
//slib::CMatrix<double> StackCoordinates(const slib::CMatrix<double>& mat);

///// expand a matrix from nxn to 3nx3n block diagonal
///**
//@verbatim
//[e  ]   [e   ]
//[ \ ] = [ e  ]
//[   ]   [  e ]
//[   \]
//@endverbatim
//*/
///// @return 3nx3n block diagonal matrix
//slib::CSparseMatrix<double> ExpandLaplacianMatrix(const slib::CSparseMatrix<double>& mat);

}// namespace hbm
