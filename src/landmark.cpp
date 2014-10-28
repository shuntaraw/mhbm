// Copyright (c) 2012-2014 Shuntaro Yamazaki (shun-yamazaki (at) aist.go.jp)
// and the National Institute of Advanced Industrial Science and Technology

#include "constant.h"

#include "landmark.h"
#include "correspondence.h"
#include "deformation.h" // transform_coordinate_matrix

/// non-rigid mesh registration
namespace hbm {

MHBMLIB_API void ImportLandmarkMap(
    const std::string& filename,
    const CMesh& mesh,
    const slib::CSparseMatrix<double>& matAllD,
    slib::CSparseMatrix<double>& matC,
    slib::CSparseMatrix<double>& matD
) {
    std::clog << "landmark <= " << filename << std::endl;
    if (!matAllD.num_rows()) {
        ThrowRuntimeError("scan landmark not loaded");
    }
    std::ifstream ifs(filename);
    if (!ifs.is_open()) {
        return;
    }
    int nlandmarks;
    ifs >> nlandmarks;
    slib::MatrixGenerator<double> genC, genD;
    for (int i = 0; i < nlandmarks; i++) {
        int lid, vid;
        ifs >> lid >> vid;
        if (ifs.fail()) {
            ThrowRuntimeError("premature end of file: " + filename);
        }
        genC.Add(i, vid - 1, 1);
        genD.Add(i, lid - 1, 1);
    }
    matC = genC.GenerateSparse(nlandmarks, mesh.vertices.size());
    matD = genD.GenerateSparse(nlandmarks, matAllD.num_rows()).MultiplyTo('N', matAllD);
    std::clog << matC.num_rows() << " landmarks"  << std::endl;
}

MHBMLIB_API slib::CSparseMatrix<double> ImportLandmarkCoordinate(
    const std::string& filename,
    const CMesh& mesh
) {
    std::clog << "landmark <= " << filename << std::endl;
    std::ifstream ifs(filename);
    if (!ifs.is_open()) {
        ThrowRuntimeError("failed to open " + filename);
    }
    std::vector<CMesh::CVertex> vertices;
    for (int i = 0;; i++) {
        std::string buf;
        getline(ifs, buf);
        if (ifs.fail()) {
            break;
        }
        std::istringstream str(buf);
        CMesh::CVertex v;
        str >> v.position[0] >> v.position[1] >> v.position[2];
        if (str.fail()) {
            ThrowRuntimeError("invalid format: " + filename);
        }
        vertices.push_back(std::move(v)); // landmark coordinate
    }
    CMesh lmpoints;
    lmpoints.vertices = std::move(vertices);
    // nearest neighbor to the mesh model
    MultidirectionalClosestPointSearch closest_point;///<closest point
    closest_point.SetParameters(&lmpoints, &mesh, std::numeric_limits<float>::max(), -1, true, SEARCH_DIRECTION::FORWARD);
    MatrixPair mat = closest_point.Find().ToPositionMatrix();
    std::clog << lmpoints.vertices.size() << " landmarks"  << std::endl;
    return mat.C.MultiplyTo('T', mat.D); ///< reset order
}

MHBMLIB_API void ExportLandmarkCoordinate(
    const std::string& filename,
    const CMesh& mesh,
    const slib::CSparseMatrix<double>& landmarks
) {
    std::clog << "landmark => " << filename << std::endl;
    if (!landmarks.num_rows()) {
        ThrowRuntimeError("no landmark");
    }
    slib::CMatrix<double> v;
    mesh.ConvertMeshToCoordinate(v);
    v = landmarks.MultiplyTo('N', v);
    std::ofstream ofs(filename);
    for (int r = 0; r < v.num_rows(); r++) {
        ofs << v(r, 0) << " " << v(r, 1) << " " << v(r, 2) << "\n";
    }
}

} // namespace hbm