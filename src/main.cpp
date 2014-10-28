// Copyright (c) 2012-2014 Shuntaro Yamazaki (shun-yamazaki (at) aist.go.jp)
// and the National Institute of Advanced Industrial Science and Technology

#include "constant.h"

#include <omp.h> // omp_set_num_threads()
#include <boost/filesystem.hpp>

#include "deformation.h"
#include "mesh.h"
#include "landmark.h"

/// entry point
class Main {
public:
    /// parse command line options
    void ParseCommandLine(int argc, char *argv[]) {
        if (argc == 7 && argv[1][0] != '-') {            // HBM-compatible format
            ParseCommandLineHBM(argv);
        } else { // Command-line swtich
            ParseCommandLineMHBM(argc, argv);
        }
    }

    /// processes a single line in a param file
    bool ParseParam() {
        char buf[BUFSIZ];
        std::istringstream str;
        int process;
        while (!param_file_.fail()) {
            param_file_.getline(buf, BUFSIZ);
            if (strrchr(buf, '#')) {
                *strrchr(buf, '#') = 0;
            }
            str.str(buf);
            str.clear();
            str >> process;
            if (!str.fail()) {
                break; // found non-empty line
            }
        }
        if (param_file_.fail()) {
            return false; // EOF
        }
        switch (process) {
        case 0: {
            // default
            char model = 'r';
            char direction = 'f';
            // read
            str >> model >> direction  ;
            std::cout << "registration: " << std::endl;
            ProcessRegistration(model, direction);
        }
        break;
        case 1: {
            // default
            float correspondence_weight = 1;
            float shape_weight = 1;
            float landmark_weight = 0;
            char model = 't';
            char direction = 'f';
            char point2point = 'p';
            float max_distance = 100;
            float max_angle = 30;
            //float rigidness_weight = 0;
            // read
            str >> correspondence_weight >> shape_weight >> landmark_weight >> model >> direction >> point2point >> max_distance >> max_angle /*>> rigidness_weight*/;
            // precompute
            correspondence_weight /= shape_weight;
            landmark_weight /= shape_weight;
            std::cout << "fitting: " << std::endl;
            ProcessDeformation(correspondence_weight, shape_weight, landmark_weight, model, direction, point2point, max_distance * max_distance, cos(M_PI / 180 * max_angle)/*, rigidness_weight*/);
        }
        break;
        case 2:
            std::cout << "approximating subdivision: " << std::endl;
            deformable_.SubdivideApproximatingSubdivision();
            break;
        case 3:
            std::cout << "interpolating subdivision: " << std::endl;
            deformable_.SubdivideInterpolatingSubdivision();
            break;
        case 4:  {
            static int dumpid = 0; // HACK: to define a local variable
            std::string filename = output_filename_;
            filename.insert(filename.find_last_of('.'), "-dump-" + std::to_string(dumpid++));    // default
            // read
            str >> filename;
            //std::cout << "dump: " << filename << std::endl;
            deformable_.mesh(). Write(filename);
        }
        break;
        case 5:
            // rebase
            std::cout << "setting current mesh as base: " << std::endl;
            deformable_.LoadMesh(deformable_.mesh());
            break;
        default:
            ThrowLogicError("undefined label");
        }
        return true;
    }

    void ExportResult() const {
        if (!output_filename_.empty()) {
            ExportModelShape(output_filename_);
        }
        if (!output_lm_filename_.empty()) {
            ExportModelLandmark(output_lm_filename_);
        }
    }

private:
    void PrintVersion() const {
        std::cout << hbm::VERSION << std::endl;
        exit(1);
    }

    void PrintUsage() const {
        std::cout << "Usage: mHBM <param> <generic> <out> <generic_landmark> <scan_landmark> <scan>" << std::endl
                  << "-or-" << std::endl
                  << "Usage: mHBM [options]" << std::endl
                  << "file options:" << std::endl
                  << "-e <model_lm>: [out] landmark of <model>" << std::endl
                  << "-g <generic>: [out] generic model" << std::endl
                  << "-l <scan_lm>: [in] landmark of <scan>" << std::endl
                  << "-m <generic_lm>: [in] landmark of <generic>" << std::endl
                  << "-p <param>: [in] parameter file" << std::endl
                  << "-r <model>: [out] deformation of <generic>" << std::endl
                  << "-w <generic_elasticity>: [in] elasticity of <generic>" << std::endl
                  << "-s <scan>: [in] scan data" << std::endl
                  << "-t <lm_pairs>: [in] pairs of valid landmark IDs" << std::endl
                  << "control options:" << std::endl
                  << "-n: disable multi-threading" << std::endl
                  << "-h: print help" << std::endl
                  << "-v: print version" << std::endl
                  ;
        exit(1);
    }

    // mHBM <param> <generic> <out> <generic_landmark> <scan_landmark> <scan>
    void ParseCommandLineHBM(char *argv[])  {
        LoadParam(argv[1]);
        LoadGenericShape(argv[2]);
        output_filename_ = argv[3];
        LoadScanShape(argv[6]);
        std::vector<std::pair<int, int>> useall;
        if (boost::filesystem::exists(argv[5])) {
            // after scandata is loaded
            LoadScanLandmark(argv[5], useall);
        }
        if (boost::filesystem::exists(argv[4])) {
            // after generic (and scan landmark) is loaded
            LoadGenericLandmark(argv[4], useall);
        }
    }

    // mHBM [options]
    void ParseCommandLineMHBM(int argc,  char *argv[]) {
        std::string generic_lm;
        std::string scan_lm;
        //std::string initshape;
        std::string elasticity;
        std::vector<std::pair<int, int>> lm_pairs;
        for (int i = 1; i < argc; i++) {
            if (argv[i][0] != '-' || strlen(argv[i]) < 2) {
                std::cerr << "invalid commandline option" << std::endl;
                PrintUsage();
            }
            switch (argv[i][1]) {
            case 'e':
                output_lm_filename_ = argv[++i];
                break;
            case 'g':
                LoadGenericShape(argv[++i]);
                break;
            //case 'i':
            //    initshape = argv[++i];
            //    break;
            case 'h':
                PrintUsage();
                break;
            case 'l':
                scan_lm = argv[++i];
                break;
            case 'm':
                generic_lm = argv[++i];
                break;
            case 'n':
                std::cout << "+ multi-threading disabled" << std::endl;
                omp_set_num_threads(1);
                mkl_set_num_threads(1);
                break;
            case 'p':
                LoadParam(argv[++i]);
                break;
            case 'r':
                output_filename_ = argv[++i];
                break;
            case 's':
                LoadScanShape(argv[++i]);
                break;
            case 't':
                LoadValidLandmarkID(argv[++i], lm_pairs);
                break;
                //case 'u':
                //    LoadUnusedLandmark(argv[++i], unused_scanlm);
                break;
            case 'v':
                PrintVersion();
                break;
            case 'w':
                elasticity = argv[++i];
                break;
            default:
                std::cerr << "invalid commandlind option" << std::endl;
                PrintUsage();
            }
        }
        if (!param_file_.is_open()) {
            std::cerr << "parameter file not found" << std::endl;
            PrintUsage();
        }
        if (!scan_lm.empty()) {
            // after scandata is loaded
            LoadScanLandmark(scan_lm, lm_pairs);
        }
        if (!generic_lm.empty()) {
            // after generic and (scan landmark) is loaded
            LoadGenericLandmark(generic_lm, lm_pairs);
        }
        if (!elasticity.empty()) {
            // after generic is loaded
            LoadGenericElasticity(elasticity);
        }
        //if (!initshape.empty()) {
        //    // after generic is loaded
        //    LoadGenericInitShape(initshape);
        //}
    }

    /// parses a parameter file
    void LoadParam(const std::string& filename) {
        if (filename.empty()) {
            return;
        }
        param_file_.open(filename);
        if (!param_file_.is_open()) {
            ThrowRuntimeError("failed to open: " + filename);
        }
        char buf[BUFSIZ];
        while (!param_file_.fail()) {
            param_file_.getline(buf, BUFSIZ);
            if (strrchr(buf, '#')) {
                *strrchr(buf, '#') = 0;
            }
            std::istringstream str(buf);
            int inputformat, outputformat; // unused
            float fit_parameter; // unused
            str >> inputformat >> outputformat >> fit_parameter;
            if (!str.fail()) {
                return;
            }
        }
        ThrowRuntimeError("premature EOF: " + filename);
    }

    /// import generic model file
    void LoadGenericShape(const std::string& filename) {
        hbm::CMesh generic;
        generic.Read(filename);
        if (generic.faces.empty()) {
            ThrowRuntimeError("empty generic mesh");
        }
        std::vector<int> selected;
        selected = generic.SelectDuplicatedFaces();
        if (!selected.empty()) {
            ThrowRuntimeError("generic mesh has duplicated face");
        }
        selected = generic.SelectNonManifoldVertices();
        if (!selected.empty()) {
            ThrowRuntimeError("generic mesh has unorientable edge");
        }
        selected = generic.SelectCollapsedFaces();
        if (!selected.empty()) {
            ThrowRuntimeError("generic mesh has collapsed face");
        }
        selected = generic.SelectUnusedVertices();
        if (!selected.empty()) {
            ThrowRuntimeError("generic mesh has unused vertex");
        }
        deformable_.LoadMesh(std::move(generic));
    }

    /// import scan data file
    void LoadScanShape(const std::string& filename) {
        scan_mesh_.Read(filename);
        if (!scan_mesh_.faces.empty()) { // triangle mesh
            std::vector<int> selected;
            selected = scan_mesh_.SelectDuplicatedFaces();
            if (!selected.empty()) {
                std::clog << "warning: scandata has duplicated faces." << std::endl;
                scan_mesh_.DeleteFaces(selected);
            }
            selected = scan_mesh_.SelectNonManifoldVertices();
            if (!selected.empty()) {
                std::clog << "warning: scandata has unorientable edges. " << std::endl;
                scan_mesh_.DeleteFaces(selected);
            }
            selected = scan_mesh_.SelectCollapsedFaces();
            if (!selected.empty()) {
                std::clog << "warning: scandata has collapsed faces. " << std::endl;
                scan_mesh_.DeleteFaces(selected);
            }
            selected = scan_mesh_.SelectUnusedVertices();
            if (!selected.empty()) {
                std::clog << "warning: scandata has unused vertex." << std::endl;
                scan_mesh_.DeleteVertices(selected);
            }
            scan_mesh_.UpdateGeometry();
        } else { // point cloud
            for (auto& v : scan_mesh_.vertices) {
                v.normal = slib::CVector <float, 3> {1, 0, 0}; // dummy non-zero vector
            }
        }
    }

    /// import generic model landmarks
    void LoadGenericLandmark(const std::string& generic_lmfile, const std::vector<std::pair<int, int>>& pairs) {
        if (deformable_.mesh().vertices.empty()) {
            ThrowRuntimeError("generic model is not loaded");
        }
        slib::CSparseMatrix<double> C;
        boost::filesystem::path p = generic_lmfile;
        if (p.extension() == ".map") {
            ImportLandmarkMap(generic_lmfile, deformable_.mesh(), scan_landmark_, C, scan_landmark_);
        } else if (p.extension() == ".smat") {
            C.Read(generic_lmfile);
        } else  {
            C = ImportLandmarkCoordinate(generic_lmfile, deformable_.mesh());
        }
        if (!pairs.empty()) {
            slib::CSparseMatrix<double> P;
            P.Resize(pairs.size(), C.num_rows());
            for (int i = 0; i < pairs.size(); i++) {
                P.Add(i, pairs[i].first, 1);
            }
            C = P.MultiplyTo('N', C);
        }
        std::clog << "using " << C.num_rows() << " landmarks " << std::endl;
        deformable_.set_landmark(std::move(C));
    }

    /// import scan landmark file
    void LoadScanLandmark(const std::string& scan_lmfile, const std::vector<std::pair<int, int>>& pairs) {
        if (scan_mesh_.vertices.empty()) {
            ThrowRuntimeError("scan data is not loaded");
        }
        scan_landmark_ = hbm::ImportLandmarkCoordinate(scan_lmfile, scan_mesh_);
        if (!pairs.empty()) {
            slib::CSparseMatrix<double> P;
            P.Resize(pairs.size(), scan_landmark_.num_rows());
            for (int i = 0; i < pairs.size(); i++) {
                P.Add(i, pairs[i].second, 1);
            }
            scan_landmark_ = P.MultiplyTo('N', scan_landmark_);
        }
        std::clog << "using " << scan_landmark_.num_rows() << " landmarks " << std::endl;
    }

    /// import generic model elasticity
    void LoadGenericElasticity(const std::string& filename) {
        if (deformable_.mesh().vertices.empty()) {
            ThrowRuntimeError("generic mesh not loaded");
        }
        std::ifstream ifs(filename);
        if (!ifs.is_open()) {
            ThrowRuntimeError("failed to open " + filename);
        }
        auto vertices = deformable_.vertex_range();
        for (auto& vertex : vertices) {
            vertex.elasticity = 1;
        }
        while (true) {
            int v;
            float w;
            ifs >> v >> w;
            if (ifs.fail()) {
                break;
            }
            if (v < 0 || v >= vertices.size() || w < 0 || w > 1) {
                ThrowRuntimeError("invalid weight at %d: elasticity = %g", v, w);
            }
            vertices [v].elasticity = w;
        }
    }

    ///// import initial deformation
    //void LoadGenericInitShape(const std::string& filename) {
    //    hbm::CMesh deformed;
    //    Read(deformed, filename);
    //    if (deformed.vertices.size() != deformable_.mesh().vertices.size()) {
    //        ThrowRuntimeError("initshape incompatible");
    //    }
    //    auto vertices = deformable_.vertex_range();
    //    for (int vid = 0; vid < vertices.size(); vid++) {
    //        vertices[vid] = deformed.vertices[vid];
    //    }
    //}

    /// export generic model landmark
    void ExportModelLandmark(const std::string& filename) const {
        if (boost::filesystem::path(filename).extension() == ".smat") {
            std::cout << "+ matrix format" << std::endl;
            deformable_.landmark().Write(filename);
        } else {
            std::cout << "+ coordinate format" << std::endl;
            ExportLandmarkCoordinate(filename, deformable_.mesh(), deformable_.landmark());
        }
        if (deformable_.landmark().num_rows() == scan_landmark_.num_rows()) {
            int nlandmarks = deformable_.landmark().num_rows();
            slib::CMatrix<double> g, s;
            deformable_.mesh().ConvertMeshToCoordinate(g);
            g = deformable_.landmark().MultiplyTo('N', g);
            scan_mesh_.ConvertMeshToCoordinate(s);
            s = scan_landmark_.MultiplyTo('N', s);
            g -= s;
            int worst_id;
            float worst_d2 = 0, average_d2 = 0;
            for (int r = 0; r < nlandmarks; r++) {
                float d2 = norm2_squared_of(make_vector_from_row(g, r));
                if (worst_d2 < d2) {
                    worst_d2 = d2;
                    worst_id = r;
                }
                average_d2 += sqrt(d2);
            }
            average_d2 /= nlandmarks;
            std::clog << "landmark distance:" << std::endl
                      << "average = " << average_d2 << ", worst = " << sqrt(worst_d2) << " (#" << worst_id << ")" << std::endl;
        }
    }

    /// export generic model shape
    void ExportModelShape(const std::string& filename) const {
        deformable_.mesh().Write(filename);
    }

    /// perform rigid registration
    void ProcessRegistration(
        char model, // {r, s, t, x, y, z}
        char direction // {f, b, m}
    ) {
        switch (model) {
        case 'r':
        case 's':
        case 't': {
            hbm::TRANSFORMATION m;
            switch (model) {
            case 'r':
                std::cout << "+ rigid transformation" << std::endl;
                m = hbm::TRANSFORMATION::RIGID;
                break;
            case 's':
                std::cout << "+ similarity " << std::endl;
                m = hbm::TRANSFORMATION::SIMILARITY;
                break;
            case 't':
                std::cout << "+ translation " << std::endl;
                m = hbm::TRANSFORMATION::TRANSLATION;
                break;
            default:
                ThrowLogicError("undefined model");
            }
            slib::CMatrix<float, 4, 4> mat;
            if (deformable_.landmark().num_rows() > 2) {
                bool enable_icp = false;
                mat = EstimateAffine(
                          deformable_.mesh(),
                          scan_mesh_,
                          &deformable_.landmark(),
                          &scan_landmark_,
                          1,
                          m,
                          0, 0, false, hbm::SEARCH_DIRECTION::FORWARD,
                          enable_icp);
            } else {
                float max_distance2 = std::numeric_limits<float>::max();
                float min_cosangle = -1;
                bool allow_border = false;
                hbm::SEARCH_DIRECTION d;
                switch (direction) {
                case 'b':
                    std::cout << "+ backward " << std::endl;
                    d = hbm::SEARCH_DIRECTION::BACKWARD;
                    break;
                case 'f':
                    std::cout << "+ forward " << std::endl;
                    d = hbm::SEARCH_DIRECTION::FORWARD;
                    break;
                case 'm':
                    std::cout << "+ bidirectional " << std::endl;
                    d = hbm::SEARCH_DIRECTION::BIDIRECTIONAL;
                    break;
                default:
                    ThrowLogicError("undefined direction");
                }
                bool enable_icp = true;
                mat = EstimateAffine(
                          deformable_.mesh(),
                          scan_mesh_,
                          0,
                          0,
                          0,
                          m,
                          max_distance2,
                          min_cosangle,
                          allow_border,
                          d,
                          enable_icp);
            }
            deformable_.ApplyTransformation(mat);
        }
        break;
        case 'x':
            std::cout << "rotation around x axis" << std::endl;
            deformable_.ApplyTransformation(slib::make_rotation_matrix(slib::CVector <float, 3> {1, 0, 0}, float(M_PI / 2)));
            break;
        case 'y':
            std::cout << "rotation around y axis" << std::endl;
            deformable_.ApplyTransformation(slib::make_rotation_matrix(slib::CVector <float, 3> {0, 1, 0}, float(M_PI / 2)));
            break;
        case 'z':
            std::cout << "rotation around z axis" << std::endl;
            deformable_.ApplyTransformation(slib::make_rotation_matrix(slib::CVector <float, 3> {0, 0, 1}, float(M_PI / 2)));
            break;
        default:
            ThrowLogicError("undefined label");
        }
    }

    // perform fitting by non-rigid deformation
    void ProcessDeformation(
        float correspondence_weight, // correspondence weight
        float shape_weight,
        float landmark_weight, // landmark weight
        char model,// {t,r,s}
        char direction, // value in {f,b}
        char point2point,//{p,s}
        float max_distance2, // maximum squared distance between corresponding points
        float min_cosangle  // minimum cos angle between corresponding points
        //float rigidness_weight // rigidness weight
    ) {
        assert(shape_weight > 0);
        correspondence_weight /= shape_weight;
        landmark_weight /= shape_weight;
        std::cout << "+ correspondence = " << correspondence_weight << std::endl;
        if (landmark_weight) {
            std::cout << "+ landmark = " << landmark_weight << std::endl;
        }

        hbm::TRANSFORMATION deformation_model;
        switch (model) {
        case 't':
            deformation_model = hbm::TRANSFORMATION::TRANSLATION;
            std::cout << "+ local translation" << std::endl;
            break;
        case 'r':
            deformation_model = hbm::TRANSFORMATION::RIGID;
            std::cout << "+ as-rigid-as-possible" << std::endl;
            break;
        case 's':
            deformation_model = hbm::TRANSFORMATION::SIMILARITY;
            std::cout << "+ as-conformal-as-possible" << std::endl;
            break;
        default:
            ThrowLogicError("undefined model");
        }

        hbm::SEARCH_DIRECTION search_direction;
        switch (direction) {
        case 'f':
            search_direction = hbm::SEARCH_DIRECTION::FORWARD;
            std::cout << "+ forward search" << std::endl;
            break;
        case 'b':
            search_direction = hbm::SEARCH_DIRECTION::BACKWARD;
            std::cout << "+ backward search" << std::endl;
            break;
        case 'm':
            search_direction = hbm::SEARCH_DIRECTION::BIDIRECTIONAL;
            std::cout << "+ bidirectional search" << std::endl;
            break;
        default:
            ThrowLogicError("undefined direction");
            break;
        }

        bool use_point_plane_distance;
        switch (point2point) {
        case 'p':
            use_point_plane_distance = false;
            break;
        case 's':
            switch (search_direction) {
            case hbm::SEARCH_DIRECTION::FORWARD:
                use_point_plane_distance = !scan_mesh_.faces.empty();
                break;
            case hbm::SEARCH_DIRECTION::BACKWARD:
                use_point_plane_distance = !deformable_.mesh().faces.empty();
                break;
            case hbm::SEARCH_DIRECTION::BIDIRECTIONAL:
                use_point_plane_distance = (!deformable_.mesh().faces.empty() && !scan_mesh_.faces.empty());
                break;
            }
            break;
        default:
            ThrowLogicError("undefined point2point");
            break;
        }
        if (use_point_plane_distance) {
            std::cout << "+ point-plane distance" << std::endl;
        } else {
            std::cout << "+ point-point distance" << std::endl;
        }

        std::cout << "+ max_distance = " << sqrt(max_distance2) << std::endl;
        std::cout << "+ max_angle = " << acos(min_cosangle) / M_PI * 180 << std::endl;
        //std::cout << "+ rigidness = " << rigidness_weight << std::endl;

        bool allow_border = false;
        auto pos = DeformNonrigid(
                       deformable_.mesh(),
                       scan_mesh_,
                       &deformable_.landmark(),
                       &scan_landmark_,
                       deformable_.org_pos(),
                       correspondence_weight,
                       landmark_weight,
                       //rigidness_weight,
                       use_point_plane_distance,
                       deformation_model,
                       max_distance2,
                       scan_mesh_.faces.size() ? min_cosangle : -1,
                       allow_border,
                       search_direction);
        deformable_.SetCoordinates(pos);
    }

    void LoadValidLandmarkID(const std::string& filename, std::vector<std::pair<int, int>>& pairs) {
        std::ifstream ifs(filename);
        if (!ifs.is_open()) {
            ThrowRuntimeError("failed to open " + filename);
        }
        pairs.clear();
        while (1) {
            std::pair<int, int> p;
            ifs >> p.first >> p.second;
            if (ifs.fail()) {
                break;
            }
            pairs.push_back(p);
        }
        //for (auto&p : pairs) {
        //    std::clog << p.first << "," << p.second<<std::endl;
        //}
    }
    //void LoadUnusedLandmark(const std::string& filename, std::vector<int>& unusedlm) {
    //    std::ifstream ifs(filename);
    //    if (!ifs.is_open()) {
    //        ThrowRuntimeError("failed to open " + filename);
    //    }
    //    unusedlm.clear();
    //    while (1) {
    //        int i;
    //        ifs >> i;
    //        if (ifs.fail()) {
    //            break;
    //        }
    //        unusedlm.push_back(i);
    //    }
    //}

    std::ifstream param_file_; // file stream
    hbm::DeformableMesh deformable_; // source mesh
    hbm::CMesh scan_mesh_; // destination mesh
    slib::CSparseMatrix<double> scan_landmark_; // scan landmarks
    std::string output_filename_; // output model filename
    std::string output_lm_filename_; // output landmark filename (optional)
};

/// the entry point
int main(int argc, char *argv[]) {
    try {
#ifdef NDEBUG
        std::clog.setstate(std::ios::failbit);
#endif
        Main m;
        m.ParseCommandLine(argc, argv);
        while (m.ParseParam()) {}
        m.ExportResult();
    } catch (std::exception& e) {
        std::cerr << "exception: " << e.what() << std::endl;
    }
    return 0;
}

