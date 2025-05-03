#ifdef WITH_CUDA
#include <cuda_runtime.h>
#endif

#include "config.hpp"
#include "field-math.hpp"
#include "optimizer.hpp"
#include "parametrizer.hpp"
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <nlohmann/json.hpp>

namespace fs = std::filesystem;
using namespace qflow;

Parametrizer field;

int main(int argc, char** argv) {
    setbuf(stdout, NULL);

#ifdef WITH_CUDA
    cudaFree(0);
#endif
    int t1, t2;
    std::string input_obj, output_obj, directory, filename, extension;
    int faces = -1;
    
    for (int i = 0; i < argc; ++i) {
        if (strcmp(argv[i], "-f") == 0) {
            sscanf(argv[i + 1], "%d", &faces);

        } else if (strcmp(argv[i], "-i") == 0) {
            input_obj = argv[i + 1];
            fs::path file_path(input_obj);
            directory = file_path.parent_path().string();
            filename = file_path.stem().string();
            extension = file_path.extension().string();
            
        } else if (strcmp(argv[i], "-o") == 0) {
            output_obj = argv[i + 1];
            
        } else if (strcmp(argv[i], "-sharp") == 0) {
            field.flag_preserve_sharp = 1;
        } else if (strcmp(argv[i], "-boundary") == 0) {
            field.flag_preserve_boundary = 1;
        } else if (strcmp(argv[i], "-adaptive") == 0) {
            field.flag_adaptive_scale = 1;
        } else if (strcmp(argv[i], "-mcf") == 0) {
            field.flag_minimum_cost_flow = 1;
        } else if (strcmp(argv[i], "-sat") == 0) {
            field.flag_aggresive_sat = 1;
        } else if (strcmp(argv[i], "-seed") == 0) {
            field.hierarchy.rng_seed = atoi(argv[i + 1]);
        }
    }

    output_obj = directory + "/" + filename + "_quad_quadriflow_" + std::to_string(faces) + extension;

    printf("%d %s %s\n", faces, input_obj.c_str(), output_obj.c_str());
    if (input_obj.size() >= 1) {
        field.Load(input_obj.c_str()); // Load mesh and normalize
        std::cout << "Extract vertices - " << field.V.rows() << " " << field.V.cols() << std::endl;
        std::cout << "Extract faces - " << field.F.rows() << " " << field.F.cols() << std::endl;


    } else {
        assert(0);
        // field.Load((std::string(DATA_PATH) + "/fertility.obj").c_str());
    }

    printf("Initialize...\n");
    t1 = GetCurrentTime64();
    field.Initialize(faces);
    t2 = GetCurrentTime64();
    printf("Use %lf seconds\n", (t2 - t1) * 1e-3);

    if (field.flag_preserve_boundary) {
        printf("Add boundary constrains...\n");
        Hierarchy& mRes = field.hierarchy;
        mRes.clearConstraints();
        for (uint32_t i = 0; i < 3 * mRes.mF.cols(); ++i) {
            if (mRes.mE2E[i] == -1) {
                uint32_t i0 = mRes.mF(i % 3, i / 3);
                uint32_t i1 = mRes.mF((i + 1) % 3, i / 3);
                Vector3d p0 = mRes.mV[0].col(i0), p1 = mRes.mV[0].col(i1);
                Vector3d edge = p1 - p0;
                if (edge.squaredNorm() > 0) {
                    edge.normalize();
                    mRes.mCO[0].col(i0) = p0;
                    mRes.mCO[0].col(i1) = p1;
                    mRes.mCQ[0].col(i0) = mRes.mCQ[0].col(i1) = edge;
                    mRes.mCQw[0][i0] = mRes.mCQw[0][i1] = mRes.mCOw[0][i0] = mRes.mCOw[0][i1] =
                        1.0;
                }
            }
        }
        mRes.propagateConstraints();
    }

    printf("Solve Orientation Field...\n");
    t1 = GetCurrentTime64();

    Optimizer::optimize_orientations(field.hierarchy);
    field.ComputeOrientationSingularities();
    t2 = GetCurrentTime64();
    printf("Use %lf seconds\n", (t2 - t1) * 1e-3);

    if (field.flag_adaptive_scale == 1) {
        printf("Estimate Slop...\n");
        t1 = GetCurrentTime64();
        field.EstimateSlope();
        t2 = GetCurrentTime64();
        printf("Use %lf seconds\n", (t2 - t1) * 1e-3);
    }
    printf("Solve for scale...\n");
    t1 = GetCurrentTime64();
    Optimizer::optimize_scale(field.hierarchy, field.rho, field.flag_adaptive_scale);
    field.flag_adaptive_scale = 1;
    t2 = GetCurrentTime64();
    printf("Use %lf seconds\n", (t2 - t1) * 1e-3);

    printf("Solve for position field...\n");
    t1 = GetCurrentTime64();
    Optimizer::optimize_positions(field.hierarchy, field.flag_adaptive_scale);

    printf("Rows and cols for faces %d, %d \n", field.hierarchy.mF.rows(), field.hierarchy.mF.cols() );
    printf("Rows and cols for vertices %d, %d \n", field.hierarchy.mV[0].rows(), field.hierarchy.mV[0].cols());
    printf("Rows and cols for orientation field %d, %d \n", field.hierarchy.mQ[0].rows(), field.hierarchy.mQ[0].cols());
    printf("Rows and cols for position field %d, %d \n", field.hierarchy.mO[0].rows(), field.hierarchy.mO[0].cols());

    field.ComputePositionSingularities();
    t2 = GetCurrentTime64();
    printf("Use %lf seconds\n", (t2 - t1) * 1e-3);
    t1 = GetCurrentTime64();
    printf("Solve index map...\n");

    std::cout << field.V.rows() << " " << field.V.cols() << std::endl;
    std::cout << field.F.rows() << " " << field.F.cols() << std::endl;

    field.ComputeIndexMap();

    std::cout << field.V.rows() << " " << field.V.cols() << std::endl;
    std::cout << field.F.rows() << " " << field.F.cols() << std::endl;

    t2 = GetCurrentTime64();
    printf("Indexmap Use %lf seconds\n", (t2 - t1) * 1e-3);
    
    printf("Rows and cols for faces %d, %d \n", field.hierarchy.mF.rows(), field.hierarchy.mF.cols());
    printf("Rows and cols for vertices %d, %d \n", field.hierarchy.mV[0].rows(), field.hierarchy.mV[0].cols());
    printf("Rows and cols for orientation field %d, %d \n", field.hierarchy.mQ[0].rows(), field.hierarchy.mQ[0].cols());
    printf("Rows and cols for position field %d, %d \n", field.hierarchy.mO[0].rows(), field.hierarchy.mO[0].cols());

    printf("Rows and cols for compact position field %d %d, %d \n", field.O_compact.size(), field.O_compact[0].rows(), field.O_compact[0].cols());
    printf("Rows and cols for compact faces %d %d, %d \n", field.F_compact.size(), field.F_compact[0].rows(), field.F_compact[0].cols());
    printf("Rows and cols for Vset %d %d \n", field.Vset.size(), field.Vset[0].size());


    std::cout << field.V.rows() << " " << field.V.cols() << std::endl;
    std::cout << field.F.rows() << " " << field.F.cols() << std::endl;


    /////////////////////////////////////// TESTING ////////////////////////////////////////////
    //int index = 13;
    //Eigen::Vector3d positionFieldSample(0, 0, 0);

    //for (int i = 0; i < field.Vset[index].size(); i++) {
    //    positionFieldSample += field.hierarchy.mO[0].col(field.Vset[index][i]);
    //    //positionFieldSample += hierarchy.mO[0].col(Vset[index][i]);

    //}
    //positionFieldSample /= field.Vset[index].size();

    //std::cout << "Vset verification - " << field.O_compact[index].transpose() << "       " << positionFieldSample.transpose() << std::endl;
    /////////////////////////////////////// TESTING ////////////////////////////////////////////

    // Export fine to coarse mappings as a CSV file
   /* std::ofstream file(directory + "/" + filename + "_quad_quadriflow_mappings_" + std::to_string(faces) + ".csv");
    for (const auto& row : field.Vset) {
        for (size_t i = 0; i < row.size(); ++i) {
            file << row[i];
            if (i != row.size() - 1) file << ",";
        }
        file << "\n";
    }
    file.close();
    std::cout << "Successfully written mappings" << std::endl;*/

    // Export as a json file
    // Create JSON and directly push the rows into it
    nlohmann::json json_output = field.Vset;

    // Write to file
    std::ofstream(directory + "/" + filename + "_quad_quadriflow_mappings_" + std::to_string(faces) + ".json")
        << json_output.dump(4); // Pretty print with indentation

    std::cout << "Successfully written mappings to JSON" << std::endl;

    std::cout << field.V.rows() << " " << field.V.cols() << std::endl;
    std::cout << field.F.rows() << " " << field.F.cols() << std::endl;


    printf("Writing the quad obj...\n");
    if (output_obj.size() < 1) {
        assert(0);
        // field.OutputMesh((std::string(DATA_PATH) + "/result.obj").c_str());
    } else {
        field.OutputMesh(output_obj.c_str());

        // Subdivided Input Mesh
        //std::ofstream os(directory + "/" + filename + "_tri_quadriflow" + extension);
        //for (int i = 0; i < field.hierarchy.mV[0].cols(); ++i) { // Vertices
        //    auto t = field.hierarchy.mV[0].col(i) * field.normalize_scale + field.normalize_offset;
        //    os << "v " << t[0] << " " << t[1] << " " << t[2] << "\n";
        //}
        //for (int i = 0; i < field.hierarchy.mF.cols(); ++i) { // Faces
        //    os << "f " << field.hierarchy.mF.col(i)[0] + 1 << " " << field.hierarchy.mF.col(i)[1] + 1
        //        << " " << field.hierarchy.mF.col(i)[2] + 1 /*<< " " << field.F.col(i)[3] + 1*/
        //        << "\n";
        //}
        //os.close();

        //

        std::ofstream os(directory + "/" + filename + "_tri_quadriflow_" + std::to_string(faces) + extension);

        // Write vertices
        for (int i = 0; i < field.hierarchy.mV[0].cols(); ++i) {
            auto t = field.hierarchy.mV[0].col(i) * field.normalize_scale + field.normalize_offset;
            os << "v " << t[0] << " " << t[1] << " " << t[2] << "\n";
        }

        // Write vertex normals
        for (int i = 0; i < field.hierarchy.mN[0].cols(); ++i) {
            auto n = field.hierarchy.mN[0].col(i);  // assume already normalized
            os << "vn " << n[0] << " " << n[1] << " " << n[2] << "\n";
        }

        // Write faces, referencing normals
        for (int i = 0; i < field.hierarchy.mF.cols(); ++i) {
            os << "f ";
            for (int j = 0; j < 3; ++j) {
                int vi = field.hierarchy.mF.col(i)[j] + 1;  // OBJ is 1-based
                os << vi << "//" << vi << " ";
            }
            os << "\n";
        }

        os.close();

    }
    printf("finish...\n");
    //	field.LoopFace(2);


    return 0;
}
