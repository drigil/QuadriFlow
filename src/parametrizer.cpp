#include "parametrizer.hpp"
#include "config.hpp"
#include "dedge.hpp"
#include "field-math.hpp"
#include "flow.hpp"
#include "localsat.hpp"
#include "optimizer.hpp"
#include "subdivide.hpp"

#include "dset.hpp"

#include <Eigen/Sparse>
#include <fstream>
#include <list>
#include <map>
#include <queue>
#include <set>

namespace qflow {

void Parametrizer::ComputeIndexMap(int with_scale) {
    // build edge info
    auto& V = hierarchy.mV[0]; // Vertices
    auto& F = hierarchy.mF; // Faces
    auto& Q = hierarchy.mQ[0]; // Orientation Field
    auto& N = hierarchy.mN[0]; // Normals
    auto& O = hierarchy.mO[0]; // Position Field
    auto& S = hierarchy.mS[0]; // ?
    // ComputeOrientationSingularities();

    BuildEdgeInfo();

    if (flag_preserve_sharp) {
        //        ComputeSharpO();
    }

    for (int i = 0; i < sharp_edges.size(); ++i) {
        if (sharp_edges[i]) {
            int e = face_edgeIds[i / 3][i % 3];
            if (edge_diff[e][0] * edge_diff[e][1] != 0) {
                Vector3d d = O.col(edge_values[e].y) - O.col(edge_values[e].x);
                Vector3d q = Q.col(edge_values[e].x);
                Vector3d n = N.col(edge_values[e].x);
                Vector3d qy = n.cross(q);
                if (abs(q.dot(d)) > qy.dot(d))
                    edge_diff[e][1] = 0;
                else
                    edge_diff[e][0] = 0;
            }
        }
    }
    std::map<int, std::pair<Vector3d, Vector3d>> sharp_constraints;
    std::set<int> sharpvert;
    for (int i = 0; i < sharp_edges.size(); ++i) {
        if (sharp_edges[i]) {
            sharpvert.insert(F(i % 3, i / 3));
            sharpvert.insert(F((i + 1) % 3, i / 3));
        }
    }

    allow_changes.resize(edge_diff.size() * 2, 1);
    for (int i = 0; i < sharp_edges.size(); ++i) {
        int e = face_edgeIds[i / 3][i % 3];
        if (sharpvert.count(edge_values[e].x) && sharpvert.count(edge_values[e].y)) {
            if (sharp_edges[i] != 0) {
                for (int k = 0; k < 2; ++k) {
                    if (edge_diff[e][k] == 0) {
                        allow_changes[e * 2 + k] = 0;
                    }
                }
            }
        }
    }
#ifdef LOG_OUTPUT
    printf("Build Integer Constraints...\n");
#endif
    BuildIntegerConstraints();

    ComputeMaxFlow();
    // potential bug
#ifdef LOG_OUTPUT
    printf("subdivide...\n");
#endif
    printf("Before subdivide_edgeDiff \n");
    printf("Rows and cols for faces %d, %d \n", hierarchy.mF.rows(), hierarchy.mF.cols());
    printf("Rows and cols for vertices %d, %d \n", hierarchy.mV[0].rows(), hierarchy.mV[0].cols());
    printf("Rows and cols for orientation field %d, %d \n", hierarchy.mQ[0].rows(), hierarchy.mQ[0].cols());
    printf("Rows and cols for position field %d, %d \n", hierarchy.mO[0].rows(), hierarchy.mO[0].cols());

    subdivide_edgeDiff(F, V, N, Q, O, &hierarchy.mS[0], V2E, hierarchy.mE2E, boundary, nonManifold,
                       edge_diff, edge_values, face_edgeOrients, face_edgeIds, sharp_edges,
                       singularities, 1);

    printf("After subdivide_edgeDiff \n");
    printf("Rows and cols for faces %d, %d \n", hierarchy.mF.rows(), hierarchy.mF.cols());
    printf("Rows and cols for vertices %d, %d \n", hierarchy.mV[0].rows(), hierarchy.mV[0].cols());
    printf("Rows and cols for orientation field %d, %d \n", hierarchy.mQ[0].rows(), hierarchy.mQ[0].cols());
    printf("Rows and cols for position field %d, %d \n", hierarchy.mO[0].rows(), hierarchy.mO[0].cols());


    allow_changes.clear();
    allow_changes.resize(edge_diff.size() * 2, 1);
    for (int i = 0; i < sharp_edges.size(); ++i) {
        if (sharp_edges[i] == 0) continue;
        int e = face_edgeIds[i / 3][i % 3];
        for (int k = 0; k < 2; ++k) {
            if (edge_diff[e][k] == 0) allow_changes[e * 2 + k] = 0;
        }
    }

#ifdef LOG_OUTPUT
    printf("Fix flip advance...\n");
    int t1 = GetCurrentTime64();
#endif
    FixFlipHierarchy();
    subdivide_edgeDiff(F, V, N, Q, O, &hierarchy.mS[0], V2E, hierarchy.mE2E, boundary, nonManifold,
                       edge_diff, edge_values, face_edgeOrients, face_edgeIds, sharp_edges,
                       singularities, 1);
    FixFlipSat();

#ifdef LOG_OUTPUT
    int t2 = GetCurrentTime64();
    printf("Flip use %lf\n", (t2 - t1) * 1e-3);
    printf("Post Linear Solver...\n");
#endif
    std::set<int> sharp_vertices;
    for (int i = 0; i < sharp_edges.size(); ++i) {
        if (sharp_edges[i] == 1) {
            sharp_vertices.insert(F(i % 3, i / 3));
            sharp_vertices.insert(F((i + 1) % 3, i / 3));
        }
    }

    printf("Before Positions Sharp \n");
    printf("Rows and cols for faces %d, %d \n", hierarchy.mF.rows(), hierarchy.mF.cols());
    printf("Rows and cols for vertices %d, %d \n", hierarchy.mV[0].rows(), hierarchy.mV[0].cols());
    printf("Rows and cols for orientation field %d, %d \n", hierarchy.mQ[0].rows(), hierarchy.mQ[0].cols());
    printf("Rows and cols for position field %d, %d \n", hierarchy.mO[0].rows(), hierarchy.mO[0].cols());


    Optimizer::optimize_positions_sharp(hierarchy, edge_values, edge_diff, sharp_edges,
                                        sharp_vertices, sharp_constraints, with_scale);

    Optimizer::optimize_positions_fixed(hierarchy, edge_values, edge_diff, sharp_vertices,
                                        sharp_constraints, flag_adaptive_scale);

    printf("Before Advanced Extract Quad \n");
    printf("Rows and cols for faces %d, %d \n", hierarchy.mF.rows(), hierarchy.mF.cols());
    printf("Rows and cols for vertices %d, %d \n", hierarchy.mV[0].rows(), hierarchy.mV[0].cols());
    printf("Rows and cols for orientation field %d, %d \n", hierarchy.mQ[0].rows(), hierarchy.mQ[0].cols());
    printf("Rows and cols for position field %d, %d \n", hierarchy.mO[0].rows(), hierarchy.mO[0].cols());

    AdvancedExtractQuad();
    printf("After Advanced Extract Quad \n");
    printf("Rows and cols for faces %d, %d \n", hierarchy.mF.rows(), hierarchy.mF.cols());
    printf("Rows and cols for vertices %d, %d \n", hierarchy.mV[0].rows(), hierarchy.mV[0].cols());
    printf("Rows and cols for orientation field %d, %d \n", hierarchy.mQ[0].rows(), hierarchy.mQ[0].cols());
    printf("Rows and cols for position field %d, %d \n", hierarchy.mO[0].rows(), hierarchy.mO[0].cols());

    printf("Rows and cols for compact position field %d %d, %d \n", O_compact.size(), O_compact[0].rows(), O_compact[0].cols());
    printf("Rows and cols for compact faces %d %d, %d \n", F_compact.size(), F_compact[0].rows(), F_compact[0].cols());
    printf("Rows and cols for Vset %d %d \n", Vset.size(), Vset[0].size());

    // Export fine to coarse mappings as a CSV file
    std::ofstream file("quad_quadriflow_mappings.csv");
    for (const auto& row : Vset) {
        for (size_t i = 0; i < row.size(); ++i) {
            file << row[i];
            if (i != row.size() - 1) file << ",";
        }
        file << "\n";
    }
    file.close();
    std::cout << "Successfully written mappings" << std::endl;
    // std::exit(0);

    FixValence();
    printf("After Fix Valence \n");
    printf("Rows and cols for compact position field %d %d, %d \n", O_compact.size(), O_compact[0].rows(), O_compact[0].cols());
    printf("Rows and cols for compact faces %d %d, %d \n", F_compact.size(), F_compact[0].rows(), F_compact[0].cols());
    printf("Rows and cols for Vset %d %d \n", Vset.size(), Vset[0].size());



    /////////////////////////////////////// TESTING ////////////////////////////////////////////
    //int index = 13;
    //Eigen::Vector3d positionFieldSample(0, 0, 0);

    //for (int i = 0; i < Vset[index].size(); i++) {
    //    positionFieldSample += O.col(Vset[index][i]);
    //    //positionFieldSample += hierarchy.mO[0].col(Vset[index][i]);

    //}
    //positionFieldSample /= Vset[index].size();

    //std::cout << "Vset verification - " << O_compact[index].transpose() << "       " << positionFieldSample.transpose() << std::endl;

    /////////////////////////////////////// TESTING ////////////////////////////////////////////


    std::vector<int> sharp_o(O_compact.size(), 0);
    std::map<int, std::pair<Vector3d, Vector3d>> compact_sharp_constraints;
    for (int i = 0; i < Vset.size(); ++i) {
        int sharpv = -1;
        for (auto& p : Vset[i]) {
            if (sharp_constraints.count(p)) {
                sharpv = p;
                sharp_o[i] = 1;
                if (compact_sharp_constraints.count(i) == 0 ||
                    compact_sharp_constraints[i].second != Vector3d::Zero()) {
                    compact_sharp_constraints[i] = sharp_constraints[sharpv];
                    O_compact[i] = O.col(sharpv);
                    compact_sharp_constraints[i].first = O_compact[i];
                }
            }
        }
    }

    std::map<std::pair<int, int>, int> o2e;
    for (int i = 0; i < F_compact.size(); ++i) {
        for (int j = 0; j < 4; ++j) {
            int v1 = F_compact[i][j];
            int v2 = F_compact[i][(j + 1) % 4];
            o2e[std::make_pair(v1, v2)] = i * 4 + j;
        }
    }
    std::vector<std::vector<int>> v2o(V.cols());
    for (int i = 0; i < Vset.size(); ++i) {
        for (auto v : Vset[i]) {
            v2o[v].push_back(i);
        }
    }
    std::vector<Vector3d> diffs(F_compact.size() * 4, Vector3d(0, 0, 0));
    std::vector<int> diff_count(F_compact.size() * 4, 0);
    for (int i = 0; i < F.cols(); ++i) {
        for (int j = 0; j < 3; ++j) {
            int v1 = F(j, i);
            int v2 = F((j + 1) % 3, i);
            if (v1 != edge_values[face_edgeIds[i][j]].x) continue;
            if (edge_diff[face_edgeIds[i][j]].array().abs().sum() != 1) continue;
            if (v2o[v1].size() > 1 || v2o[v2].size() > 1) continue;
            for (auto o1 : v2o[v1]) {
                for (auto o2 : v2o[v2]) {
                    auto key = std::make_pair(o1, o2);
                    if (o2e.count(key)) {
                        int dedge = o2e[key];
                        Vector3d q_1 = Q.col(v1);
                        Vector3d q_2 = Q.col(v2);
                        Vector3d n_1 = N.col(v1);
                        Vector3d n_2 = N.col(v2);
                        Vector3d q_1_y = n_1.cross(q_1);
                        Vector3d q_2_y = n_2.cross(q_2);
                        auto index = compat_orientation_extrinsic_index_4(q_1, n_1, q_2, n_2);
                        double s_x1 = S(0, v1), s_y1 = S(1, v1);
                        double s_x2 = S(0, v2), s_y2 = S(1, v2);
                        int rank_diff = (index.second + 4 - index.first) % 4;
                        if (rank_diff % 2 == 1) std::swap(s_x2, s_y2);
                        Vector3d qd_x = 0.5 * (rotate90_by(q_2, n_2, rank_diff) + q_1);
                        Vector3d qd_y = 0.5 * (rotate90_by(q_2_y, n_2, rank_diff) + q_1_y);
                        double scale_x = (with_scale ? 0.5 * (s_x1 + s_x2) : 1) * hierarchy.mScale;
                        double scale_y = (with_scale ? 0.5 * (s_y1 + s_y2) : 1) * hierarchy.mScale;
                        Vector2i diff = edge_diff[face_edgeIds[i][j]];
                        Vector3d C = diff[0] * scale_x * qd_x + diff[1] * scale_y * qd_y;

                        diff_count[dedge] += 1;
                        diffs[dedge] += C;
                        auto key = std::make_pair(o2, o1);
                        if (o2e.count(key)) {
                            int dedge = o2e[key];
                            diff_count[dedge] += 1;
                            diffs[dedge] -= C;
                        }
                    }
                }
            }
        }
    }

    for (int i = 0; i < F.cols(); ++i) {
        Vector2i d1 = rshift90(edge_diff[face_edgeIds[i][0]], face_edgeOrients[i][0]);
        Vector2i d2 = rshift90(edge_diff[face_edgeIds[i][1]], face_edgeOrients[i][1]);
        if (d1[0] * d2[1] - d1[1] * d2[0] < 0) {
            for (int j = 0; j < 3; ++j) {
                int v1 = F(j, i);
                int v2 = F((j + 1) % 3, i);
                for (auto o1 : v2o[v1]) {
                    for (auto o2 : v2o[v2]) {
                        auto key = std::make_pair(o1, o2);
                        if (o2e.count(key)) {
                            int dedge = o2e[key];
                            diff_count[dedge] = 0;
                            diffs[dedge] = Vector3d(0, 0, 0);
                        }
                    }
                }
            }
        }
    }

    for (int i = 0; i < diff_count.size(); ++i) {
        if (diff_count[i] != 0) {
            diffs[i] /= diff_count[i];
            diff_count[i] = 1;
        }
    }

    Optimizer::optimize_positions_dynamic(F, V, N, Q, Vset, O_compact, F_compact, V2E_compact,
                                          E2E_compact, sqrt(surface_area / F_compact.size()),
                                          diffs, diff_count, o2e, sharp_o,
                                          compact_sharp_constraints, flag_adaptive_scale);

    //    optimize_quad_positions(O_compact, N_compact, Q_compact, F_compact, V2E_compact,
    //    E2E_compact,
    //                            V, N, Q, O, F, V2E, hierarchy.mE2E, disajoint_tree,
    //                            hierarchy.mScale, false);
}

}  // namespace qflow
