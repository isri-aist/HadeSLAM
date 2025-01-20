#ifndef STELLA_VSLAM_OPTIMIZER_G2O_SE3_EQUIRECTANGULAR_REPROJ_EDGE_H
#define STELLA_VSLAM_OPTIMIZER_G2O_SE3_EQUIRECTANGULAR_REPROJ_EDGE_H

#include "stella_vslam/type.h"
#include "stella_vslam/optimize/internal/landmark_vertex.h"
#include "stella_vslam/optimize/internal/se3/shot_vertex.h"

#include <g2o/core/base_binary_edge.h>

namespace stella_vslam {
namespace optimize {
namespace internal {
namespace se3 {

class equirectangular_reproj_edge final : public g2o::BaseBinaryEdge<2, Vec2_t, landmark_vertex, shot_vertex> {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    equirectangular_reproj_edge();

    bool read(std::istream& is) override;

    bool write(std::ostream& os) const override;

    void computeError() override;

    void linearizeOplus() override;

    Vec2_t cam_project(const Vec3_t& pos_c) const;

    double cols_, rows_;
};

inline equirectangular_reproj_edge::equirectangular_reproj_edge()
    : g2o::BaseBinaryEdge<2, Vec2_t, landmark_vertex, shot_vertex>() {}

inline bool equirectangular_reproj_edge::read(std::istream& is) {
    for (unsigned int i = 0; i < 2; ++i) {
        is >> _measurement(i);
    }
    for (int i = 0; i < information().rows(); ++i) {
        for (int j = i; j < information().cols(); ++j) {
            is >> information()(i, j);
            if (i != j) {
                information()(j, i) = information()(i, j);
            }
        }
    }
    return true;
}

inline bool equirectangular_reproj_edge::write(std::ostream& os) const {
    for (unsigned int i = 0; i < 2; ++i) {
        os << measurement()(i) << " ";
    }
    for (int i = 0; i < information().rows(); ++i) {
        for (int j = i; j < information().cols(); ++j) {
            os << " " << information()(i, j);
        }
    }
    return os.good();
}

inline void equirectangular_reproj_edge::computeError() {
    const auto v1 = static_cast<const shot_vertex*>(_vertices.at(1));
    const auto v2 = static_cast<const landmark_vertex*>(_vertices.at(0));
    const Vec2_t obs(_measurement);
    _error = obs - cam_project(v1->estimate().map(v2->estimate()));
}

inline void equirectangular_reproj_edge::linearizeOplus() {
    auto vj = static_cast<shot_vertex*>(_vertices.at(1));
    const g2o::SE3Quat& cam_pose_cw = vj->shot_vertex::estimate();
    const Mat33_t rot_cw = cam_pose_cw.rotation().toRotationMatrix();

    auto vi = static_cast<landmark_vertex*>(_vertices.at(0));
    const Vec3_t& pos_w = vi->landmark_vertex::estimate();
    const Vec3_t pos_c = cam_pose_cw.map(pos_w);

    const auto pcx = pos_c(0);
    const auto pcy = pos_c(1);
    const auto pcz = pos_c(2);
    const auto L = pos_c.norm();

    // 回転に対する微分
    const Vec3_t d_pc_d_rx(0, -pcz, pcy);
    const Vec3_t d_pc_d_ry(pcz, 0, -pcx);
    const Vec3_t d_pc_d_rz(-pcy, pcx, 0);
    // 並進に対する微分
    const Vec3_t d_pc_d_tx(1, 0, 0);
    const Vec3_t d_pc_d_ty(0, 1, 0);
    const Vec3_t d_pc_d_tz(0, 0, 1);
    // 3次元点に対する微分
    const Vec3_t d_pc_d_pwx = rot_cw.block<3, 1>(0, 0);
    const Vec3_t d_pc_d_pwy = rot_cw.block<3, 1>(0, 1);
    const Vec3_t d_pc_d_pwz = rot_cw.block<3, 1>(0, 2);

    // 状態ベクトルを x = [rx, ry, rz, tx, ty, tz, pwx, pwy, pwz] として，
    // 導関数ベクトル d_pcx_d_x, d_pcy_d_x, d_pcz_d_x を作成
    VecR_t<9> d_pcx_d_x;
    d_pcx_d_x << d_pc_d_rx(0), d_pc_d_ry(0), d_pc_d_rz(0),
        d_pc_d_tx(0), d_pc_d_ty(0), d_pc_d_tz(0),
        d_pc_d_pwx(0), d_pc_d_pwy(0), d_pc_d_pwz(0);
    VecR_t<9> d_pcy_d_x;
    d_pcy_d_x << d_pc_d_rx(1), d_pc_d_ry(1), d_pc_d_rz(1),
        d_pc_d_tx(1), d_pc_d_ty(1), d_pc_d_tz(1),
        d_pc_d_pwx(1), d_pc_d_pwy(1), d_pc_d_pwz(1);
    VecR_t<9> d_pcz_d_x;
    d_pcz_d_x << d_pc_d_rx(2), d_pc_d_ry(2), d_pc_d_rz(2),
        d_pc_d_tx(2), d_pc_d_ty(2), d_pc_d_tz(2),
        d_pc_d_pwx(2), d_pc_d_pwy(2), d_pc_d_pwz(2);

    // 導関数ベクトル d_L_d_x を作成
    const VecR_t<9> d_L_d_x = (1.0 / L) * (pcx * d_pcx_d_x + pcy * d_pcy_d_x + pcz * d_pcz_d_x);

    // ヤコビ行列を作成
    MatRC_t<2, 9> jacobian = MatRC_t<2, 9>::Zero();
    jacobian.block<1, 9>(0, 0) = -(cols_ / (2 * M_PI)) * (1.0 / (pcx * pcx + pcz * pcz))
                                 * (pcz * d_pcx_d_x - pcx * d_pcz_d_x);
    jacobian.block<1, 9>(1, 0) = -(rows_ / M_PI) * (1.0 / (L * std::sqrt(pcx * pcx + pcz * pcz)))
                                 * (L * d_pcy_d_x - pcy * d_L_d_x);

    // g2oの変数にセット
    // 3次元点に対する微分
    _jacobianOplusXi = jacobian.block<2, 3>(0, 6);
    // 姿勢に対する微分
    _jacobianOplusXj = jacobian.block<2, 6>(0, 0);
}

inline Vec2_t equirectangular_reproj_edge::cam_project(const Vec3_t& pos_c) const {
    const double theta = std::atan2(pos_c(0), pos_c(2));
    const double phi = -std::asin(pos_c(1) / pos_c.norm());
    return {cols_ * (0.5 + theta / (2 * M_PI)), rows_ * (0.5 - phi / M_PI)};
}

//G. Caron, June 2024
class stereo_equirectangular_reproj_edge final : public g2o::BaseBinaryEdge<3, Vec3_t, landmark_vertex, shot_vertex> {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    stereo_equirectangular_reproj_edge();

    bool read(std::istream& is) override;

    bool write(std::ostream& os) const override;

    void computeError() override;

    void linearizeOplus() override;

    bool depth_is_positive() const;

    Vec3_t cam_project(const Vec3_t& pos_c) const;

    double cols_, rows_, true_baseline_;
};

inline stereo_equirectangular_reproj_edge::stereo_equirectangular_reproj_edge()
    : g2o::BaseBinaryEdge<3, Vec3_t, landmark_vertex, shot_vertex>() {}

inline bool stereo_equirectangular_reproj_edge::read(std::istream& is) {
    for (unsigned int i = 0; i < 3; ++i) {
        is >> _measurement(i);
    }
    for (unsigned int i = 0; i < 3; ++i) {
        for (unsigned int j = i; j < 3; ++j) {
            is >> information()(i, j);
            if (i != j) {
                information()(j, i) = information()(i, j);
            }
        }
    }
    return true;
}

inline bool stereo_equirectangular_reproj_edge::write(std::ostream& os) const {
    for (unsigned int i = 0; i < 3; ++i) {
        os << measurement()(i) << " ";
    }
    for (unsigned int i = 0; i < 3; ++i) {
        for (unsigned int j = i; j < 3; ++j) {
            os << " " << information()(i, j);
        }
    }
    return os.good();
}

inline void stereo_equirectangular_reproj_edge::computeError() {
    const auto v1 = static_cast<const shot_vertex*>(_vertices.at(1));
    const auto v2 = static_cast<const landmark_vertex*>(_vertices.at(0));
    const Vec3_t obs(_measurement);
    _error = obs - cam_project(v1->estimate().map(v2->estimate()));
}

inline void stereo_equirectangular_reproj_edge::linearizeOplus() {
auto vj = static_cast<shot_vertex*>(_vertices.at(1));
    const g2o::SE3Quat& cam_pose_cw = vj->shot_vertex::estimate();
    const Mat33_t rot_cw = cam_pose_cw.rotation().toRotationMatrix();

    auto vi = static_cast<landmark_vertex*>(_vertices.at(0));
    const Vec3_t& pos_w = vi->landmark_vertex::estimate();
    const Vec3_t pos_c = cam_pose_cw.map(pos_w);

    const auto pcx = pos_c(0);
    const auto pcy = pos_c(1);
    const auto pcz = pos_c(2);
    const auto L = pos_c.norm();

    // 回転に対する微分
    const Vec3_t d_pc_d_rx(0, -pcz, pcy);
    const Vec3_t d_pc_d_ry(pcz, 0, -pcx);
    const Vec3_t d_pc_d_rz(-pcy, pcx, 0);
    // 並進に対する微分
    const Vec3_t d_pc_d_tx(1, 0, 0);
    const Vec3_t d_pc_d_ty(0, 1, 0);
    const Vec3_t d_pc_d_tz(0, 0, 1);
    // 3次元点に対する微分
    const Vec3_t d_pc_d_pwx = rot_cw.block<3, 1>(0, 0);
    const Vec3_t d_pc_d_pwy = rot_cw.block<3, 1>(0, 1);
    const Vec3_t d_pc_d_pwz = rot_cw.block<3, 1>(0, 2);

    // 状態ベクトルを x = [rx, ry, rz, tx, ty, tz, pwx, pwy, pwz] として，
    // 導関数ベクトル d_pcx_d_x, d_pcy_d_x, d_pcz_d_x を作成
    VecR_t<9> d_pcx_d_x;
    d_pcx_d_x << d_pc_d_rx(0), d_pc_d_ry(0), d_pc_d_rz(0),
        d_pc_d_tx(0), d_pc_d_ty(0), d_pc_d_tz(0),
        d_pc_d_pwx(0), d_pc_d_pwy(0), d_pc_d_pwz(0);
    VecR_t<9> d_pcy_d_x;
    d_pcy_d_x << d_pc_d_rx(1), d_pc_d_ry(1), d_pc_d_rz(1),
        d_pc_d_tx(1), d_pc_d_ty(1), d_pc_d_tz(1),
        d_pc_d_pwx(1), d_pc_d_pwy(1), d_pc_d_pwz(1);
    VecR_t<9> d_pcz_d_x;
    d_pcz_d_x << d_pc_d_rx(2), d_pc_d_ry(2), d_pc_d_rz(2),
        d_pc_d_tx(2), d_pc_d_ty(2), d_pc_d_tz(2),
        d_pc_d_pwx(2), d_pc_d_pwy(2), d_pc_d_pwz(2);

    // 導関数ベクトル d_L_d_x を作成
    const VecR_t<9> d_L_d_x = (1.0 / L) * (pcx * d_pcx_d_x + pcy * d_pcy_d_x + pcz * d_pcz_d_x);

    // ヤコビ行列を作成
    MatRC_t<2, 9> jacobian = MatRC_t<2, 9>::Zero();
    jacobian.block<1, 9>(0, 0) = -(cols_ / (2 * M_PI)) * (1.0 / (pcx * pcx + pcz * pcz))
                                 * (pcz * d_pcx_d_x - pcx * d_pcz_d_x);
    jacobian.block<1, 9>(1, 0) = -(rows_ / M_PI) * (1.0 / (L * std::sqrt(pcx * pcx + pcz * pcz)))
                                 * (L * d_pcy_d_x - pcy * d_L_d_x);

    //additional part for u_right
    const auto pcx_right = pos_c(0)+true_baseline_;
    const Vec3_t d_pc_right_d_ry(pcz, 0, -pcx_right);
    const Vec3_t d_pc_right_d_rz(-pcy, pcx_right, 0);
    VecR_t<9> d_pcx_right_d_x;
    d_pcx_right_d_x << d_pc_d_rx(0), d_pc_right_d_ry(0), d_pc_right_d_rz(0),
        d_pc_d_tx(0), d_pc_d_ty(0), d_pc_d_tz(0),
        d_pc_d_pwx(0), d_pc_d_pwy(0), d_pc_d_pwz(0);
    VecR_t<9> d_pcz_right_d_x;
    d_pcz_right_d_x << d_pc_d_rx(2), d_pc_right_d_ry(2), d_pc_right_d_rz(2),
        d_pc_d_tx(2), d_pc_d_ty(2), d_pc_d_tz(2),
        d_pc_d_pwx(2), d_pc_d_pwy(2), d_pc_d_pwz(2);
    MatRC_t<1, 9> jacobian_right = MatRC_t<1, 9>::Zero();
    jacobian_right = -(cols_ / (2 * M_PI)) * (1.0 / (pcx_right * pcx_right + pcz * pcz))
                                 * (pcz * d_pcx_right_d_x - pcx * d_pcz_right_d_x);


    // g2oの変数にセット
    // 3次元点に対する微分
    _jacobianOplusXi.block<2,3>(0,0) = jacobian.block<2, 3>(0, 6);
    // to do
    _jacobianOplusXi.block<1,3>(2,0) = jacobian_right.block<1,3>(0, 6); // derivative of u_right with respect to 3D point
    // 姿勢に対する微分
    _jacobianOplusXj.block<2,6>(0,0) = jacobian.block<2, 6>(0, 0);
    //to do
    _jacobianOplusXj.block<1,6>(2,0) = jacobian_right.block<1,6>(0, 0); // derivative of u_right with respect to the pose
}

inline bool stereo_equirectangular_reproj_edge::depth_is_positive() const {
    return true;
}

inline Vec3_t stereo_equirectangular_reproj_edge::cam_project(const Vec3_t& pos_c) const {
    const double theta = std::atan2(pos_c(0), pos_c(2));
    const double phi = -std::asin(pos_c(1) / pos_c.norm());
    return {cols_ * (0.5 + theta / (2 * M_PI)), rows_ * (0.5 - phi / M_PI), cols_ * (0.5 + std::atan2(pos_c(0) - true_baseline_, pos_c(2)) / (2.0 * M_PI))};
}

} // namespace se3
} // namespace internal
} // namespace optimize
} // namespace stella_vslam

#endif // STELLA_VSLAM_OPTIMIZER_G2O_SE3_EQUIRECTANGULAR_REPROJ_EDGE_H
