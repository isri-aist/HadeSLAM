//G. Caron, Mar. 2021, rebase Feb. 2023

#ifndef STELLA_VSLAM_OPTIMIZER_G2O_SE3_OMNIDIRECTIONAL_REPROJ_EDGE_H
#define STELLA_VSLAM_OPTIMIZER_G2O_SE3_OMNIDIRECTIONAL_REPROJ_EDGE_H

#include "stella_vslam/type.h"
#include "stella_vslam/optimize/internal/landmark_vertex.h"
#include "stella_vslam/optimize/internal/se3/shot_vertex.h"

#include <g2o/core/base_binary_edge.h>

namespace stella_vslam {
namespace optimize {
namespace internal {
namespace se3 {

class omnidirectional_reproj_edge final : public g2o::BaseBinaryEdge<2, Vec2_t, landmark_vertex, shot_vertex> {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    omnidirectional_reproj_edge();

    bool read(std::istream& is) override;

    bool write(std::ostream& os) const override;

    void computeError() override;

    void linearizeOplus() override;

    Vec2_t cam_project(const Vec3_t& pos_c) const;

    double fx_, fy_, cx_, cy_, xi_;
};

inline omnidirectional_reproj_edge::omnidirectional_reproj_edge()
    : g2o::BaseBinaryEdge<2, Vec2_t, landmark_vertex, shot_vertex>() {}

inline bool omnidirectional_reproj_edge::read(std::istream& is) {
    for (unsigned int i = 0; i < 2; ++i) {
        is >> _measurement(i);
    }
    for (unsigned int i = 0; i < 2; ++i) {
        for (unsigned int j = i; j < 2; ++j) {
            is >> information()(i, j);
            if (i != j) {
                information()(j, i) = information()(i, j);
            }
        }
    }
    return true;
}

inline bool omnidirectional_reproj_edge::write(std::ostream& os) const {
    for (unsigned int i = 0; i < 2; ++i) {
        os << measurement()(i) << " ";
    }
    for (unsigned int i = 0; i < 2; ++i) {
        for (unsigned int j = i; j < 2; ++j) {
            os << " " << information()(i, j);
        }
    }
    return os.good();
}

inline void omnidirectional_reproj_edge::computeError() {
    const auto v1 = static_cast<const shot_vertex*>(_vertices.at(1));
    const auto v2 = static_cast<const landmark_vertex*>(_vertices.at(0));
    const Vec2_t obs(_measurement);
    _error = obs - cam_project(v1->estimate().map(v2->estimate()));
}

inline void omnidirectional_reproj_edge::linearizeOplus() {
    auto vj = static_cast<shot_vertex*>(_vertices.at(1)); //camera pose
    const ::g2o::SE3Quat& cam_pose_cw = vj->shot_vertex::estimate();

    auto vi = static_cast<landmark_vertex*>(_vertices.at(0)); //3D point
    const Vec3_t& pos_w = vi->landmark_vertex::estimate(); //3D point in the world frame
    const Vec3_t pos_c = cam_pose_cw.map(pos_w); //3D point in the camera frame

    const auto X = pos_c(0);
    const auto Y = pos_c(1);
    const auto Z = pos_c(2);

    const double rho = sqrt(X*X+Y*Y+Z*Z);
    const double den = Z + xi_ * rho;
    const double den2 = den*den;

    const auto x = X / den;
    const auto x2 = x*x;
    const auto y = Y / den;
    const auto y2 = y*y;

		auto fact = 1+(1-xi_*xi_)*(x2+y2);
		if(fact < 0)
			fact = 0;
    const auto alpha = std::sqrt(fact); //be careful if xi_ > 1

    const Mat33_t rot_cw = cam_pose_cw.rotation().toRotationMatrix();

    // g2oの変数にセット / Set to g2o variable
    // 3次元点に対する微分 /  Differentiation for 3-dimensional points: du / dwX
    _jacobianOplusXi(0, 0) = -fx_ * rot_cw(0, 0) / den + fx_ * X * rot_cw(2, 0) / den2; //-du / dwX
    _jacobianOplusXi(0, 1) = -fx_ * rot_cw(0, 1) / den + fx_ * X * rot_cw(2, 1) / den2; //-du / dwY
    _jacobianOplusXi(0, 2) = -fx_ * rot_cw(0, 2) / den + fx_ * X * rot_cw(2, 2) / den2; //-du / dwZ

    _jacobianOplusXi(1, 0) = -fy_ * rot_cw(1, 0) / den + fy_ * Y * rot_cw(2, 0) / den2; //-dv / dwX
    _jacobianOplusXi(1, 1) = -fy_ * rot_cw(1, 1) / den + fy_ * Y * rot_cw(2, 1) / den2; //-dv / dwY
    _jacobianOplusXi(1, 2) = -fy_ * rot_cw(1, 2) / den + fy_ * Y * rot_cw(2, 2) / den2; //-dv / dwZ

    // 姿勢に対する微分 / Differentiation for pose: du / dp
    _jacobianOplusXj(0, 0) = fx_ * x*y; // du / dr_X
    _jacobianOplusXj(0, 1) = -fx_ * ((1+x2)*alpha-xi_*y2)/(alpha + xi_) ; // du / dr_Y
    _jacobianOplusXj(0, 2) = fx_ * y; //du / dr_Z
    _jacobianOplusXj(0, 3) = -fx_*(1+x2*(1-xi_*(alpha+xi_))+y2)/(rho*(alpha+xi_)); // du / dt_X
    _jacobianOplusXj(0, 4) = fx_ * xi_*x*y / rho; //du / dt_Y
    _jacobianOplusXj(0, 5) = fx_ * alpha * x / rho; //du / dt_Z

    _jacobianOplusXj(1, 0) = fy_ * ((1+y2)*alpha-xi_*x2)/(alpha + xi_); // dv / dr_X
    _jacobianOplusXj(1, 1) = -fy_ * x*y; // dv / dr_Y
    _jacobianOplusXj(1, 2) = -fy_ * x; //dv / dr_Z
    _jacobianOplusXj(1, 3) = fy_ * xi_*x*y / rho; // dv / dt_X
    _jacobianOplusXj(1, 4) = -fy_*(1+y2*(1-xi_*(alpha+xi_))+x2)/(rho*(alpha+xi_)); //dv / dt_Y
    _jacobianOplusXj(1, 5) = fy_ * alpha * y / rho; //dv / dt_Z
}

inline Vec2_t omnidirectional_reproj_edge::cam_project(const Vec3_t& pos_c) const {      
    // reproject onto the image
    // sphere or 3D to normalized image plane
    const double den = pos_c(2) + xi_ * sqrt(pos_c(0)*pos_c(0)+pos_c(1)*pos_c(1)+pos_c(2)*pos_c(2));

    return {fx_ * pos_c(0) / den + cx_, fy_ * pos_c(1) / den + cy_};
}

} // namespace se3
} // namespace internal
} // namespace optimize
} // namespace stella_vslam

#endif // STELLA_VSLAM_OPTIMIZER_G2O_SE3_OMNIDIRECTIONAL_REPROJ_EDGE_H
