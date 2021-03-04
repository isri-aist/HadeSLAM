//G. Caron, Mar. 2021
#include "openvslam/optimize/g2o/se3/omnidirectional_pose_opt_edge.h"

namespace openvslam {
namespace optimize {
namespace g2o {
namespace se3 {

omnidirectional_pose_opt_edge::omnidirectional_pose_opt_edge()
    : ::g2o::BaseUnaryEdge<2, Vec2_t, shot_vertex>() {}

bool omnidirectional_pose_opt_edge::read(std::istream& is) {
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

bool omnidirectional_pose_opt_edge::write(std::ostream& os) const {
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

void omnidirectional_pose_opt_edge::linearizeOplus() {
    auto vi = static_cast<shot_vertex*>(_vertices.at(0));  //3D point
    const ::g2o::SE3Quat& cam_pose_cw = vi->shot_vertex::estimate(); //camera pose
    const Vec3_t pos_c = cam_pose_cw.map(pos_w_); //3D point in the world frame

    const auto X = pos_c(0);
    const auto Y = pos_c(1);
    const auto Z = pos_c(2);

    const double rho = sqrt(X*X+Y*Y+Z*Z);
    const double den = Z + xi_ * rho;

    const auto x = X / den;
    const auto x2 = x*x;
    const auto y = Y / den;
    const auto y2 = y*y;

    const auto alpha = std::sqrt(1+(1-xi_*xi_)*(x2+y2)); //be careful if xi_ > 1

    // 姿勢に対する微分 / Differentiation for pose: du / dp
    _jacobianOplusXi(0, 0) = fx_ * x*y; // du / dr_X
    _jacobianOplusXi(0, 1) = -fx_ * ((1+x2)*alpha-xi_*y2)/(alpha + xi_) ; // du / dr_Y
    _jacobianOplusXi(0, 2) = fx_ * y; //du / dr_Z
    _jacobianOplusXi(0, 3) = -fx_*(1+x2*(1-xi_*(alpha+xi_))+y2)/(rho*(alpha+xi_)); // du / dt_X
    _jacobianOplusXi(0, 4) = fx_ * xi_*x*y / rho; //du / dt_Y
    _jacobianOplusXi(0, 5) = fx_ * alpha * x / rho; //du / dt_Z

    _jacobianOplusXi(1, 0) = fy_ * ((1+y2)*alpha-xi_*x2)/(alpha + xi_); // dv / dr_X
    _jacobianOplusXi(1, 1) = -fy_ * x*y; // dv / dr_Y
    _jacobianOplusXi(1, 2) = -fy_ * x; //dv / dr_Z
    _jacobianOplusXi(1, 3) = fy_ * xi_*x*y / rho; // dv / dt_X
    _jacobianOplusXi(1, 4) = -fy_*(1+y2*(1-xi_*(alpha+xi_))+x2)/(rho*(alpha+xi_)); //dv / dt_Y
    _jacobianOplusXi(1, 5) = fy_ * alpha * y / rho; //dv / dt_Z

}

} // namespace se3
} // namespace g2o
} // namespace optimize
} // namespace openvslam
