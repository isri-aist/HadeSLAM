//G. Caron, Mar. 2021

#ifndef OPENVSLAM_OPTIMIZER_G2O_SE3_OMNIDIRECTIONAL_REPROJ_EDGE_H
#define OPENVSLAM_OPTIMIZER_G2O_SE3_OMNIDIRECTIONAL_REPROJ_EDGE_H

#include "openvslam/type.h"
#include "openvslam/optimize/g2o/landmark_vertex.h"
#include "openvslam/optimize/g2o/se3/shot_vertex.h"

#include <g2o/core/base_binary_edge.h>

namespace openvslam {
namespace optimize {
namespace g2o {
namespace se3 {

class omnidirectional_reproj_edge final : public ::g2o::BaseBinaryEdge<2, Vec2_t, landmark_vertex, shot_vertex> {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    omnidirectional_reproj_edge();

    bool read(std::istream& is) override;

    bool write(std::ostream& os) const override;

    void computeError() override {
        const auto v1 = static_cast<const shot_vertex*>(_vertices.at(1));
        const auto v2 = static_cast<const landmark_vertex*>(_vertices.at(0));
        const Vec2_t obs(_measurement);
        _error = obs - cam_project(v1->estimate().map(v2->estimate()));
    }

    void linearizeOplus() override;

    bool depth_is_positive() const {
        const auto v1 = static_cast<const shot_vertex*>(_vertices.at(1));
        const auto v2 = static_cast<const landmark_vertex*>(_vertices.at(0));
        return 0.0 < (v1->estimate().map(v2->estimate()))(2);
    }

    inline Vec2_t cam_project(const Vec3_t& pos_c) const {      
        // reproject onto the image
        // sphere or 3D to normalized image plane
        const double den = pos_c(2) + xi_ * sqrt(pos_c(0)*pos_c(0)+pos_c(1)*pos_c(1)+pos_c(2)*pos_c(2));

        return {fx_ * pos_c(0) / den + cx_, fy_ * pos_c(1) / den + cy_};
    }

    double fx_, fy_, cx_, cy_, xi_;
};

} // namespace se3
} // namespace g2o
} // namespace optimize
} // namespace openvslam

#endif // OPENVSLAM_OPTIMIZER_G2O_SE3_OMNIDIRECTIONAL_REPROJ_EDGE_H
