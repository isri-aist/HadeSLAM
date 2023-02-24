//G. Caron, Mar. 2021, rebase Feb. 2023

#ifndef STELLA_VSLAM_CAMERA_OMNIDIRECTIONAL_H
#define STELLA_VSLAM_CAMERA_OMNIDIRECTIONAL_H

#include "stella_vslam/camera/base.h"

namespace stella_vslam {
namespace camera {

class omnidirectional final : public base {
public:
 //   EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    omnidirectional(const std::string& name, const color_order_t& color_order,
                    const unsigned int cols, const unsigned int rows, const double fps,
                    const double fx, const double fy, const double cx, const double cy,
                    const double xi);

    omnidirectional(const YAML::Node& yaml_node);

    ~omnidirectional() override;

    void show_parameters() const override final;

    image_bounds compute_image_bounds() const override final;

		cv::Point2f undistort_point(const cv::Point2f& dist_pt) const override final;

    Vec3_t convert_point_to_bearing(const cv::Point2f& undist_pt) const override final;

    cv::Point2f convert_bearing_to_point(const Vec3_t& bearing) const override final;

    bool reproject_to_image(const Mat33_t& rot_cw, const Vec3_t& trans_cw, const Vec3_t& pos_w, Vec2_t& reproj, float& x_right) const override final;

    bool reproject_to_bearing(const Mat33_t& rot_cw, const Vec3_t& trans_cw, const Vec3_t& pos_w, Vec3_t& reproj) const override final;

    nlohmann::json to_json() const override final;

//! Override for optimization
    void undistort_points(const std::vector<cv::Point2f>& dist_pts, std::vector<cv::Point2f>& undist_pts) const override final;
    void undistort_keypoints(const std::vector<cv::KeyPoint>& dist_keypts, std::vector<cv::KeyPoint>& undist_keypts) const override final;

    //-------------------------
    // Parameters specific to this model

    //! pinhole-like params
    const double fx_;
    const double fy_;
    const double cx_;
    const double cy_;
    const double fx_inv_;
    const double fy_inv_;

    //! Barreto's model params
    const double xi_;

};

std::ostream& operator<<(std::ostream& os, const omnidirectional& params);

} // namespace camera
} // namespace stella_vslam

#endif // STELLA_VSLAM_CAMERA_OMNIDIRECTIONAL_H
