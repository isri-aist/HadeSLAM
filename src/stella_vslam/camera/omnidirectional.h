//G. Caron, Mar. 2021, rebase Feb. 2023

#ifndef STELLA_VSLAM_CAMERA_OMNIDIRECTIONAL_H
#define STELLA_VSLAM_CAMERA_OMNIDIRECTIONAL_H

#include "stellavslam/camera/base.h"

namespace stellavslam {
namespace camera {

class omnidirectional final : public base {
public:
    omnidirectional(const std::string& name, const color_order_t& color_order,
                    const unsigned int cols, const unsigned int rows, const double fps,
                    const double fx, const double fy, const double cx, const double cy,
                    const double xi);

    omnidirectional(const YAML::Node& yaml_node);

    ~omnidirectional() override;

    void show_parameters() const override final;

    image_bounds compute_image_bounds() const override final;

    cv::KeyPoint undistort_keypoint(const cv::KeyPoint& distorted) const override final {
        return distorted;
    }

    void undistort_keypoints(const std::vector<cv::KeyPoint>& dist_keypts, std::vector<cv::KeyPoint>& undist_keypts) const override final;

    Vec3_t convert_keypoint_to_bearing(const cv::KeyPoint& undist_keypt) const override final {
        // Barreto's model
        Vec3_t bearing;

        // digital (pixel) to normalized
        const auto x_normalized = (undist_keypt.pt.x - cx_) / fx_;
        const auto y_normalized = (undist_keypt.pt.y - cy_) / fy_;
        auto fact = 1.0+(1.0-xi_*xi_)*(x_normalized*x_normalized + y_normalized*y_normalized);

        if(fact >= 0.0)
        {
          fact = (xi_ + std::sqrt(fact)) / (x_normalized*x_normalized + y_normalized*y_normalized + 1.0);

          // convert to spherical coordinates
          bearing[0] = fact * x_normalized;
          bearing[1] = fact * y_normalized;
          bearing[2] = fact - xi_;
        }
        else
        {
        	// Rough
          //bearing[0] = bearing[1] = 0.0;
          // should project to the closest on sphere
          bearing[0] = x_normalized;
          bearing[1] = y_normalized;
          bearing[2] = 1.0;
        }

        return bearing;
    }

    void convert_keypoints_to_bearings(const std::vector<cv::KeyPoint>& undist_keypts, eigen_alloc_vector<Vec3_t>& bearings) const override final;

    cv::KeyPoint convert_bearing_to_keypoint(const Vec3_t& bearing) const override final {
        cv::KeyPoint undist_keypt;

        // sphere or 3D to normalized image plane
        const auto den = bearing(2) + xi_ * sqrt(bearing(0)*bearing(0)+bearing(1)*bearing(1)+bearing(2)*bearing(2));
        const auto x_normalized = bearing(0) / den;
        const auto y_normalized = bearing(1) / den;

        // convert to pixel image coordinated
        undist_keypt.pt.x = fx_ * x_normalized + cx_;
        undist_keypt.pt.y = fy_ * y_normalized + cy_;

        return undist_keypt;
    }

    void convert_bearings_to_keypoints(const eigen_alloc_vector<Vec3_t>& bearings, std::vector<cv::KeyPoint>& undist_keypts) const override final;

    bool reproject_to_image(const Mat33_t& rot_cw, const Vec3_t& trans_cw, const Vec3_t& pos_w, Vec2_t& reproj, float& x_right) const override final;

    bool reproject_to_bearing(const Mat33_t& rot_cw, const Vec3_t& trans_cw, const Vec3_t& pos_w, Vec3_t& reproj) const override final;

    nlohmann::json to_json() const override final;

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

} // namespace camera
} // namespace stellavslam

#endif // STELLA_VSLAM_CAMERA_OMNIDIRECTIONAL_H
