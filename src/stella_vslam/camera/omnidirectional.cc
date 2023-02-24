//G. Caron, Mar. 2021, rebase Feb. 2023
#include "stella_vslam/camera/omnidirectional.h"

#include <iostream>

#include <spdlog/spdlog.h>
#include <nlohmann/json.hpp>

namespace stella_vslam {
namespace camera {

omnidirectional::omnidirectional(const std::string& name, const color_order_t& color_order,
                         const unsigned int cols, const unsigned int rows, const double fps,
                         const double fx, const double fy, const double cx, const double cy,
                         const double xi)
    : base(name, setup_type_t::Monocular, model_type_t::Omnidirectional, color_order, cols, rows, fps, 0.0, 0.0, 0.0),
      fx_(fx), fy_(fy), cx_(cx), cy_(cy), fx_inv_(1.0 / fx), fy_inv_(1.0 / fy), xi_(xi) {
    spdlog::debug("CONSTRUCT: camera::omnidirectional");

    img_bounds_ = compute_image_bounds();

    inv_cell_width_ = static_cast<double>(num_grid_cols_) / (img_bounds_.max_x_ - img_bounds_.min_x_);
    inv_cell_height_ = static_cast<double>(num_grid_rows_) / (img_bounds_.max_y_ - img_bounds_.min_y_);
}

omnidirectional::omnidirectional(const YAML::Node& yaml_node)
    : omnidirectional(yaml_node["name"].as<std::string>(),
                  load_color_order(yaml_node),
                  yaml_node["cols"].as<unsigned int>(),
                  yaml_node["rows"].as<unsigned int>(),
                  yaml_node["fps"].as<double>(),
                  yaml_node["fx"].as<double>(),
                  yaml_node["fy"].as<double>(),
                  yaml_node["cx"].as<double>(),
                  yaml_node["cy"].as<double>(),
                  yaml_node["xi"].as<double>()) {}

omnidirectional::~omnidirectional() {
    spdlog::debug("DESTRUCT: camera::omnidirectional");
}

void omnidirectional::show_parameters() const {
    show_common_parameters();
    std::cout << "  - fx: " << fx_ << std::endl;
    std::cout << "  - fy: " << fy_ << std::endl;
    std::cout << "  - cx: " << cx_ << std::endl;
    std::cout << "  - cy: " << cy_ << std::endl;
    std::cout << "  - xi: " << xi_ << std::endl;
    std::cout << "  - min x: " << img_bounds_.min_x_ << std::endl;
    std::cout << "  - max x: " << img_bounds_.max_x_ << std::endl;
    std::cout << "  - min y: " << img_bounds_.min_y_ << std::endl;
    std::cout << "  - max y: " << img_bounds_.max_y_ << std::endl;
}

image_bounds omnidirectional::compute_image_bounds() const {
    spdlog::debug("compute image bounds");

        return image_bounds{0.0, cols_, 0.0, rows_};

}

cv::Point2f omnidirectional::undistort_point(const cv::Point2f& dist_pt) const {
    return dist_pt;
}

Vec3_t omnidirectional::convert_point_to_bearing(const cv::Point2f& undist_pt) const {
		// digital (pixel) to normalized
    const auto x_normalized = (undist_pt.x - cx_) / fx_;
    const auto y_normalized = (undist_pt.y - cy_) / fy_;
    auto fact = 1.0+(1.0-xi_*xi_)*(x_normalized*x_normalized + y_normalized*y_normalized);

    if(fact >= 0.0)
    {
      fact = (xi_ + std::sqrt(fact)) / (x_normalized*x_normalized + y_normalized*y_normalized + 1.0);
      // convert to spherical coordinates
      return Vec3_t{fact * x_normalized, fact * y_normalized, fact - xi_};
    }
    else
    {
      //bearings.at(idx) = Vec3_t{0, 0, 1}; //rough
      return Vec3_t{x_normalized, y_normalized, 1}; // should project to the closest on sphere
      //bearings.erase(bearings.begin()+idx); //could be could if discarding undist_keypts too
    }

}

cv::Point2f omnidirectional::convert_bearing_to_point(const Vec3_t& bearing) const {
		// sphere or 3D to normalized image plane
    /*const*/ auto den = bearing[2] + xi_ * std::sqrt(bearing[0]*bearing[0]+bearing[1]*bearing[1]+bearing[2]*bearing[2]);
    if(fabs(den) < 1e-10)
    	den = 1.;
    	 
    const auto x_normalized = bearing[0] / den;
    const auto y_normalized = bearing[1] / den;

    // convert to pixel image coordinated
		return cv::Point2f(fx_ * x_normalized + cx_, fy_ * y_normalized + cy_);

}

bool omnidirectional::reproject_to_image(const Mat33_t& rot_cw, const Vec3_t& trans_cw, const Vec3_t& pos_w, Vec2_t& reproj, float& x_right) const {
    // convert to camera-coordinates
    const Vec3_t pos_c = rot_cw * pos_w + trans_cw;

    /*
    // check if the point is visible
    if (pos_c(2) <= 0.0) {
        return false;
    }
    */

    // reproject onto the image
    // sphere or 3D to normalized image plane
    const auto den = pos_c(2) + xi_ * sqrt(pos_c(0)*pos_c(0)+pos_c(1)*pos_c(1)+pos_c(2)*pos_c(2));
    reproj(0) = fx_ * pos_c(0) / den + cx_;
    reproj(1) = fy_ * pos_c(1) / den + cy_;

    //x_right = reproj(0) - focal_x_baseline_ * z_inv;
    x_right = 0.0;

    // // check if the point is visible
    //return (img_bounds_.min_x_ < reproj(0) && reproj(0) < img_bounds_.max_x_
    //        && img_bounds_.min_y_ < reproj(1) && reproj(1) < img_bounds_.max_y_);
    return true;
}

bool omnidirectional::reproject_to_bearing(const Mat33_t& rot_cw, const Vec3_t& trans_cw, const Vec3_t& pos_w, Vec3_t& reproj) const {
    // convert to camera-coordinates
    reproj = (rot_cw * pos_w + trans_cw).normalized();

    return true;
}

nlohmann::json omnidirectional::to_json() const {
    return {{"model_type", get_model_type_string()},
            {"setup_type", get_setup_type_string()},
            {"color_order", get_color_order_string()},
            {"cols", cols_},
            {"rows", rows_},
            {"fps", fps_},
            {"focal_x_baseline", focal_x_baseline_},
            {"num_grid_cols", num_grid_cols_},
            {"num_grid_rows", num_grid_rows_},
            {"fx", fx_},
            {"fy", fy_},
            {"cx", cx_},
            {"cy", cy_},
            {"xi", xi_}};
}

std::ostream& operator<<(std::ostream& os, const omnidirectional& params) {
    os << "- name: " << params.name_ << std::endl;
    os << "- setup: " << params.get_setup_type_string() << std::endl;
    os << "- fps: " << params.fps_ << std::endl;
    os << "- cols: " << params.cols_ << std::endl;
    os << "- rows: " << params.rows_ << std::endl;
    os << "- color: " << params.get_color_order_string() << std::endl;
    os << "- model: " << params.get_model_type_string() << std::endl;
    os << "- focal x baseline: " << params.focal_x_baseline_ << std::endl;
    os << "- num grid cols: " << params.num_grid_cols_ << std::endl;
    os << "- num grid rows: " << params.num_grid_rows_ << std::endl;
    os << "  - fx: " << params.fx_ << std::endl;
    os << "  - fy: " << params.fy_ << std::endl;
    os << "  - cx: " << params.cx_ << std::endl;
    os << "  - cy: " << params.cy_ << std::endl;
    os << "  - xi: " << params.xi_ << std::endl;
    return os;
}

//! Override for optimization

void omnidirectional::undistort_points(const std::vector<cv::Point2f>& dist_pts, std::vector<cv::Point2f>& undist_pts) const {
    undist_pts = dist_pts;
}

void omnidirectional::undistort_keypoints(const std::vector<cv::KeyPoint>& dist_keypts, std::vector<cv::KeyPoint>& undist_keypts) const {
    undist_keypts = dist_keypts;
}

} // namespace camera
} // namespace stella_vslam
