//G. Caron, Mar. 2021, rebase Feb. 2023
#include "stellavslam/camera/omnidirectional.h"

#include <iostream>

#include <spdlog/spdlog.h>
#include <nlohmann/json.hpp>

namespace stellavslam {
namespace camera {

omnidirectional::omnidirectional(const std::string& name, const color_order_t& color_order,
                         const unsigned int cols, const unsigned int rows, const double fps,
                         const double fx, const double fy, const double cx, const double cy,
                         const double xi)
    : base(name, setup_type_t::Monocular, model_type_t::Omnidirectional, color_order, cols, rows, fps, 0, 0),
      fx_(fx), fy_(fy), cx_(cx), cy_(cy), fx_inv_(1.0 / fx), fy_inv_(1.0 / fy),
      xi_(xi) {
    spdlog::debug("CONSTRUCT: camera::omnidirectional");

    img_bounds_ = compute_image_bounds();

    inv_cell_width_ = static_cast<double>(num_grid_cols_) / (img_bounds_.max_x_ - img_bounds_.min_x_);
    inv_cell_height_ = static_cast<double>(num_grid_rows_) / (img_bounds_.max_y_ - img_bounds_.min_y_);
}

omnidirectional::omnidirectional(const YAML::Node& yaml_node)
    : omnidirectional(yaml_node["Camera.name"].as<std::string>(),
                  load_color_order(yaml_node),
                  yaml_node["Camera.cols"].as<unsigned int>(),
                  yaml_node["Camera.rows"].as<unsigned int>(),
                  yaml_node["Camera.fps"].as<double>(),
                  yaml_node["Camera.fx"].as<double>(),
                  yaml_node["Camera.fy"].as<double>(),
                  yaml_node["Camera.cx"].as<double>(),
                  yaml_node["Camera.cy"].as<double>(),
                  yaml_node["Camera.xi"].as<double>()) {}

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

void omnidirectional::undistort_keypoints(const std::vector<cv::KeyPoint>& dist_keypts, std::vector<cv::KeyPoint>& undist_keypts) const {
    undist_keypts = dist_keypts;
}

void omnidirectional::convert_keypoints_to_bearings(const std::vector<cv::KeyPoint>& undist_keypts, eigen_alloc_vector<Vec3_t>& bearings) const {
    bearings.resize(undist_keypts.size());
    for (unsigned long idx = 0; idx < undist_keypts.size(); ++idx) {
        // Barreto's model
        Vec3_t bearing;

        // digital (pixel) to normalized
        const auto x_normalized = (undist_keypts.at(idx).pt.x - cx_) / fx_;
        const auto y_normalized = (undist_keypts.at(idx).pt.y - cy_) / fy_;
        auto fact = 1.0+(1.0-xi_*xi_)*(x_normalized*x_normalized + y_normalized*y_normalized);

        if(fact >= 0.0)
        {
          fact = (xi_ + std::sqrt(fact)) / (x_normalized*x_normalized + y_normalized*y_normalized + 1.0);
          // convert to spherical coordinates
          bearings.at(idx) = Vec3_t{fact * x_normalized, fact * y_normalized, fact - xi_};
        }
        else
        {
          //bearings.at(idx) = Vec3_t{0, 0, 1}; //rough
          bearings.at(idx) = Vec3_t{x_normalized, y_normalized, 1}; // should project to the closest on sphere
          //bearings.erase(bearings.begin()+idx); //could be could if discarding undist_keypts too
        }
    }
}

void omnidirectional::convert_bearings_to_keypoints(const eigen_alloc_vector<Vec3_t>& bearings, std::vector<cv::KeyPoint>& undist_keypts) const {
    undist_keypts.resize(bearings.size());
    for (unsigned long idx = 0; idx < bearings.size(); ++idx) {
        // sphere or 3D to normalized image plane
        /*const*/ auto den = bearings.at(idx)(2) + xi_ * sqrt(bearings.at(idx)(0)*bearings.at(idx)(0)+bearings.at(idx)(1)*bearings.at(idx)(1)+bearings.at(idx)(2)*bearings.at(idx)(2));
        if(fabs(den) < 1e-10)
        	den = 1.;
        	 
        const auto x_normalized = bearings.at(idx)(0) / den;
        const auto y_normalized = bearings.at(idx)(1) / den;

        // convert to pixel image coordinated
        undist_keypts.at(idx).pt.x = fx_ * x_normalized + cx_;
        undist_keypts.at(idx).pt.y = fy_ * y_normalized + cy_;
    }
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

} // namespace camera
} // namespace stellavslam
