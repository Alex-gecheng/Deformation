#pragma once
    
#include <vector>
#include <string>
#include <Eigen/Dense>

typedef std::pair<int, Eigen::Vector3d> p_pair;

int readPointCSV(const std::string& filename, std::vector<p_pair>& points);

