#include "readcsv.h"
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <Eigen/Dense>


int readPointCSV(const std::string& filename, std::vector<p_pair>& points) {
    std::ifstream file(filename);
    std::string line;

    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return 0;
    }

    //std::getline(file, line);

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string token;
        std::vector<std::string> columns;

        while (std::getline(ss, token, ',')) {
            columns.push_back(token);
        }

        const size_t sz = columns.size();
        if (sz < 6) {
            std::cerr << "Warning: Not enough columns in line: " << line << std::endl;
            continue;
        }

        try {
            double x1 = std::stod(columns[sz - 6]);
            double y1 = std::stod(columns[sz - 5]);
            double z1 = std::stod(columns[sz - 4]);

            double x2 = std::stod(columns[sz - 3]);
            double y2 = std::stod(columns[sz - 2]);
            double z2 = std::stod(columns[sz - 1]);

            points.emplace_back(
                Eigen::Vector3d(x1, y1, z1),
                Eigen::Vector3d(x2, y2, z2)
            );
        }
        catch (const std::exception& e) {
            std::cerr << "Warning: Invalid number format in line: " << line << std::endl;
            continue;
        }
    }

    std::cout << "Read " << points.size() << " points from " << filename << std::endl;
    return 1;
}