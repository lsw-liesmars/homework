/**
 * PnP算法测试程序
 * 包含多个测试用例，验证不同场景下的求解精度
 */

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <random>

using namespace std;

// 复用主程序中的类定义
#include "pnp_solver.cpp"

// 添加噪声
void addNoise(vector<vector<double>>& points, double noise_level) {
    random_device rd;
    mt19937 gen(rd());
    normal_distribution<> d(0, noise_level);
    
    for (auto& pt : points) {
        for (auto& coord : pt) {
            coord += d(gen);
        }
    }
}

// 测试用例1：立方体顶点（无噪声）
void test1_cube_perfect() {
    cout << "\n" << string(80, '=') << endl;
    cout << "测试用例1：立方体顶点（理想情况，无噪声）" << endl;
    cout << string(80, '=') << endl;
    
    vector<vector<double>> world_points = {
        {0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0},
        {0, 0, 1}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1}
    };
    
    Matrix K(3, 3, {800, 0, 320, 0, 800, 240, 0, 0, 1});
    
    // 设置真实外参
    double angle = 25.0 * M_PI / 180.0;
    Matrix R_true(3, 3, {
        cos(angle), -sin(angle), 0,
        sin(angle), cos(angle), 0,
        0, 0, 1
    });
    Matrix T_true(3, 1, {{-1.5}, {0.8}, {4.0}});
    
    // 生成图像点
    vector<vector<double>> image_points;
    Matrix RT(3, 4);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) RT.data[i][j] = R_true.data[i][j];
        RT.data[i][3] = T_true.data[i][0];
    }
    
    for (const auto& wp : world_points) {
        Matrix X(4, 1, {wp[0], wp[1], wp[2], 1.0});
        Matrix x = K * RT * X;
        image_points.push_back({x.data[0][0] / x.data[2][0], x.data[1][0] / x.data[2][0]});
    }
    
    // 求解
    Matrix R_solved, T_solved;
    if (PnPSolver::solve(image_points, world_points, K, R_solved, T_solved)) {
        verifyResult(image_points, world_points, K, R_solved, T_solved);
    }
}

// 测试用例2：平面上的点
void test2_planar_points() {
    cout << "\n" << string(80, '=') << endl;
    cout << "测试用例2：平面点配置" << endl;
    cout << string(80, '=') << endl;
    
    vector<vector<double>> world_points = {
        {0, 0, 0}, {2, 0, 0}, {4, 0, 0},
        {0, 2, 0}, {2, 2, 0}, {4, 2, 0},
        {0, 4, 0}, {2, 4, 0}, {4, 4, 0}
    };
    
    Matrix K(3, 3, {1000, 0, 512, 0, 1000, 384, 0, 0, 1});
    
    double angle_x = 15.0 * M_PI / 180.0;
    double angle_z = 45.0 * M_PI / 180.0;
    Matrix R_true(3, 3, {
        cos(angle_z), -sin(angle_z), 0,
        sin(angle_z) * cos(angle_x), cos(angle_z) * cos(angle_x), -sin(angle_x),
        sin(angle_z) * sin(angle_x), cos(angle_z) * sin(angle_x), cos(angle_x)
    });
    Matrix T_true(3, 1, {{-2.0}, {1.5}, {8.0}});
    
    // 生成图像点
    vector<vector<double>> image_points;
    Matrix RT(3, 4);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) RT.data[i][j] = R_true.data[i][j];
        RT.data[i][3] = T_true.data[i][0];
    }
    
    for (const auto& wp : world_points) {
        Matrix X(4, 1, {wp[0], wp[1], wp[2], 1.0});
        Matrix x = K * RT * X;
        image_points.push_back({x.data[0][0] / x.data[2][0], x.data[1][0] / x.data[2][0]});
    }
    
    // 求解
    Matrix R_solved, T_solved;
    if (PnPSolver::solve(image_points, world_points, K, R_solved, T_solved)) {
        verifyResult(image_points, world_points, K, R_solved, T_solved);
    }
}

// 测试用例3：带噪声的数据
void test3_with_noise() {
    cout << "\n" << string(80, '=') << endl;
    cout << "测试用例3：带噪声的观测数据（噪声标准差=1.0像素）" << endl;
    cout << string(80, '=') << endl;
    
    vector<vector<double>> world_points = {
        {0, 0, 0}, {3, 0, 0}, {0, 3, 0}, {3, 3, 0},
        {0, 0, 3}, {3, 0, 3}, {0, 3, 3}, {3, 3, 3},
        {1.5, 1.5, 0}, {1.5, 1.5, 3}
    };
    
    Matrix K(3, 3, {900, 0, 400, 0, 900, 300, 0, 0, 1});
    
    double angle = 20.0 * M_PI / 180.0;
    Matrix R_true(3, 3, {
        cos(angle), 0, sin(angle),
        0, 1, 0,
        -sin(angle), 0, cos(angle)
    });
    Matrix T_true(3, 1, {{-1.0}, {0.5}, {6.0}});
    
    // 生成图像点
    vector<vector<double>> image_points;
    Matrix RT(3, 4);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) RT.data[i][j] = R_true.data[i][j];
        RT.data[i][3] = T_true.data[i][0];
    }
    
    for (const auto& wp : world_points) {
        Matrix X(4, 1, {wp[0], wp[1], wp[2], 1.0});
        Matrix x = K * RT * X;
        image_points.push_back({x.data[0][0] / x.data[2][0], x.data[1][0] / x.data[2][0]});
    }
    
    // 添加噪声
    addNoise(image_points, 1.0);
    
    // 求解
    Matrix R_solved, T_solved;
    if (PnPSolver::solve(image_points, world_points, K, R_solved, T_solved)) {
        verifyResult(image_points, world_points, K, R_solved, T_solved);
    }
}

// 测试用例4：最少点数（6个点）
void test4_minimum_points() {
    cout << "\n" << string(80, '=') << endl;
    cout << "测试用例4：最少点数配置（6个点）" << endl;
    cout << string(80, '=') << endl;
    
    vector<vector<double>> world_points = {
        {0, 0, 0}, {5, 0, 0}, {0, 5, 0},
        {5, 5, 0}, {2.5, 2.5, 5}, {0, 0, 5}
    };
    
    Matrix K(3, 3, {850, 0, 320, 0, 850, 240, 0, 0, 1});
    
    double angle = 35.0 * M_PI / 180.0;
    Matrix R_true(3, 3, {
        1, 0, 0,
        0, cos(angle), -sin(angle),
        0, sin(angle), cos(angle)
    });
    Matrix T_true(3, 1, {{0.0}, {-2.0}, {10.0}});
    
    // 生成图像点
    vector<vector<double>> image_points;
    Matrix RT(3, 4);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) RT.data[i][j] = R_true.data[i][j];
        RT.data[i][3] = T_true.data[i][0];
    }
    
    for (const auto& wp : world_points) {
        Matrix X(4, 1, {wp[0], wp[1], wp[2], 1.0});
        Matrix x = K * RT * X;
        image_points.push_back({x.data[0][0] / x.data[2][0], x.data[1][0] / x.data[2][0]});
    }
    
    // 求解
    Matrix R_solved, T_solved;
    if (PnPSolver::solve(image_points, world_points, K, R_solved, T_solved)) {
        verifyResult(image_points, world_points, K, R_solved, T_solved);
    }
}

// 测试用例5：大角度旋转
void test5_large_rotation() {
    cout << "\n" << string(80, '=') << endl;
    cout << "测试用例5：大角度旋转（60度）" << endl;
    cout << string(80, '=') << endl;
    
    vector<vector<double>> world_points = {
        {-2, -2, 0}, {2, -2, 0}, {-2, 2, 0}, {2, 2, 0},
        {-2, -2, 2}, {2, -2, 2}, {-2, 2, 2}, {2, 2, 2}
    };
    
    Matrix K(3, 3, {750, 0, 320, 0, 750, 240, 0, 0, 1});
    
    double angle = 60.0 * M_PI / 180.0;
    double c = cos(angle), s = sin(angle);
    Matrix R_true(3, 3, {
        c, -s, 0,
        s, c, 0,
        0, 0, 1
    });
    Matrix T_true(3, 1, {{1.0}, {-1.0}, {7.0}});
    
    // 生成图像点
    vector<vector<double>> image_points;
    Matrix RT(3, 4);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) RT.data[i][j] = R_true.data[i][j];
        RT.data[i][3] = T_true.data[i][0];
    }
    
    for (const auto& wp : world_points) {
        Matrix X(4, 1, {wp[0], wp[1], wp[2], 1.0});
        Matrix x = K * RT * X;
        image_points.push_back({x.data[0][0] / x.data[2][0], x.data[1][0] / x.data[2][0]});
    }
    
    // 求解
    Matrix R_solved, T_solved;
    if (PnPSolver::solve(image_points, world_points, K, R_solved, T_solved)) {
        verifyResult(image_points, world_points, K, R_solved, T_solved);
    }
}

int main() {
    cout << "\n";
    cout << "##############################################################################" << endl;
    cout << "#                                                                            #" << endl;
    cout << "#          无人机摄影测量 PnP 算法测试套件                                  #" << endl;
    cout << "#          测试 x = K[R|T]X 的求解算法                                      #" << endl;
    cout << "#                                                                            #" << endl;
    cout << "##############################################################################" << endl;
    
    // 运行所有测试
    test1_cube_perfect();
    test2_planar_points();
    test3_with_noise();
    test4_minimum_points();
    test5_large_rotation();
    
    cout << "\n";
    cout << "##############################################################################" << endl;
    cout << "#                          所有测试完成！                                    #" << endl;
    cout << "##############################################################################" << endl;
    cout << endl;
    
    return 0;
}
