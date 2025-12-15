/**
 * 无人机摄影测量作业：空间点求相机参数 (改进版)
 * 求解方程：x = K[R|T]X
 * 
 * 改进：
 * 1. 使用更稳定的DLT算法实现
 * 2. 改进SVD实现
 * 3. 更好的正交化处理
 */

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>

using namespace std;

// 矩阵类
class Matrix {
public:
    int rows, cols;
    vector<vector<double>> data;
    
    Matrix() : rows(0), cols(0) {}
    
    Matrix(int r, int c) : rows(r), cols(c) {
        data.resize(r, vector<double>(c, 0.0));
    }
    
    Matrix(int r, int c, const vector<double>& values) : rows(r), cols(c) {
        data.resize(r, vector<double>(c));
        for (int i = 0; i < r; i++) {
            for (int j = 0; j < c; j++) {
                data[i][j] = values[i * c + j];
            }
        }
    }
    
    double& operator()(int i, int j) { return data[i][j]; }
    const double& operator()(int i, int j) const { return data[i][j]; }
    
    Matrix operator*(const Matrix& other) const {
        Matrix result(rows, other.cols);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < other.cols; j++) {
                for (int k = 0; k < cols; k++) {
                    result(i, j) += data[i][k] * other(k, j);
                }
            }
        }
        return result;
    }
    
    Matrix transpose() const {
        Matrix result(cols, rows);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                result(j, i) = data[i][j];
            }
        }
        return result;
    }
    
    Matrix inverse() const {
        if (rows != cols) throw runtime_error("只能对方阵求逆");
        int n = rows;
        Matrix aug(n, 2 * n);
        
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) aug(i, j) = data[i][j];
            aug(i, i + n) = 1.0;
        }
        
        for (int i = 0; i < n; i++) {
            int maxRow = i;
            for (int k = i + 1; k < n; k++) {
                if (abs(aug(k, i)) > abs(aug(maxRow, i))) maxRow = k;
            }
            swap(aug.data[i], aug.data[maxRow]);
            
            double pivot = aug(i, i);
            if (abs(pivot) < 1e-10) throw runtime_error("矩阵不可逆");
            
            for (int j = 0; j < 2 * n; j++) aug(i, j) /= pivot;
            
            for (int k = 0; k < n; k++) {
                if (k != i) {
                    double factor = aug(k, i);
                    for (int j = 0; j < 2 * n; j++) {
                        aug(k, j) -= factor * aug(i, j);
                    }
                }
            }
        }
        
        Matrix result(n, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                result(i, j) = aug(i, j + n);
            }
        }
        return result;
    }
    
    void print(const string& name = "") const {
        if (!name.empty()) cout << name << ":" << endl;
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                cout << setw(12) << setprecision(6) << data[i][j] << " ";
            }
            cout << endl;
        }
        cout << endl;
    }
    
    double norm() const {
        double sum = 0;
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                sum += data[i][j] * data[i][j];
            }
        }
        return sqrt(sum);
    }
};

// 简化的特征值分解（用于对称矩阵）
void eigenSymmetric(const Matrix& A, Matrix& V, vector<double>& eigenvalues) {
    int n = A.rows;
    Matrix B = A;
    V = Matrix(n, n);
    for (int i = 0; i < n; i++) V(i, i) = 1.0;
    
    const int maxIter = 100;
    const double tol = 1e-10;
    
    for (int iter = 0; iter < maxIter; iter++) {
        double maxVal = 0;
        int p = 0, q = 1;
        for (int i = 0; i < n - 1; i++) {
            for (int j = i + 1; j < n; j++) {
                if (abs(B(i, j)) > maxVal) {
                    maxVal = abs(B(i, j));
                    p = i;
                    q = j;
                }
            }
        }
        
        if (maxVal < tol) break;
        
        double theta;
        if (abs(B(p, p) - B(q, q)) < 1e-10) {
            theta = M_PI / 4.0;
        } else {
            theta = 0.5 * atan(2.0 * B(p, q) / (B(p, p) - B(q, q)));
        }
        
        double c = cos(theta), s = sin(theta);
        
        // 应用Givens旋转
        for (int i = 0; i < n; i++) {
            double Bip = B(i, p), Biq = B(i, q);
            B(i, p) = c * Bip - s * Biq;
            B(i, q) = s * Bip + c * Biq;
        }
        for (int j = 0; j < n; j++) {
            double Bpj = B(p, j), Bqj = B(q, j);
            B(p, j) = c * Bpj - s * Bqj;
            B(q, j) = s * Bpj + c * Bqj;
        }
        for (int i = 0; i < n; i++) {
            double Vip = V(i, p), Viq = V(i, q);
            V(i, p) = c * Vip - s * Viq;
            V(i, q) = s * Vip + c * Viq;
        }
    }
    
    eigenvalues.resize(n);
    for (int i = 0; i < n; i++) eigenvalues[i] = B(i, i);
    
    // 排序（从大到小）
    for (int i = 0; i < n - 1; i++) {
        for (int j = i + 1; j < n; j++) {
            if (eigenvalues[j] > eigenvalues[i]) {
                swap(eigenvalues[i], eigenvalues[j]);
                for (int k = 0; k < n; k++) {
                    swap(V(k, i), V(k, j));
                }
            }
        }
    }
}

// SVD分解（针对矩阵求最小奇异向量）
Matrix svdMinVector(const Matrix& A) {
    Matrix AtA = A.transpose() * A;
    Matrix V;
    vector<double> eigenvalues;
    eigenSymmetric(AtA, V, eigenvalues);
    
    // 返回最小特征值对应的特征向量（最后一列）
    int n = V.rows;
    Matrix minVec(n, 1);
    for (int i = 0; i < n; i++) {
        minVec(i, 0) = V(i, n - 1);
    }
    return minVec;
}

// 强制正交化（Gram-Schmidt）
Matrix orthogonalize(const Matrix& R_approx) {
    Matrix R(3, 3);
    
    // 第一行：归一化
    double norm1 = 0;
    for (int j = 0; j < 3; j++) norm1 += R_approx(0, j) * R_approx(0, j);
    norm1 = sqrt(norm1);
    for (int j = 0; j < 3; j++) R(0, j) = R_approx(0, j) / norm1;
    
    // 第二行：正交化并归一化
    double dot = 0;
    for (int j = 0; j < 3; j++) dot += R_approx(1, j) * R(0, j);
    for (int j = 0; j < 3; j++) R(1, j) = R_approx(1, j) - dot * R(0, j);
    double norm2 = 0;
    for (int j = 0; j < 3; j++) norm2 += R(1, j) * R(1, j);
    norm2 = sqrt(norm2);
    for (int j = 0; j < 3; j++) R(1, j) /= norm2;
    
    // 第三行：叉乘
    R(2, 0) = R(0, 1) * R(1, 2) - R(0, 2) * R(1, 1);
    R(2, 1) = R(0, 2) * R(1, 0) - R(0, 0) * R(1, 2);
    R(2, 2) = R(0, 0) * R(1, 1) - R(0, 1) * R(1, 0);
    
    return R;
}

// PnP求解
bool solvePnP(const vector<vector<double>>& img_pts,
              const vector<vector<double>>& world_pts,
              const Matrix& K,
              Matrix& R,
              Matrix& T) {
    
    int n = img_pts.size();
    if (n < 6) {
        cerr << "至少需要6个点!" << endl;
        return false;
    }
    
    // 归一化图像坐标
    Matrix K_inv = K.inverse();
    vector<vector<double>> norm_pts(n, vector<double>(3));
    
    for (int i = 0; i < n; i++) {
        Matrix img(3, 1);
        img(0, 0) = img_pts[i][0];
        img(1, 0) = img_pts[i][1];
        img(2, 0) = 1.0;
        Matrix norm = K_inv * img;
        for (int j = 0; j < 3; j++) norm_pts[i][j] = norm(j, 0);
    }
    
    // 构建DLT方程 A * m = 0
    Matrix A(2 * n, 12);
    for (int i = 0; i < n; i++) {
        double x = norm_pts[i][0];
        double y = norm_pts[i][1];
        double X = world_pts[i][0];
        double Y = world_pts[i][1];
        double Z = world_pts[i][2];
        
        // 第一个方程
        A(2*i, 0) = X; A(2*i, 1) = Y; A(2*i, 2) = Z; A(2*i, 3) = 1;
        A(2*i, 4) = 0; A(2*i, 5) = 0; A(2*i, 6) = 0; A(2*i, 7) = 0;
        A(2*i, 8) = -x*X; A(2*i, 9) = -x*Y; A(2*i, 10) = -x*Z; A(2*i, 11) = -x;
        
        // 第二个方程
        A(2*i+1, 0) = 0; A(2*i+1, 1) = 0; A(2*i+1, 2) = 0; A(2*i+1, 3) = 0;
        A(2*i+1, 4) = X; A(2*i+1, 5) = Y; A(2*i+1, 6) = Z; A(2*i+1, 7) = 1;
        A(2*i+1, 8) = -y*X; A(2*i+1, 9) = -y*Y; A(2*i+1, 10) = -y*Z; A(2*i+1, 11) = -y;
    }
    
    // SVD求解
    Matrix sol = svdMinVector(A);
    
    // 重建外参矩阵
    Matrix M(3, 4);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 4; j++) {
            M(i, j) = sol(i * 4 + j, 0);
        }
    }
    
    // 提取R和T
    Matrix R_approx(3, 3);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            R_approx(i, j) = M(i, j);
        }
    }
    
    // 计算尺度因子
    double scale = R_approx.norm() / sqrt(3.0);
    
    // 正交化
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            R_approx(i, j) /= scale;
        }
    }
    R = orthogonalize(R_approx);
    
    // 提取平移
    T = Matrix(3, 1);
    for (int i = 0; i < 3; i++) {
        T(i, 0) = M(i, 3) / scale;
    }
    
    return true;
}

// 验证结果
void verifyResult(const vector<vector<double>>& img_pts,
                  const vector<vector<double>>& world_pts,
                  const Matrix& K,
                  const Matrix& R,
                  const Matrix& T) {
    cout << "========== 验证结果 ==========" << endl;
    
    double total_error = 0;
    int n = img_pts.size();
    
    for (int i = 0; i < n; i++) {
        Matrix X(4, 1);
        X(0, 0) = world_pts[i][0];
        X(1, 0) = world_pts[i][1];
        X(2, 0) = world_pts[i][2];
        X(3, 0) = 1.0;
        
        Matrix RT(3, 4);
        for (int r = 0; r < 3; r++) {
            for (int c = 0; c < 3; c++) {
                RT(r, c) = R(r, c);
            }
            RT(r, 3) = T(r, 0);
        }
        
        Matrix x_proj = K * RT * X;
        double u = x_proj(0, 0) / x_proj(2, 0);
        double v = x_proj(1, 0) / x_proj(2, 0);
        
        double error = sqrt((u - img_pts[i][0]) * (u - img_pts[i][0]) + 
                           (v - img_pts[i][1]) * (v - img_pts[i][1]));
        total_error += error;
        
        cout << "点" << i << ": 原始(" << img_pts[i][0] << ", " << img_pts[i][1] << ") "
             << "重投影(" << u << ", " << v << ") 误差=" << error << endl;
    }
    
    cout << "平均重投影误差: " << total_error / n << " 像素" << endl << endl;
}

int main() {
    cout << "========== 无人机摄影测量：PnP相机参数求解 (改进版) ==========" << endl << endl;
    
    // 测试用例
    vector<vector<double>> world_points = {
        {0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0},
        {0, 0, 1}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1}
    };
    
    Matrix K(3, 3, {800, 0, 320, 0, 800, 240, 0, 0, 1});
    cout << "相机内参矩阵 K:" << endl;
    K.print();
    
    // 真实外参
    double angle = 30.0 * M_PI / 180.0;
    Matrix R_true(3, 3, {
        cos(angle), 0, sin(angle),
        0, 1, 0,
        -sin(angle), 0, cos(angle)
    });
    Matrix T_true(3, 1, {{-2.0}, {1.0}, {5.0}});
    
    cout << "真实旋转矩阵 R:" << endl;
    R_true.print();
    cout << "真实平移向量 T:" << endl;
    T_true.print();
    
    // 生成图像点
    vector<vector<double>> image_points;
    Matrix RT_true(3, 4);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) RT_true(i, j) = R_true(i, j);
        RT_true(i, 3) = T_true(i, 0);
    }
    
    cout << "生成的图像点:" << endl;
    for (const auto& wp : world_points) {
        Matrix X(4, 1, {wp[0], wp[1], wp[2], 1.0});
        Matrix x = K * RT_true * X;
        double u = x(0, 0) / x(2, 0);
        double v = x(1, 0) / x(2, 0);
        image_points.push_back({u, v});
        cout << "  (" << u << ", " << v << ")" << endl;
    }
    cout << endl;
    
    // 求解
    Matrix R_solved, T_solved;
    cout << "========== 开始求解PnP ==========" << endl;
    if (solvePnP(image_points, world_points, K, R_solved, T_solved)) {
        cout << "求解成功！" << endl << endl;
        
        cout << "求解得到的旋转矩阵 R:" << endl;
        R_solved.print();
        cout << "求解得到的平移向量 T:" << endl;
        T_solved.print();
        
        verifyResult(image_points, world_points, K, R_solved, T_solved);
        
        cout << "========== 与真实值对比 ==========" << endl;
        cout << "旋转矩阵差异:" << endl;
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                cout << setw(12) << R_solved(i, j) - R_true(i, j) << " ";
            }
            cout << endl;
        }
        cout << endl;
        
        cout << "平移向量差异:" << endl;
        for (int i = 0; i < 3; i++) {
            cout << setw(12) << T_solved(i, 0) - T_true(i, 0) << endl;
        }
    }
    
    return 0;
}
