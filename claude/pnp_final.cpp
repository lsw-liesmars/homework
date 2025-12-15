/**
 * 无人机摄影测量作业：PnP算法最终完整版
 * 求解: x = K[R|T]X
 * 基于标准DLT + QR分解方法
 */

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>

using namespace std;

const double EPS = 1e-10;

// 3D向量
struct Vec3 {
    double x, y, z;
    Vec3(double x = 0, double y = 0, double z = 0) : x(x), y(y), z(z) {}
    Vec3 operator-(const Vec3& v) const { return Vec3(x - v.x, y - v.y, z - v.z); }
    Vec3 operator+(const Vec3& v) const { return Vec3(x + v.x, y + v.y, z + v.z); }
    Vec3 operator*(double s) const { return Vec3(x * s, y * s, z * s); }
    double dot(const Vec3& v) const { return x * v.x + y * v.y + z * v.z; }
    Vec3 cross(const Vec3& v) const {
        return Vec3(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);
    }
    double norm() const { return sqrt(x * x + y * y + z * z); }
    Vec3 normalized() const { double n = norm(); return Vec3(x / n, y / n, z / n); }
};

// 矩阵类
class Matrix {
public:
    int rows, cols;
    vector<double> data;
    
    Matrix() : rows(0), cols(0) {}
    Matrix(int r, int c, double val = 0) : rows(r), cols(c), data(r * c, val) {}
    Matrix(int r, int c, const vector<double>& d) : rows(r), cols(c), data(d) {}
    
    double& operator()(int i, int j) { return data[i * cols + j]; }
    const double& operator()(int i, int j) const { return data[i * cols + j]; }
    
    Matrix operator*(const Matrix& B) const {
        Matrix C(rows, B.cols);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < B.cols; j++) {
                for (int k = 0; k < cols; k++) {
                    C(i, j) += (*this)(i, k) * B(k, j);
                }
            }
        }
        return C;
    }
    
    Matrix transpose() const {
        Matrix T(cols, rows);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                T(j, i) = (*this)(i, j);
            }
        }
        return T;
    }
    
    Matrix inverse() const {
        if (rows != cols) throw runtime_error("只能对方阵求逆");
        int n = rows;
        vector<double> aug(n * 2 * n);
        
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                aug[i * 2 * n + j] = (*this)(i, j);
            }
            aug[i * 2 * n + i + n] = 1.0;
        }
        
        for (int i = 0; i < n; i++) {
            int pivot = i;
            for (int k = i + 1; k < n; k++) {
                if (abs(aug[k * 2 * n + i]) > abs(aug[pivot * 2 * n + i])) {
                    pivot = k;
                }
            }
            if (pivot != i) {
                for (int j = 0; j < 2 * n; j++) {
                    swap(aug[i * 2 * n + j], aug[pivot * 2 * n + j]);
                }
            }
            
            double d = aug[i * 2 * n + i];
            if (abs(d) < EPS) throw runtime_error("矩阵奇异");
            
            for (int j = 0; j < 2 * n; j++) {
                aug[i * 2 * n + j] /= d;
            }
            
            for (int k = 0; k < n; k++) {
                if (k != i) {
                    double f = aug[k * 2 * n + i];
                    for (int j = 0; j < 2 * n; j++) {
                        aug[k * 2 * n + j] -= f * aug[i * 2 * n + j];
                    }
                }
            }
        }
        
        Matrix inv(n, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                inv(i, j) = aug[i * 2 * n + j + n];
            }
        }
        return inv;
    }
    
    void print(const string& name = "") const {
        if (!name.empty()) cout << name << ":\n";
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                cout << setw(12) << setprecision(6) << (*this)(i, j) << " ";
            }
            cout << "\n";
        }
        cout << "\n";
    }
};

// 求对称矩阵的特征值和特征向量
void eigenDecomposition(const Matrix& A, vector<double>& eigenvalues, Matrix& eigenvectors) {
    int n = A.rows;
    Matrix B = A;
    eigenvectors = Matrix(n, n);
    for (int i = 0; i < n; i++) eigenvectors(i, i) = 1.0;
    
    for (int iter = 0; iter < 100; iter++) {
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
        
        if (maxVal < EPS) break;
        
        double theta;
        if (abs(B(p, p) - B(q, q)) < EPS) {
            theta = M_PI / 4.0;
        } else {
            theta = 0.5 * atan(2.0 * B(p, q) / (B(p, p) - B(q, q)));
        }
        
        double c = cos(theta), s = sin(theta);
        
        Matrix R = Matrix(n, n);
        for (int i = 0; i < n; i++) R(i, i) = 1.0;
        R(p, p) = c; R(q, q) = c;
        R(p, q) = -s; R(q, p) = s;
        
        B = R.transpose() * B * R;
        eigenvectors = eigenvectors * R;
    }
    
    eigenvalues.resize(n);
    for (int i = 0; i < n; i++) eigenvalues[i] = B(i, i);
    
    // 排序
    for (int i = 0; i < n - 1; i++) {
        for (int j = i + 1; j < n; j++) {
            if (eigenvalues[j] > eigenvalues[i]) {
                swap(eigenvalues[i], eigenvalues[j]);
                for (int k = 0; k < n; k++) {
                    swap(eigenvectors(k, i), eigenvectors(k, j));
                }
            }
        }
    }
}

// PnP求解器
bool solvePnP(const vector<pair<double, double>>& imagePoints,
              const vector<Vec3>& worldPoints,
              const Matrix& K,
              Matrix& R,
              Vec3& t) {
    
    int n = imagePoints.size();
    if (n < 6) {
        cerr << "需要至少6个点!" << endl;
        return false;
    }
    
    // 步骤1：将图像点转到归一化平面
    Matrix K_inv = K.inverse();
    vector<pair<double, double>> normPoints(n);
    
    for (int i = 0; i < n; i++) {
        double u = imagePoints[i].first;
        double v = imagePoints[i].second;
        double x = K_inv(0, 0) * u + K_inv(0, 1) * v + K_inv(0, 2);
        double y = K_inv(1, 0) * u + K_inv(1, 1) * v + K_inv(1, 2);
        double z = K_inv(2, 0) * u + K_inv(2, 1) * v + K_inv(2, 2);
        normPoints[i] = {x / z, y / z};
    }
    
    // 步骤2：构建DLT方程组
    Matrix A(2 * n, 12);
    for (int i = 0; i < n; i++) {
        double x = normPoints[i].first;
        double y = normPoints[i].second;
        double X = worldPoints[i].x;
        double Y = worldPoints[i].y;
        double Z = worldPoints[i].z;
        
        A(2*i, 0) = X;  A(2*i, 1) = Y;  A(2*i, 2) = Z;  A(2*i, 3) = 1;
        A(2*i, 4) = 0;  A(2*i, 5) = 0;  A(2*i, 6) = 0;  A(2*i, 7) = 0;
        A(2*i, 8) = -x*X; A(2*i, 9) = -x*Y; A(2*i, 10) = -x*Z; A(2*i, 11) = -x;
        
        A(2*i+1, 0) = 0;  A(2*i+1, 1) = 0;  A(2*i+1, 2) = 0;  A(2*i+1, 3) = 0;
        A(2*i+1, 4) = X;  A(2*i+1, 5) = Y;  A(2*i+1, 6) = Z;  A(2*i+1, 7) = 1;
        A(2*i+1, 8) = -y*X; A(2*i+1, 9) = -y*Y; A(2*i+1, 10) = -y*Z; A(2*i+1, 11) = -y;
    }
    
    // 步骤3：SVD求解（求A^T*A的最小特征向量）
    Matrix AtA = A.transpose() * A;
    vector<double> eigenvalues;
    Matrix eigenvectors;
    eigenDecomposition(AtA, eigenvalues, eigenvectors);
    
    // 最小特征值对应的特征向量
    vector<double> P(12);
    for (int i = 0; i < 12; i++) {
        P[i] = eigenvectors(i, 11);  // 最后一列
    }
    
    // 步骤4：重建投影矩阵并分解
    // P = [r1^T; r2^T; r3^T]，其中ri是行向量
    Vec3 r1(P[0], P[1], P[2]);
    Vec3 r2(P[4], P[5], P[6]);
    Vec3 r3(P[8], P[9], P[10]);
    Vec3 t_vec(P[3], P[7], P[11]);
    
    // 归一化尺度（使R的第三行norm接近1）
    double scale = r3.norm();
    if (scale < EPS) {
        cerr << "尺度因子过小!" << endl;
        return false;
    }
    
    r1 = r1 * (1.0 / scale);
    r2 = r2 * (1.0 / scale);
    r3 = r3 * (1.0 / scale);
    t_vec = t_vec * (1.0 / scale);
    
    // 步骤5：正交化R矩阵（使用Gram-Schmidt）
    r1 = r1.normalized();
    r2 = (r2 - r1 * r2.dot(r1)).normalized();
    r3 = r1.cross(r2);  // 叉乘保证正交
    
    // 检查并调整方向（确保det(R) > 0）
    if (r1.cross(r2).dot(r3) < 0) {
        r1 = r1 * (-1);
        r2 = r2 * (-1);
        r3 = r3 * (-1);
        t_vec = t_vec * (-1);
    }
    
    R = Matrix(3, 3);
    R(0, 0) = r1.x; R(0, 1) = r1.y; R(0, 2) = r1.z;
    R(1, 0) = r2.x; R(1, 1) = r2.y; R(1, 2) = r2.z;
    R(2, 0) = r3.x; R(2, 1) = r3.y; R(2, 2) = r3.z;
    
    t = t_vec;
    
    return true;
}

// 验证结果
void verifyResult(const vector<pair<double, double>>& imagePoints,
                  const vector<Vec3>& worldPoints,
                  const Matrix& K,
                  const Matrix& R,
                  const Vec3& t) {
    cout << "========== 验证结果 ==========\n";
    
    double totalError = 0;
    for (size_t i = 0; i < imagePoints.size(); i++) {
        // 投影: x = K * [R|t] * X
        Vec3 X = worldPoints[i];
        double xc = R(0, 0) * X.x + R(0, 1) * X.y + R(0, 2) * X.z + t.x;
        double yc = R(1, 0) * X.x + R(1, 1) * X.y + R(1, 2) * X.z + t.y;
        double zc = R(2, 0) * X.x + R(2, 1) * X.y + R(2, 2) * X.z + t.z;
        
        double u_proj = (K(0, 0) * xc + K(0, 2) * zc) / zc;
        double v_proj = (K(1, 1) * yc + K(1, 2) * zc) / zc;
        
        double error = sqrt(pow(u_proj - imagePoints[i].first, 2) + 
                           pow(v_proj - imagePoints[i].second, 2));
        totalError += error;
        
        cout << "点" << i << ": 原始(" << imagePoints[i].first << ", " << imagePoints[i].second << ") "
             << "重投影(" << u_proj << ", " << v_proj << ") "
             << "误差=" << error << "\n";
    }
    
    cout << "平均重投影误差: " << totalError / imagePoints.size() << " 像素\n\n";
}

int main() {
    cout << "========== 无人机摄影测量：PnP相机参数求解 (完整版) ==========\n\n";
    
    // 测试数据
    vector<Vec3> worldPoints = {
        Vec3(0, 0, 0), Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(1, 1, 0),
        Vec3(0, 0, 1), Vec3(1, 0, 1), Vec3(0, 1, 1), Vec3(1, 1, 1)
    };
    
    Matrix K(3, 3);
    K(0, 0) = 800; K(0, 1) = 0; K(0, 2) = 320;
    K(1, 0) = 0; K(1, 1) = 800; K(1, 2) = 240;
    K(2, 0) = 0; K(2, 1) = 0; K(2, 2) = 1;
    K.print("相机内参矩阵 K");
    
    // 真实外参
    double angle = 30.0 * M_PI / 180.0;
    Matrix R_true(3, 3);
    R_true(0, 0) = cos(angle); R_true(0, 1) = 0; R_true(0, 2) = sin(angle);
    R_true(1, 0) = 0; R_true(1, 1) = 1; R_true(1, 2) = 0;
    R_true(2, 0) = -sin(angle); R_true(2, 1) = 0; R_true(2, 2) = cos(angle);
    Vec3 t_true(-2.0, 1.0, 5.0);
    
    R_true.print("真实旋转矩阵 R");
    cout << "真实平移向量 T: (" << t_true.x << ", " << t_true.y << ", " << t_true.z << ")\n\n";
    
    // 生成图像点
    vector<pair<double, double>> imagePoints;
    cout << "生成的图像点:\n";
    for (const auto& wp : worldPoints) {
        double xc = R_true(0, 0) * wp.x + R_true(0, 1) * wp.y + R_true(0, 2) * wp.z + t_true.x;
        double yc = R_true(1, 0) * wp.x + R_true(1, 1) * wp.y + R_true(1, 2) * wp.z + t_true.y;
        double zc = R_true(2, 0) * wp.x + R_true(2, 1) * wp.y + R_true(2, 2) * wp.z + t_true.z;
        
        double u = (K(0, 0) * xc + K(0, 2) * zc) / zc;
        double v = (K(1, 1) * yc + K(1, 2) * zc) / zc;
        imagePoints.push_back({u, v});
        cout << "  (" << u << ", " << v << ")\n";
    }
    cout << "\n";
    
    // 求解PnP
    Matrix R_solved;
    Vec3 t_solved;
    cout << "========== 开始求解PnP ==========\n";
    if (solvePnP(imagePoints, worldPoints, K, R_solved, t_solved)) {
        cout << "求解成功！\n\n";
        
        R_solved.print("求解得到的旋转矩阵 R");
        cout << "求解得到的平移向量 T: (" << t_solved.x << ", " << t_solved.y << ", " << t_solved.z << ")\n\n";
        
        verifyResult(imagePoints, worldPoints, K, R_solved, t_solved);
        
        cout << "========== 与真实值对比 ==========\n";
        cout << "旋转矩阵差异:\n";
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                cout << setw(12) << R_solved(i, j) - R_true(i, j) << " ";
            }
            cout << "\n";
        }
        cout << "\n平移向量差异: ("
             << t_solved.x - t_true.x << ", "
             << t_solved.y - t_true.y << ", "
             << t_solved.z - t_true.z << ")\n";
    } else {
        cout << "求解失败!\n";
    }
    
    return 0;
}
