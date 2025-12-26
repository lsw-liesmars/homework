/**
 * 无人机摄影测量作业：空间点求相机参数
 * 求解方程：x = K[R|T]X
 * 其中 x(图像点), K(内参矩阵), X(空间点) 已知
 * 求解 R(旋转矩阵) 和 T(平移向量)
 * 
 * 算法步骤：
 * 1. 使用DLT方法求解线性方程，得到投影矩阵M = K[R|T]
 * 2. 通过矩阵分解从M中提取R和T
 */

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

// 矩阵类定义
class Matrix {
public:
    int rows, cols;
    vector<vector<double>> data;
    
    // 默认构造函数
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
    
    // 矩阵乘法
    Matrix operator*(const Matrix& other) const {
        if (cols != other.rows) {
            throw runtime_error("矩阵维度不匹配！");
        }
        Matrix result(rows, other.cols);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < other.cols; j++) {
                for (int k = 0; k < cols; k++) {
                    result.data[i][j] += data[i][k] * other.data[k][j];
                }
            }
        }
        return result;
    }
    
    // 转置
    Matrix transpose() const {
        Matrix result(cols, rows);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                result.data[j][i] = data[i][j];
            }
        }
        return result;
    }
    
    // 矩阵求逆（使用高斯消元法）
    Matrix inverse() const {
        if (rows != cols) {
            throw runtime_error("只能对方阵求逆！");
        }
        int n = rows;
        Matrix augmented(n, 2 * n);
        
        // 构造增广矩阵 [A | I]
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                augmented.data[i][j] = data[i][j];
            }
            augmented.data[i][i + n] = 1.0;
        }
        
        // 高斯-约当消元
        for (int i = 0; i < n; i++) {
            // 选主元
            int maxRow = i;
            for (int k = i + 1; k < n; k++) {
                if (abs(augmented.data[k][i]) > abs(augmented.data[maxRow][i])) {
                    maxRow = k;
                }
            }
            swap(augmented.data[i], augmented.data[maxRow]);
            
            // 归一化
            double pivot = augmented.data[i][i];
            if (abs(pivot) < 1e-10) {
                throw runtime_error("矩阵不可逆！");
            }
            for (int j = 0; j < 2 * n; j++) {
                augmented.data[i][j] /= pivot;
            }
            
            // 消元
            for (int k = 0; k < n; k++) {
                if (k != i) {
                    double factor = augmented.data[k][i];
                    for (int j = 0; j < 2 * n; j++) {
                        augmented.data[k][j] -= factor * augmented.data[i][j];
                    }
                }
            }
        }
        
        // 提取逆矩阵
        Matrix result(n, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                result.data[i][j] = augmented.data[i][j + n];
            }
        }
        return result;
    }
    
    // 打印矩阵
    void print(const string& name = "") const {
        if (!name.empty()) {
            cout << name << ":" << endl;
        }
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                cout << setw(12) << setprecision(6) << data[i][j] << " ";
            }
            cout << endl;
        }
        cout << endl;
    }
};

// SVD分解（使用Jacobi方法求解对称矩阵的特征值和特征向量）
class SVD {
public:
    Matrix U, S, V;
    
    SVD(const Matrix& A) : U(A.rows, A.rows), S(A.rows, A.cols), V(A.cols, A.cols) {
        // 计算 A^T * A
        Matrix AtA = A.transpose() * A;
        
        // 对 A^T * A 进行特征分解得到 V 和奇异值平方
        eigenDecomposition(AtA, V, S);
        
        // 计算奇异值（对角矩阵S的对角元素开方）
        for (int i = 0; i < S.rows && i < S.cols; i++) {
            S.data[i][i] = sqrt(max(0.0, S.data[i][i]));
        }
        
        // 计算 U = A * V * S^-1
        Matrix A_V = A * V;
        for (int i = 0; i < U.rows; i++) {
            for (int j = 0; j < U.cols && j < S.cols; j++) {
                if (S.data[j][j] > 1e-10) {
                    U.data[i][j] = A_V.data[i][j] / S.data[j][j];
                } else {
                    U.data[i][j] = 0;
                }
            }
        }
    }
    
private:
    // Jacobi方法求对称矩阵特征值和特征向量
    void eigenDecomposition(const Matrix& A, Matrix& eigenvectors, Matrix& eigenvalues) {
        int n = A.rows;
        Matrix B = A;
        eigenvectors = Matrix(n, n);
        
        // 初始化为单位矩阵
        for (int i = 0; i < n; i++) {
            eigenvectors.data[i][i] = 1.0;
        }
        
        // Jacobi迭代
        const int maxIter = 100;
        const double tolerance = 1e-10;
        
        for (int iter = 0; iter < maxIter; iter++) {
            // 找最大非对角元素
            double maxVal = 0;
            int p = 0, q = 1;
            for (int i = 0; i < n; i++) {
                for (int j = i + 1; j < n; j++) {
                    if (abs(B.data[i][j]) > maxVal) {
                        maxVal = abs(B.data[i][j]);
                        p = i;
                        q = j;
                    }
                }
            }
            
            if (maxVal < tolerance) break;
            
            // 计算旋转角度
            double theta;
            if (abs(B.data[p][p] - B.data[q][q]) < 1e-10) {
                theta = M_PI / 4.0;
            } else {
                theta = 0.5 * atan(2.0 * B.data[p][q] / (B.data[p][p] - B.data[q][q]));
            }
            
            double c = cos(theta);
            double s = sin(theta);
            
            // 构造旋转矩阵并应用
            Matrix R = Matrix(n, n);
            for (int i = 0; i < n; i++) {
                R.data[i][i] = 1.0;
            }
            R.data[p][p] = c;
            R.data[q][q] = c;
            R.data[p][q] = -s;
            R.data[q][p] = s;
            
            B = R.transpose() * B * R;
            eigenvectors = eigenvectors * R;
        }
        
        // 提取特征值（对角元素）
        eigenvalues = Matrix(n, n);
        for (int i = 0; i < n; i++) {
            eigenvalues.data[i][i] = B.data[i][i];
        }
        
        // 按特征值从大到小排序
        for (int i = 0; i < n - 1; i++) {
            for (int j = i + 1; j < n; j++) {
                if (eigenvalues.data[j][j] > eigenvalues.data[i][i]) {
                    swap(eigenvalues.data[i][i], eigenvalues.data[j][j]);
                    for (int k = 0; k < n; k++) {
                        swap(eigenvectors.data[k][i], eigenvectors.data[k][j]);
                    }
                }
            }
        }
    }
};

// PnP求解器
class PnPSolver {
public:
    /**
     * 求解相机外参数
     * @param image_points 图像点坐标 (u, v, 1)，每行一个点
     * @param world_points 世界坐标系中的点 (X, Y, Z, 1)，每行一个点
     * @param K 相机内参矩阵 3x3
     * @param R 输出：旋转矩阵 3x3
     * @param T 输出：平移向量 3x1
     */
    static bool solve(const vector<vector<double>>& image_points,
                     const vector<vector<double>>& world_points,
                     const Matrix& K,
                     Matrix& R,
                     Matrix& T) {
        
        int n = image_points.size();
        if (n < 6) {
            cerr << "至少需要6个点才能求解PnP问题！" << endl;
            return false;
        }
        
        // 步骤1：将图像点归一化（去除内参影响）
        Matrix K_inv = K.inverse();
        vector<vector<double>> normalized_points(n, vector<double>(3));
        
        for (int i = 0; i < n; i++) {
            Matrix img_pt(3, 1);
            img_pt.data[0][0] = image_points[i][0];
            img_pt.data[1][0] = image_points[i][1];
            img_pt.data[2][0] = 1.0;
            
            Matrix norm_pt = K_inv * img_pt;
            normalized_points[i][0] = norm_pt.data[0][0];
            normalized_points[i][1] = norm_pt.data[1][0];
            normalized_points[i][2] = norm_pt.data[2][0];
        }
        
        // 步骤2：构建DLT线性方程组 A * m = 0
        // 其中 m 是 [R|T] 矩阵展开的12维向量
        Matrix A(2 * n, 12);
        
        for (int i = 0; i < n; i++) {
            double x = normalized_points[i][0];
            double y = normalized_points[i][1];
            double X = world_points[i][0];
            double Y = world_points[i][1];
            double Z = world_points[i][2];
            
            // 第一个方程（x方向）
            A.data[2*i][0] = X;
            A.data[2*i][1] = Y;
            A.data[2*i][2] = Z;
            A.data[2*i][3] = 1;
            A.data[2*i][4] = 0;
            A.data[2*i][5] = 0;
            A.data[2*i][6] = 0;
            A.data[2*i][7] = 0;
            A.data[2*i][8] = -x * X;
            A.data[2*i][9] = -x * Y;
            A.data[2*i][10] = -x * Z;
            A.data[2*i][11] = -x;
            
            // 第二个方程（y方向）
            A.data[2*i+1][0] = 0;
            A.data[2*i+1][1] = 0;
            A.data[2*i+1][2] = 0;
            A.data[2*i+1][3] = 0;
            A.data[2*i+1][4] = X;
            A.data[2*i+1][5] = Y;
            A.data[2*i+1][6] = Z;
            A.data[2*i+1][7] = 1;
            A.data[2*i+1][8] = -y * X;
            A.data[2*i+1][9] = -y * Y;
            A.data[2*i+1][10] = -y * Z;
            A.data[2*i+1][11] = -y;
        }
        
        // 步骤3：使用SVD求解 A * m = 0（m是A的右奇异向量，对应最小奇异值）
        Matrix AtA = A.transpose() * A;
        SVD svd(AtA);
        
        // 最小特征值对应的特征向量（V的最后一列）
        Matrix M(3, 4);
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 4; j++) {
                M.data[i][j] = svd.V.data[i * 4 + j][11];  // 最后一列
            }
        }
        
        // 步骤4：从M中提取R和T
        // M = [R_raw | T_raw]
        Matrix R_raw(3, 3);
        Matrix T_raw(3, 1);

        for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                        R_raw.data[i][j] = M.data[i][j];
                }
                T_raw.data[i][0] = M.data[i][3];
        }

        // 步骤5：正交化与尺度恢复
        // 对 R_raw 进行 SVD 分解: R_raw = U * S * V^T
        SVD svd_r(R_raw);
        Matrix U = svd_r.U;
        Matrix V = svd_r.V;
        Matrix R_ortho = U * V.transpose(); // 得到最接近的旋转矩阵

        // 恢复尺度，通过比较 R_raw 和 R_ortho 的范数来计算尺度
        double norm_raw = 0, norm_ortho = 0;
        for(int i=0; i<3; ++i) 
                for(int j=0; j<3; ++j) {
                        norm_raw += R_raw.data[i][j] * R_raw.data[i][j];
                        norm_ortho += R_ortho.data[i][j] * R_ortho.data[i][j]; // 其实这里是3
                }
        double scale = sqrt(norm_raw / norm_ortho);
        Matrix T_ortho(3, 1);
        // T 也需要除以尺度
        for(int i=0; i<3; ++i) T_ortho.data[i][0] = T_raw.data[i][0] / scale;

        // 符号修正
        if (determinant3x3(R_ortho) < 0) {
                // 整个 R 取反变成 +1
                for(int i=0; i<3; ++i)
                        for(int j=0; j<3; ++j)
                                R_ortho.data[i][j] = -R_ortho.data[i][j];
                
                // T 也要跟随取反！(保持 R 和 T 的相对关系)
                for(int i=0; i<3; ++i)
                        T_ortho.data[i][0] = -T_ortho.data[i][0];
        }

        // 手性约束
        // 变换第一个点，看 Z 是否大于 0
        double x = world_points[0][0];
        double y = world_points[0][1];
        double z = world_points[0][2];
        // P_cam = R * P_world + T
        double z_cam = R_ortho.data[2][0] * x + R_ortho.data[2][1] * y + R_ortho.data[2][2] * z + T_ortho.data[2][0];

        if (z_cam < 0) {
                // 如果点在相机背面，说明 R 和 T 都反了
                for(int i=0; i<3; ++i) {
                        for(int j=0; j<3; ++j) R_ortho.data[i][j] = -R_ortho.data[i][j];
                        T_ortho.data[i][0] = -T_ortho.data[i][0];
                }
        }

        // 赋值回输出变量
        R = R_ortho;
        T = T_ortho;

        return true;
    }
    
private:
    // 计算3x3矩阵的行列式
    static double determinant3x3(const Matrix& M) {
        return M.data[0][0] * (M.data[1][1] * M.data[2][2] - M.data[1][2] * M.data[2][1])
             - M.data[0][1] * (M.data[1][0] * M.data[2][2] - M.data[1][2] * M.data[2][0])
             + M.data[0][2] * (M.data[1][0] * M.data[2][1] - M.data[1][1] * M.data[2][0]);
    }
};

// 辅助函数：验证结果
void verifyResult(const vector<vector<double>>& image_points,
                  const vector<vector<double>>& world_points,
                  const Matrix& K,
                  const Matrix& R,
                  const Matrix& T) {
    cout << "========== 验证结果 ==========" << endl;
    cout << "重投影误差：" << endl;
    
    double total_error = 0;
    int n = image_points.size();
    
    for (int i = 0; i < n; i++) {
        // 世界坐标
        Matrix X(4, 1);
        X.data[0][0] = world_points[i][0];
        X.data[1][0] = world_points[i][1];
        X.data[2][0] = world_points[i][2];
        X.data[3][0] = 1.0;
        
        // 构造 [R|T]
        Matrix RT(3, 4);
        for (int r = 0; r < 3; r++) {
            for (int c = 0; c < 3; c++) {
                RT.data[r][c] = R.data[r][c];
            }
            RT.data[r][3] = T.data[r][0];
        }
        
        // 计算投影：x_proj = K * [R|T] * X
        Matrix x_proj = K * RT * X;
        
        // 归一化
        double u_proj = x_proj.data[0][0] / x_proj.data[2][0];
        double v_proj = x_proj.data[1][0] / x_proj.data[2][0];
        
        // 原始图像点
        double u_orig = image_points[i][0];
        double v_orig = image_points[i][1];
        
        // 计算误差
        double error = sqrt((u_proj - u_orig) * (u_proj - u_orig) + 
                           (v_proj - v_orig) * (v_proj - v_orig));
        total_error += error;
        
        cout << "点 " << i << ": 原始(" << u_orig << ", " << v_orig << ") "
             << "重投影(" << u_proj << ", " << v_proj << ") "
             << "误差=" << error << endl;
    }
    
    cout << "平均重投影误差: " << total_error / n << endl;
    cout << endl;
}

int main() {
    cout << "========== 无人机摄影测量：PnP相机参数求解 ==========" << endl << endl;
    
    // 测试用例：已知的世界坐标点
    vector<vector<double>> world_points = {
        {0, 0, 0},
        {1, 0, 0},
        {0, 1, 0},
        {1, 1, 0},
        {0, 0, 1},
        {1, 0, 1},
        {0, 1, 1},
        {1, 1, 1}
    };
    
    // 已知的相机内参矩阵
    Matrix K(3, 3, {
        800, 0, 320,
        0, 800, 240,
        0, 0, 1
    });
    
    cout << "相机内参矩阵 K:" << endl;
    K.print();
    
    // 真实的旋转和平移（用于生成测试数据）
    double angle = 30.0 * M_PI / 180.0;  // 30度
    Matrix R_true(3, 3, {
        cos(angle), 0, sin(angle),
        0, 1, 0,
        -sin(angle), 0, cos(angle)
    });
    
    Matrix T_true(3, 1, {
        {-2.0},
        {1.0},
        {5.0}
    });
    
    cout << "真实旋转矩阵 R:" << endl;
    R_true.print();
    cout << "真实平移向量 T:" << endl;
    T_true.print();
    
    // 生成图像点坐标
    vector<vector<double>> image_points;
    Matrix RT_true(3, 4);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            RT_true.data[i][j] = R_true.data[i][j];
        }
        RT_true.data[i][3] = T_true.data[i][0];
    }
    
    cout << "生成的图像点坐标:" << endl;
    for (const auto& wp : world_points) {
        Matrix X(4, 1, {wp[0], wp[1], wp[2], 1.0});
        Matrix x = K * RT_true * X;
        
        double u = x.data[0][0] / x.data[2][0];
        double v = x.data[1][0] / x.data[2][0];
        
        image_points.push_back({u, v});
        cout << "  (" << u << ", " << v << ")" << endl;
    }
    cout << endl;
    
    // 求解PnP问题
    cout << "========== 开始求解PnP ==========" << endl;
    Matrix R_solved, T_solved;
    
    if (PnPSolver::solve(image_points, world_points, K, R_solved, T_solved)) {
        cout << "求解成功！" << endl << endl;
        
        cout << "求解得到的旋转矩阵 R:" << endl;
        R_solved.print();
        
        cout << "求解得到的平移向量 T:" << endl;
        T_solved.print();
        
        // 验证结果
        verifyResult(image_points, world_points, K, R_solved, T_solved);
        
        // 比较真实值和求解值
        cout << "========== 与真实值对比 ==========" << endl;
        cout << "旋转矩阵差异:" << endl;
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                cout << setw(12) << R_solved.data[i][j] - R_true.data[i][j] << " ";
            }
            cout << endl;
        }
        cout << endl;
        
        cout << "平移向量差异:" << endl;
        for (int i = 0; i < 3; i++) {
            cout << setw(12) << T_solved.data[i][0] - T_true.data[i][0] << endl;
        }
    } else {
        cout << "求解失败！" << endl;
    }
    
    return 0;
}
