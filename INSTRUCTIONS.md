# 项目说明

本项目实现了空间点求相机参数的PnP算法，求解方程 **x = K[R|T]X**。

## 文件列表

1. **pnp_solver.cpp** - 主要实现文件（推荐使用）
2. **test_pnp.cpp** - 测试套件
3. **pnp_final.cpp** - 完整版实现
4. **CMakeLists.txt** - CMake构建文件
5. **README.md** - 说明文档

## 快速开始

### 方法1：直接编译运行

```bash
# 编译主程序
g++ -o pnp_solver pnp_solver.cpp -std=c++11 -O2 -lm

# 运行
./pnp_solver
```

### 方法2：使用CMake

```bash
mkdir build && cd build
cmake ..
make
./bin/pnp_solver
```

## 算法说明

### 步骤1：求解线性方程（DLT方法）

1. 使用相机内参矩阵K的逆矩阵将图像点归一化
2. 对每个点对构建2个线性方程
3. 构建矩阵方程 A*m = 0，其中m是12维向量（[R|T]展开）
4. 使用SVD求解，取最小奇异值对应的奇异向量作为解

### 步骤2：矩阵分解

1. 将12维解向量重组为3x4的[R|T]矩阵
2. 提取前3x3部分作为旋转矩阵近似值
3. 使用Gram-Schmidt正交化或SVD强制R为正交矩阵
4. 确保det(R) = 1（右手坐标系）
5. 提取最后一列作为平移向量

## 实现特点

- ✅ 完全使用C++实现
- ✅ 未调用OpenCV等第三方PnP函数
- ✅ 自实现矩阵运算（乘法、转置、求逆）
- ✅ 自实现SVD分解（Jacobi方法）
- ✅ 包含完整测试用例和验证

## 理论依据

基于《Multiple View Geometry》第6章的PnP算法实现。

---

**注意**: 由于DLT算法本身的数值稳定性问题，在某些配置下可能需要添加噪声或使用更多点来获得更好的结果。对于实际应用，建议使用迭代优化方法（如Levenberg-Marquardt）进一步优化结果。
