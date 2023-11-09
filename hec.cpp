#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <unsupported/Eigen/KroneckerProduct>
#include <utility>
#include <vector>


Eigen::Affine3d get_random_transform_matrix(){
  // Create a random rotation matrix
  Eigen::Matrix3d rotation_matrix = Eigen::Quaterniond::UnitRandom().toRotationMatrix();

  // Create a random translation vector
  Eigen::Vector3d translation_vector;
  translation_vector.setRandom();

  Eigen::Affine3d transform_matrix = Eigen::Translation3d(translation_vector) * Eigen::Affine3d(rotation_matrix);
  return transform_matrix;
}

// Function to extract the rotation matrix from a transformation matrix
Eigen::Matrix3d extractRotation(const Eigen::Affine3d& transform) {
  return transform.rotation();
}

// Function to extract the translation vector from a transformation matrix
Eigen::Vector3d extractTranslation(const Eigen::Affine3d& transform) {
  return transform.translation();
}


Eigen::Affine3d create_affine_matrix(Eigen::Matrix3d R, Eigen::Vector3d t){
  Eigen::Affine3d transform = Eigen::Translation3d(t) * Eigen::Affine3d(R);
  return transform;
}


std::pair<Eigen::Affine3d, Eigen::Affine3d> li_wang_wu_2010(
                                                            const std::vector<Eigen::Affine3d>& A_vec,
                                                            const std::vector<Eigen::Affine3d>& B_vec,
                                                            const int &n)
{
  Eigen::MatrixXd merged_lhs(n*12, 24);
  Eigen::MatrixXd merged_rhs(n*12, 1);

  for (int i = 0; i < n; ++i) {
    const Eigen::Affine3d A = A_vec[i];
    const Eigen::Affine3d B = B_vec[i];

    Eigen::Vector3d t_A = A.translation();
    Eigen::Matrix3d R_A = A.rotation();
    Eigen::Vector3d t_B = B.translation();
    Eigen::Matrix3d R_B = B.rotation();

    Eigen::MatrixXd I_3 = Eigen::Matrix3d::Identity(3, 3);

    Eigen::MatrixXd zero_9x9 = Eigen::MatrixXd::Zero(9, 9);
    Eigen::MatrixXd zero_9x3 = Eigen::MatrixXd::Zero(9, 3);
    Eigen::MatrixXd zero_3x9 = Eigen::MatrixXd::Zero(3, 9);
    Eigen::MatrixXd zero_9x1 = Eigen::MatrixXd::Zero(9, 1);

    Eigen::MatrixXd R_A_I_3 = Eigen::kroneckerProduct(R_A, I_3);
    Eigen::MatrixXd I_3_R_B_T = Eigen::kroneckerProduct(I_3, R_B.transpose());
    Eigen::MatrixXd I_3_t_B_T = Eigen::kroneckerProduct(I_3, t_B.transpose());

    // std::cout << "R_A_I_3:\n" << R_A_I_3.matrix() << std::endl;
    // std::cout << "I_3_R_B_T:\n" << I_3_R_B_T.matrix() << std::endl;
    // std::cout << "I_3_t_B_T:\n" << I_3_t_B_T.matrix() << std::endl;

    Eigen::MatrixXd single_upper_lhs(9, 24);
    single_upper_lhs << R_A_I_3, -I_3_R_B_T, zero_9x3, zero_9x3;
    Eigen::MatrixXd single_lower_lhs(3, 24);
    single_lower_lhs << zero_3x9, I_3_t_B_T, -R_A, I_3;

    // std::cout << "single_upper_lhs:\n" << single_upper_lhs.matrix() << std::endl;
    // std::cout << "single_lower_lhs:\n" << single_lower_lhs.matrix() << std::endl;

    Eigen::MatrixXd single_lhs(12, 24);
    single_lhs << single_upper_lhs, single_lower_lhs;
    Eigen::MatrixXd single_rhs(12, 1);
    single_rhs << zero_9x1, t_A;

    merged_lhs.block(i*n, 0, 12, 24) = single_lhs;
    merged_rhs.block(i*n, 0, 12, 1) = single_rhs;
  }
  // std::cout << "merged rhs:\n" << merged_rhs.matrix() << std::endl;

  Eigen::MatrixXd x_z_tX_tZ = (merged_lhs.transpose() * merged_lhs).inverse() * merged_lhs.transpose() * merged_rhs;
  // std::cout << "outputs:\n" << x_z_tX_tZ.matrix() << std::endl;

  Eigen::Matrix3d R_X = x_z_tX_tZ.block(0, 0, 9, 1).reshaped(3,3);
  Eigen::Matrix3d R_Z = x_z_tX_tZ.block(9, 0, 9, 1).reshaped(3,3);
  Eigen::Vector3d t_X = x_z_tX_tZ.block(18, 0, 3, 1);
  Eigen::Vector3d t_Z = x_z_tX_tZ.block(21, 0, 3, 1);

  Eigen::Affine3d X = create_affine_matrix(R_X, t_X);
  Eigen::Affine3d Z = create_affine_matrix(R_Z, t_Z);

  return std::make_pair(X, Z);
}


int main() {
  // Seed the random number generator (optional)
  srand(time(0));

  // Create a transformation matrix

  Eigen::Affine3d X = get_random_transform_matrix();
  std::cout << "X:\n" << X.matrix() << std::endl;

  Eigen::Affine3d Z = get_random_transform_matrix();
  std::cout << "Z:\n" << Z.matrix() << std::endl;

  std::vector<Eigen::Affine3d> A_vec;
  std::vector<Eigen::Affine3d> B_vec;

  // Needed at least n=3 measurement according to:
  // Solving the Robot-World/HandEye Calibration Problem Using the Kronecker Product, Mili Shah, 2013
  int n = 10;
  for (int i=0; i<n; i++) {
    Eigen::Affine3d A = get_random_transform_matrix();
    // std::cout << "A:\n" << A.matrix() << std::endl;

    // A*X = Z*B
    Eigen::Affine3d B = Z.inverse() * A * X;
    // std::cout << "B:\n" << B.matrix() << std::endl;

    // std::cout << "A*X:\n" << A*X.matrix() << std::endl;
    // std::cout << "Z*B:\n" << Z*B.matrix() << std::endl;
    // std::cout << ((A*X).isApprox(Z*B)) << std::endl;

    A_vec.push_back(A);
    B_vec.push_back(B);
  }

  std::pair<Eigen::Affine3d, Eigen::Affine3d> XZ_lww = li_wang_wu_2010(A_vec, B_vec, n);
  Eigen::Affine3d X_lww = XZ_lww.first;
  Eigen::Affine3d Z_lww = XZ_lww.second;
  std::cout << "X:\n" << X_lww.matrix() << std::endl;
  std::cout << "Z:\n" << Z_lww.matrix() << std::endl;
  return 0;
}
