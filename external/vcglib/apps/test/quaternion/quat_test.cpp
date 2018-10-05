#include <iostream>
using namespace std;
#include <vcg/math/quaternion.h>
#include <vcg/math/matrix33.h>


using namespace vcg;



ostream &operator<<(ostream &o, Point3f &q) {
  o.precision(6) ;
  o.setf( ios::fixed, ios::floatfield ) ;
  o << "[" << q[0] << " " << q[1] << " " << q[2]  << "]";
  return o;
}

ostream &operator<<(ostream &o, Quaternionf &q) {
  o.precision(6) ;
  o.setf( ios::fixed, ios::floatfield ) ;
  o << "[" << q[0] << " " << q[1] << " " << q[2] << " " << q[3] << "]";
  return o;
}

ostream &operator<<(ostream &o, Matrix33f &m) {
  o.precision(6) ;
  o.setf( ios::fixed, ios::floatfield ) ;
  for(int i = 0; i < 3; i++) {
    o << "[";
    for(int j = 0; j < 3; j++) {
      o << m[i][j] << " ";
    }
    o << "] ";
  }
  return o;
}
ostream &operator<<(ostream &o, Matrix44f &m) {
  o.precision(6) ;
  o.setf( ios::fixed, ios::floatfield ) ;
  for(int i = 0; i < 4; i++) {
    o << "[";
    for(int j = 0; j < 4; j++) {
      o << m[i][j] << " ";
    }
    o << "] ";
  }
  return o;
}

bool verify(Quaternionf &q){
  cout << "Quaternion: " << q << endl;
  Matrix33f m;
  q.ToMatrix(m);
  cout << "To Matrix: " << m << endl;
  cout << "Row norms: " << m.GetRow(0).Norm() << " "
                        << m.GetRow(1).Norm() << " "
                        << m.GetRow(2).Norm() << endl;
  Point3f p(3, 4, 5);
  Point3f qp = q.Rotate(p);
  Point3f mp = m*p;
  cout << "Rotating p: " << p << endl;
  cout << "q*p = " << qp << " m*p: " << mp << endl;
  q.FromMatrix(m);
  cout << "Back to: " << q << endl << endl;
  return true;
}

bool verify(Matrix33f &m) {
  cout << "Matrix: " << m << endl;
  cout << "Det: " << m.Determinant() << endl;
  cout << "Row norms: " << m.GetRow(0).Norm() << " "
                        << m.GetRow(1).Norm() << " "
                        << m.GetRow(2).Norm() << endl;
  cout << "Column norms: " << m.GetColumn(0).Norm() << " "
                        << m.GetColumn(1).Norm() << " "
                        << m.GetColumn(2).Norm() << endl;
  Matrix33f im = m.transpose() * m;
  im = im*m;
  cout << "Check ortonormality: " <<  im << endl;
  Quaternionf q;
  q.FromMatrix(m);
  cout << "To Quaternion: " << q << endl;

  float alpha = 2*acos(q[0]);
  Point3f axis(q[1], q[2], q[3]);
  cout << "Norm: " << axis.SquaredNorm() + q[0]*q[0] << endl;
  axis.Normalize();
  cout << "angle: " << 2*acos(q[0]) << " Axis: " << axis << endl;

  Point3f p(3, 4, 5);
  Point3f qp = q.Rotate(p);
  Point3f mp = m*p;
  cout << "Rotating p: " << p << endl;
  cout << "q*p = " << qp << " m*p: " << mp << endl;
  q.ToMatrix(m);
  cout << "Back to: " << m << endl<< endl;
}

int main() {

  Quaternionf q(1, 0, 0, 0); /*identity*/
//  cout << "Verify identity: " << endl;
  verify(q);

  q.FromAxis(M_PI/2, Point3f(0, 0, 1));
//  cout << "Verify 90 degrees\n";
  verify(q);

  q.FromAxis(M_PI/4, Point3f(0, 1, 1));
//  cout << "Verify 45 degrees\n";
  verify(q);

  Matrix33f m;
  m[0][0] =  0.70145;  m[0][1] = 0.372035; m[0][2] = 0.607913;
  m[1][0] = -0.628023; m[1][1] = 0.725922; m[1][2] = 0.2804;
  m[2][0] = -0.336978; m[2][1] = -0.57847; m[2][2] =  0.742845;

  cout << "verify matrix: " << endl;
  verify(m);
  verify(m);

  q.FromAxis(0.7, Point3f(-0.20, -0.42, -0.83));
  cout << "verify from axis: " << endl;
  verify(q);

//  Quaternionf q(0.648947, -0.114828, -0.104375, -0.385262);
  Quaternionf iq= Inverse(q);
  //Matrix33f m;

  cout << "matrix: " <<  m << endl;

Point3f p(3, 4, 5);
  //test this matrix:
  cout << "norms: " << m.GetRow(0).Norm() << " " << m.GetRow(1).Norm() << " " << m.GetRow(2).Norm() << endl;
  q.FromMatrix(m);
  q.ToMatrix(m);
  cout << "quaternion: " << q[0] << " " << q[1] << " " << q[2] << " " << q[3] <<endl;
  cout << "matrix: " << endl;
	for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      cout << m[i][j] << " ";
    }
    cout << endl;
  }
  cout << endl;

  cout << "Point: " << p[0] << " " << p[1] << " " << p[2] << endl;

  Point3f r = q.Rotate(p);
  cout << "Rotated by q: " << r[0] << " " << r[1] << " " << r[2] << endl;
  r = m*p;
  cout << "Rotated by m: " << r[0] << " " << r[1] << " " << r[2] << endl;

  q.FromMatrix(m);
  cout << "quaternion: " << q[0] << " " << q[1] << " " << q[2] << " " << q[3] <<endl;
  return 0;
}
