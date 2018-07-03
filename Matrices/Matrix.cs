/*
 *  Name    Matrices
 *  Date    1 July 2018
 *  Author  Hugh Connor
 *  Purpose
 *      Computer graphics homogeneous matrix routines.
 *  Use a vector to implement matrix-vector multiplication.
 *  Member(s):  mat (matrix values)
 *  Methods:
 *      Initialization:
 *          Matrix()            //  Identity
 *          Matrix(Matrix Mx)   //  Copy another Matrix.
 *          GetElement(r, c)    //  get mat[r, c]
 *          SetElement(r, c, v) //  mat[r, c] = v
 *          Identity()          //  Reset to an identity matrix.
 *      Transformations:
 *          Translate(d, v)     //  translate in direction d for a distance v.
 *          Translate(dx, dy, dz)   //  Translate in all three canonical directions at once.
 *          Rotate(d, v)        //  Rotate about canonical axis d for v radians.
 *          Scale(d, v)         //  Scale in the canonical direction d for size v.
 *          Scale(dx, dy, dz)   //  Scale in all three canonical directions at once.
 *          Shear(c, s, v)      //  Shear the given coordinate s by shearing coordinate s by a factor v.
 *          Shear(d, v1, v2)    //  Shear in the direction d with the other two values v1, v2.
 *      Operations:
 *          Multiply(Matrix mMult)  //  Post-multiply the given matrix to this matrix.
 *          MultiplyScalar(s)   //  Multiply each element by th scalar s.
 *          MultiplyVector(Vector v)    //  Post-multiply this matrix with the given vector, returning a vector.
 *          Transpose()         //  Transpose this matrix across the main diagonal.
 *          Determinant()       //  Return the determinant of this 4x4 matrix.  Recurse in method Det.
 *          Det(s, v)           //  Calculate the determinant of the s x s matrix, recursively call itself for smaller matrices.
 *          Trace()             //  Calculate the trace of this matrix.
 *          Adjoint(Matrix mA)  //  Calculate this matrix as the adjoint of the given matrix.
 *          Inverse()           //  Calculate the inverse of this matrix or its adjoint.
 *      Print:
 *          Display()           //  Return a formatted string of this matrix.
 *          
 * 
 * YYYYMMDD Who     What
 * -------- ---     -----------------------------------------------
 * 20180701 HFC     Initial deployment
 */

using System;

using Vectors;

namespace Matrices
{
    public class Matrix
    {
        //  Initialize as an identity matrix.

        private double[,] mat = new double [4,4] { 
                                    {1.0, 0.0, 0.0, 0.0}, 
                                    {0.0, 1.0, 0.0, 0.0},
                                    {0.0, 0.0, 1.0, 0.0},
                                    {0.0, 0.0, 0.0, 1.0} };

        public Matrix()
        {
        }

        public Matrix(Matrix mX)
        {
            for (int i = 0; i < 4; i++)
                for (int j = 0; j < 4; j++)
                    this.mat[i, j] = mX.mat[i, j];
        }

        public double GetElement(int r, int c)
        {
            return this.mat[r, c];
        }

        public void SetElement(int r, int c, double v)
        {
            this.mat[r, c] = v;
        }

        public void Identity()
        {
            for (int i = 0; i < 4; i++)
                for (int j = 0; j < 4; j++)
                    this.mat[i, j] = (i == j) ? 1.0 : 0.0;
        }

        public void Translate(char d, double v)
        {
            switch (d)
            {
                case 'x':
                case 'X':
                    this.mat[0, 3] = v;
                    break;
                case 'y':
                case 'Y':
                    this.mat[1, 3] = v;
                    break;
                case 'z':
                case 'Z':
                    this.mat[2, 3] = v;
                    break;
                default:
                    break;
            }
        }

        public void Translate(double dx, double dy, double dz)
        {
            this.mat[0, 3] = dx;
            this.mat[1, 3] = dy;
            this.mat[2, 3] = dz;
        }

        public void Rotate(char d, double v)
        {
            double c = Math.Cos(v);
            double s = Math.Sin(v);

            switch (d)
            {
                case 'x':
                case 'X':
                    this.mat[1, 1] = c;
                    this.mat[1, 2] = -s;
                    this.mat[2, 1] = s;
                    this.mat[2, 2] = c;
                    break;
                case 'y':
                case 'Y':
                    this.mat[0, 0] = c;
                    this.mat[0, 2] = s;
                    this.mat[2, 0] = -s;
                    this.mat[2, 2] = c;
                    break;
                case 'z':
                case 'Z':
                    this.mat[0, 0] = c;
                    this.mat[0, 1] = -s;
                    this.mat[1, 0] = s;
                    this.mat[1, 1] = c;
                    break;
                default:
                    break;
            }
        }

        public void Scale(char d, double v)
        {
            switch (d)
            {
                case 'x':
                case 'X':
                    this.mat[0, 0] = v;
                    break;
                case 'y':
                case 'Y':
                    this.mat[1, 1] = v;
                    break;
                case 'z':
                case 'Z':
                    this.mat[2, 2] = v;
                    break;
                default:
                    break;
            }
        }

        public void Scale(double dx, double dy, double dz)
        {
            this.mat[0, 0] = dx;
            this.mat[1, 1] = dy;
            this.mat[2, 2] = dz;
        }

        public void Shear(char c, char s, double v)
        {
            //  Shear the coordinate c using the coordinate s.

            if ("xyz".IndexOf(c) >= 0 && "xyz".IndexOf(s) >= 0 && c != s)
            {
                int row = (c == 'x') ? 0 : (c == 'y') ? 1 : 2;
                int col = (s == 'x') ? 0 : (s == 'y') ? 1 : 2;
                this.mat[row, col] = v;
            }
        }
        public void Shear(char d, double v1, double v2)
        {
            //  Shear the coordinate d using the given both shearing values.

            switch (d)
            {
                case 'x':
                case 'X':
                    this.mat[1, 0] = v1;
                    this.mat[2, 0] = v2;
                    break;
                case 'y':
                case 'Y':
                    this.mat[0, 1] = v1;
                    this.mat[2, 1] = v2;
                    break;
                case 'z':
                case 'Z':
                    this.mat[0, 2] = v1;
                    this.mat[1, 2] = v2;
                    break;
                default:
                    break;
            }
        }

        public void Multiply(Matrix mMult)
        {
            //  Post-multiply this by mMult.

            Matrix temp = new Matrix(this); //  Create a working copy.

            for (int i = 0; i < 4; i++)
                for (int j = 0; j < 4; j++)
                {
                    double v = 0.0;
                    for (int k = 0; k < 4; k++)
                        v += mMult.mat[k, i] * temp.mat[j, k];
                    this.mat[j, i] = v;
                }
        }

        public void MultiplyScalar(double s)
        {
            for (int i = 0; i < 4; i++)
                for (int j = 0; j < 4; j++)
                    this.mat[i, j] *= s;
        }

        public Vector MultiplyVector(Vector v)
        {
            Vector mv = new Vector();

            for (int i = 0; i < 4; i++)
            {
                double vi = 0;
                for (int j = 0; j < 4; j++)
                    vi += v.GetElement(j) * this.mat[i, j];
                mv.SetElement(i, vi);
            }
            return mv;
        }

        public void Transpose()
        {
            for (int i = 0; i < 4; i++)
                for (int j = 0; j < 4; j++)
                    if (i != j)
                    {
                        double temp = this.mat[i, j];
                        this.mat[i, j] = this.mat[j, i];
                        this.mat[j, i] = temp;
                    }
        }

        public double Determinant()
        {
            double[] v = new double[16];
            int k = 0;
            for (int i = 0; i < 4; i++)
                for (int j = 0; j < 4; j++)
                    v[k++] = this.mat[i, j];

            double val = Det(4, v);

            return val;
        }

        private double Det(int s, double[] v)
        {
            double val = 0.0;
            if (s == 2)
                val = v[0] * v[3] - v[1] * v[2];
            else if (s == 1)
                val = v[0];
            else
            {
                //  Allocate space for the next smaller determinant.
                double[] vs = new double[(s - 1) * (s - 1)];
                for (int i = 0; i < s; i++) //  for each value in row[0]...
                    if (v[i] != 0.0)    //  Skip zero factors.
                    {
                        int k = 0;  //  column index to skip the ith column
                        int n = 0;  //  next determinant array index
                        for (int j = s; j < (s * s); j++)
                        {
                            if (k != i) //  the ith column is skipped
                            {
                                vs[n++] = v[j];
                            }
                            k++;
                            if (k == s) k = 0;  //  Reset column index
                        }
                        double det = this.Det(s - 1, vs);   //  recurse to next smaller determinant
                        if ((i & 1) == 0)   //  Alternate row[0] +/-
                            val += v[i] * det;
                        else
                            val -= v[i] * det;
                    }
            }
            return val;
        }

        public double Trace()
        {
            return this.mat[0, 0] + this.mat[1, 1] + this.mat[2, 2] + this.mat[3, 3]; 
        }

        public void Adjoint(Matrix mA)
        {
            Matrix temp = new Matrix(this); //  Create a working copy.

            double[] vs = new double[9];    //  3x3 matrices for each adjoint minor.

            for (int i = 0; i < 4; i++)     //  row
                for (int j = 0; j < 4; j++) //  column
                {
                    int n = 0;  //  sub-matrix insert counter
                    for (int k = 0; k < 4; k++)     //  row
                        for (int l = 0; l < 4; l++) //  column
                            if (k != i && l != j)
                                vs[n++] = temp.mat[k, l];
                    double det = temp.Det(3, vs);
                    if (((i + j) & 1) == 0)
                        this.mat[j, i] = det;
                    else
                        this.mat[j, i] = -det;
                }
        }

        public void Inverse()
        {
            const double epsilon = 1.00e-6;
            double det = this.Determinant();
            this.Adjoint(this);
            if (Math.Abs(det) > epsilon)
                this.MultiplyScalar(1.00 / det);
        }

        public string Display()
        {
            string xx = "";

            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 4; j++)
                    xx += this.mat[i, j].ToString() + " ";
                xx += "\n";
            }
            return xx;
        }
    }
}

