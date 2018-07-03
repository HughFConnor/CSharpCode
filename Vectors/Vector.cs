/*
 *  Name    Vectors
 *  Date    2 July 2018(
 *  Author  Hugh Connor
 *  Purpose
 *      Computer graphics homogeneous vector routines.
 *  Member(s):  vec (vector values)
 *  Methods:
 *      Initialization:
 *          Vector()            //  0
 *          Vector(Vector Mx)   //  Copy another Vector.
 *          Vector(a, b, c, d)  //  Create a Vector from values.
 *          Initialize(a, b, c, d)  //  Initialize a Vector from values.
 *          GetElement(i)       //  get vec[i]
 *          SetElement(i, v)    //  vec[i] = v
 *      Operations:
 *          Dot(Vector vDot)    //  Returns the dot product of vDot and this vector.
 *          Cross(Vector vDot)  //  Returns the vector cross product of vectors a and b.
 *      Print:
 *          Display()           //  Return a formatted string of this matrix.
 * 
 * YYYYMMDD Who     What
 * -------- ---     -----------------------------------------------
 * 20180702 HFC     Initial deployment
 */
using System;

namespace Vectors
{
    public class Vector
    {
        private double [] vec = new double[4] {0.0, 0.0, 0.0, 0.0};

        public Vector()
        {
        }

        public Vector(Vector vX)
        {
            for (int i = 0; i < 4; i++)
                this.vec[i] = vX.vec[i];
        }

        public Vector(double a, double b, double c, double d)
        {
            this.vec[0] = a;
            this.vec[1] = b;
            this.vec[2] = c;
            this.vec[3] = d;
        }

        public void Initialize(double a, double b, double c, double d)
        {
            this.vec[0] = a;
            this.vec[1] = b;
            this.vec[2] = c;
            this.vec[3] = d;
        }

        public double GetElement(int i)
        {
            return vec[i];
        }

        public void SetElement(int i, double v)
        {
            vec[i] = v;
        }

        public double Dot(Vector vDot)
        {
            //  non-homogeneous vector dot product 
            return this.vec[0] * vDot.vec[0] + 
                    this.vec[1] * vDot.vec[1] +
                    this.vec[2] * vDot.vec[2] +
                    this.vec[3] * vDot.vec[3];
        }

        public void Cross(Vector a, Vector b)
        {
            //  Create copies of the input, then calculate the cross product.

            Vector A = new Vector(a);
            Vector B = new Vector(b);

            this.vec[0] = A.vec[1] * B.vec[2] - A.vec[2] * B.vec[1];
            this.vec[1] = A.vec[2] * B.vec[0] - A.vec[0] * B.vec[2];
            this.vec[2] = A.vec[0] * B.vec[1] - A.vec[1] * B.vec[0];
            this.vec[3] = 0.0;
        }

        public string Display()
        {
            string xx = "";

            for (int i = 0; i < 4; i++)
                xx += this.vec[i] + " ";

            return xx;
        }
    }
}
