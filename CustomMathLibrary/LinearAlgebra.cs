using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Factorization;

namespace CustomMathLibrary
{
    public class LinearAlgebra
    {
        /// <summary>
        /// For a given plane move a copy of it along its normal. Determine the new value "d" and return the coefficients for the new plane equation.
        /// </summary>
        /// <param name="normal">The normal from the old plane (also the normal from the moved plane)</param>
        /// <param name="pointOnOriginalPlane">A point from the original plane</param>
        /// <param name="factor">The factor on how far the move the copied plane</param>
        /// <returns>The coefficients d, a, b, c of the new plane equation in this order</returns>
        public static float[] CopyPlaneAndMoveAlongNormal(float[] normal, float[] pointOnOriginalPlane, int factor)
        {
            CustomVector3 convertedNormal = CustomVector3Converter.FloatArrayToCustomVector3(normal);
            CustomVector3 convertedPointOnOriginalPlane =
                CustomVector3Converter.FloatArrayToCustomVector3(pointOnOriginalPlane);

            CustomVector3 pointOnCopiedPlane = convertedPointOnOriginalPlane + convertedNormal * factor;
            float dNew = -(convertedNormal.x * pointOnCopiedPlane.x + convertedNormal.y * pointOnCopiedPlane.y +
                           convertedNormal.z * pointOnCopiedPlane.z);

            return new[] {dNew, convertedNormal.x, convertedNormal.y, convertedNormal.z};
        }


        /// <summary>
        ///
        /// Determines the direction vector for a symmetry line
        ///  
        /// </summary>
        /// <param name="normal"></param>
        /// <param name="direction"></param>
        /// <returns>The direction vector</returns>
        public static float[] DetermineDirectionVectorForSymmetryLine(float[] normal, float[] direction)
        {
            /*
                 *	Fuer Symmetrieebene/gerade muss folgendes gelten:
                 *
                 *  (https://de.wikipedia.org/wiki/Orthogonalprojektion)
                 * 
                 * 	Oberkieferebene Parameter (a, b, c, d)
                 * 	g: Richtungsvektor zwischen den beiden Direktionspunkten 
                 * 	r: Richtungsvektor entlang der Symmetrieebene
                 *
                 * 	g * r = 0 (Skalarprodukt)
                 * 	a*r1 + b*r2 + c*r3 = 0
                 *
                 * 	Unterbestimmtes LGS - > QR mit 3x3 Matrix
                 * 
                 */

            Matrix<float> a = Matrix<float>.Build.DenseOfArray(new float[,]
            {
                {direction[0], direction[1], direction[2]},
                {normal[0], normal[1], normal[2]},
                {0, 0, 0}
            });

            return SolveUnderdeterminedLinearEquation(a);

            /*
            
                        OLD VERSION
            
            QR<float> qr = a.QR();
            Matrix<float> q = qr.Q;
            Matrix<float> r = qr.R;

            
            Vector<float> vPortion = Vector<float>.Build.DenseOfArray(new[] {0f, 0, 0});
            Vector<float> wPortion = Vector<float>.Build.DenseOfArray(new[] {0f, 0, 1});

            for (int i = r.RowCount - 2; i >= 0; i--)
            {
                // Von rechts ohne das Diagonalelement j > i
                for (int j = r.ColumnCount - 1; j > i; j--)
                {
                    vPortion[i] -= r[i, j] * vPortion[j];
                    wPortion[i] -= r[i, j] * wPortion[j];
                }

                vPortion[i] /= r[i, i];
                wPortion[i] /= r[i, i];
            }
            
            Vector<float> rV = q.Transpose() * vPortion;
            Vector<float> rW = q.Transpose() * wPortion;


            //  r = (r1,r2,r3) wobei bspw. ri = vi + lambda * wi ist
            float r1 = rV[0] + 1 * rW[0];
            float r2 = rV[1] + 1 * rW[1];
            float r3 = rV[2] + 1 * rW[2];

            CustomVector3 tempResult = new CustomVector3(r1, r2, r3);

            float[] result = CustomVector3Converter.CustomVector3ToFloatArray(tempResult.normalized);

            return result;
            */
        }

        /// <summary>
        ///  Solves an underdetermined linear equation system with 2 equations and 3 unknown.
        /// </summary>
        /// <param name="matrix">The 3x3 matrix, which represents the linear equation system</param>
        /// <returns>One possible solution, where the parameter is set as 1.</returns>
        public static float[] SolveUnderdeterminedLinearEquation(Matrix<float> matrix)
        {
            Matrix<float>[] qr = QRDecomposition(matrix);
            Vector<float>[] portions = GaussianEliminationUnderdeterminedEq(qr[1]);
            Vector<float> rV = qr[0].Transpose() * portions[0];
            Vector<float> rW = qr[0].Transpose() * portions[1];


            //  r = (r1,r2,r3) wobei bspw. ri = vi + lambda * wi ist
            float r1 = rV[0] + 1 * rW[0];
            float r2 = rV[1] + 1 * rW[1];
            float r3 = rV[2] + 1 * rW[2];

            CustomVector3 tempResult = new CustomVector3(r1, r2, r3);
            float[] result = CustomVector3Converter.CustomVector3ToFloatArray(tempResult.normalized);

            return result;
        }

        /// <summary>
        /// Solves a given 3x3 triangular matrix using the gaussian elimination.
        /// The matrix has to represent an underdetermined linear equation system. The last row has to be filled with zeros.
        /// </summary>
        /// <param name="rMatrix">The matrix to solve</param>
        /// <returns>The two parts of the solution</returns>
        public static Vector<float>[] GaussianEliminationUnderdeterminedEq(Matrix<float> rMatrix)
        {
            Vector<float> vPortion = Vector<float>.Build.DenseOfArray(new[] {0f, 0, 0});
            Vector<float> wPortion = Vector<float>.Build.DenseOfArray(new[] {0f, 0, 1});

            for (int i = rMatrix.RowCount - 2; i >= 0; i--)
            {
                // Von rechts ohne das Diagonalelement j > i
                for (int j = rMatrix.ColumnCount - 1; j > i; j--)
                {
                    vPortion[i] -= rMatrix[i, j] * vPortion[j];
                    wPortion[i] -= rMatrix[i, j] * wPortion[j];
                }

                vPortion[i] /= rMatrix[i, i];
                wPortion[i] /= rMatrix[i, i];
            }

            return new[] {vPortion, wPortion};
        }


        /// <summary>
        /// Computes the QR-Decomposition of a given m x n matrix (uses the MathNet.Numerics library).
        /// </summary>
        /// <param name="matrix">The input matrix</param>
        /// <returns>The Q and R matrix</returns>
        public static Matrix<float>[] QRDecomposition(Matrix<float> matrix)
        {
            QR<float> qr = matrix.QR();
            Matrix<float> q = qr.Q;
            Matrix<float> r = qr.R;

            return new[] {q, r};
        }
    }
}