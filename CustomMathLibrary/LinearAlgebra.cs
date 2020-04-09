using System;
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

        /// <summary>
        /// Checks whether two 3 dimensional vectors are scalar to each other with a small epsilon
        /// </summary>
        /// <param name="first">First vector</param>
        /// <param name="second">Second vector</param>
        /// <returns>bool value</returns>
        public static bool CheckDotProduct(double[] first, double[] second)
        {
            double result = first[0] * second[0] + first[1] * second[1] +
                            first[2] * second[2];
            double epsilon = 1e-16;
            return Math.Abs(result) < epsilon;
        }

        /// <summary>
        /// Calculate whether a straight line defined by its local vector and direction vector has a certain intersection with the given plane (in coordinate form)
        /// using the lambda of the straight line
        /// 
        ///   http://www.songho.ca/math/line/line.html
        ///   https://de.serlo.org/mathe/geometrie/analytische-geometrie/lagebeziehung-punkten-geraden-ebenen/lagebeziehung-geraden-einer-ebene/lagebeziehungen-geraden-ebenen
        /// 
        /// </summary>
        /// <param name="planeCoords">The coefficients a, b, c, and d of the plane</param>
        /// <param name="localVec">The local vector of the straight line</param>
        /// <param name="dirVecPlus">+ The direction vector of the straight line</param>
        /// <param name="dirVecMinus">- The direction vector of the straight line</param>
        /// <returns>The intersection point as a CustomVector3</returns>
        private CustomVector3 DetermineIntersectionForLineAndPlane(float[] planeCoords, CustomVector3 localVec,
            CustomVector3 dirVecPlus,
            CustomVector3 dirVecMinus)
        {
            float nominator = -(planeCoords[0] * localVec.x + planeCoords[1] * localVec.y +
                                planeCoords[2] * localVec.z +
                                planeCoords[3]);
            float denominatorPlus = planeCoords[0] * dirVecPlus.x + planeCoords[1] * dirVecPlus.y +
                                    planeCoords[2] * dirVecPlus.z;
            float denominatorMinus = planeCoords[0] * dirVecMinus.x + planeCoords[1] * dirVecMinus.y +
                                     planeCoords[2] * dirVecMinus.z;


            

            if (!float.IsNaN(nominator / denominatorPlus))
            {
                float lambda =
                    nominator / denominatorPlus;
                return localVec + lambda * dirVecPlus;
            }

            if (!float.IsNaN(nominator / denominatorMinus))
            {
                float lambda =
                    nominator / denominatorMinus;
                return localVec + lambda * dirVecPlus;
            }

            return CustomVector3.positiveInfinity;
        }
        
        /// <summary>
        ///  Converts the parameter form of a plane equation to its coordinate equivalent
        /// </summary>
        /// <param name="localVec">The local vector of the plane</param>
        /// <param name="dirVecOne">The first direction vector of the plane</param>
        /// <param name="dirVecTwo">The second direction vector of the plane</param>
        /// <returns>The parameters a, b, c, and d of this plane</returns>
        private float[] ConvertParameterTooCoordinateForm(CustomVector3 localVec, CustomVector3 dirVecOne, CustomVector3 dirVecTwo)
        {
            CustomVector3 normal = CustomVector3.Cross(dirVecOne, dirVecTwo);
            return new[]
            {
                normal.x, normal.y, normal.z, -(normal.x * localVec.x + normal.y * localVec.y + normal.z * localVec.z)
            };
        }
        
        /// <summary>
        /// Converts the coordinate form of a given plane equation to its parametric representation.
        /// Coordinate form
        ///                 P: a*x + b*y + c*z + d = 0
        /// Parametric form
        ///                 P: r0 + lambda * v + mu * u where r0, v and u are vectors
        ///
        ///  https://de.wikipedia.org/wiki/Parameterform
        /// 
        /// </summary>
        /// <param name="parameters">The parameters d, a, b, c from plane in this order</param>
        /// <returns>The location r0 and the two distance vectors v and u</returns>
        public float[][] ConvertCoordinateFormToParameterForm(float[] parameters)
        {
            float a = parameters[1];
            float b = parameters[2];
            float c = parameters[3];
            float d = parameters[0];

            double epsilon = 1e-16;
            float[] dirVecU = {-b, a, 0};
            if (Math.Abs(dirVecU[0] + dirVecU[1] + dirVecU[2]) < epsilon)
            {
                
            }
            
            dirVecU[0] = 0;
            dirVecU[1] = -c;
            dirVecU[2] = b;
            if (Math.Abs(dirVecU[0] + dirVecU[1] + dirVecU[2]) < epsilon)
            {
                dirVecU[0] = -c;
                dirVecU[1] = 0;
                dirVecU[2] = a;
            }


            float[] vPortion = {0, -c, b};
            float[] pPortion = {0f, 0, 0};

            if (!float.IsNaN(d / a))
            {
                pPortion = new[] {d / a, 0, 0};
            }
            else if (!float.IsNaN(d / b))
            {
                pPortion = new[] {0, d / b, 0};
            }
            else if (!float.IsNaN(d / c))
            {
                pPortion = new[] {0, 0, d / c};
            }

            return new[] {pPortion, dirVecU, vPortion};
        }
        
    }
}