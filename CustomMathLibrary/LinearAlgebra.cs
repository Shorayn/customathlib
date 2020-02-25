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
            CustomVector3 convertedPointOnOriginalPlane = CustomVector3Converter.FloatArrayToCustomVector3(pointOnOriginalPlane);

            CustomVector3 pointOnCopiedPlane = convertedPointOnOriginalPlane + convertedNormal * factor;
            float dNew = -(convertedNormal.x * pointOnCopiedPlane.x + convertedNormal.y * pointOnCopiedPlane.y +
                           convertedNormal.z * pointOnCopiedPlane.z);

            return new[] {dNew, convertedNormal.x, convertedNormal.y, convertedNormal.z};
        }
        
        
        private float[] DetermineDirectionVectorForSymmetryLine(float[] normal, float[] direction)
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

            Matrix<float> A = Matrix<float>.Build.DenseOfArray(new float[,]
            {
                {direction[0], direction[1], direction[2]},
                {normal[0], normal[1], normal[2]},
                {0, 0, 0}
            });
            
            QR<float> qr = A.QR();
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


            //Debug.Log(wPortion);

            Vector<float> rV = q.Transpose() * vPortion;
            Vector<float> rW = q.Transpose() * wPortion;


            //  r = (r1,r2,r3) wobei bspw. xi = vi + lambda * wi ist
            float r1 = rV[0] + 1 * rW[0];
            float r2 = rV[1] + 1 * rW[1];
            float r3 = rV[2] + 1 * rW[2];

            CustomVector3 result = new CustomVector3(r1, r2, r3);

            return null;
            //return result.normalized;
        }
        
        
       
    }
}