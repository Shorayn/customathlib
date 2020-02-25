using System;

namespace CustomMathLibrary
{
    public struct CustomVector3
    {
        

        public float x;
        public float y;
        public float z;

        public CustomVector3(float x, float y, float z)
        {
            this.x = x;
            this.y = y;
            this.z = z;
        }

        public static CustomVector3 operator +(CustomVector3 a, CustomVector3 b)
        {
            return new CustomVector3(a.x + b.x, a.y + b.y, a.z + b.z);
        }
        
        public static CustomVector3 operator -(CustomVector3 a, CustomVector3 b)
        {
            return new CustomVector3(a.x - b.x, a.y - b.y, a.z - b.z);
        }
        
        public static CustomVector3 operator *(CustomVector3 a, float d)
        {
            return new CustomVector3(a.x * d, a.y * d, a.z * d);
        }
        
        public static CustomVector3 operator /(CustomVector3 a, float d)
        {
            return new CustomVector3(a.x / d, a.y / d, a.z / d);
        }
        
        public static bool operator ==(CustomVector3 lhs, CustomVector3 rhs)
        {
            float num1 = lhs.x - rhs.x;
            float num2 = lhs.y - rhs.y;
            float num3 = lhs.z - rhs.z;
            return (double) num1 * (double) num1 + (double) num2 * (double) num2 + (double) num3 * (double) num3 < 9.99999943962493E-11;
        }

        public static bool operator !=(CustomVector3 lhs, CustomVector3 rhs)
        {
            return !(lhs == rhs);
        }

        public override string ToString()
        {
            return $"X-Value: {this.x}, Y-Value: {this.y}, Z-Value: {this.z}";
        }
        
        public bool Equals(CustomVector3 other)
        {
            return x.Equals(other.x) && y.Equals(other.y) && z.Equals(other.z);
        }

        public override bool Equals(object obj)
        {
            return obj is CustomVector3 other && Equals(other);
        }

        public override int GetHashCode()
        {
            unchecked
            {
                var hashCode = x.GetHashCode();
                hashCode = (hashCode * 397) ^ y.GetHashCode();
                hashCode = (hashCode * 397) ^ z.GetHashCode();
                return hashCode;
            }
        }

        // TODO: CONTINUE HERE
        public static CustomVector3 Normalize(CustomVector3 a)
        {
            return new CustomVector3();
        }

        /// <summary>
        /// Return the length of the vector
        /// </summary>
        public float magnitude => (float) Math.Sqrt((double) this.x * (double) this.x + (double) this.y * (double) this.y + (double) this.z * (double) this.z);
    }
}