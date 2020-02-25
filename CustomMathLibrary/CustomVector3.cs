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
        
        private static readonly CustomVector3 zeroCusVec = new CustomVector3(0.0f,0.0f,0.0f);

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

        public void Normalize()
        {
            float num = CustomVector3.Magnitude(this);
            if ((double) num > 9.99999974737875E-06)
                this = this / num;
            else
                this = CustomVector3.zero;
        }

        public static CustomVector3 Normalize(CustomVector3 vector)
        {
            float number = CustomVector3.Magnitude(vector);
            return (double) number > 9.99999974737875E-06 ? vector / number : CustomVector3.zero;
        }

        public CustomVector3 normalized => CustomVector3.Normalize(this);

        /// <summary>
        /// Returns the length of the vector
        /// </summary>
        public float magnitude => (float) Math.Sqrt((double) this.x * (double) this.x + (double) this.y * (double) this.y + (double) this.z * (double) this.z);

        public static float Magnitude(CustomVector3 a)
        {
            return (float) Math.Sqrt((double) a.x * (double) a.x + (double) a.y * (double) a.y + (double) a.z * (double) a.z);
        }

        public static CustomVector3 zero => CustomVector3.zeroCusVec;


        public static CustomVector3 Cross(CustomVector3 lhs, CustomVector3 rhs)
        {
            return new CustomVector3((float) ((double) lhs.y * (double) rhs.z - (double) lhs.z * (double) rhs.y), (float) ((double) lhs.z * (double) rhs.x - (double) lhs.x * (double) rhs.z), (float) ((double) lhs.x * (double) rhs.y - (double) lhs.y * (double) rhs.x));
        }

        public static float Dot(CustomVector3 lhs, CustomVector3 rhs)
        {
            return (float) ((double) lhs.x * (double) rhs.x + (double) lhs.y * (double) rhs.y +(double) lhs.z * (double) rhs.z);
        }
    }
}