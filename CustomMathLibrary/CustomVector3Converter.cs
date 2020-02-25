namespace CustomMathLibrary
{
    public static class CustomVector3Converter
    {
        /// <summary>
		/// Converts a float array to a CustomVector3
		/// </summary>
		/// <param name="vectorData"></param>
		/// <returns>Vector3 with 3 entries</returns>
		public static CustomVector3 FloatArrayToCustomVector3(float[] vectorData)
		{
			return new CustomVector3(vectorData[0], vectorData[1], vectorData[2]);
		}
   
   
		/// <summary>
		/// Converts CustomVector3 to float array
		/// </summary>
		/// <param name="customVector3"></param>
		/// <returns>Float array with x, y, z value</returns>
		public static float[] CustomVector3ToFloatArray(CustomVector3 customVector3)
		{
			return new float[]{customVector3.x, customVector3.y, customVector3.z};
		}
    
		public static CustomVector3 DoubleArrayToCustomVector3(double[] values)
		{
			return new CustomVector3((float) values[0], (float) values[1], (float) values[2]);
		}
   
		public static double[] CustomVector3ToDoubleArray(CustomVector3 customVector3)
		{
			return new double[] {customVector3.x, customVector3.y, customVector3.z};
		}

		
		/// <summary>
		/// Converts a 2 dimensional float array to a CustomVector3 array
		/// </summary>
		/// <param name="data"></param>
		/// <returns>Vector3 array</returns>
		public static CustomVector3[] CustomVector3ArrayParser(float[,] data)
		{
			CustomVector3[] tmpVertices = new CustomVector3[data.GetLength(0)];

			for (int i = 0; i < data.GetLength(0); i++)
			{
				CustomVector3 current = new CustomVector3(data[i,0], data[i,1], data[i,2]);
				tmpVertices[i] = current;

			}

			return tmpVertices;
		}

		/// <summary>
		/// Converts a CustomVector3 array to a 2 dimensional float array
		/// </summary>
		/// <param name="tempVertices"></param>
		/// <returns>2 dimensional float array</returns>
		public static float[,] TwoDimFloatArrayVectorParser(CustomVector3[] tempVertices)
		{
			float[,] tmp = new float[tempVertices.Length,3];

			for (int i = 0; i < tempVertices.Length; i++)
			{
				float[] current = CustomVector3ToFloatArray(tempVertices[i]);
				tmp[i,0] = current[0];
				tmp[i,1] = current[1];
				tmp[i,2] = current[2];
            
			}

			return tmp;
		}

		
		public static double[][] TwoDimDoubleArrayVectorParser(CustomVector3[] tempVertices)
		{
			double[][] tmp = new double[tempVertices.Length][];

			for (int i = 0; i < tempVertices.Length; i++)
			{
				double[] current = CustomVector3ToDoubleArray(tempVertices[i]);
				tmp[i] = current;
			}

			return tmp;
		}
        
        
        
        
    }
}