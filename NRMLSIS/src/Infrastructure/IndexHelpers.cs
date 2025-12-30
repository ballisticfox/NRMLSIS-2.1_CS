// #######################################################################
// MSIS® (NRL-SOF-014-1) SOFTWARE
// NRLMSIS® empirical atmospheric model software. Use is governed by the
// Open Source Academic Research License Agreement contained in the file
// nrlmsis2.1_license.txt, which is part of this software package. BY
// USING OR MODIFYING THIS SOFTWARE, YOU ARE AGREEING TO THE TERMS AND
// CONDITIONS OF THE LICENSE.
// #######################################################################

using System;

namespace NRLMSIS.Infrastructure
{
    /// <summary>
    /// Constants for B-spline calculations, making index mapping explicit and type-safe.
    /// Fortran MSIS uses negative indices for splines (e.g., -5:0), which we map to C# 0-based arrays.
    /// </summary>
    public static class SplineConstants
    {
        /// <summary>Minimum Fortran spline index (-5)</summary>
        public const int FortranMinIndex = -5;
        
        /// <summary>Maximum Fortran spline index (0)</summary>
        public const int FortranMaxIndex = 0;
        
        /// <summary>Offset to convert Fortran spline index to C# array index (+5)</summary>
        public const int IndexOffset = 5;
        
        /// <summary>Minimum spline order (2 = linear)</summary>
        public const int MinOrder = 2;
        
        /// <summary>Maximum spline order (6 = quintic)</summary>
        public const int MaxOrder = 6;
        
        /// <summary>Offset to convert spline order to C# array index (-2)</summary>
        public const int OrderOffset = 2;
        
        /// <summary>Number of spline indices (6, from -5 to 0)</summary>
        public const int IndexCount = 6;
        
        /// <summary>Number of spline orders (5, from order 2 to 6)</summary>
        public const int OrderCount = 5;
    }
    
    /// <summary>
    /// Enumeration of B-spline orders for type safety.
    /// </summary>
    public enum SplineOrder
    {
        /// <summary>Linear spline (order 2)</summary>
        Linear = 2,
        
        /// <summary>Quadratic spline (order 3)</summary>
        Quadratic = 3,
        
        /// <summary>Cubic spline (order 4)</summary>
        Cubic = 4,
        
        /// <summary>Quartic spline (order 5)</summary>
        Quartic = 5,
        
        /// <summary>Quintic spline (order 6)</summary>
        Quintic = 6
    }
    
    /// <summary>
    /// Type-safe wrapper for B-spline weight arrays with Fortran-compatible indexing.
    /// Handles the conversion between Fortran's Sz(-5:0, 2:6) and C# Sz[0:5, 0:4].
    /// </summary>
    /// <remarks>
    /// In Fortran, the spline array is indexed as Sz(i, k) where:
    /// - i ranges from -5 to 0 (spline index relative to reference point)
    /// - k ranges from 2 to 6 (spline order)
    /// 
    /// In C#, we use a [6,5] array where:
    /// - First dimension [0:5] maps to Fortran i=-5:0 via index = i + 5
    /// - Second dimension [0:4] maps to Fortran k=2:6 via index = k - 2
    /// </remarks>
    public readonly struct SplineWeights
    {
        private readonly double[,] _data;
        
        /// <summary>
        /// Creates a new SplineWeights wrapper around an existing array.
        /// </summary>
        /// <param name="data">The underlying [6,5] spline weight array</param>
        /// <exception cref="ArgumentException">
        /// Thrown if the array is not [6,5] in size
        /// </exception>
        public SplineWeights(double[,] data)
        {
            if (data == null)
                throw new ArgumentNullException(nameof(data));
            
            if (data.GetLength(0) != SplineConstants.IndexCount || 
                data.GetLength(1) != SplineConstants.OrderCount)
            {
                throw new ArgumentException(
                    $"Spline weight array must be [{SplineConstants.IndexCount},{SplineConstants.OrderCount}] " +
                    $"but received [{data.GetLength(0)},{data.GetLength(1)}]",
                    nameof(data));
            }
            
            _data = data;
        }
        
        /// <summary>
        /// Access spline weight using Fortran-style indices.
        /// </summary>
        /// <param name="fortranIndex">Spline index in Fortran convention (-5 to 0)</param>
        /// <param name="order">Spline order (2 to 6)</param>
        /// <returns>The spline weight value</returns>
        /// <exception cref="ArgumentOutOfRangeException">
        /// Thrown if indices are out of valid range
        /// </exception>
        public double this[int fortranIndex, SplineOrder order]
        {
            get
            {
                ValidateIndices(fortranIndex, order);
                int i = fortranIndex + SplineConstants.IndexOffset;
                int k = (int)order - SplineConstants.OrderOffset;
                return _data[i, k];
            }
            set
            {
                ValidateIndices(fortranIndex, order);
                int i = fortranIndex + SplineConstants.IndexOffset;
                int k = (int)order - SplineConstants.OrderOffset;
                _data[i, k] = value;
            }
        }
        
        /// <summary>
        /// Access spline weight using C#-style 0-based indices (direct array access).
        /// </summary>
        /// <param name="csharpIndex">C# array index (0 to 5)</param>
        /// <param name="csharpOrder">C# array index for order (0 to 4)</param>
        public double GetDirect(int csharpIndex, int csharpOrder)
        {
            return _data[csharpIndex, csharpOrder];
        }
        
        /// <summary>
        /// Set spline weight using C#-style 0-based indices (direct array access).
        /// </summary>
        public void SetDirect(int csharpIndex, int csharpOrder, double value)
        {
            _data[csharpIndex, csharpOrder] = value;
        }
        
        /// <summary>
        /// Gets the underlying data array (for performance-critical code).
        /// </summary>
        public double[,] RawData => _data;
        
        /// <summary>
        /// Extracts spline weights for a specific order into a pre-allocated buffer.
        /// This is used to get the 4 weights needed for temperature calculation.
        /// </summary>
        /// <param name="order">The spline order to extract</param>
        /// <param name="startIndex">Starting Fortran index (-3 for typical use)</param>
        /// <param name="count">Number of weights to extract (typically 4)</param>
        /// <param name="buffer">Pre-allocated buffer to receive the weights</param>
        public void ExtractWeights(SplineOrder order, int startIndex, int count, Span<double> buffer)
        {
            if (buffer.Length < count)
                throw new ArgumentException($"Buffer must have at least {count} elements", nameof(buffer));
            
            int orderIndex = (int)order - SplineConstants.OrderOffset;
            
            for (int i = 0; i < count; i++)
            {
                int fortranIdx = startIndex + i;
                int csharpIdx = fortranIdx + SplineConstants.IndexOffset;
                buffer[i] = _data[csharpIdx, orderIndex];
            }
        }
        
        private static void ValidateIndices(int fortranIndex, SplineOrder order)
        {
            if (fortranIndex < SplineConstants.FortranMinIndex || 
                fortranIndex > SplineConstants.FortranMaxIndex)
            {
                throw new ArgumentOutOfRangeException(
                    nameof(fortranIndex),
                    fortranIndex,
                    $"Fortran spline index must be in range [{SplineConstants.FortranMinIndex}, {SplineConstants.FortranMaxIndex}]");
            }
            
            if ((int)order < SplineConstants.MinOrder || 
                (int)order > SplineConstants.MaxOrder)
            {
                throw new ArgumentOutOfRangeException(
                    nameof(order),
                    order,
                    $"Spline order must be in range [{SplineConstants.MinOrder}, {SplineConstants.MaxOrder}]");
            }
        }
    }
    
    /// <summary>
    /// Helper methods for working with MSIS arrays and indices.
    /// </summary>
    public static class IndexHelper
    {
        /// <summary>
        /// Converts a Fortran spline index to C# array index.
        /// </summary>
        /// <param name="fortranIndex">Fortran index (-5 to 0)</param>
        /// <returns>C# array index (0 to 5)</returns>
        public static int FortranToCSharpSplineIndex(int fortranIndex)
        {
            if (fortranIndex < SplineConstants.FortranMinIndex || 
                fortranIndex > SplineConstants.FortranMaxIndex)
            {
                throw new ArgumentOutOfRangeException(
                    nameof(fortranIndex),
                    fortranIndex,
                    $"Fortran index must be in range [{SplineConstants.FortranMinIndex}, {SplineConstants.FortranMaxIndex}]");
            }
            
            return fortranIndex + SplineConstants.IndexOffset;
        }
        
        /// <summary>
        /// Converts a C# array index to Fortran spline index.
        /// </summary>
        /// <param name="csharpIndex">C# array index (0 to 5)</param>
        /// <returns>Fortran index (-5 to 0)</returns>
        public static int CSharpToFortranSplineIndex(int csharpIndex)
        {
            if (csharpIndex < 0 || csharpIndex >= SplineConstants.IndexCount)
            {
                throw new ArgumentOutOfRangeException(
                    nameof(csharpIndex),
                    csharpIndex,
                    $"C# index must be in range [0, {SplineConstants.IndexCount - 1}]");
            }
            
            return csharpIndex - SplineConstants.IndexOffset;
        }
        
        /// <summary>
        /// Converts a spline order to C# array index.
        /// </summary>
        /// <param name="order">Spline order (2 to 6)</param>
        /// <returns>C# array index (0 to 4)</returns>
        public static int SplineOrderToIndex(SplineOrder order)
        {
            return (int)order - SplineConstants.OrderOffset;
        }
        
        /// <summary>
        /// Converts a C# array index to spline order.
        /// </summary>
        /// <param name="index">C# array index (0 to 4)</param>
        /// <returns>Spline order (2 to 6)</returns>
        public static SplineOrder IndexToSplineOrder(int index)
        {
            if (index < 0 || index >= SplineConstants.OrderCount)
            {
                throw new ArgumentOutOfRangeException(
                    nameof(index),
                    index,
                    $"Index must be in range [0, {SplineConstants.OrderCount - 1}]");
            }
            
            return (SplineOrder)(index + SplineConstants.OrderOffset);
        }
    }
    
    /// <summary>
    /// Constants for altitude regimes in the atmosphere.
    /// Makes the code more readable than using magic numbers.
    /// </summary>
    public static class AltitudeRegimes
    {
        /// <summary>
        /// Fully mixed region boundary (70 km).
        /// Below this altitude, the atmosphere is well-mixed with constant mixing ratios.
        /// </summary>
        public const double ZetaF = 70.0;
        
        /// <summary>
        /// Bates profile boundary (122.5 km).
        /// Above this altitude, temperature follows the Bates profile.
        /// Below this altitude, temperature is determined by splines.
        /// </summary>
        public const double ZetaB = 122.5;
        
        /// <summary>
        /// Reference altitude for active minor species (85 km).
        /// </summary>
        public const double ZetaA = 85.0;
        
        /// <summary>
        /// Reference height for chemical/dynamical taper (100 km).
        /// </summary>
        public const double ZetaGamma = 100.0;
        
        /// <summary>
        /// Inverse scale height of chemical/dynamical taper (1/30 km^-1).
        /// </summary>
        public const double HGamma = 1.0 / 30.0;
        
        /// <summary>
        /// Checks if altitude is in the well-mixed region (below 70 km).
        /// </summary>
        public static bool IsFullyMixed(double altitude) => altitude < ZetaF;
        
        /// <summary>
        /// Checks if altitude is in the spline region (70-122.5 km).
        /// </summary>
        public static bool IsSplineRegion(double altitude) => 
            altitude >= ZetaF && altitude < ZetaB;
        
        /// <summary>
        /// Checks if altitude is in the Bates profile region (>122.5 km).
        /// </summary>
        public static bool IsBatesRegion(double altitude) => altitude >= ZetaB;
        
        /// <summary>
        /// Gets a human-readable description of the altitude regime.
        /// </summary>
        public static string GetRegionName(double altitude)
        {
            if (IsFullyMixed(altitude))
                return "Fully Mixed (Homosphere)";
            else if (IsSplineRegion(altitude))
                return "Transition Region (Spline)";
            else
                return "Upper Thermosphere (Bates Profile)";
        }
    }
}