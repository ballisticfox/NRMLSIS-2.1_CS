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
    // ============================================================================
    //  SPLINE INDEXING
    // ============================================================================

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

        /// <summary>Minimum Fortran weight index (-4)</summary>
        public const int FortranMinWeightIndex = -4;

        /// <summary>Offset for weight array (w[-4:0] maps to w[0:4])</summary>
        public const int WeightOffset = 4;
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

    // ============================================================================
    //  DENSITY PROFILE COLUMN INDICES
    // ============================================================================

    /// <summary>
    /// Named constants for Beta array columns in density profile calculations.
    /// Makes parameter extraction self-documenting instead of using magic numbers.
    /// </summary>
    public static class DensityProfileColumns
    {
        // Column indices (before Bl offset is applied)

        /// <summary>Column 0: Log of reference density or mixing ratio</summary>
        public const int LogReferenceDensity = 0;

        /// <summary>Column 1: Turbopause height (ZetaM)</summary>
        public const int TurbopauseHeight = 1;

        /// <summary>Column 2: Scale height of lower portion of effective mass profile (HML)</summary>
        public const int LowerScaleHeight = 2;

        /// <summary>Column 3: Scale height of upper portion of effective mass profile (HMU)</summary>
        public const int UpperScaleHeight = 3;

        /// <summary>Column 4: Chapman term coefficient (C)</summary>
        public const int ChapmanCoefficient = 4;

        /// <summary>Column 5: Chapman term reference height (ZetaC)</summary>
        public const int ChapmanReferenceHeight = 5;

        /// <summary>Column 6: Chapman term scale height (HC)</summary>
        public const int ChapmanScaleHeight = 6;

        /// <summary>Column 7: Chemical/dynamical term coefficient (R)</summary>
        public const int ChemicalDynamicCoefficient = 7;

        /// <summary>Column 8: Chemical/dynamical term reference height (ZetaR)</summary>
        public const int ChemicalDynamicReferenceHeight = 8;

        /// <summary>Column 9: Chemical/dynamical term scale height (HR)</summary>
        public const int ChemicalDynamicScaleHeight = 9;
    }

    // ============================================================================
    //  BASIS FUNCTION ARRAY SIZES
    // ============================================================================

    /// <summary>
    /// Constants for basis function array extraction.
    /// Documents the size of various basis function subsets.
    /// </summary>
    public static class BasisFunctionSizes
    {
        /// <summary>Number of primary geomagnetic basis functions (13)</summary>
        public const int GeomagneticPrimary = 13;

        /// <summary>Number of secondary geomagnetic basis functions per component (7)</summary>
        public const int GeomagneticSecondary = 7;

        /// <summary>Number of UT-dependent basis functions (9)</summary>
        public const int UniversalTimeDependent = 9;

        /// <summary>Total geomagnetic basis functions (13 + 2*7 = 27)</summary>
        public const int GeomagneticTotal = GeomagneticPrimary + 2 * GeomagneticSecondary;
    }

    // ============================================================================
    //  SPLINE WEIGHT ARRAY HELPERS
    // ============================================================================

    /// <summary>
    /// Helper methods for accessing B-spline weight arrays with Fortran-style indexing.
    /// Provides cleaner syntax than manual index arithmetic.
    /// </summary>
    public static class SplineWeightHelpers
    {
        /// <summary>
        /// Gets a weight value using Fortran-style negative indexing.
        /// w[-4:0] maps to w[0:4] in C#
        /// </summary>
        /// <param name="weights">Weight array (size 5)</param>
        /// <param name="fortranIndex">Fortran index (-4 to 0)</param>
        public static double GetWeight(double[] weights, int fortranIndex)
        {
            if (fortranIndex < -4 || fortranIndex > 0)
                throw new ArgumentOutOfRangeException(nameof(fortranIndex),
                    $"Fortran weight index must be in range [-4, 0], got {fortranIndex}");

            return weights[fortranIndex + SplineConstants.WeightOffset];
        }

        /// <summary>
        /// Sets a weight value using Fortran-style negative indexing.
        /// w[-4:0] maps to w[0:4] in C#
        /// </summary>
        /// <param name="weights">Weight array (size 5)</param>
        /// <param name="fortranIndex">Fortran index (-4 to 0)</param>
        /// <param name="value">Value to set</param>
        public static void SetWeight(double[] weights, int fortranIndex, double value)
        {
            if (fortranIndex < -4 || fortranIndex > 0)
                throw new ArgumentOutOfRangeException(nameof(fortranIndex),
                    $"Fortran weight index must be in range [-4, 0], got {fortranIndex}");

            weights[fortranIndex + SplineConstants.WeightOffset] = value;
        }

        /// <summary>
        /// Gets a spline value using Fortran-style indexing.
        /// S[-5:0, 2:6] maps to S[0:5, 0:4] in C#
        /// </summary>
        /// <param name="splines">Spline array [6,5]</param>
        /// <param name="fortranSplineIndex">Fortran spline index (-5 to 0)</param>
        /// <param name="order">Spline order (2 to 6)</param>
        public static double GetSpline(double[,] splines, int fortranSplineIndex, int order)
        {
            if (fortranSplineIndex < -5 || fortranSplineIndex > 0)
                throw new ArgumentOutOfRangeException(nameof(fortranSplineIndex),
                    $"Fortran spline index must be in range [-5, 0], got {fortranSplineIndex}");

            if (order < 2 || order > 6)
                throw new ArgumentOutOfRangeException(nameof(order),
                    $"Spline order must be in range [2, 6], got {order}");

            int i = fortranSplineIndex + SplineConstants.IndexOffset;
            int k = order - SplineConstants.OrderOffset;
            return splines[i, k];
        }

        /// <summary>
        /// Sets a spline value using Fortran-style indexing.
        /// S[-5:0, 2:6] maps to S[0:5, 0:4] in C#
        /// </summary>
        /// <param name="splines">Spline array [6,5]</param>
        /// <param name="fortranSplineIndex">Fortran spline index (-5 to 0)</param>
        /// <param name="order">Spline order (2 to 6)</param>
        /// <param name="value">Value to set</param>
        public static void SetSpline(double[,] splines, int fortranSplineIndex, int order, double value)
        {
            if (fortranSplineIndex < -5 || fortranSplineIndex > 0)
                throw new ArgumentOutOfRangeException(nameof(fortranSplineIndex),
                    $"Fortran spline index must be in range [-5, 0], got {fortranSplineIndex}");

            if (order < 2 || order > 6)
                throw new ArgumentOutOfRangeException(nameof(order),
                    $"Spline order must be in range [2, 6], got {order}");

            int i = fortranSplineIndex + SplineConstants.IndexOffset;
            int k = order - SplineConstants.OrderOffset;
            splines[i, k] = value;
        }
    }

    // ============================================================================
    //  BASIS FUNCTION EXTRACTION HELPERS
    // ============================================================================

    /// <summary>
    /// Helper methods for extracting basis function subsets from the global basis function array.
    /// Makes array slicing operations self-documenting.
    /// </summary>
    public static class BasisFunctionExtractor
    {
        /// <summary>
        /// Extracts primary geomagnetic basis functions (13 elements starting at CMag).
        /// </summary>
        /// <param name="globalBasisFunctions">Global basis function array</param>
        /// <returns>Array of 13 primary geomagnetic basis functions</returns>
        public static double[] ExtractGeomagneticPrimary(double[] globalBasisFunctions)
        {
            double[] result = new double[BasisFunctionSizes.GeomagneticPrimary];
            Array.Copy(globalBasisFunctions, Constants.CMag, result, 0, BasisFunctionSizes.GeomagneticPrimary);
            return result;
        }

        /// <summary>
        /// Extracts secondary geomagnetic basis functions (7x2 elements after primary).
        /// </summary>
        /// <param name="globalBasisFunctions">Global basis function array</param>
        /// <returns>Array [7,2] of secondary geomagnetic basis functions</returns>
        public static double[,] ExtractGeomagneticSecondary(double[] globalBasisFunctions)
        {
            double[,] result = new double[BasisFunctionSizes.GeomagneticSecondary, 2];
            int offset = Constants.CMag + BasisFunctionSizes.GeomagneticPrimary;

            for (int i = 0; i < BasisFunctionSizes.GeomagneticSecondary; i++)
            {
                result[i, 0] = globalBasisFunctions[offset + i];
                result[i, 1] = globalBasisFunctions[offset + BasisFunctionSizes.GeomagneticSecondary + i];
            }

            return result;
        }

        /// <summary>
        /// Extracts UT-dependent basis functions (9 elements starting at CUt).
        /// </summary>
        /// <param name="globalBasisFunctions">Global basis function array</param>
        /// <returns>Array of 9 UT-dependent basis functions</returns>
        public static double[] ExtractUniversalTimeDependent(double[] globalBasisFunctions)
        {
            double[] result = new double[BasisFunctionSizes.UniversalTimeDependent];
            Array.Copy(globalBasisFunctions, Constants.CUt, result, 0, BasisFunctionSizes.UniversalTimeDependent);
            return result;
        }

        /// <summary>
        /// Extracts geomagnetic parameters from a BasisSubset Beta array.
        /// </summary>
        /// <param name="subset">The basis subset containing Beta parameters</param>
        /// <param name="column">Column index (before Bl offset)</param>
        /// <returns>Array of geomagnetic parameters</returns>
        public static double[] ExtractGeomagneticParameters(BasisSubset subset, int column)
        {
            double[] parms = new double[Constants.NMag];
            int columnIndex = column - subset.Bl;

            for (int i = 0; i < Constants.NMag; i++)
            {
                parms[i] = subset.Beta[Constants.CMag + i, columnIndex];
            }

            return parms;
        }
    }

    // ============================================================================
    //  BOUNDS CHECKING HELPER
    // ============================================================================

    /// <summary>
    /// Helper methods for common bounds checking patterns.
    /// </summary>
    public static class BoundsHelper
    {
        /// <summary>
        /// Checks if an index is within valid bounds [min, max).
        /// </summary>
        /// <param name="index">Index to check</param>
        /// <param name="min">Minimum valid value (inclusive)</param>
        /// <param name="max">Maximum valid value (exclusive)</param>
        /// <returns>True if index is in [min, max)</returns>
        public static bool IsInBounds(int index, int min, int max)
        {
            return index >= min && index < max;
        }

        /// <summary>
        /// Checks if an index plus offset is within valid bounds [min, max).
        /// Common pattern in spline calculations.
        /// </summary>
        /// <param name="baseIndex">Base index</param>
        /// <param name="offset">Offset to add</param>
        /// <param name="min">Minimum valid value (inclusive)</param>
        /// <param name="max">Maximum valid value (exclusive)</param>
        /// <returns>True if (baseIndex + offset) is in [min, max)</returns>
        public static bool IsInBounds(int baseIndex, int offset, int min, int max)
        {
            int index = baseIndex + offset;
            return index >= min && index < max;
        }
    }

    /// <summary>
    /// Type-safe wrapper for B-spline weight arrays with Fortran-compatible indexing.
    /// Handles the conversion between Fortran's Sz(-5:0, 2:6) and C# Sz[0:5, 0:4].
    /// </summary>
    public readonly struct SplineWeights
    {
        private readonly double[,] _data;

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

        public double GetDirect(int csharpIndex, int csharpOrder)
        {
            return _data[csharpIndex, csharpOrder];
        }

        public void SetDirect(int csharpIndex, int csharpOrder, double value)
        {
            _data[csharpIndex, csharpOrder] = value;
        }

        public double[,] RawData => _data;

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

        public static int SplineOrderToIndex(SplineOrder order)
        {
            return (int)order - SplineConstants.OrderOffset;
        }

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
        public const double ZetaF = 70.0;
        public const double ZetaB = 122.5;
        public const double ZetaA = 85.0;
        public const double ZetaGamma = 100.0;
        public const double HGamma = 1.0 / 30.0;

        public static bool IsFullyMixed(double altitude) => altitude < ZetaF;
        public static bool IsSplineRegion(double altitude) => altitude >= ZetaF && altitude < ZetaB;
        public static bool IsBatesRegion(double altitude) => altitude >= ZetaB;

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