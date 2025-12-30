// #######################################################################
// MSIS® (NRL-SOF-014-1) SOFTWARE
// NRLMSIS® empirical atmospheric model software. Use is governed by the
// Open Source Academic Research License Agreement contained in the file
// nrlmsis2.1_license.txt, which is part of this software package. BY
// USING OR MODIFYING THIS SOFTWARE, YOU ARE AGREEING TO THE TERMS AND
// CONDITIONS OF THE LICENSE.
// #######################################################################

using System;

namespace NRLMSIS
{
    /// <summary>
    /// Computation context that caches intermediate results for performance optimization.
    /// This class is NOT thread-safe by design - create one instance per thread.
    /// Reusing a context instance for multiple calculations improves performance by
    /// avoiding redundant computations when input parameters haven't changed.
    /// </summary>
    /// <remarks>
    /// Performance benefits:
    /// - Avoids recalculating basis functions when location/time unchanged
    /// - Avoids recalculating splines when altitude unchanged
    /// - Reduces memory allocations by reusing buffers
    /// - Improves cache locality by grouping related data
    ///
    /// Thread safety:
    /// - Each thread should have its own Context instance
    /// - Do NOT share a single context across multiple threads
    /// - For parallel calculations, create one context per thread
    /// </remarks>
    public sealed class Context : IDisposable
    {
        // ===== Cached Input Parameters =====

        /// <summary>Last calculated day of year (NaN if uncached)</summary>
        private double _lastDay = double.NaN;

        /// <summary>Last calculated universal time in seconds (NaN if uncached)</summary>
        private double _lastUtsec = double.NaN;

        /// <summary>Last calculated latitude in degrees (NaN if uncached)</summary>
        private double _lastLat = double.NaN;

        /// <summary>Last calculated longitude in degrees (NaN if uncached)</summary>
        private double _lastLon = double.NaN;

        /// <summary>Last calculated altitude in km (NaN if uncached)</summary>
        private double _lastAltitude = double.NaN;

        /// <summary>Last calculated F10.7 flux (NaN if uncached)</summary>
        private double _lastSflux = double.NaN;

        /// <summary>Last calculated F10.7 average flux (NaN if uncached)</summary>
        private double _lastSfluxavg = double.NaN;

        /// <summary>Last calculated Ap geomagnetic indices</summary>
        private readonly double[] _lastAp = new double[7];

        // ===== Cached Computed Values =====

        /// <summary>Cached horizontal and temporal basis functions (gf array)</summary>
        private readonly double[] _basisFunctions = new double[Constants.MaxNbf];

        /// <summary>Cached B-spline weights for current altitude</summary>
        private readonly double[,] _splineWeights = new double[6, 5];

        /// <summary>Cached B-spline reference index</summary>
        private int _splineIndex;

        /// <summary>Cached temperature profile parameters</summary>
        private TemperatureProfile _temperatureProfile = new TemperatureProfile();

        /// <summary>Cached density profile parameters for each species</summary>
        private readonly DensityParameters[] _densityProfiles = new DensityParameters[Constants.NSpec];

        // ===== Reusable Buffers (to avoid allocations) =====

        /// <summary>Reusable buffer for density output (10 elements)</summary>
        private readonly double[] _densityBuffer = new double[10];

        /// <summary>Reusable buffer for spline weights extraction (4 elements)</summary>
        private readonly double[] _weightBuffer = new double[4];

        /// <summary>Reusable buffer for spline work array (5 elements)</summary>
        private readonly double[] _splineWorkBuffer = new double[5];

        // ===== Initialization =====

        /// <summary>
        /// Creates a new computation context with cleared cache.
        /// </summary>
        public Context()
        {
            // Initialize density profiles array
            for (int i = 0; i < _densityProfiles.Length; i++)
            {
                _densityProfiles[i] = new DensityParameters();
            }

            // Initialize lastAp with sentinel values
            Array.Fill(_lastAp, double.NaN);
        }

        // ===== Cache Validation Methods =====

        /// <summary>
        /// Checks if the location and time parameters have changed since last calculation.
        /// </summary>
        /// <returns>True if parameters changed and cache is invalid, false if cache is still valid</returns>
        public bool IsLocationTimeChanged(double day, double utsec, double lat, double lon,
                                         double sflux, double sfluxavg, double[] ap)
        {
            // Check if any parameter has changed
            bool changed = day != _lastDay ||
                          utsec != _lastUtsec ||
                          lat != _lastLat ||
                          lon != _lastLon ||
                          sflux != _lastSflux ||
                          sfluxavg != _lastSfluxavg;

            // Check Ap array
            if (!changed)
            {
                for (int i = 0; i < 7; i++)
                {
                    if (ap[i] != _lastAp[i])
                    {
                        changed = true;
                        break;
                    }
                }
            }

            return changed;
        }

        /// <summary>
        /// Checks if the altitude parameter has changed since last calculation.
        /// </summary>
        public bool IsAltitudeChanged(double altitude)
        {
            return altitude != _lastAltitude;
        }

        // ===== Cache Update Methods =====

        /// <summary>
        /// Updates the cached location and time parameters.
        /// Call this after recalculating basis functions.
        /// </summary>
        public void UpdateLocationTimeCache(double day, double utsec, double lat, double lon,
                                           double sflux, double sfluxavg, double[] ap)
        {
            _lastDay = day;
            _lastUtsec = utsec;
            _lastLat = lat;
            _lastLon = lon;
            _lastSflux = sflux;
            _lastSfluxavg = sfluxavg;
            Array.Copy(ap, _lastAp, 7);
        }

        /// <summary>
        /// Updates the cached altitude parameter.
        /// Call this after recalculating spline weights.
        /// </summary>
        public void UpdateAltitudeCache(double altitude)
        {
            _lastAltitude = altitude;
        }

        // ===== Accessors for Cached Values =====

        /// <summary>Gets the cached basis functions array (read-write access)</summary>
        public double[] BasisFunctions => _basisFunctions;

        /// <summary>Gets the cached spline weights array (read-write access)</summary>
        public double[,] SplineWeights => _splineWeights;

        /// <summary>Gets or sets the cached spline index</summary>
        public int SplineIndex
        {
            get => _splineIndex;
            set => _splineIndex = value;
        }

        /// <summary>Gets the cached temperature profile parameters (read-write access)</summary>
        public TemperatureProfile TemperatureProfile => _temperatureProfile;

        /// <summary>Gets the cached density profile parameters array (read-write access)</summary>
        public DensityParameters[] DensityProfiles => _densityProfiles;

        // ===== Accessors for Reusable Buffers =====

        /// <summary>
        /// Gets the reusable density buffer as a span for zero-copy operations.
        /// The buffer contains 10 elements for density output.
        /// </summary>
        public Span<double> DensityBuffer => _densityBuffer;

        /// <summary>
        /// Gets the reusable weight buffer as a span.
        /// The buffer contains 4 elements for spline weight extraction.
        /// </summary>
        public Span<double> WeightBuffer => _weightBuffer;

        /// <summary>
        /// Gets the reusable spline work buffer as a span.
        /// The buffer contains 5 elements for spline calculations.
        /// </summary>
        public Span<double> SplineWorkBuffer => _splineWorkBuffer;

        // ===== Cache Management =====

        /// <summary>
        /// Resets all cached values to their initial state.
        /// Use this if you want to force recalculation of all parameters.
        /// </summary>
        public void Reset()
        {
            _lastDay = double.NaN;
            _lastUtsec = double.NaN;
            _lastLat = double.NaN;
            _lastLon = double.NaN;
            _lastAltitude = double.NaN;
            _lastSflux = double.NaN;
            _lastSfluxavg = double.NaN;
            Array.Fill(_lastAp, double.NaN);

            Array.Clear(_basisFunctions, 0, _basisFunctions.Length);
            Array.Clear(_splineWeights, 0, _splineWeights.Length);
            _splineIndex = 0;

            _temperatureProfile = new TemperatureProfile();
            for (int i = 0; i < _densityProfiles.Length; i++)
            {
                _densityProfiles[i] = new DensityParameters();
            }
        }

        /// <summary>
        /// Disposes resources held by this context.
        /// Currently this class doesn't hold unmanaged resources, but implementing
        /// IDisposable allows for future enhancements and follows best practices.
        /// </summary>
        public void Dispose()
        {
            // Currently nothing to dispose, but good practice for future
            // This allows for future enhancements like pooled memory, etc.
        }
    }
}