// Vulcan L0: in-process Thermo RAW -> columnar batches.
// One call per batch of scans returns flat primitive arrays (no per-spectrum Python).
// Semantics are documented in vulcan/io/STORE_SPEC.md (section "Thermo direct").
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Text;
using ThermoFisher.CommonCore.Data;
using ThermoFisher.CommonCore.Data.Business;
using ThermoFisher.CommonCore.Data.FilterEnums;
using ThermoFisher.CommonCore.Data.Interfaces;
using ThermoFisher.CommonCore.RawFileReader;

namespace Vulcan.IO
{
    public sealed class Batch
    {
        public int N;
        public int[] Scan; public byte[] MsLevel; public double[] RtMin; public float[] InjTimeMs;
        public double[] PrecMz; public double[] IsoTargetMz; public double[] MonoMz;
        public float[] IsoLowerOffset; public float[] IsoUpperOffset;
        public short[] Charge; public int[] MasterScan;
        public byte[] Activation; public float[] Nce; public byte[] SuppActivation; public float[] SuppEnergy;
        public byte[] Analyzer; public float[] Resolution; public float[] FaimsCv;
        public double[] Tic; public double[] BasePeakMz; public double[] BasePeakInt;
        public float[] ScanLow; public float[] ScanHigh; public sbyte[] Polarity; public byte[] PeakSource;
        public byte[] FilterUtf8; public int[] FilterOffsets;
        public double[] Sps; public int[] SpsOffsets;
        public long NPeaks; public float[] Mz; public float[] Intensity; public long[] PeakOffsets;
    }

    public sealed class RawExtractor : IDisposable
    {
        // activation codes (keep in sync with vulcan/io/schema.py ACTIVATION)
        public const byte ACT_NONE = 0, ACT_CID = 1, ACT_HCD = 2, ACT_ETD = 3, ACT_ECD = 4, ACT_MPD = 5,
            ACT_PQD = 6, ACT_UVPD = 7, ACT_NETD = 8, ACT_PTR = 9, ACT_OTHER = 255;
        // analyzer codes (schema.py ANALYZER)
        public const byte AN_UNKNOWN = 0, AN_OT = 1, AN_IT = 2, AN_TOF = 3, AN_ASTRAL = 4, AN_QUAD = 5, AN_SECTOR = 6;

        readonly IRawDataPlus raw;
        public readonly int FirstScan, LastScan;
        public readonly string InstrumentModel, InstrumentName, InstrumentSerial, SoftwareVersion, CreationDate;
        readonly int tInj = -1, tMaster = -1, tCharge = -1, tMono = -1, tRes = -1, tFaims = -1, tSps = -1, tSps2 = -1, tFaimsOn = -1;
        readonly int[] lastScanOfLevel = new int[16];
        public bool IncludeReferencePeaks = true;   // TRFP default (lock-mass/exception peaks kept)
        public bool IncludeFilterString = true;
        // M2b: one GetTrailerExtraValues(scan) call per scan instead of one GetTrailerExtraValue(scan, i) per field
        // (same typed objects; store bit-identity checked by scripts/m2b/ingest_identity.py).
        public bool FastTrailer = true;
        // M2b profiling: Stopwatch ticks per section (0 stats, 1 event+reactions, 2 trailer, 3 filter object,
        // 4 filter string, 5 peak fetch, 6 peak copy, 7 other); ProfTicksPerSecond = Stopwatch.Frequency.
        public bool Profile = false;
        public readonly long[] ProfTicks = new long[8];
        public static readonly long ProfTicksPerSecond = System.Diagnostics.Stopwatch.Frequency;
        object[] trv;   // current scan's trailer values (FastTrailer)

        public RawExtractor(string path)
        {
            raw = RawFileReaderAdapter.FileFactory(path);
            if (raw == null || raw.IsError) throw new Exception("cannot open RAW: " + path + " " + (raw?.FileError?.ErrorMessage ?? ""));
            if (raw.InAcquisition) throw new Exception("RAW file still in acquisition: " + path);
            raw.SelectInstrument(Device.MS, 1);
            FirstScan = raw.RunHeaderEx.FirstSpectrum;
            LastScan = raw.RunHeaderEx.LastSpectrum;
            var inst = raw.GetInstrumentData();
            InstrumentModel = inst.Model ?? ""; InstrumentName = inst.Name ?? ""; InstrumentSerial = inst.SerialNumber ?? "";
            SoftwareVersion = inst.SoftwareVersion ?? "";
            CreationDate = raw.FileHeader.CreationDate.ToString("o", CultureInfo.InvariantCulture);
            var hdr = raw.GetTrailerExtraHeaderInformation();
            for (int i = 0; i < hdr.Length; i++)
            {
                string l = hdr[i].Label.Trim().TrimEnd(':').Trim();
                switch (l)
                {
                    case "Ion Injection Time (ms)": tInj = i; break;
                    case "Master Scan Number": tMaster = i; break;
                    case "Charge State": tCharge = i; break;
                    case "Monoisotopic M/Z": tMono = i; break;
                    case "FT Resolution": case "Orbitrap Resolution": if (tRes < 0) tRes = i; break;
                    case "FAIMS CV": tFaims = i; break;
                    case "FAIMS Voltage On": tFaimsOn = i; break;
                    case "SPS Masses": tSps = i; break;
                    case "SPS Masses Continued": tSps2 = i; break;
                }
            }
            for (int i = 0; i < lastScanOfLevel.Length; i++) lastScanOfLevel[i] = -1;
        }

        public void Dispose() { raw?.Dispose(); }

        static byte MapActivation(ActivationType a)
        {
            switch (a)
            {
                case ActivationType.CollisionInducedDissociation: return ACT_CID;
                case ActivationType.HigherEnergyCollisionalDissociation: return ACT_HCD;
                case ActivationType.ElectronTransferDissociation: return ACT_ETD;
                case ActivationType.ElectronCaptureDissociation: return ACT_ECD;
                case ActivationType.MultiPhotonDissociation: return ACT_MPD;
                case ActivationType.PQD: return ACT_PQD;
                case ActivationType.UltraVioletPhotoDissociation: return ACT_UVPD;
                case ActivationType.NegativeElectronTransferDissociation: return ACT_NETD;
                case ActivationType.ProtonTransferReaction: return ACT_PTR;
                default: return ACT_OTHER;
            }
        }

        static byte MapAnalyzer(MassAnalyzerType m)
        {
            switch (m)
            {
                case MassAnalyzerType.MassAnalyzerFTMS: return AN_OT;
                case MassAnalyzerType.MassAnalyzerITMS: return AN_IT;
                case MassAnalyzerType.MassAnalyzerTOFMS: return AN_TOF;
                case MassAnalyzerType.MassAnalyzerASTMS: return AN_ASTRAL;
                case MassAnalyzerType.MassAnalyzerSQMS: case MassAnalyzerType.MassAnalyzerTQMS: return AN_QUAD;
                case MassAnalyzerType.MassAnalyzerSector: return AN_SECTOR;
                default: return AN_UNKNOWN;
            }
        }

        static double ParseD(object o, double dflt)
        {
            if (o == null) return dflt;
            if (o is double d) return d;
            if (o is float f) return f;
            if (o is int i) return i;
            if (o is short s) return s;
            if (o is long l) return l;
            if (o is byte b) return b;
            if (o is bool bo) return bo ? 1 : 0;
            var str = o.ToString().Trim();
            return double.TryParse(str, NumberStyles.Float, CultureInfo.InvariantCulture, out var v) ? v : dflt;
        }

        object Tr(int scan, int idx)
        {
            if (idx < 0) return null;
            if (FastTrailer) return trv != null && idx < trv.Length ? trv[idx] : null;
            return raw.GetTrailerExtraValue(scan, idx);
        }

        long tp;  // profiling timestamp
        void Lap(int sec)
        {
            if (!Profile) return;
            long t = System.Diagnostics.Stopwatch.GetTimestamp();
            ProfTicks[sec] += t - tp;
            tp = t;
        }

        static void Grow<T>(ref T[] a, long need)
        {
            if (need <= a.Length) return;
            long n = Math.Max(need, (long)(a.Length * 1.5) + 1024);
            Array.Resize(ref a, (int)n);
        }

        public Batch Read(int first, int last)
        {
            if (first < FirstScan) first = FirstScan;
            if (last > LastScan) last = LastScan;
            int n = Math.Max(0, last - first + 1);
            var b = new Batch
            {
                N = n, Scan = new int[n], MsLevel = new byte[n], RtMin = new double[n], InjTimeMs = new float[n],
                PrecMz = new double[n], IsoTargetMz = new double[n], MonoMz = new double[n],
                IsoLowerOffset = new float[n], IsoUpperOffset = new float[n], Charge = new short[n], MasterScan = new int[n],
                Activation = new byte[n], Nce = new float[n], SuppActivation = new byte[n], SuppEnergy = new float[n],
                Analyzer = new byte[n], Resolution = new float[n], FaimsCv = new float[n], Tic = new double[n],
                BasePeakMz = new double[n], BasePeakInt = new double[n], ScanLow = new float[n], ScanHigh = new float[n],
                Polarity = new sbyte[n], PeakSource = new byte[n], FilterOffsets = new int[n + 1], SpsOffsets = new int[n + 1],
                PeakOffsets = new long[n + 1]
            };
            var mz = new float[Math.Max(1024, n * 256)]; var it = new float[mz.Length];
            var filt = new byte[Math.Max(1024, n * 64)]; int filtLen = 0;
            var sps = new List<double>();
            long np = 0;
            var utf8 = Encoding.UTF8;
            for (int k = 0; k < n; k++)
            {
                int s = first + k;
                b.Scan[k] = s;
                if (Profile) tp = System.Diagnostics.Stopwatch.GetTimestamp();
                var stats = raw.GetScanStatsForScanNumber(s);
                Lap(0);
                var ev = raw.GetScanEventForScanNumber(s);
                Lap(1);
                trv = FastTrailer ? raw.GetTrailerExtraValues(s) : null;
                int lvl = (int)ev.MSOrder; if (lvl < 1) lvl = 1; // MSOrderType: Ms=1, Ms2=2 ... (negative values = special scans)
                b.MsLevel[k] = (byte)lvl;
                b.RtMin[k] = stats.StartTime;
                b.Tic[k] = stats.TIC; b.BasePeakMz[k] = stats.BasePeakMass; b.BasePeakInt[k] = stats.BasePeakIntensity;
                b.ScanLow[k] = (float)stats.LowMass; b.ScanHigh[k] = (float)stats.HighMass;
                b.Polarity[k] = (sbyte)(ev.Polarity == PolarityType.Negative ? -1 : 1);
                b.Analyzer[k] = MapAnalyzer(ev.MassAnalyzer);
                b.InjTimeMs[k] = (float)ParseD(Tr(s, tInj), double.NaN);
                b.Resolution[k] = (float)ParseD(Tr(s, tRes), double.NaN);
                Lap(2);
                // FAIMS: trailer value if present and FAIMS on; else scan-event CV
                float cv = float.NaN;
                var flt = raw.GetFilterForScanNumber(s);
                if (flt.CompensationVoltage == TriState.On && flt.CompensationVoltageCount > 0)
                    cv = (float)flt.CompensationVoltageValue(0);
                else if (tFaims >= 0)
                {
                    bool on = tFaimsOn < 0 || ParseD(Tr(s, tFaimsOn), 0) != 0 || (Tr(s, tFaimsOn)?.ToString().Trim().ToLower() == "on");
                    if (on) { double v = ParseD(Tr(s, tFaims), double.NaN); if (!double.IsNaN(v)) cv = (float)v; }
                }
                b.FaimsCv[k] = cv;
                Lap(3);

                b.PrecMz[k] = double.NaN; b.IsoTargetMz[k] = double.NaN; b.MonoMz[k] = double.NaN;
                b.IsoLowerOffset[k] = float.NaN; b.IsoUpperOffset[k] = float.NaN;
                b.Nce[k] = float.NaN; b.SuppEnergy[k] = float.NaN; b.MasterScan[k] = -1; b.Charge[k] = 0;
                if (lvl >= 2)
                {
                    int nr = ev.MassCount;
                    int pi = Math.Min(lvl - 2, nr - 1);
                    if (pi >= 0)
                    {
                        var r = ev.GetReaction(pi);
                        double w = r.IsolationWidth;
                        double target = r.PrecursorMass;
                        b.IsoTargetMz[k] = target;
                        if (w > 0) { b.IsoLowerOffset[k] = (float)(w / 2 - r.IsolationWidthOffset); b.IsoUpperOffset[k] = (float)(w / 2 + r.IsolationWidthOffset); }
                        b.Activation[k] = MapActivation(r.ActivationType);
                        b.Nce[k] = (float)r.CollisionEnergy;
                        if (ev.SupplementalActivation == TriState.On && nr > pi + 1)
                        {
                            var r2 = ev.GetReaction(pi + 1);
                            b.SuppActivation[k] = MapActivation(r2.ActivationType);
                            b.SuppEnergy[k] = (float)r2.CollisionEnergy;
                        }
                    }
                    double mono = ParseD(Tr(s, tMono), 0);
                    if (mono > 0) b.MonoMz[k] = mono;
                    b.PrecMz[k] = mono > 0 ? mono : b.IsoTargetMz[k];
                    b.Charge[k] = (short)ParseD(Tr(s, tCharge), 0);
                    int ms = tMaster >= 0 ? (int)ParseD(Tr(s, tMaster), -1) : -1;
                    if (ms <= 0) ms = lastScanOfLevel[Math.Min(lvl - 1, 15)];
                    b.MasterScan[k] = ms;
                    if (tSps >= 0)
                    {
                        string sp = (Tr(s, tSps)?.ToString() ?? "") + (tSps2 >= 0 ? "," + (Tr(s, tSps2)?.ToString() ?? "") : "");
                        foreach (var tok in sp.Split(new[] { ',', ' ' }, StringSplitOptions.RemoveEmptyEntries))
                            if (double.TryParse(tok, NumberStyles.Float, CultureInfo.InvariantCulture, out var v) && v > 0) sps.Add(v);
                    }
                }
                b.SpsOffsets[k + 1] = sps.Count;
                lastScanOfLevel[Math.Min(lvl, 15)] = s;
                Lap(1);

                if (IncludeFilterString)
                {
                    var fs = flt.ToString();
                    int bc = utf8.GetByteCount(fs);
                    Grow(ref filt, filtLen + bc);
                    filtLen += utf8.GetBytes(fs, 0, fs.Length, filt, filtLen);
                }
                b.FilterOffsets[k + 1] = filtLen;
                Lap(4);

                // peaks: centroid stream if present (FTMS/Astral), else centroided segmented scan, else Thermo centroider on profile
                double[] pm = null, pin = null; int L = 0; byte src = 0;
                var cs = raw.GetCentroidStream(s, IncludeReferencePeaks);
                if (cs != null && cs.Length > 0) { pm = cs.Masses; pin = cs.Intensities; L = cs.Length; src = 0; }
                else if (stats.IsCentroidScan)
                {
                    var seg = raw.GetSegmentedScanFromScanNumber(s, stats);
                    pm = seg.Positions; pin = seg.Intensities; L = pm.Length; src = 1;
                }
                else
                {
                    var sc = raw.GetSimplifiedCentroids(s);
                    pm = sc.Masses; pin = sc.Intensities; L = pm.Length; src = 2;
                }
                b.PeakSource[k] = src;
                Lap(5);
                Grow(ref mz, np + L); Grow(ref it, np + L);
                for (int j = 0; j < L; j++) { mz[np + j] = (float)pm[j]; it[np + j] = (float)pin[j]; }
                np += L;
                b.PeakOffsets[k + 1] = np;
                Lap(6);
            }
            b.NPeaks = np; b.Mz = mz; b.Intensity = it; // arrays may be longer than NPeaks
            b.FilterUtf8 = filt; b.Sps = sps.ToArray();
            return b;
        }
    }
}
