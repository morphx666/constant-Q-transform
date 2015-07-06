using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using InoueLab;
using MySignalProcessing;
using System.Numerics;
using CsvFileIO;

namespace ConstantQtransform
{
    class Program
    {
        public static double fs = 44100;
        public static double Q = 23;
        public static double fmin = 60;
        public static double fmax = 6000;
        public static double fratio = 1.0 / 24; //1オクターブの周波数bin数
        public static double q_rate = 20.0 * fratio; //Q比

        static void Main(string[] args)
        {
            double[] freqs = get_freqs(fmin, get_freq_num(fmin, fmax, fratio), fratio);
            Complex[,] kernel = calcKernel_matrix(freqs, 1024);
            double[,] kernel_pow = new double[kernel.GetLength(0), kernel.GetLength(1)];
            for (int i = 0; i < kernel.GetLength(0); i++)
            {
                for (int j = 0; j < kernel.GetLength(1); j++)
                {
                    kernel_pow[i, j] = kernel[i, j].Magnitude;
                }
            }
            CsvFileIO.CsvFileIO.WriteData("kernel.csv", kernel_pow);
        }
        //--------------------------------------------------------------------------
        //対数周波数bin数の計算
        static int get_freq_num(double fmin, double fmax, double fratio)
        {
            return (int)Math.Round(Math.Log(fmax / fmin, 2) / fratio) + 1;
        }
        //--------------------------------------------------------------------------
        //各対数周波数binに対する周波数を計算
        static double[] get_freqs(double fmin, int nfreq, double fratio)
        {
            double[] freqs = new double[nfreq];
            for (int i = 0; i < freqs.Length; i++)
                freqs[i] = fmin * (Math.Pow(2, i * fratio));
            return freqs;
        }

        //--------------------------------------------------------------------------
        //Kernel行列の計算
        static Complex[,] calcKernel_matrix(double[] freqs, int fftlen)
        {
            Complex[,] kernel = new Complex[freqs.Length, fftlen];
            Complex[] tmp_kernel = new Complex[fftlen];

            for (int k = 0; k < freqs.Length; k++)
            {
                double freq = freqs[k];
                int N_k = (int)(fs * Q / freq);
                int start_win = (fftlen - N_k) / 2; //FFT窓の中心を解析部分に合わせる
                double[] hamming = new double[N_k];
                for (int i = 0; i < N_k; i++)
                    hamming[i] = 1;
                Nm.Windowing(hamming, Nm.DataWindowType.Hamming);
                for (int i = start_win; i < start_win + N_k; i++)
                    tmp_kernel[i] = hamming[i - start_win] / N_k * Complex.Exp(2 * Math.PI * Complex.ImaginaryOne * Q * (i - start_win) / N_k);
                tmp_kernel = Nm.FastFourierTransform(tmp_kernel, false); //sw==falseでfftlenで割る
                for (int i = 0; i < fftlen; i++)
                {
                    if (tmp_kernel[i].Magnitude <= 0.0054)
                        kernel[k, i] = Complex.Zero;
                    else
                        kernel[k, i] = Complex.Conjugate(tmp_kernel[i]);
                }
            }

            return kernel;
        }
    }
}
