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
        public static double fs;
        public static double fmin = 60;
        public static double fmax = 6000;
        public static double fratio = 1.0 / 24; //1オクターブの周波数bin数
        public static double q_rate = 20.0 * fratio; //Q比
        public static double Q = (int)(1.0 / (Math.Pow(2, fratio) - 1) * q_rate);

        static void Main(string[] args)
        {
            //double[] freqs = get_freqs(fmin, get_freq_num(fmin, fmax, fratio), fratio);
            //Complex[,] kernel = calcKernel_matrix(freqs);
            //double[,] kernel_pow = new double[kernel.GetLength(0), kernel.GetLength(1)];
            //for (int i = 0; i < kernel.GetLength(0); i++)
            //{
            //    for (int j = 0; j < kernel.GetLength(1); j++)
            //    {
            //        kernel_pow[i, j] = kernel[i, j].Magnitude;
            //    }
            //}
            Tuple<double[][], int> tu = WaveFile.Load(@"C:\Users\優\Desktop\音素材\mix3.wav");
            double[] wavdata = tu.Item1[0];
            fs = tu.Item2;
            double[,] output = constantQ_transformP(wavdata);

            CsvFileIO.CsvFileIO.WriteData("output.csv", output);
        }
        //--------------------------------------------------------------------------
        //対数周波数bin数の計算
        static int get_freq_num(double fmin, double fmax, double fratio)
        {
            return (int)Math.Round(Math.Log((double)fmax / fmin, 2) / fratio) + 1;
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
        static Complex[,] calcKernel_matrix(double[] freqs)
        {
            //窓幅N[k]の最大値
            int fftlen = (int)(Math.Pow(2, Math.Ceiling((Math.Log(fs * Q / freqs[0], 2)))));
            Complex[,] kernel = new Complex[freqs.Length, fftlen];
            Complex[] tmp_kernel = new Complex[fftlen];

            for (int k = 0; k < freqs.Length; k++)
            {
                double freq = freqs[k];
                int N_k = (int)((double)(fs * Q) / freq);     //窓幅
                int start_win = (fftlen - N_k) / 2; //FFT窓の中心を解析部分に合わせる
                double[] hamming = new double[N_k];
                for (int i = 0; i < N_k; i++)
                    hamming[i] = 1;
                Nm.Windowing(hamming, Nm.DataWindowType.Hamming);
                for (int i = start_win; i < start_win + N_k; i++)
                {
                    tmp_kernel[i] = hamming[i - start_win] / N_k * Complex.Exp(2 * Math.PI * Complex.ImaginaryOne * Q * (i - start_win) / N_k);
                }
                tmp_kernel = Nm.FastFourierTransform(tmp_kernel, false); //sw==falseでfftlenで割る
                double[] d = new double[tmp_kernel.Length];
                for (int i = 0; i < tmp_kernel.Length; i++)
                    d[i] = tmp_kernel[i].Magnitude;
                double max = d.Max();
                for (int i = 0; i < fftlen; i++)
                {
                    //if (tmp_kernel[i].Magnitude <= 0.0054)
                    //    kernel[k, i] = Complex.Zero;
                    //else
                    kernel[k, i] = Complex.Conjugate(tmp_kernel[i]);
                }
            }

            return kernel;
        }
        //--------------------------------------------------------------------------
        //constant Q transform
        static Complex[,] constantQ_transform(double[] wavdata)
        {
            double[] new_wavdata;
            Complex[,] kernel;
            Complex[,] output;
            int wav_length = wavdata.Length;
            int T_step;
            int F_step;
            int fft_step;
            int fft_start_point;
            int hop;

            double[] freqs = get_freqs(fmin, get_freq_num(fmin, fmax, fratio), fratio);
            kernel = calcKernel_matrix(freqs);

            hop = (int)(Math.Round(0.01 * fs)); //100 frames per 1 sec
            T_step = wav_length / hop;
            F_step = freqs.Length;
            fft_step = (int)(Math.Pow(2, Math.Ceiling((Math.Log(fs * Q / freqs[0], 2))))); //greater than Max of N[k]
            new_wavdata = new double[wav_length + fft_step];
            output = new Complex[F_step, T_step];

            for (int i = 0; i < wav_length; i++)
                new_wavdata[i + fft_step / 2] = wavdata[i];

            for (int t = 0; t < T_step; t++)
            {
                fft_start_point = hop * t;
                Complex[] partial_data = new Complex[fft_step];
                for (int i = 0; i < fft_step; i++)
                    partial_data[i] = new Complex(new_wavdata[fft_start_point + i], 0);

                Complex[] fft = Nm.FastFourierTransform(partial_data, true);
                Complex[] cq = new Complex[F_step];
                for (int i = 0; i < F_step; i++)
                    for (int j = 0; j < fft.Length; j++)         //fft.Length = fft_step
                        cq[i] += fft[j] * kernel[i, j];

                for (int f = 0; f < F_step; f++)
                    output[f, t] = cq[f];
            }
            return output;
        }
        //--------------------------------------------------------------------------
        //constant Q transform(power)
        static double[,] constantQ_transformP(double[] wavdata)
        {
            double[,] output;
            Complex[,] c = constantQ_transform(wavdata);

            int I = c.GetLength(0), J = c.GetLength(1);
            output = new double[I, J];
            for (int i = 0; i < I; i++)
            {
                for (int j = 0; j < J; j++)
                {
                    output[i, j] = c[i, j].Magnitude * c[i, j].Magnitude;
                }
            }
            return output;
        }

    }
}
