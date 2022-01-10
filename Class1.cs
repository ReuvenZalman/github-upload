using System;
using System.IO;

namespace clMatlab
{
    /// <summary>
    ///  This class performs an important function.
    /// </summary>
    public class clMatlab
    {
        /// <summary>
        /// this will be the tooltip
        /// </summary>
        /// <param name="path">Big brain time.</param>
        public static int[,] getMatrixFromPath(string path, string sp,int rowOffset = 0)
        {
            int[,] M = null;
            string[] rows = File.ReadAllLines(path);
            string[] row = rows[rowOffset].Split(sp);
            row = new string[row.Length];
            M = new int[rows.Length - rowOffset, row.Length - 1];
            for (int i = rowOffset; i < rows.Length; i++)
            {
                row = rows[i].Split(sp);
                for (int j = 0; j < row.Length - 1; j++) 
                {
                    M[i - rowOffset, j] = Convert.ToInt32(row[j]);
                }
            }
            return M;
        }
        public static void printMatrix(int[] M)
        {
            int H = M.Length;
            Console.WriteLine("Printing Matrix Length = {0}", H);
            for (int i = 0; i < H; i++)
            {
                Console.Write("{0,5}: ", i);
                Console.Write("{0} ", M[i]);
                Console.WriteLine("");
            }
        }
        public static void printMatrix(int[,] M)
        {
            int H = M.GetLength(0);
            int W = M.GetLength(1);
            Console.WriteLine("Printing Matrix Size = {0} x {1} = {2}", H, W, H * W);
            for (int i = 0; i < H; i++)
            {
                Console.Write("{0,5}: ", i);
                int j = 0;
                for (; j < W - 1; j++) 
                {
                    Console.Write("{0} ",M[i, j]);
                }
                Console.WriteLine("{0}", M[i, j]);
            }
            Console.WriteLine("");
        }
        public static void printMatrix(double[] M)
        {
            int W = M.Length;
            Console.WriteLine("Printing Vector of length = {0}",W);
            int i = 0;
            for (; i < W-1; i++)
            {
                Console.Write("{0} ", M[i]);
            }
            Console.WriteLine("{0}\n", M[i]);
        }
        public static void printMatrix(double[,] M)
        {
            int H = M.GetLength(0);
            int W = M.GetLength(1);
            Console.WriteLine("Printing Matrix Size = {0} x {1} = {2}", H, W, H * W);
            for (int i = 0; i < H; i++)
            {
                Console.Write("{0,5}: ", i);
                int j = 0;
                for (; j < W - 1; j++)
                {
                    Console.Write("{0} ", M[i, j]);
                }
                Console.WriteLine("{0}", M[i, j]);
            }
        }
        public static void printMatrixLean(int[,] M)
        {
            int H = M.GetLength(0);
            int W = M.GetLength(1);
            int R = 5;
            int C = 5;

            Console.WriteLine("Printing Matrix Size = {0} x {1} = {2}\n", H, W, H * W);
            int i = 0;
            for (; i < R + 1; i++)
            {
                Console.Write("{0,10}: ", i);
                int j = 0;
                for (; j < C + 1; j++)
                {
                    Console.Write("{0,10} ", M[i, j]);
                }
                Console.Write("... ");
                j = W - 1 - C;
                for (; j < W - 1; j++)
                {
                    Console.Write("{0,10} ", M[i, j]);
                }
                Console.WriteLine("{0,10}", M[i, j]);
            }
            for (int k = 0; k < 3; k++)
            {
                Console.WriteLine("{0,22}{0,11}{0,11}", '.');
            }
            i = H - 1 - R;
            for (; i < H; i++)
            {
                Console.Write("{0,10}: ", i);
                int j = 0;
                for (; j < C + 1; j++)
                {
                    Console.Write("{0,10} ", M[i, j]);
                }
                Console.Write("... ");
                j = W - 1 - C;
                for (; j < W - 1; j++)
                {
                    Console.Write("{0,10} ", M[i, j]);
                }
                Console.WriteLine("{0,10}", M[i, j]);
            }
            Console.WriteLine("\n");
        }
        public static void printMatrixLean(double[,] M)
        {
            int H = M.GetLength(0);
            int W = M.GetLength(1);
            if (H <= 10 || W <= 10) 
            {
                clMatlab.printMatrix(M);
                return;
            }
            int R = 5;
            int C = 5;

            Console.WriteLine("Printing Matrix Size = {0} x {1} = {2}\n", H, W, H * W); int i = 0;
            for (; i < R + 1; i++)
            {
                Console.Write("{0,10}: ", i);
                int j = 0;
                for (; j < C + 1; j++)
                {
                    Console.Write("{0,10} ", M[i, j]);
                }
                Console.Write("... ");
                j = W - 1 - C;
                for (; j < W - 1; j++)
                {
                    Console.Write("{0,10} ", M[i, j]);
                }
                Console.WriteLine("{0,10}", M[i, j]);
            }
            for (int k = 0; k < 3; k++)
            {
                Console.WriteLine("{0,22}{0,11}{0,11}", '.');
            }
            i = H - 1 - R;
            for (; i < H; i++)
            {
                Console.Write("{0,10}: ", i);
                int j = 0;
                for (; j < C + 1; j++)
                {
                    Console.Write("{0,10} ", M[i, j]);
                }
                Console.Write("... ");
                j = W - 1 - C;
                for (; j < W - 1; j++)
                {
                    Console.Write("{0,10} ", M[i, j]);
                }
                Console.WriteLine("{0,10}", M[i, j]);
            }
            Console.WriteLine("\n");
        }
        public static double MatrixMean(int[,] M)
        {
            double mm = 0;
            int H = M.GetLength(0);
            int W = M.GetLength(1);
            for (int i = 0; i < H; i++)
            {
                for (int j = 0; j < W; j++)
                {
                    mm += M[i, j];
                }
            }
            mm /= (W*H);
            return mm;
        }
        public static double MatrixMean(double[,] M)
        {
            double mm = 0;
            int H = M.GetLength(0);
            int W = M.GetLength(1);
            for (int i = 0; i < H; i++)
            {
                for (int j = 0; j < W; j++)
                {
                    mm += M[i, j];
                }
            }
            mm /= Convert.ToDouble(W * H);
            return mm;
        }
        public static int MatrixIdxOfMin(double[] M)
        {
            //Needs generalization for any n!
            int W = M.Length;
            int minIdx = 0;
            for (int j = 1; j < W; j++)
            {
                if (M[j] < M[minIdx])
                {
                    minIdx = j;
                }
            }
            return minIdx;
        }
        public static int[] MatrixIdxOfMin(double[] M, int n)
        {
            //Needs generalization for any n!
            int W = M.Length;
            int[] minIdx = new int[n];
            minIdx[0] = 0;
            for (int j = 1; j < W; j++)
            {
                if (M[j] < M[minIdx[0]])
                {
                    minIdx[0] = j;
                }
            }
            return minIdx;
        }
        public static int[,] MatrixIdxOfMin(double[,] M, int n)
        {
            //Needs generalization for any n!
            int H = M.GetLength(0);
            int W = M.GetLength(1);
            int[,] minIdx = new int[n, 2];
            minIdx[0, 0] = 0;
            minIdx[0, 1] = 0;
            for (int i = 0; i < H; i++)
            {
                for (int j = 1; j < W; j++)
                {
                    if (M[i, j] < M[minIdx[0, 0], minIdx[0, 1]])
                    {
                        minIdx[0, 0] = i;
                        minIdx[0, 1] = j;
                    }
                }
            }
            return minIdx;
        }
        public static double[,] ConvertIntArrayToDouble(int[,] M)
        {
            int H = M.GetLength(0);
            int W = M.GetLength(1);
            double[,] D = new double[H, W];
            for (int i = 0; i < H; i++)
            {
                for (int j = 0; j < W; j++)
                {
                    D[i, j] = Convert.ToDouble(M[i, j]);
                }
            }
            return D;
        }
        public static void MatrixAdd(double[] M, double val)
        {
            int W = M.Length;
            for (int j = 0; j < W; j++)
            {
                M[j] += val;
            }
            return;
        }
        public static void MatrixAdd(double[,] M, double val)
        {
            int H = M.GetLength(0);
            int W = M.GetLength(1);
            for (int i = 0; i < H; i++)
            {
                for (int j = 0; j < W; j++)
                {
                    M[i,j] += val;
                }
            }
            return;
        }
        public static void MatrixMultiply(int[,] M, int val)
        {
            int H = M.GetLength(0);
            int W = M.GetLength(1);
            for (int i = 0; i < H; i++)
            {
                for (int j = 0; j < W; j++)
                {
                    M[i, j] *= val;
                }
            }
            return;
        }
        public static void MatrixMultiply(double[,] M, double val)
        {
            int H = M.GetLength(0);
            int W = M.GetLength(1);
            for (int i = 0; i < H; i++)
            {
                for (int j = 0; j < W; j++)
                {
                    M[i, j] *= val;
                }
            }
            return;
        }
        public static void MatrixMultiply(double[] M, double val)
        {
            int W = M.Length;
            for (int j = 0; j < W; j++)
            {
                M[j] *= val;
            }
            return;
        }
        public static void ZeroBelowThreshold(int[,] M, int th = 0)
        {
            int H = M.GetLength(0);
            int W = M.GetLength(1);
            for (int i = 0; i < H; i++)
            {
                for (int j = 0; j < W; j++)
                {
                    M[i, j] = (M[i, j] <= th) ? 0 : M[i, j];
                }
            }
            return;
        }
        public static void ZeroBelowThreshold(double[,] M, double th = 0)
        {
            int H = M.GetLength(0);
            int W = M.GetLength(1);
            for (int i = 0; i < H; i++)
            {
                for (int j = 0; j < W; j++)
                {
                    M[i, j] = (M[i, j] <= th) ? 0 : M[i, j];
                }
            }
            return;
        }
        public static int[,] FindAboveThreshold(double[,] M, double th = 0)
        {
            int H = M.GetLength(0);
            int W = M.GetLength(1);
            int k = 0;
            int[,] ijs = new int[H * W, 2];
            for (int i = 0; i < H; i++)
            {
                for (int j = 0; j < W; j++)
                {
                    if (M[i, j] > th) 
                    {
                        ijs[k, 0] = i;
                        ijs[k, 1] = j;
                        k++;
                    }
                }
            }
            int[,] ijs0 = new int[k, 2];
            for (int i = 0; i < k; i++)
            {
                ijs0[i, 0] = ijs[i, 0];
                ijs0[i, 1] = ijs[i, 1];
            }
            return ijs0;
        }
        public static int[,] GetPartOfMatrix(int[,] M, int ii, int ie, int ji, int je)
        {
            int[,] Mout = new int[ie - ii + 1, je - ji + 1];
            for (int i = 0; i < Mout.GetLength(0); i++)
            {
                for (int j = 0; j < Mout.GetLength(1); j++)
                {
                    Mout[i, j] = M[ii + i, ji + j];
                }
            }
            return Mout;
        }
        public static double[,] GetPartOfMatrix(double[,] M, int ii,int ie, int ji, int je)
        {
            double[,] Mout = new double[ie - ii + 1, je - ji + 1];
            for (int i = 0; i < Mout.GetLength(0); i++)
            {
                for (int j = 0; j < Mout.GetLength(1); j++)
                {
                    Mout[i, j] = M[ii + i, ji + j];
                }
            }
            return Mout;
        }
        public static int[] sumOfRows(int[,] M)
        {
            int[] hProf = new int[M.GetLength(1)];
            for (int j = 0; j < hProf.Length; j++)
            {
                for (int i = 0; i < M.GetLength(0); i++)
                {
                    hProf[j] += M[i, j];
                }
            }
            return hProf;
        }
        public static double[] sumOfRows(double[,] M)
        {
            double[] hProf = new double[M.GetLength(1)];
            for (int j = 0; j < hProf.Length; j++)
            {
                for (int i = 0; i < M.GetLength(0); i++)
                {
                    hProf[j] += M[i, j];
                }
            }
            return hProf;
        }
        public static int[] sumOfCols(int[,] M)
        {
            int[] vProf = new int[M.GetLength(0)];
            for (int i = 0; i < vProf.Length; i++)
            {
                for (int j = 0; j < M.GetLength(1); j++)
                {
                    vProf[i] += M[i, j];
                }
            }
            return vProf;
        }
        public static double[] sumOfCols(double[,] M)
        {
            double[] vProf = new double[M.GetLength(0)];
            for (int i = 0; i < vProf.Length; i++)
            {
                for (int j = 0; j < M.GetLength(1); j++)
                {
                    vProf[i] += M[i, j];
                }
            }
            return vProf;
        }
        public static double[] rcFilter(double[] V, double k = 0.5)
        {
            double[] Vout = new double[V.Length];
            Vout[0] = (1 - k) * V[0];
            for (int i = 1; i < V.Length; i++)
            {
                Vout[i] = k * Vout[i - 1] + (1 - k) * V[i];
            }
            return Vout;
        }
        public static double[] flipArray(double[] V)
        {
            double[] Vout = new double[V.Length];
            for (int i = 0; i < V.Length; i++)
            {
                Vout[i] = V[V.Length - i - 1];
            }
            return Vout;
        }
        public static int centroid(double[] V)
        {
            double loc = 0;
            double sum = 0;
            for (int i = 0; i < V.Length; i++)
            {
                loc += (i + 1) * V[i];
                sum += V[i];
            }
            return Convert.ToInt32(Math.Round(loc / sum) - 1);
        }
        public static int[,] maxk(int[,] M, int n = 1)
        {
            int[] maxVals = new int[n];
            int[,] maxIdxs = new int[n, 2];
            int H = M.GetLength(0);
            int W = M.GetLength(1);
            int val = 0;
            for (int i = 0; i < H; i++)
            {
                for (int j = 0; j < W; j++)
                {
                    val = M[i, j];
                    if (val > maxVals[n - 1])
                    {
                        int k = 0;
                        for (; k < n; k++)
                        {
                            if (val > maxVals[k])
                            {
                                break;
                            }
                        }
                        for (int s = n - 1; s > k; s--)
                        {
                            maxVals[s] = maxVals[s - 1];
                            maxIdxs[s, 0] = maxIdxs[s - 1, 0];
                            maxIdxs[s, 1] = maxIdxs[s - 1, 1];
                        }
                        maxVals[k] = val;
                        maxIdxs[k, 0] = i;
                        maxIdxs[k, 1] = j;
                    }
                }
            }
            return maxIdxs;
        }
        public static int[] maxk(double[] M, int n = 1)
        {
            double[] maxVals = new double[n];
            int[] maxIdxs = new int[n];
            int H = M.Length;
            double val;
            for (int i = 0; i < H; i++)
            {
                val = M[i];
                if (val > maxVals[n - 1])
                {
                    int k = 0;
                    for (; k < n; k++)
                    {
                        if (val > maxVals[k])
                        {
                            break;
                        }
                    }
                    for (int s = n - 1; s > k; s--)
                    {
                        maxVals[s] = maxVals[s - 1];
                        maxIdxs[s] = maxIdxs[s - 1];
                    }
                    maxVals[k] = val;
                    maxIdxs[k] = i;
                }
            }
            return maxIdxs;
        }
        public static int[,] maxk(double[,] M, int n = 1)
        {
            double[] maxVals = new double[n];
            int[,] maxIdxs = new int[n, 2];
            int H = M.GetLength(0);
            int W = M.GetLength(1);
            double val;
            for (int i = 0; i < H; i++)
            {
                for (int j = 0; j < W; j++)
                {
                    val = M[i, j];
                    if (val > maxVals[n - 1]) 
                    {
                        int k = 0;
                        for (; k < n; k++)
                        {
                            if (val > maxVals[k])
                            {
                                break;
                            }
                        }
                        for (int s = n-1; s > k; s--)
                        {
                            maxVals[s] = maxVals[s-1];
                            maxIdxs[s, 0] = maxIdxs[s - 1, 0];
                            maxIdxs[s, 1] = maxIdxs[s - 1, 1];
                        }
                        maxVals[k] = val;
                        maxIdxs[k, 0] = i;
                        maxIdxs[k, 1] = j;
                    }
                }
            }
            return maxIdxs;
        }
        public static int[,] MatrixAddition(int[,] A, int[,] B)
        {
            int H = A.GetLength(0);
            int W = A.GetLength(1);
            int[,] C = new int[H, W];
            if (H!=B.GetLength(0) || W != B.GetLength(1))
            {
                Console.WriteLine("Matrix dimensions aren't equal.");
                return C;
            }
            for (int i = 0; i < H; i++)
            {
                for (int j = 0; j < W; j++)
                {
                    C[i, j] = A[i, j] + B[i, j];
                }
            }
            return C;
                
        }
        public static int[] Unique(int[] M ,int del)
        {
            int N = M.Length;
            int[] idxs = new int[N];
            int[] P = new int[N];
            int k = 0;
            P[k] = M[0];
            idxs[k] = 0;
            for (int i = 1; i < N; i++)
            {
                for (int j = 0; j < k + 1; j++)
                {
                    if (Math.Abs(P[j] - M[i]) < del) 
                    {
                        break;
                    }
                    if (j == k) 
                    {
                        k++;
                        P[k] = M[i];
                        idxs[k] = i;
                    }
                }
            }
            int[] Pout = new int[k + 1];
            for (int i = 0; i < Pout.Length; i++)
            {
                Pout[i] = M[idxs[i]];
            }
            return Pout;
        }
        public static int[,] Unique(int[,] M, int del)
        {
            int N = M.GetLength(0);
            int[] idxs = new int[N];
            int[,] P = new int[N,2];
            int k = 0;
            P[k, 0] = M[0, 0];
            P[k, 1] = M[0, 1];
            idxs[k] = 0;
            for (int i = 1; i < N; i++)
            {
                for (int j = 0; j < k + 1; j++)
                {
                    if (Math.Abs(P[j,0] - M[i,0]) < del && Math.Abs(P[j, 1] - M[i, 1]) < del)
                    {
                        break;
                    }
                    if (j == k)
                    {
                        k++;
                        P[k, 0] = M[i, 0];
                        P[k, 1] = M[i, 1];
                        idxs[k] = i;
                    }
                }
            }
            int[,] Pout = new int[k + 1,2];
            for (int i = 0; i < Pout.GetLength(0); i++)
            {
                Pout[i, 0] = M[idxs[i], 0];
                Pout[i, 1] = M[idxs[i], 1];
            }
            return Pout;
        }
        public static int[,] getPartOfMatrix(int[,] M, int[] centCoo, int d)
        {
            int dd = d / 2;
            int T = centCoo[0] - dd + (1 - d % 2);
            int B = centCoo[0] + dd;
            int L = centCoo[1] - dd + (1 - d % 2);
            int R = centCoo[1] + dd;
            int[,] K = new int[B - T + 1, R - L + 1];
            for (int i = 0; i < K.GetLength(0); i++)
            {
                for (int j = 0; j < K.GetLength(1); j++)
                {
                    K[i, j] = M[T + i, L + j];
                }
            }
            return K;
        }
    }
}
