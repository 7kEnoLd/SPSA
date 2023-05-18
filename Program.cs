using static System.Formats.Asn1.AsnWriter;

namespace SPSA
{
    internal class Program
    {
        static void Main(string[] args)
        {
            int[] count = new int[7];
            double[] x11 = GetX(scope);
            for (int l = 0; l < iter; l++)
            {
                for (int i = 0; i < scope + 2; i++)
                {
                    x11[i] = BCD.FuncDescend(i, scope, x11);
                }
                if (GetTol(scope, x11))
                {
                    count[0] = l; break;
                }
            }

            double[] x22 = GetX(scope);
            for (int l = 0; l < iter; l++)
            {
                x22 = SPSA.CalSPSA(scope, x22, l);
                if (GetTol(scope, x22))
                {
                    count[1] = l; break;
                }
            }

            double[] x33 = GetX(scope);
            for (int l = 0; l < iter; l++)
            {
                if (l == 0)
                {
                    SPSA.df1 = Program.GetFuncDen(scope, x33);
                    SPSA.x_1 = x33;
                    x33 = SPSA.CalSPSA(scope, x33, l);
                }
                else
                {
                    x33 = SPSA.CalSPSABB(scope, x33, l);
                }
                if (GetTol(scope, x33))
                {
                    count[2] = l; break;
                }
            }

            double[] x44 = GetX(scope);
            for (int l = 0; l < iter; l++)
            {
                SPSA_E.xMin = x44;
                SPSA_E.fMin = Program.GetFuncNum(scope, x44);

                x44 = SPSA_E.CalSPSA(scope, x44, l);
                double f = Program.GetFuncNum(scope, x44);
                if (f < SPSA_E.fMin)
                {
                    SPSA_E.fMin = f;
                    SPSA_E.xMin = x44;
                }
                SPSA_E.m1++;
                SPSA_E.m2++;
                if (SPSA_E.m2 == SPSA_E.m2_Max)
                {
                    x44 = SPSA_E.xMin;
                    SPSA_E.m2 = 0;
                }
                if (GetTol(scope, x44))
                {
                    count[3] = l; break;
                }
            }

            SPSA_E.m1 = 0;
            SPSA_E.m2 = 0;

            double[] x55 = GetX(scope);
            for (int l = 0; l < iter; l++)
            {
                if (l == 0)
                {
                    SPSA_E.xMin = x55;
                    SPSA_E.fMin = Program.GetFuncNum(scope, x55);
                    SPSA_E.df1 = Program.GetFuncDen(scope, x55);
                    SPSA_E.x_1 = x55;
                    x55 = SPSA_E.CalSPSA(scope, x55, l);
                }
                else
                {
                    x55 = SPSA_E.CalSPSABB(scope, x55, l);
                    double f = Program.GetFuncNum(scope, x55);
                    if (f < SPSA_E.fMin)
                    {
                        SPSA_E.fMin = f;
                        SPSA_E.xMin = x55;
                    }
                    SPSA_E.m1++;
                    SPSA_E.m2++;
                    if (SPSA_E.m2 == SPSA_E.m2_Max)
                    {
                        x55 = SPSA_E.xMin;
                        SPSA_E.m2 = 0;
                    }
                }
                if (GetTol(scope, x55))
                {
                    count[4] = l; break;
                }
            }

            //阻尼牛顿法
            double[] x66 = GetX(scope);
            for (int l = 0; l < iter; l++)
            {
                if (l == 0)
                {
                    x66 = ZNewt.CalZNewt(scope, x66);
                }
                else
                {
                    x66 = ZNewt.CalZNewt(scope, x66);
                }
                if (GetTol(scope, x66))
                {
                    count[5] = l; break;
                }
            }


            //共轭梯度法
            double[] x77 = GetX(scope);
            for (int l = 0; l < iter; l++)
            {
                if (l == 0)
                {
                    x77 = ConjGrads.CalConjGrads(scope, x77, l);
                }
                else
                {
                    x77 = ConjGrads.CalConjGrads(scope, x77, l);
                }
                if (GetTol(scope, x77))
                {
                    count[6] = l; break;
                }
            }

            for (int i = 0; i < fig.Count; i++)
            {
                if (fig[i] > tol)
                {
                    Console.WriteLine(fig[i]);
                }
                else
                {
                    Console.WriteLine(fig[i]);
                    Console.WriteLine();
                }
            }
        }

        #region 底层函数
        public static int scope = 1000;

        static List<double> fig = new List<double>();

        static int iter = 100000;

        static double tol = 0.01;

        public static double GetFuncNum(int scope, double[] x)
        {
            double f = 0;
            for (int i = 0; i < scope; i++)
            {
                f += Math.Pow(-x[i] * x[i] + x[i + 1] * x[i + 1] + x[i + 2] * x[i + 2], 2) + Math.Pow(x[i] * x[i] - x[i + 1] * x[i + 1] + x[i + 2] * x[i + 2], 2) + Math.Pow(x[i] * x[i] + x[i + 1] * x[i + 1] - x[i + 2] * x[i + 2], 2);
            }
            return f;
        }

        public static double[] GetFuncDen(int scope, double[] x)
        {
            double[] f = new double[scope + 2];
            for (int i = 0; i < scope; i++)
            {
                f[i] += 4 * x[i] * (3 * x[i] * x[i] - x[i + 1] * x[i + 1] - x[i + 2] * x[i + 2]);
                f[i + 1] += 4 * x[i + 1] * (3 * x[i + 1] * x[i + 1] - x[i] * x[i] - x[i + 2] * x[i + 2]);
                f[i + 2] += 4 * x[i + 2] * (3 * x[i + 2] * x[i + 2] - x[i + 1] * x[i + 1] - x[i] * x[i]);
            }
            return f;
        }

        public static double[,] GetFuncDen2(int scope, double[] x)
        {
            double[,] f = new double[scope + 2, scope + 2];
            for (int i = 0; i < scope; i++)
            {
                f[i, i] += 36 * x[i] * x[i] - 4 * x[i + 1] * x[i + 1] - 4 * x[i + 2] * x[i + 2];
                f[i, i + 1] += -8 * x[i] * x[i + 1];
                f[i, i + 2] += -8 * x[i] * x[i + 2];
                f[i + 1, i] += -8 * x[i] * x[i + 1];
                f[i + 1, i + 1] += 36 * x[i + 1] * x[i + 1] - 4 * x[i] * x[i] - 4 * x[i + 2] * x[i + 2];
                f[i + 1, i + 2] += -8 * x[i + 1] * x[i + 2];
                f[i + 2, i] += -8 * x[i] * x[i + 2]; ;
                f[i + 2, i + 1] += -8 * x[i + 1] * x[i + 2];
                f[i + 2, i + 2] += 36 * x[i + 2] * x[i + 2] - 4 * x[i + 1] * x[i + 1] - 4 * x[i] * x[i];
            }
            return f;
        }

        static double[] GetX(int scope)
        {
            double[] x = new double[scope + 2];
            if (scope % 3 == 1)
            {
                int count = 0;
                for (int i = 0; i < scope + 2; i++)
                {
                    switch (count)
                    {
                        case 0:
                            x[i] = 1; break;
                        case 1:
                            x[i] = 2; break;
                        default:
                            x[i] = 1; break;
                    }
                    ++count;
                    if (count == 3) count = 0;
                }
                return x;
            }
            else
            {
                Console.WriteLine("非法变量数");
                return x;
            }
        }

        static bool GetTol(int scope, double[] x)
        {
            List<double> f = new List<double>();
            double tol_cur = 0;
            double[] tol_den = GetFuncDen(scope, x);
            for (int i = 0; i < scope + 2; i++)
            {
                tol_cur += Math.Pow(tol_den[i], 2);
            }
            tol_cur = Math.Sqrt(tol_cur);
            fig.Add(tol_cur);
            return tol_cur < tol;
        }


        //只改变单个值的函数值
        public static double GetFuncNumSingle(int scope, double[] x_ori, double x1, int loca)
        {
            double[] x = new double[scope + 2];
            for (int i = 0; i < x.Length; i++)
            {
                x[i] = x_ori[i];
            }
            x[loca] = x1;
            double f = 0;
            for (int i = 0; i < scope; i++)
            {
                f += Math.Pow(-x[i] * x[i] + x[i + 1] * x[i + 1] + x[i + 2] * x[i + 2], 2) + Math.Pow(x[i] * x[i] - x[i + 1] * x[i + 1] + x[i + 2] * x[i + 2], 2) + Math.Pow(x[i] * x[i] + x[i + 1] * x[i + 1] - x[i + 2] * x[i + 2], 2);
            }
            return f;
        }
        #endregion
    }
}