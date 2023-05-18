using MathNet.Numerics.Distributions;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SPSA
{
    public class SPSA_E
    {
        private static double alpha = 0.602;
        private static double gamma = 0.101;
        private static double a = 20;
        private static double c = 1.9;
        private static double A = 5000;
        private static int search = 10000;
        private static int scope = Program.scope;

        public static int m1 = 0;
        public static int m2 = 0;
        public static int m1_Max = 1000;
        public static int m2_Max = 100;

        public static double[] xMin = new double[scope + 2];
        public static double fMin = 10000000;

        private static double[] BernoulliGenerator(int scope)
        {
            double[] bernoulli = new double[scope + 2];
            Bernoulli bern = new Bernoulli(0.5);
            for (int i = 0; i < bernoulli.Length; i++)
            {
                bernoulli[i] = bern.Sample();
            }
            return bernoulli;
        }



        public static double[] CalSPSA(int scope, double[] x, int l)
        {
            double[] x1 = new double[scope + 2];
            double[] x2 = new double[scope + 2];
            double[] delta = BernoulliGenerator(scope);
            for (int i = 0; i < scope + 2; i++)
            {
                if (delta[i] == 0)
                {
                    delta[i] = -1;
                }
                x1[i] = x[i] + c / Math.Pow(l + 1, gamma) * delta[i];
                x2[i] = x[i] - c / Math.Pow(l + 1, gamma) * delta[i];
            }

            double[] f1 = new double[scope + 2];
            double[] f2 = new double[scope + 2];
            for (int i = 0; i < scope + 2; i++)
            {
                f1[i] = Program.GetFuncNumSingle(scope, x, x1[i], i) - Program.GetFuncNum(scope, x);
                f2[i] = Program.GetFuncNumSingle(scope, x, x2[i], i) - Program.GetFuncNum(scope, x);
            }
            double[] g = new double[scope + 2];
            for (int i = 0; i < scope + 2; i++)
            {
                g[i] = (f1[i] - f2[i]) / (2 * c / Math.Pow(l + 1, gamma) * delta[i]);
            }
            double[] xm = new double[scope + 2];
            double g_Max = g.Max();
            for (int i = 0; i < scope + 2; i++)
            {
                g[i] = g[i] / g_Max;
                xm[i] = x[i] + a / Math.Pow(l + 1 + A, alpha) * g[i];
            }
            double f_des1;
            double f_des = Program.GetFuncNum(scope, xm);
            for (int i = 0; i < search; i++)
            {
                for (int j = 0; j < scope + 2; j++)
                {
                    xm[j] = x[j] - a / Math.Pow(l + 1 + A, alpha) * g[j] / 2;
                }
                f_des1 = Program.GetFuncNum(scope, xm);
                if (f_des1 < f_des)
                {
                    break;
                }
            }
            if (Program.GetFuncNum(scope, x1) < Program.GetFuncNum(scope, x2) && Program.GetFuncNum(scope, x1) < Program.GetFuncNum(scope, xm))
            {
                return x1;
            }
            else if (Program.GetFuncNum(scope, x2) < Program.GetFuncNum(scope, x1) && Program.GetFuncNum(scope, x2) < Program.GetFuncNum(scope, xm))
            {
                return x2;
            }
            else
            {
                return xm;
            }
        }



        public static double[] CalSPSABB(int scope, double[] x, int l)
        {
            double[] x1 = new double[scope + 2];
            double[] x2 = new double[scope + 2];
            double[] delta = BernoulliGenerator(scope);
            for (int i = 0; i < scope + 2; i++)
            {
                if (delta[i] == 0)
                {
                    delta[i] = -1;
                }
                x1[i] = x[i] + c / Math.Pow(l + 1, gamma) * delta[i];
                x2[i] = x[i] - c / Math.Pow(l + 1, gamma) * delta[i];
            }

            double[] f1 = new double[scope + 2];
            double[] f2 = new double[scope + 2];
            for (int i = 0; i < scope + 2; i++)
            {
                f1[i] = Program.GetFuncNumSingle(scope, x, x1[i], i) - Program.GetFuncNum(scope, x);
                f2[i] = Program.GetFuncNumSingle(scope, x, x2[i], i) - Program.GetFuncNum(scope, x);
            }
            double[] g = new double[scope + 2];
            for (int i = 0; i < scope + 2; i++)
            {
                g[i] = (f1[i] - f2[i]) / (2 * c / Math.Pow(l + 1, gamma) * delta[i]);
            }
            double[] xm = new double[scope + 2];

            //BB


            double[] df = Program.GetFuncDen(scope, x);
            double p = 0; double q = 0;
            for (int i = 0; i < scope + 2; i++)
            {
                p += (x[i] - x_1[i]) * (df[i] - df1[i]);
                q += Math.Pow(df[i] - df1[i], 2);
            }
            p /= q;
            for (int i = 0; i < scope + 2; i++)
            {
                x_1[i] = x[i];
            }
            for (int i = 0; i < scope + 2; i++)
            {
                df1[i] = df[i];
            }



            for (int i = 0; i < scope + 2; i++)
            {
                xm[i] = x[i] - p * g[i];
            }
            double f_des1;
            double f_des = Program.GetFuncNum(scope, xm);
            for (int i = 0; i < search; i++)
            {
                for (int j = 0; j < scope + 2; j++)
                {
                    xm[j] = x[j] - p * g[j] / 2;
                }
                f_des1 = Program.GetFuncNum(scope, xm);
                if (f_des1 < f_des)
                {
                    break;
                }
            }


            if (Program.GetFuncNum(scope, x1) < Program.GetFuncNum(scope, x2) && Program.GetFuncNum(scope, x1) < Program.GetFuncNum(scope, xm))
            {
                return x1;
            }
            else if (Program.GetFuncNum(scope, x2) < Program.GetFuncNum(scope, x1) && Program.GetFuncNum(scope, x2) < Program.GetFuncNum(scope, xm))
            {
                return x2;
            }
            else
            {
                return xm;
            }
        }

        public static double[] df1 = new double[scope + 2];
        public static double[] x_1 = new double[scope + 2];

    }
}
