using MathNet.Numerics.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using static System.Formats.Asn1.AsnWriter;

namespace SPSA
{
    public class ZNewt
    {
        static int iterZ = 100;
        //移动倍率
        static int scopeZ = 2;
        //黄金分割率
        static double ratioZ = 0.618;
        //精度
        static double tolZ = 0.0001;
        public static double[] CalZNewt(int scope, double[] x)
        {
            double[] den = Program.GetFuncDen(scope, x);
            double[,] den2 = Program.GetFuncDen2(scope, x);
            Matrix<double> hessian = Matrix<double>.Build.DenseOfArray(den2);
            Matrix<double> hessian_inverse = hessian.Inverse();
            double[,] den2_inverse = hessian_inverse.ToArray();
            double[] xm = new double[scope + 2];
            for (int i = 0; i < scope + 2; i++)
            {
                for (int j = 0; j < scope + 2; j++)
                {
                    xm[i] += den2_inverse[i, j] * den[j];
                }
            }
            return Z(xm, x, scope);
        }

        //黄金分割法求区间
        static double[] Z(double[] xm, double[] x, int scope)
        {
            double[] a = new double[xm.Length];
            double[] b = new double[xm.Length];
            double[] x1 = new double[xm.Length];
            double[] x2 = new double[xm.Length];
            for (int i = 0; i < xm.Length; i++)
            {
                a[i] = x[i] - scopeZ * xm[i];
                b[i] = x[i];
                x1[i] = b[i] - ratioZ * (b[i] - a[i]);
                x2[i] = a[i] + ratioZ * (b[i] - a[i]);
            }
            double f1 = Program.GetFuncNum(scope, x1);
            double f2 = Program.GetFuncNum(scope, x2);

            for (int i = 0; i < iterZ; i++)
            {
                if (f1 > f2)
                {
                    for (int k = 0; k < xm.Length; k++)
                    {
                        a[k] = x1[k]; x1[k] = x2[k];
                        x2[k] = a[k] + ratioZ * (b[k] - a[k]);
                    }
                    f1 = f2;
                    f2 = Program.GetFuncNum(scope, x2);
                }
                else
                {
                    for (int k = 0; k < xm.Length; k++)
                    {
                        b[k] = x2[k]; x2[k] = x1[k];
                        x1[k] = b[k] - ratioZ * (b[k] - a[k]);
                    }
                    f2 = f1;
                    f1 = Program.GetFuncNum(scope, x1);
                }
                double temp = 0;
                for (int j = 0; j < xm.Length; j++)
                {
                    if (Math.Abs(b[j] - a[j]) >= temp)
                    {
                        temp = Math.Abs(b[j] - a[j]);
                    }
                }
                if (temp < tolZ)
                {
                    break;
                }
            }

            double[] lambda = new double[xm.Length];
            for (int i = 0; i < xm.Length; i++)
            {
                lambda[i] = 0.5 * (a[i] + b[i]);
            }
            return lambda;
        }
    }
}
