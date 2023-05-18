using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using static System.Formats.Asn1.AsnWriter;

namespace SPSA
{
    public class ConjGrads
    {
        static double g0 = 0;
        static double[] d0 = new double[Program.scope + 2];

        //此函数用于迭代共轭梯度法中的x,选用的是Fletcher-Reeves公式
        public static double[] CalConjGrads(int scope, double[] x, int l)
        {
            double[] xm = new double[scope + 2];
            double[] g = new double[scope + 2];
            double[] d = new double[scope + 2];
            double[][] G = GetG(scope);
            for (int i = 0; i < scope + 2; i++)
            {
                for (int j = 0; j < scope + 2; j++)
                {
                    g[i] += G[i][j] * x[j];
                    if (l == 0)
                    {
                        d[i] -= G[i][j] * x[j];
                    }
                }
            }
            double p = 0; double q = 0; double[] q1 = new double[scope + 2];
            if (l != 0)
            {
                double b = 0;
                for (int i = 0; i < scope + 2; i++)
                {
                    p += g[i] * g[i];
                }
                b = p / g0; g0 = p;
                for (int i = 0; i < scope + 2; i++)
                {
                    d[i] = -g[i] + b * d0[i];
                }
                for (int i = 0; i < scope + 2; i++)
                {
                    for (int j = 0; j < scope + 2; j++)
                    {
                        q1[i] += d[j] * G[j][i];
                    }
                    q += d[i] * q1[i];
                }
                p /= q;
                for (int i = 0; i < scope + 2; i++)
                {
                    d0[i] = d[i];
                }
                for (int i = 0; i < scope + 2; i++)
                {
                    xm[i] = x[i] + p * d[i];
                }
            }
            else
            {
                for (int i = 0; i < scope + 2; i++)
                {
                    p += g[i] * g[i];
                    for (int j = 0; j < scope + 2; j++)
                    {
                        q1[i] += d[j] * G[j][i];
                    }
                    q += d[i] * q1[i];
                }
                g0 = p;
                p /= q;
                for (int i = 0; i < scope + 2; i++)
                {
                    d0[i] = d[i];
                }
                for (int i = 0; i < scope + 2; i++)
                {
                    xm[i] = x[i] + p * d[i];
                }
            }
            return xm;
        }

        static double[][] GetG(int scope)
        {
            double[][] G = new double[scope + 2][];
            //构造矩阵G
            for (int i = 0; i < scope + 2; i++)
            {
                G[i] = new double[scope + 2];
                for (int j = 0; j < scope + 2; j++)
                {
                    if (i == 0)
                    {
                        G[i][0] = 6; G[i][1] = -2; G[i][2] = -2;
                    }
                    else if (i == 1)
                    {
                        G[i][0] = -2; G[i][1] = 12; G[i][2] = -4; G[i][2] = -2;
                    }
                    else if (i == scope)
                    {
                        G[i][scope - 2] = -2; G[i][scope - 1] = -4; G[i][scope] = 12; G[i][scope + 1] = -2;
                    }
                    else if (i == scope + 1)
                    {
                        G[i][scope - 1] = -2; G[i][scope] = -2; G[i][scope + 1] = 6;
                    }
                    else
                    {
                        G[i][i] = 18; G[i][i - 1] = -4; G[i][i - 2] = -2; G[i][i + 1] = -4; G[i][i + 2] = -2;
                    }
                }
            }
            return G;
        }
    }
}
